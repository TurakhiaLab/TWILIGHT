#include "subtree-refine.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace msa {
namespace accurate {

namespace {

/**
 * Compute the self-alignment score of a sequence using localScoringMatrix.
 * Matches MAFFT's selfscore denominator calculation:
 *   selfscore[i] = sum_{k=0}^{L-1} matrix[c_k][c_k]
 */
float computeSequenceSelfScore(
    const SequenceDB::SequenceInfo* sInfo,
    char type,
    const Params& param)
{
    if (!sInfo) return 1.0f;
    const char* raw = sInfo->alnStorage[sInfo->storage];
    int len = sInfo->len;
    float score = 0.0f;
    for (int k = 0; k < len; ++k) {
        char c = raw[k];
        if (c == '-' || c == '.') continue;
        int idx = letterIdx(type, c);
        if (idx >= 0 && idx < param.matrixSize) {
            score += param.localScoringMatrix[idx][idx];
        }
    }
    return (score > 0.0f) ? score : static_cast<float>(std::max(1, len));
}

struct UpgmaCluster {
    phylogeny::Node* node = nullptr;
    int size = 1;
    float height = 0.0f;
    bool active = true;
};

/**
 * Build guide subtree using UPGMA (Average Linkage).
 */
phylogeny::Node* buildSubtreeUPGMA(
    int N,
    const std::string& rootID,
    std::vector<UpgmaCluster>& clusters,
    std::vector<std::vector<float>>& distMatrix)
{
    int numActive = N;
    int internalCounter = 1;

    while (numActive > 1) {
        // Find closest pair of active clusters (u, v) with u < v
        int bestU = -1, bestV = -1;
        float minDist = std::numeric_limits<float>::infinity();

        for (int i = 0; i < N; ++i) {
            if (!clusters[i].active) continue;
            for (int j = i + 1; j < N; ++j) {
                if (!clusters[j].active) continue;
                if (distMatrix[i][j] < minDist) {
                    minDist = distMatrix[i][j];
                    bestU = i;
                    bestV = j;
                }
            }
        }

        if (bestU == -1 || bestV == -1) break;

        // Determine node identifier: the final merged node takes rootID
        std::string nodeID = (numActive == 2)
            ? rootID
            : ("node_rf_" + std::to_string(internalCounter++));

        phylogeny::Node* parentNode = new phylogeny::Node(nodeID, 0.0f);
        parentNode->grpID = -1;

        float mergeHeight = minDist * 0.5f;
        // Ensure ultrametric monotonicity
        if (mergeHeight < clusters[bestU].height) mergeHeight = clusters[bestU].height;
        if (mergeHeight < clusters[bestV].height) mergeHeight = clusters[bestV].height;

        clusters[bestU].node->parent = parentNode;
        clusters[bestV].node->parent = parentNode;
        clusters[bestU].node->branchLength = std::max(1e-5f, mergeHeight - clusters[bestU].height);
        clusters[bestV].node->branchLength = std::max(1e-5f, mergeHeight - clusters[bestV].height);

        parentNode->children.push_back(clusters[bestU].node);
        parentNode->children.push_back(clusters[bestV].node);

        if (numActive == 2) {
            clusters[bestU].node = parentNode;
            clusters[bestV].active = false;
            break;
        }

        // Update average distances from new cluster (placed at bestU) to other active clusters
        int sizeU = clusters[bestU].size;
        int sizeV = clusters[bestV].size;
        int newSize = sizeU + sizeV;

        for (int k = 0; k < N; ++k) {
            if (!clusters[k].active || k == bestU || k == bestV) continue;
            float newDist = (sizeU * distMatrix[bestU][k] + sizeV * distMatrix[bestV][k]) / static_cast<float>(newSize);
            distMatrix[bestU][k] = distMatrix[k][bestU] = newDist;
        }

        clusters[bestU].node = parentNode;
        clusters[bestU].size = newSize;
        clusters[bestU].height = mergeHeight;
        clusters[bestV].active = false;

        --numActive;
    }

    // Locate the root node
    for (int i = 0; i < N; ++i) {
        if (clusters[i].active) {
            clusters[i].node->identifier = rootID;
            clusters[i].node->branchLength = 0.0f;
            clusters[i].node->parent = nullptr;
            return clusters[i].node;
        }
    }
    return nullptr;
}

struct NjCluster {
    phylogeny::Node* node = nullptr;
    bool active = true;
};

/**
 * Build guide subtree using Neighbor-Joining (NJ).
 */
phylogeny::Node* buildSubtreeNJ(
    int N,
    const std::string& rootID,
    std::vector<NjCluster>& clusters,
    std::vector<std::vector<float>>& distMatrix)
{
    int numActive = N;
    int internalCounter = 1;

    while (numActive > 2) {
        // Calculate net divergence R(i) for each active cluster
        std::vector<float> R(N, 0.0f);
        for (int i = 0; i < N; ++i) {
            if (!clusters[i].active) continue;
            float sumD = 0.0f;
            for (int j = 0; j < N; ++j) {
                if (i != j && clusters[j].active) {
                    sumD += distMatrix[i][j];
                }
            }
            R[i] = sumD / static_cast<float>(numActive - 2);
        }

        // Find pair (u, v) minimizing Q(u, v) = D(u, v) - R(u) - R(v)
        int bestU = -1, bestV = -1;
        float minQ = std::numeric_limits<float>::infinity();

        for (int i = 0; i < N; ++i) {
            if (!clusters[i].active) continue;
            for (int j = i + 1; j < N; ++j) {
                if (!clusters[j].active) continue;
                float Q = distMatrix[i][j] - R[i] - R[j];
                if (Q < minQ) {
                    minQ = Q;
                    bestU = i;
                    bestV = j;
                }
            }
        }

        if (bestU == -1 || bestV == -1) break;

        std::string nodeID = "node_rf_" + std::to_string(internalCounter++);
        phylogeny::Node* parentNode = new phylogeny::Node(nodeID, 0.0f);
        parentNode->grpID = -1;

        float dUV = distMatrix[bestU][bestV];
        float brU = 0.5f * dUV + 0.5f * (R[bestU] - R[bestV]);
        float brV = dUV - brU;
        if (brU < 1e-5f) brU = 1e-5f;
        if (brV < 1e-5f) brV = 1e-5f;

        clusters[bestU].node->parent = parentNode;
        clusters[bestV].node->parent = parentNode;
        clusters[bestU].node->branchLength = brU;
        clusters[bestV].node->branchLength = brV;

        parentNode->children.push_back(clusters[bestU].node);
        parentNode->children.push_back(clusters[bestV].node);

        // Update distance matrix for parentNode placed at bestU
        for (int k = 0; k < N; ++k) {
            if (!clusters[k].active || k == bestU || k == bestV) continue;
            float newDist = 0.5f * (distMatrix[bestU][k] + distMatrix[bestV][k] - dUV);
            if (newDist < 1e-5f) newDist = 1e-5f;
            distMatrix[bestU][k] = distMatrix[k][bestU] = newDist;
        }

        clusters[bestU].node = parentNode;
        clusters[bestV].active = false;
        --numActive;
    }

    // Join last two active clusters to root
    int u = -1, v = -1;
    for (int i = 0; i < N; ++i) {
        if (clusters[i].active) {
            if (u == -1) u = i;
            else { v = i; break; }
        }
    }

    if (u != -1 && v != -1) {
        phylogeny::Node* rootNode = new phylogeny::Node(rootID, 0.0f);
        rootNode->grpID = -1;
        float dUV = distMatrix[u][v];
        float brU = std::max(1e-5f, dUV * 0.5f);
        float brV = std::max(1e-5f, dUV * 0.5f);

        clusters[u].node->parent = rootNode;
        clusters[v].node->parent = rootNode;
        clusters[u].node->branchLength = brU;
        clusters[v].node->branchLength = brV;

        rootNode->children.push_back(clusters[u].node);
        rootNode->children.push_back(clusters[v].node);
        return rootNode;
    }

    return (u != -1) ? clusters[u].node : nullptr;
}

} // anonymous namespace

phylogeny::Tree* refineGroupSubtreeFromHomTable(
    HierarchyGroup& group,
    const LocalHomTable& localHomTable,
    phylogeny::Tree* oldSubT,
    SequenceDB* database,
    Option* option,
    const Params& param)
{
    if (!oldSubT || !oldSubT->root) {
        return oldSubT;
    }

    int N = localHomTable.getNumSeqs();
    if (N <= 1) {
        return oldSubT;
    }

    // 1. Compute selfScore for each sequence in group
    std::vector<float> selfScores(N, 1.0f);
    for (int i = 0; i < N; ++i) {
        int gIdx = localHomTable.getGlobalIndex(i);
        const auto* sInfo = (gIdx >= 0 && gIdx < static_cast<int>(database->sequences.size()))
            ? database->sequences[gIdx] : nullptr;
        selfScores[i] = computeSequenceSelfScore(sInfo, option->type, param);
    }

    // 2. Compute MAFFT score2dist distance matrix
    //    D(i, j) = 1.0 - Score(i, j) / min(selfScore_i, selfScore_j)
    std::vector<std::vector<float>> distMatrix(N, std::vector<float>(N, 0.0f));
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            float minSelf = std::min(selfScores[i], selfScores[j]);
            float rawScore = static_cast<float>(localHomTable.getPairResult(i, j).score);
            float sim = (minSelf > 0.0f) ? (rawScore / minSelf) : 0.0f;
            if (sim > 1.0f) sim = 1.0f;
            if (sim < 0.0f) sim = 0.0f;
            float dist = 1.0f - sim;
            distMatrix[i][j] = distMatrix[j][i] = dist;
        }
    }

    // Preserve original root identifier for seamless connection to upper hierarchy
    std::string rootID = oldSubT->root->identifier;
    if (rootID.empty()) {
        rootID = group.rootIdentifier.empty() ? "node_root" : group.rootIdentifier;
    }

    // 3. Initialize leaf nodes for each sequence
    phylogeny::Node* newRootNode = nullptr;

    if (option->refineMethod == "nj") {
        std::vector<NjCluster> clusters(N);
        for (int i = 0; i < N; ++i) {
            int gIdx = localHomTable.getGlobalIndex(i);
            const auto* sInfo = database->sequences[gIdx];
            phylogeny::Node* leaf = new phylogeny::Node(sInfo->name, 0.0f);
            leaf->grpID = -1;
            leaf->placed = false;
            leaf->seqsIncluded = { gIdx };
            leaf->alnLen = sInfo->len;
            leaf->alnNum = 1;
            leaf->alnWeight = sInfo->weight;
            clusters[i].node = leaf;
            clusters[i].active = true;
        }
        newRootNode = buildSubtreeNJ(N, rootID, clusters, distMatrix);
    } else {
        // Default: UPGMA
        std::vector<UpgmaCluster> clusters(N);
        for (int i = 0; i < N; ++i) {
            int gIdx = localHomTable.getGlobalIndex(i);
            const auto* sInfo = database->sequences[gIdx];
            phylogeny::Node* leaf = new phylogeny::Node(sInfo->name, 0.0f);
            leaf->grpID = -1;
            leaf->placed = false;
            leaf->seqsIncluded = { gIdx };
            leaf->alnLen = sInfo->len;
            leaf->alnNum = 1;
            leaf->alnWeight = sInfo->weight;
            clusters[i].node = leaf;
            clusters[i].size = 1;
            clusters[i].height = 0.0f;
            clusters[i].active = true;
        }
        newRootNode = buildSubtreeUPGMA(N, rootID, clusters, distMatrix);
    }

    if (!newRootNode) {
        return oldSubT;
    }

    // 4. Create new Tree and populate allNodes
    phylogeny::Tree* newTree = new phylogeny::Tree();
    newTree->root = newRootNode;
    newTree->m_numLeaves = N;

    // Collect all nodes into allNodes map
    std::vector<phylogeny::Node*> stack = { newTree->root };
    while (!stack.empty()) {
        phylogeny::Node* curr = stack.back();
        stack.pop_back();
        newTree->allNodes[curr->identifier] = curr;
        for (phylogeny::Node* child : curr->children) {
            stack.push_back(child);
        }
    }

    // Update tree metrics and sequence weights
    newTree->calLeafNum();
    phylogeny::updateLevels(newTree->root, 1);
    newTree->calSeqWeight();

    // 5. Clean up old subtree and return refined tree
    delete oldSubT;
    return newTree;
}

} // namespace accurate
} // namespace msa


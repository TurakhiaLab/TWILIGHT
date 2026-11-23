#include "subtree-refine.hpp"
#include "align-method.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <chrono>
#include <cctype>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_arena.h>
#include <tbb/enumerable_thread_specific.h>

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

namespace {

static void collectSubtreeLeaves(phylogeny::Node* node, std::vector<int>& outLeaves, SequenceDB* database) {
    if (!node) return;
    if (node->children.empty()) {
        if (!node->seqsIncluded.empty()) {
            outLeaves.push_back(node->seqsIncluded[0]);
        } else {
            auto it = database->name_map.find(node->identifier);
            if (it != database->name_map.end()) {
                outLeaves.push_back(it->second->id);
            }
        }
        return;
    }
    for (auto* ch : node->children) {
        collectSubtreeLeaves(ch, outLeaves, database);
    }
}

struct Bipartition {
    std::vector<int> groupA;
    std::vector<int> groupB;
};

static void collectCandidateBipartitions(
    phylogeny::Node* root,
    const std::vector<int>& allLeaves,
    SequenceDB* database,
    std::vector<Bipartition>& outBips)
{
    if (!root || allLeaves.size() <= 2) return;

    std::vector<phylogeny::Node*> stack = { root };
    std::vector<Bipartition> rawBips;

    while (!stack.empty()) {
        phylogeny::Node* curr = stack.back();
        stack.pop_back();

        if (curr != root) {
            std::vector<int> cladeLeaves;
            collectSubtreeLeaves(curr, cladeLeaves, database);
            std::sort(cladeLeaves.begin(), cladeLeaves.end());
            cladeLeaves.erase(std::unique(cladeLeaves.begin(), cladeLeaves.end()), cladeLeaves.end());

            if (!cladeLeaves.empty() && cladeLeaves.size() < allLeaves.size()) {
                std::unordered_set<int> aSet(cladeLeaves.begin(), cladeLeaves.end());
                std::vector<int> bLeaves;
                bLeaves.reserve(allLeaves.size() - cladeLeaves.size());
                for (int s : allLeaves) {
                    if (aSet.find(s) == aSet.end()) {
                        bLeaves.push_back(s);
                    }
                }

                if (!bLeaves.empty()) {
                    rawBips.push_back({std::move(cladeLeaves), std::move(bLeaves)});
                }
            }
        }

        for (auto* ch : curr->children) {
            stack.push_back(ch);
        }
    }

    // Filter duplicate partitions (bipartition {A, B} is identical to {B, A})
    std::set<std::vector<int>> seenA;
    for (auto& bip : rawBips) {
        std::sort(bip.groupA.begin(), bip.groupA.end());
        std::sort(bip.groupB.begin(), bip.groupB.end());
        // Canonicalize: smaller size or lexicographically smaller first
        if (bip.groupA.size() > bip.groupB.size() || 
           (bip.groupA.size() == bip.groupB.size() && bip.groupA > bip.groupB)) {
            std::swap(bip.groupA, bip.groupB);
        }
        if (seenA.insert(bip.groupA).second) {
            outBips.push_back(std::move(bip));
        }
    }
}

static IntPairVec findGroupGappyColumns(const std::vector<int>& seqs, int alnLen, SequenceDB* database) {
    IntPairVec gappy;
    int runStart = -1;
    for (int col = 0; col < alnLen; ++col) {
        bool allGap = true;
        for (int sIdx : seqs) {
            char c = database->sequences[sIdx]->alnStorage[database->sequences[sIdx]->storage][col];
            if (c != '-' && c != '.' && c != ' ') {
                allGap = false;
                break;
            }
        }
        if (allGap) {
            if (runStart == -1) runStart = col;
        } else {
            if (runStart != -1) {
                gappy.push_back({runStart, col - runStart});
                runStart = -1;
            }
        }
    }
    if (runStart != -1) {
        gappy.push_back({runStart, alnLen - runStart});
    }
    return gappy;
}

static float computeCurrentBipartitionWSP(
    const std::vector<int>& seqsA,
    const std::vector<int>& seqsB,
    int alnLen,
    SequenceDB* database,
    char type,
    const Params& param)
{
    if (alnLen <= 0 || seqsA.empty() || seqsB.empty()) return 0.0f;

    const int maxMatrixDim = (type == 'p') ? 21 : 5;
    const float gapOpen = param.gapOpen;
    const float gapExtend = param.gapExtend;

    std::vector<const char*> rawA(seqsA.size());
    for (size_t i = 0; i < seqsA.size(); ++i) {
        auto* s = database->sequences[seqsA[i]];
        rawA[i] = s->alnStorage[s->storage];
    }
    std::vector<const char*> rawB(seqsB.size());
    for (size_t i = 0; i < seqsB.size(); ++i) {
        auto* s = database->sequences[seqsB[i]];
        rawB[i] = s->alnStorage[s->storage];
    }

    std::vector<char> prevCharA(seqsA.size(), 'A');
    std::vector<char> prevCharB(seqsB.size(), 'A');

    float totalWSP = 0.0f;
    int countsA[21];
    int countsB[21];

    for (int col = 0; col < alnLen; ++col) {
        std::fill(countsA, countsA + maxMatrixDim, 0);
        std::fill(countsB, countsB + maxMatrixDim, 0);
        int numResA = 0, numResB = 0;
        int numGapOpenA = 0, numGapOpenB = 0;
        int numGapExtA = 0, numGapExtB = 0;

        for (size_t i = 0; i < seqsA.size(); ++i) {
            char ch = rawA[i][col];
            bool isGap = (ch == '-' || ch == '.' || ch == ' ');
            if (isGap) {
                bool prevIsGap = (col == 0) ? false : (prevCharA[i] == '-' || prevCharA[i] == '.' || prevCharA[i] == ' ');
                if (!prevIsGap) numGapOpenA++;
                else numGapExtA++;
            } else {
                numResA++;
                int aIdx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                if (aIdx >= 0 && aIdx < maxMatrixDim) countsA[aIdx]++;
            }
            prevCharA[i] = ch;
        }

        for (size_t i = 0; i < seqsB.size(); ++i) {
            char ch = rawB[i][col];
            bool isGap = (ch == '-' || ch == '.' || ch == ' ');
            if (isGap) {
                bool prevIsGap = (col == 0) ? false : (prevCharB[i] == '-' || prevCharB[i] == '.' || prevCharB[i] == ' ');
                if (!prevIsGap) numGapOpenB++;
                else numGapExtB++;
            } else {
                numResB++;
                int aIdx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                if (aIdx >= 0 && aIdx < maxMatrixDim) countsB[aIdx]++;
            }
            prevCharB[i] = ch;
        }

        // Inter-group substitution score
        float colSub = 0.0f;
        if (param.scoringMatrix) {
            for (int a = 0; a < maxMatrixDim; ++a) {
                if (countsA[a] == 0) continue;
                for (int b = 0; b < maxMatrixDim; ++b) {
                    if (countsB[b] > 0) {
                        colSub += static_cast<float>(countsA[a]) * countsB[b] * param.scoringMatrix[a][b];
                    }
                }
            }
        }

        // Inter-group gap score
        float colGap = static_cast<float>(numResA) * (numGapOpenB * gapOpen + numGapExtB * gapExtend) +
                       static_cast<float>(numResB) * (numGapOpenA * gapOpen + numGapExtA * gapExtend);

        totalWSP += (colSub + colGap);
    }

    return totalWSP;
}

static float computeBipartitionWSP(
    const std::vector<int>& seqsA,
    const std::vector<int>& seqsB,
    const std::vector<int>& colsA,
    const std::vector<int>& colsB,
    const alnPath& path,
    SequenceDB* database,
    char type,
    const Params& param)
{
    const int newLen = static_cast<int>(path.size());
    if (newLen <= 0 || seqsA.empty() || seqsB.empty()) return 0.0f;

    const int maxMatrixDim = (type == 'p') ? 21 : 5;
    const float gapOpen = param.gapOpen;
    const float gapExtend = param.gapExtend;

    std::vector<const char*> rawA(seqsA.size());
    for (size_t i = 0; i < seqsA.size(); ++i) {
        auto* s = database->sequences[seqsA[i]];
        rawA[i] = s->alnStorage[s->storage];
    }
    std::vector<const char*> rawB(seqsB.size());
    for (size_t i = 0; i < seqsB.size(); ++i) {
        auto* s = database->sequences[seqsB[i]];
        rawB[i] = s->alnStorage[s->storage];
    }

    int cA = 0;
    int cB = 0;
    std::vector<char> prevCharA(seqsA.size(), 'A');
    std::vector<char> prevCharB(seqsB.size(), 'A');

    float totalWSP = 0.0f;
    int countsA[21];
    int countsB[21];

    for (int k = 0; k < newLen; ++k) {
        int8_t code = path[k];
        bool hasA = (code == 0 || code == 2);
        bool hasB = (code == 0 || code == 1);

        int origColA = hasA ? colsA[cA++] : -1;
        int origColB = hasB ? colsB[cB++] : -1;

        std::fill(countsA, countsA + maxMatrixDim, 0);
        std::fill(countsB, countsB + maxMatrixDim, 0);
        int numResA = 0, numResB = 0;
        int numGapOpenA = 0, numGapOpenB = 0;
        int numGapExtA = 0, numGapExtB = 0;

        for (size_t i = 0; i < seqsA.size(); ++i) {
            char ch = (origColA >= 0) ? rawA[i][origColA] : '-';
            bool isGap = (ch == '-' || ch == '.' || ch == ' ');
            if (isGap) {
                bool prevIsGap = (k == 0) ? false : (prevCharA[i] == '-' || prevCharA[i] == '.' || prevCharA[i] == ' ');
                if (!prevIsGap) numGapOpenA++;
                else numGapExtA++;
            } else {
                numResA++;
                int aIdx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                if (aIdx >= 0 && aIdx < maxMatrixDim) countsA[aIdx]++;
            }
            prevCharA[i] = ch;
        }

        for (size_t i = 0; i < seqsB.size(); ++i) {
            char ch = (origColB >= 0) ? rawB[i][origColB] : '-';
            bool isGap = (ch == '-' || ch == '.' || ch == ' ');
            if (isGap) {
                bool prevIsGap = (k == 0) ? false : (prevCharB[i] == '-' || prevCharB[i] == '.' || prevCharB[i] == ' ');
                if (!prevIsGap) numGapOpenB++;
                else numGapExtB++;
            } else {
                numResB++;
                int aIdx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                if (aIdx >= 0 && aIdx < maxMatrixDim) countsB[aIdx]++;
            }
            prevCharB[i] = ch;
        }

        // Inter-group substitution score
        float colSub = 0.0f;
        if (param.scoringMatrix) {
            for (int a = 0; a < maxMatrixDim; ++a) {
                if (countsA[a] == 0) continue;
                for (int b = 0; b < maxMatrixDim; ++b) {
                    if (countsB[b] > 0) {
                        colSub += static_cast<float>(countsA[a]) * countsB[b] * param.scoringMatrix[a][b];
                    }
                }
            }
        }

        // Inter-group gap score
        float colGap = static_cast<float>(numResA) * (numGapOpenB * gapOpen + numGapExtB * gapExtend) +
                       static_cast<float>(numResB) * (numGapOpenA * gapOpen + numGapExtA * gapExtend);

        totalWSP += (colSub + colGap);
    }

    return totalWSP;
}

static void applyRealignedBipartition(
    const std::vector<int>& groupA,
    const std::vector<int>& groupB,
    const std::vector<int>& colsA,
    const std::vector<int>& colsB,
    const alnPath& path,
    SequenceDB* database)
{
    const int newLen = static_cast<int>(path.size());

    // Update Group A sequences
    for (int sIdx : groupA) {
        auto* seq = database->sequences[sIdx];
        int storeFrom = seq->storage;
        int storeTo = 1 - storeFrom;
        seq->memCheck(newLen + 1);

        int cA = 0;
        for (int k = 0; k < newLen; ++k) {
            int8_t code = path[k];
            if (code == 0 || code == 2) {
                int origCol = colsA[cA++];
                seq->alnStorage[storeTo][k] = seq->alnStorage[storeFrom][origCol];
            } else {
                seq->alnStorage[storeTo][k] = '-';
            }
        }
        seq->alnStorage[storeTo][newLen] = '\0';
        seq->len = newLen;
        seq->changeStorage();
    }

    // Update Group B sequences
    for (int sIdx : groupB) {
        auto* seq = database->sequences[sIdx];
        int storeFrom = seq->storage;
        int storeTo = 1 - storeFrom;
        seq->memCheck(newLen + 1);

        int cB = 0;
        for (int k = 0; k < newLen; ++k) {
            int8_t code = path[k];
            if (code == 0 || code == 1) {
                int origCol = colsB[cB++];
                seq->alnStorage[storeTo][k] = seq->alnStorage[storeFrom][origCol];
            } else {
                seq->alnStorage[storeTo][k] = '-';
            }
        }
        seq->alnStorage[storeTo][newLen] = '\0';
        seq->len = newLen;
        seq->changeStorage();
    }
}

struct RefineScratch {
    std::vector<float> flatConsistency;
    std::vector<int> colToProfA;
    std::vector<int> colToProfB;
    std::vector<int> refPosToCol;
    std::vector<int> qryPosToCol;
    std::vector<int> touchedRef;
    std::vector<int> touchedQry;
    int maxSeqsAllocated = 0;
    int maxResAllocated = 0;

    void ensureCapacity(int numSeqs, int maxResLen, int currentAlnLen, int refLen, int qryLen) {
        size_t matrixSize = static_cast<size_t>(refLen) * qryLen;
        if (flatConsistency.size() < matrixSize) {
            flatConsistency.resize(matrixSize);
        }
        std::fill_n(flatConsistency.data(), matrixSize, 0.0f);

        if (colToProfA.size() < static_cast<size_t>(currentAlnLen + 64)) {
            colToProfA.assign(currentAlnLen + 64, -1);
            colToProfB.assign(currentAlnLen + 64, -1);
        }

        if (numSeqs > maxSeqsAllocated || maxResLen > maxResAllocated) {
            maxSeqsAllocated = std::max(maxSeqsAllocated, numSeqs + 32);
            maxResAllocated = std::max(maxResAllocated, maxResLen + 256);
            size_t mapSize = static_cast<size_t>(maxSeqsAllocated) * maxResAllocated;
            refPosToCol.assign(mapSize, -1);
            qryPosToCol.assign(mapSize, -1);
            touchedRef.clear();
            touchedQry.clear();
        }
    }
};

static tbb::enumerable_thread_specific<RefineScratch> refineScratchPool;

} // anonymous namespace

static double computeSubtreeInitialScore(
    const std::vector<int>& allLeaves,
    int alnLen,
    SequenceDB* database,
    char type,
    const Params& param,
    bool hasCons,
    float cWeight,
    std::shared_ptr<LocalHomTable> localHomTable,
    bool isTCS)
{
    if (alnLen <= 0 || allLeaves.empty()) return 0.0;

    const int maxMatrixDim = (type == 'p') ? 21 : 5;
    const float gapOpen = param.gapOpen;
    const float gapExtend = param.gapExtend;
    const size_t numSeqs = allLeaves.size();

    std::vector<const char*> rawSeq(numSeqs);
    for (size_t i = 0; i < numSeqs; ++i) {
        auto* s = database->sequences[allLeaves[i]];
        rawSeq[i] = s->alnStorage[s->storage];
    }

    std::vector<char> prevChar(numSeqs, 'A');
    double totalWSP = 0.0;
    int counts[21];

    for (int col = 0; col < alnLen; ++col) {
        std::fill(counts, counts + maxMatrixDim, 0);
        int numRes = 0;
        int numGapOpen = 0;
        int numGapExt = 0;

        for (size_t i = 0; i < numSeqs; ++i) {
            char ch = rawSeq[i][col];
            bool isGap = (ch == '-' || ch == '.' || ch == ' ');
            if (isGap) {
                bool prevIsGap = (col == 0) ? false : (prevChar[i] == '-' || prevChar[i] == '.' || prevChar[i] == ' ');
                if (!prevIsGap) numGapOpen++;
                else numGapExt++;
            } else {
                numRes++;
                int aIdx = letterIdx(type, static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
                if (aIdx >= 0 && aIdx < maxMatrixDim) counts[aIdx]++;
            }
            prevChar[i] = ch;
        }

        // Intra-subtree substitution score across all pairs (i < j)
        double colSub = 0.0;
        if (param.scoringMatrix) {
            for (int a = 0; a < maxMatrixDim; ++a) {
                if (counts[a] == 0) continue;
                if (counts[a] > 1) {
                    colSub += ((static_cast<double>(counts[a]) * (counts[a] - 1)) / 2.0) * param.scoringMatrix[a][a];
                }
                for (int b = a + 1; b < maxMatrixDim; ++b) {
                    if (counts[b] > 0) {
                        colSub += (static_cast<double>(counts[a]) * counts[b]) * param.scoringMatrix[a][b];
                    }
                }
            }
        }

        // Intra-subtree gap score across all pairs (i < j)
        double colGap = static_cast<double>(numRes) * (numGapOpen * gapOpen + numGapExt * gapExtend);
        totalWSP += (colSub + colGap);
    }

    double totalCons = 0.0;
    if (hasCons && localHomTable) {
        // Map residue positions to alignment columns for all sequences in allLeaves
        std::vector<std::vector<int>> resToCol(numSeqs);
        for (size_t i = 0; i < numSeqs; ++i) {
            int localIdx = localHomTable->getLocalIndex(allLeaves[i]);
            if (localIdx < 0) continue;
            resToCol[i].reserve(database->sequences[allLeaves[i]]->len);
            for (int col = 0; col < alnLen; ++col) {
                char ch = rawSeq[i][col];
                if (ch != '-' && ch != '.' && ch != ' ') {
                    resToCol[i].push_back(col);
                }
            }
        }

        for (size_t i = 0; i < numSeqs; ++i) {
            int localA = localHomTable->getLocalIndex(allLeaves[i]);
            if (localA < 0) continue;
            const auto& rMapA = resToCol[i];

            for (size_t j = i + 1; j < numSeqs; ++j) {
                int localB = localHomTable->getLocalIndex(allLeaves[j]);
                if (localB < 0) continue;
                const auto& rMapB = resToCol[j];

                if (isTCS) {
                    const auto& tcsRes = localHomTable->getTCSPairResult(localA, localB);
                    for (const auto& entry : tcsRes.entries) {
                        int posA = (localA < localB) ? entry.pos1 : entry.pos2;
                        int posB = (localA < localB) ? entry.pos2 : entry.pos1;
                        if (posA < static_cast<int>(rMapA.size()) && posB < static_cast<int>(rMapB.size())) {
                            if (rMapA[posA] == rMapB[posB]) {
                                totalCons += entry.score;
                            }
                        }
                    }
                } else {
                    const auto& pairRes = localHomTable->getPairResult(localA, localB);
                    for (const auto& seg : pairRes.segments) {
                        if (seg.importance <= 0.0f) continue;
                        int segLen = seg.end1 - seg.start1 + 1;
                        for (int k = 0; k < segLen; ++k) {
                            int posA = (localA < localB) ? (seg.start1 + k) : (seg.start2 + k);
                            int posB = (localA < localB) ? (seg.start2 + k) : (seg.start1 + k);
                            if (posA < static_cast<int>(rMapA.size()) && posB < static_cast<int>(rMapB.size())) {
                                if (rMapA[posA] == rMapB[posB]) {
                                    totalCons += seg.importance;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return totalWSP + (hasCons ? (cWeight * totalCons) : 0.0);
}

void iterativeRefineGroupSubtree(
    HierarchyGroup& group,
    phylogeny::Tree* subT,
    std::shared_ptr<LocalHomTable> localHomTable,
    SequenceDB* database,
    Option* option,
    const Params& param,
    int teamThreads)
{
    if (!subT || !subT->root || subT->root->numLeaves <= 2) return;

    std::vector<int> allLeaves;
    collectSubtreeLeaves(subT->root, allLeaves, database);
    std::sort(allLeaves.begin(), allLeaves.end());
    allLeaves.erase(std::unique(allLeaves.begin(), allLeaves.end()), allLeaves.end());
    if (allLeaves.size() <= 2) return;

    int currentAlnLen = subT->root->alnLen;
    if (currentAlnLen <= 0) {
        currentAlnLen = database->sequences[allLeaves[0]]->len;
        subT->root->alnLen = currentAlnLen;
    }

    // 0. Initial cleanup of any columns that are entirely gap across all sequences in this subtree
    IntPairVec subtreeAllGaps = findGroupGappyColumns(allLeaves, currentAlnLen, database);
    if (!subtreeAllGaps.empty()) {
        std::vector<bool> isAllGapCol(currentAlnLen, false);
        for (const auto& g : subtreeAllGaps) {
            for (int k = 0; k < g.second; ++k) {
                isAllGapCol[g.first + k] = true;
            }
        }
        int newAlnLen = 0;
        for (int col = 0; col < currentAlnLen; ++col) {
            if (!isAllGapCol[col]) newAlnLen++;
        }
        for (int sIdx : allLeaves) {
            auto* seq = database->sequences[sIdx];
            int storeFrom = seq->storage;
            int storeTo = 1 - storeFrom;
            seq->memCheck(newAlnLen + 1);
            int wCol = 0;
            for (int col = 0; col < currentAlnLen; ++col) {
                if (!isAllGapCol[col]) {
                    seq->alnStorage[storeTo][wCol++] = seq->alnStorage[storeFrom][col];
                }
            }
            seq->alnStorage[storeTo][newAlnLen] = '\0';
            seq->len = newAlnLen;
            seq->changeStorage();
        }
        currentAlnLen = newAlnLen;
        subT->root->alnLen = currentAlnLen;
    }

    std::vector<Bipartition> candidateBips;
    collectCandidateBipartitions(subT->root, allLeaves, database, candidateBips);
    if (candidateBips.empty()) return;

    // Sort candidate bipartitions by balance descending (major balanced backbone cuts first, leaf twigs last)
    std::sort(candidateBips.begin(), candidateBips.end(), [](const Bipartition& a, const Bipartition& b) {
        int balA = std::min(static_cast<int>(a.groupA.size()), static_cast<int>(a.groupB.size()));
        int balB = std::min(static_cast<int>(b.groupA.size()), static_cast<int>(b.groupB.size()));
        if (balA != balB) return balA > balB;
        return a.groupA.size() > b.groupA.size();
    });

    int maxSubtreeSeqs = localHomTable ? localHomTable->getNumSeqs() : static_cast<int>(allLeaves.size());
    int maxSubtreeRes = 0;
    for (int sIdx : allLeaves) {
        if (database->sequences[sIdx]->len > maxSubtreeRes) {
            maxSubtreeRes = database->sequences[sIdx]->len;
        }
    }

    int maxIters = (option->refineIters > 0) ? option->refineIters : 1;
    int batchSize = (option->refineBatchSize > 0) ? option->refineBatchSize : ((teamThreads > 0) ? teamThreads : 4);

    struct CandidateResult {
        int branchIdx = -1;
        float deltaScore = 0.0f;
        alnPath path;
        std::vector<int> colsA;
        std::vector<int> colsB;
    };

    auto refineStart = std::chrono::high_resolution_clock::now();
    bool hasCons = (localHomTable && option->accurate);
    bool isTCS = (option->accuracyLevel == 3 || option->accurateMethod == "tcs");
    float cWeight = (hasCons ? option->consistencyWeight : 0.0f);
    double initialSubtreeScore = computeSubtreeInitialScore(
        allLeaves, currentAlnLen, database, option->type, param, hasCons, cWeight, localHomTable, isTCS
    );

    int actualIters = 0;
    int totalUpdatesApplied = 0;
    float cumulativeImprovement = 0.0f;

    std::atomic<int64_t> totalConsUs{0};
    std::atomic<int64_t> totalDpUs{0};
    std::atomic<int64_t> totalWspUs{0};

    for (int iter = 0; iter < maxIters; ++iter) {
        actualIters++;
        int updatesInThisIter = 0;
        float iterImprovement = 0.0f;

        for (size_t bStart = 0; bStart < candidateBips.size(); bStart += batchSize) {
            size_t bEnd = std::min(candidateBips.size(), bStart + batchSize);
            size_t curBatchCount = bEnd - bStart;
            std::vector<CandidateResult> batchResults(curBatchCount);

            tbb::this_task_arena::isolate([&] {
                tbb::parallel_for(tbb::blocked_range<size_t>(0, curBatchCount), [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t idx = range.begin(); idx < range.end(); ++idx) {
                        size_t b = bStart + idx;
                        const auto& bip = candidateBips[b];
                        IntPairVec gappyA = findGroupGappyColumns(bip.groupA, currentAlnLen, database);
                        IntPairVec gappyB = findGroupGappyColumns(bip.groupB, currentAlnLen, database);

                        int totalGappyA = 0;
                        for (const auto& g : gappyA) totalGappyA += g.second;
                        int totalGappyB = 0;
                        for (const auto& g : gappyB) totalGappyB += g.second;

                        IntPair lens = {currentAlnLen - totalGappyA, currentAlnLen - totalGappyB};
                        if (lens.first <= 0 || lens.second <= 0) continue;

                        // Build column mappings for Group A and Group B
                        std::vector<int> colsA;
                        colsA.reserve(lens.first);
                        int gcIdxA = 0;
                        for (int col = 0; col < currentAlnLen; ++col) {
                            if (gcIdxA < static_cast<int>(gappyA.size()) && col == gappyA[gcIdxA].first) {
                                col += gappyA[gcIdxA].second - 1;
                                gcIdxA++;
                            } else {
                                colsA.push_back(col);
                            }
                        }

                        std::vector<int> colsB;
                        colsB.reserve(lens.second);
                        int gcIdxB = 0;
                        for (int col = 0; col < currentAlnLen; ++col) {
                            if (gcIdxB < static_cast<int>(gappyB.size()) && col == gappyB[gcIdxB].first) {
                                col += gappyB[gcIdxB].second - 1;
                                gcIdxB++;
                            } else {
                                colsB.push_back(col);
                            }
                        }

                        phylogeny::Node nodeA("groupA", 0.0f);
                        nodeA.seqsIncluded = bip.groupA;
                        nodeA.alnNum = bip.groupA.size();
                        nodeA.alnLen = currentAlnLen;
                        nodeA.alnWeight = static_cast<float>(bip.groupA.size());

                        phylogeny::Node nodeB("groupB", 0.0f);
                        nodeB.seqsIncluded = bip.groupB;
                        nodeB.alnNum = bip.groupB.size();
                        nodeB.alnLen = currentAlnLen;
                        nodeB.alnWeight = static_cast<float>(bip.groupB.size());

                        NodePair nodePair(&nodeA, &nodeB);

                        auto tC0 = std::chrono::high_resolution_clock::now();
                        RefineScratch& scratch = refineScratchPool.local();
                        bool hasCons = (localHomTable && option->accurate);
                        float curCons = 0.0f;
                        if (hasCons) {
                            scratch.ensureCapacity(maxSubtreeSeqs, maxSubtreeRes, currentAlnLen, lens.first, lens.second);

                            for (int cA = 0; cA < lens.first; ++cA) scratch.colToProfA[colsA[cA]] = cA;
                            for (int cB = 0; cB < lens.second; ++cB) scratch.colToProfB[colsB[cB]] = cB;

                            for (int sIdx : bip.groupA) {
                                int localA = localHomTable->getLocalIndex(sIdx);
                                if (localA < 0) continue;
                                const auto* seq = database->sequences[sIdx];
                                const char* str = seq->alnStorage[seq->storage];
                                int resIdx = 0;
                                for (int col = 0; col < currentAlnLen; ++col) {
                                    char ch = str[col];
                                    if (ch != '-' && ch != '.') {
                                        int cA = scratch.colToProfA[col];
                                        if (cA >= 0) {
                                            int mIdx = localA * scratch.maxResAllocated + resIdx;
                                            scratch.refPosToCol[mIdx] = cA;
                                            scratch.touchedRef.push_back(mIdx);
                                        }
                                        resIdx++;
                                    }
                                }
                            }

                            for (int sIdx : bip.groupB) {
                                int localB = localHomTable->getLocalIndex(sIdx);
                                if (localB < 0) continue;
                                const auto* seq = database->sequences[sIdx];
                                const char* str = seq->alnStorage[seq->storage];
                                int resIdx = 0;
                                for (int col = 0; col < currentAlnLen; ++col) {
                                    char ch = str[col];
                                    if (ch != '-' && ch != '.') {
                                        int cB = scratch.colToProfB[col];
                                        if (cB >= 0) {
                                            int mIdx = localB * scratch.maxResAllocated + resIdx;
                                            scratch.qryPosToCol[mIdx] = cB;
                                            scratch.touchedQry.push_back(mIdx);
                                        }
                                        resIdx++;
                                    }
                                }
                            }

                            bool isTCS = (option->accuracyLevel == 3 || option->accurateMethod == "tcs");
                            for (int sA : bip.groupA) {
                                int localA = localHomTable->getLocalIndex(sA);
                                if (localA < 0) continue;
                                for (int sB : bip.groupB) {
                                    int localB = localHomTable->getLocalIndex(sB);
                                    if (localB < 0) continue;

                                    if (isTCS) {
                                        const auto& tcsRes = localHomTable->getTCSPairResult(localA, localB);
                                        if (tcsRes.entries.empty()) continue;

                                        if (localA < localB) {
                                            for (const auto& entry : tcsRes.entries) {
                                                int posA = entry.pos1;
                                                int posB = entry.pos2;
                                                int idxA = localA * scratch.maxResAllocated + posA;
                                                int idxB = localB * scratch.maxResAllocated + posB;
                                                int cA = scratch.refPosToCol[idxA];
                                                int cB = scratch.qryPosToCol[idxB];
                                                if (cA >= 0 && cB >= 0) {
                                                    scratch.flatConsistency[cA * lens.second + cB] += entry.score;
                                                }
                                            }
                                        } else {
                                            for (const auto& entry : tcsRes.entries) {
                                                int posB = entry.pos1;
                                                int posA = entry.pos2;
                                                int idxA = localA * scratch.maxResAllocated + posA;
                                                int idxB = localB * scratch.maxResAllocated + posB;
                                                int cA = scratch.refPosToCol[idxA];
                                                int cB = scratch.qryPosToCol[idxB];
                                                if (cA >= 0 && cB >= 0) {
                                                    scratch.flatConsistency[cA * lens.second + cB] += entry.score;
                                                }
                                            }
                                        }
                                    } else {
                                        const auto& pairRes = localHomTable->getPairResult(localA, localB);
                                        if (pairRes.segments.empty()) continue;

                                        if (localA < localB) {
                                            for (const auto& seg : pairRes.segments) {
                                                float imp = seg.importance;
                                                if (imp <= 0.0f) continue;
                                                int segLen = seg.end1 - seg.start1 + 1;
                                                for (int k = 0; k < segLen; ++k) {
                                                    int posA = seg.start1 + k;
                                                    int posB = seg.start2 + k;
                                                    int idxA = localA * scratch.maxResAllocated + posA;
                                                    int idxB = localB * scratch.maxResAllocated + posB;
                                                    int cA = scratch.refPosToCol[idxA];
                                                    int cB = scratch.qryPosToCol[idxB];
                                                    if (cA >= 0 && cB >= 0) {
                                                        scratch.flatConsistency[cA * lens.second + cB] += imp;
                                                    }
                                                }
                                            }
                                        } else {
                                            for (const auto& seg : pairRes.segments) {
                                                float imp = seg.importance;
                                                if (imp <= 0.0f) continue;
                                                int segLen = seg.end1 - seg.start1 + 1;
                                                for (int k = 0; k < segLen; ++k) {
                                                    int posB = seg.start1 + k;
                                                    int posA = seg.start2 + k;
                                                    int idxA = localA * scratch.maxResAllocated + posA;
                                                    int idxB = localB * scratch.maxResAllocated + posB;
                                                    int cA = scratch.refPosToCol[idxA];
                                                    int cB = scratch.qryPosToCol[idxB];
                                                    if (cA >= 0 && cB >= 0) {
                                                        scratch.flatConsistency[cA * lens.second + cB] += imp;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            // Compute consistency score of the current bipartition alignment
                            for (int col = 0; col < currentAlnLen; ++col) {
                                int cA = scratch.colToProfA[col];
                                int cB = scratch.colToProfB[col];
                                if (cA >= 0 && cB >= 0) {
                                    curCons += scratch.flatConsistency[cA * lens.second + cB];
                                }
                            }

                            for (int cA = 0; cA < lens.first; ++cA) scratch.colToProfA[colsA[cA]] = -1;
                            for (int cB = 0; cB < lens.second; ++cB) scratch.colToProfB[colsB[cB]] = -1;
                            for (int mIdx : scratch.touchedRef) scratch.refPosToCol[mIdx] = -1;
                            scratch.touchedRef.clear();
                            for (int mIdx : scratch.touchedQry) scratch.qryPosToCol[mIdx] = -1;
                            scratch.touchedQry.clear();
                        }
                        auto tC1 = std::chrono::high_resolution_clock::now();
                        totalConsUs.fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(tC1 - tC0).count());

                        auto tDp0 = std::chrono::high_resolution_clock::now();
                        msa::AlnProfile profRef, profQry;
                        std::pair<IntPairVec, IntPairVec> gappyCols = {gappyA, gappyB};
                        alignment_helper::calculatePSGP(profRef, profQry, nodePair, database, option, param, lens, gappyCols);

                        const float* flatConsPtr = hasCons ? scratch.flatConsistency.data() : nullptr;
                        float cWeight = option->consistencyWeight;
                        float dpScore = 0.0f;
                        std::vector<int8_t> aln_wo_gc = alignProfile(profRef, profQry, param, nullptr, cWeight, -1, &dpScore, flatConsPtr, lens.second);
                        auto tDp1 = std::chrono::high_resolution_clock::now();
                        totalDpUs.fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(tDp1 - tDp0).count());

                        if (aln_wo_gc.empty()) continue;

                        auto tW0 = std::chrono::high_resolution_clock::now();
                        float curWSP = computeCurrentBipartitionWSP(bip.groupA, bip.groupB, currentAlnLen, database, option->type, param);
                        float candWSP = computeBipartitionWSP(bip.groupA, bip.groupB, colsA, colsB, aln_wo_gc, database, option->type, param);
                        float deltaWSP = candWSP - curWSP;

                        // Compute consistency score of the candidate realignment
                        float candCons = 0.0f;
                        if (hasCons) {
                            int cA = 0;
                            int cB = 0;
                            for (int8_t code : aln_wo_gc) {
                                if (code == 0) {
                                    if (cA < lens.first && cB < lens.second) {
                                        candCons += scratch.flatConsistency[cA * lens.second + cB];
                                    }
                                }
                                if (code == 0 || code == 2) cA++;
                                if (code == 0 || code == 1) cB++;
                            }
                        }
                        float deltaCons = candCons - curCons;

                        // Unified composite acceptance score: deltaWSP + consistencyWeight * deltaCons
                        float delta = deltaWSP + (hasCons ? (cWeight * deltaCons) : 0.0f);
                        auto tW1 = std::chrono::high_resolution_clock::now();
                        totalWspUs.fetch_add(std::chrono::duration_cast<std::chrono::microseconds>(tW1 - tW0).count());

                        if (delta > 0.0f) {
                            batchResults[idx] = {static_cast<int>(b), delta, std::move(aln_wo_gc), std::move(colsA), std::move(colsB)};
                        }
                    }
                });
            });

            int bestIdx = -1;
            float maxDelta = 0.0f;
            for (size_t idx = 0; idx < curBatchCount; ++idx) {
                if (batchResults[idx].deltaScore > maxDelta) {
                    maxDelta = batchResults[idx].deltaScore;
                    bestIdx = idx;
                }
            }

            if (bestIdx >= 0 && maxDelta > 1e-3f) {
                size_t winningBip = bStart + bestIdx;
                applyRealignedBipartition(
                    candidateBips[winningBip].groupA,
                    candidateBips[winningBip].groupB,
                    batchResults[bestIdx].colsA,
                    batchResults[bestIdx].colsB,
                    batchResults[bestIdx].path,
                    database
                );

                currentAlnLen = batchResults[bestIdx].path.size();
                subT->root->alnLen = currentAlnLen;
                iterImprovement += maxDelta;
                cumulativeImprovement += maxDelta;
                updatesInThisIter++;
                totalUpdatesApplied++;
            }
        }

        if (updatesInThisIter == 0) {
            break; // Early stop if an entire pass over all branches yielded no further improvements
        }
    }

    // Recompute subT->root->msaFreq from the refined sequences
    int profileSize = (option->type == 'n') ? 6 : 22;
    subT->root->msaFreq.assign(currentAlnLen, std::vector<float>(profileSize, 0.0f));
    for (int sIdx : allLeaves) {
        float w = database->sequences[sIdx]->weight;
        int storage = database->sequences[sIdx]->storage;
        for (int t = 0; t < currentAlnLen; ++t) {
            int letterIndex = letterIdx(option->type, toupper(database->sequences[sIdx]->alnStorage[storage][t]));
            if (letterIndex >= 0 && letterIndex < profileSize) {
                subT->root->msaFreq[t][letterIndex] += w;
            }
        }
    }
    subT->root->alnLen = currentAlnLen;

    auto refineEnd = std::chrono::high_resolution_clock::now();
    auto refineMs = std::chrono::duration_cast<std::chrono::milliseconds>(refineEnd - refineStart).count();

    double percentImp = 0.0;
    if (std::abs(initialSubtreeScore) > 1e-4) {
        percentImp = (cumulativeImprovement / std::abs(initialSubtreeScore)) * 100.0;
    }

    if (option->printDetail || cumulativeImprovement > 0.0f) {
        std::cerr << "[Iterative Refine Group " << group.groupID << "] Seqs=" << allLeaves.size() 
                  << " | Time=" << refineMs << " ms"
                  << " | Updates=" << totalUpdatesApplied
                  << " | Score Improvement=+" << std::fixed << std::setprecision(2) << percentImp << "%\n";
    }
}

} // namespace accurate
} // namespace msa


#include "panmanWriter.hpp"
#include "panmanCodec.hpp"
#include "panmanStructs.hpp"
#include "panman.capnp.h"

#include <capnp/message.h>
#include <capnp/serialize.h>
#include <kj/std/iostream.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/lzma.hpp>

#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cassert>
#include <mutex>
#include <chrono>

namespace mga {
namespace panman {

namespace {

struct MafGappedConsensus {
    std::string gapFree;
    std::vector<std::pair<uint32_t, uint32_t>> gaps; // (nucPosition, gapLength)
    std::vector<std::pair<int, int>> oldToNew; // oldPos -> (nucPosition, gapPos)
};

MafGappedConsensus parseMafGappedConsensus(const std::string& gapped) {
    MafGappedConsensus out;
    out.oldToNew.resize(gapped.size());
    size_t pendingGap = 0;
    for (size_t oldPos = 0; oldPos < gapped.size(); oldPos++) {
        if (gapped[oldPos] == '-') {
            pendingGap++;
            out.oldToNew[oldPos] = {(int)out.gapFree.size(), (int)pendingGap - 1};
        } else {
            if (pendingGap > 0) {
                out.gaps.emplace_back((uint32_t)out.gapFree.size(), (uint32_t)pendingGap);
                pendingGap = 0;
            }
            out.oldToNew[oldPos] = {(int)out.gapFree.size(), -1};
            out.gapFree.push_back(gapped[oldPos]);
        }
    }
    if (pendingGap > 0) {
        out.gaps.emplace_back((uint32_t)out.gapFree.size(), (uint32_t)pendingGap);
    }
    return out;
}

void collectNodesPreorder(phylogeny::Node* node, std::vector<phylogeny::Node*>& nodes) {
    if (!node) return;
    nodes.push_back(node);
    for (const auto& child : node->children) {
        collectNodesPreorder(child.get(), nodes);
    }
}

void collectNodesPostorder(phylogeny::Node* node, std::vector<phylogeny::Node*>& nodes) {
    if (!node) return;
    for (const auto& child : node->children) {
        collectNodesPostorder(child.get(), nodes);
    }
    nodes.push_back(node);
}

std::vector<::mga::panman::NucMutData> consolidateNucMutations(const std::vector<::mga::panman::NucMutData>& muts) {
    if (muts.empty()) return {};

    std::vector<::mga::panman::NucMutData> sortedMuts = muts;
    std::sort(sortedMuts.begin(), sortedMuts.end(), [](const ::mga::panman::NucMutData& a, const ::mga::panman::NucMutData& b) {
        if (a.primaryBlockId != b.primaryBlockId) return a.primaryBlockId < b.primaryBlockId;
        if (a.nucPosition != b.nucPosition) return a.nucPosition < b.nucPosition;
        if (a.nucGapExist != b.nucGapExist) return a.nucGapExist < b.nucGapExist;
        return a.nucGapPosition < b.nucGapPosition;
    });

    std::vector<::mga::panman::NucMutData> out;
    out.reserve(sortedMuts.size());

    size_t i = 0;
    while (i < sortedMuts.size()) {
        const auto& first = sortedMuts[i];
        uint8_t firstType = first.mutInfo & 0x07;
        size_t j = i + 1;

        // Group up to 6 contiguous mutations with identical type, block, and gap status
        while (j < sortedMuts.size() && (j - i) < 6) {
            const auto& curr = sortedMuts[j];
            uint8_t currType = curr.mutInfo & 0x07;

            if (curr.primaryBlockId != first.primaryBlockId || currType != firstType || curr.nucGapExist != first.nucGapExist) {
                break;
            }

            if (!first.nucGapExist) {
                if (curr.nucPosition != first.nucPosition + static_cast<int32_t>(j - i)) {
                    break;
                }
            } else {
                if (curr.nucPosition != first.nucPosition || 
                    curr.nucGapPosition != first.nucGapPosition + static_cast<int32_t>(j - i)) {
                    break;
                }
            }
            j++;
        }

        size_t L = j - i;
        uint32_t payload = 0;
        for (size_t k = 0; k < L; ++k) {
            uint32_t code = (sortedMuts[i + k].mutInfo >> 8) & 0x0F;
            payload |= (code << (8 + 4 * k));
        }

        ::mga::panman::NucMutData consolidated;
        consolidated.primaryBlockId = first.primaryBlockId;
        consolidated.nucPosition = first.nucPosition;
        consolidated.nucGapPosition = first.nucGapPosition;
        consolidated.nucGapExist = first.nucGapExist;
        consolidated.mutInfo = payload | (static_cast<uint32_t>(L) << 4) | (firstType & 0x07);

        out.push_back(consolidated);
        i = j;
    }

    return out;
}

} // namespace

bool PanmanWriter::writePanMAN(BlockSet* rootBlockSet, 
                            Tree* tree, 
                            const std::string& newickString, 
                            const std::string& outputPath) {
    if (!rootBlockSet || !tree || !tree->root) {
        std::cerr << "[PanmanWriter Error] Null rootBlockSet or tree provided.\n";
        return false;
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    std::cout << "[PanmanWriter] Building PanMAN backbone from linearGraph & running TBB Fitch parsimony passes..." << std::endl;

    // 1. Get Linear Block Sequence (Backbone) & Audit VBlocks
    const auto& linearBlocks = rootBlockSet->getLinearizeBlocks();
    size_t numInstances = linearBlocks.size();
    std::cout << "[PanmanWriter] Backbone length (numInstances): " << numInstances << std::endl;

    // --- VBlock Audit: Compare getAllBlocks() vs linearBlocks ---
    std::unordered_map<std::string, size_t> seqTotalBpInAllBlocks;
    {
        std::set<VBlockID> linearSet(linearBlocks.begin(), linearBlocks.end());
        std::set<VBlockID> allVBlocksInGraph;
        std::map<VBlockID, std::vector<std::string>> vblockSequences;
        std::map<VBlockID, size_t> vblockLengths;
        std::map<VBlockID, bool> vblockDistant;
        std::unordered_map<std::string, size_t> seqTotalBpInLinear;

        for (auto& blkWeak : rootBlockSet->getAllBlocks()) {
            auto blk = blkWeak.lock();
            if (!blk) continue;

            int maxCopy = blk->getMaxCopy();
            for (int c = 0; c <= maxCopy; ++c) {
                bool hasSegments = false;
                std::vector<std::string> seqsWithCopy;
                for (const auto& seqPair : blk->getSequences()) {
                    for (const auto& segPair : seqPair.second.getSegments()) {
                        if (segPair.second.getCopyCount() == c) {
                            hasSegments = true;
                            seqsWithCopy.push_back(seqPair.first);
                            break;
                        }
                    }
                }
                if (hasSegments) {
                    VBlockID vid = {blk->getId(), c};
                    allVBlocksInGraph.insert(vid);
                    vblockSequences[vid] = seqsWithCopy;
                    vblockLengths[vid] = blk->getConsensusString().size();
                    vblockDistant[vid] = blk->isDistant(c);
                }
            }

            for (const auto& seqPair : blk->getSequences()) {
                for (const auto& segPair : seqPair.second.getSegments()) {
                    seqTotalBpInAllBlocks[seqPair.first] += std::abs(segPair.second.getEnd() - segPair.second.getStart());
                }
            }
        }

        for (size_t i = 0; i < numInstances; i++) {
            const auto& vblk = linearBlocks[i];
            auto blk = rootBlockSet->getBlock(vblk.first);
            if (!blk) continue;
            for (const auto& seqPair : blk->getSequences()) {
                for (const auto& segPair : seqPair.second.getSegments()) {
                    if (segPair.second.getCopyCount() == vblk.second) {
                        seqTotalBpInLinear[seqPair.first] += std::abs(segPair.second.getEnd() - segPair.second.getStart());
                    }
                }
            }
        }

        std::vector<VBlockID> missingVBlocks;
        size_t totalMissingBp = 0;
        for (const auto& vid : allVBlocksInGraph) {
            if (!linearSet.count(vid)) {
                missingVBlocks.push_back(vid);
                totalMissingBp += vblockLengths[vid];
            }
        }

        std::cout << "\n=======================================================\n";
        std::cout << "[PanmanWriter VBlock Audit]\n";
        std::cout << "  Total VBlocks in getAllBlocks(): " << allVBlocksInGraph.size() << "\n";
        std::cout << "  Total VBlocks in LinearGraph:    " << linearBlocks.size() << " (unique: " << linearSet.size() << ")\n";
        if (missingVBlocks.empty()) {
            std::cout << "  Status: ALL VBlocks are INCLUDED in LinearGraph! (0 missing)\n";
        } else {
            std::cout << "  Status: MISSING " << missingVBlocks.size() << " VBlocks! (Total missing consensus bp: " << totalMissingBp << ")\n";
            std::cout << "  -----------------------------------------------------\n";
            std::cout << "  Missing VBlocks List:\n";
            for (const auto& vid : missingVBlocks) {
                std::cout << "    - VBlock (" << vid.first << ", copy " << vid.second << ")"
                          << " | Length: " << vblockLengths[vid] << " bp"
                          << " | isDistant: " << (vblockDistant[vid] ? "YES" : "NO")
                          << " | SeqCount: " << vblockSequences[vid].size()
                          << " | Seqs: [";
                for (size_t s = 0; s < vblockSequences[vid].size(); ++s) {
                    std::cout << vblockSequences[vid][s] << (s + 1 < vblockSequences[vid].size() ? ", " : "");
                }
                std::cout << "]\n";
            }
            std::cout << "  -----------------------------------------------------\n";
        }

        std::cout << "  Per-Sequence Length Audit (AllBlocks vs LinearGraph):\n";
        auto seqNames = rootBlockSet->getSequences();
        for (const auto& seq : seqNames) {
            size_t total = seqTotalBpInAllBlocks[seq];
            size_t lin = seqTotalBpInLinear[seq];
            long long diff = static_cast<long long>(lin) - static_cast<long long>(total);
            std::cout << "    - " << seq << ": AllBlocks = " << total << " bp, Linear = " << lin << " bp (Diff: " << diff << " bp)\n";
        }
        std::cout << "=======================================================\n\n";
    }

    // Collect pre-order & post-order node lists and map to flat integer indices
    std::vector<phylogeny::Node*> preorderNodes;
    collectNodesPreorder(tree->root.get(), preorderNodes);
    size_t numNodes = preorderNodes.size();

    std::unordered_map<uintptr_t, int> nodeToIdxMap;
    std::unordered_map<std::string, int> seqNameToLeafIdx;

    for (size_t idx = 0; idx < numNodes; idx++) {
        nodeToIdxMap[reinterpret_cast<uintptr_t>(preorderNodes[idx])] = static_cast<int>(idx);
        if (preorderNodes[idx]->is_leaf()) {
            seqNameToLeafIdx[preorderNodes[idx]->identifier] = static_cast<int>(idx);
        }
    }

    std::vector<int> parentIdx(numNodes, -1);
    std::vector<std::vector<int>> childrenIdx(numNodes);

    for (size_t idx = 0; idx < numNodes; idx++) {
        phylogeny::Node* pNode = preorderNodes[idx]->parent;
        if (pNode) {
            parentIdx[idx] = nodeToIdxMap[reinterpret_cast<uintptr_t>(pNode)];
        }
        for (const auto& child : preorderNodes[idx]->children) {
            childrenIdx[idx].push_back(nodeToIdxMap[reinterpret_cast<uintptr_t>(child.get())]);
        }
    }

    std::vector<phylogeny::Node*> postorderNodes;
    collectNodesPostorder(tree->root.get(), postorderNodes);
    std::vector<int> postorderIndices(numNodes);
    for (size_t i = 0; i < numNodes; i++) {
        postorderIndices[i] = nodeToIdxMap[reinterpret_cast<uintptr_t>(postorderNodes[i])];
    }

    // 2. Pre-parse Consensus & Map Sequences to Backbone Columns
    std::vector<MafGappedConsensus> blockParsed(numInstances);
    std::vector<std::unordered_map<std::string, std::string>> columnRows(numInstances);
    std::unordered_map<std::string, std::vector<int>> alignedSequences;
    std::unordered_map<std::string, std::vector<int>> alignedStrandSequences;

    auto seqNames = rootBlockSet->getSequences();
    for (const auto& seq : seqNames) {
        alignedSequences[seq].assign(numInstances, -1);
        alignedStrandSequences[seq].assign(numInstances, -1);
    }

    std::unordered_map<std::string, size_t> seqWrittenBp;
    std::unordered_map<std::string, int> seqOverwrittenSegments;
    std::unordered_map<std::string, std::set<std::pair<BlockID, int>>> writtenSegments;

    for (size_t i = 0; i < numInstances; i++) {
        const auto& vblk = linearBlocks[i];
        auto blk = rootBlockSet->getBlock(vblk.first);
        if (!blk) continue;

        std::string consStr = blk->getConsensusString();
        blockParsed[i] = parseMafGappedConsensus(consStr);

        auto& seqs = blk->getSequences();
        for (const auto& seqName : seqNames) {
            auto it = seqs.find(seqName);
            if (it != seqs.end()) {
                const auto& segmentsMap = it->second.getSegments();
                int matchCount = 0;
                for (const auto& segPair : segmentsMap) {
                    if (segPair.second.getCopyCount() == vblk.second) {
                        matchCount++;
                        if (matchCount > 1) {
                            seqOverwrittenSegments[seqName]++;
                            std::cout << "    💥 [Overwrite Detail] Seq: " << seqName 
                                      << " in Block_" << blk->getId() << " (vblk.copy=" << vblk.second << ")"
                                      << " | Seg start=" << segPair.first << ", end=" << segPair.second.getEnd()
                                      << ", len=" << std::abs(segPair.second.getEnd() - segPair.second.getStart()) << " bp\n";
                        }
                        writtenSegments[seqName].insert({blk->getId(), segPair.first});

                        bool segIsRev = segPair.second.isReverse();
                        alignedSequences[seqName][i] = static_cast<int>(i);
                        alignedStrandSequences[seqName][i] = segIsRev ? 0 : 1;

                        std::string gappedRow = consStr;
                        for (const auto& var : segPair.second.getVariants()) {
                            if (var.getType() == VariantType::SNV && var.getStart() < (int)gappedRow.size()) {
                                gappedRow[var.getStart()] = var.getAlt();
                            } else if (var.getType() == VariantType::GAP) {
                                for (int p = var.getStart(); p < var.getEnd() && p < (int)gappedRow.size(); p++) {
                                    gappedRow[p] = '-';
                                }
                            }
                        }
                        columnRows[i][seqName] = gappedRow;
                    }
                }
            }
        }
    }

    // Calculate actual non-gap base count written to PanMAN for each sequence
    for (size_t i = 0; i < numInstances; i++) {
        for (const auto& kv : columnRows[i]) {
            size_t nonGap = 0;
            for (char c : kv.second) {
                if (c != '-' && c != 'x' && c != 'X') nonGap++;
            }
            seqWrittenBp[kv.first] += nonGap;
        }
    }

    std::cout << "\n[PanmanWriter Segment & Base Audit]\n";
    for (const auto& seq : seqNames) {
        size_t totalBp = seqTotalBpInAllBlocks[seq];
        size_t writtenBp = seqWrittenBp[seq];
        long long diff = static_cast<long long>(writtenBp) - static_cast<long long>(totalBp);
        int over = seqOverwrittenSegments[seq];
        std::cout << "  - " << seq << ": AllBlocks = " << totalBp << " bp | Written = " << writtenBp 
                  << " bp (Diff: " << diff << " bp, Overwritten Segments: " << over << ")\n";
    }

    // Check for any unwritten segments
    for (auto& blkWeak : rootBlockSet->getAllBlocks()) {
        auto blk = blkWeak.lock();
        if (!blk) continue;
        for (const auto& seqPair : blk->getSequences()) {
            const auto& seqName = seqPair.first;
            for (const auto& segPair : seqPair.second.getSegments()) {
                if (!writtenSegments[seqName].count({blk->getId(), segPair.first})) {
                    std::cout << "    ⚠️ Unwritten Segment in " << seqName 
                              << ": Block_" << blk->getId() << " (start " << segPair.second.getStart() 
                              << ", end " << segPair.second.getEnd() << ", copy " << segPair.second.getCopyCount() 
                              << ", len " << std::abs(segPair.second.getEnd() - segPair.second.getStart()) << " bp)\n";
                }
            }
        }
    }
    std::cout << "=======================================================\n\n";

    // 2.5 Compute Rotation Indexes and Sequence Inversion for circular/inverted genomes
    std::unordered_map<std::string, int> presentCounter, minRank, maxRank;
    std::unordered_map<std::string, int64_t> minStart, maxStart;
    std::unordered_map<std::string, int> fwdCount, revCount;

    for (size_t r = 0; r < numInstances; r++) {
        const auto& vblk = linearBlocks[r];
        auto blk = rootBlockSet->getBlock(vblk.first);
        if (!blk) continue;

        for (const auto& seqName : seqNames) {
            if (alignedSequences[seqName][r] != -1) {
                int rank = presentCounter[seqName]++;
                bool isRev = (alignedStrandSequences[seqName][r] == 0);
                if (isRev) revCount[seqName]++; else fwdCount[seqName]++;

                int startCoord = -1;
                auto it = blk->getSequences().find(seqName);
                if (it != blk->getSequences().end()) {
                    for (const auto& segPair : it->second.getSegments()) {
                        if (segPair.second.getCopyCount() == vblk.second) {
                            startCoord = std::min(segPair.second.getStart(), segPair.second.getEnd());
                            break;
                        }
                    }
                }

                if (startCoord != -1) {
                    auto mit = minStart.find(seqName);
                    if (mit == minStart.end() || startCoord < mit->second) {
                        minStart[seqName] = startCoord;
                        minRank[seqName] = rank;
                    }
                    auto Mit = maxStart.find(seqName);
                    if (Mit == maxStart.end() || startCoord > Mit->second) {
                        maxStart[seqName] = startCoord;
                        maxRank[seqName] = rank;
                    }
                }
            }
        }
    }

    std::unordered_map<std::string, int> rotationIndexes;
    std::unordered_map<std::string, bool> sequenceInverted;
    for (const auto& seq : seqNames) {
        bool inv = (revCount[seq] > fwdCount[seq]);
        sequenceInverted[seq] = inv;
        int rot = 0;
        if (inv) {
            if (maxRank.count(seq)) rot = maxRank[seq];
        } else {
            if (minRank.count(seq)) rot = minRank[seq];
        }
        rotationIndexes[seq] = rot;
    }

    // 3. Parallel TBB Fitch Parsimony Pass (Block & Nucleotide Mutations)
    std::vector<std::vector<::mga::panman::BlockMutData>> perNodeBlockMuts(numNodes);
    std::vector<std::vector<::mga::panman::NucMutData>> perNodeNucMuts(numNodes);
    std::vector<std::mutex> nodeMutexes(numNodes);

    tbb::parallel_for(size_t(0), numInstances, [&](size_t i) {
        const auto& parsed = blockParsed[i];

        // --- (A) Block Fitch Pass ---
        std::vector<int> bStates(numNodes, 1); // 1 = absent, 2 = forward, 4 = reverse
        for (const auto& seq : seqNames) {
            auto sIt = seqNameToLeafIdx.find(seq);
            if (sIt != seqNameToLeafIdx.end()) {
                int leafIdx = sIt->second;
                if (alignedSequences[seq][i] == -1) {
                    bStates[leafIdx] = 1;
                } else if (alignedStrandSequences[seq][i] == 1) {
                    bStates[leafIdx] = 2;
                } else {
                    bStates[leafIdx] = 4;
                }
            }
        }

        // Post-order Forward Pass for Block
        for (int idx : postorderIndices) {
            if (!childrenIdx[idx].empty()) {
                int orSt = 0, andSt = 0;
                bool first = true;
                for (int cIdx : childrenIdx[idx]) {
                    int cSt = bStates[cIdx];
                    if (first) { andSt = cSt; first = false; }
                    else { andSt &= cSt; }
                    orSt |= cSt;
                }
                bStates[idx] = (andSt != 0) ? andSt : orSt;
            }
        }

        // Pre-order Backward Pass for Block (default state before root is 1)
        for (size_t idx = 0; idx < numNodes; idx++) {
            int pIdx = parentIdx[idx];
            if (pIdx == -1) {
                int currSt = 1;
                while (!(bStates[idx] & currSt) && currSt < 16) { currSt <<= 1; }
                bStates[idx] = currSt;
            } else {
                int pSt = bStates[pIdx];
                if (pSt & bStates[idx]) {
                    bStates[idx] = pSt;
                } else {
                    int currSt = 1;
                    while (!(bStates[idx] & currSt) && currSt < 16) { currSt <<= 1; }
                    bStates[idx] = currSt;
                }
            }
        }

        // Assign Block Mutations (including Root)
        for (size_t idx = 0; idx < numNodes; idx++) {
            int pIdx = parentIdx[idx];
            int pSt = (pIdx == -1) ? 1 : bStates[pIdx];
            int cSt = bStates[idx];
            if (pSt != cSt) {
                ::mga::panman::BlockMutData bmd;
                bmd.primaryBlockId = static_cast<int32_t>(i);
                bmd.secondaryBlockId = -1;
                bmd.chrIdx = 0;

                if (pSt == 1 && cSt != 1) {
                    bmd.blockMutInfo = 1; // BI (Insertion) = 1
                    bmd.inversion = (cSt == 4);
                    std::lock_guard<std::mutex> lock(nodeMutexes[idx]);
                    perNodeBlockMuts[idx].push_back(bmd);
                } else if (pSt != 1 && cSt == 1) {
                    bmd.blockMutInfo = 0; // BD (Deletion) = 0
                    bmd.inversion = false;
                    std::lock_guard<std::mutex> lock(nodeMutexes[idx]);
                    perNodeBlockMuts[idx].push_back(bmd);
                } else if (pSt != 1 && cSt != 1 && pSt != cSt) {
                    bmd.blockMutInfo = 0; // BD with inversion = true is Block Inversion
                    bmd.inversion = true;
                    std::lock_guard<std::mutex> lock(nodeMutexes[idx]);
                    perNodeBlockMuts[idx].push_back(bmd);
                }
            }
        }

        // --- (B) Build 2D Sequence Structure for Leaf Nodes ---
        // Template sequence: parsed.gapFree.size() + 1 positions
        std::vector<std::pair<char, std::vector<char>>> templateSequence(parsed.gapFree.size() + 1, {'-', {}});
        for (size_t j = 0; j < parsed.gapFree.size(); j++) {
            templateSequence[j].first = parsed.gapFree[j];
        }
        for (const auto& gapEntry : parsed.gaps) {
            if (gapEntry.first < templateSequence.size()) {
                templateSequence[gapEntry.first].second.resize(gapEntry.second, '-');
            }
        }

        std::unordered_map<std::string, std::vector<std::pair<char, std::vector<char>>>> seqSequences;
        for (const auto& seq : seqNames) {
            if (alignedSequences[seq][i] == -1 || !columnRows[i].count(seq)) {
                continue;
            }
            auto currSeq = templateSequence;
            const std::string& gRow = columnRows[i][seq];
            for (size_t oldPos = 0; oldPos < gRow.size(); oldPos++) {
                char c = gRow[oldPos];
                int newPos = parsed.oldToNew[oldPos].first;
                int gapPos = parsed.oldToNew[oldPos].second;
                if (gapPos >= 0) {
                    if (newPos < (int)currSeq.size() && gapPos < (int)currSeq[newPos].second.size()) {
                        currSeq[newPos].second[gapPos] = c;
                    }
                } else {
                    if (newPos < (int)currSeq.size()) {
                        currSeq[newPos].first = c;
                    }
                }
            }
            seqSequences[seq] = std::move(currSeq);
        }

        // --- (C) Non-Gap Nucleotide Fitch Pass ---
        for (size_t j = 0; j < parsed.gapFree.size(); j++) {
            int consSt = (1 << PanmanCodec::symbolToCode(parsed.gapFree[j]));
            std::vector<int> nStates(numNodes, 0); // 0 = Block absent
            for (const auto& seq : seqNames) {
                auto sIt = seqNameToLeafIdx.find(seq);
                if (sIt != seqNameToLeafIdx.end()) {
                    int leafIdx = sIt->second;
                    if (alignedSequences[seq][i] != -1 && seqSequences.count(seq)) {
                        char c = seqSequences[seq][j].first;
                        if (c != '-') {
                            nStates[leafIdx] = (1 << PanmanCodec::symbolToCode(c));
                        } else {
                            nStates[leafIdx] = 1; // Gap / deletion state (code 0)
                        }
                    }
                }
            }

            // Post-order Forward Pass
            for (int idx : postorderIndices) {
                if (!childrenIdx[idx].empty()) {
                    int orSt = 0, andSt = 0;
                    bool first = true;
                    for (int cIdx : childrenIdx[idx]) {
                        int cSt = nStates[cIdx];
                        if (cSt != 0) {
                            if (first) { andSt = cSt; first = false; }
                            else { andSt &= cSt; }
                            orSt |= cSt;
                        }
                    }
                    if (!first) {
                        nStates[idx] = (andSt != 0) ? andSt : orSt;
                    } else {
                        nStates[idx] = 0;
                    }
                }
            }

            // Pre-order Backward Pass
            for (size_t idx = 0; idx < numNodes; idx++) {
                int pIdx = parentIdx[idx];
                if (pIdx == -1) {
                    if (nStates[idx] != 0) {
                        if (nStates[idx] & consSt) {
                            nStates[idx] = consSt;
                        } else {
                            int currSt = 1;
                            while (!(nStates[idx] & currSt) && currSt < (1 << 16)) { currSt <<= 1; }
                            nStates[idx] = currSt;
                        }
                    }
                } else {
                    int pSt = nStates[pIdx];
                    if (pSt == 0 || nStates[idx] == 0) {
                        nStates[idx] = 0;
                    } else {
                        if (pSt & nStates[idx]) {
                            nStates[idx] = pSt;
                        } else {
                            int currSt = 1;
                            while (!(nStates[idx] & currSt) && currSt < (1 << 16)) { currSt <<= 1; }
                            nStates[idx] = currSt;
                        }
                    }
                }
            }

            // Assign Non-Gap Mutations
            for (size_t idx = 0; idx < numNodes; idx++) {
                if (nStates[idx] == 0) continue;
                int pIdx = parentIdx[idx];
                int pSt = (pIdx == -1 || nStates[pIdx] == 0) ? consSt : nStates[pIdx];
                int cSt = nStates[idx];
                if (pSt != cSt) {
                    uint8_t mutType = 0;
                    if (pSt == 1 && cSt != 1) {
                        mutType = 2; // Insertion (NI = 2)
                    } else if (pSt != 1 && cSt == 1) {
                        mutType = 1; // Deletion (ND = 1)
                    } else {
                        mutType = 0; // Substitution (NS = 0)
                    }

                    uint8_t code = 0;
                    if (cSt > 1) {
                        int temp = cSt;
                        while (temp > 1) { temp >>= 1; code++; }
                    }

                    ::mga::panman::NucMutData nmd;
                    nmd.primaryBlockId = static_cast<int32_t>(i);
                    nmd.nucPosition = static_cast<int32_t>(j);
                    nmd.nucGapPosition = -1;
                    nmd.nucGapExist = false;
                    nmd.mutInfo = (static_cast<uint32_t>(code) << 8) | (1u << 4) | (mutType & 0x07);

                    std::lock_guard<std::mutex> lock(nodeMutexes[idx]);
                    perNodeNucMuts[idx].push_back(nmd);
                }
            }
        }

        // --- (D) Gap Nucleotide Fitch Pass ---
        for (size_t j = 0; j < templateSequence.size(); j++) {
            for (size_t k = 0; k < templateSequence[j].second.size(); k++) {
                int consSt = 1; // Default consensus is gap (code 0, state 1)
                std::vector<int> nStates(numNodes, 0); // 0 = Block absent
                for (const auto& seq : seqNames) {
                    auto sIt = seqNameToLeafIdx.find(seq);
                    if (sIt != seqNameToLeafIdx.end()) {
                        int leafIdx = sIt->second;
                        if (alignedSequences[seq][i] != -1 && seqSequences.count(seq)) {
                            char c = seqSequences[seq][j].second[k];
                            if (c != '-') {
                                nStates[leafIdx] = (1 << PanmanCodec::symbolToCode(c));
                            } else {
                                nStates[leafIdx] = 1; // Gap state
                            }
                        }
                    }
                }

                // Post-order Forward Pass
                for (int idx : postorderIndices) {
                    if (!childrenIdx[idx].empty()) {
                        int orSt = 0, andSt = 0;
                        bool first = true;
                        for (int cIdx : childrenIdx[idx]) {
                            int cSt = nStates[cIdx];
                            if (cSt != 0) {
                                if (first) { andSt = cSt; first = false; }
                                else { andSt &= cSt; }
                                orSt |= cSt;
                            }
                        }
                        if (!first) {
                            nStates[idx] = (andSt != 0) ? andSt : orSt;
                        } else {
                            nStates[idx] = 0;
                        }
                    }
                }

                // Pre-order Backward Pass
                for (size_t idx = 0; idx < numNodes; idx++) {
                    int pIdx = parentIdx[idx];
                    if (pIdx == -1) {
                        if (nStates[idx] != 0) {
                            if (nStates[idx] & consSt) {
                                nStates[idx] = consSt;
                            } else {
                                int currSt = 1;
                                while (!(nStates[idx] & currSt) && currSt < (1 << 16)) { currSt <<= 1; }
                                nStates[idx] = currSt;
                            }
                        }
                    } else {
                        int pSt = nStates[pIdx];
                        if (pSt == 0 || nStates[idx] == 0) {
                            nStates[idx] = 0;
                        } else {
                            if (pSt & nStates[idx]) {
                                nStates[idx] = pSt;
                            } else {
                                int currSt = 1;
                                while (!(nStates[idx] & currSt) && currSt < (1 << 16)) { currSt <<= 1; }
                                nStates[idx] = currSt;
                            }
                        }
                    }
                }

                // Assign Gap Mutations
                for (size_t idx = 0; idx < numNodes; idx++) {
                    if (nStates[idx] == 0) continue;
                    int pIdx = parentIdx[idx];
                    int pSt = (pIdx == -1 || nStates[pIdx] == 0) ? consSt : nStates[pIdx];
                    int cSt = nStates[idx];
                    if (pSt != cSt) {
                        uint8_t mutType = 0;
                        if (pSt == 1 && cSt != 1) {
                            mutType = 2; // Insertion (NI = 2)
                        } else if (pSt != 1 && cSt == 1) {
                            mutType = 1; // Deletion (ND = 1)
                        } else {
                            mutType = 0; // Substitution (NS = 0)
                        }

                        uint8_t code = 0;
                        if (cSt > 1) {
                            int temp = cSt;
                            while (temp > 1) { temp >>= 1; code++; }
                        }

                        ::mga::panman::NucMutData nmd;
                        nmd.primaryBlockId = static_cast<int32_t>(i);
                        nmd.nucPosition = static_cast<int32_t>(j);
                        nmd.nucGapPosition = static_cast<int32_t>(k);
                        nmd.nucGapExist = true;
                        nmd.mutInfo = (static_cast<uint32_t>(code) << 8) | (1u << 4) | (mutType & 0x07);

                        std::lock_guard<std::mutex> lock(nodeMutexes[idx]);
                        perNodeNucMuts[idx].push_back(nmd);
                    }
                }
            }
        }
    });

    // 4. Serialize Cap'n Proto TreeGroup Message
    capnp::MallocMessageBuilder message;
    ::panman::TreeGroup::Builder treeGroupBuilder = message.initRoot<::panman::TreeGroup>();

    auto treesBuilder = treeGroupBuilder.initTrees(1);
    ::panman::Tree::Builder treeBuilder = treesBuilder[0];
    treeBuilder.setAlphabet(::panman::Alphabet::DNA);

    // Set Newick String
    std::string nwkStr = newickString;
    if (nwkStr.empty() && tree && tree->root) {
        nwkStr = tree->getNewickString();
    }
    if (!nwkStr.empty()) {
        auto newickList = treeBuilder.initNewick(1);
        newickList.set(0, nwkStr);
    } else {
        treeBuilder.initNewick(0);
    }

    // Build Consensus Sequence Map (ConsensusSeqToBlockIds)
    std::map<std::vector<uint32_t>, std::vector<std::pair<int64_t, bool>>> consensusMap;
    std::map<int64_t, int64_t> blockIdToLength;
    for (size_t i = 0; i < numInstances; i++) {
        int64_t blockId = static_cast<int64_t>(i) << 32;
        std::vector<uint32_t> packedCons = PanmanCodec::packConsensusSequence(blockParsed[i].gapFree);
        consensusMap[packedCons].push_back(std::make_pair(blockId, false));
        blockIdToLength[i] = blockParsed[i].gapFree.size();
    }

    auto consensusSeqMapBuilder = treeBuilder.initConsensusSeqMap(consensusMap.size());
    size_t cIdx = 0;
    for (const auto& kv : consensusMap) {
        ::panman::ConsensusSeqToBlockIds::Builder cItem = consensusSeqMapBuilder[cIdx++];
        cItem.setBlockLength(blockIdToLength[kv.second[0].first >> 32]);
        auto blockIdBuilder = cItem.initBlockId(kv.second.size());
        auto blockGapExistBuilder = cItem.initBlockGapExist(kv.second.size());
        auto conSeqBuilder = cItem.initConsensusSeq(kv.first.size());

        for (size_t i = 0; i < kv.second.size(); ++i) {
            blockIdBuilder.set(i, kv.second[i].first);
            blockGapExistBuilder.set(i, kv.second[i].second);
        }
        for (size_t i = 0; i < kv.first.size(); ++i) {
            conSeqBuilder.set(i, kv.first[i]);
        }
    }

    // Build GapList per block
    auto gapsBuilder = treeBuilder.initGaps(numInstances);
    for (size_t i = 0; i < numInstances; i++) {
        ::panman::GapList::Builder gl = gapsBuilder[i];
        gl.setBlockId(static_cast<int64_t>(i) << 32);
        gl.setBlockGapExist(false);

        auto posBuilder = gl.initNucPosition(blockParsed[i].gaps.size());
        auto lenBuilder = gl.initNucGapLength(blockParsed[i].gaps.size());
        for (size_t j = 0; j < blockParsed[i].gaps.size(); j++) {
            posBuilder.set(j, blockParsed[i].gaps[j].first);
            lenBuilder.set(j, blockParsed[i].gaps[j].second);
        }
    }

    // Build Tree Nodes & Mutations
    auto nodesBuilder = treeBuilder.initNodes(numNodes);
    for (size_t nIdx = 0; nIdx < numNodes; ++nIdx) {
        ::panman::Node::Builder nodeBuilder = nodesBuilder[nIdx];

        std::map<int32_t, std::vector<::mga::panman::NucMutData>> bNucMuts;
        std::map<int32_t, ::mga::panman::BlockMutData> bBlockMuts;

        auto consolidatedNucMuts = consolidateNucMutations(perNodeNucMuts[nIdx]);
        for (const auto& nm : consolidatedNucMuts) {
            bNucMuts[nm.primaryBlockId].push_back(nm);
        }
        for (const auto& bm : perNodeBlockMuts[nIdx]) {
            bBlockMuts[bm.primaryBlockId] = bm;
        }

        std::set<int32_t> affectedBlocks;
        for (const auto& kv : bNucMuts) affectedBlocks.insert(kv.first);
        for (const auto& kv : bBlockMuts) affectedBlocks.insert(kv.first);

        auto mutationsBuilder = nodeBuilder.initMutations(affectedBlocks.size());
        size_t mIdx = 0;
        for (int32_t bId : affectedBlocks) {
            ::panman::Mutation::Builder mutBuilder = mutationsBuilder[mIdx++];
            int64_t fullBlockId = static_cast<int64_t>(bId) << 32;
            mutBuilder.setBlockId(fullBlockId);
            mutBuilder.setBlockGapExist(false);

            if (bBlockMuts.count(bId)) {
                mutBuilder.setBlockMutExist(true);
                mutBuilder.setBlockMutInfo(bBlockMuts[bId].blockMutInfo);
                mutBuilder.setBlockInversion(bBlockMuts[bId].inversion);
                mutBuilder.setChrIdx(bBlockMuts[bId].chrIdx);
            } else {
                mutBuilder.setBlockMutExist(false);
                mutBuilder.setBlockMutInfo(2); // Present / no-change
                mutBuilder.setBlockInversion(false);
                mutBuilder.setChrIdx(0);
            }

            if (bNucMuts.count(bId)) {
                const auto& nucList = bNucMuts[bId];
                auto nucMutListBuilder = mutBuilder.initNucMutation(nucList.size());
                for (size_t i = 0; i < nucList.size(); ++i) {
                    nucMutListBuilder[i].setNucPosition(nucList[i].nucPosition);
                    nucMutListBuilder[i].setNucGapPosition(nucList[i].nucGapPosition);
                    nucMutListBuilder[i].setNucGapExist(nucList[i].nucGapExist);
                    nucMutListBuilder[i].setMutInfo(nucList[i].mutInfo);
                }
            } else {
                mutBuilder.initNucMutation(0);
            }
        }

        nodeBuilder.initAnnotations(0);
    }

    // Initialize circular / rotation / inverted structures
    treeBuilder.initCircularSequences(0);
    treeBuilder.initRotationIndexes(0);
    treeBuilder.initSequencesInverted(0);
    treeBuilder.initChrLists(0);

    // 5. Write Compressed Message to File
    std::ofstream outputFile(outputPath, std::ios::binary);
    if (!outputFile.is_open()) {
        std::cerr << "[PanmanWriter Error] Cannot open file for writing: " << outputPath << "\n";
        return false;
    }

    boost::iostreams::filtering_streambuf<boost::iostreams::output> outBuffer;
    boost::iostreams::lzma_params params;
    params.level = 9;
    outBuffer.push(boost::iostreams::lzma_compressor(params));
    outBuffer.push(outputFile);

    std::ostream outStream(&outBuffer);
    kj::std::StdOutputStream kjOutputStream(outStream);

    ::capnp::writeMessage(kjOutputStream, message);

    boost::iostreams::close(outBuffer);
    outputFile.close();

    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = endTime - startTime;
    std::cout << "[PanmanWriter] Successfully wrote PanMAN file in " << diff.count() << " seconds: " << outputPath << std::endl;
    return true;
}

} // namespace panman
} // namespace mga

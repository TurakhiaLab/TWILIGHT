#include "mga.hpp"
#include "kseq.h"
#include "block_manager.hpp"
#include "phylogeny.hpp"

#include <unordered_set>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <climits>
#include <limits>
#include <memory>
#include "message_manager.hpp"

#include <zlib.h>
#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>

namespace fs = boost::filesystem;


KSEQ_INIT2(, gzFile, gzread);

// Gzip compression function
std::string gzip_compress(const std::string &input) {
    z_stream zs{};
    deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY); // 15+16 = gzip header
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(input.data()));
    zs.avail_in = input.size();
    std::string output;
    output.resize(compressBound(input.size())); 
    zs.next_out = reinterpret_cast<Bytef*>(&output[0]);
    zs.avail_out = output.size();
    int ret = deflate(&zs, Z_FINISH);
    if (ret != Z_STREAM_END) {
        mga::msg::debug("Compression error: %d\n", ret);
        deflateEnd(&zs);
        return {};
    }
    output.resize(zs.total_out);
    deflateEnd(&zs);
    return output;
}

// The function signature in mga.hpp must match this
std::unique_ptr<BlockManager> mga::io::readSequences(std::string& fileName, Option& option, phylogeny::Tree& tree) {
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    
    std::string seqFileName = fileName;
    mga::msg::debug("Reading sequences from %s...\n", seqFileName.c_str());
    
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        mga::msg::error("Cannot open sequence file: %s\n", seqFileName.c_str());
        exit(1);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);

    auto blockManager = std::make_unique<BlockManager>();

    int seqNum = 0;
    uint64_t totalLen = 0;
    int maxLen = 0;
    int minLen = INT_MAX;

    std::unordered_set<std::string> loadedNames;

    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        std::string seqName_full = kseq_rd->name.s;
        
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";
        phylogeny::Node* seqNode = nullptr;

        auto it_full = tree.allNodes.find(seqName_full);
        if (it_full != tree.allNodes.end()) {
            seqName = seqName_full;
            seqNode = it_full->second;
        } else {
            auto it_noblank = tree.allNodes.find(seqName_noblank);
            if (it_noblank != tree.allNodes.end()) {
                seqName = seqName_noblank;
                seqNode = it_noblank->second;
            }
        }

        if (seqNode != nullptr) {
            if (loadedNames.count(seqName)) {
                mga::msg::warn("Duplicate leaf name found in sequence file: %s. Only first occurrence will be kept.\n", seqName.c_str());
            } else {
                std::string seqContent = std::string(kseq_rd->seq.s, seqLen);

                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
                if (seqLen == 0) mga::msg::debug("Null sequence encountered: %s\n", seqName.c_str());
                totalLen += seqLen;
                
                BlockSet* blockSet = blockManager->createBlockSet(seqNode->identifier);
                blockSet->addSequenceName(seqName);
                
                std::shared_ptr<Block> newBlock = blockSet->createBlock(seqContent);

                Sequence info(seqName);
                info.addSegment(0, seqLen);

                // Add sequence to block. CIGAR is empty, so no variations will be generated.
                newBlock->addSequence(info);

                loadedNames.insert(seqName);

                blockManager->addSequenceLength(seqName, seqLen);
                blockManager->addSequence(seqName, seqContent);
                seqNum++;
            }
        }
    }
    
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    if (seqNum == 0) {
        mga::msg::error("No sequences were read from the input that match any node in the tree.\n");
    }

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> seqReadDuration = seqReadEnd - seqReadStart;

    uint32_t avgLen = (seqNum > 0) ? totalLen / seqNum : 0;

    mga::msg::info("Read %d sequences (Total: %lu bp, Avg: %u bp, Range: [%d, %d] bp) in %.2f ms\n",
                   seqNum, totalLen, avgLen, (seqNum > 0 ? minLen : 0), maxLen, seqReadDuration.count());

    return blockManager;
}


void mga::io::writeAlignment(std::string fileName, StringPairs& seqs, bool compressed, bool append) {
    auto writeStart = std::chrono::high_resolution_clock::now();

    if (compressed) {
        fileName += ".gz";
        std::vector<std::string> compressed_chunks(seqs.size());
        tbb::parallel_for(size_t(0), seqs.size(), [&](size_t i) {
            std::string fasta_chunk = ">" + seqs[i].first + "\n" + seqs[i].second + "\n";
            compressed_chunks[i] = gzip_compress(fasta_chunk);
            std::string().swap(seqs[i].second);
        });
        std::ofstream outFile;
        if (append) outFile.open(fileName, std::ios::binary | std::ios::app);
        else        outFile.open(fileName, std::ios::binary);
        if (!outFile) {
            mga::msg::error("Cannot open file for writing: %s\n", fileName.c_str());
            exit(1);
        }
        for (auto &c : compressed_chunks) {
            outFile.write(c.data(), c.size());
        }
        outFile.close();
    }
    else {
        std::ofstream outFile;
        if (append) outFile.open(fileName, std::ios::app);
        else        outFile.open(fileName);
        if (!outFile) {
            mga::msg::error("Cannot open file for writing: %s\n", fileName.c_str());
            exit(1);
        }
        for (auto it = seqs.begin(); it != seqs.end(); ++it) {
            outFile << ('>' + it->first + '\n');
            outFile << (it->second + '\n');
        }
        outFile.close();
    }

    auto writeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> writeDuration = writeEnd - writeStart;
    mga::msg::info("Alignment written to %s in %.2f ms\n", fileName.c_str(), writeDuration.count());
    return;
}


std::string generateAlignmentString(const std::string& consensus, const Variants& variations) {
    std::string aliSeq = consensus;
    
    for (auto& var : variations) {
        if (var.getType() == VariantType::SNV) {
            int pos = var.getStart();
            if (pos >= 0 && pos < aliSeq.length()) {
                aliSeq[pos] = var.getAlt();
            }
        } else if (var.getType() == VariantType::GAP) {
            for (int i = var.getStart(); i < var.getEnd(); ++i) {
                if (i >= 0 && i < aliSeq.length()) {
                    aliSeq[i] = '-';
                }
            }
        }
    }
    
    return aliSeq;
}

// Helper function to write a single MAF block (either homology-merged or single-copy)
static void writeSingleMAFBlock(std::ofstream& mafFile,
                                const std::shared_ptr<Block>& blk,
                                int target_copy,
                                const std::string& blockLabel,
                                const std::unordered_map<std::string, int>& seqLengths) {
    if (!blk) return;

    const std::string& consensus = blk->getConsensus().getConsensusString();
    int consLen = consensus.length();

    int segmentCount = 0;
    int totalVarLen = 0;

    // Recalculate score considering only eligible segments
    for (auto& seqEntry : blk->getSequences()) {
        for (auto& segPair : seqEntry.second.getSegments()) {
            Segment& seg = segPair.second;
            if (target_copy != -1 && seg.getCopyCount() != target_copy) continue;

            segmentCount++;
            for (auto& var : seg.getVariants()) {
                totalVarLen += (var.getEnd() - var.getStart());
            }
        }
    }

    if (segmentCount == 0) return;

    double score = static_cast<double>(consLen * segmentCount - totalVarLen);

    // 1. Output 'a' line
    mafFile << "a score=" << std::fixed << std::setprecision(1) << score 
            << " blockID=" << blockLabel << "\n";

    // 2. Output 's' lines for each sequence segment
    for (auto& seqEntry : blk->getSequences()) {
        std::string seqID = seqEntry.first;
        long srcSize = 0;
        auto itLen = seqLengths.find(seqID);
        if (itLen != seqLengths.end()) {
            srcSize = itLen->second;
        }

        for (auto& segPair : seqEntry.second.getSegments()) {
            Segment& seg = segPair.second;
            if (target_copy != -1 && seg.getCopyCount() != target_copy) continue;

            std::string aliString = generateAlignmentString(consensus, seg.getVariants());

            // Count ungapped base size to strictly adhere to UCSC MAF specification
            int size = 0;
            for (char c : aliString) {
                if (c != '-' && c != '.') size++;
            }

            int fwdStart = std::min(seg.getStart(), seg.getEnd());
            int fwdEnd   = std::max(seg.getStart(), seg.getEnd());
            char strand  = seg.isReverse() ? '-' : '+';

            // In UCSC MAF, minus strand start is measured from the 5' end of the reverse complement:
            // mafStart = srcSize - fwdEnd
            int mafStart = fwdStart;
            if (strand == '-') {
                mafStart = (srcSize > 0) ? (srcSize - fwdEnd) : fwdStart;
            }

            mafFile << "s " << std::left << std::setw(20) << seqID << " "
                    << std::right << std::setw(10) << mafStart << " "
                    << std::right << std::setw(8) << size << " "
                    << strand << " "
                    << std::right << std::setw(10) << srcSize << " "
                    << aliString << "\n";
        }
    }
    mafFile << "\n";
}

void mga::io::writeMAF(BlockSet* blockSet, const std::string& outputFileName, const std::unordered_map<std::string, int>& seqLengths) {
    if (!blockSet) return;
    auto start = std::chrono::high_resolution_clock::now();
    mga::msg::debug("Writing MAF to %s...\n", outputFileName.c_str());

    std::ofstream mafFile(outputFileName);
    if (!mafFile.is_open()) {
        mga::msg::error("Could not open file %s for writing MAF.\n", outputFileName.c_str());
        return;
    }

    // 1. MAF Header
    mafFile << "##maf version=1\n";
    mafFile << "# Generated from BlockSet ID: " << blockSet->getId() << "\n";
    mafFile << "# Mode: Homology-Merged\n\n";

    // 2. Output physical blocks in the order of appearance on linear backbone first, then unvisited blocks
    std::unordered_set<BlockID> seen_blocks;
    auto vblocks = blockSet->getLinearizeBlocks();
    for (const auto& vid : vblocks) {
        if (seen_blocks.insert(vid.first).second) {
            if (auto blk = blockSet->getBlock(vid.first)) {
                std::string blockLabel = std::to_string(blk->getId());
                writeSingleMAFBlock(mafFile, blk, -1, blockLabel, seqLengths);
            }
        }
    }

    auto all_blocks = blockSet->getAllBlocks();
    for (auto& weak_blk : all_blocks) {
        if (auto blk = weak_blk.lock()) {
            if (seen_blocks.insert(blk->getId()).second) {
                std::string blockLabel = std::to_string(blk->getId());
                writeSingleMAFBlock(mafFile, blk, -1, blockLabel, seqLengths);
            }
        }
    }

    mafFile.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    mga::msg::info("MAF written to %s in %.2f ms\n", outputFileName.c_str(), duration.count());
}

void mga::io::writeLinearMAF(BlockSet* blockSet, const std::string& outputFileName, const std::unordered_map<std::string, int>& seqLengths) {
    if (!blockSet) return;
    auto start = std::chrono::high_resolution_clock::now();
    mga::msg::debug("Writing linear MAF to %s...\n", outputFileName.c_str());

    std::ofstream mafFile(outputFileName);
    if (!mafFile.is_open()) {
        mga::msg::error("Could not open file %s for writing Linear MAF.\n", outputFileName.c_str());
        return;
    }

    // 1. MAF Header
    mafFile << "##maf version=1\n";
    mafFile << "# Generated from BlockSet ID: " << blockSet->getId() << "\n";
    mafFile << "# Mode: Linearized\n\n";

    // 2. Traverse linear backbone and output blocks with blockID=<BlockID>.<copy num>
    auto vblocks = blockSet->getLinearizeBlocks();
    for (const auto& vid : vblocks) {
        if (auto blk = blockSet->getBlock(vid.first)) {
            std::string blockLabel = std::to_string(vid.first) + "." + std::to_string(vid.second);
            writeSingleMAFBlock(mafFile, blk, vid.second, blockLabel, seqLengths);
        }
    }

    mafFile.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    mga::msg::info("Linear MAF written to %s in %.2f ms\n", outputFileName.c_str(), duration.count());
}

void mga::io::writeConcatMSA(BlockSet* blockSet, const std::string& outputFileName, double threshold, int minBlockLen) {
    auto start = std::chrono::high_resolution_clock::now();
    mga::msg::debug("Generating concatenated MSA to %s...\n", outputFileName.c_str());

    if (!blockSet) {
        mga::msg::error("Null BlockSet provided to writeConcatMSA.\n");
        return;
    }

    // 1. Gather all unique sequence names in deterministic order
    std::vector<std::string> allSeqNames(blockSet->getSequences().begin(), blockSet->getSequences().end());
    if (allSeqNames.empty()) {
        mga::msg::warn("No sequences in BlockSet for writeConcatMSA.\n");
        return;
    }
    std::sort(allSeqNames.begin(), allSeqNames.end());
    size_t totalSeqs = allSeqNames.size();

    // Normalize threshold (accept e.g. 95 or 0.95)
    double normThreshold = (threshold > 1.0) ? (threshold / 100.0) : threshold;

    // 2. Collect candidate blocks in linear backbone order (Option A)
    std::unordered_set<BlockID> seen_blocks;
    std::vector<std::shared_ptr<Block>> candidate_blocks;

    auto vblocks = blockSet->getLinearizeBlocks();
    for (const auto& vid : vblocks) {
        if (seen_blocks.insert(vid.first).second) {
            if (auto blk = blockSet->getBlock(vid.first)) {
                candidate_blocks.push_back(blk);
            }
        }
    }

    // Append any remaining distant / unvisited blocks not on the linear backbone
    auto all_blocks = blockSet->getAllBlocks();
    for (auto& weak_blk : all_blocks) {
        if (auto blk = weak_blk.lock()) {
            if (seen_blocks.insert(blk->getId()).second) {
                candidate_blocks.push_back(blk);
            }
        }
    }

    // 3. Filter blocks exceeding occupancy threshold and minimum length
    std::vector<std::shared_ptr<Block>> eligible_blocks;
    for (auto& blk : candidate_blocks) {
        if (!blk) continue;
        if (static_cast<int>(blk->getConsensus().getConsensusString().length()) < minBlockLen) {
            continue;
        }
        size_t presentSeqCount = 0;
        const auto& seqMap = blk->getSequences();
        for (const auto& seqName : allSeqNames) {
            auto it = seqMap.find(seqName);
            if (it != seqMap.end() && !it->second.getSegments().empty()) {
                presentSeqCount++;
            }
        }
        double occupancy = static_cast<double>(presentSeqCount) / static_cast<double>(totalSeqs);
        if (occupancy >= normThreshold) {
            eligible_blocks.push_back(blk);
        }
    }

    mga::msg::debug("writeConcatMSA: %zu / %zu blocks passed occupancy (>= %.2f%%) and length (>= %d bp) filters. Generating concatenated MSA...\n",
                    eligible_blocks.size(), candidate_blocks.size(), normThreshold * 100.0, minBlockLen);

    // 4. Calculate total MSA length and pre-allocate string buffers
    size_t totalMsaLen = 0;
    for (const auto& blk : eligible_blocks) {
        totalMsaLen += blk->getConsensus().getConsensusString().length();
    }

    std::unordered_map<std::string, std::string> msaMap;
    for (const auto& seqName : allSeqNames) {
        msaMap[seqName].reserve(totalMsaLen);
    }

    // 5. Build MSA strings block by block
    for (const auto& blk : eligible_blocks) {
        const std::string& consensus = blk->getConsensus().getConsensusString();
        size_t consLen = consensus.length();
        const auto& seqMap = blk->getSequences();

        for (const auto& seqName : allSeqNames) {
            auto seqIt = seqMap.find(seqName);
            const Segment* bestSeg = nullptr;

            if (seqIt != seqMap.end() && !seqIt->second.getSegments().empty()) {
                int bestTotalVar = std::numeric_limits<int>::max();
                int bestGapLen = std::numeric_limits<int>::max();
                int bestCopyCount = std::numeric_limits<int>::max();

                for (const auto& segPair : seqIt->second.getSegments()) {
                    const Segment& seg = segPair.second;
                    int snpCount = 0;
                    int gapLen = 0;
                    for (const auto& var : seg.getVariants()) {
                        if (var.getType() == VariantType::SNV) {
                            snpCount++;
                        } else if (var.getType() == VariantType::GAP) {
                            gapLen += (var.getEnd() - var.getStart());
                        }
                    }
                    int totalVar = snpCount + gapLen;
                    int copyCount = seg.getCopyCount();

                    bool isBetter = false;
                    if (bestSeg == nullptr) {
                        isBetter = true;
                    } else if (totalVar < bestTotalVar) {
                        isBetter = true;
                    } else if (totalVar == bestTotalVar) {
                        // Tie-breaker 1: choose fewer gaps (gaps are uninformative, SNPs carry phylogenetic signal)
                        if (gapLen < bestGapLen) {
                            isBetter = true;
                        } else if (gapLen == bestGapLen) {
                            // Tie-breaker 2: choose smaller copy number
                            if (copyCount < bestCopyCount) {
                                isBetter = true;
                            } else if (copyCount == bestCopyCount) {
                                // Tie-breaker 3: deterministic start coordinate order
                                if (seg.getStart() < bestSeg->getStart()) {
                                    isBetter = true;
                                }
                            }
                        }
                    }

                    if (isBetter) {
                        bestSeg = &seg;
                        bestTotalVar = totalVar;
                        bestGapLen = gapLen;
                        bestCopyCount = copyCount;
                    }
                }
            }

            if (bestSeg != nullptr) {
                std::string aliStr = generateAlignmentString(consensus, bestSeg->getVariants());
                if (aliStr.length() < consLen) {
                    aliStr.append(consLen - aliStr.length(), '-');
                }
                msaMap[seqName] += aliStr;
            } else {
                // Sequence absent in this block: fill with gaps
                msaMap[seqName].append(consLen, '-');
            }
        }
    }

    // 6. Write FASTA format output
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        mga::msg::error("Could not open file %s for writing concat MSA.\n", outputFileName.c_str());
        return;
    }

    for (const auto& seqName : allSeqNames) {
        outFile << ">" << seqName << "\n";
        outFile << msaMap[seqName] << "\n";
    }

    outFile.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    mga::msg::info("Concatenated MSA written to %s (Length: %zu bp) in %.2f ms (%.3f s)\n",
                   outputFileName.c_str(), totalMsaLen, duration.count(), duration.count() / 1000.0);
}
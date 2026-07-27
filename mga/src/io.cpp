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
#include <memory>

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
        std::cerr << "Compression error: " << ret << "\n";
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
    
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
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
                printf("WARNING: duplicate leaf names found in the sequence file! Leaf name: %s. Only the first occurrence will be kept.\n", seqName.c_str());
            } else {
                std::string seqContent = std::string(kseq_rd->seq.s, seqLen);

                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
                if (seqLen == 0) std::cerr << "Null sequences, " << seqName << '\n';
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
        std::cerr << "Error: no sequences were read from the input that are also in the tree.\n";
    }

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;

    uint32_t avgLen = (seqNum > 0) ? totalLen / seqNum : 0;

    // blockManager->updateLongestSequences();

    std::cerr << "===== Sequence Summary =====\n";
    std::cerr << "Number of sequences read and found in tree: " << seqNum << '\n';
    std::cerr << "Max. Length: " << maxLen << '\n';
    std::cerr << "Min. Length: " << minLen << '\n';
    std::cerr << "Avg. Length: " << avgLen << '\n';
    std::cerr << "Sequences read in " <<  seqReadTime.count() / 1000000 << " ms\n";

    return blockManager;
}


void mga::io::writeAlignment(std::string fileName, StringPairs& seqs, bool compressed, bool append) {
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
            fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
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
            fprintf(stderr, "ERROR: cant open file: %s\n", fileName.c_str());
            exit(1);
        }
        for (auto it = seqs.begin(); it != seqs.end(); ++it) {
            outFile << ('>' + it->first + '\n');
            outFile << (it->second + '\n');
        }
        outFile.close();
    }
    return;
}


std::string generateAlignmentString(const std::string& consensus, Variants& variations) {
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

void mga::io::writeMAF(BlockSet* blockSet, const std::string& outputFileName, bool grouped) {
    std::ofstream mafFile(outputFileName);
    if (!mafFile.is_open()) {
        std::cerr << "Error: Could not open file " << outputFileName << " for writing MAF.\n";
        return;
    }

    // 1. MAF Header
    mafFile << "##maf version=1\n";
    mafFile << "# Generated from BlockSet ID: " << blockSet->getId() << "\n";
    mafFile << "# Mode: " << (grouped ? "Homology-Merged" : "Linearized") << "\n\n";
    
    // 🌟 定義一個小結構來統一後續的處理流程
    struct IterItem {
        std::shared_ptr<Block> blk;
        int target_copy; // -1 代表 Grouped 模式 (全印)
    };
    std::vector<IterItem> items_to_write;

    // 🌟 2. 根據模式準備要輸出的 Block 與 Copy 目標
    if (grouped) {
        // [Grouped 模式]：直接走訪所有實體 Block，target_copy 設為 -1
        auto all_blocks = blockSet->getAllBlocks();
        for (auto& weak_blk : all_blocks) {
            if (auto blk = weak_blk.lock()) {
                items_to_write.push_back({blk, -1});
            }
        }
    } else {
        // [Linearized 模式]：走訪 1D 骨幹，根據 VBlockID 決定要印哪個 Copy
        // 假設你的 getLinearizeBlocks() 現在回傳的是 std::vector<VBlockID> (即 pair<BlockID, int>)
        auto vblocks = blockSet->getLinearizeBlocks(); 
        for (auto& vid : vblocks) {
            if (auto blk = blockSet->getBlock(vid.first)) {
                items_to_write.push_back({blk, vid.second});
            }
        }
    }

    // 3. 執行統一的輸出邏輯
    for (auto& item : items_to_write) {
        auto blk = item.blk;
        int target_copy = item.target_copy;

        BlockID currentBlockId = blk->getId();
        const std::string& consensus = blk->getConsensus().getConsensusString();
        int consLen = consensus.length();
        
        int segmentCount = 0;
        int totalVarLen = 0;
        
        // 🌟 重新計算 Score，只算符合條件的 Segments
        for (auto& seqEntry : blk->getSequences()) {
            for (auto& segPair : seqEntry.second.getSegments()) {
                Segment& seg = segPair.second;
                
                // [過濾機制] 如果不是 grouped 模式，且 Copy 號碼不對，就跳過！
                if (target_copy != -1 && seg.getCopyCount() != target_copy) continue;

                segmentCount++;
                for (auto& var : seg.getVariants()) {
                    totalVarLen += (var.getEnd() - var.getStart());
                }
            }
        }
        
        if (segmentCount == 0) continue;

        double score = static_cast<double>(consLen * segmentCount - totalVarLen);

        // 🌟 決定 a 行的 blockID 標籤
        std::string blockLabel = std::to_string(currentBlockId);
        if (target_copy != -1) {
            blockLabel += "_C" + std::to_string(target_copy); // 產生像 "5_C0", "5_C1" 這樣的標籤
        }

        // 3. 'a' line
        mafFile << "a score=" << std::fixed << std::setprecision(1) << score 
                << " blockID=" << blockLabel << "\n";

        // 4. 輸出各個 Segment 的 's' 行
        for (auto& seqEntry : blk->getSequences()) {
            std::string seqID = seqEntry.first;
            
            for (auto& segPair : seqEntry.second.getSegments()) {
                Segment& seg = segPair.second;
                
                // 🌟 同樣的過濾機制
                if (target_copy != -1 && seg.getCopyCount() != target_copy) continue;
                
                std::string aliString = generateAlignmentString(consensus, seg.getVariants());
                
                int start = seg.getStart();
                int size = std::abs(seg.getEnd() - seg.getStart()); 
                char strand = seg.isReverse() ? '-' : '+';
                long srcSize = 0; 

                mafFile << "s " << std::left << std::setw(20) << seqID << " "
                        << std::right << std::setw(10) << start << " "
                        << std::right << std::setw(8) << size << " "
                        << strand << " "
                        << std::right << std::setw(10) << srcSize << " "
                        << aliString << "\n";
            }
        }
        mafFile << "\n"; 
    }

    mafFile.close();
}
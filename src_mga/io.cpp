#include "mga.hpp"
#include "kseq.h"
#include "block.hpp"
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

void mga::io::writeMAF(BlockSet* blockSet, const std::string& outputFileName) {
    std::ofstream mafFile(outputFileName);
    if (!mafFile.is_open()) {
        std::cerr << "Error: Could not open file " << outputFileName << " for writing MAF.\n";
        return;
    }

    // 1. MAF Header
    mafFile << "##maf version=1\n";
    mafFile << "# Generated from BlockSet ID: " << blockSet->getId() << "\n\n";
    BlockWeakPtrs blocks_;
    
    auto blocks_ids = blockSet->getLinearizeBlocks();

    for (auto& id : blocks_ids) blocks_.push_back(blockSet->getBlock(id));
    
    // 🚨 新增：用來追蹤每個 BlockID 寫入的次數 (Occurrence tracker)
    std::unordered_map<BlockID, int> blockCounter; 

    // 2. Blocks
    for (auto& block_ptr : blocks_) {
        auto blk = block_ptr.lock();
        if (!blk) continue; // 防呆：如果弱指標失效則跳過

        BlockID currentBlockId = blk->getId();
        const std::string& consensus = blk->getConsensus();
        int consLen = consensus.length();
        
        int segmentCount = 0;
        int totalVarLen = 0;
        
        for (auto& seqEntry : blk->getSequences()) {
            for (auto& segPair : seqEntry.second.getSegments()) {
                segmentCount++;
                for (auto& var : segPair.second.getVariants()) {
                    totalVarLen += (var.getEnd() - var.getStart());
                }
            }
        }
        
        if (segmentCount == 0) continue;

        double score = static_cast<double>(consLen * segmentCount - totalVarLen);

        // 🚨 新增：取得並增加這個 BlockID 的出現次數
        int currentCount = blockCounter[currentBlockId]++;
        
        // 組合成你要的 extension 格式，例如 "3_0", "3_1"
        std::string customBlockID = std::to_string(currentBlockId) + "_" + std::to_string(currentCount);

        // 3. 'a' line (加上自訂的 blockID)
        mafFile << "a score=" << std::fixed << std::setprecision(1) << score 
                << " blockID=" << customBlockID << "\n";

        // 4. 輸出各個 Segment 的 's' 行
        for (auto& seqEntry : blk->getSequences()) {
            std::string seqID = seqEntry.first;
            
            for (auto& segPair : seqEntry.second.getSegments()) {
                Segment& seg = segPair.second;
                
                // 產生包含 SNV 與 Gap 的比對序列字串
                std::string aliString = generateAlignmentString(consensus, seg.getVariants());
                
                // 準備 s 行參數
                int start = seg.getStart();
                int size = std::abs(seg.getEnd() - seg.getStart()); // 實際佔用的鹼基數量
                char strand = seg.isReverse() ? '-' : '+';
                
                // 由於目前的資料結構未直接存放整條 src 染色體的總長度，這裡先填 0 (這不影響純序列分析，但若是上傳 genome browser 需後處理修正)
                long srcSize = 0; 

                // 格式化輸出以保持整齊
                mafFile << "s " << std::left << std::setw(20) << seqID << " "
                        << std::right << std::setw(10) << start << " "
                        << std::right << std::setw(8) << size << " "
                        << strand << " "
                        << std::right << std::setw(10) << srcSize << " "
                        << aliString << "\n";
            }
        }
        mafFile << "\n"; // 每個 Block 結束後空一行
    }

    mafFile.close();
}
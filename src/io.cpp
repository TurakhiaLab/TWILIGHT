#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include "score.hpp"
#include "kseq.h"
#include <zlib.h>
#include <filesystem>
#include <tbb/parallel_for.h>
#include <tbb/spin_rw_mutex.h>
#include <chrono>

namespace fs = std::filesystem;

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
        msg::error("Compression error: %d\n", ret);
        deflateEnd(&zs);
        return {};
    }
    deflateEnd(&zs);
    // Return a string with exact size and capacity to avoid retaining oversized buffers
    return std::string(output.data(), zs.total_out);
}


void msa::io::readSequenceNames(std::string seqFile, std::unordered_set<std::string>& seqNames)
{
    gzFile f_rd = gzopen(seqFile.c_str(), "r");
    if (!f_rd) {
        throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Failed to open file: " + seqFile);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    while (kseq_read(kseq_rd) >= 0) {
        std::string seqName_full = kseq_rd->name.s;
        seqNames.insert(seqName_full);
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);
    return;
}

void msa::io::readSequences(std::string fileName, SequenceDB* database, Option* option, Tree*& tree, int subtree) {
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = fileName;
    bool placed = (option->alnMode == 3 && seqFileName == option->seqFile);
    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Can't open file: " + seqFileName);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    
    int seqNum_init = database->sequences.size(), maxLen = 0, minLen = INT_MAX;
    int seqNum = seqNum_init;
    uint64_t totalLen = 0; // For average length calculation
    std::vector<int> seqsLens; // For median length calculation
    
    // Read sequences
    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        std::string seqName_full = kseq_rd->name.s;
        std::string seqName_noblank = seqName_full.substr(0, seqName_full.find(' '));
        std::string seqName = "";
        if (tree->allNodes.find(seqName_full) != tree->allNodes.end()) seqName = seqName_full;
        else if (tree->allNodes.find(seqName_noblank) != tree->allNodes.end()) seqName = seqName_noblank;
        if (seqName != "") {
            if (database->name_map.find(seqName) != database->name_map.end()) {
                msg::warn("Duplicate leaf names found in the sequence file! Leaf name: %s. Only the first occurrence will be kept.\n", seqName.c_str());
            }
            else {
                int subtreeIdx = tree->allNodes[seqName]->grpID;
                std::string seq = std::string(kseq_rd->seq.s, seqLen);
                // seqs[seqName] = std::make_pair(seq, subtreeIdx);
                if (seqLen > maxLen) maxLen = seqLen;
                if (seqLen < minLen) minLen = seqLen;
                if (seqLen == 0) msg::warn("Null sequence: %s\n", seqName.c_str());
                totalLen += seqLen;
                database->addSequence(seqNum, seqName, seq, subtreeIdx, tree->allNodes[seqName]->weight, option->debug, option->alnMode);
                if (option->alnMode == PLACE_WO_TREE) database->subtreeAln[database->name_map[seqName]->id] = alnPath (seqLen, 0);
                tree->allNodes[seqName]->placed = placed; 
                ++seqNum;
                seqsLens.push_back(seqLen);
            }
        }
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);

    // Prune Tree if necessary
    if (tree->m_numLeaves != seqNum && (option->alnMode == DEFAULT_ALN || option->alnMode == ACCURATE_ALN)) {
        msg::warn("Mismatch between the number of leaves and the number of sequences, (%lu != %d)\n", tree->m_numLeaves, seqNum);
        int kk = 0;
        for (auto node: tree->allNodes) {
            if (node.second->is_leaf()) {
                if (database->name_map.find(node.second->identifier) == database->name_map.end()) {
                    msg::info("Missing %s\n", node.second->identifier.c_str());
                    ++kk;
                }
            }
        }
        msg::info("Prune the tree according to the existing sequences.\n");
        std::unordered_set<std::string> seqNames;
        for (auto it = database->name_map.begin(); it != database->name_map.end(); ++it) seqNames.insert(it->first);
        phylogeny::pruneTree(tree, seqNames);
    }

    // Identify low quality sequences

    if (seqNum == seqNum_init) {
        throw msa::TwilightError(msa::ErrorCode::EMPTY_INPUT, "No sequences were read from the input.");
    }

    std::sort(seqsLens.begin(), seqsLens.end());

    uint32_t avgLen = totalLen/(seqNum - seqNum_init);
    uint32_t medLen = seqsLens[(seqNum - seqNum_init)/2];

    bool min_max = (option->minLen > 0 || option->maxLen < INT32_MAX);

    int minLenTh = (min_max) ? option->minLen : (option->lenDev > 0) ? static_cast<int>(medLen*(1-option->lenDev)) : option->minLen;
    int maxLenTh = (min_max) ? option->maxLen : (option->lenDev > 0) ? static_cast<int>(medLen*(1+option->lenDev)) : option->maxLen;
    std::atomic<int> numLowQ;
    tbb::spin_rw_mutex lowQMutex;
    stringPairVec lowQSeqs;
    numLowQ.store(0);
    if (option->alnMode != 3 || placed) {
        tbb::parallel_for(tbb::blocked_range<int>(0, seqNum), [&](tbb::blocked_range<int> range){ 
        for (int i = range.begin(); i < range.end(); ++i) {
            auto seq = database->sequences[i];
            if (option->alnMode == 3 && !tree->allNodes[seq->name]->placed) continue;
            int ambig = (option->type == 'n') ? 4 : 20;
            seq->lowQuality = (seq->len > maxLenTh || seq->len < minLenTh);
            if (!seq->lowQuality) {
                int ambigCount = 0;
                for (int j = 0; j < seq->len; ++j) {
                    if (letterIdx(option->type, toupper(seq->alnStorage[0][j])) == ambig) ambigCount++;
                }
                seq->lowQuality = (ambigCount > (seq->len * option->maxAmbig));
            }
            if (seq->lowQuality) {
                numLowQ.fetch_add(1);
                if ((!option->noFilter && option->writeFiltered)) {
                    {
                        tbb::spin_rw_mutex::scoped_lock lock(lowQMutex);
                        lowQSeqs.push_back({seq->name, std::string(seq->alnStorage[0],seq->len)});
                    }
                }
                if (!option->noFilter) database->sequences[i]->len = 0;
            }
        }
        });
    }
    if (!lowQSeqs.empty()) {
        fs::path o(option->outFile);
        fs::path p(option->seqFile);
        std::string outPath = (o.parent_path().string() == "") ? "." : o.parent_path().string(); 
        std::string lowQFileName = (subtree != -1) ? 
            outPath + "/subtree-" + std::to_string(subtree) + ".filtered.fasta" :
            outPath + "/" + p.stem().string() + ".filtered.fasta";
        writeAlignment(lowQFileName, lowQSeqs, option->compressed, false);
    }
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    
    // Print summary
    int numRead = (seqNum - seqNum_init);
    if (option->accurate) {
        option->resolveAccurateParams(numRead);
        if (option->printDetail) {
            std::string modeDesc;
            if (option->sampleMode == 0) {
                modeDesc = "0 (auto -> mode " + std::to_string(option->resolvedAccMode) + ")";
            } else if (option->sampleMode == 1) {
                modeDesc = "1 (all-to-all)";
            } else if (option->sampleMode == 2) {
                modeDesc = "2 (high sample rate)";
            } else {
                modeDesc = "3 (low sample rate)";
            }
            msg::debug("[Accurate Mode] Resolved parameters for %d sequences: sample-mode=%s, sample-max-group=%d, sample-rate=%g, sample-top-centrality=%g\n",
                       numRead, modeDesc.c_str(), option->accMaxGroup, option->accSampleRate, option->accSampleTop);
        }
    }
    if (option->alnMode != 3 || placed) {
        msg::info("===== Sequence Summary =====\n");
        msg::info("Number : %d\n", numRead);
        msg::info("Max. Length: %d\n", maxLen);
        msg::info("Min. Length: %d\n", minLen);
        msg::info("Avg. Length: %u\n", avgLen);
        msg::info("Med. Length: %u\n", medLen);
        if (option->noFilter)
            msg::info("Deferred sequences: %d\n", numLowQ.load());
        else
            msg::info("Excluded sequences: %d\n", numLowQ.load());
        msg::info("Sequences read in %lld ms\n", static_cast<long long>(seqReadTime.count() / 1000000));
    }
    else {
        msg::info("==== Backbone Alignment ====\n");
        msg::info("Number : %d\n", numRead);
        msg::info("Length:  %u\n", avgLen);
        msg::info("Backbone alignment read in %lld ms\n", static_cast<long long>(seqReadTime.count() / 1000000));
    }
}

void msa::io::readAlignment(std::string msaFileName, SequenceDB* database, Option* option, Node*& node) {
    int profileSize = option->type == 'n' ? 6 : 22;
    gzFile f_rd = gzopen(msaFileName.c_str(), "r");
    if (!f_rd) {
        throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Can't open file: " + msaFileName);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    int seqNum = 0, msaLen = 0;
    std::string seqName;
    while (kseq_read(kseq_rd) >= 0) {
        int seqLen = kseq_rd->seq.l;
        bool calFreq = true;
        if (seqNum == 0) {
            msaLen = seqLen;
            seqName = kseq_rd->name.s;
            node->msaFreq = Profile (msaLen, std::vector<float> (profileSize, 0.0));
        }
        else {
            if (seqLen != msaLen) {
                msg::warn("Length of \"%s\" (%d) does not match in %s (%d)\n", (kseq_rd->name.s), seqLen, msaFileName.c_str(), msaLen);
                calFreq = false;
            }
        }
        if (!calFreq) continue;
        std::string seq = std::string(kseq_rd->seq.s, seqLen);
        tbb::parallel_for(0, seqLen, [&](int j) {
            int letterIndex = letterIdx(option->type, toupper(seq[j]));
            node->msaFreq[j][letterIndex] += 1.0;
        });
        ++seqNum;
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);
    node->alnNum = seqNum;
    node->alnLen = msaLen;
    node->alnWeight = static_cast<float>(seqNum);
    return;
}

phylogeny::Tree* msa::io::readAlignments_and_buildTree(SequenceDB *database, Option *option) {
    auto alnReadStart = std::chrono::high_resolution_clock::now();
    std::string path = option->msaDir;
    std::vector<std::string> files;
    int subtreeIdx = 0;
    int totalFile = 0;
    // Read alignment file names
    std::error_code ec;
    fs::path msaFolder(path);
    for (fs::recursive_directory_iterator it(msaFolder, ec), eit; it != eit; it.increment(ec)) {
        if (ec) {
            it.pop();
            continue;
        }
        if (!fs::is_directory(it->path())) {
            files.push_back(it->path().string());
            totalFile += 1;
        }
    }
    // Read alignments and calculate profiles
    msg::info("====== Alignment Info ======\n");
    std::sort(files.begin(), files.end());
    std::vector<Node*> nodes;
    int profileSize = option->type == 'n' ? 6 : 22;
    for (auto msaFileName: files) {
        std::string nodeName = "node_" + std::to_string(subtreeIdx+1);
        Node* node = new Node(nodeName, 1.0);
        node->grpID = 0;
        node->seqsIncluded.push_back(subtreeIdx);
        readAlignment(msaFileName, database, option, node);
        database->subtreeAln[subtreeIdx] = alnPath (node->msaFreq.size(), 0);
        nodes.push_back(node);
        database->subAlnFiles.push_back({msaFileName, subtreeIdx});
        ++subtreeIdx;
        msg::info("[%d/%zu] %s (Count: %d, Length: %d)\n",
                  subtreeIdx, files.size(), fs::path(msaFileName).filename().string().c_str(),
                  node->alnNum, node->alnLen);
    }
    // Sort alignments by sequence count
    std::sort(nodes.begin(), nodes.end(), [](Node* a, Node* b) {
        return a->alnNum > b->alnNum;  // descending order
    });
    // Build a Tree
    Tree* T = new Tree();
    T->root = nodes.front();
    T->allNodes[nodes.front()->identifier] = nodes.front();
    for (int i = 1; i < nodes.size(); ++i) {
        nodes[i]->parent = T->root;
        T->root->children.push_back(nodes[i]);
        T->allNodes[nodes[i]->identifier] = nodes[i];
    }
    phylogeny::updateLevels(T->root, 1);
    auto alnReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnReadTime = alnReadEnd - alnReadStart;
    
        
    msg::info("Read %zu MSA files in %g ms.\n", T->allNodes.size(), alnReadTime.count() / 1e6);
    return T;
}

void msa::io::readBackboneAlignment(Tree* T, SequenceDB *database, Option *option) {
    auto alnReadStart = std::chrono::high_resolution_clock::now();
    int subtreeIdx = -1;
    readAlignment(option->backboneAlnFile, database, option, T->root);
    database->subtreeAln[subtreeIdx] = alnPath (T->root->msaFreq.size(), 0);
    T->root->seqsIncluded.push_back(subtreeIdx);
    auto alnReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnReadTime = alnReadEnd - alnReadStart;
    msg::info("Backbone alignment read in %lld ms. (Count: %d, Length: %zu)\n",
              static_cast<long long>(alnReadTime.count() / 1000000), T->root->alnNum, T->root->msaFreq.size());
    return;
}

void msa::io::writePrunedTree(Tree* T, msa::Option* option) {
    fs::path o(option->outFile);
    fs::path t(option->treeFile);
    std::string prunedTreeFileName = (o.parent_path().string() == "") ? t.filename().string() + ".pruned.nwk" : o.parent_path().string() + "/" + t.filename().string() + ".pruned.nwk";
    std::string out_str = T->getNewickString() + "\n";
	std::ofstream outFile(prunedTreeFileName);
    if (!outFile) {
        throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Can't open file: " + prunedTreeFileName);
    }
	outFile << out_str;
	outFile.close();
    return;
}

void msa::io::writeSubtrees(Tree* tree, PartitionInfo* partition, Option* option) {
    std::string tempDir = option->tempDir;
    for (auto subroot: partition->partitionsRoot) {
        int subtreeIdx = tree->allNodes[subroot.first]->grpID;
        Tree* subT = new Tree(subroot.second.first, false);
        std::string subtreeName = "subtree-" + std::to_string(subtreeIdx);
        std::string subtreeTreeFile = option->tempDir + '/' + subtreeName + ".nwk";
        std::string out_str = subT->getNewickString() + "\n";
        std::ofstream outFile(subtreeTreeFile);
        if (!outFile) {
            delete subT;
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Failed to open file: " + subtreeTreeFile);
        }
	    outFile << out_str;
	    outFile.close();
        delete subT;
    }
    return;
}

void msa::io::writeSubAlignments(SequenceDB* database, Option* option, int subtreeIdx, int alnLen) {
    std::string fileName = option->tempDir + "/subtree-" + std::to_string(subtreeIdx) + ".subalignment.aln";
    database->subAlnFiles.push_back({fileName, subtreeIdx}); // Store file name and subtreeIdx
    std::vector<std::pair<std::string, std::string>> seqs;
    msg::info("Subalignment Length: %d\n", alnLen);
    writeAlignment(fileName, database, alnLen, option->compressed);
    return;
}

int msa::io::update_and_writeAlignment(SequenceDB* database, Option* option, std::string fileName, int subtreeIdx) {
    const int outBuffSize = 10000;
    int totalSeqs = 0;
    bool nochange = false;
    char gapType = (option->alnMode == PLACE_WO_TREE) ? '.' : '-';
    stringPairVec seqs_before, seqs_after;
    if ((option->alnMode == DEFAULT_ALN || option->alnMode == ACCURATE_ALN) && option->compressed) fileName += ".gz";
    fs::path p(fileName);
    std::string filename = p.stem().string();
    std::string finalAlnFileName = option->tempDir + "/" + filename + ".final.aln";
    gzFile f_rd = gzopen(fileName.c_str(), "r");
    if (!f_rd) {
        throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Failed to open file: " + fileName);
    }
    kseq_t* kseq_rd = kseq_init(f_rd);
    while (kseq_read(kseq_rd) >= 0) {
        std::string seqName = kseq_rd->name.s;
        size_t seqLen = kseq_rd->seq.l;
        if (seqLen == database->subtreeAln[subtreeIdx].size()) {
            // No need to update the alignment since no insertions are introduced
            nochange = true;
            break;
        }
        seqs_before.push_back({seqName, std::string(kseq_rd->seq.s, seqLen)});
        if (seqs_before.size() == outBuffSize) { // Update Sequences
            seqs_after.resize(seqs_before.size());
            tbb::parallel_for(tbb::blocked_range<int>(0, seqs_before.size()), [&](tbb::blocked_range<int> range){ 
            for (int n = range.begin(); n < range.end(); ++n) {
                seqs_after[n].first = seqs_before[n].first;
                seqs_after[n].second.reserve(database->subtreeAln[subtreeIdx].size());
                size_t r = 0;         
                for (size_t j = 0; j < database->subtreeAln[subtreeIdx].size(); ++j) {
                    if ((database->subtreeAln[subtreeIdx][j] & 0xFFFF) == 0) {
                        seqs_after[n].second += seqs_before[n].second[r];
                        ++r;
                    }
                    else {
                        seqs_after[n].second += gapType;
                    }
                }
                if (seqs_after[n].second.size() != database->subtreeAln[subtreeIdx].size()) {
                    msg::error("%s length not match. alnSeq (%zu) != msaAln (%zu)\n", seqName.c_str(), seqs_after[n].second.size(), database->subtreeAln[subtreeIdx].size());
                }
            }
            });
            if (totalSeqs == 0) writeAlignment(finalAlnFileName, seqs_after, option->compressed, false);  
            else                writeAlignment(finalAlnFileName, seqs_after, option->compressed, true);   
            totalSeqs += seqs_after.size();
            seqs_after.clear();
            seqs_before.clear();
        }
    }
    kseq_destroy(kseq_rd);
    gzclose(f_rd);
    // Update the remaining sequences
    if (!nochange) { 
        seqs_after.resize(seqs_before.size());
        tbb::parallel_for(tbb::blocked_range<int>(0, seqs_before.size()), [&](tbb::blocked_range<int> range){ 
        for (int n = range.begin(); n < range.end(); ++n) {
            seqs_after[n].first = seqs_before[n].first;
            seqs_after[n].second.reserve(database->subtreeAln[subtreeIdx].size());
            size_t r = 0;         
            for (size_t j = 0; j < database->subtreeAln[subtreeIdx].size(); ++j) {
                if ((database->subtreeAln[subtreeIdx][j] & 0xFFFF) == 0) {
                    seqs_after[n].second += seqs_before[n].second[r];
                    ++r;
                }
                else {
                    seqs_after[n].second += gapType;
                }
            }
            if (seqs_after[n].second.size() != database->subtreeAln[subtreeIdx].size()) {
                msg::error("%s length not match. alnSeq (%zu) != msaAln (%zu)\n", seqs_after[n].first.c_str(), seqs_after[n].second.size(), database->subtreeAln[subtreeIdx].size());
            }
        }
        });
        if (totalSeqs == 0) writeAlignment(finalAlnFileName, seqs_after, option->compressed, false);  
        else                writeAlignment(finalAlnFileName, seqs_after, option->compressed, true);   
        totalSeqs += seqs_after.size();
        seqs_after.clear();
        seqs_before.clear();
    }
    if (option->alnMode == PLACE_WO_TREE) msg::info("Final Alignment Length: %zu\n", database->subtreeAln[subtreeIdx].size());
    database->subtreeAln[subtreeIdx].clear();
    if (nochange) {
        std::string command = "cp " + fileName + " " + finalAlnFileName;
        int delResult = system(command.c_str());
        if (delResult != 0) msg::error("Unable to copy %s to %s.\n", fileName.c_str(), finalAlnFileName.c_str());
    }
    if (option->deleteTemp && (option->alnMode == 0 || option->alnMode == ACCURATE_ALN)) {
        std::string command = "rm " + fileName;
        int delResult = system(command.c_str());
        if (delResult != 0) msg::error("Unable to delete %s.\n", fileName.c_str());
    }
    return totalSeqs;
}

void msa::io::update_and_writeAlignments(SequenceDB* database, Option* option, int& totalSeqs) {
    int proceeded = 0;   
    totalSeqs = 0;
    for (auto subAln: database->subAlnFiles) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        totalSeqs += update_and_writeAlignment(database, option, subAln.first, subAln.second);
        ++proceeded;
        auto subtreeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
        msg::info("Finish writing alignment of subtree No. %d. (%d/%zu) in %lld ms.\n",
                  subAln.second, proceeded, database->subAlnFiles.size(), static_cast<long long>(subtreeTime.count() / 1000000));
    }
    return;
}

void msa::io::writeFinalMSA(SequenceDB* database, Option* option, int alnLen, const Params* param) {
    std::string fileName = option->outFile;
    if (database->currentTask == 2) {
        if (option->compressed) fileName += ".gz";
        std::string command = "cat " + option->tempDir + "/*.final.aln* > " + fileName;
        int catResult = system(command.c_str());
        if (catResult != 0) {
            msg::error("Unable to concatenate subtree alignments. Please check the available space on your device and merge files with \"cat <temp-dir>/*.final.aln <output-file>\". Temporary files will not be deleted.\n");
        }
        if (option->scoreMode > 0) {
            msa::score::evaluateAndOutputScore(database, option, param, alnLen);
        }
        if (catResult == 0 && option->deleteTemp) {
            command = "rm -rf " + option->tempDir;
            int delResult = system(command.c_str());
            if (delResult != 0) msg::error("Unable to delete temporary files.\n");
        }
        return;
    }
    if (option->alnMode != 2) {
        std::vector<std::pair<std::string,std::string>> seqs;
        msg::info("Final Alignment Length: %d\n", alnLen);
        writeAlignment(fileName, database, alnLen, option->compressed);
        if (option->scoreMode > 0) {
            msa::score::evaluateAndOutputScore(database, option, param, alnLen);
        }
        if (option->deleteTemp && option->tempDir != "") {
            std::string command = "rm -rf " + option->tempDir;
            int delResult = system(command.c_str());
            if (delResult != 0) msg::error("Unable to delete temporary files.\n");
        }
    }
    return;
}

void msa::io::writeAlignment(std::string fileName, SequenceDB* database, int alnLen, bool compressed) {
    if (compressed) {
        fileName += ".gz";
        std::ofstream outFile(fileName, std::ios::binary);
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Failed to open file: " + fileName);
        }
        // Batch chunking: compress sequences in parallel batches to keep threads busy while bounding memory usage
        const size_t totalSeqs = database->sequences.size();
        const size_t batchSize = 1000;
        for (size_t bStart = 0; bStart < totalSeqs; bStart += batchSize) {
            size_t bEnd = std::min(bStart + batchSize, totalSeqs);
            size_t curBatchSize = bEnd - bStart;
            std::vector<std::string> compressed_chunks(curBatchSize);
            tbb::parallel_for(size_t(0), curBatchSize, [&](size_t idx) {
                size_t s = bStart + idx;
                auto seq = database->sequences[s];
                if (seq && !seq->lowQuality) {
                    std::string fasta_chunk = ">" + seq->name + "\n" + std::string(&seq->alnStorage[seq->storage][0], alnLen) + "\n";
                    compressed_chunks[idx] = gzip_compress(fasta_chunk);
                }
            });
            for (const auto &c : compressed_chunks) {
                if (!c.empty()) outFile.write(c.data(), c.size());
            }
        }
        outFile.close();
    }
    else {
        std::ofstream outFile(fileName);
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Failed to open file: " + fileName);
        }
        for (int s = 0; s < database->sequences.size(); ++s) {
            auto seq = database->sequences[s];
            if (!seq->lowQuality) {
                outFile << '>' << seq->name << "\n";
                outFile.write(&seq->alnStorage[seq->storage][0], alnLen);
                outFile << '\n';
            }
        }
        outFile.close();
    }
    return;
}

void msa::io::writeAlignment(std::string fileName, stringPairVec& seqs, bool compressed, bool append) {
    if (compressed) {
        fileName += ".gz";
        std::ofstream outFile;
        if (append) outFile.open(fileName, std::ios::binary | std::ios::app);
        else        outFile.open(fileName, std::ios::binary);
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Can't open file: " + fileName);
        }
        // Batch chunking: compress sequences in parallel batches to keep threads busy while bounding memory usage
        const size_t totalSeqs = seqs.size();
        const size_t batchSize = 1000;
        for (size_t bStart = 0; bStart < totalSeqs; bStart += batchSize) {
            size_t bEnd = std::min(bStart + batchSize, totalSeqs);
            size_t curBatchSize = bEnd - bStart;
            std::vector<std::string> compressed_chunks(curBatchSize);
            tbb::parallel_for(size_t(0), curBatchSize, [&](size_t idx) {
                size_t i = bStart + idx;
                std::string fasta_chunk = ">" + seqs[i].first + "\n" + seqs[i].second + "\n";
                compressed_chunks[idx] = gzip_compress(fasta_chunk);
                std::string().swap(seqs[i].second);
            });
            for (const auto &c : compressed_chunks) {
                if (!c.empty()) outFile.write(c.data(), c.size());
            }
        }
        outFile.close();
    }
    else {
        std::ofstream outFile;
        if (append) outFile.open(fileName, std::ios::app);
        else        outFile.open(fileName);
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR, "Can't open file: " + fileName);
        }
        for (auto it = seqs.begin(); it != seqs.end(); ++it) {
            outFile << ('>' + it->first + '\n');
            outFile << (it->second + '\n');
        }
        outFile.close();
    }
    return;
}
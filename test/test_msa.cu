
#include <fstream>
#include <string>
#include <boost/program_options.hpp> 
#include "../src/kseq.h"
#include "zlib.h"

#ifndef TREE_HPP
#include "../src/tree.hpp"
#endif

#ifndef MSA_HPP
#include "../src/msa.hpp"
#endif

#ifndef ALIGN_HPP
#include "../src/align.cuh"
#endif

#ifndef ALIGN_TALCO_HPP
#include "../src/TALCO-XDrop.hpp"
#endif


namespace po = boost::program_options;

#include "../src/treePartition.cpp"



KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("MSA Command Line Arguments");


void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,s", po::value<std::string>()->required(), "Tip sequences - Fasta format (required)")
        ("machine,m",  po::value<std::string>()->default_value("gpu"), "Run on gpu or cpu")
        ("max_leaves,l",  po::value<int>()->default_value(200), "Maximum number of leaves per subtree")
        ("help,h", "Print help messages");

}


void printTree(Node* node)
{
    if (node->parent == nullptr)
        std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << "ROOT\t" << node->grpID << std::endl;
    else
        std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << "\t" << node->grpID << std::endl;  

    if (node->children.size() == 0) return;
    // std::cout << "Print children\n";
    for (auto &c: node->children) printTree(c);
}

void printLeaves(Node* node)
{
    if (node->children.size() == 0) {
        std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << '\t' << node->grpID << std::endl;
        return;
    }
    // else
    //     std::cout << std::setw(7) << node->identifier << ": " << node->branchLength << "\t" << node->parent->identifier << "\t" << node->grpID << std::endl;  

    // if (node->children.size() == 0) return;
    // std::cout << "Print children\n";
    for (auto &c: node->children) printLeaves(c);
}

Tree* readNewick(po::variables_map& vm)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::string treeFileName = vm["tree"].as<std::string>();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Can't open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; inputStream >> newick;
    Tree *T = new Tree(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Newick string read in: " <<  treeBuiltTime.count() << " ns\n";

    // printTree(T->root);

    return T;
}

void readSequences(po::variables_map& vm, msa::utility* util)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName = vm["sequences"].as<std::string>();

    gzFile f_rd = gzopen(seqFileName.c_str(), "r");
    if (!f_rd) {
        fprintf(stderr, "ERROR: cant open file: %s\n", seqFileName.c_str());
        exit(1);
    }

    kseq_t* kseq_rd = kseq_init(f_rd);

    while (kseq_read(kseq_rd) >= 0) {
        size_t seqLen = kseq_rd->seq.l;
        util->seqs[kseq_rd->name.s] = std::string(kseq_rd->seq.s, seqLen);
    }
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in: " <<  seqReadTime.count() << " ns\n";
}

void checkAlignment(std::vector<std::string>& ref)
{
    size_t len = 0;
    bool set = false;

    for (auto &r: ref)
    {
        if (!set) len = r.size();
        if (r.size() != len)
        {
            fprintf(stderr, "Error: Alignment Size do not match\n");
        }
    }
}

/*
void msaPostOrderTraversal(Node* node, msa::utility* util, Params& param, int grpID)
{
    // std::cout << node->identifier << '\n';
    // if (node->children.size() == 0) std::cout << node->identifier << '\n';
    if (!(node->grpID==-1 || node->grpID==grpID)) return;

    if (node->children.size()==0) 
    {
        node->msa.push_back(util->seqs[node->identifier]);
        return;
    }
    // for (auto& child: node->children) msaPostOrderTraversal(child, util, param, grpID);
    std::pair<std::vector<std::string>, std::vector<std::string>> alignments;
    std::vector<std::string> ref;
    
    size_t childIndex = 0;
    for (childIndex=0; childIndex<node->children.size(); childIndex++)
    {
        if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->msa.size() > 0))
        {
            ref = node->children[childIndex]->msa;
            break;
        }
    }

    


    
    // std::cout << node->identifier << '\n';
    if (childIndex == node->children.size() - 1) {
        node->msa = ref;
        return;
    }

    // if (node->identifier == "node_13") {
    for (size_t i=childIndex+1; i<node->children.size(); i++)
    {
        if (!(node->children[i]->grpID == -1 || node->children[i]->grpID == grpID))
        {
            continue;
        }
        
        // if (node->identifier == "node_11") {
        auto alnStart = std::chrono::high_resolution_clock::now();
        std::vector<std::string> query = node->children[i]->msa;
        if (ref.size() == 0 || query.size() == 0) continue;
        // Malloc 
        int32_t seqLen = ref[0].size() + query[0].size();
        int32_t seqNum = ref.size() + query.size();
        char* hostRef = (char*)malloc(seqLen * seqNum * sizeof(char));
        char* hostQry = (char*)malloc(seqLen * seqNum * sizeof(char));
        char* hostAln = (char*)malloc(seqLen * seqNum * sizeof(char));
        int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
        int16_t* hostParam = (int16_t*)malloc(6 * sizeof(int16_t)); 
        // float* hostFreqRef = (float*)malloc(seqLen * 5 * sizeof(float));
        // float* hostFreqQry = (float*)malloc(seqLen * 5 * sizeof(float));
        // int32_t* hostH     = (int32_t*)malloc(3 * seqLen * 2 * sizeof(int32_t));
        // int32_t* hostI     = (int32_t*)malloc(2 * seqLen * 2 * sizeof(int32_t));
        // int32_t* hostD     = (int32_t*)malloc(2 * seqLen * 2 * sizeof(int32_t));
        // int32_t* hostWfLL  = (int32_t*)malloc(seqLen * 2 * sizeof(int32_t));
        // int32_t* hostWfLen = (int32_t*)malloc(seqLen * 2 * sizeof(int32_t));
        // int8_t* hostTB = (int8_t*)malloc(seqLen*seqLen * sizeof(int8_t)); 

        int seqCount = 0;
        for (int j = 0; j < seqLen*seqNum; ++j) { 
            hostAln[j] = 0;
            if (seqCount < ref.size()) {
                if (j%seqLen < ref[seqCount].size()) {
                    hostRef[j] = ref[seqCount][j%seqLen];
                }
                else hostRef[j] = 0;
            } 
            else hostRef[j] = 0;
            if (seqCount < query.size()) {
                if (j%seqLen < query[seqCount].size()) {
                    hostQry[j] = query[seqCount][j%seqLen];
                }
                else hostQry[j] = 0;
            }
            else hostQry[j] = 0;
            if (j%seqLen == seqLen-1) ++seqCount;
        }
        
        

        int numBlocks = 1; 
        int blockSize = 256;
        
        hostSeqInfo[0] = seqLen;
        hostSeqInfo[1] = ref[0].size();
        hostSeqInfo[2] = query[0].size();
        hostSeqInfo[3] = ref.size();
        hostSeqInfo[4] = query.size();
        hostSeqInfo[5] = numBlocks;
        hostSeqInfo[6] = blockSize;
        
        hostParam[0] = param.match;
        hostParam[1] = param.mismatch;
        hostParam[2] = param.gapOpen;
        hostParam[3] = param.gapExtend;
        hostParam[4] = param.xdrop;
        hostParam[5] = param.marker;
        // Cuda Malloc
        char* deviceRef;
        char* deviceQry;
        char* deviceAln;
        int32_t* deviceSeqInfo;
        int16_t* deviceParam;

        cudaMalloc((void**)&deviceRef, seqLen * seqNum * sizeof(char));
        cudaMalloc((void**)&deviceQry, seqLen * seqNum * sizeof(char));
        cudaMalloc((void**)&deviceAln, seqLen * seqNum * sizeof(char));
        cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
        cudaMalloc((void**)&deviceParam, 6 * sizeof(int16_t));

        // Copy to Device
        cudaMemcpy(deviceRef, hostRef, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceQry, hostQry, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceAln, hostAln, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceParam, hostParam, 6 * sizeof(int16_t), cudaMemcpyHostToDevice);
    
        
    
        // printf("Before kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        printf("refSize: %d, qrySize: %d\n", ref.size(), query.size());
        printf("refLen: %d, qryLen: %d\n", ref[0].size(), query[0].size());
        alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
            deviceRef, 
            deviceQry, 
            deviceParam, 
            deviceAln, 
            deviceSeqInfo
        );

        // alignGrpToGrp_cuda<<<numBlocks, blockSize>>>(
        //     deviceRef, 
        //     deviceQry, 
        //     deviceParam, 
        //     deviceAln, 
        //     deviceSeqInfo
        // );
        cudaDeviceSynchronize();
        // printf("After kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        cudaError_t err;
        // err = cudaMemcpy(hostRef, deviceRef, seqLen * seqNum * sizeof(char), cudaMemcpyDeviceToHost);
        // cudaGetErrorString(err);
        // err = cudaMemcpy(hostQry, deviceQry, seqLen * seqNum * sizeof(char), cudaMemcpyDeviceToHost);
        // cudaGetErrorString(err);
        err = cudaMemcpy(hostAln, deviceAln, seqLen * seqNum * sizeof(char), cudaMemcpyDeviceToHost);
        cudaGetErrorString(err);
        cudaDeviceSynchronize();
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Aln "<<  node->identifier << "," << node->children[i]->identifier <<" in: " <<  alnTime.count() << " ns\n";
        ref.clear();


        
        for (int j = 0; j < (hostSeqInfo[3] + hostSeqInfo[4]); ++j) {
            std::string s = "";
            size_t alnIdx = j * seqLen;
            // std::cout << j << '\n'; 
            while (true) {
                if (hostAln[alnIdx] != 'A' && hostAln[alnIdx] != 'C' && hostAln[alnIdx] != 'G' && 
                    hostAln[alnIdx] != 'T' && hostAln[alnIdx] != '-' && hostAln[alnIdx] != 'N') break;
                // std::cout << hostAln[alnIdx];
                s += hostAln[alnIdx];
                ++alnIdx;
            }
            // std::reverse(s.begin(), s.end()); 
            // std::cout << s << '\n';
            ref.push_back(s);
        }        
        // for (int i = 0; i < ref[0].size(); ++i) {
        //     bool same = true;
        //     char firstchar = ref[0][i];
        //     for (int j = 1; j < ref.size(); ++j) {
        //         if ( ref[j][i] != firstchar) same = false;
        //     }
        //     if (!same || i < 100) {
        //         std::cout << "===POSITION  "<< i <<"  =================\n";
        //         for (int j = 0; j < ref.size(); ++j) {
        //             std::cout << ref[j].substr(max(0, i - 10), 20) << '\n';
        //         }
                
        //     }
        // }
        std::cout << "Aln Seq Length: " << ref[0].size() << '\n';
        
        // for (auto &s: alignments.first) ref.push_back(s);
        // for (auto &s: alignments.second) ref.push_back(s);
        // alignments.first.clear(); alignments.second.clear();

        // if (node->identifier == "node_6") {
        //     for (int i = 0; i < 64; ++i) {
        //         std::cout << hostRef[i];
        //     }
        //     std::cout << '\n';
        //     for (int i = 0; i < 64; ++i) {
        //         std::cout << hostQry[i];
        //     }
        //     std::cout << '\n'; 
        //     // for (int i = 0; i < 20; ++i) {
        //     //     std::cout << hostQry[seqLen+i];
        //     // }
        //     // std::cout << '\n'; 
        // }

        // free device memory
        cudaFree(deviceRef);
        cudaFree(deviceQry);
        cudaFree(deviceAln);
        cudaFree(deviceParam);
        cudaFree(deviceSeqInfo);

        cudaDeviceSynchronize();
        // free host memory
        free(hostRef);
        free(hostQry);
        free(hostAln);
        free(hostParam);
        free(hostSeqInfo);

        


        // alignGrpToGrp(ref, query, param, alignments);
        // ref.clear();
        // for (auto &s: alignments.first) ref.push_back(s);
        // for (auto &s: alignments.second) ref.push_back(s);
        // alignments.first.clear(); alignments.second.clear();
        // }
    // }
    }
    // if (node->identifier == "node_6") {
    //     std::cout << "RESULT...\n"; 
    // for (int j = 0; j < ref.size(); ++j) {
    //     for (int k = 0; k < 64; ++k) {
    //         std::cout << ref[j][k];
    //     }
    //     std::cout << '\n';
    // }    
    // }
    node->msa = ref;

    return;
}
*/

/*
void msaPostOrderTraversal_gpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{
    // std::cout << node->identifier << '\n';
    // if (node->children.size() == 0) std::cout << node->identifier << '\n';

    // assign msa to all nodes
    for (auto n: nodes) {
        // std::cout << n.first->identifier << '\n';
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msa.push_back(util->seqs[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
        // std::cout << n.second->identifier << '\n';
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msa.push_back(util->seqs[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
    }

    // store all sequences to array
    int32_t seqNum = 0;
    int32_t seqLen = 0;
    int32_t pairNum = nodes.size();
    std::vector<std::string> seqs;
    std::vector<std::vector<uint8_t>> freq;
    std::vector<std::pair<int32_t, int32_t>> seqIdx;
    std::vector<std::pair<int32_t, int32_t>> len;
    // get maximum sequence/profile length 
    for (auto n: nodes) {
        int32_t qryLen = tree->allNodes[n.second->identifier]->msa[0].size();
        int32_t refLen = tree->allNodes[n.first->identifier]->msa[0].size();
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    // store info to array 
    for (auto n: nodes) {
        int32_t qryIdx = 0;
        int32_t refIdx = 0;
        int32_t qryLen = tree->allNodes[n.second->identifier]->msa[0].size();
        int32_t refLen = tree->allNodes[n.first->identifier]->msa[0].size();
        refIdx = seqNum;
        std::vector<uint8_t> temp;
        for (int i = 0; i < 12*seqLen; ++i) temp.push_back(0);
        assert(temp.size() == 12*seqLen);
        
        for (auto seq: tree->allNodes[n.first->identifier]->msa) {
            for (int s = 0; s < refLen; ++s) {
                if      (seq[s] == 'A' || seq[s] == 'a') temp[6*s+0]+=1;
                else if (seq[s] == 'C' || seq[s] == 'c') temp[6*s+1]+=1;
                else if (seq[s] == 'G' || seq[s] == 'g') temp[6*s+2]+=1;
                else if (seq[s] == 'T' || seq[s] == 't') temp[6*s+3]+=1;
                else if (seq[s] == 'N' || seq[s] == 'n') temp[6*s+4]+=1;
                else                                     temp[6*s+5]+=1;
            }
            ++seqNum;
            seqs.push_back(seq);
        }
        qryIdx = seqNum;
        for (auto seq: tree->allNodes[n.second->identifier]->msa) {
            for (int s = 0; s < qryLen; ++s) {
                if      (seq[s] == 'A' || seq[s] == 'a') temp[6*(seqLen+s)+0]+=1;
                else if (seq[s] == 'C' || seq[s] == 'c') temp[6*(seqLen+s)+1]+=1;
                else if (seq[s] == 'G' || seq[s] == 'g') temp[6*(seqLen+s)+2]+=1;
                else if (seq[s] == 'T' || seq[s] == 't') temp[6*(seqLen+s)+3]+=1;
                else if (seq[s] == 'N' || seq[s] == 'n') temp[6*(seqLen+s)+4]+=1;
                else                                     temp[6*(seqLen+s)+5]+=1;
            }
            ++seqNum;
            seqs.push_back(seq);
        }
        seqIdx.push_back(std::make_pair(refIdx, qryIdx));
        len.push_back(std::make_pair(refLen, qryLen));
        freq.push_back(temp);
        
    }

    // Print alignment info
    // for (int i = 0; i < pairNum; ++i) {
    //     printf("No.%d [%s, %s] (%d, %d)\n", i, nodes[i].first->identifier.c_str(), nodes[i].second->identifier.c_str(), len[i].first, len[i].second);    
    // }
    
    // Malloc
    // char* hostSeqs = (char*)malloc(seqLen * seqNum * sizeof(char));
    uint8_t* hostFreq = (uint8_t*)malloc(12*seqLen * pairNum * sizeof(uint8_t));
    int8_t* hostAln = (int8_t*)malloc(2*seqLen * pairNum * sizeof(int8_t));
    // int32_t* hostIdx = (int32_t*)malloc(2*pairNum * sizeof(int32_t));
    int32_t* hostLen = (int32_t*)malloc(2*pairNum * sizeof(int32_t));
    int32_t* hostAlnLen = (int32_t*)malloc(pairNum * sizeof(int32_t));
    int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
    int16_t* hostParam = (int16_t*)malloc(6 * sizeof(int16_t)); 
    
    // copy info from array host mem
    // int seqCount = 0;
    // for (int j = 0; j < seqLen*seqNum; ++j) { 
    //     if (seqCount < seqs.size()) {
    //         if (j%seqLen < seqs[seqCount].size()) {
    //             hostSeqs[j] = seqs[seqCount][j%seqLen];
    //         }
    //         else hostSeqs[j] = 0;
    //     }
    //     if (j%seqLen == seqLen-1) ++seqCount;
    // }
    for (int j = 0; j < 2*pairNum; ++j) { 
        // if (j%2 == 0) hostIdx[j] = seqIdx[j/2].first;
        // else          hostIdx[j] = seqIdx[j/2].second;
        if (j%2 == 0) hostLen[j] = len[j/2].first;
        else          hostLen[j] = len[j/2].second;
    }
    for (int j = 0; j < pairNum; ++j) {
        for (int l = 0; l < 12*seqLen; ++l) {
            hostFreq[12*seqLen*j+l] = freq[j][l];
        }
    }
    // for (int i = 0; i < 2*pairNum; ++i) {
    //     std::cout << hostIdx[i] << ',';    
    // }
    // std::cout << '\n';
    // for (int i = 0; i < 2*pairNum; ++i) {
    //     std::cout << hostLen[i] << ',';    
    // }
    // std::cout << '\n';
    // for (int i = 0; i < 70; ++i) std::cout << tree->allNodes[nodes[0].first->identifier]->msa[0][i] << ',';
    // std::cout << '\n';
    // for (int i = 0; i < 70; ++i) std::cout << (hostFreq[6*i+3] & 0xFFFF) << ',';
    // std::cout << '\n';


    for (int j = 0; j < 2*seqLen*pairNum; ++j) { 
        hostAln[j] = 0;
    }
    for (int j = 0; j < pairNum; ++j) { 
        hostAlnLen[j] = 0;
    }

    int numBlocks = 1024; 
    int blockSize = 512;

    hostSeqInfo[0] = seqLen;
    hostSeqInfo[1] = seqNum;
    hostSeqInfo[2] = pairNum;
    hostSeqInfo[3] = numBlocks;
    hostSeqInfo[4] = blockSize;
    hostSeqInfo[5] = numBlocks;
    hostSeqInfo[6] = blockSize;
        
    hostParam[0] = param.match;
    hostParam[1] = param.mismatch;
    hostParam[2] = param.gapOpen;
    hostParam[3] = param.gapExtend;
    hostParam[4] = param.xdrop;
    hostParam[5] = param.marker;

    // Cuda Malloc
    // char* deviceSeqs;
    uint8_t* deviceFreq;
    int8_t* deviceAln;
    // int32_t* deviceIdx;
    int32_t* deviceLen;
    int32_t* deviceAlnLen;
    int32_t* deviceSeqInfo;
    int16_t* deviceParam;
    auto kernelStart = std::chrono::high_resolution_clock::now();
    // cudaMalloc((void**)&deviceSeqs, seqLen * seqNum * sizeof(char));
    cudaMalloc((void**)&deviceFreq, 12*seqLen * pairNum * sizeof(uint8_t));
    cudaMalloc((void**)&deviceAln, 2*seqLen * pairNum * sizeof(int8_t));
    // cudaMalloc((void**)&deviceIdx, 2*pairNum * sizeof(int32_t));
    cudaMalloc((void**)&deviceLen, 2*pairNum * sizeof(int32_t));
    cudaMalloc((void**)&deviceAlnLen, pairNum * sizeof(int32_t));
    cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
    cudaMalloc((void**)&deviceParam, 6 * sizeof(int16_t));

    // Copy to device
    // cudaMemcpy(deviceSeqs, hostSeqs, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceFreq, hostFreq, 12*seqLen * pairNum * sizeof(uint8_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceAln, hostAln, 2*seqLen * pairNum * sizeof(int8_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(deviceIdx, hostIdx, 2*pairNum * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceLen, hostLen, 2*pairNum * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceAlnLen, hostAlnLen, pairNum * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceParam, hostParam, 6 * sizeof(int16_t), cudaMemcpyHostToDevice);
    
    printf("Before kernel %s\n", cudaGetErrorString(cudaGetLastError()));
    alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
        // deviceSeqs, 
        deviceFreq,
        deviceAln, 
        // deviceIdx, 
        deviceLen,
        deviceAlnLen,
        deviceSeqInfo, 
        deviceParam
    );
    
    cudaDeviceSynchronize();
    printf("After kernel %s\n", cudaGetErrorString(cudaGetLastError()));

    cudaMemcpy(hostAln, deviceAln, 2*seqLen * pairNum * sizeof(int8_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(hostAlnLen, deviceAlnLen, pairNum * sizeof(int32_t), cudaMemcpyDeviceToHost);
    auto kernelEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    std::cout << "KernelTime "<< kernelTime.count() / 1000000<< " ms\n";

    // for (int k = 0; k < pairNum; ++k) {
    //     std::cout << hostAlnLen[k] << ','; 
    // }
    // std::cout << '\n';

    for (int k = 0; k < pairNum; ++k) {
        std::vector<std::string> alignment;
        int32_t refNum = seqIdx[k].second - seqIdx[k].first;
        int32_t qryNum = (k != pairNum-1) ? seqIdx[k+1].first - seqIdx[k].second : seqNum - seqIdx[k].second;
        int32_t refStart = seqIdx[k].first;
        int32_t qryStart = seqIdx[k].second;
        int32_t refIndex = 0;
        int32_t qryIndex = 0;
        // printf("k: %d, refNum: %d, qryNum: %d\n", k, refNum, qryNum);
        for (int j = 0; j < qryNum + refNum; ++j) alignment.push_back("");

        // for (int j = hostAlnLen[k] - 1; j >= 0; --j) {
        for (int j = 0; j < hostAlnLen[k]; ++j) {
            // if (nodes[k].first->identifier == "node_13") {
            //     // std::cout << (hostAln[k*2*seqLen+j] & 0xFFFF) << ',';
            //     if ((hostAln[k*2*seqLen+j] & 0xFFFF) != 0) 
            //         std::cout << "TB: " << (hostAln[k*2*seqLen+j] & 0xFFFF) << "at "<< j << '\n';
            // }
            if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 0) {
                for (size_t i=0; i<refNum; i++) alignment[i] += seqs[refStart+i][refIndex]; 
                for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += seqs[qryStart+i][qryIndex];
                qryIndex++;refIndex++;
            }
            else if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 2) {
                for (size_t i=0; i<refNum; i++) alignment[i] += seqs[refStart+i][refIndex];  
                for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += "-"; 
                refIndex++;
            }
            else {
                for (size_t i=0; i<refNum; i++) alignment[i] += "-"; 
                for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += seqs[qryStart+i][qryIndex]; 
                qryIndex++;
            }
        }
        // std::cout << '\n';
        // std::cout << nodes[k].first->identifier << ',' << nodes[k].second->identifier << '\n';
        // if (nodes[k].first->identifier == "node_13") {   
        // for (int i = 0; i < hostAlnLen[k]; ++i) {
        //     bool same = true;
        //     char firstchar = alignment[0][i];
        //     for (int j = 1; j < qryNum+refNum; ++j) {
        //         if ( alignment[j][i] != firstchar) same = false;
        //     }
        //     if (!same) {
        //         std::cout << "==========^= "<< i <<"  =================\n";
        //         for (int j = 0; j < qryNum+refNum; ++j) {
        //             std::cout << alignment[j].substr(max(0, i - 10), 20) << '\n';
        //         }
         
        //     }
        // }
        // }
        tree->allNodes[nodes[k].first->identifier]->msa = alignment;
    }

    // free device memory
    // cudaFree(deviceSeqs);
    cudaFree(deviceFreq);
    cudaFree(deviceAlnLen);
    cudaFree(deviceAln);
    // cudaFree(deviceIdx);
    cudaFree(deviceParam);
    cudaFree(deviceSeqInfo);
    cudaDeviceSynchronize();
    // free host memory
    // free(hostSeqs);
    free(hostFreq);
    // free(hostIdx);
    free(hostAlnLen);
    free(hostAln);
    free(hostParam);
    free(hostSeqInfo);

    return;
}
*/

void msaPostOrderTraversal_gpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{

    // assign msa to all nodes
    for (auto n: nodes) {
        // std::cout << n.first->identifier << '\n';
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msa.push_back(util->seqs[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
        // std::cout << n.second->identifier << '\n';
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msa.push_back(util->seqs[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
    }

    int numBlocks = 1024; 
    int blockSize = 512;

    
    // int alignSize = nodes.size() < numBlocks ? nodes.size() : numBlocks;
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        int32_t qryLen = tree->allNodes[n.second->identifier]->msa[0].size();
        int32_t refLen = tree->allNodes[n.first->identifier]->msa[0].size();
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    int round = nodes.size() / numBlocks + 1;
    for (int r = 0; r < round; ++r) {
        int alignSize = (nodes.size() - r*numBlocks) < numBlocks ? (nodes.size() - r*numBlocks) : numBlocks;
        if (alignSize == 0) break;
        // store all sequences to array
        int32_t seqNum = 0;
        int32_t pairNum = alignSize;
        std::vector<std::string> seqs;
        std::vector<std::vector<uint8_t>> freq;
        std::vector<std::pair<int32_t, int32_t>> seqIdx;
        std::vector<std::pair<int32_t, int32_t>> len;
        // store info to array 
        for (int n = 0; n < alignSize; ++n) {
            int32_t nIdx = n + r*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t qryLen = tree->allNodes[nodes[nIdx].second->identifier]->msa[0].size();
            int32_t refLen = tree->allNodes[nodes[nIdx].first->identifier]->msa[0].size();
            refIdx = seqNum;
            std::vector<uint8_t> temp;
            for (int i = 0; i < 12*seqLen; ++i) temp.push_back(0);
            
            assert(temp.size() == 12*seqLen);

            for (auto seq: tree->allNodes[nodes[nIdx].first->identifier]->msa) {
                for (int s = 0; s < refLen; ++s) {
                    if      (seq[s] == 'A' || seq[s] == 'a') temp[6*s+0]+=1;
                    else if (seq[s] == 'C' || seq[s] == 'c') temp[6*s+1]+=1;
                    else if (seq[s] == 'G' || seq[s] == 'g') temp[6*s+2]+=1;
                    else if (seq[s] == 'T' || seq[s] == 't') temp[6*s+3]+=1;
                    else if (seq[s] == 'N' || seq[s] == 'n') temp[6*s+4]+=1;
                    else                                     temp[6*s+5]+=1;
                }
                ++seqNum;
                seqs.push_back(seq);
            }
            qryIdx = seqNum;
            for (auto seq: tree->allNodes[nodes[nIdx].second->identifier]->msa) {
                for (int s = 0; s < qryLen; ++s) {
                    if      (seq[s] == 'A' || seq[s] == 'a') temp[6*(seqLen+s)+0]+=1;
                    else if (seq[s] == 'C' || seq[s] == 'c') temp[6*(seqLen+s)+1]+=1;
                    else if (seq[s] == 'G' || seq[s] == 'g') temp[6*(seqLen+s)+2]+=1;
                    else if (seq[s] == 'T' || seq[s] == 't') temp[6*(seqLen+s)+3]+=1;
                    else if (seq[s] == 'N' || seq[s] == 'n') temp[6*(seqLen+s)+4]+=1;
                    else                                     temp[6*(seqLen+s)+5]+=1;
                }
                ++seqNum;
                seqs.push_back(seq);
            }
            // printf("len: (%d, %d), num: (%d, %d)\n", refLen, qryLen, refIdx, qryIdx);
            seqIdx.push_back(std::make_pair(refIdx, qryIdx));
            len.push_back(std::make_pair(refLen, qryLen));
            freq.push_back(temp);
        }
        // Malloc
        uint8_t* hostFreq = (uint8_t*)malloc(12*seqLen * pairNum * sizeof(uint8_t));
        int8_t* hostAln = (int8_t*)malloc(2*seqLen * pairNum * sizeof(int8_t));
        int32_t* hostLen = (int32_t*)malloc(2*pairNum * sizeof(int32_t));
        int32_t* hostAlnLen = (int32_t*)malloc(pairNum * sizeof(int32_t));
        int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
        int16_t* hostParam = (int16_t*)malloc(6 * sizeof(int16_t)); 
        // Store Info to host mem
        for (int j = 0; j < 2*pairNum; ++j) { 
            if (j%2 == 0) hostLen[j] = len[j/2].first;
            else          hostLen[j] = len[j/2].second;
        }
        for (int j = 0; j < pairNum; ++j) {
            for (int l = 0; l < 12*seqLen; ++l) {
                hostFreq[12*seqLen*j+l] = freq[j][l];
            }
        }
        for (int j = 0; j < 2*seqLen*pairNum; ++j) { 
            hostAln[j] = 0;
        }
        for (int j = 0; j < pairNum; ++j) { 
            hostAlnLen[j] = 0;
        }
        hostSeqInfo[0] = seqLen;
        hostSeqInfo[1] = seqNum;
        hostSeqInfo[2] = pairNum;
        hostSeqInfo[3] = numBlocks;
        hostSeqInfo[4] = blockSize;
        hostSeqInfo[5] = numBlocks;
        hostSeqInfo[6] = blockSize;
        hostParam[0] = param.match;
        hostParam[1] = param.mismatch;
        hostParam[2] = param.gapOpen;
        hostParam[3] = param.gapExtend;
        hostParam[4] = param.xdrop;
        hostParam[5] = param.marker;

        // Cuda Malloc
        uint8_t* deviceFreq;
        int8_t* deviceAln;
        int32_t* deviceLen;
        int32_t* deviceAlnLen;
        int32_t* deviceSeqInfo;
        int16_t* deviceParam;
        auto kernelStart = std::chrono::high_resolution_clock::now();
        cudaMalloc((void**)&deviceFreq, 12*seqLen * pairNum * sizeof(uint8_t));
        cudaMalloc((void**)&deviceAln, 2*seqLen * pairNum * sizeof(int8_t));
        cudaMalloc((void**)&deviceLen, 2*pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceAlnLen, pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
        cudaMalloc((void**)&deviceParam, 6 * sizeof(int16_t));
        // Copy to device
        cudaMemcpy(deviceFreq, hostFreq, 12*seqLen * pairNum * sizeof(uint8_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceAln, hostAln, 2*seqLen * pairNum * sizeof(int8_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceLen, hostLen, 2*pairNum * sizeof(int32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceAlnLen, hostAlnLen, pairNum * sizeof(int32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(deviceParam, hostParam, 6 * sizeof(int16_t), cudaMemcpyHostToDevice);

        // printf("Before kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
            deviceFreq,
            deviceAln, 
            deviceLen,
            deviceAlnLen,
            deviceSeqInfo, 
            deviceParam
        );

        cudaDeviceSynchronize();
        // printf("After kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        // Copy to host
        cudaMemcpy(hostAln, deviceAln, 2*seqLen * pairNum * sizeof(int8_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(hostAlnLen, deviceAlnLen, pairNum * sizeof(int32_t), cudaMemcpyDeviceToHost);
        auto kernelEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
        if (round > 1) {
            printf("Round. %d align %d pairs. KernelTime: %d ms\n", r, alignSize, kernelTime.count() / 1000000);
        }
        else {
            std::cout << "KernelTime "<< kernelTime.count() / 1000000<< " ms\n";
        }

        for (int k = 0; k < pairNum; ++k) {
            std::vector<std::string> alignment;
            int32_t refNum = seqIdx[k].second - seqIdx[k].first;
            int32_t qryNum = (k != pairNum-1) ? seqIdx[k+1].first - seqIdx[k].second : seqNum - seqIdx[k].second;
            int32_t refStart = seqIdx[k].first;
            int32_t qryStart = seqIdx[k].second;
            int32_t refIndex = 0;
            int32_t qryIndex = 0;
            // printf("k: %d, refNum: %d, qryNum: %d\n", k, refNum, qryNum);
            // printf("k: %d, length: %d\n", k, hostAlnLen[k]);
            for (int j = 0; j < qryNum + refNum; ++j) alignment.push_back("");
            int nIdx = k + r*numBlocks;
            // printf("k: %d, length: %d, %s\n", k, hostAlnLen[k], nodes[nIdx].first->identifier.c_str());
            // for (int j = hostAlnLen[k] - 1; j >= 0; --j) {
            if (hostAlnLen[k] <= 0) {
                std::vector<std::string> reference, query;
                std::vector<int8_t> aln;
                for (auto s: tree->allNodes[nodes[nIdx].first->identifier]->msa) reference.push_back(s);
                for (auto s: tree->allNodes[nodes[nIdx].second->identifier]->msa) query.push_back(s);
                Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, param.xdrop, param.marker);
                Talco_xdrop::Align (
                    talco_params,
                    reference,
                    query,
                    aln
                );
                for (int j = 0; j < aln.size(); ++j) {
                    // std::cout << j << ',' << refIndex << ',' << qryIndex << '\n';
                    if ((aln[j] & 0xFFFF) == 0) {
                        for (size_t i=0; i<refNum; i++) alignment[i]        += reference[i][refIndex]; 
                        for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += query[i][qryIndex];
                        qryIndex++;refIndex++;
                    }
                    else if ((aln[j] & 0xFFFF) == 2) {
                        for (size_t i=0; i<refNum; i++) alignment[i]        += reference[i][refIndex]; 
                        for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += '-';
                            refIndex++;
                        }
                    else {
                        for (size_t i=0; i<refNum; i++) alignment[i]        += '-'; 
                        for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += query[i][qryIndex];
                        qryIndex++;
                    }
                }
                printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                // std::cout << "Start: " << tree->allNodes[nodes[nIdx].first->identifier]->identifier << '\n';
                // std::cout << "CPU fallback on No. " << k << ", Alignment Length: "<< aln.size() << '\n';
                // tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
                // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;
                
            }
            else {
                for (int j = 0; j < hostAlnLen[k]; ++j) {
                    if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 0) {
                        for (size_t i=0; i<refNum; i++) alignment[i] += seqs[refStart+i][refIndex]; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += seqs[qryStart+i][qryIndex];
                        qryIndex++;refIndex++;
                    }
                    else if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 2) {
                        for (size_t i=0; i<refNum; i++) alignment[i] += seqs[refStart+i][refIndex];  
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += "-"; 
                        refIndex++;
                    }
                    else {
                        for (size_t i=0; i<refNum; i++) alignment[i] += "-"; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)] += seqs[qryStart+i][qryIndex]; 
                        qryIndex++;
                    }
                }
                // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;
            }
            tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
            tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
            tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;
            // printf("(refNum, qryNum) = (%d, %d)\n", refNum, qryNum);
            // std::cout << "Finish: " << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ',' << tree->allNodes[nodes[nIdx].second->identifier]->identifier << '\n';
        }     
        // free memory
        cudaFree(deviceFreq);
        cudaFree(deviceAlnLen);
        cudaFree(deviceAln);
        cudaFree(deviceParam);
        cudaFree(deviceSeqInfo);
        cudaDeviceSynchronize();
        free(hostFreq);
        free(hostAlnLen);
        free(hostAln);
        free(hostParam);
        free(hostSeqInfo);
        
    }
    return;
}

void msaPostOrderTraversal_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{
    // assign msa to all nodes
    for (auto n: nodes) {
        // std::cout << n.first->identifier << '\n';
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msa.push_back(util->seqs[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
        // std::cout << n.second->identifier << '\n';
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msa.push_back(util->seqs[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
    }


    Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, param.xdrop, param.marker);

    auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        std::vector<std::string> reference, query;
        std::vector<int8_t> aln;
        for (auto s: tree->allNodes[n.first->identifier]->msa) reference.push_back(s);
        for (auto s: tree->allNodes[n.second->identifier]->msa) query.push_back(s);
        Talco_xdrop::Align (
            talco_params,
            reference,
            query,
            aln
        );
        
        int32_t refIndex = 0;
        int32_t qryIndex = 0;
        int32_t refNum = reference.size();
        int32_t qryNum = query.size();
        std::vector<std::string> alignment;
        for (int j = 0; j < refNum+qryNum; ++j) alignment.push_back("");
        for (int j = 0; j < aln.size(); ++j) {
            // std::cout << j << ',' << refIndex << ',' << qryIndex << '\n';
            if ((aln[j] & 0xFFFF) == 0) {
                for (size_t i=0; i<refNum; i++) alignment[i]        += reference[i][refIndex]; 
                for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += query[i][qryIndex];
                qryIndex++;refIndex++;
            }
            else if ((aln[j] & 0xFFFF) == 2) {
                for (size_t i=0; i<refNum; i++) alignment[i]        += reference[i][refIndex]; 
                for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += '-';
                refIndex++;
            }
            else {
                for (size_t i=0; i<refNum; i++) alignment[i]        += '-'; 
                for (size_t i=0; i<qryNum; i++) alignment[i+refNum] += query[i][qryIndex];
                qryIndex++;
            }
        }
        tree->allNodes[n.first->identifier]->refStartPos = refNum;
        tree->allNodes[n.first->identifier]->msa.clear();
        tree->allNodes[n.first->identifier]->msa = alignment;
        // std::cout << "Finish: " << tree->allNodes[n.first->identifier]->identifier << " RefEnd " << tree->allNodes[n.first->identifier]->refEndPos << '\n';
    }
    
    auto kernelEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    std::cout << "RunTime "<< kernelTime.count() / 1000000<< " ms\n";

    return;
}


void transitivityMerge_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util)
{
    // assign msa to all nodes
    /*
    for (auto n: nodes) {
        // std::cout << n.first->identifier << '\n';
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msa.push_back(util->seqs[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
        // std::cout << n.second->identifier << '\n';
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msa.push_back(util->seqs[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msa.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msa.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msa.push_back(util->seqs[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msa = node->children[childIndex]->msa;
                        break;
                    }
                }
            }
        }
    }    
    */
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        // std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        std::vector<bool> allGapsR, allGapsQ;
        std::vector<std::string> reference, query;
        std::vector<std::string> alignment;
        // std::vector<int8_t> aln;
        for (auto s: tree->allNodes[n.first->identifier]->msa) reference.push_back(s);
        for (auto s: tree->allNodes[n.second->identifier]->msa) query.push_back(s);
        int32_t refLen = reference[0].size();
        int32_t qryLen = query[0].size();
        int32_t refNum = reference.size();
        int32_t qryNum = query.size();
        int32_t seqLen = max(refLen, qryLen);
        int32_t refStart = tree->allNodes[n.first->identifier]->refStartPos;
        int32_t qryStart = tree->allNodes[n.second->identifier]->refStartPos;
        int32_t parentNumRef = refNum - refStart;
        int32_t parentNumQry = qryNum - qryStart;
        // std::cout << refNum << ',' << qryNum << ',' << refLen << ',' << qryLen << '\n';
        // std::cout << refStart << ',' << qryStart << ',' << parentNumRef << ',' << parentNumQry << '\n';
        if ((parentNumRef == qryNum || parentNumQry == refNum) && tree->allNodes[n.first->identifier]->parent != nullptr) continue; 
        assert(parentNumRef == parentNumQry);

        for (int i = 0; i < refNum + qryStart; ++i) alignment.push_back("");
        // for (int i = 0; i < refLen; ++i) allGapsR.push_back(true);
        // for (int i = 0; i < qryLen; ++i) allGapsQ.push_back(true);
        // for (int j = refStart; j < refNum; ++j) {
        //     for (int r = 0; r < refLen; ++r) {
        //         if (!allGapsR[r]) continue;
        //         if (reference[j][r] != '-') allGapsR[r] = false;
        //     }
        // }
        // for (int j = qryStart; j < qryNum; ++j) {
        //     for (int q = 0; q < qryLen; ++q) {
        //         if (!allGapsQ[q]) continue;
        //         if (query[j][q] != '-') allGapsQ[q] = false;
        //     }
        // }

        for (int i = 0; i < seqLen; ++i) {
            if (i < refLen) {
                bool allGaps = true;
                for (int j = refStart; j < refNum; ++j) {
                    if (reference[j][i] != '-') {
                        allGaps = false;
                        break;
                    }
                }
                allGapsR.push_back(allGaps);
            }
            if (i < qryLen) {
                bool allGaps = true;
                for (int j = qryStart; j < qryNum; ++j) {
                    if (query[j][i] != '-') {
                        allGaps = false;
                        break;
                    }
                }
                allGapsQ.push_back(allGaps);
            }
        }
        
        int32_t rIdx = 0, qIdx = 0;
        assert(allGapsR.size() == refLen);
        assert(allGapsQ.size() == qryLen);
        while (rIdx < refLen && qIdx < qryLen) {
            // if (tree->allNodes[n.first->identifier]->identifier == "node_102") printf("(%d,%d)\n", rIdx, qIdx);
            // if (rIdx > 6400 && rIdx < 6700) std::cout << allGapsR[rIdx] << ',' << allGapsQ[qIdx] << '\n';
            if (allGapsR[rIdx] == false && allGapsQ[qIdx] == false) {
                for (size_t i=0; i<refStart; i++)      alignment[i]          += reference[i][rIdx]; 
                for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += query[i][qIdx];
                for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += reference[i][rIdx];
                qIdx++;rIdx++;
            }
            else if (allGapsR[rIdx] == true && allGapsQ[qIdx] == false) {
                int consecGap = 0;
                int k = rIdx;
                while (allGapsR[k] && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refStart; i++)      alignment[i]          += reference[i][rIdx]; 
                    for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += '-';
                    for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
                    rIdx += 1;
                }
                // if (k == refLen - 1) break;
            }
            else if (allGapsR[rIdx] == false && allGapsQ[qIdx] == true) {
                
                int consecGap = 0;
                int k = qIdx;
                while (allGapsQ[k] && k < qryLen) {
                    ++consecGap;
                    ++k;
                }
                // std::cout << "Q:" << qIdx << "consecGap: " << consecGap << '\n';
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refStart; i++)      alignment[i]          += '-'; 
                    for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += query[i][qIdx];
                    for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
                    qIdx += 1;
                }
                // if (k == qryLen - 1) break;
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (allGapsR[rIdx] && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refStart; i++)      alignment[i]          += reference[i][rIdx]; 
                    for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += '-';
                    for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
                    rIdx += 1;
                }
                consecGap = 0;
                while (allGapsQ[qIdx] && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refStart; i++)      alignment[i]          += '-'; 
                    for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += query[i][qIdx];
                    for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
                    qIdx += 1;
                }
            }
        }
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                for (size_t i=0; i<refStart; i++)      alignment[i]          += reference[i][g]; 
                for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += '-';
                for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                for (size_t i=0; i<refStart; i++)      alignment[i]          += '-'; 
                for (size_t i=0; i<qryStart; i++)      alignment[i+refStart] += query[i][g];;
                for (size_t i=refStart; i<refNum; i++) alignment[i+qryStart] += '-';
            }
        }
        
        tree->allNodes[n.first->identifier]->refStartPos = refStart+qryStart;
        tree->allNodes[n.first->identifier]->msa.clear();
        tree->allNodes[n.first->identifier]->msa = alignment;
        // std::cout << "Check (Length, SeqNum) = ("<< alignment[0].size() << ", " << alignment.size() << ')' << '\n';
    }
    
    // auto kernelEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    // std::cout << "RunTime "<< kernelTime.count() / 1000<< " us\n";

    return;
}

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& hier, std::stack<Node*> msaStack, int grpID, int mode) {
    int hierIdx = 0;
    // mode 0: msa, 1: merge
    while(!msaStack.empty()) {
        Node* node = msaStack.top();
        if (!(node->grpID==-1 || node->grpID==grpID)) {
            msaStack.pop();
            continue;
        };
        if (node->children.size()==0) {
            msaStack.pop();
            continue;
        }
        size_t childIndex = 0;
        if (mode == 0) {
            for (childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID))
                {
                    break;
                }
            }
            // std::cout << node->identifier << '\n';
            if (childIndex == node->children.size() - 1) {
                msaStack.pop();
                continue;
            }
            for (size_t i=childIndex+1; i<node->children.size(); i++)
            {
                if (!(node->children[i]->grpID == -1 || node->children[i]->grpID == grpID))
                {
                    continue;
                }
                hier.push_back(std::make_pair(std::make_pair(node, node->children[i]),hierIdx));
                ++hierIdx;
            }
            msaStack.pop();
        }
        else {
            for (size_t i=0; i<node->children.size(); i++) {
                if (!(node->children[i]->grpID == -1 || node->children[i]->grpID == grpID))
                {
                    continue;
                }
                hier.push_back(std::make_pair(std::make_pair(node, node->children[i]),hierIdx));
                ++hierIdx;
            }
            msaStack.pop();
        }
    }
    // if (mode == 1) for (auto h: hier) {
    //     std::cout << h.first.first->identifier << ',' << h.first.second->identifier << ',' << h.first.first->level << ',' << h.second << '\n';
    // }
    hierIdx = 0;
    std::stack<std::pair<Node*, int>> hierStack; 
    // Node* tempRoot = hier[0].first.first->parent;
    Node* preNode = hier[0].first.first;
    size_t prelevel = hier[0].first.first->level;
    hier[0].second = hierIdx;
    for (int k = 1; k < hier.size(); ++k) {
        if (!hierStack.empty()) {
            if (hier[k].first.first->identifier == hierStack.top().first->identifier) {
                hierIdx = max(hierIdx+1, hierStack.top().second);
                hier[k].second = hierIdx; 
                prelevel = hier[k].first.first->level;
                hierStack.pop();
            }
            else {
                if (mode == 0) {
                    if (hier[k].first.first->level <= prelevel) {
                        hier[k].second = ++hierIdx;
                        prelevel = hier[k].first.first->level;
                    }
                    else {
                        hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                        hier[k].second = 0;
                        hierIdx = 0;
                        prelevel = hier[k].first.first->level;
                    }
                }
                else {
                    if (hier[k].first.first->level < prelevel || (hier[k].first.first->level == prelevel && hier[k].first.first->identifier == preNode->identifier)) {
                        hier[k].second = ++hierIdx;
                        prelevel = hier[k].first.first->level;
                    }
                    else {
                        if (preNode->parent->identifier == hierStack.top().first->identifier) {
                            hierStack.top().second = max(hierIdx+1, hierStack.top().second);
                        }
                        else {
                            hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                        }
                        
                        hier[k].second = 0;
                        hierIdx = 0;
                        prelevel = hier[k].first.first->level;
                    }
                }
            }
        }
        else {
            if (mode == 0) {
                if (hier[k].first.first->level <= prelevel) {
                    hier[k].second = ++hierIdx;
                    prelevel = hier[k].first.first->level;
                }
                else {
                    hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                    hier[k].second = 0;
                    hierIdx = 0;
                    prelevel = hier[k].first.first->level;
                }
            }
            else {
                if (hier[k].first.first->level < prelevel || (hier[k].first.first->level == prelevel && hier[k].first.first->identifier == preNode->identifier)) {
                    hier[k].second = ++hierIdx;
                    prelevel = hier[k].first.first->level;
                }
                else {
                    hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                    hier[k].second = 0;
                    hierIdx = 0;
                    prelevel = hier[k].first.first->level;
                }
            }
        }
        preNode = hier[k].first.first;
    }
    // for (int k = 0; k < hier.size(); ++k) {
    //     std::cout << hier[k].first.first->identifier << ',' << hier[k].first.second->identifier << ',' << hier[k].first.first->level << ',' << hier[k].second << '\n';
    // }
}


void getPostOrderList(Node* node, std::stack<Node*>& msaStack) {
    std::stack<Node*> s1;
    std::vector<std::vector<Node*>> hier;
    size_t lowestLevel = 0; 
    
    s1.push(node); 
    Node* current; 
  
    while (!s1.empty()) { 
        current = s1.top(); 
        msaStack.push(current); 
        if (current->level > lowestLevel) lowestLevel = current->level;
        s1.pop(); 
        for (auto ch: current->children) {
            if (ch->grpID == current->grpID) {
                s1.push(ch);
            }      
        }
    } 
    return;
}


int main(int argc, char** argv) {

    auto mainStart = std::chrono::high_resolution_clock::now();
    std::string seqFileName;
    

    parseArguments(argc, argv);

    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << mainDesc << std::endl;
        // Return with error code 1 unless the user specifies help
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    // Define MSA utility
    msa::utility* util = new msa::utility;

    // Read Tree (Newick String)
    // std::cout << "Start read tree...\n";
    Tree* T = readNewick(vm);
    
    printTree(T->root);

    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    // paritionInfo_t * partition = new paritionInfo_t(2,0, "longest"); /*Starting with zero partition*/
    
    int maxSubtreeSize = vm["max_leaves"].as<int>();
    paritionInfo_t * partition = new paritionInfo_t(maxSubtreeSize, 0, "centroid");
    // std::cout << "Start Partition ..... \n";
    partitionTree(T->root, partition);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Partition the tree in: " <<  treeBuiltTime.count() << " ns\n";
    // printLeaves(T->root);

    // Read Input Sequences (Fasta format)
    readSequences(vm, util);

    // Test alignSeqToSeq
    // Node* root = T->root;
    int marker = 128;
    int xdrop = 200;
    Params param(2,-1,-2,-1,xdrop,marker);

    for (auto p: partition->partitionsRoot) {
        std::cout << p.first << std::setw(5) << p.second.second << '\n';
    }

    Tree* newT = reconsturctTree(T->root, partition->partitionsRoot);
    printTree(newT->root);
    // for (auto n: newT->allNodes) {
    //     std::cout << n.first << ',' << n.second->identifier << "\nChildren\n";
    //     for (auto ch: n.second->children) {
    //         std::cout << ch->identifier << ',';
    //     }
    //     std::cout << '\n';
    // }

    
    // printTree(T->root);
    // size_t maxlen = 0;
    // for (auto s: util->seqs) {
    //     if (s.second.size() > maxlen) maxlen = s.second.size();
    // }


    
    std::cout << "Start MSA...\n";
    auto msaStart = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    for (auto &p: partition->partitionsRoot)
    {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
        int grpID = p.second.first->grpID;
        getMsaHierachy(subhier, msaStack, grpID, 0);
        for (auto h: subhier) {
            while (hier.size() < h.second+1) {
                std::vector<std::pair<Node*, Node*>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
    }
    std::string machine = vm["machine"].as<std::string>();
    int level = 0;
    for (auto m: hier) {
        std::cout << "Aln level: " << level << '\n';
        auto alnStart = std::chrono::high_resolution_clock::now();
        if (machine == "cpu" || machine == "CPU" || machine == "Cpu") {
            msaPostOrderTraversal_cpu(T, m, util, param);
        }
        else if (machine == "gpu" || machine == "GPU" || machine == "Gpu") {
            msaPostOrderTraversal_gpu(T, m, util, param);
        }
        else {
            fprintf(stderr, "Error: Unrecognized machine type: %s\n", machine.c_str()); 
            exit(1);
        }
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (m.size() > 1) std::cout << "Aln "<<  m.size() <<" pairs in " <<  alnTime.count() / 1000000 << " ms\n";
        else              std::cout << "Aln "<<  m.size() <<" pair in " <<  alnTime.count() / 1000000 << " ms\n";
        ++level;
    }
    
    auto msaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds msaTime = msaEnd - msaStart;
    std::cout << "MSA in " <<  msaTime.count() / 1000000000 << " s\n";

    // Push MSA to roots 
    for (auto p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<Node*> msaArray;
        while (!msaStack.empty()) {
            msaArray.push_back(msaStack.top());
            msaStack.pop();
        }
        // for (auto m: msaArray) std::cout << m->identifier << '(' << m->msa.size() << "),";
        // std::cout << '\n';
        if (msaArray.back()->msa.size() == 0 && msaArray.size() > 1) {
            if (msaArray.size() == 2) {
                T->allNodes[msaArray.back()->identifier]->msa = msaArray[0]->msa;
                break;
            }
            for (int m = msaArray.size()-2; m >=0; --m) {
                if (msaArray[m]->msa.size()>0) {
                    T->allNodes[msaArray.back()->identifier]->msa = msaArray[m]->msa;
                    break;
                }
            }
        }
        // for (auto m: msaArray) std::cout << m->identifier << '(' << m->msa.size() << "),";
        // std::cout << '\n';
    }
    if (partition->partitionsRoot.size() == 1) {
        auto mainEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mainTime = mainEnd - mainStart;
        std::cout << "Total Execution in " <<  mainTime.count() / 1000000000 << " s\n";
        return;
    }
    // type 1 alignment
    // for (auto n: newT->allNodes) {
    //     std::cout << n.first << ':' << T->allNodes[n.first]->msa.size() << '\n';
    //     // std::cout << '\n';
    // }
    std::cout << "Start MSA between parent and children\n";
    std::vector<std::pair<Node*, Node*>> type1Aln;
    for (auto n: newT->allNodes) {
        for (auto m: n.second->children) {
            type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
        }
    }
    // for (auto n: type1Aln) {
    //     std::cout << n.first->identifier << ',' << n.second->identifier << '\n';
    // }
    auto alnStart = std::chrono::high_resolution_clock::now();
    if (machine == "cpu" || machine == "CPU" || machine == "Cpu") {
        msaPostOrderTraversal_cpu(T, type1Aln, util, param);
    }
    else if (machine == "gpu" || machine == "GPU" || machine == "Gpu") {
        msaPostOrderTraversal_gpu(T, type1Aln, util, param);
    }
    else {
        fprintf(stderr, "Error: Unrecognized machine type: %s\n", machine.c_str()); 
        exit(1);
    }
    auto alnEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    if (type1Aln.size() > 1) std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pairs in " <<  alnTime.count() / 1000000 << " ms\n";
    else                     std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pair in "  <<  alnTime.count() / 1000000 << " ms\n";
    
    

    paritionInfo_t * newPartition = new paritionInfo_t(std::numeric_limits<size_t>::max(), 0, "centroid");
    // std::cout << "Start Partition ..... \n";
    partitionTree(newT->root, newPartition);
    hier.clear();
    for (auto &p: newPartition->partitionsRoot)
    {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
        int grpID = p.second.first->grpID;
        getMsaHierachy(subhier, msaStack, grpID, 1);
        while (!msaStack.empty()) {
            // std::cout << msaStack.top()->identifier << ',' << msaStack.top()->level << ',';
            msaStack.pop();
        }
        // std::cout << '\n';
        for (auto h: subhier) {
            while (hier.size() < h.second+1) {
                std::vector<std::pair<Node*, Node*>> temp;
                hier.push_back(temp);
            }
            hier[h.second].push_back(h.first);
        }
    }
    level = 0;
    auto mergeStart = std::chrono::high_resolution_clock::now();
    for (auto m: hier) {
        std::cout << "Transtivity merge level: " << level << '\n';
        auto submergeStart = std::chrono::high_resolution_clock::now();
        transitivityMerge_cpu(T, m, util);
        // if (machine == "cpu" || machine == "CPU" || machine == "Cpu") {
        //     msaPostOrderTraversal_cpu(T, m, util, param);
        // }
        // else if (machine == "gpu" || machine == "GPU" || machine == "Gpu") {
        //     msaPostOrderTraversal_gpu(T, m, util, param);
        // }
        // else {
        //     fprintf(stderr, "Error: Unrecognized machine type: %s\n", machine.c_str()); 
        //     exit(1);
        // }
        auto submergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds submergeTime = submergeEnd - submergeStart;
        if (m.size() > 1) std::cout << "Merge "<<  m.size() <<" subtree pairs in " << submergeTime.count() / 1000000 << " ms\n";
        else              std::cout << "Merge "<<  m.size() <<" subtree pair in " <<  submergeTime.count() / 1000000 << " ms\n";
        ++level;
    }
    auto mergeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
    std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << T->root->msa.size() << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    // for (auto &p: partition->partitionsRoot)
    // {
    //     if (p.second.first->partitionParent == nullptr) continue;
    //     std::pair<std::vector<std::string>, std::vector<std::string>> alignments;
    //     std::vector<std::string> ref = p.second.first->partitionParent->msa;
    //     std::vector<std::string> query = p.second.first->msa;
    //     alignGrpToGrp(ref, query, param, alignments);
    //     ref.clear();
    //     for (auto &s: alignments.first) ref.push_back(s);
    //     for (auto &s: alignments.second) ref.push_back(s);
    //     p.second.first->msa = ref;
    // }
    
    // for (auto &p: partition->partitionsRoot)
    // {
    //     std::cout << p.first << std::endl;
    //     for (auto& s: p.second.first->msa)
    //     {
    //         std::cout << s << "\n";
    //     }
    // }
    
    


    // std::vector<std::string> refs, querys;
    // std::pair<std::vector<std::string>,std::vector<std::string>> alignments;
    // refs.push_back("-ATGACCAGCTTGAAACGGTCGCAGACAGAAAGGCCCGT-TGCCACCGAAAGGGCCTCCGTTGTGGGAACAGATGGCACTC----CCAAAGTGCACTCTGATGACTTCTACATGC--GAC-GCTTCAGGTCCCAGAA-TGGCAGTTTAGGATCCTCAGTCATGG--CTCCTGTGGGGCCCCCTCGAAGTGAAGGCCCTCACCATATAACCTCAACCCCCGGAGTCCCAAAGATGGGGGTTAGGGCA-------AGAATTGCAGATTGGCCCCCAAGA-AAGGACAACGTGAAA-GAATCT---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCATGTCCTCCAAAAGCAGTCCTGTGAGTCAGGGAAGTTCTGTTAGCCTCAATTCCAATGACTCAGCTATGTTAAAAAGCATACAGAACAC--GCTGAAAAACAAGACAAGACCGTCGGAGAACATGGACTCCAGATTTCTCATGCCTGAAGCCTATCCCAGCTCCCCTAGGAAAGCTCTTCGCAGAATACGGCAGCGGAGC-AACAGCGACATCACCATAAGTGAACTTGATGTGGATAGCTTTGATGAATGTATCTCACCCACATATAAA-ACTGGACCATCACTGCACAGGGAATACGGTAGCACATCTTCAATTGATAAGC-AGGGAACGTCCGGAGAAGGCTTTTTTGACTTGTTAAAGGGCTACAAAGATGACAAACCTGACC---GAGGTCCAACCCCAACCAAGCTCAGTGACTTTCTTATCGCCGGTGGGGGCAAGGGGTCTGGTTTCTC--CTTGGATGTTATAGA---TGGTGCCATTTCACAAAGGGAGAACCTCAGGCTGTTTAAGGAGAGGGAAAAACCAC-TCAAGCGACGTTCCAAGTCTGAAACTGGAGACTCATCTATTTTTCGTAAACTACGCAATGCCAAAGGTGAA---GAAC---TTGGG---AAA---TCATCGGATCTTGAAGAT---CGATCAGA------TTCTGTC---CCTTGGACGTGCCCTAAGTGC--TTTGCTCACTACGATGTCCAGAGTANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNN-NNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAACACTGGAACTATTTTGGGGCTGATGAAAATCT-TGGTCCGGTGGCTGTGAGCATTCGAAGGGAAAAACCAGAAGAAATGAAAGA-AAATGGATCTCCATACAACTACCGGATAATTTTTAGAACTAGTGAGCTCATGACACTG-AGAGGGTCAGTCCTGGAGGATGCCA-TCCCCTCAACAGCAAAGCACTCA-ACAGCCAGGGGGC-TTCCCCTGAAAGAGGTGCTGGAGCATGTGATTCCCGAGCT-------CAATGTCCAGTGCCTGCGGTTGGC--------------------------------------------------CTTCAACA--------------------------CCCCCAAAGTCA---------------CAGAACA-------------------G-CTGATGAAACTGGATGAACAAGG----------------------------------------------------------------------------------------------------------------------GCTGAACTATCAGCAGAAAGTAGGCATCATGTACTGCAAAGCTGGACAGAGCACAGAAGAAGAGATGTACAACAATGAATCAGCCAGCCCAGCCTTTGAGGAATTCCTTCAGCTGTTGGGAGAACGAGTTCGGCTCAAGGGATTTGAGAAGTACCGTGCACAGCTTGACAACAAAACTGACTCCACTGGAACTCATTCTCTGTATACGACATACAAAGACTATGAAATCATGTTCCATGTTTCTACCATG-CTGCCCTACACACCTAACAACAAGCAACAG---------------------------NNNNN-NNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNN------NNNNNNNNNNNNNN---NNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCAGTCACCAGGTCCAGAGATGTGCCTTCC---TTTGGGCCTCCCATCCCC-AAAGGAGTCACT-TTCCCTAAATCGAATGTGTTCAGAGACTTCCTTTTGGCCAA--AGTGATTAATGCAGAGAATGCTGCTCA-TAAATC-AGAAAAGTTCCGGGCCA--TGGCTACTCGGAC-CCGCCAG---GAATACCTGAAAGATCTGGCAG-AAAAGAATGTCACCAACACCCCTATTGATCCTTCTGGCAAGTTCCCATTCATCTCTCTGGCCTC-CAAGAAGAAGGAGAAGTCCAAGCCATATCC--AGGAGCTGAGCTCAGTAGC-ATGGGTGCCATTGTGT--GGGCAATTCGGGCCAAAGA-CTACAACCAGGCTATGGA-AATCGACTGTCTTTTGGGGATCTCCAATGAGTTCATT-GTCCTCAT-TGAACAGGAAACAAAGAGCGTGGTTTTCAATTGCTCCTGCAG-AGATGTGATAGGGTGGACTTCAACTGACAACAGCC-TCAAAATTTTCTATGAAC-GAGGAGAATGTGTTTCTGTGGAGAGTTTCATAAACAA---TGAGGATATCAAAGAGAT-TGTCAAAAGGTTACAGNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNCTTGGCTTCCATGTCAACTAT-GAGGGCATCGTGGCTG-ATGTAGAGCCCTATGGCTATGCCTGGC-AGGCAGGCCTAAGGCAGGGCAGCCGGCTGGTAGAGATCTGCAAAGTG--GCAGTGGCCACTCTGA-GCCATGAGCAGATGATC-GATCTCCTGAGAACATCCGTCACAGTGAAGGTGGTCATTATCCCCCCC--CATGATGACTGCACCCCACGGAGN------NNNNNN-NNNNNNNN--NNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-----------------------NNNNNNNNNNNNNNN-NNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAC--TGTCTCCTGGTTCAGACATCTATGTGACGGTCTCATCCATCGCTTTAGCAAGATCCCAGCAGTGCCGTAACTCCCCTAGCAACCTATCTTCATCCAGCGAGACTGGCTCTGGTGGGGGCACTTACAGACAGAAGTCCATGCCAGAA-GGGTTTGGAGTGAGCCGCAGATCCCCAGCCTCCATTGACAGGCAGAATACCCAGTCAGATAT---CGG-CGGCAGCGGAAAATCCACACCCAGCTGGCAAAGGAGT-GAGGACAGCATTGCTGACCAGATGGCTTACAGTTATAGAGGACCTC--AGGATTTCAATTCT---TTTGTCCTCGAGCAGCATGAATAT-ACAGAGCCAACGTGCCACCTCCCAGCAGTATCAAAGGTACTGCCAGCTTTCCGAGAGAGCCCCAGTGGGAGATTGATGCGGCAGGATCCAGTGGTTCATTTGTCTCCAAACAAGCAAGGGCATTCTGA-TAGCCACTACTCGAGCCACTCCAGT--AGCAATACCCTCTCCAGCAATGCATCGAGCGCC---CACAGTGACGAGAAGTGGTATGA---TGGGGACCGCACTGAATCCGAACTCAGCAGC-TATAACTATCTG-CAGGGCACGTCTGCGGACAGCGGCATCGAC-ACCACCTCCTA-CGGCCCCGGCCAC-GGCAGCACAGCCT-CCCTGGGGGCTG-CCACATCATCGCCTCGCT-CA---GGGCCTG-GCAAGGAGAAGGTGGCACCCCTGTGGCACAGTTCCAGTGAAGTCATC---TCCCTGGCAGACCGGACTTTAGAGGCGGAGAGCCACGGCATGGACCGGAAAGCCGAGTCTTCACTGAGCTTGGATATCCACAGCAAGAGCCAAGGAGGCTCCGGCCCTCTGACGAGGGAGAACAGCACTTTCAGTGTCAGTGATACTGC-CT-CCCACACAAGTACCATGAGCTCCCGACACTCG-GCCAGCCCGGT-CGTTTTCACCAGTG-CCAGGAGTTCC-CCTAAAGAAGAGCTTCATCCTGCCACCTCTTCGCAGCTCGTACCTTCATT------CTCTTCTTCCTCCTCCTCCTCGTCTGGTCCCAGGACTTTCTACCCTCGCCAGGGTGCT-ACGAGCAAGTACCTGATTGGGTGGAAAAAACCAGAAGGAACCATAAATTCTGTGGGATTTATGGACACAAGAAAGCGTCATCAGAGTGATGGCAATGAAATAACCCACACCAGGCTACAGGCCTCAACCAGGGACCTCCGAGCATCTCCTAAGCCAGCTTCCAAGTCCACCATTG--AAGAGGATCTAAAGAAACTGATAGACCTTGAAAGCCCAACTCCTGAATCACAGAAAAATTTT-AAG------TTCCACGCACTCTCCTCTCCTCAGTCTCCTTTCCCTCCCACCCCCACCTCAAGGCGGGC-CTTGCACAGAACGCTGTCGGATGAGAGTATTTACAGTGGCCAGAGAGAGCACTTCTTCACCTC--TAGG-GCCTCACTTCTGGACCAAGCC---CTGCCCAACGATGTCCTCTTCAGCAGCACGTACCCTTCTCTCCCCAAGTCTCTTCCGCTGCGGAGACCTTCTTAC-ACCTTAGGGATGAAGTCATTGCAT-GGAGAGTTCTCGGCTTCAGACAGCTCCCTTACTGACATCCAGGAGACCCGA---AGGCAGCCCATGCCCGACCCCGGCCTGATGCCCTTGCCTGATGCTGCT-GCAGA-TTTGGATTGGTCCAACTTGGTAGATGCTGCCAAAGCCTACGAG------GTCCAAAGAGCGTCGTTCTTTGCTGCTAGTGATGAAAATCATCGCCCCTTGAGTGCTGCATCCAACAGCGATCAG-CTTGAGGACCAGGCTCTGGCCCAGATGAAGTCTTACAGCAGCAGTAAAGATTCCTCTCCC---ACTCTGGCTTCCAAAGTGGACCAGCTGGAAGGCATGTTGAAGATGCTTCGGGAAGA-TTTGAAG---AAG-------------GAAAAAGAAGATAAAGCCCACCTTCAGGCCGA--AGTTCAGCACTTGCGTGAGGACAACCTGAGGCTACAAGAGGAGTCCCAGAATGCCTCGGACAAGCTGAAGAAGTTCACAGAGTGGGTCTTCAATACCATAGACATGAGCTA");
    // refs.push_back("--------------------------------------------------------------------------------------------------------------------------------GTC-------------------------------------------------------------------------------------------------AAGATGGGGGTTAGGGAAGCGA-ATAGAAGAGAGGATTGGCCCCCAAGA-AAGGAAAATGTGAAG-GAATCCATGAACGAGTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCCTGTCCTCCAAGAGCAGTCCTGTGAGTCAGGGAAGTTCTGTTAGCCTCAATTCCAATGACTCAGCTATGCTGAAAAGCATACAGAACAC--ACTGAAGAACAAGACAAGACCTTCGGAGAACATGGACTCCAGATTTCTCATGCCGGAAGCCTATCCCAGCTCCCCTAGAAAAGCTCTCCGTAGAATACGGCAGCGAAGC-AACAGTGATATCACCATAAGTGAACTTGATGTGGATAGCTTTGATGAATGTATCTCACCTACATATAAG-ACGGGACCATCACTGCACAGGGAATATGGTAGCACATCTTCAATTGATAAAC-AGGGAACATCTGGAGAAAACTTTTTTGACTTGTTAAAGGGCTACAAAGATGACAAATCTGACC---GAGGTCCGACCCCAACGAAGCTCAGTGACTTTCTTATTGCTGG------CAAAGGGTCTGGTTTCTC--CTTGGATGTTTTGGATGGTGGGTCCATCTCACAAAGGGAGAATCTCAGGCTCTTTAAGGAAAGGGAAAAACCAC-TCAAGCGGCGTTCCAAGTCTGAAACAGGAGACTCATCTATTTTTCGTAAACTACGCAATGCCAAAGGTGAA---GAAC---TTGGA---AAG---TCATCAGATCTTGAAGATAACCGGTCAGAAGA---TTCTGTCAGGCCTTGGACATGCCCTAAGTGT--TTTGCCCACTATGATGTCCAGAGTATATTATTTGACTTGAACGAGGCAATTATGAATAGGC---ACAATGTTATTAAGAGGAGAAA-CACCACCACAGGAGCTTCGGCTGCAGCCGTGGCATCCTT-GGTCTCTGGACCTTTGTCCCATTCAGCCAGTTTTAGCTCTCCCAT-GGGCAGTACGGAGGACCTGAATTCCAAGGGAAG------CCTTGGTATGGACCAGGGAGATGATAAAAGCAATGATCTTGT-AATGAGCTGTCCATATTTTCGGAATGAGATTGGTGGAGAAGGTGAAAGGAAGATAAGCCTTTCGAAATCAAAT--TCTGGTTCTTTTAGTGGATATGAAAGTACCTCCTTTGAGTCTACCCTTACCTCCCATTGCACCAACGCT-GGAGTGGCAGTTCTTGAAGTACCCAAGGAGAACTTGGTGTTGC-ATCTAGACA-GAGT-AAAAAGATACATTGTGGAACATGTGGATCTTGGTGCATACTACTATAGGAAATTTTTCTATCAGAAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTCATGANNNNN-NNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------NNNNNNNNNNNNNNNNNNNNNNNN--------------------------------------------------NNNNNNNN--------------------------NNNNNNNNNNNN---------------NNNNNNN-------------------N-NNNNNNNNNNNNNNNNNNNNNNN----------------------------------------------------------------------------------------------------------------------NCTGAACTATCAGCAGAAAGTAGGCATCATGTACTGCAAAGCTGGACAGAGCACCGAGGAAGAAATGTACAACAATGAGTCAGCCAGCCCAGCCTTTGAGGAGTTCCTTCAGCTGTTGGGAGAACGCATTAGGCTCAAAGGATTTGAGAAATATCGTGCACAGCTTGACACCAAGACTGACTCCACTGGAACTCATTCTCTGTACACAACATATAAAGACTATGAAATCATGTTTCATGTTTCTACCATG-CTGCCCTACACACCTAACAACAAGCAACAG---------------------------CTCCT-CCGCAAGCGGCAC--ATTGGAAATGATATCGTAACAATAGTTTTCCAAGAACCTGGAGCACAGCCATTCAGCC--CAAAAAACATTCGTTC------CCACTTTCAGCACG---TTTTCG--TCATTGTCAGGGCTCACAACCCCTGCACGGACAGTGTCTGTTACAGTGTGGCAGTCACCAGATCCAGAGATGTGCCTTCC---T------------TCCCC-AAAGGAGTCACT-TTCCTTAAATCAAATGTGTTCAGGGACTTCTCTTGGCCCAA--AGTGATCAATGCAGAAAATGCTGCTCA-TAAGTC-AGAAAAGTTTCGGGCCA--TGGCTACTCGGAC-CCGCCAG---GAATAC---AAAGATCTGGCCG-AAAAGAATGTCACCAACACTCCTGTTGATCCTTCTGGCAAGTTCCCGTTCATCTCTCTGGCCTC-CAAGAAGAAAGAGAAGTCCAAGCCGTATCC--AGGAGCTGAGCTCAGTAGC-ATGGGTGCCATTGTGT--GGGCAATCCGGGCCAAAGA-CTACAACCAGGCTGTGGA-GATCGACTGTCTTCTAGGGATCTCCAATGAGTTCGTC-GTCCTCAT-CGAACAGGAGACAAAGAGCGTGGTTTTCAACTGCTCCTGCAG-AGATGTGATAGGGTGGACTTCGACTGACACCAGCC-TCAAGATCTTCTATGAAC-GAGGAGAATGCGTTTCTGTGGAGAGTTTCATAAACAA---TGAGGATATCAAAGAGAT-TGTCAGAAGGTTACAGTTTGTGTCAAAAGGTTGTGAA-TCAGTGGAGATGACTCTGCGAAGGAATGGGCTA-GGACAGCTTGGCTTCCATGTCAACTAC-GAGGGCATCGTGGCTG-ATGTAGAGCCCTATGGCTATGCCTGGC-AGGCAGGCCTGAGGCAGGGCAGCCGGTTGGTAGAGATCTGCAAAGTA--GCAGTGGCCACTCTGA-GCCATGAGCAAATGATC-GATCTCCTCAGAACATCTGTCACGGTGAAGGTGGTCATCATCCCCCCT--CATGATGACTGCACCCCACGGAGG------AGTTGT-TCTGAAAC--CTACCGCATGCCAGTGA--T-------------GGAGTACAAGATGAATGAAGGGGTTTCGTATGAATTC--AAATTTCCTTTCCGAAATAACAACAAATGGCAGAGGAATGCCAGTAAAGGG------GCACATTCGCCTCAAG-TCCCAGCCCAGGTGC-AGA-GTCCCATGACCTCAAGGATGAATGCTGGGAAAGGAGATGGCAAGATGCCTCCTCCGGAAAGAGCTACCAACATTCCTCGAAGCATCTCTAGTGATGGGCGTCCCCTGGAGA---GGCGGC--TGTCTCCTGGTTCAGACATCTATGTGACGGTCTCATCCATTGCTTTAGCGAGATCCCAGCAGTGCCGTAACTCCCCTAGCAACCTGTCGTCCTCCAGTGAAACTGGCTCCGGTGGGGGCACTTACAGACAGAAGTCCATGCCTGAA-GGATTTGGAGTAAGCCGCAGATCCCCAGCCTCCATTGACAGGCAGAACACCCAGTCAGATAT---TGG-TGGCAGCGGAAAATCCACACCCAGCTGGCAAAGGAGT-GAGGACAGCATTGTTGACCAAATGGCTTACAGTTATAGAGGACCTC--AGGATTTCAATTCT---TTTGTCCTCGAGCAGCATGAATAT-ACAGAGCCAACGTGCCACCTTCCAGCAGTATCAAAGGTACTGCCAGCTTTCCGAGAGAGTCCCAGTGGGAGGTTGATGCGGCAGGATCCAGTGGTGCATTTGTCTCCGAACAAGCAAGGGCATTCTGA-TAGCCACTACTCGAGCCACTCCAGT--AGCAATACCCTCTCCAGCAACGCCTCGAGTGTC---CACAGTGATGAGAAGTGGTATGA---TGGGGACCGCACTGAATCCGAACTCAGCAGC-TATAACTATCTG-CAGGGCACCTCTGCGGACAGTGGCATTGAC-ACCACCTCCTA-TGGCCCCGGCCAC-GGCAGCACAGCCT-CCCTGGGGGCAG-CCACATCGTCGCCTCGCT-CA---GGGCCAG-ACAAGGAGAAGGTGGCACCTTTATGGCAAAGCTCCAGTGAAGTGATC---TCCTGTGCCGACCGGACTTTAGAGGCAGACGGTCATGGCATGGACCGGAAGGCAGAGTCT---CTGAGCTTGGATCTCCACAGCAAGAGCCAAGTAGGCTCGGGCCCTCTGACAAGGGAGAACAGCGCTTTCAGTATTAGTGACGCTGC-AT-CCCACACAAGTACTATGAGCTCCCGACACTCG-GCCAGCCCCGT-GGTTTTCACCAGTG-CCAGGAGTTCC-CCTAAAGAAGAGCTTCATCCTGCCACCCCCTCGCAGCTCACACCTTCGTT------------CTCCTCGTCTTCCTCAAGTGGTCCCAGGACTTTCTACCCTCGCCAGGGTGCT-ACAAGCAAGTACCTGATTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNN------TTCCACGCGCTCTCGTCTCCTCAGTCTCCTTTCCCTCCCACGCCAGCCTCTAGGCGTGC-CTTGCACAGAACACTGTCGGATGAGAGTATTTACAGTGGCCAGAGGGAGCACTTCTTCACCTC--GAGG-GCTTCGCTTCTGGACCAAGCC---CTGCCCAATGATGTCCTCTTCAGTAGCACATATCCTTCTCTGCCCAAGTCTCTTCCACTGCGGAGACCCTCTTAC-ACCTTAGGGATGAAATCATTGCAT-GGAGAGTTCTCTGCTTCAGACAGCTCCCTCACGGACATCCAGGAGACCCGA---AGGCAGCCCATGCCTGACCCTGGCCTGATGCCCTTGCCTGATGCTGCT-GCAGA-CTTGGATTGGTCCAACTTGGTGGATGCTGCCAAGGCCTATGAG------GTCCAGAGAGCCTCATTCTTTGCTGCTAGTGATGAAAACCATCGCCCCTTAAGTGCTGCATCCAACAGCGATCAG-CTTGAGGACCAGGCTCTGGTTCAGATGAAGTCATACAGCAGCAGTAAAGATTCTTCTCCC---ACTCTGGCTTCTAAAGTGGACCAACTGGAAGGCATGCTGAAGATGCTTCGGGAAGA-TTTGAAG---AAG--------------------------------------------------------------------------------------------------------------------------------------------------------------");
    // refs.push_back("GATCGAGTACAGATAGCATAACTTCTGGGTAAGGCACTTGAAGAAGTCGAACAGCCTCCGTAAGACCCTGAGGAGAACGTCGGAGTCCAACAGGAGTGCGTTCACGACTTGAAGCCGGACTTCCACCCGAAACAGAAGAAAAAG-------------GAA---GAACTCTAGAAGGGCCTCGTCGAAGTCGTATGGAAGGTCAACCAGGTGAAATCTTCGGTCTCACCCTCTCCTTAGAAATGACNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNGACCAGGAGTTCGACTAGTGACAACCTGCGTC-GTGAGTCCTCGGCCACCAAAAGT--AGTAATCGTCGTTTCTTGCTCCGAGAGACCTG------GAGTATCCGAAACCGACGTAGGTGGTTCAACCTGGTCAGGTCTAGACGTCGACGTAGCCCGTTCCCGTAGTCCGGCCCCAGCCCGTACCCGACGGCAGCCCAGAGGAC-GTACAGTCACTCCCTCGACAGCCTCCGGCTCTTAAGGGGTACGTCCCTGAAGTAGGGATTCCACATGCTTCCGGAGGCGTCACCTTCTCTGAAC--CCCTCCCTTCCCATGCACGATCACTTCTCCTGCAGTAACCCGTCCCGAACCAGGTCTTCACTGCG-GGAC--CTCCCCTTCTTCCCAAGGGAGAC-CGGTGACATCT-ACGAGAGCAGGCTGTCGCAGGCCACGTTCCGGGGGGAACTCCACCCCCACGACCCCTTTCCCCTGACTCCTCTCCTC-TCCCGTACCTT------GAATTTTAAGAAGACACTGAGTCCCCAACCCGAAAGGTCCAGGTAGTCAAAGAACTCCGGGAGGACCTGTCACCTGAACCTGCAACCGACGCCCCTACGGGCGTCCAGGGACCACCTCCGGACGTCCGCCCACACCCGGTAGAGTGACGGTAGCGAGACTACTGCGAAAGAGCACAGGTATTTAGGGTGCCTCAACTACCAGGGGAGCCCAAAGAAGGTAGGTTAGTCCA-TGAACGAGCACCGTGGGACCGCGCCCATCTTGCAGGAGTCCGGCCTCCTCCTCCTGCTCCTCCTCCTCCTC---TTCCTCCCGCGGTCGACGCTCCTCGACCGGCCTACCTCGAGAAGAAACCCCCTCGAGGACCGTGACGACTTTTGGTGGCCCGACCGT--CTCACAGCCCTCGAGTACCATGAACACACCCTTCGCCGTAGCAA-CTACGACTTTCACGACAAGAGGGCGCAGTCTCCTGGGCGCGGGTCGACCGAGAACGACACCTATAGGTTCGAGTCCCTTCT-GAGACGAAAGGCCAGGTCCGGCACCGAGAGACGGAGATTTCAGGCCA-GACGGTCCCTCTAGTGAAGTGACCTCGACACGGTGTCCCCCCGGTGAAAGAGGAACGGCCCGGG---CCTCGCGCCCCTGCGGCACCGTCGGGGCTCGCTCCGACACGACGGCACCGACCCCGGCATCCTCCACCACAGCTACGGCGACAGGCGTCTCCACGGGACGTCTATCAATGTCGACGACTCAAGTCTGAGTCCCACCAGGGGT---AGTATGGTGAACAGCAGTGACACCCGCGAGCTGCGCAACGACCTCTCCCATAACGATGACCTCACCGAGCTCATCACCGATAGTCTTACCGGAACGGACAAACCTCTGTTCACGTGGTGACCCAGGACGGCGCAGTCAGAGGGTGACCCCGAGAGAGCCTTTCGGCCGTCGTGGAAACTATGACGCCCCTCCACCGTGCAGCCGAGACACATAAGTACGACGAGCTCCTGTTTTCTTAACTTTAGGACTTCAGGAGATA-TTGACATTCGGTAGACCAGTCGTTACGACAGAAGCGGGGAAACGGTCGACCCGCACCTAAAGGGCGACGGTGAC---TCCAGACTGACGCACAAGACGGCCAGCTACCTCCGGCCCCTGGACGCCGAGTGAGGCTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAG---AGAGGTCCCCTGCGGGT-AGTGACCTCTACGAAGAGCCTTACAACCGCCGAGAAAGGCCCCCTCCGTAGAAGGGCAGAGGAAAGGGTCGCAAGTCGGCTCGCCAGTAC-CCCGAGACGTGAACCCTCCCCC-GGACTCCGCTCACGCC------GGGGAATGACCGTGAGGAGACGGTAAACAACAATAAAGCCTTTCCTTTAAACTTAAGCATCCTTTGGGGAAGCAAGTAGAACATAAGG-------------CGGTGACCGTACGCCATCCGAAGTCTCGTTGA------GGAGGCTCCCCACGTCAGAAGTACCCCCCCCTACTACTGGTGGAAGTGGCACTGGCTGCAAGAGTCCTCCAGTTAGTAGACGAGCACCGAGTCTCACCGGTGCCGGTGAAACGTCTAGAGGTGGTCGGCCGACGGGACGGAGTCCGGACGGACGGTCCGTATCGGCATCCCGAGGTGCAGTCGGTGCTACGGGAGCATTAACTGCACCTTCGGCTCGACCGGATCGGGCAAGGAGGCGTCCCAGTAAAGGTAGCTGAGCGTC-GGAAAACTGTGT-TCGACATTGGCGAACTGTTAGAGAAACTATAGGAGC---AACGACTCCTTTGAGAGGTACCTTTGTGTAA-GAGGAG-CGAGCATCTTCTAAAAC--TCCGACCGCAGCCACCTCCAGGTAGGATAGTGTAGAGCCGTCCTCGTCAACTTTTGGTGCGAGAAACAAAGGACGAGCTACTCTTGTTACTTGAGTAACCTCTAAGGGTTCTCTGTCAGATAAAGGTAGCGGACCAACAT-C--AGAAACCGGGC-CTGACGG-GTGTGTTCCCGTGGGT-ACGATGACTCGAGTCGAGGACCTATGCCGAACCTGAAAAGGAAGAAGAATCTCCGGTTTCTCTACTTGCCCTTGAACGGTCTTCCT-AGCTACCCCCACAACCACTGTAAGAAAAGACGGTCCAGAAAGTCCATAAGGACCGCCCAGGCCCAACGGTACCGGGCCTTGAAAAGTCTAAATACTCGGCGCAAAAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGACTTATGTAAACTAAATCCCTTGCACTGAGGAAACCCCTACCCACCGGGTTTCCT-TCCGTGTAGAGACCTGGACCACTGACGGTGTGACATTGTGTGTGGCAGACACGTCCCCGACACCTGGGACTGCTACTGCTTGTGCACGACTTTCACTCTCGCCTAC-AAAAAGCCCGACTTACCGACGCGAGGTCCAAGGACCTTTTGATAACACTGCTATAGTAAGGGTTACACGGCGAACGCCTCCTC---------------------------GACAACGAACAACAATCCGCAC--ATCCCGTCGTACCACCTTTGCACCTTGTATTAAAGTATCAGAAACATCCAACACATGTCTCTCACCCAAGGTCACCTCAGGCAGAACCACAGCTCGACCCGGGCCATAAAGAGCTTCGGAAACTCAGCTTACGCAAGGGGGTC-GTCGACGTCCTTGAGGAGTTTCCGACTCGACCGCCTAAGCAACAACATGTAGAGGAGGAGCCACGAGACGGGCCGAGACGTCATGTACTACGGCTGGAAGACGA-CTATCAAGTCG-----------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCGAAGTAC----TCG-------------------ACGAGGC---------------ACTGGAACCCCC---------------------------ACAACT--TC--------------------------------------------------CGGTTCGCGTCCGTTACGTGTAAC-------TCGAGCCCTTAGTGTACGAGGTCGTGAAGGAAGTCCCCTTCCGGAGACCGGCTCCTC--ACGA---ACCGCCA---ACTCCCCTACCGTAGGA-GGTTCT-GACTAGGGGA-GTCACAGTACTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNN--NGGAAGACTATCTTC-TTAAAGGAT--ATTATC---ATACGTGGCTCTAGGT-GTACGA-----GGTGCTACATGGAAAAGTGAGACAGGTCCACGTTATGGTCCG-AGAGGAACCCCTGAAGTTCCTGGCGGTGCG-GCCGTAACC-A-CGTTACCCTCGACTCCCAGCTGAGCTT-CCTCCGTGAAAGTGTGGGTGACTTTCTTGGT-CTCGAACTGAAGCTGTCCGACTAGAAGGAAAGGGGAAGGGGTGGTTAGAGTAACGCCTTC---ATCCCCGTCGAGTACTGGTCAAGCAACGAGAATAGCAGCGGGACCAGGTAGGGCTCC------GACGGAAACCTTAAGTCCAG-GAGACACGACGGGTACCCCCTCGACTTCGACCGGCTCACCCTGTTTCCCGGGCTCTGGTCCCTCCGATGTCGCCGCCGACTCCGCGAGCACCACCACAAAG-AGGAGAATTAC-TGTAACA---CGGA-CAAGTACTACCGGAGCAAGTCCAGTTTATTATATGAGACCTGTAGCATCAC-CCGTTTCGTGAAACCCGTGCAGGTCCCGGACTGGTTC---AGAAGCCTCG-CCAAC-AGTAGGTCCAGGCTACT---GAAGGGCT---CGAG---GAGGGGAAACCGTAACGCGTCAAATGCCTTTTACCTGCTCAGGGGCCAAAGGCTGAACCTTGCG-GCGAAGTCCCCG-AAGAGGGCA-AGGAATTTCTCGGACTCCAACAGGGAGACCCTCTACCCGGGT---AGCTATTGTAGGTCCCTCTTCGGCCTGGGGAACGAGGGGGGCCGCTACTCCTTCAGTGACTCGAACCAGCCCCAACCTGGAG---C-CAGCCTAAACAGTAGGAAC-ATCGGGAACTTGTTCAGCTTCTTCGAAAGGGGCCTGCAGGGGACGAACAGCTACCTCCTCCACGACGGTATGAGGGACACGTCCCTGCCGGGGCAGAACA-TCCACCCCCTC--TGCG-TGAGTAGTT-TCGACAGGTGC---AGTTCGAGCGACTACCACTACAGTGACAACGAGGC--GACGGCCTAGGACGCGTCCCGAAAGGCGCCCCTC-GAC-CCCATCCGAAGCCCGT-ACTCCTT-G--GCCCTC-AGGTACGAGAGGCGGCCGGCGCGGA-ACGACGGCTCGCACA---AGACCTACGAAGACTCGTACCGACTCAGTAACCTTAACTCCGAGTGAC-TC-GA-AGGGACTGAGTGTCCCGACGAACACCTCCTGTACGAGAGTTCCGTTGACCTCCAAAGGTACAGGACCGAGCTTGCCGA---CCGAAGAAAGTGCAACAGGAAGGACCCCCCGGTTAGACGCTAG-GA------GCGGGATTGGGGGTAGAAACCCTGAGGTCCCCACCTCCACTACACCACTCCCGGAAGTGACGCTCCCCCGGGTTGCCCGCGGTACTGTCTCCT-AGGATTCGACGGCAAGACCCTGGACTTCGCGGCGTACATCTTCAGCAGACACACCTGAAACC---CTCGAGGTAGACGCGGGTGGGGACTCCGGGCGAGCGACCGGTGCCCGGAAAGGCAGACGCTGGCGAAGTCCGACCAGTA---------------------------------------------");
    // refs.push_back("GATCGAGTACAGATAGCATCCCTTCTGGGTAAGTCACTTGAAGAAGTCGAACAGACTCCGTAAGACCCTAAGGAGGACGTCGGAGTCCAACAGAAGGGCGTTCACGACTTGAAGTCGGACGTCCACCCGAAACAGAAGGAAAAG-------------NNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACGACATTCTGAAGTAGACCCGGTATCG-GACCAGGAGTTCGACTAGCGACAACCTACGTCGT-GAGTTCCCCGCTACCAAAAGT--AGTGATCGTCGTTTCATGCTCCGAGAGACCTG------GAGTATCCGAAACCGCCGTAGGTGGTTCAACCTGGTTAGGTCTAGGCGTCGCCGTAGCCCGTCCCCGTAGTCCGGTCCCAGCCCGTCCCCGACGGAAGCCCAGAGAACC-TACAGTCACTCCCTCGACAGGCTCCGACTCTTGAGAGGTACGTTGCTGAAGTAAGGATTCCACATTCTTCCGGAAGCGTTCCCCTCTCTGAAC--CCCTCTCTTCCCATGCACGATGACTTTTCCTGTAGCAACCCGTCCCGAACCAGGTCTTCGCTACG-GGA--GCTCCCCTTCTTCACGAGGGAGAC-CGACGACATTT-ATGAGAGTAGGCTGTCCCAAGACACGTTCCAGGCGGAACTCCACCCACACCCTCCCTTCCCTCTGACTCCTCTCCTC-TCGCGCACCTT------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--TACGGGAACAAACAAACCTCTGTTTACTTGGTGACCTAGGACGGCGTAATTAGAGGGTGACCCCGAGAGAGCTTTTCGGCCGTCATGGAAACTATGACGACCCTCTACCGTGCAACCGAGACATATAAGTACGACGAGCTCCTGTTTTCTTAACTTTAGGACTCCAGGAGATA-TTGACATTCGGTAGACCAGTCGCTACGACAGGAGTGAAGAGACGGTCGACCCGCACCTAAAAGGTGACGGCGGT---TACAGACTGACCCACAAGACGGACAGTTACCTCCGACCCCTAGACGCCGAGTGAGGTTTGGGAAGACTGTACTTAAAAACAGACATTCACGGGG-GGTGTCTCGGTCAGAGTGACCTAC---TATCCAACGATCCTCTCAATGCCGTGACGACCCTAGAACGATTTCGTTACCTACTCTGGCAGTGTATCTACAGGCTTGGTTCTCTGTCGNNNN---NNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NGAGGCACCCCATGTCAGTAGCACCCCTCCTTACTACTGCTGGAAGTGACACTGTCTACAAGAGTCCTCTAGTTAGTAGACGAGTACCGAGTCTCACCGGTGACGGTGAAACGTCTAGAGGTGGTCGGACGACGGGACTGCGTCGGGACGGACGGTCCGTATCGGTATCCCGAGGTGCAGTCGGTGCTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATCGGGTAAGGAGGCGTCTCAGTAGAGATGACTAAGTGTT-GGAAACCAGTGT-TTGACGTTAGAAAACTGTTAGAGAAACTATAGGAGT------AATTACTTTGAGAGATGTCTTTATGTAA-GGGGAG-CGAGTATCTTCTAAAAC--TCCGACCACAGTCAACTTCAGGTGGGTTAGTGTAGAGACGTCCTTGTTAACTTTTGGTGCGAGAAACAAAGAACAAGTTACTCCTGCTATTTGAGTAACCTCTAAGGGTTTTCCGTCAGTTAAAGGTACCGAAACAACATC---AGAAACCGGGCTTA-ACGG-GTGTGTTACCGAGGAC-ACGATGACTCGAGTCGAGGACCTATACCAAACCTGAAGAGGAAGAAGAAACTCCGGTCTCTTTACTTACCCTTGACCCGTCTTCCT-AGCTATCCCCACAACCACTGTAAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACATTGTGTGTGACAGTCACGTCCCCGATACTTGGGCCTGTTACTGCTTTTGTACGACTTTCACGCTTGCCTAC-AAAAAACCCGACTTACCGACACGAGGGCCAAGAACCTTTTGATAACAATGCTATAGTAAAGGTTACACGGCGAACGCATCCTC---------------------------GACAACGAATAACAATCCACAC--ATGCCTTCGTACCATCTTTGTACCTTGTATTAAAGTATCAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGTCAAAACCACAGTTCGACACGTGCTATGAAGAGTTTAGGAAACTCGGCTTGCGCGAGAGGGTT-GTCGACGTCCTTAAGAAGTTTACGACCCGGCCGACTAAGTAACAACATGTAGAGAAGAAGTCACGAGACGGGTCGAAACGTCATGTACTACGGATGAAAGACGA-CTATCAAGTCG-----------------------------------------------------------------------------------------------------------------------GGAACGAGTAGTTCAAAGTAC----TCA-------------------ACGAGAC---------------ACTGAAACCCTC---------------------------ACAACT--TC--------------------------------------------------CGGTTGGCGTCCGTGACCTGTAAC-------TCGAGCCCTTAACGCACAAGGTCGTGAAGAAAGTCCCCTTCGGGAGACCGACAGCTC--ACGA---AACGACA---ACTCCCTTACCGTAGGA-GGTCCT-GGCTTGGGGA-GTCACAGTACTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNN-NNNNNNNGAT--ATTATC---ATACGTGGTTCTAGGTGC-ACAA-----GGTGTTATATAGAAAAGTGAGATAGATCTACGTTGTGGTTCA-AAAGGAACCCGTGAAGTTCATGGCGGTGAG-GACGTAAAC-A-CGTTACCCTCGATTCCCGTCTGAGTTT-CCTCCGTGAAAGTGTCGGTGATTTCCTTGGT-CTTAAACTAAAACTTTCCGATTAAAAGGAGAGCGGAAGAGGTGGATAGAGTAAGGCTTTT---ATGCCTGTCGAGTAATGTTCAAGCAACGAAAATAGTAGAGGGACCAGGTACGAATCC------GAAGGAAACCTTAAGTCCAG-GAGACACGACGGGTACCCCCTCGACTTCGACCGACTTACCCTGTTTCCAGGTCTCTGGTTACTACGGTGCCGTCGACGGCTTCGGGGGCACCACCATAAAG-AGGAGAATTAT-TGTAACA---CGGA-CAAGTACTAACGGAGTAAGTTTAGTTTATTATATGAGACCTGTAGTATTAC-CCGTTTCGTGAAACCTGTTCAGGTCCCAGACTGTCTC---AGAAGACTAG-CCAA-TAGAAGTTCTAGACTGCT---AAAGGGTT---CAAG---TAGTGGAAACCGTAACGCATTAAATGCTTTTTACCTACTCAGAGGTCAAAGTCTGAACCTTGC-AGAGAATTCACCAAAA-AGGGAA-AGGAATTTCTCGGACTCCAAGAGGGAGACACTCTACCCGGGT---AGATATTGTAGGTTTCTCTTTGGTCTTGGGAACGGAGGTGGTGGCTATTCTTTCAGTGACTCGAACCAACCTCAACCTGGAG---C-TAGTCTAAACAGTAGAAAC-ATCGGGAAATTATTTAGTTTTTTCGAAAGAGGTCTGCAAGGGACGAATAGTTAACTTCTACACGATGGTATGAGGGACACGTTACTACCAGGTCAGAATA-TACACCCACTC--TATG-TAAGTAGTT-TCGATAGGTGT---AGGTCAAGTGAATACCACTATAGTGACAACGAGGC--GACGGCATAAGACGCTTCTCGAAAGGACCCTCTC-GAC-CCTATCCGAAGGCCGT-ACTCTTT-AGA--CCTC-AGGTACGAGAGGCTGCCAGAACAGA-ACAAAAAGTCACACA---AGACATACGAAAAGTTGTACCGACTCAGTGACCTCAACTCCGATTGTC-TT-GA-AGGGACTGAGTGTCCCGACGAAAACCTCCTGTACGAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATGCAAAAGAAAGGAACCCCCGGTTAGACGTTAG-GA------ACGGGATTGGGGGTAGAAACCCTGAGGTCCCCAGCTCCAATATACCACTCTTGGGAGTAAAGCTCCCCCGGGATGTCCTCGGTAGTGACTCCT-AGGATTCGACGGTAAAACCCTGGACTTTGCGGCGTACATCTTCAGTAGCCACACCTGAAACC---CTCACGGTAGAACAGGCTGCTACCTCCGGGAGAGGCACCTTTGTCCGGAAAGGCAGACGCTGGAAAAGTTCGACCAGTA---------------------------------------------");
    // refs.push_back("GATCGAGTACAGATAGCATAACTTCTGGGTAAGACACTTGAAGAAGTCGAACAGACTCCGTAAGACCCTGAGGAGGACGTCGGAGTTCAACAGGAGGGCGTTCACGACTTGAAGTCGGACCTCCACCCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGGGCTTCGTAGAAGTCGTACGGAAGGTTGACCAGGTGAAATCTTCGGTCCCACCCTCTCCTTAGAAATGACGACGACATTCTGAAGTAGACCCGGTATCA-GACCAGGAGCTCGACTAGCGACAACCTACGTCGT-GAGTTCCCCACTACCAAAAGT--AGTGATCGTCGTTTCTTACTCCGAGAGACCTG------GAGTATCCGGAACCGTCGTAGGTGGTTCAACCTGGTTAGGTCTAGACGTCGTCGTAGTCCGTCCCCGTAGTCCGGTCCCAGCCCGTACCCGACGGAAGCCCAGAGAACC-TACAGCCATTCCCTCGACAGACTCCGACTCTTGAGAGGCACGTTACTGAAGTAAGGATTCCACATTCTTCCGGAGGCGTTCCCTTCTCTAAAC--CCCTCTCTTCCCATGCACGATGACTTCTCCTGCAGCAACCCGTTCCGAACCAGGTCTTCACTACA-GGA--ACTCCCCTTCTTCACGAGGGAGAC-CGACGACATTT-ATGAGAGTAGGCTGTCGCAAGACACGTTCCAGGCGGAACTCCACCCCCACCATCCCTTCCCTCTGACTCCTCTCCTC-TCACGTACCTT------GAATTTTAAGAAGACACTAAGCCCTCAACCTGAAAGTTCCAGCTAGTCAAAGAAATCTAGGAGAAGTTACCACCTGAACCTCAAACCGAATCCCCTACGGGGTTCCAGGGACCAACTCCATGCGTCGGACCACACGCTATAAAGTAACGGTAGTGAGACTACCGCGAAAGCACACAGGTATTTAGGGTGTCTCAAATACCAAGGAAGTCCGAAGAAGGTAGGTTAGTCCATGAA-CGATCATCGCGGGACCGCGCCCATCTTGCAAGACCCGGGTCTACTCCTACTCCTCCTCCTCCTC------TTCCTCCCCTGCTCAACACGTCTTCACCGTCCTACTTCGAGGAGAAACCCACTTGAAGACCGTGACCACTTTTGTTGACCCGACCGT-CTCACAGCTCTCGAGTACCATGAACACACCCTTCGTCGTAGTAAA-TATGACTTTCACGACAAGAGGGAACAGTCTCCCAATTACGGACGAACCGAGAACGACACCTATAGGTTCGAGTCCCTTCTG-AGACAAAAGGCCAGGTTCGGCACCGAGAGACAGAGATTTCAGGCTAGACGGTTACTCTAGTGAAGTGACCTCGACACGGTGTCCCCTCAGTGAAAGAGGAACCGACCGGG---GCTCGCTCCACT---ACATCGACGGGGGTCCCTCCGACACGACGACACCGACCCGGGTATCCTCCACCACAGCTACGGTGACAGTCGTCTCCACGGAACGTCTATTAATATCGACAACTCGAGTCTAAGTCACGCCAGAGGT---AGTATGGTGAAGAGCAGTGACACCCGTGAGCTACGTAACGACCTCTCCCATAACGATGACCTCACCGAGCTCATCACCGATAGTCT--TACGGGAACAAACAAACCTCTGTTTACTTGGTGACCTAGGACGGCGTAATTAGAGGGTGACCCCGAGAGAGCCTTTCGACCGTCATGGAAACTATGACGACCCTCTACCGTACAACCGAGACATATAAGTACGACGAGCTCCTGTTTTCTTAACTTTAGGACTCCAGGAGATA-TTGACATTCGGTAGACCAGTCGTTACGACAGGAGTGAAGAAACGGTCGACCGACACCTAAAAGGCGACGGTGGT---TACAGACTGACCCACAAGACGGACAGTTACCTTCGACCCCTAGACGCCGAGTGAGGTTTGGGAAGACCGTACCTGAAAACAGACATTCACGGGG-GTGGTCTCGGTCAGAGTGACCTACTTCTGTCCAACGATCCCCTCAATGCCGTGACGACCCTAGAACGATTTCGTTACCTACTCTGACAGTGTATCTACAGGCTTGGTCCTCTTTCGGCGG---AGAGATCACCCGCGGGT-AGTGACCTCTACGAAGCTCCTTACAACCGTCGAGAAAGACCTCCTCCGTAGAAGGGTAGAGGAAAGGGTCATAAGTCGGCGCTCCAGTACCCT-GAGACGTGGACCCGACCTT-AGACTCCCCTTACTCA------GGGGAATAACCGTAAGGAGACGGTGAACAATAATAAAGCCTTCCCTTTGAACTTAAGCATCCTTTGAGGAAGAAGTAAAAACATGAGG-------------TAGTGACCGTACGCCATCTAAAGTCTCGTTGA------GGAAGCACCCCACGTCAGTAGTACCCCTCCTTACTACTGGTGGAAGTGACACTGTCTACAAGAGTCCTCTAGTTAGTAGACGAGTACCGAGTCTCACCGGTGACGGTGGAACGTCTAGAGATGGTCAGACGACGGGACGGCGTCGGGGCGGACGGTCCGCATCGGTATCCCGAGGTGCAGTCGGTGTTACGGGAGCATCAACTGTACCTTCGGTTCGACGGGATCAGGTAAGGAGGCGTCTCAGTAGAGGTGACTAAGTGTT-GGAAAACTCTGT-TTGACGTTAGAAAACTGTTAGAGAAACTATAGGAGT------AATTACTTTGAGAGGTGTCTTTATGTAA-GAGGAG-CGAGTATCTTCTAAAAC--TCCGACAACAGTCAACTCCAGGTGGGTTAGTGTAGAGACGTCCTTGTTAACTTTTGATGTGAGAAACAAAGAACAAGTTACTCTTGCTACTTGAGTAACCTCTGAGGGTTTTCCGTCAGTTAAAGGTACCGGAACAACATC---AGAAACCGGGCTTA-ATGG-GTGTGTTACCGGGGGT-ACGATGACTCAAGCCGAGGACCTATACCAAACCTAAAGAGGAAGAAGAACCTCCGGTCTCTCTACTTACCCTTGAACGGTCTTCCT-AGTTATCCCCACAACCACTGTAAGAAAAGACGGTCCAGAAAGTCCATAAGGACCGCCCAGGCTCAACGGTACCGGGCCTTGAAAAGACTAAATACTCGTCGTAAAAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGACTTGTGTAAACTGAATCCNNNNNNNNNNNN----------------------------------------------------------------------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACCTTTTGATAACAATGTTATAGTAAGGGTTACACAGCGAACGCATC------------------------------GACAACGAATAACAATCCACAC--ATACCTTCGTACCATCTTTGTACCTTGTATTAAAGTATTAGAAACATACAACATATGTCTCTTACCCAAGGTCACCTCAGTCAAAACCACAGTTCGACACGTGCTATGAAGAGTTTAGGAAACTCAGCTTGTGCAAGAGGATT-GTTGACTTCCTTAAGGAGTTTCCGACCTGGCCGACTAAGTAACAACATGTAGAGAAGAAGGCACGAGACGGGTCGAAATGTCATGTACTACGGATGAAAGACGA-CTATCAAGTCG-----------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTAC----TCG-------------------ACGAGAC---------------ACTGAAACCCTC---------------------------ACAACT--TC--------------------------------------------------CGGTCGGCGTCCGTGACCTGTAAC-------TCGAGTCCTTAGTGCACGAGGTCGTGAAGAAAGTCCCCTTCGGGAGACCGACAGCTC--ACGA---AACGACA---ACTCCCTTACCGTAGGA-GGTCCT-GGCGTGGGGA-GTCACAGTACTCGAGTGATCAAGCTTTTTAATAAGCCATCAACATGCCACAAGGTAAAAGAAAGTAAAGAAGACCGAAAAGGGAAGATTACGAGTGTCGGTGACCTGGTTCTAAGAG-TAGTCGGGGTTTTATCAAGGTCACA--AGAAAGACTATCTT-TTTGAAGGAT--ATTATC---ATACGTGGTTCTAGGTGC-ACAA-----GGTGTTACATAGAAAAGTGAGATAGATCTACGTTGTGGTTCA-AAAGGAACCCGTGAAGTTCATGACGGTGAG-GACGTAAAC-A-CGTTACTCTCGATTCCCGTCTGAGTTT-CCTCCGTGAAAGTGTCGGTGATTTCCTTGGT-CTTAAACTAAAACTGTCCGAATAAAAGGAAAGTGGAAGAGGTGGATAGAGTAAGGCTTTT---ATGCCTGTCGAGTAATGTTCAAGTAACGAAAATAGTAGAGGGACCAGTTACACCACCTTAAGCGAAGGAAACCTTAAGTCCAG-GAGTCACGACGGGTACCCCCTCGATTTCGACCGACTTACCCTGTTTCCAGGTCTCTGGTTCCTACGGTGTCGTCGACGGCTTCGAGGGCACCACCACAAAG-AGGAGAATTAC-TGTAACA---CGGA-CAAGTGTTAACGGAGTAAGTTTAGTTTATTATATGAGACCTGTAGCATTAC-CCGTTTCGTGAAACCTGTGCAGGTGCCGGCCTGTCTC---AGAAGACTAG-CCAA-TAGAAGTTCTAGACTACT---AAAGGGTT---CAAG---TAGTGGAAACCGCAACGCATTAAATGCTTTTTACCTACTCAGAGGTCAAAGTCTGAATCTTGC-CGCGAACTCACCAAAA-AGGGAA-AGAAATTTCTCGGACTCCAAGAGGGAGACACTCTACCCGGGT---AGATATTGTAGGTCCCTCTTTGGTCTTGGGAACGAAGGTGACCGCTACTCTTTCAGTGACTCGAACCAGCCTCAACCTGGAG---C-TAGTCTAAACAGTAGAAAC-ATGGGGAAATTGTCTAGTTTTTTCGAAAGAGGCCTTCAAGGGACGAATAGCTAACTTCTACACGATGGTATAAGGGACACGTCGCTACCAGGTCAGAATA-TGCACCCACTC--TATG-TAAGCAGTT-TCGATAGGTGC---AGGTCAAGTGAATACCACTATAGCGACAACGAGGC--GACGGCATAAGACGCTTCTCGAAAGGACCCTCTC-CAC-CCTATCCGAAGTCCGT-ACTCTTT-GGA--CCTC-AGGTACCAGAGGCTCCCAGAGCAGA-ACGAAAAGTCGCACA---AGACATACGAAAAGTTGTACCGACTCAGTGGCCTCAACTCCGATTGTC-TT-GA-AGGGACTGAGTGTCCCGACGAAAACCTCCTATTCGAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAACTGTAAAAGGAAGGAACCCCCGGTTAGACGTTAG-GA------ACGGGATTGGGGGTAGAAACCCTGAGGTCCCCAACTCCAATATACCACTCTTGGAAGTGAAGCTCCCCCGGGATGTCCTCGGTACTGACTCCT-AGGGTTCGACGGTAAGACCCTGGACTTCGCGGCGTACATCTTCAGTAGCCACACCTGAAACC---ACCACGGTAGACAAGGCTGTTACCTCCGGGAGAGCGACCTTTGTCCGGAAAGACAGACCCTAGAAAAATTCGACCAGTA---------------------------------------------");
    // refs.push_back("ATGACCAGCCTGAAGCGGTCCCAGACCGAGAGGCCCGTGGCCACGGAGAGGGTGTCGCTGGTGGGCCCGGAGGGGCCGTCGGCCAAGGTCCACACTGACGACTTCTACATGAGGCGCTTCCGCTCGCA-GAACGGCAGCCTGGGCCCGCCGGTCATGGCCCCCGTGGGCCTCCCCCGCAGCGAGAGCACCCACCACGTCACCTCCACCCCGGGGGTCCCCAAGATGGGGGTCCGGGCG------AGGATCGCCGACTGGCCGCCGA-GGAAGGACACC-----GTCAAGGAGT---CG---GG-TCG-GCCGGGCCAGGAGGGCGAGCCGGCCACCGGCCCCGAGAGCCTGTCGTCCAAGAGCAGCCCCGGGAGCCAGGGCAGCTCGGTCAGCCTGAACTCGAGCGACGCCGCCATGTTGAAGAGCCTCCACAACAC--CCTGAAGAGCAAGACCAGGCCAGTGGAGCGGATGGACTCCCGCTTCCTCATGCCCGAGGCCTACCCGGCCTCCCCGAGGAAGGCCCTGCGCCGGATCCGGCAGCGGA-GCAACAGCGACATCACCATCAGCGAGCTGGACGCCGACAGCTTCGACGAGTGCATCTCGCCCACGTACAAG-AGCGGGCCGTCGCTCCACCGGGAGTACGGCAGCACGTCGTCCATAGACAAGC-AGGGGACCTCGGGGGAGAGTTTCTTCGACCTGTTGAAGGGCTTCAAGGGGGACCACCCGGAGCCCCGGGGAGGGACGCCCACCCCACTGGGCGACCTGCTGATGGCCAG-CGGGAGCAAGGGCTCCGGCT-TCTCCCTGGACATCGTCGA---CGGGCCCGCCGTCCCCCGGGAGAACCTCCGACTCTTCAAGGAGAGGGAGAAGCCGC-CCAAGAGACGCTGCAAGTCGGAGACGGGGGACTCGTCCATCTTCCGGAAGCTCCGGAATGCCAAGGGCGAAGGGGAGCCGGCGGGG--AAA---C-CGGCGGACCAGGACGACGGCCGGTCGGAGGA---CGGCGTCCGGCCCTGGACCTGCCCCAGGTGC--TTCGCCCATTACGACGTCCAGAGCGTCCTATTTGACCTGAACGAGGCGGTGGTGAACAGAC---ACAACGTCATCAAGAGGCGGAA-CACGACCACGGGGGCCTCGGCCGCCGCCGTGGCCTCCTTGGTTTCGGGCCCGCCGGCCC-ACTCGGCCAGCTTCGGCTCCCCCATGGGC-AGCACCGAAGACCTGAACTCCAAAGGGAG------CCTCGGGGTGGACCAGGGGGACGACAAGAGCAACGAGCTCGTCA-TGAGCTGCCCCTATTTCCGGAACGAGATGGGCGGGGAGGGGGAGAGGAAGATCAGCCTATCCAGGTCCA--ACTCGGGGACCTTCGGCGGGGGCGAAAACGCCGCGTTCGAGTCCACGCTTAGCTCCCACTGCACCAACG-CCGGAGTGGCCGTCCTCGAGGTGCCCAAGGACAGT---TTGGTTTTGCACTTGGATCGGGTCAAGCGGTACATCGTGGAACACGTGGATCTGGGCGCCTACTATTACCGGAAATTCTTCTATCAGAAGGAACACTGGAACTATTTTGGAGCAGATGAGAA-CCTCGGCCCTGTAGCCGTGAGCATCCGAAGAGAAAAACCCGAGGAAATGAAAG-AAAATGGGCCGCCTCACAACTATCGGATAATTTTCAGAACTAGTGAGATTATTACCTTACACACAGTAAGCACTCAAGAA----A-TACCATCGATCTCTCTTGAGCGG-ACGAGTACAGAGC-TGCCACCAAAGAAAAGGGGGGAATCAGTCAGCCCAGAGATAACACAAAAATGTTTGGCCCTTGCATTTTCCGTACCTAAGAATAAAGATGAGGAAATGGAATTAGTCAGTCAGTTAATACACTATAACAGTACTGATGGACTGATTCATTGTACTCTCCCAAGAGCACCATCAGAAACCATTCAAAACACAAAAAAGGAA-ATCAAAAGATTCATGTGTCTCAGTGAGTTTGG------------------------------------------------------------------------------------------------------------------AAAGCTCA-GCTACCAGCTGAAAGTGGGGATCATGTACAGCAAAGCCGGGCAGAGCACGGAAGAAGAGATGTACAATAACGAATCCGCCGGGCCAGGTTTCGAAGAG--TTCCTTCAGCTTTTGGG--AGAACGAGTCCGTCTCAAAGGCTTCGAGAAGTA-TCGC----GCCCAGCTCGATACCAAGACTGATTCAACTGGCACCTACTCTCTCTACACTACTTACAAG---GACTACGAAATCATGTTCCATGTTTCCACGAT-GCTACCCTACACGCCAAA-CAACAAGCAACAG---------------------------CTCTTAAGGA-AGCGGCATATCGGCAACGATATTGTGACGATCG-TTTTCCAAGAACCTGGTTCCCAGCCCTTCAGCCCAAAGAACATTCGCTCTCACTTCCAGCACGTTTTTGTGATCGTTCGGGTGCATAACCCTTGCACCGAGAATG-TGTGTTACAGCGTGGCCGTCACCAGGTCCAGAGACGTTCCTTCC--TTCG-GCCCCCCGATCCCAA-AGGGGGTCACC-TTCCCCAAATCCAACGTGTTCCGGGACTTCCTGCTGGCC--AAGGTCATCAACGCGGAGAACGCCGCCCACAAATCG--GAAAAGTTCCGGGCCA---TGGCCACCCGGACGCG---CCAGGAGTACCTG-AAGGACCTGGCGGAAAAGAACGTCACCAACACCCCCATCGACCCTTCCGGAAAGTTCCCCTTCATCTCCCTG-GCCTCCAAGAAGAAAGAGAAGTCGAAGCCGTACCCGGGGG--CCGAGCTCCACAGC-GCGGGGGCT-ATC-GTGTGGGGCGTCCACGCCAAAGACTACC-ACAAGGGGATGG-AGATCGACTGCCTGCTGGCCATCTCCAACGAGT-TCCTGGTCCTCATC-GAACAGGAAACCAAGAGC-GTG--GTGTTCAACTGTTCCTGCAGGGACGTCATAGGGTGGACCTCAACCGACACCAGCCTCAAAATATTTTAC-GAACGGGGGGAGTGCGTGTTCGTGGAGAGCTTCCAGAACAATCTGGAGGACGTCAGGGAGATTGTCAAAAGGCTGGAGTTTGTGACGAAAGGCTGTGA-GTCGGTGGAGATGACGCTGCGGCGGAACGGGCT-GGGCCAGCTGGGCTTCCACGTCAACTAT-GAGGGCATCGTGGC-GGACGTGGAGCCGTACGGCTATGCCTGGCAGGCCGGC-CTGAGGCAGGGCAGCCGGCTCGTGGAGATCTGCACCGTG--GCCGT-GGCCACCCTCAGCCACGAGCAGATGATC-GACCTCCTGAGGACCTCGGTCACCGTCAAGGTCGTCATCGTCCCCCCGCAC--GAGGACTGCACCCCGCGAAGGTTAGGCAAGCCCTCTAGACATCTTTTTTTGTCAGAGATCCCTGGTTGGGAGGGAGGGGAGGATGAGTA--------------TACTGAAATTTGTGAGTTTGGATTGACGA-----GGCAGAGAGA-------------------CCTCACCTCCTCCAGGAGGC-CTTCCCAGCTTGCGCACCCTTT-CACAGCAAGGCCAGGTGGCGCGGGGGGGTGGATAAAAATGGGGCCAAAATGTAAAATTGA--AAGAGCACC-AAAGGCTGTCAAAGAAGACAGTAGAACAAGTGAT----A---GATGG--TTGTCTCCTGGTTCGGACCATTACGTGACCGTGTCATCCATGGCTTTGGGAAGATCCCCGCAGTGCCGTAATTCTCCTAGTAGCCTGTCCTCCTCTAGCGACACTGGTTCTGGTGGAGGGACTTACCGACAGAAGTCCATGCCCGAAG-GGTTTGGAATGAGTCGCAGGTCCCCGGCTTCCATCGACAGGCAGAGTACCCAAACGGATAT---T-GGCGGAAGCGGGAAATCCACGCCCAGTTGGCAGAGGAG-TCAGGACAGCATTGGCGACCACTTG---------------------------------------------------------------------GAACCGACGTGCCATATCCCCGCAGTTTCCAAGGTCCTGCCATCGTTCAGAGACAGCCCTGGCAGCAGACTAATAAGACAGGGACCGGTGGTTCATTTGCCTGCAAGCAAGCAGGGGCACTCTGAC-AGCCATTACTCGAGCCACTCGAGCAG--CAACACCCTCTCCAGCAACGGCTCCAGCGCC---CACAGCGACGAGAAATGGTACGAGAGCGGCGAGCGAGGCGATGCGGAGCCCAGCGGC-TACACTT-TCCTCCAGGGCACCTCGGCCGACAGCGGCGTTGACA-CCACCTGCTA-CGGC-CCCTCCCACGGCAGCACGTCCT--CCCTCGGGG---CCACCTCGTCCCCCCGG-ACGGCCGGGCCCGGCAAAGACAAGGG-GGCCTCCCTCTGGCACAGCTCCAGCGAGGTCATCT---CCATCGCCGACCGGACGTTAGAGAAGGAGAGCCACGGGATGGAGAGGAAGGCGGAATCCTCCCTCAGCCTCGACCTCCACAGCAAGAGCCAACCGGGCTCCAACCCCTTAACGAGAGAGGGCAGCACTTATAGCATCAACGATGCCGCTTCGCACGCCAGTACCATGACTTCCCGACATTC---TGCCAACCCAGT-TGTTTTCTGCAGTG-CCAGAAGTTCGCCTAAAGA-GGAGCTTCATATTGCCACTTCCCTGCAACTCGCGCC------------------------------------TGGGCCTAGGACTTTCTATCCTCGCCAGGGGGCTACT-AGCAAGTACCTGATTGGATGGAAAAAACCCGAAGGGACAATAAACTCTGTGGGCTTTATGGATACCAGAAAACGCCATCAGAGTGATGGGAACGAAATAGGCCCCTCCAGGTTGCAGGCGTCTGTGAGAGAACTTCGGTCATCTCCCAAGCAGGCCGCTAAAGCCACCATCG--AAGAGGAGCTCAAGAGACTGATTGACCTCGAGAGCCAAACGTTGGAGTCCCAGAAGAATTTTAAGAAATGTTTTAGAGTATGTATCG-CTCCCCAACCTCCATTCCCTACCACCCCGACCACGAGGCGGACACTG-CACCGGACCCTGTCGGATGAGAGCATTTACAGCGGTCAGAGGGAAAAGTATTTCAACTCC--CGCA-CCTCCCTTCTGGACCAAGCCC---TGCCCAACGATGTCCTGTTTAGTAGCACCTACCCTTCTCTGCCCAAGTCTTTACCTTTGCGGAGGCCGTCTTAT-ACTCTAGGGATGAAGTCCTTACACAGT-GAGTTTTCGGCATCAGACAGCTCTCTCCCCGACGTCCAGGAGAACCGA---AGTCAGCCAATGCAGGACCCGGGCCTCATGCCCCTGCCTGACTCTGCT-TCTGA-CCTGGACTGGTCAAACTTGGTCGATGCTGCCAAAGCCTTTGAGGTCACAGTCCAACGAGCCTCATTCTTTGCTGCTAATGAAGAAAACCACCGACCGCTTAGCACTGCTTCCAATATTGATCAATTAGAGGATCAGTCTGCCGCATCAC-TGAAATCTTATGCAGGCAGTAAAGATTCATCCCC---TACTCTGGCTTCTAAAGTGGACCAATTGGAAGGTATGCTGAAGATGCTTCAAGAGGA-CTTGACG---AAG-------------GAGAAAGAAGACAAGGCTCACCTCCAAGCGGAAGTTCAGCACCTGAGGGAAGACAACCTCCGTTTGCAAGAGGAGTCCCAGAGCGCCTCGGAC-AAGCTGAAAAAATTCACCGAGTGGGTCTTCAACACCATCGACATGAGCTAA");
    // refs.push_back("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ACATCCAGTTGCCTTGATAGCATGTCCTCCAAAGGCAGTCCTGGAAGTCAAGGTAGTGCAGTTAACCTCAATTCCAGTGACTCTGCCATGTTGAAAAGCATACAGAACAC--TCTTAAAAGCAAAACAAGACAATCAGAGAACATGGATTCCAGATTTCTCATGCCTGAAGCTTATCCTAGCTCCCCTAGGAAAGCCCTTCGCAGAATAAGGCAACGTA-GTAACAGTGATATTACTATAAGTGAACTTGATGTGGATAGCTTTGATGAATGTATTTCACCTACATATAAA-ACTGGTCCATCACTGCACAGGGAATATGGTAGCACATCATCAATAGATAAGC-AGGGTACTTCAGGAGAAAGTTTTTTTGATTTGTTGAAAGGGTACAAAGAAGATAAATGTGATCAGAGAGGCCCAACTCCAACAAAACTGAGTGAATTTTTGATAGCCAG-TGGAAGCAAGGGTTCTGGCT-TTTCTCTGGATGTGATAGA---TGGACCTATTACACAGAGAGAGAACCTTAGATTATTTAAAGAAAGGGAGAAACCAC-TCAAAAGACGTTCTAAGTCAGAAACAGGGGACTCATCCATTTTTCGTAAGTTGCGTAATGCCAAATGTGAAGGGGAAC---TTGGA--AAA---T-CATCGGACCTTGAAGACAACCGATCAGAGGA---TTCTGTCAGGCCTTGGACTTGTCCAAAGTGT--TTTGGTCACTATGATGTCCAGAGCATATTGTTTGACTTGAATGAGGCAATTATAAATAGGC---ACAACGTTATAAAACGGAGAAA-TACGACCACAGGAGCATCTGCAGCAGCTGTTGCTTCTTTAGTTTCTGGACCTTTATCTC-ACTCTGCAAGTTTCAGTTCCCCTATGGGC-AGTACTGAGGACTTGAATTCTAAAGGAAG------CCTCAGCATGGACCAGGGAGATGATAAGAGTAACGACCTTGTAA-TGAGCTGTCCATACTTTCGAAATGAGATAGGTGGAGAAGGAGAAAGGAAGATCAGTTTATCCAAATCAA--ATTCTGGTTCCTTCAGTA---GTGAAAGTGCCTCATTTGAATCTACACTTAGCTCTCATTGCACAAATG-CAGGAGTGGCTGTTCTTGAAGTGCCCAAGGAAAAT---TTGGTTTTGCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAACACTGGAACTATTTTGGAGCAGATGAGAA-CCTAGGCCCAGTGGCTGTGAGCATACGAAGAGAAAAGCCCGAAGAATTGAAAG-AAAATGGACCACCATACAACTATCGAATAATATTCAGAACTAGTGAGCTCCTGACATTA-AGAGGTGCTGTGCTAGAAGATGCCA-TCCCTTCAACTGCAAAGCATTGT-ACAGCTCGAGGAC-TACCTTTGAAGGAAGTGCTTGAGCATGTGATTCCTGAGCT-------CAATGTTCAGTGTCTGAGGTTGGC--------------------------------------------------CTTTAACA--------------------------CCCCCAAAGTCA---------------CAGAAC-A-------------------ACTCATGAAACTAGACGAGCAAGG---------------------------------------------------------------------------------------------------------------------GNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NN-NN---NNNNNNNNNNNNNNNNNNNCTGATTCAACTGGTACCCATTCTCTTTATACAACATAC---AAAGACTATGAAATCATGTTCCATGTATCTACTCT-GCTGCCATATACACCTAA-CAACAAGCAACAG---------------------------NNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNTGTGGCTGTAACCAGGTCCAGAGATGTGCCATCC--TTTG-GACCACCAATCCCCA-AAGGAGTCACT-TTCCCTAAGTCAAATGTGTTCAGGGACTTCCTTTTGGCC--AAAGTTATTAATGCGGAGAATGCTGCTCA-TAAGTCT-GAGAAGTTCCGTGCCA---TGGCTACCCGGACCCG---CCAGGAGTACTTA-AAAGATTTAGCAGAAAAGAATGTAACCAATACTCCCATTGATCCTTCTGGCAAATTCCCCTTTATATCTCTT-GCCTCCAAGAAAAAAGAGAAATCCAAACCTTATCCTG--GAGCAGAGTTCAGTAGCA-TGGGTGC-CAT-TGTGTGGGGTGTCCATGCCAAAGACTACA-ATAAAGCTATGG-AGATAGACTGCCTTCTGGGCATCTCCAATGAGT-TCATAGTCCTCAT-TGAA---GAGACAAAGA---GTGTGGTGTTCAATTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTGTGACTAAAGGATGTGA-GTCAGTGGAGATGACCTTGAGGAGAAATGGGCT-CGGACAACTTGGCTTCCATGTCAACTAC-GAGGGCATTGTGGC-AGATGTAGAACCTTATGGCTATGCCTGGCA-GGCAGGACTGAGGCAAGGCAGCCGCCTAGTGGAAATCTGTAAAGTG--GCTGT-TGCCACCCTAAGCCATGAGCAAATGATT-GACCTTCTGAGAACATCAGTCACAGTAAAAGTGGTCATTATACCACCCCAT--GATGATTGTACCCCAAGAAGN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-----------------------NNN-NNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNN-NNNNNN-NNNNNNNNNNNNN-NNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNG--TTGTCTCCTGGTTCTGACATCTATGTGACAGTTTCATCCATTGCTTTAGCGAGATCACAGCAATGTCGTAACTCCCCTAGCAACCTGTCTTCATCTAGTGAGACTGGCTCTGTTGGGGGCACTTACAGACAAAAATCCATGCCTGAAG-GGTTTGGGATGAGCCGCAGATCCCCAGCCTCCATTGATCGTCAGAACACTCAGACAGATAT---C-AGTGGCAGTGGGAAGGCCACGCCCAGCTGGCAGAGAAG-TGAGGAGAGCATTGCTGACCAGATGGCTTACAGTTATAGAGGACCTC--AGGATTTCAATTCT---TTTGTCCTCGAGCAGCATGAATAT-ACAGAGCCAACATGCCATCTCCCAGCTGTATCAAAGGTGTTGCCATCTTTCCGAGAGAGTCCTAGTGGGAGATTAATGCGGCAGGATCCAGTTGTTCATTTATCACCAAACAAACAAGGGCATTCTGAT-AGTCATTACTCAAGCCACTCCAGCAG--CAATACTCTCTCCAGCAATGCATCGA-GTGCT--CACAGTGATGAGAAGTGGTATGA---AGGTGATCGCACAGAGTCGGAGCTCAATAGC-TACAACTAT-CTTCAAGGCACCTCTGCGGATAGTGGCATAG-ATACCACTTCATA-CGGC-CCCAGCCATGGCAGCACAGCAT--CTCTGGGGGCTGCCACATCATCCCCTCGCT-CA---------GACAAAGGAAAAGT-GGCACCACTGTGGCACAGCTCCAGTGAAGTTGTCT---CCATGGCAGATAGGACCTTAGA---AGAAAGTCATGGAATAGACCGGAAAACAGAGTCTTCACTGAGCCTGGATATCCATAGCAAGAGCCAACCCACCTCAAATCCTTTAACAAGGGAGAACAGTACTTTTAGTATTAATGATGCTGCTTCTCACACAAGNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNN-NNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGCCACCAGAGTGATGGCAATGAAATAGCCCATACCAGGTTGCGTGCCTCTGCCAGGGATCTCCGGGCATCCCCAAAGCAAACATCTAAGTCAACAATTG--AGGAGGATCTAAAAAAACTCATTGACTTTGAAAGCCCAACTCCCGAATCACAGAAAAATTTTA-AG------TTCCATGGGCTTTCCTCTCCACAGTCTCCTTTCCCATCCACTCCTACCACACGGCGGA-CCTTACATAGAACACTGTCAGATGAAAGTATTTACAGTGGTCAAAGAGAGCACTACTTCAACTCA--AGGA-CCTCGATTCTGGACCAAGCCT---TGCCCAATGATGTTCTCTTCAGCAGCACATATCCCTCTCTCCCAAAGAGTCTACCTTTGCGGAGGCCATCATAT-ACATTAGGAATGAAGTCATTACATGNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAGTAAAGATTCTTCTCC---TACACTGGCTTCTAAAGTGGACCAATTGGAAGGGATGCTGAAGATGCTCCAGGAAGA-TTTGAAG---AAG--------------------------------------------------------------------------------------------------------------------------------------------------------------");
    // refs.push_back("ATGACCAGCTTGAAGCGGTCCCAGACTGAAAGACCTGTTGCAACTGAAAGAGCATCAATGGTTGGAACAGATGGGATTT---CCAAAGTTCACACTGATGACTTCTACATGCGCCGCTTCAGGTCTCA-GAATGGCAGTTTGGGTTCTTCTGTCATGGCTCCTGTGGGCCCTCCTCGAAGTGAAGGTTCTCACCATATTACCTCCACTCCCGGAGTCCCAAAGATGGGGGTTAGGGCA------AGGATTGCAGATTGGCCACC-AAGGAAGGACAAT-----GTTAAAGAAT---CA---AG-CCGA-TCAAGTCAAGAAATGGAGACGTCAAGTTGCCTTGACAGCATGTCCTCCAAAGGCAGTCCTGGAAGTCAAGGTAGTTCAGTTAACCTCAACTCCAGTGACTCTGCGATGTTGAAAAGCATACAGAATAC--TCTAAAAAGCAAAACCAGACAATCAGAGAACATGGATTCCAGATTTCTCATGCCGGAAGCTTATCCTAGCTCCCCTAGGAAAGCCCTTCGCAGAATAAGACAACGTA-GTAACAGTGACATTACTATAAGTGAACTTGATGTGGATAGCTTTGACGAATGTATCTCACCCACATACAAA-ACTGGTCCATCACTACACAGGGAGTATGGTAGCACATCATCAATTGATAAGC-AGGGAACTTCAGGAGAAAGCTTTTTTGATTTGTTGAAAGGATACAAAGATGATAAATCTGATCAGAGGGGCCCAACCCCAACAAAACTGAGTGAATTTTTGATAGCCGG-TGGAAGCAAGGGTTCAGGCT-TTTCTTTGGATGTGATAGA---TGGACCTATCACACAGAGAGAGAACCTCAGACTTTTTAAGGAAAGGGAGAAACCAC-TCAAGAGACGGTCTAAGTCAGAAACGGGGGATTCATCCATTTTTCGTAAATTACGTAATGCCAAATGTGAAGGGGAAC---TTGGA--AAA---T-CATCAGACCTTGAAGATAACCGATCAGAGGA---CTCTGTCAGACCTTGGACTTGTCCTAAGTGT--TTTGGTCATTATGATGTCCAGAGCATATTGTTTGACTTGAATGAGGCAGTTATCAATAGAC---ACAATGTAATAAAAAGGAGAAA-CACAACCACAGGAGCATCTGCAGCTGCTGTTGCTTCTTTAGTTTCTGGGCCTTTATCAC-ATTCTGCAAGTTTCAGTTCCCCTATGGGA-AGTACTGAGGACTTGAATTCTAAGGGAAG------CCTCAGCATGGACCAGGGAGATGATAAGAGTAATGACCTTGTAA-TGAGCTGTCCGTATTTTCGAAATGAGATTGGTGGAGAAGGAGAAAGGAAGATCAGTTTATCCAAATCGA--ATTCTGGTTCCTTCAGTA---GTGAAAGTGCCTCATTTGAATCAACTCTTAGCTCTCATTGCACAAATG-CAGGAGTGGCTGTACTTGAAGTGCCCAAGGAAAAT---TTGGTATTGCATTTAGATAGAGTAAAGAGATACATTGTGGAACATGTGGATCTTGGTGCATACTATTACAGGAAATTTTTTTATCAGAAGGAACACTGGAACTATTTTGGGGCAGATGAGAA-CCTGGGTCCAGTGGCTATAAGCATTCGAAGAGAAAAACCAGAAGAATTGAAAG-AAAATGGACCTCCATACAACTATCGAATAATATTCAGAACTAGTGAGCTCCTGACATTA-AGAGGTGCTGTACTAGAAGATGCCA-TCCCTTCAACTGCAAAGCATTGT-ACAGCTAGAGGAT-TGCCCTTGAAAGAAGTACTGGAGCATGTGATTCCTGAGCT-------AAATGTTCAGTGCCTGAGGTTGGC--------------------------------------------------CTTCAACA--------------------------CCCCCAAAGTCA---------------CAGAGC-A-------------------ACTCATGAAACTGGACGAACAAGGGGGGAACTTAAGACACTGTCAGACACTGGGGAAAAGGAAGAGAGTTGGAGGAGTGATGCCAAGGGGTTTTCAAGTGAAGAGAAAATTCCAGGTCTTCAGAGAGATCCAAGAGGAAGAGTTAA-GCTATCAATTGAAAGTGGGAATCATGTATTGCAAATCTGGACAGAGCACCGAAGAAGAGATGTACAACAATGAGTCAGCAGGCCCTGCTTTTGAAGAA--TTCCTCCAACTTTTGGG--AGAACGAGTTCGTCTTAAAGGATTTGAGAAGTA-TC-GT---GCTCAACTTGATACTAAAACTGATTCAACTGGTACCCACTCTCTTTACACAACTTAT---AAAGACTATGAAATTATGTTTCATGTATCTACTTT-GCTGCCATATACACCTAA-CAACAAGCAACAG---------------------------CTCCTACGGA-AGAGGCACATTGGAAATGACATTGTTACAATAG-TTTTCCAAGAACCAGGAGCCCAGCCCTTCAGCCCAAAGAACATTCGCTCTCACTTCCAGCATGTTTTTGTCATTGTCCGAGTTCATAATCCCTGCACTGACAGT-GTCTGTTACAGTGTGGCTGTAACCAGGTCCAGAGACGTGCCATCC--TTTG-GCCCACCAATCCCCA-AAGGAGTCACT-TTCCCTAAGTCCAATGTGTTCAGGGACTTCCTTTTGGCC--AAAGTAATTAACGCGGAGAATGCTGCTCA-TAAGTCT-GAGAAATTCCGTGCCA---TGGCTACCCGGACCCG---TCAAGAATACCTA-AAAGATTTGGCAGAGAAGAATGTAACCAATACTCCCATTGATCCTTCTGGCAAATTCCCCTTCATCTCTCTT-GCCTCCAAGAAAAAAGAGAAATCCAAGCCTTATCCTG--GAGCAGAATTCAGTAGCA-TGGGTGC-CAT-TGTATGGAGTGTCCATGCCAAAGACTACA-ATAAGGCAATGG-AGATAGACTGTCTTCTGGGCATCTCAAATGAAT-TTATAGTTCTCAT-TGAGCAAGAAACCAAAA---GTGTGGTATTTAATTGTTCCTGCAGGGATGTAATAGGATGGACTTCAACTGACACCAGCATCAAAATCTTCTA-TGAGAGAGGGGAATGTGTTTCTGTGGAGAGCTTCAATAACAATAGTGAGGATATCAAAGAGATTGTCAAAAGATTGCAGTTTGTGACTAAAGGATGTGA-GTCAGTGGAGATGACTTTGAGGAGAAATGGGCT-TGGACAGCTTGGCTTCCATGTAAACTAT-GAGGGTATTGTGGC-AGATGTGGAACCTTATGGTTATGCCTGGCA-GGCAGGATTAAGGCAAGGCAGCCGCCTGGTGGAAATTTGTAAAGTG--GCTGT-TGCTACCTTAAGCCATGAACAAATGATT-GACCTTCTGAGAACATCAGTCACAGTGAAAGTGGTCATTATACCACCCCAT--GATGATTGTACTCCAAGAAGG------AGCTGCTCTGAAACATACCGCATGCCAGTGAT-------------GGACTACAAAATGAATGAAGGGGTTTCCTATGAATTCAAGTTCCCTTTCCGAAATAACAACAAATGGCAGAGAAATGCCAACAAAGGGCAAGGACCTCATGGGTCCCAAG-TAC-CCTCCCAGGTTCAG---AGTCCAGTGATTTCTCGAATGACTGCTGGAAAAGGAGATGGAAA-GATGCAG-CCTCCA-GAACGAGCTGCTA-ATAT-TCCTCGAAGCATCTCTAGTGATGGTCGCCCCCTGGAAA---GGAGG--TTGTCTCCTGGCTCTGACATCTATGTGACAGTTTCATCCATTGCTTTAGCGAGATCACAGCAATGTCGTAACTCCCCCAGCAACCTGTCTTCATCTAGTGAGACTGGCTCTGGTGGGGGCACTTACAGACAAAAATCCATGCCTGAAG-GGTTTGGAATGAGCCGAAGATCCCCAGCCTCCATTGACCGGCAGAACACTCAGACAGATAT---T-GGTGGCAGTGGGAAAGCCACACCCAGCTGGCAGAGAAG-TGAGGAGAGCATTGCTGACCAAATG---------------------------------------------------------------------GAGCCAACATGCCACCTCCCAGCTGTATCAAAGGTATTGCCATCTTTCCGAGAGAGCCCTAGTGGGAGATTAATGCGGCAAGATCCAGTTGTTCATTTGTCACCAAACAAACAAGGGCATTCTGAT-AGTCATTACTCCAGCCACTCCAGTAG--CAATACTCTGTCCAGCAATGCATCTA-GTGCT--CACAGTGATGAGAAGTGGTATGA---TGGAGATCGCACAGAGTCTGAACTCAATAGC-TACAACTAT-CTTCAAGGCACCTCTGCAGATAGTGGCATCG-ATACCACCTCCTA-CGGC-CCCAGCCACGGGAGCACAGCAT--CTCTGGGGGCTGCTACATCATCCCCACGCT-CA---GGGCCTGGAAAGGAAAAGGT-GGCCCCTCTGTGGCACAGCTCCAGCGAAGTGGTCT---CCATGGCAGATAGGACCTTAGAAAAAGAGAGTCATGGTGTGGACCGGAAGACTGAATCTTCACTGAGCCTGGATATCCATAGCAAGAGCCAACCCACCTCTAATCCTTTAACAAGGGAGAACAGTACTTTTAGCATTAACGATGCTGCATCTCATACAAGTACCATGAGTTCCCGACACTC---TGCCAGCCCAGT-TGTTTTCACCAGTG-CCAGAAGTTCACCCAAAGAA-GAGCTTCATTCTGCCACTTCCCCGCAGCTCGCGCCTTCTTT------CTCCTCCTCCTCTTCCTCCTCTTCTGGTCCTAGAACTTTCTACCCTCGGCAGGGTGCTACT-AGCAAGTACCTGATTGGATGGAAAAAACCTGAAGGGACAATAAACTCTGTGGGGTTTATGGATACAAGAAAACGCCACCAAAGTGATGGCAATGAAATAGCCCACACCAGGTTGCGTGCCTCTGCCAGGGACCTCCGGGCATCTCCAAAGCAAACATCTAAGTCAACAATTG--AAGAGGATCTAAAAAAACTAATCGACTTTGAAAGCCCAACTCCTGAATCACAGAAAAATTTTA-AG------TTCCATGGACTTTCCTCTCCACAATCTCCTTTCCCTTCCACTCCTACCACAAGGCGGA-CCTTGCACAGAACGTTGTCGGATGAAAGTATTTACAGCAGTCAAAGAGAACACTACTTCCCCTCG--CGGA-ACTCGATTCTGGACCAAGCCT---TGCCCAATGATGTCCTCTTCAGTAGCACCTACCCTTCTCTCCCAAAGTCTCTCCCTTTACGGAGGCCATCGTAC-ACATTAGGAATGAAGTCATTGCACGGA-GAGTTCTCGGCCTCAGACAGCTCCCTCCCTGATATCCAGGAGCCAAGA---AGGCAGCCCATGCCAGACCCAGGTCTCATGCCCCTGCCAGACTCTGCT-TCTGA-CCTGGACTGGTCAAACTTGGTAGATGCTGCCAAAGCTTTTGAG------GTCCAGCGTGCCTCATTCTTTGCTGCTAGTGATGAAAACCACCGACCCCTGAGTGCTGCATCCAACAGCGACCAGCTAGAGGA-CCAAGCTATGGCCCAGATGAAGTCTTACGGAAGCAGTAAAGATTCTTCTCC---TACACTGGCTTCTAAAGTGGACCAATTGGAAGGTATGCTGAAGATGCTTCAAGAAGA-TTTGAAG---AAG-------------GAAAAAGAAGACAAAGCTCACCTTCAAGCTGAAGTCCAGCATCTTCGGGAAGACAACTTGCGTTTACAGGAAGAGTCCCAGAGTGCATCGGA-CAAGCTGAAAAAATTCACTGAATGGGTCTTTAACACAATTGATATGAGCTA-");
    // refs.push_back("A-GTCAAGTATAGCTAACACAACTTTTGGGTGAGCCACTTAAAGAAGTCGAAGAGTCAGCGTAAGACCCTA-AGGA-GGACATTGGCGTTCAACAGGAGGGCGTTTACA-ACGTGAAGTC-GGACTTC-TAACCGAAACAGGAGAAAAAG-------------G-AACCCGGAATTTAGGAGAACGTCGTAGAAGTCGTTTGGAAGTTCAACCAG-GTGAAATCTCCGCTCTCAACCACCTCTGAGAAATG--GAC-ATATTCCGAAATAGACACGTC-----ATTCGACCAG---GAGTCGG-ACTAGTGGTGACCTCCGTCGTGAGTGTCCGGACACCAAAAGTAGGCGTTGTCGTTTCTTACTTC-GAGCGACTTG------AAGTT-TCCGAAACCGTCGTAGATGGTCCAAACTGGTCAG-GTCTAGCCTCC--GTCTCAGTCCGTTACCGTATTCTGGTCCCAGACCG-TACCCGACGGCAGAAACGAGAA-CCTATAGACACTCCCTCGACAGACTCCGTCTTTTGA-GAGGTACATCA-CTAA-AGTAAGGGTTTCAT-ATACTACCAGA-GGCATTACCGTCGCTAAAACCCTCTCTTCC-CATACACG-ACCACTTGTCCTGTA-GTAAACCGTACCGAACCAGGTTCTCACTCCGGGCGCTCCGCTTCTTCACCAGAGAAACTGGTGACATTTATGAGAGTAGACTATTTCAGGAAACGTTCCGTGCGGAACACCATCCTCACGACCCTTTCCCTCTGAC---TGTCCTG-TCTCACACACT------GAA-TTT-TAAGAAAACCCTAAGTCCACAACCTGATAGTTCTAGATAG----TCAAAGAAGTCCAGAAGAAGTTGACAC--C-TGAAGCTACGAGCAAAACCACTACGAGCGTCTAGAGAGCGACTCCGTGCGT--CG-GCTCCCACCCG-GTAGAGTAACGGTAGTGA--GACCACTGCAAAAGAACATAGGTATTTGGGGTGACTCAAATAGCAGGGAAGCCCGAAGAAGGTAGGCTAGTCT--ATGAACGATCATCGTGGGACCGCACCCATCTTTCAGGATCCTGGTCTACTCCTCCTTCTCCTCCTCCTC------TTTCTTCCACG-CTCGACC-CCTCTTCATCGTCTTACTTCGAGAAGAAATCCACTCGAAGACCGTGACCACTTTTGTTGGCCCGACCGTCTCACAGCGCTCGAGTACCATGAACATACACTCTGCCGTAGTGATTCTGATATTCGT-GACAAGAGGGAGCAATCTCCTAACCTCCGACCAACTGAGAACGACA-CTTATAGGTCTGAG-TCACTTCTGAGTCGAAAAGCCAGAT-AC---ACCGAGAGAAA--A----AGAT--TACAGGA---CAGACGTTACCTCTGGTGAAGTGAC-CTCGACAC-GGTGTCACCACAG-TGAAAAAGGAAGCAATCG-GG---GCAA-GCTCCT---CT---ACAACGACGAGGGTCTCTTCGTCACGACGGTACCGATCCAG-GTATCCTTCACCATA--GATACGGCGATAGTCGGCTTCAG--GG--GACTTCTATCAACATCGACAAGTCGAG-A-CGAAGTCACGCGAGAGGTGACAGTATG--GTGAAGAGTAGTGACACTCGTGAC-CTCCGTAACGACCTCTCTCATAACGACGACCTTACCGACCTCATCACCGACAGTC-TTACTGGGACG----AATA-AACCTCTGT-TCACTTGGTGT-CCTAGGACAGA-GTAGTCAGAAGGT-GACCC-TGAGAG-AG-ACTTCCTACCGTCAT--GAAAACTATGACGACCCTCTACTGTGC---AACCAAG---------------------------------------------------------------GTAGACCAGTCGCTATGATAGGAGTGAAGAAACGGTCGACCCACACCTAAACGGCGATGGTGGG---CATAGACTAACGCACGAGAC-GGAC-AGTTACCTTCGTCCCCTAGATGCCGATTAAGGTTTGGG--AAGACCGTATCTAAAAACAGACATACACGAGGGTCGTCTCGGACAGAGTGAACTCCAACTATTCAACGAACCTCTTAACGCTGTGACAACACTGGCACGATTACGTTACCTGCTGTGGCAGTGCATCTGTAGTCTTGGTCCTCTGTTGGA-GGCCGAGAGGTTGCCCGCAGGTAGTGATCTCTACGAAGCCCCCTACAACCGACGAGCAAGACCCCCACCGTAAAAGGGAAGAGGAAAGGGACTCGGTTAGGATCTACCGTGACCTGAAACTTGAACCCTACCAT-GGACTCGGTGTACTCCAGGGACAGGAAACAAACGTAAAGCGACGGTAAACAATAATAAAG-ACTTCCCTTTAAACTTAAGTATGCTGTAAGGAAGTAAGTGAAACATAAGG-------------TAGTGACCGTATGCCATCAAAAGACTTGTTGA------AGACGCACCACACGTCAGGAGTATACCTCCCTGTTACTGCTGGAAGTGGCAATGACTACAGGAGTCTTCTAGTTAGTAAACGAGTACCGAGTCTCACCGATGTCGGTGGAATGTTTAAAGGTGGTCGGACGATGGGACGGCGTCAGGACGGAC--GGTGCGCATGGGCATCCCGAGTTGTAGACGGTGTTAAGGAAGTATCAAATGTACCTTCGGTTCGACTGGGTC---TGGTAAGGATGCGTCTCAGTAGAGGTGACAAAGTGTAGGGAATCAGTGTTCGAGGTCGGAAAACTGCTAGAG-AAACTACAGTAGTTACAACGACTACTTCGAGAGGTTTCTTTGCGTAAGTGGAGAGAGTATCTTATAAAAATACAAACACAGGCAACTTCAGGTTGGATAATGTAGAGCTGTCCTCGTTAATTTGTGGTGTGAAAAACACAGGACAAGTTACTCCTGATGTTTGAGTAATCTCTGTGGGTCCTCTGTCAGATAAAGGTATGG--GAACGACATCAGAAAGCGTACCTGTGAGGTGTGTTATCGTGGTCTTGATATCTAGAGACGAGGACC-TATTCCGAACCTGAAAAGAAAGAAAAAACTTCGTTCACTCTACTTCCCCTTAAACGGTCTTCCCAGTTATCCACACAACCACTGTAAGAAGAGACGGTTCAGAAAGTTCATAAGAACTGCCCAG-GACCACCGGTAACGTGCTTTGAAAAGACTAAATACTCGTCGTAAGAGACGTAATTACTGAAACCGATCATTCTTCAGGGACT-TGTGCAAACT-AAATCCTTTACACTGCG-GAAAACCTTATCCA-CCAGGTTTCCTACCATGTAGAGACCTAGAACAGTGTCGTTGTGACATTGTTTGTGACAGCCATGTCCCTGGCACGTG--AGAGT-GTTACTGTTTCTGCACGACCTTCACTCT-AGCCTACAAGAAACCCGACTTAC-CAACACGAGGTCCGAGGACCTTTTGATAACAGTGTTATAGTAAAGGTTATACAGCGAAAGAGTCT-TC---------------------------GACAACAAACAACAAACCTCATATACCGTTGTATCATCTTTGTACCTTGTACTAAAGTATCAGGAACATACAACATATGTCTCTCACTCATGGACAACTTAGTCAGAAGCACAGGTCGACTCGTGCTATGAAAAGTTTAGGAAACTCGGCTTGAGCAAGAGGGTTCTCGACTTCTTTGAGAA-GTTTCCGACCTGGACGGCTGAGCAATAACATGTAAAGGAGAAGACACGAAACGGGTCGAAACGTCATGTACTACGGGTGGAATTCAACCATTGAATTG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCGAAGTATTCG-------------------ACAAGAC---------------ACTGAAATCCTC---------------------------ACAACTTC---------------------------------------------------CGGTCGGAGTCCGTTACTTATAAG-------TCGAGTCCTTGGTGTATGAGTTCGTGAAGAAAATTCCCTTCAGGAGACCGACACGTCACGAAACGTCAACTTCCTTATCGTAGAAGGTCGTGTCTCGGAGA-ATTACAGTATTCGAGTGACCAAGACTTTTACTAGGCCATCAACATGACTCCAGGCAAAAGAAACTAAAGGAGACCAAAAAGAGAAGAATACGAGTGTCGATGACCTGGTTCCAAGAGTAGACGGGGTTTTATCAAGGTCACAAGAAAGACTATCTTTTTAAAAGC-CATTATTATACGTGGTTCCAGATGTACAAGTTGTTACATAGAAAAGTGAGACAGATCTACATTTTAGTTTAAGAGGAACCCTTGAAGATCTTGACGGTGAGGTCGTAAGCATGTTACTCTTGAGTCTCAACTGAGTTTGCTCCGTGAAAGTGTAGGTGACTTCCTTGGACTTAATCTAAACCTTTCCGAATAGAAGGAGAGGGGAAGGGGTGGTTAGAGTAATGCTTTTATTCCCGTC-GAGTAGTGTTCAAGCAACGAAAATAGTAGGGGGACCAGGTATAACTCT------GAAGGAAAACTTAAGTCTAGAAGTCACGACGGGTAACCTCTTAACTTTGAACGTCTTACCCTATTCCCAGGTCTCTGCTCACTGCGATGTCGTCGTCGTCTACGAGGACAACAACACAAGGAAGAAAATTATTGTAATA---CGGACAACTAACGACGAAGTAAATTTAGTTTATTTTATGAGACATGTAGTATTACTCGTTTTGTGAAACCTGTTCAGGTTCCGAAC-TGCGAT---AGAAGACTGGACAATAGAAGTTCTAGACTTCT---AAAGGGTT---CGAGGGGAAGTGGAAACCGTAATGCGTCGAATGCTTTTTACCTA-CTTAGAGGACAAAGACTAAAACTCGCAGAGAACTCACCGAAGAGAGAAAG-GAATTTCTCGGAGTCCAAGAG--GGAAACTCATTATCCAGGT---AGCTAGTGTAGATCTCTCTTTGGACTTGGGAACGAAGGTGACCGTTAGTCTTCAAGTGACTCGAATCAACACGGAC-CCGAAGAGA-CTAGTCTGAAAAGTAGAAATATCGGGAAATCATTCAGTTTTTTTGAAAGGGGACTTCAAGGGACGAATAGATAACTTCTACACGACGGTATGAGAGATACGTCTCTACCAGGTCTAAATTTGCACCCACTTTAT-GTAAGTAGTTTTGATAGGTGTAGTTC-GAGTGAATACCATTACAGCGACAACGATGCGACAGAATAAGATGCTTCTCGAAAGGAACCCCTTGACCCTATCGGTAGTCCGCATTCTTTGGACCTCAGATCTAAGAGTCTAACAGACCAGAACGAAAATTCACACAAGACCTATGAGAAATTGTACTGTCTTAGTGATCGCAAGTCCAATTGCCTCGAAGGGACTGAGGGTCCTGACGAAAATCTTCTGTACGAC-AGTTCTGTTGACCGTCAAAGTTATAGAACTGA--TCTGGACGA---ACAAAGGAACTACAAGAGGAAGGATCCCCCGGTTAGTCGTTAAGA------ACGGGACTGGGGGTAAAACCCCTGA---GGTCCTCATCTCCATTGCACTACCCTTGAAAGTAACGCTCCTCCAGGTTGACCTCGG-TACTGACGTCTAG----GATCCGACGGTAAGACACTGGATTTCGCGGAGTACATCTTCAGT-AGGCTCACCTGAAAAC---ACCACGGTAGACAAG---GTTGACTGTGGG---GTCACCGTTGTCGGGAAAGCCAGACCCTGGCGAAGTTCGACCAGTA-----------");

    // querys.push_back("-GATCGAGTACAGTTACCACAACTTCTGGGTAAGGCACTTGAAGAAGTCGAACAGACTCCGTAAGACCCTGAGGAGGAC-GTCGGAGTCCAACAGGAGAGCGTCCACGACGTGAAGGCGGACTTCCACGCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGGGCTTCGTAGAAGTCGTATGGAAGGTCGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTCAGAAATGACGACGACATTCCGAAGTAGACCCGGTCTCGGACAAGGAGGTCGACTAGTGACAACCTACGTCGTGAGTTCCCCGCTACCAAAAGTAGTGATC--GTCGTTTTTTTCTCCGAGAGACCTG------GAGTATCCGAAACCGTCGTA--GATGGTCCAACCTGGTTAGGTTCAGACGTCGTCGTAGTCCGTCCCCGTAGTCCGGTCC-CAGGCCGTATCCGACGGACGCGCAGAGGACCTACAGCCACTCCCTCG-ACAGGCTCCGGCTCTTGAGAGGTACGTTGCTGAAGTAAGGATTCCA-TATTCTTCCGGAAGAGTTACCCTCTCTGAACCCCTCTCTTCCCATGCACGATGACT-TCTC-CTGCAGCAACCCGTTCCGAACCAGGTCTTCACTTCGGGAGCTCCACTTTTTCA-CGAGGGAGACCGGTGACATTTACGAGAGCAGGCTGTCACAAGACACGTTCCGGGCGGCGCTCCACCCCCAC--GACCCCTTTCCTCTGACTCCTCTCCT-CTCACGCACCTT------GAATTTTAAGAAGACATTTA-GTCCTCAACCCGAAAGTTCCAGCTAGTCAAAGAAATCTAGAAGAAGTTACCACCTGAACCTCCA-ACCGAATCCTCTACGGGCCTCCAGAGAC---CAAC-TCCGTGCGTCGGACCATACCCGATAAA-GCAACGGTAGCGA-GACCACTGCGAAAGAGCACAGGTATTTAGGGTGTC-TCAAATACCAAGGA-AGTC--CAAAAAAGGTAGGTTAATCCATGAAC-GATCATCGCGGGACCGCTC--CCATTTTTGAAGATCCTGGTCTCCTCCTCCTCCT-CCTCCTCCTC------TTCCTGCCACGCTCGACGCTCCTCCGCCGACCTACTTCGAGAAGAAATCCACTTGAAGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGAACACACCCTTCGTCGTAGTAAATACGACTTCCACGATAAGAGGGAACAGTCTCCCGAACTCTGTCG-AACCGAGAACCACACATACAGATTCGA---GTCCCTTCTGAGACAAAAGGCCAGGTCCGGCACC-GAGAGACTGAGGTTTCAGGCTAGACGGTAC-CT--CTAATGAAGTGACCTCGACACGGTGTCCCCACCGTGAAAGAGGAACGGACCGGG---ACTCGCTCCACTGCTACACCGTCGGGGGTCCCTCCGACACGACGGCACCGACCCCGGTATTCTCCACCA----------------CAGTTACGGTGACAGTCGTCTCCACGGAACGTCTATCAATATCGACAACTCAAGTCTAAGGCACGCCAGGGGT---AGTATGGTGAAG-AGTAGT---GATAC-CCGTGAGCTGCGTAACGACCTCT-CTCATAACGATGACCTCACCGAACTCATCACCGATAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGTAATTAGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTACAACC----GAGACATATAA-GTACGACGAGC-TCCTGTTTTCTTAACTTTAGGACTCCAGGA-GATATTGACATTCGGTAAACCAGTCGTTACGATAGGAGTGAAGAAACGGTCGATCCATAACTAAAAGGCGACGGTGGC---TATAGACTGACCCAC-----AAGACGGACAGTTACCTCCGGCCCCT-----AGACGCCGACTGAGGTTTGGGAAGCCCGTACCTGAAGACGGACATCCACGGGGGGTGTCTTGGCCATAGTGACCTACTTCTGTCCAACGATCCTCTCAAGGCTGT---GACCCTAGAACGATTTCGGTACCTACTCTGGCAGTGTATCTACAGACTTGGTCCTCT-GTCGGCGG---AGAGATCACCCGCGGGTAGTGACCTCTACGAAGCTCCCTACAACCGTCGAGA-A-AGACCTCCTCCGTAGAAGGGTAGAGGAAAAGGTCGTAAGTCGGCACTCCAGTACCCTGAG-ACGTGGACCCTACCCT-GAACTCCACT-TACTCC------GGGGAACGACCGCAAGGAGAC-GGTAAACAATA--ATAAAGCCTTCCCTTTGAACTTAAGCATACTTTGTGGAAGTAAGTAAAACATGAGG-------------TATTGACCGTAC-GCCATCCAAAGCCTCGTCGA------GGAGGCACCCCACGTCAGTAGTACCCCTCCTTACTACTGTTGGAAGTG--GCACTGTCTACAAGAGTCCTCTAGCTAGTAGACGAGTACCGATTCTCACCGATGTCGGTGGAACGTCTAGAGGT-GGTCCG-CTGACGGGACGGAGTCGGGAC-GGACGGTCCGTATTGGCATCCCGAGGTGTAG-GCGG-TGCTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATTGGGTAAAGAAGCGTCTCAGTAGAGGTGACTAA--GTGTCGGAAAACT-TTGTTTGACGTTGGAAAACTGTTAGAGAAACTAGAGGAGC---AACAATTATTTTGAGAGGTGACTTTGTGTAAGAGGAGCAAGTATCTTCTAAAAC-TCCGACCACAGTCACCTTCAGGTGGGATAGT-GTAGAGATGTCCTTGTTAACTTCTGGTGCGAGAAGCAAAGGACAAGTTACTCGTGCTACTTGAGTAACCTCTAGGGATTTTCCGTCAGATCAAGG-TACCGGAATGACA--TCAGAAACCGGGCCTGGTGGGTGTGTTACCGGGG-GTACGACGACTCGAGTCGA-GGACCCATACCGAA-TCTGAAAAGGAAGAAGAACCTTC-GGTCTCTCTACTTACCTTTGAACGGTCTTCCC-AGCTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGCTCAACGGTACCGGGCCTTGAAGAGGCTAAATACTCGTCG-TAAAAGACGCAAATAGTGAAACCGGTTTTCCTTCAGGGAC--TTGTGTAA--ACTAAATCCCTTTCACT--GGGGAAATCCTTACCCTC-CGGGTTTCCTTCCGTGTAGAGACCTGGACCATTGT--CGG-TGTGATATTGTCTGTGACAGTCTCGTCCCTAACAC-GTGGGACTGCTACTGCTTTTGCACGACCTTCACCCTAGCC-TACAAAAAACCCGACTTGCCGACA----CGAGGTCCGAGAACCTTTTGCTAACAATGTTATAGTAAAGGTTACACGGCGA-AGGAATCCT--C---------------------------GACAACGAACAACAACCCACACATACCGTCGTA----CCATCTTTGTACCTTGTATTAAAGTATTAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGTCAAA-ACCATAGTTCGACACGCGCTATGAAGAGTTTAGGAAACTCGGCTTGAGCGAGAGGGTTATCAACTTCCTTAAGAAGTTTCCGACCCGGTCGACTGAGTAACAACATGTAGAGAAGAAGTCACGAGACAGGTCGAAACGTCATGTACTACGGATGAAAGACGACCATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTACTCG-------------------ACGAGAC---------------ACTGGAACCCGC--------------------------ACAAC-TTC--------------------------------------------------CGGTTAGCGTCCGTGA-CCT------GTAAC-------TCGAGTCCTTGGTGCACGAGGTCGTGAAGAAACTCTCCGT---CGGGAGACCGGCAGCTCACGAACC--GACAGCTGCCTTAC-CGTAGGAGGTCCTGACTTGGAG--A-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGACATCAACATAC-CTC-TAGGTAAAAGAAAG-TAAAGAAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTTACAAGGAAGACCATCTTTTTAAAAGATATTATCATACGTGGGTC-TAG-GTGCACA--AGGTGCTACATAGAAAAGTGAGATAGATCCACGTTGTGGTTCAAGAGGAACCCGTGAAGTT-CATGACGGTGAGGACGTAAACACGTTACCCTTGATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGTCTTAAACT-----AAAACTTTCCGACTAAAAGGAGAGGGGAAGAGGTGGATAGAGTAAGGCTTTTATGCCTGTCGAGTAATGTTCGAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---AAGGCACGACGGGT-AACCCCTCGATTTTGACCGACTTACTCTGTCTCC-CGGTCTCTGGTTCCTACGGTGTCGTCGACGTCTTCGAGGTCACCACCACAAAGAAGAGAATTAGTG---TAACA----CGGACAAGTATTAGCGGAGTAAGTTTAGT--TTATTATATGAG-ACCTGTA-GTAT-CACCCGTTTCGTGAAACCTGTGCAGGTCCCGGACTG-TCTC---AGAAGACT-AGCCAATAGAAGCTCTAGACT-ACT---GAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--ATGCTTTTTATCTACTCAGAGGTC-AAAGGCTAAAACTTGC-AGCGAACTCACCAAAAA--GGGAAAGGAATTTTTCGGACTCCAAGAGA--------G-AGACACTCTATCCGGGC------A-GATATTGCAGGTTTCTCTTTGGTCTTG----GGAACGGTGGTGGTCATTACTCTTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CTAGTCTAAACAGTAGAAACATCGGGAAATCGTTTAGTTTTTTCG-AA-AGAGG-CC---TGCAAGGGACAAATAGCTAACTTCTACACGATGGTATAAGGGACACGTCCCTACCAGGTCAGAACATACATCCACTCTATGTAAGTAGTTTCGATAGGTGTAGTT-CAAGTGAATACCACTATAGTGACAACGAAGCGACGGCATAAGACGCTTCTCGAAAAGACCCC-CTTGACCCCATCCGAAGTCCGTACTCCTTAGACCTCAGGTACAAGAGGCTGCCAGAACAGAACAAAAAGTC--GCACAAGACATACGAAAAGTCG--TACCG--ACTCAGTAACCTCAACTCCGATTG-TCTTG---A--AGGGACTGAGTGTCCTG-AC--GAAAACCTCCTGTCC-GATAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATACAAAAGGAAAGAACCCCCGGTT-AGACGTTAGGA------ACGGGAATGGG----GGTAGAAACCCTGAGGTCCCCAACTCCAATATACCACTCTTGGAAGTGAAGCCCCCCCAGGATGTCCTCGGTATTGACTACTAGGATTCGACGGTAA---AACCCTGGCCTTTGCGGCGTACATCTTTAGTAGTCACACCTGAAATC---CCCACGGTA-GACACGGTTGTTGTCTCCAGGACAGTCACC--GTTCTCCGGAAAGACAGACACTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-AAAAAA----AATGACCAGCTTGAAACGGTCACAGACAGAAAGGCCCCTTGCCA-CTGACAGGGCCTCTGTTGTTGGC-ACAGACGGCACCC---CCAAAGTCCACACTGATGATTTCTACATGCGGCGCTTCCGGTCCCAAAATGGCAGCTTAG---GATCATCAGTTATGGCTCCTGTAGGACCCCCCCGAAGTGAAGGTTCTCACCATATAACCTCAACCCCTGGAGTCCCAAAGATGGGGGTAAGAGCA------AGGATTGCAGATTGGCCCCCAAGAAAG-GAAAACATAAAAGAATCT---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGATAGCCTGTCCTCCAAAAGCAGTCCTGTGAGTCAGGGAAGTTCTGTTAGCCTCAATTCCAATGACTCAGCCATGCTGAAAAGCATACAGAACACGCTGAAAAACAAGACAAGACCGTCGGAGAACATGGACTCCAGATTTCTCATGCCTG-AAGCCTACCCCAGCTCCCCCAGAAAAGCTCTTCGCAGAATACGGCAGCGAAGCAACAGTGATATCACCATAAGTGAACTTGA-TGTGGATAGCTTTGATGAATGTATCTCACCTACATACAAGACTGGACCATCACTGCACAGGGAATATGGTAGCACATCTTCAATTGATAAACAGGGAACGTCTG-GAGAAAGCTTTTTTGATTTGT-TAAAGGGCTACAAAGATGACAAATCTGATA---GAGGTCCAACTCCAACCAAGCTCAGTGACTTTCTCATTACTGGTGGTAGCAAGGGTTCTGGTTTCTCTTTGGATGTTATAGA----CGGGCCTATCTCACAGA-GAGAGAACCTCAGGCTTTTTAAGGAAAGGGAAAAACCACTCAAGCGACGTTCAAAATCTGAAACTGGAGACTCATCTATTTTTCGTAAATTGCGCAATGCCAAAGGTGAA---GAAC---TTGGGAAG---TCATCAGATCTTGAAGATAACC--GATCAGAAGA---CTCTGTCAGGCCCTGGACATGTCCAAAGTGCTTTGCCCACTATGATGTCCAGAGTATATTATTTGATTT-GAATGAGGCAATTATGAACAGGC---ACAATGTTATTAAGAGGAGAAACACCACCACTGGAGCTTCCGCAGCTGCC-GTGGCATCCTT-GGTCTCTGGACCTCTG-TCTCATTCAGCCAGTTTTAGCTCCCCAATGGGCAGCACAGAGGACCTGAATTCCAAAGGAAG------CCTCAGCATGGACCAGGGAG-ATGATAAAAGCAATGAGCTTGTAATGAGCTGTCCATATTTTCGGAATGAGATAGGTGGAGAAGGGGAGAGGAAAATCAGCCTTTCAAAATCAAATTCTGGCTCCTTTAGTGGATGTGAAAGTGCCTCCTTTGAGTCTACCCTTAGTTCCCATTGCACAAATGCAGGAGTGGCAGTACTTGAAGTGCCCAAGGAGAACTTGGTGTTGCACCTAGATAGAGTGAAA-A----------------GATACATCGTGGAACACGTAGATCTGGGT-GCATACTATTATAGAAAATTTTTCTACCAGAAGGAACACTGGAACTATTTTGGGGCTG-ATGAGA---ATCTT-GGTCCAGTGGCTGTGAGCATTCG-AAGGGAAAAACCAGATGAAATGAAAGAAAATGGATCTCCATACAACTACCGGATAA--TTTTTAGAACTAGTGAGCT-CATGACACTG-AGAGGTT--CGGT--CCTGGAGGATGCCATTCCGTCAAC-AGCCAAGCACT-CAACGGCCAGA-GGCCTGCCTCT----CAAAGAAGTGCTGGAGCACGTGGTTCCCGAGCT-------CAATGTCCAGTGCCTGCGGTTGGC--------------------------------------------------CTTCAACA--------------------------CGC-----CCAA-GGTCA---------------C-----AGA-GCA-------------------GCT-CATGAAACTGGATGAACAAGG----------------------------------------------------------------------------------------------------------------------GCTGAACTACCAGCAGAAAGTAGGCATCATGTACTGCAAAGCTGGACAGAGCACTGAAGA-AGAGATGTACAACAACGAGTCAGCTGGCCCAGCCTTTGAAGAATTCCTCCAACTATTGGGAG-AGCGAGTTCGGCTCAAAGGATTTGAGA-AGTATCGAGCACAGCTTGATACCAAAACT-GAC---TCCACTGGA--ACCCATTCTCTGTACACAACATA--CAAAGATTATGAAATTATGTTCCATGTTTCTACCATGCTGCCA-TACACACCCAAC-AACAAGCAACAG---------------------------CTCCTGAGGAAGCGGC---ACA----TTGGAAATG--ATATTGTA-ACAATTGTTTTCCAAGAGCCTGGAGCACAGCCATTCAGCCCAAAAAATATCCGATCCC-ACTTCC-AGCATG-TTTTCGTCATCGTTAGGGTGC-ACAATCCCTGCTCTGACAGTGTCTGTTATAGTGTGGCTGTCACCAGGTCCAGAGATGTGCCTTCCTTTGGACCTCCCATTCCTAA-AGGGGTCACTTTCCCTAAGTCAAA-TGTGTTCAGGGACT-TCCTTTTGGCGAAAGTGATTAATGCAGAAAATGCTGCTCATAAATCGGAGAAGTTCCGGGCC-ATGGCAACTCGGACCCGCCAGGAATACCTGAAAGATCTGGCAGAAAAGAATGTCACCAACACCCCT-ATCGACCCTTC-TGGCAAGTTTCCATTCATCTCTCTGGCTTCCAAGAAGAAGGAAAAGTCTAAGCCATATCCAGGAGCCGAGCTCAGCAGCATGG-GGGCC--ATTGTG--TGGGCAGTCCGGGCTGAAGACTACAACAAGGCCATGGAGCTAGACTGCCTTTTAGGGATCTCCAATGAG-TTCATTGTGCTCATTGAACAGGA-AACAAAGAGCGTGATCTTCAATTGTTCCTGTAGAGATGTGATAGGGTGGACTTCAACTGACACCAGCCTCAAAA-TCTTCTATGAACGAGGAGAATGTGTTTCAGTGGGGAGTTTTATTAACAA---CGAGGAGATCAAAGAAATTGTCAAAAGGTTGCAGTTTGTTTCCAAAGG--CTGTG-AA--TCGGTGGAGATGACTCT--GCGAAGAAATGGGCTAGGACAGCTTGGCTTCCATGTCAACTATGAGGGCATCGTGGCAGATGTGGAGCCCTACGGTTATGCCTGGCAGGCAGGGCTGAGGCAGGGCAGTCGCCTGGTGGAGATCTGCAAGGTGGCGGTTGCCACTTTGAGCCATGAGCAGATGATCGACCTCCTGAGAACATCTGTCACGGTGAAGGTTGTCATCATTCCCCCGCATGATGACTGCACCCCGCGGAGG------AGTTGCTCTGAAACCTACCGCATGCCAGTGAT-------------G----GAGTACAAAATGAATGAAGGTGTTTCATACGAATTCAAGTTTCCCTTCCGAAATAATAACAAATGGCAGAGGAACACCAGCAAGGGG------CCTCATTCACCTCAAG-TCCCGTCCCAGATGCAGAGTCCCATGACCTCGCGGCTGAATGCTGGAAAAGGAGATGGGAAGATGCCTCCTCCAGA-AAGAGCTGCCAACATCCCTCGAAGCATCTCCAGTGACGGGCGCCCACTAGAGA---GGCGGCTGTCTCCTGGTTCGGACATCTATGTGACGGTCTCATCCATGGCTTTAGCAAGATCCCAG---TGTCGGAACTCCCCTAGCAACCTGTCTTCATCCAGTGACACTGGTTCTGTGGGGGGCACTTACAGGCAGAAGTCCATGCCTGAAGGGTTTGGAGTGAGCCGTAGATCCCCGGCCTCCATTGACAGACAGAACACCCAGTCAGATAT---CGGTGGCAGCGGA-AAATCCACACCTAGCTGGCAAAGAAGTGAGGATAGCATTGCTGACCAGATGGCTTACAGTTATAGAGGA-CCT----CAGGATTTCAATTCTTTTGTCCTCGAGCAGCATGAATATACAGAGCCAACATG----CCATCTCCCAGCAG-TATCAAAGG--TACTGCCAGCTTTC-CGAGAGAGCCCCAGTGGGAGATTAATGCGGCAGGA-TCCAGTGGTTCATTTGTCTCCAAACAAACAAGGGCATTCTGATAGCCACTACTCGA-GCCACTCCAGTAGCAATACTCTCTCCAGCAATGCGTCGAGTGCCCACAGTG-----ATGAGAAGTGGTACGA---TGGGGACCGCACAGAATCCGAACTCAACAGCTATAACTATCTGCAAGGC-ACCTCTGCTGA--CAGTGGCATTGACACCACCTCTTATGGCCCCAGCCATGGCAGCACAGCCTCCCTGGGGGCTGCCACATCATCACCTCGCTCA---GGGCCAGGCAAGGAGAAAGTGGCACCCCTATG--GCA--CAGCTCCAGTG--AAGTAATCTCCATGGCAGA-----TCGGACTTTGGAG--ACAG-AGAGCC-ACGGCCTGGACCGGAAAGCAGAGTCTTCCCTGA-GCTTAGACATACAC-AGC-AAGAGCCAAGCCGGCTCGAGCCCTCTGACAAGGGAGAACAGCACCTTCAGTATAA---ACGATGCTGCTTCC-CACACAAGTACCATGAGCTCCCGACACTCTG-CC-AGCCCGGTGGTTTTCACCAGTGCCCGAAGTTCACCTAAAGAAGAGCTTCATCCAGCCGCCCCCTCA---CAGCTC-GCACCGTCCTT------CTCCTCCTCCTCCTCCTCCTCCTCTGGT-CCTAGGA-GTTT-TTACCCT--CGCCAGGGCGCTACTAGCAAGTACCTA-ATTGGATGGAAAAA--ACCTGAAGGAACCA-TAAACT-CCGTGGGATTTATGGA-CACG-AGAAAGCGTCAT-CAG-AG--TGATGGCA--ATGAAATAGCCCACACCAG-GCTG-CGTGCCTCAACC-----AGAGACCTCCGGGCATC--TCCTAAGCCAACCTCCA-AGTCCACCATT--------GAAGAAGATCTAAAGAAACT---AATCGACCTTGAAAGCCCAACTCCTGAATCA----CAGAAGA-GTTTTAAG------TTCCACGC----ACTCTCCTCT-CCTCAGTCTCCTTTCCCCAGCACCCCCACCTCACGGCGGGCCTTGCACAGAACACTGTCGG-ACGAGAGCATT---TACAGTAGCCAGAGGGAGCACTTTTTCACCTCCAGG---GCATCACTTCTGGACCAAGCCCTGCCCAACGATGTC-CTCT-TCAGTAGCACGTACCCTTCTCTCCCCAAGTCACTCCCATTGAGGAGGCCTTCTTACACCTTAGGAATGAAATCGCTGCATGGAGAGTTCTCGGCCTCCGAC-AGCTCCCT-CACTGACATCCAGGAGACCCGCAGGCAGCCTATGCCCGACCCTG-GCC--TGATGCCCCTGCCTGACACTGC--TGCGG--ACTTGGATTGGTCCAACCTGGTAGATGCTGCCAA--AGCCTACGAG------G-TCCAGAGAGCCTCGTTTTTTGCTGCTAGTGATGAAAACCATCGCCCCTTGAGTGCTGCATCCAACAGTGATCAGCTGGAGGACCAGGCTCTGGCCCAGATG-AAGCCTTTCAGCAGCAGTAAAGATTCCT----C-TCCCACTCTGGCTTCTAAAGTGGACCAGCTGGAAGGTATGCTGAAGATGCTTCGGGAAGATTTGAAG---AAG-------------GAAAAAGAAGACAAAGCCCACCTTCAGGCGGAGGTGCAGCACCTGCGAGAGGACAACCTGAGGCTACAGGAGGAGTCCCAGAACGCCTCGGACAAGCTGAAGAAGTTCACGGAATGGGTCTTCAACACCATAGACATGAGCTAGGG");
    // querys.push_back("--AAAA-------TGACCAGCTTGAAACGGTCACAGACAGAAAGGCCTCTTGCCA-CTGACAGGGCCTCTGTTGTTGGC-ACAGACNGCACCC---CCAAAGTCCACACTGATGATTTCTACATGCGGCGCTTCCGGTCCCAAAATGGCAGCTTAG---GATCATCAGTTATGGCTCCTGTAGGACCCCCCCGAAGTGAAGGTTCTCACCATATAACCTCAACCCCCGGAGTCCCAAAAATGGGGGTAAGGGCA------AGGATTGCAGATTGGCCCCCAAGAAAG-GAAAACATAAAAGAATCT---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGATAGCCTGTCCTCCAAAAGCAGTCCTGTGAGTCAGGGAAGTTCTGTTAGCCTCAATTCCAATGACTCAGCCATGCTGAAAAGCATACAGAACACGCTGAAAAACAAGACAAGACCGTCGGAGAACATGGACTCCAGATTTCTCATGCCTG-AAGCCTACCCCAGCTCCCCCAGAAAAGCTCTTCGCAGAATACGCCAGCGAAGCAACAGTGATATCACCATAAGTGAACTTGA-TGTGGATAGCTTTGATGAATGTATCTCACCTACATACAAGACTGGACCATCACTGCACAGGGAATATGGTAGCACATCTTCAATTGATAAACAGGGAACATCCG-GAGAAAGCTTTTTTGATTTGT-TAAAGGGCTACAAAGATGACAAATCTGATC---GAGGTCCAACTCCAACCAAGCTCAGTGACTTTCTCATTACTGGTGGTGGCAAGGGTTCTGGTTTCTCTTTGGATGTAATAGA----CGGGCCTATCTCACAGA-GAGAGAACCTCAGGCTTTTTAAGGAAAGGGAAAAACCACTCAAGCGACGTTCAAAATCTGAAACTGGAGACTCATCTATTTTTCGTAAATTGCGCAATGCCAAAGGTGAA---GAAC---TTGGGAAG---TCATCAGATCTTGAAGATAACC--GATCAGAAGA---CTCTGTCAGGCCCTGGACATGTCCAAAGTGCTTTGCCCACTATGATGTCCAGAGTATATTATTTGATTT-GAATGAGGCAATTATGAACAGGC---ACAATGTTATTAAGAGGAGAAACACCACCACTGGAGCTTCCGCAGCTGCC-GTGGCATCCTT-GGTCTCTGGACCTCTG-TCTCATTCAGCCAGTTTTAGCTCCCCAATGGGCAGCACAGAGGACCTGAATTCCAAAGGAAG------CCTCAGCATGGACCAGGGAG-ATGATAAAAGCAATGAGCTTGTAATGAGCTGTCCATATTTTCGGAATGAGATAGGTGGAGAAGGGGAGAGGAAAATCAGCCTTTCGAAATCAAATTCTGGCTCCTTTAGTGGATGTGAAAGTGCCTCCTTTGAGTCTACCCTTAGTTCCCATTGCACAAATGCAGGAGTGGCAGNACTTGAAGTGCCCAAGGAGAACTTGGTGTTGCACCTAGATAGAGTGAAA-A----------------GATACATCGTGGAACACGTAGATCTGGGT-GCATACTATTATAGAAAATTTTTCTACCAGAAGGAACACTGGAACTATTTTGGGGCTG-ATGAGA---ATCTT-GGTCCAGTGGCTGTGAGCATTCG-AAGGGAAAAACCAGATGAAATGAAAGAAAATGGATCTCCGTACAACTACCGAATAA--TTTTTAGAACTAGTGAGCT-CATGACACTG-AGAGGTT--CGGT--CCTGGAGGACGCCATTCCGTCGAC-AGCCAAGCACT-CGACAGCCAGA-GGCCTGCCTCT----CAAAGAAGTGCTGGAGCACGTGGTTCCTGAGCT-------CAATGTCCAGTGCCTGCGGTTGGC--------------------------------------------------CTTCAACA--------------------------CAC-----CCAA-GGTCA---------------C-----AGA-GCA-------------------GCT-CATGAAACTGGATGAACAAGG----------------------------------------------------------------------------------------------------------------------GCTGAACTACCAGCAGAAAGTAGGCATCATGTACTGCAAAGCTGGACAGAGCACTGAAGA-AGAGATGTACAACAATGAGTCAGCTGGCCCAGCCTTTGAAGAATTCCTTCAACTATTGGGAG-AGCGAGTTCGGCTCAAAGGATTTGAGA-AGTATCGAGCACAGCTTGATACCAAAACT-GAC---TCCACTGGA--ACCCATTCTCTGTACACAACATA--CAAAGATTATGAAATTATGTTCCATGTTTCTACCATGCTGCCA-TACACACCCAAC-AACAAGCAACAG---------------------------CTCCTGAGGAAGCGGC---ACA----TTGGAAATG--ATATCGTA-ACAATTGTTTTCCAAGAGCCTGGAGCACAGCCATTCAGCCCAAAAAACATCCGATCCC-ACTTCC-AGCACG-TTTTTGTCATCGTCAGGGTGC-ACAATCCGTGCTCTGACAGTGTCTGTTATAGTGTGGCTGTTACCAGGTCCAGAGATGTGCCTTCCTTTGGGCCTCCCATTCCTAA-AGGGGTCACTTTCCCTAAGTCAAA-TGTGTTCAGGGACT-TCCTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNTTTGTTTCAAAAGG--CTGTG-AA--TCGGTGGAGATGACTCT--GCGAAGAAATGGGCTAGGACAGCTTGGCTTCCATGTCAACTACGAGGGCATCGTGGCGGATGTGGAGCCCTACGGTTATGCCTGGCAGGCAGGGCTGAGGCAGGGCAGTCGCCTGGTGGAGATCTGCAAGGTGGCGGTAGCCACTCTGAGCCATGAGCAGATGATCGACCTCCTGAGAACATCTGTCACGGTGAAGGTTGTCATCATTCCCCCCCATGATGACTGCACCCCGCGGAGN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------------N----NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNN------NNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNGCTGTCTCCTGGTTCGGACATCTATGTGACGGTCTCATCCATGGCTTTAGCAAGATCCCAG---TGTCGGAACTCTCCTAGCAACTTGTCTTCATCCAGTGATACTGGTTCTGTGGGGGGCACTTACAGGCAGAAGTCCATGCCCGAAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTTACAGTTATAGAGGA-CCT----CAGGATTTCAATTCTTTTGTCCTCGAGCAGCATGAATATACAGNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNT-CTGATAGCCACTACTCGA-GCCACTCCAGTAGCAATACTCTCTCCAGCAATGCGTCGAGTGCCCATAGTG-----ATGAGAAGTGGTACGA---TGGGGACCGCACAGAATCCGAACTCAACAGCTATAACTATCTGCAAGGC-ACCTCTGCTGA--CAGTGGCATTGACACCACCTCTTATGGCCCCAGCCACGGCAGCACAGCCTCGCTGGGGGCTGCCACATCGTCACCTCGCTCA---GGGCCAGGCAAGGAGAAGGTGGCACCCCTATG--GCA--CAGCTCCAGTG--AAGTAATCTCCATGGCAGA-----TCGGACTTTGGAG--ACAG-AGAGCC-ACGGCCTGGACCGGAAAACAGAGTCTTCCCTGA-GCTTAGACATACAC-AGC-AAGAGCCAAGCCGGCTCGAGTCCTCTGACAAGGGAGAACAGCACCTTCAGTATAA---ACGATGCTGCTTCC-CACACAAGTACCATGAGCTCCCGACACTCTG-CC-AGCCCAGTGGTTTTCACCAGTGCCCGGAGTTCACCTAAAGAAGAGCTTCATCCAGCTGCCCCCTCA---CAGCTC-GCACCATCCTT------CTCCTCCTCTTCCTCCTCCTCCTCTGGT-CCTAGGA-GTTT-TTACCCT--CGCCAGGGCGCTACTAGCAAGTACCTG-ATTGGATGGAAAAA--ACCCGAAGGAACCA-TAAACT-CCGTGGGATTTATGGA-CACG-AGAAAGCGTCAT-CAG-AG--CGATGGCA--ATGAAATAGCCCACACCAG-GCTG-CGTGCCTCAACC-----AGAGACCTCCGGGCATC--TCCTAAGCCAACCTCCA-AGTCCACCATT--------GAAGAAGATCTAAAGAAACT---AATTGATCTTGAAAGCCCAACTCCTGAATCA----CAGAAGA-GTTTTAAG------TTCCACGC----ACTCTCCTCT-CCTCAGTCTCCTTTCCCCAGCACCCCCACCTCACGGCGGGCCTTGCACAGAACACTGTCGG-ACGAGAGCATT---TACAGTAGCCAGAGGGAGCACTTTTTCACCTCCAGG---GCGTCACTTCTGGACCAAGCCCTGCCCAACGACGTC-CTCT-TCAGTAGCACGTACCCTTCTCTCCCCAAGTCGCTCCCGTTGAGGAGGCCTTCTTACACCTTAGGAATGAAATCGCTGCATGGAGAGTTCTCGGCCTCAGAC-AGCTCCCT-CACTGACATCCAGGAGACCCGCAGGCAGCCTATGCCCGACCCTG-GCC--TGATGCCCCTGCCTGACACTGC--TGCAG--ACTTGGATTGGTCCAACCTGGTAGATGCTGCCAA--AGCCTATGAG------G-TCCAGAGAGCCTCGTTTTTTGCTGCTAGTGATGAAAACCATCGCCCCTTGAGTGCTGCATCCAACAGTGATCAGCTGGAGGACCAGGCTCTGGCCCAGATG-AAGCCTTACAGCAGCAGTAAAGACTCCT----C-TCCCACTCTGGCTTCTAAAGTGGACCAGCTGGAAGGTATGCTGAAGATGCTTCGGGAAGATTTGAAG---AAG-------------GAAAAAGAAGATAAAGCTCACCTTCAGGCGGAGGTGCAGCACCTGCGAGAGGACAACCTGAGGCTACAGGAGGAGTCCCAGAACGCCTCGGACAAGCTGAAGAAGTTCACAGAATGGGTCTTCAACACCATAGACATGAGCTAGGG");
    // querys.push_back("-------------ATACCACAACTTCTGGGTAAGACACTTGAAGAAGTCGAACAGGCTCCGCAAGACCCTGAGGAGGAC-ATCGGAGTCCAACAGGAGAGCGTCCACGACGTGGAGGCGGACTTCCACTCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGGGCTTCGTAGAAGTCGTATGGAAGGTCGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTCAGAAATGACGACGACATTCCGAAGTAGACCCGGTCTCGGACCAGGAGGTCGACTAGTGACAACCTACGTCGTGAGTTCCCCGCTACCAAAAGTAGTGATC--GTCGTTTTTTACTCCGAGAGACCTG------GAGTATCCGAAACCGTCGTA--GATGGTCCAACCTGGTTAGGTTCAGACGTCGTCACAGTCCGTCCCCGTAGTCCGGTCC-CAGCCCGTATCCGACGGACGCCCAGAGGACCTACAGTCACTCCCTCG-ACAGGCTCCGACTCTTGAGAGGTACGTCGCTAAAGTAAGGATTCCA-CATTCTTCCGGAGGAGTTGCCCTCGCTGAACCCCTCTCTTCCCATGCACGATGACT-TCTC-CTGCAGCAACCCGTCCCGAACCAGGTCTTCACTGCGGGACCTCCACTTTTTCA-CGAGGGAGACCGATAACATTTACGAGAGCAGGCTGTCACAAGACACGTTCCGGGCGGCACTCCACCCCCAC--GACCCCTTTCCTCTGACTCCTCTCCT-CTCACGCACCTT------GAATTTTGAGAAGACACTAA-GTCCTCAACCCGAAAGTTCTAGTTAATCAAAGAAATCTAGAAGAAGTTACCACCTGAACCTCCA-ACCGAATCCTCTACGGGCCTCCAGAGAC---CAAC-TCCGTGCGTCGGACCACACCCGATAAA-GTAACGGTAGCGA-GACTACTGCGAAAGAGCACAGGTATTTAGGGTGCC-TCAAATACCAAGGA-AGCC--CAAAAAAGGTAGGTTAGTCCATGAAC-GATCATCGCGGGACCGCTC--CCATTTTTGAGGATCCTGGTCTCCTCCTCCTCCT-TCTCCTCCTC------TTCCTACCACGCTCGACACTCCCCCGTCGACCTACTTCGAGAAGAAATCCACTTGAGGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGAACACACCCTTCGTCGTAGCAAATATGACTTCCACGACAAGAGGGAACAGTCTCCCCAGCTCGGCCG-AACCGAGAACGACACATACAGATTCGAG---TCCCTTCTGAGACAAAAGGCCAGGTCCGGCACCG-AGAGACAGAGGTTTCAGGCTAGACGGTAC-CT--CTAATGAAGTGACCTCGACACGGTATCCCCACGGTGAAAGAGGAACGGACCGGG---ACTCGCTCCACTGCTACACCGTCGGGGGTCGCTCCGACACGACGGCACCGACCCCGGTATTCTCCACCA----------------CAGTTACGGTGACAGTCGTCTCCACGGAACGTCTATCAATATCGACAACTCAAGCCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGTAGT---GATAC-CCGTGAACTGCGTAACGACCTCT-CTCATAACGATGACCTCACCGAGCTCATCACCGATAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGTAATTAGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTACAACC----GAGACATATAA-GTACGACGAGCT-CCTGTTTTCTTAACTTTAGGACTCCAGGA-GATATTGACATTCGGTAGACCAGTCGTTACGATAGGAGTGAAGAAACGGTCGATCCGCACCTAAAAGGCGACGGTGGT---TATAGACTGACCCAC-----AAGACGGACAGTTACCTCCGACCCCT-----AGATGCCGAGTGAGGTTTGGGAAGCCCGTACCTGAAGACGGACATTCACGGGGGGTGTCTTGGTCATAGTGACCTACTTCTGTTCAACGATCCTCTCAAGGCTGT---GACCCTAGAACGATTTCGGTACCTACTCTGGCAGTGTATCTACAGGCTTGGTCCTCT-GTCGGCGG---AGAGATCACCCGCGGGCAGTGACCTCTACGAAGCTCCCTACAACCGCCGAGA-A-AGACCTCCTCCGTAGAAGGGTAGAGGAAAAGGTCGTAAGTCGGCGCTCCAGTACCCTGAG-ACGTGGACCCTGCCCT-GAACTCCACT-TACTCC------GGGGAACGACCGCAAGGAGAC-GGTGAACAATA--ATAAAGCCTTCCCTTTGAACTTAAGCATACTTTGTGGAAGTAAGTAAAACATGAGG-------------TAGTGACCGTAC-GCCATCCAAAGTCTCGTTGA------GGAGGCGCCCCACGTCAGTAGTACGCCCCCTTACTACTGTTGGAAGTG--GCACTGTCTACAAGAGTCCTCCAGCTAGTAGACGAGTACCGAGTCTCACCGATGGCGGTGGAACGTCTAGAGGT-GGTCCG-CTGACGGGACGGAGTCGGGAC-GGACGGTCCGTATTGGCATCCCGAGGTGTAG-GCGG-TGTTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATCGGGTAAAGAAGCGTCTCAGTAGAGGTGGCTAA--GTGTCGGAAAACT-TTGTTTGACGTTGGAAAACTGTTAGAGAAACTAGAGGAGT---TACAATTATTTTGATGGGTGACTTTGTGTAAGAGGAGCAAGTATCTTCTAAAAC-TCCGACCACAGTCAACTTCAGGTGGGATAGT-GTAGAGATGTCCTTGTTAACTTCTGGTGCGAGAAACAAAGGACAAGTTACTCGTGTTACTTGAGTAACCTCTAGGGATTTTCCGTCAGATCAAGG-TACCGGAACAACA--TCAGAAGTCGGGCCTGACGGGTATGTTACCGGGG-GTACGACGACTCGAGCCGAGG-ACCTATACCGAA-TCTGAAAAGGAAGAAGAACCTTC-GGTCTCTCTACTTGCCTTTGAACGGTCTTCCC-AGCTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGCTCAACGGTACCGGGCTTTGAAGAGGCTAAATACTCGTCG-TAAAAGACGTAATTAGTGAAAGCGGTTTTCCTTCAGGGAC--TTGTGTAA--ACTGAATCCCTTTCACT--GGGGAAATCCTTACCCTC-CGGGTTTCCTTCCGTGTAGAGACCTGGACCATTGT--CGG-TGTGATATTGTCTGTGACAGTCTCGTGCCTAACAC-GTGGGACTGCTACTGCTTTTGCACGACCTTCACCCTAGCC-TACAAAAAACCCGACTTACCGACA----CGAGGTCCGAGAACCTTTTGTTAACAATGCTATAGTAAAGGTTACACGGCGA-AGGAGTCCTC-----------------------------GACAACAAACAACAACCCACACATACCGTCGTA----CCATCTTTGTACCTTGTATTAAAGTATTAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGTCAAA-ACCATAGTTCGACACGAGCTATGAAGAGTTTAGGAAACTCGGCTTGAGCGAGAGGGTTATCAACTTCTTTAAGAAGTTTCCGACCCGGTCGACTGAGTAACAACATGTAGAGAAGAAGTCACGAGACAGGTCGAAACGTCATGTACTACGGATGAAAGACGACCATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTACTCG-------------------ACGAGAC---------------ACTGGAACCCAC--------------------------ACAAC-TTC--------------------------------------------------CGGTTGGCGTCCGTGA-CCT------GTAAC-------TCGAGTCCTTGGTGCACGAGGTCGTGAAGAAACTCTCCGT---CCGGAGACCGACAGCTCACGAACC--GACAGCTGCCTTAC-CGCAGGAGGTCCTGGCTTGGAG--A-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAACATGC-CTC-TAGGTAAAAGAAAG-TAAAGTAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTCACAAGGAAGACCATCTTTTTAAAAGATATTATCATACGTGGGTC-TAGA-TGCACA--AGGTGCTACATAGAAAAGTGAGATAGATCCACGTTGTGGTTCAAGAGGAACCCGTGAAGTT-CATGACGGTGAGGACGTAAACACGTTACCCTTGATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGTCTTAAACT-----AAAACTTTCCGACTAAAAGGAGAGGGGAAGAGGTGGATAGAGTAAGGCTTTTATACCTGTCGAGTAATGTTCGAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---GAGACACGACGGGT-AACCCCTCGATTTTGACCGACTTACTCTGTCTCC-AGGTCTCTGGTTCCTACGGTGCCGTCGACGCCTTCGAGGTCACCACCACAAAGAGGAGAATTATTG---TAACA----CGGACAAGTATTAACGGAGTAAGTTTAGT--TTATTATATGAG-ACCTGTA-GTAT-CACCCGTTTCGTGAAACCTGTACAGGTCCCGGACTG-TCTC---AGAAGACT-AGCCAATAGAAGTTCTAGACT-ACT---GAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--ATGCTTTTTATCTACTCAGAGGTC-AAAGTCTAAAACTTGC-AGCGAACTCACCAAAAA--GGGAAAGGAATTTTTCGGACTCCAAGAGA--------GA-GACACTCTATCCGGGC------A-GATAATGTAGGTTTCTCTTTGGTCTTG----GGAACGGTGGTGGTCATTACTCTTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CTAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTTTTCG-AA-AGAGG-TC---TACAAGGGACAAATAGTTAACTTCTACACGATGGTATAAGGGACACGTCACTACCAGGTCAGAACATACATCCACTCTATGTAAGTAGTTTCGATAGGTGTAGTTC-AAGTGAATACCACTATAGTGACAACGAAGCGACCGCATAAGACGCTTCTCGAAAAGACCCC-CTCGACCCCATCCGAAGTCCGTACTCTTTAGACCTCAGGTACAAGAGGCTGCCAGAACAGAACAAAAAGTC--GCACAAGACATACGAAAAGTCG--TACCG--ACTCAGTAACCTTAACTCCGATTG-TCTTG---A--AGGGACTGAGTGTCCTG-AC--GAAAACCTCCTGTCC-GATAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATACAAAAGGAAAGAACCCCCGGTT-AGACGTTAGGA------ACGGGAATGGG----GGTAAAAACCCTGAGGCCCCCAACTCCAATATACCACTCTTGGAAGTGAAGCCCCCCCAGGATGTCCTCGGTATTGACTACTAGGATTCGACGGTAA---AACCCTGGCCTTCGCGGCGTACATCTTTAGTAGTCACACCTGAAACC---CCCACGGCAG-ACACGGTTGTTGTCTCCGGGACAGTCACC--GTTCTCCGGAAAGACAGACACTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-------------ATACCACAACTTCTGGGTAAGGCACTTGAAGAAGTCGAACAGGCTCCGTAAGACCCTGAGGAGGAC-ATCGGAGTCCAACAGGAGAGCGTCCACGACGTGGAGGCGGACTTCCACTCGAAACAGAAGAAAAAG-------------GAA---GAAGTCTAGAAGGGCTTCGTAGAAGTCGTATGGAAGGTCGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTCAGAAATGACGACGACATTCCGAAGTAGACCCGGTCTCGGACCAGGAGGTCGACTAGTGACAACCTACGTCGTGAGTTCCCCGCTACCAAAAGTAGTGATC--GTCGTTTTTTACTCCGAGAGACCTG------GAGTATCCGAAACCGTCGTA--GATGGTCCAACCTGGTTAGGTTCAGACGTCGTCACAGTCCGTCCCCGTAGTCCGGTCC-CAGCCCGTATCCGACGGACGCCCAGAGGACCTACAGTCACTCCCTCG-ACAGGCTCCGGCTCTTGAGAGGTACGTCGCTAAAGTAAGGATTCCA-CATTCTTCCGGAGGAGTTACCCTCGCTGAACCCCTCTCTTCCCATGCACGATGACT-TCTC-CTGCAGCAACCCGTCCCGAACCAGGTCTTCACTGCGGGACCTCCACTTTTTCA-CGAGGGAGACCGATGACATTTACGAGAGCAGGCTGTCACAAGACACGTTCCGGGCGGCACTCCACCCCCAC--GACCCCTTTCCTCTGACTCCTCTCCT-CTCACGCACCTT------GAATTTTGAGAAGACACTAA-GTCCTCAACCCGAAAGTTCTAGTTAATCAAAGAAATCTAGAAGAAGTTACCACCTGAACCTCCA-ACCGAATCCTCTACGGGCCTCCAGAGAC---CAAC-TCCGTGCGTCGGACCACACCCGATAAA-GTAACGGTAGCGA-GACTACTGCGAAAGAGCACAGGTATTTAGGGTGCC-TCAAATACCAAGGA-AGCC--CAAAAAAGGTAGGTTAGTCCATGAAC-GATCATCGCGGGACCGCTC--CCATTTTTGAGGATCCTGGTCTCCTCCTCCTCCT-TCTCCTCCTC------TTCCTACCACGCTCGACACTCCCCCGTCGACCTACTTCGAGAAGAAATCCACTTGAGGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGAACACACCCTTCGTCGTAGCAAATATGACTTCCACGACAAGAGGGAACAGTCTCCCGAGCTCGGCCG-AACCGAGAACGACACATACAGATTCGAG---TCCCTTCTGAGACAAAAGGCCAGGTCCGGCACCG-AGAGACAGAGGTTTCAGGCTAGACGGTAC-CT--CTAATGAAGTGACCTCGACACGGTATCCCCACGGTGAAAGAGGAACGGACTGGG---ACTCGCTCCACTGCTACACCGTCGGGGGTCACTCCGACACGACGGCACCGACCCCGGTATTCTCCACCA----------------CAGTTACGGTGACAGTCGTCTCCACGGAACGTCTATCAATATCGACAACTCAAGCCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGTAGT---GATAC-CCGTGAGCTGCGTAACGACCTCT-CTCATAACGATGACCTCACCGAGCTCATCACCGATAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGTAATTAGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTACAACC----GAGACATATAA-GTACGACGAGCT-CCTGTTTTCTTAACTTTAGGACTCCAGGA-GATATTGACATTCGGTAGACCAGTCGTTACGATAGGAGTGAAGAAACGGTCGATCCGCACCTAAAAGGCGACGGTGGT---TATAGACTGACCCAC-----AAGACGGACAGTTACCTCCGACCCCT-----AGATGCCGAGTGAGGTTTGGGAAGCCCGTACCTGAAGACGGACATTCACGGGGGGTGTCTTGGTCATAGTGACCTACTTCTGTTCAACGATCCTCTCAAGGCTGT---GACCCTAGAACGATTTCGGTACCTACTCTGGCAGTGTATCTACAGGCTTGGTCCTCT-GTCGGCGG---AGAGATCACCCGCGGGCAGTGACCTCTACGAAGCTCCCTACAACCGCCGAGA-A-AGACCTCCTCCGTAGAAGGGTAGAGGAAAAGGTCGTAAGTCGGCGCTCCAGTACCCTGAG-ACGTGGACCCTGCCCT-GAACTCCACT-TACTCC------GGGGAACGACCGCAAGGAGAC-GGTGAACAATA--ATAAAGCCTTCCCTTTGAACTTAAGTATACTTTGTGGAAGTAAGTAAAACATGAGG-------------TAGTGACCGTAC-GCCATCCAAAGTCTCGTTGA------GGAGGCGCCCCACGTCAGTAGTACGCCCCCTTACTACTGTTGGAAGTG--GCACTGTCTACAAGAGTCCTCCAGTTAGTAGACGAGTACCGAGTCTCACCGATGGCGGTGGAACGTCTAGAGAT-GGTCCG-CTGACGGGACGGAGTCGGGAC-GGACGGTCCGTATTGGCATCCCGAGGTGTAG-GCGG-TGCTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATCGGGTAAAGAAGCGTCTCAGTAGAGGTGGCTAA--GTGTCGGAAAACT-TTGTTTGACGTTGGAAAACTGTTAGAGAAACTAGAGGAGT---TACAATTATTTTGAGGGGTGACTTTGTGTAAGAGGAGCAAGTATCTTCTAAAAC-TCCGACCACAGTCAACTTCAGGTGGGATAGT-GTAGAGATGTCCTTGTTAACTTCTGGTGCGAGAAACAAAGGACAAGTTACTCGTGTTACTTGAGTAACCTCTAGGGATTTTCCGTCAGATCAAGG-TACCGGAACAACA--TCAGAAGTCGGGCCTGACGGGTATGTTACCGGGG-GTACGACGACTCGAGCCGAGG-ACCTATACCGAA-TCTGAAAAGGAAGAAGAACCTTC-GGTCTCTCTACTTACCTTTGAACGGTCTTCCC-AGCTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGCTCAACGGTACCGGGCCTTGAAGAGGCTAAATACTCGTCG-TAAAAGACGTAATTAGTGAAAGCGGTTTTCCTTCAGGGAC--TTGTGTAA--ACTGAATCCCTTTCACT--GGGGAAATCCTTACCCTC-CGGGTTTCCTTCCGTGTAGAGACCTGGACCATTGT--CGG-TGTGATATTGTCTGTGACAGTCTCGTGCCTAACAC-GTGGGACTGCTACTGCTTTTGCACGACCTTCACCCTAGCC-TACAAAAAACCCGACTTACCGACA----CGAGGTCCGAGAACCTTTTGTTAACAATGCTATAGTAAAGGTTACACGGCGA-AGGAGTCCTC-----------------------------GACAACAAACAACAATCCACACATACCGTCGTA----CCATCTTTGTACCTTGTATTAAAGTATTAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGTCAAA-ACCATAGTTCGACACGAGCTATGAAGAGTTTAGGAAACTCGGCTTGAGCGAGAGGGTTATCAACTTCCTTAAGAAGTTTCCGACCCGGTCGACTGAGTAACAACATGTAGAGAAGAAGTCACGAGACAGGTCGAAACGTCATGTACTACGGATGAAAGACGACCATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTACTCG-------------------ACGAGAC---------------ACTGGAACCCAC--------------------------ACAAC-TTC--------------------------------------------------CGGTTGGCGTCCGTGA-CCT------GTAAC-------TCGAGTCCTTGGTGCACGAGGTCGTGAAGAAACTCTCCGT---CCGGAGACCGACAGCTCACGAACC--GACAGCTGCCTTAC-CGCAGGAGGTCCTGCCTTGGAG--A-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAACATGC-CTC-TAGGTAAAAGAAAG-TAAAGTAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTCACAAGGAAGACCATCTTTTTAAAAGATATTATCATACGTGGGTC-TAGA-TGCACA--AGGTGCTACATAGAAAAGTGAGATAGATCCACGTTGTGGTTCAAGAGGAACCCGTGAAGTT-CATGACGGTGAGGACGTAAACACGTTACCCTTGATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGTCTTAAACT-----AAAGCTTTCCGACTAAAAGGAGAGGGGAAGAGGTGGATAGAGTAAGGCTTTTATACCTGTCGAGTAATGTTCGAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---GAGACACGACGGGT-AACCCCTCGATTTTGACCGACTTACTCTGTCTCC-AGGTCTCTGGTTCCTACGGTGCCGTCGACGCCTTCGAGGTCACCACCACAAAGAGGAGAATTATTG---TAACA----CGGACAAGTATTAACGGAGTAAGTTTAGT--TTATTATATGAG-ACCTGTA-GTAT-CACCCGTTTCGTGAAACCTGTACAGGTCCCGGACTG-TCTC---AGAAGGCT-AGCCAATAGAAGTTCTAGACT-ACT---GAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--ATGCTTTTTATCTACTCAGAGGTC-AAAGTCTAAAACTTGC-AGCGAACTCACCAAAAA--GGGAAAGGAATTTTTCGGACTCCAAGAGA--------GA-GACACTCTATCCGGGC------A-GATAATGTAGGTTTCTCTTTGGTCTTG----GGAACGGTGGTGGTCATTACTCTTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CTAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTTTTCG-AA-AGAGG-CC---TACAAGGGACAAATAGTTAACTTCTACACGATGGTATAAGGGACACGTCACTACCAGGTCAGAACATACATCCACTCTATGTAAGTAGTTTCGATAGGTGTAGTTC-AAGTGAATACCACTATAGTGACAACGAAGCGACCGCATAAGACGCTTCTCGAAAAGACCCC-CTCGACCCCATCCGAAGTCCGTACTCTTTAGACCTCAGGTACAAGAGGCTGCCAGAACAGAACAAAAAGTC--GCACAAGACATACGAAAAGTCG--TACCG--ACTCAGTAACCTTAACTCCGATTG-TCTTG---A--AGGGACTGAGTGTCCTG-AC--GAAAACCTCCTGTCC-GATAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATACAAAAGGAAAGAACCCCCGGTT-AGACGTTAGGA------ACGGGAATGGG----GGTAAAAACCCTGAGGCCCCCAACTCCAATATACCACTCTTGGAAGTGAAGCCCCCCCAGGATGTCCTCGGTATTGACTACTAGGATTCGACGGTAA---AACCCTGGCCTTCGCGGCGTACATCTTTAGTAGTCACACCTGAAACC---CCCACGGCAG-ACACGGTTGTTGTCTCCGGGACAGTCACC--GTTCTCCGGAAAGACAGACACTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-TACAG-------ATACCACAACTTCTGGGTAAGGCACTTGAAGAAGTCGAACAGGCTTCGCAAGACCCTGAGGAGGAC-ATCGGAGTCCAACAGGAGAGCGTCCACGACATGGAGACGGACTTCCACCCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGGGCTTCGTAGAAGTCGTATGGAAGGTCGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTCAGAAATGACGACGACATTCCGAAGTAGACCCGGTCTCGGACCAGGAGATCGACTAGTGACAACCTACGTCGTGAGTTCCCCGCTACCAAAAGTAGTGATC--GTCGTTTTTTGCTCCGAGAGACCTG------GAGCATCCGAAACCGTCGTA--GATGGTCCAACCTGGTTAGGTTCAGACGTCGTCACAGTCCGTCCCCGTAGTCCGGTCC-CAGTCCGTATCCGACGGACGCCCAGAGGACCTACAGTCACTCCCTCG-ACAGGCTCCGGCTCTTGAGAGGTACGTCGCTAAAGTAAGGATTCCA-CATTCTTCCGGAGGAGTTGCCCTCGCTGAACCCCTCTCTTCCCATGCACGATGACT-TCTC-CTGCAGCAACCCGTCCCGAACCAGGTCTTCACTGCGGGACCTCCACTTTTTCA-CGAGGGAGACCGATGACATTTACGAGAGCAGGCTGTCACAAGACACGTTCCGGGCGGCACTCCACCCCCAC--GACCCCTTTCCTCTGACTCCTCTCCT-CTCACGCACCTT------GAATTTTGAGAAGACACTAA-GTCCTCAACCCGAAAGTTCTAGTTAATCAAAGAAATCTAGAAGAAGTTACCACCTGAACCTCCA-ACCGAATCCTCTACGGGCCTCCAGAGAC---CAAC-TCCGTGCGTCGGACCACACCCGATAAA-GTAACGGTAGCGA-GACTACTGCGAAAGAGCACAGGTATTTAGGGTGCC-TCAAATACCAAGGA-AGCC--CAAAAAAGGTAGGTTAGTCCATGAAC-GATCATCGTGGGACCGCTC--CCATTTTTGAGGATCCTGGTCTCCTCCTCCTCCT-TCTCCTCCTC------TTCCTACCACGCTCGACACTCCCCCGTCGGCCTACTTCGAGAAGAAATCCACTTGAAGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGAACACACCCTTCGTCGTAGCAAATATGACTTCCACGACAAGAGGGAACAGTCTCCCGAGCTCGGCCG-AACCGAGAACGACACATACAGATTCGA---GTCCCTTCTGAGACAAAAGGCCAGGTCCGGCACC-GAGAGACAGAGGTTTCAGGCTAGACGGTAC-CT--CTACTGAAGTGACCTCGATACGGTATCCCCACGGTGAAAGAGGAACGGACCGGG---ACTCGCTCCACTGCTACACCGTCGGGGGTCGCTCCGACACGACGGCACCGACCCCGGTATTCTCCACCA----------------CAGTTACGGTGACAGTCGTCTCCACGGAACGTCTATCAATATCGACAACTCAAGCCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGTAGT---GATAC-CCGTGAGCTGCGTAACGACCTCT-CTCATAACGATGACCTCACCGAGCTCATCACCGATAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGTAATTAGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTACAACC----GAGACATATAA-GTACGACGAGC-TCCTGTTTTCTTAACTTTAGGACTCCAGGA-GATATTGACATTCGGTAGACCAGTCGTTACGATAGGAGTGAAGAAACGGTTGATCCGCACCTAAAAGGCGACGGTGGT---TATAGACTGACCCAC-----AAGACGGACAGTTACCTCCGACCCCT-----AGATGCCGAGTGAGGTTTGGGAAGCCCGTACCTGAAGACGGACATTCACGGGGGGTGTCTTGGTCATAGTGACCTACTTCTGTTCAACGATCCTCTCAAGGCTGT---GACCCTAGAACGATTTCGGTACCTACTCTGGCAGTGTATCTACAGGCTTGGTCCTCT-GTCGGCGG---AGAGATCACCCGCGGGCAGTGACCTCTACGAAGCTCCCTACAACCGTCGAGA-A-AGACCTCCTCCGTAGAAGGGTAGAGGAAAAGGTCGTAAGTCGGCGCTCCAGTACCCTGAG-ACGTGGACCCTGCCCT-GAACTCCACT-TACTCC------GGGGAACGACCGCAAGGAGAC-GGTAAACAATA--ATAAAGCCTTCCCTTTGAACTTAAGCATACTTTGTGGAAGTAAGTAAAACATGAGG-------------TAGTGACCGTAC-GCCATCCAAAGTCTCGTTGA------GGAGGCGCCCCACGTCAGTAGTACCCCCCCTTACTACTGTTGGAAGTG--GCACTGTCTACAAGAGTCCTCCAGCTAGTAGACGAGTACCGAGTCTCACCGATGACGGTGGAACGTCTAGAGGT-GGTCCG-CTGACGGGACGGAGTCGGGAC-GGACGGTCCGTATTGGCATCCCGAGGTGTAG-ACGG-TGCTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATCGGGTAAAGAAGCGTCTCAGTAGAGGTGGCTAA--GTGTCGGAAAACT-TTGTTTGACGTTGGAAAACTGTTAGAGAAACTAGAGGAGT---AACAATTATTTTGAGGGGTGACTTTGTGTAAGAGGAGCAAGTATCTTCTAAAAC-TCCGACCACAGTCAACTTCAGGTGGGATAGT-GTAGAGATGTCCTTGTTAACTTCTGGTGCGAGAAACAAAGGACAAGTTACTCGTGTTACTTGAGTAACCTCTAGGGATTTTCCGTCAGATCAAGG-TACCGGAACAACA--TCAGAAGTCGAGCCTGACGGGTATGTTACCGGGG-GTACGACGACTCGAGCCGAGG-ACCTATACCGAA-TCTGAAAAGGAAGAAGAACCTTC-GGTCTCTCTACTTACCCTTGAACGGTCTTCCC-AGCTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGCTCAACGGTACCGGGCCTTGAAGAGGCTAAATACTCGTCG-TAAAAGACGTAATTAGTGAAAGCGGTTTTCCTTCAGGGAC--TTGTGTAA--ACTGAATCCCTTTCACT--GGGGAAATCCTTAACCTC-CGGGTTTCCTTCCGTGTAGAGACCTGGACCATTGT--CGG-TGTGATATTGTCTGTGACAGTCTCGTACCTAACAC-GTGGGACTGCTACTGCTTTTGCACGACTTTCACCCTAGCC-TACAAAAAACCCGACTTGCCGAC----GCGAGGTCCGAGAACCTTTTGTTAACAATGTTATAGTAAAGGTTACACGGCGA-AGGAGTCCT--C---------------------------GACAACGAACAACAACCCACACATACCGTCGTA----CCATCTTTGTACCTTGTATTAAAGTATTAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGTCAAA-ACCATAGTTCGACACGAGCTATGAAGAGTTTAGGAAACTCGGCTTGAGCGAGAGGGTTATCAACTTCCTTAAGAAGTTTCCGACCCGGTCGACTGAGTAACAACATGTAGAGAAGAAGTCACGAGACAGGTCGAAACGTCATGTACTACGGATGAAAGACGACCATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTACTCG-------------------ACGAGAC---------------ACTGGAACCCGC--------------------------ACAAC-TTC--------------------------------------------------CGGTTGGCGTCCGTGA-CCT------GTAAC-------TCGAGCCCTTGGTGTACGAGGTCGTGAAGAAACTCTCCGT---CCGGAGACCGGCAGCTCACGAACC--GACAGCTGCCTTAC-CGCAGGAGGTCCTGGCTTGGAG--A-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAACATGC-CTC-TAGGTAAAAGAAAG-TAAAGTAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTCACAAGGAAGACCATTTTTTTAAAAGATATTATCATACGTGGGTC-TAGA-TGCACA--AGGTGCTACATAGAAAAGTGAGATAGATCCACGTTGTGGTTCAAGAGGAACCCGTGAAGTT-CATGACGGTGAGGACGTAAACACGTTACCCTTGATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGTCTTAAACT-----AAAACTTTCCGACTAAAAGGAGAGGGGAAGAGGTGGATAGAGTAAGGCTTTTATGCCTGTCGAGTAATGTTCGAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---GAGACACGACGGGT-AACCCCTCGATTTTGACCGACTTACTCTGTCTCC-AGGTCTCTGGTTCCTACGGTGCCGTCGACGCCTTCGAGGTCACCACCACAAAGAGGAGAATTATTG---TAACA----CGGACAAGTATTAACGGAGTAAGTTTAGT--TTATTATATGAG-ACCTGTA-GTAT-CACCCGTTTCGTGAAACCTGTACGGGTCCCGGACTG-TCTT---AGAAGACT-AGCCAATAGAAGTTCTAGACT-ACT---GAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--ATGCTTTTTATCTACTCAGAGGTC-AAAGTCTAAAACTTGC-AGCGAACTCACCAAAAA--GGGAAAGGAATTTTTCGGACTCCAAGAGA--------G-AGACACTCTATCCGGGC------A-GATATTGTAGGTTTCTCTTTGGTCTTG----GGAACGGTGGTGGTCATTACTCTTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CTAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTTTTCG-AA-AGAGG-CC---TACAAGGGACAAATAGTTAACTTCTACACGATGGTATAAGGGACACGTCACTACCAGGTCAGAACATACATCCACTCTATGTAAGTAGTTTCGATAGGTGTAGTT-CAAGTGAATACCACTATAGTGACAACGAAGCGACCGCATAAGACGCTTCTCGAAAAGACCCC-CTCGACCCCATCCGAAGTCCGTACTCTTTAGACCTCAGGTACAAGAGGCTGCCAGAACAGAACAAAAAGTC--ACACAAGACATACGAAAAGTCG--TACCG--ACTCAGTAACCTTAACTCCGATTG-TCTTG---A--AGGGACTGAGTGTCCTG-AC--GAAAACCTCCTGTCC-GATAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATACAAAAGGAAAGAACCCCCGGTT-AGACGTTAGGA------ACGGGAATGGG----GGTAGAAACCCTGAGGCCCCCAACTCCAATATACCACTCTTGGAAGTGAAGCCCTCCCAGGATGTCCTCGGTTTTGACTACTAGGGTTCGACGGTAA---AACCCTGGCCTTCGCGGCGTACATCTTTAGTAGTCACACCTGAAACC---CCCACGGCA-GACACGGTTGTTGTCTCCGGGACAGTCACC--GTTCTCCGGAAAGACAGACACTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-GATCGAGTACAGATACCACAACTTCTGGGTAAGACACTTGAAGAAGTCGAACAGACTCCGTAAGACCCTGAGGAGGAC-ATCGGAGTCCAACAGGAGAGCGTTCACGACTTGAAGACGGACCTCCACCCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGGGCTTCGTAGAAGTCGTATGGAAGATTGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTTAGAAATGACGACGACATTCTGAAGTAGACCCGGTCTCGGACCAGGAGTTCGACTAGTGGCGACCTACGTCGTGAGTTCCCTGCTACCAAAAGTAGTGATC--GTCGTTTTTTACTCCGAGAGACCTG------GAGTATCCGAAACCGTCGTA--GATGGTCCAACCTGGTTAGGTCTAGACGTCGTCACAGCCCGTCCCCGTAGTCCGGCCC-CAGTCCGTACCCGACGGAAGCCCAGAGGACCTGCAGTCACTCCCTTG-ACAGGCTCCGACTCTTGAGCGGTACGTCGCTGAAGTAGGGGTCCCA-CATTCTGCCGGAGGCGTTGCCGTCTCTGAACCCCTCTCTTCCCATGCACGACGACTT-CTCC-TGTAGCAACCCGTCCCGAACCAGGTCTTCGCTTCGGGAGCTCCCCTTTTTCAC-GAGGGAGACCGACGACATCTACGAGAGCAGGCTGTCGCAGGACACGTTCCGGGCGGAGCTCCACCCCCAACA--CCCCTTTCCTCAGACCCCTCTACTG-TCGCGCACCTT------GAATTTTAAGAAGACACTAA-GTCCTCAACCCGAAAGTTCCAGCTAGTCAAAGAAATCCAGGAGGAGTTACCACCTGAACCTCCA-GCCGAATCCCCTGCGGGCGTCCAGGGAC---CACC-TCCGGGCGTCGGACCACACCCGATAGA-GCAACGGCAGTGA-GACTACTGCGAAAGAACACAGGTATTTAGGGTGCC-TCAAATACCAAGGA-AGAC--CGAAAAAGGTAGGTTAATCCATAAAC-GATCATCGCGGGACCGCTC--CCATCTTTCAGGATCCTGGCCTTCTCCTCCTCCT-CCTCCTCCTC------TTCCTGCCACGCTCGACACTCCGCCACCGTCCCACTTCGAGAAGAAATCCACTTGAAGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGACCACACCCTTCGTCACAGCAACTATGACTTTCACGACAAGAGGGAACAGTCGCCCAAGCTCGGGTG-GACCGAGAACGACACCTCTAGGTTCGAG---TCCCTCCTGAGACAAAAGGCCAGGTATGGTACTGG-GAGACAAAGATTTCAGGCTAGACGGTTC-CT--CTAGTGAAGTGACCTCGATACGGTGTCTCCACGGTGAAAGAGGAACGGACCCGG---GCTCGCACCACTGCTGCACCGTCGGGGGTCCCTCCGGCACGACGGCACCGACCCCGGTATCCTCCACCA----------------CAGCTACGGTGACAGTCGCCTCCACGGAACGTCTATCAATATCGACAACTCGAGCCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGCAGT---GACAC-CCGTGAGCTTCGCAACGACCTCT-CCCATAACGATGACCTCACCGAGCTCATCACCGATAGGCTTACCGGAACAAACAAACCTCTGTCTACTTGGTGCCCT--AGGACGGCGTAATTAGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTGCAACC----GAGACATATAAG-TACGACGAGCTCC-TGTTTTCTTAACTTTAGGACTCCAGGAG-ATATTGACATTCGGTAGACCAGCCGTTACGATAGGAGCGAGGAAACGGTCGACCCACACCTAAAGGGTGACGGCGGCGGCTACAGACGGACCCAC-----GAGACGGACAGTTACCTCCGGCCCCT-----AGACGCCGAGTGAGGTTTGGGAAGCCCGTACCTGAAGACGGACATCCACGGGGGGGGTCTCGGTCACAGTGACCTGCTTCTTTCTAACGACCCCCTCAAAGCTGT---GACCCTAGAACGATTGCGGTACCTACTCTGACAGTGTATCTACAGGCTCGGTCCTCT-GTCGGCGG---AGAGCTCCCCCGCGGGCAGTGACCTCTACGAAGCTCCTTATAACCGCCGAGA-A-AGGCCTCCGCCGTAGAAGGGCAGAGGAAAGGGTCGTAAGTCGGCGCTCCAGTACCCTGAG-ACGTGAACCCTACCCT-GAACTCCGCT-TACTCC------AGGGAACAACCGTAAGGAGAC-GGTAAACAATA--ATAAAGCCTTCCCTTTGAACTTAAGCATACTTTGCGGAAGTAAGTAAAACATGAGG-------------TAGTGACCGTAC-GCCATCCAAAGACTTGTCGA------GGAGGCGCCCCACGTCAGCAGTACCCCACCTTACTACTGCTGGAAGTG--GCACTGCCTCCAGGAGTCCTCCAGCTAGTAGACGAGCACCGAGTCTCACCGGTGGCGGTGGAACGTCTAGAGGT-GGTCCG-CCGACGGGACGGAGTCGGGAC-GGACGGTCCGCATCGGCATCCCGAGGTGCAG-CCGG-TGCTACGGGAGCATCAACTGCACCTTCGGTTCGACGGGATTGGGTAAAGAAGCGTCTCAGTAGAGGTGACTAA--GTGTTGGAAAACT-TTGTTTGACGTTGGAAAACTGCTAGAGAAACTATAGGAGT---AACAATTACTTTGAGAGGTGACTTTATGTAAGAGGAGCAAGCATCTTCTAAAAC-TCCGACCACAGTCAACTTCAGGTGGGATAGT-GTAGAGACGTCCTTGTTAACTTTTGGTGCGAGAAACAAAGGACAAGCTACTCCTGTTACTTGAGTAACCTCTAGGGGTTTTCTCTCAGATAAAGG-TACCGGAACAACA--TCAGAAACCGGGCCTGCCGGGTGTGCTACCGAGG-GTACGATGACTCGAGTCGAGGA-CCCATACCGAA-TCTGAAAAGGAAGAAGAACCTCC-GGTCTCTCTACTTACCTTTGAACGGTCTTCCC-AGTTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGCTCAGCGGTACCGGGCCTTGAAAAGACTGAATACTCGTCG-TAAAAGGCGTAATTAGTGAAACCGGTTTTCCTTCAGGGAC--TTGTGCAA--ACTGAATCCCTTTCACT--GGGGAAATCCCTACCCTC-CGGGTTTTCTTCCGTGTAGAGACCTGGATCACTGC--CGG-TGTGACATTGTCTGTGGAAGACTCGTCCCCAATAC-GCGGGACTGTTACTGCTTTTGCACGACTTTCACGCTGGCC-TACAAAAAACCCGACTTGCCGACA----CGAGGTCCGAGAACCTTTTGCTAACAATGTTATAGTAAAGGCTACACGGCGA-AGGCATCCTC---------------------------GACAACGAACAACAATCCACACATACCGTCGTACC----ATCTTTGTACCTTGTATTAAAGTATCAGAAACATACAACACATGTCTCTTACCCAAGGTCACCTCAGCCAAAACCAC---AG---GACACGTGCTATAAAGAGTTAGGGATCCTGAGC---------AGGGTCATCGACTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNG---------------------------------------------------------------------------------------------------------------------GGAACGAGTAGGTCAAAGTACTCG-------------------ACGAGAC---------------ACTGGAACCCTC--------------------------ACAAC-TTC--------------------------------------------------CGGTTCGCGTCCGTGA-CGT------GCAAC-------TCGAGCCCTTAGTGCACGAGCTCGTGAAGAAAGTCTCCGT---CGGGAGACCGACAACTCACGAACC--GGCAACTCCCTTAC-CGTAGGAGGTTCTGTCTTGGAG--A-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAATATAC-CTC-TAGGTAAAAGAAAG-TAAAGAAGACCAAAAAGGGAAGCTTACGAGTGCCGGTGGCCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTCACAAGGAAGACCATCTTTTTAAAGGATATCATCATCCGCGGTTC-TAG-GTGCACA--AGGTGTTACATAGAAAAGTGAGATAGATCTACGTTGTGGTTTAAAAGGAACCCGTGAAGTT-CATGACGGTGAGGACGTAAACACGTTACCCTCGATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGTCTTAACCT-----AAACCTTTCCGACTAAAAGGAAAGCGGAAGAGGTGGATAGAGTAAGGCTTTTATCCCTGTCGAGTAATGTTCGAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---GAGCCACGACGGGT-AGCCCCTCGATTTCGACCGCCTTACTCTGTTTCC-GGGCCTCTGGTCCCTACGGTGCCGCCGGCGGCTGCGGGGCCACCACCACAAGGAGGAGAACTAGTG---TAACA----CGGACAAGTATTAGCGGAGCAAGTCCAGC--TTATTATACGAG-ACCTGTA-GTAT-CACCCGTTTCGTGAACCCCGTGCAGGTCCCGGACTG-CCTCCTCAGAAGGCT-GGCCAACAGAAGGTCCAGGCT-GCT---GAAGGGCT---CCAG----GAGCGGAAA-CCGTAA--CGCGTTGA--ACGCCTTCTACCTCCTCAGCGGTC-AAAGGCTGAACCTCGC-CGCGAACTCACCAAACA--GGGAGAGGAATTTCTCGGACTCCAAGAGA--------GA-GACACTCTCTCCGGGC------A-GTTAGTGCAGGTTCCTCTTTGGCCTTG----GGAACGGCGGCGGCCGCTACTCCTTCAGCG----ACTCGAACCAACCGCAACCCGGGG---CTAGCCCAAACAGCAGAAACATCGGGAAATTGTTCAGCTTTTTCG-AG-AGGGG-CC---TGCAGGGGACGAACAGCTACCTGCTGCACGACGGCATGAGGGACACGTCCCTGCCCGGGCGGAACATGCACCCACTCTATGTAAGTAGCTTCGACAGGCGCAGTTC-GAGTGAATACCACTACAGCGACAACGACGCGACGGCATAGGACGCGTCGCGGAAGGCCCCG-CTCGACCCCATCCGGAGCCCGTACTCCTTGGCCCTCAGGTACGAGAGGCGGCCGGCGCGGAACGAGAACTC--GCACAAGACCTACGAGAAGTCG--TACCG--GCTCAGCGAGCTCAACTCCGATTG-TCTCG---A--CGGGACTGAGTGTCCCG-AC--GAAAACCTCCTGTCC-GAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGACCTTGCCGA---TCTAAGGAAATACAAAAGGAACGCACCCCCGGTC-AGACGCTAGGA------ACGGGAGTGCG----GGTAAAACCCCTGAGGCCCCCAGCTCCACTACACCACGCTCGGAAGTGAAGCACCCCCAGGGTGGCCTCGGTACTGCCTCCTAGGATTCGACGGCAAAAC---CCTGGCCTTCGCAGCGTACATCTTCAGCAGTCACACGTGAAACC---CCCGCGGTAG-ACACGGCTGTTGGCTCCGGGACAGCCACC--GTTGTCCGGAAAGGCAGACGCTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-GATCGAGTACAGATACCACAACTTCTGGGTGAGACACTTGAAGAAGTCGAACAGGCTCCGTAAGACCCTAAGGAGGAC-GTCGGAGTCCAACAGGAGAGCGTCCACGACTTGAAGACGGACTTCTACCCGAAACAGAAGAAAAAG-------------GAA---GAAGTTTAGAAGTGCTTCGTAGAAGTCGTATGGAAGGTCGACCAGGTGAAATCTTCGGTCTCATCCTCTCCTTAGAAATGACGACGACATTCTGAAGTAGACCCGGTCTCGGACCCGGAGCTCGACTAGCGACAACCTACGTCGTGAGTTCCCTGCTACCAAAAGTAGTGATC--GTCGTTTTTTACTCCGAGAGACCTG------GAGTATCCGAAACCGTCGCA--GATGGTCCAACCTGGTTAGATCTAGGCGTCGTCACAGTCCGTCCCCGTAGTCCGGCCC-CAGGCCGTACCCGACGGAAGCCCAGAGGACCTGTAGTCACTCTCTCG-ACAGGCTCCGACTCTTAAGAGGTACGTCGCTGAAGTAGGGCTCGCA-CATCCTTCCGGACGCGTCGCCGTCTCTGAACCCCTCTCTGCCCATNNNNNNNNNNNN-NNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNN------GAATTTTAAGAAGACACTAA-GTCCTCAACCCGAAAGTTTCAGCTAGTCAAAGAAATCTAGAAGGAGTTACCACCTGAACCTCCA-ACCGAATCCGCTCCGGGCCTCCAGGGAC---CACC-TCCGTGCGTCGGAACACACCCGATAAA-GTAACGGTAGTGA-GACTACTGCGAAAGAACACAGGTACTTAGGGTGCC-TCAAATAACAAGGA-AGAC--CAAAAAAGGTAGGTTAATCCATAAAC-GATCATCGCGGGACCGCTC--CCATCTTTCAGGACCCGGGTCTCCTTCTCCTCCT-CCTCCTCCTC------TTCCTGCCACGCTCGACACTCCTTCGCCGTCCCACTTCAAGAAGAAACCCACTTGAAGCCCGTGACCACTTTTGGTGACCCGACCG-TC-TCACAGCCCTCGAGTACCATGACCACACTCTTCGTCGTAGCAAATATGACTTGCACGACAAGAGGGAGCAGTCTCCCAAGCTCGGTCG-GCCCGAGAACGACACCTATAGGTTCGAG---TCCCTCCTGAGACGAAAGGCCAGGTACGGCACCGA-GAGACGCAGGTCTCAGGCTAGACGGTAC-CT--CTAGTGGAGTGACCTCGACACGGTGTCCCCCCGGTGGAAGAGGAACGGACCCGG---ACTCGCCCCACTGCTGCACCGCCGGGGGTCCCTCCGGCACGACGGCACCGACCCCGGCATCCTCCACCA----------------CAGCTACGGTGACAGTCGCCTCCACGGGACGTCTATCAACATCGACAACTCGAGCCTAAGGCACGCCAGGGGT---AGCATGGTGAAG-AGCAGT---GACAC-CCGCGAGCTCCGCAACGACCTCT-CCCATAACGACGACCTCACCGAGCTCATCACCGACAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGTAATTGGA--GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTAT-GACGACCCTCT-ACCGTACAGCC----GANNNNNNNNNN-NNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNGGTAGACCAGCCGTTACGACAGGAGTGAAGAAACGGTCGAACCACACCTGAAAGGTGACGGTGAC---TACAGGCGGACCCAC-----GAGACGGACAGCTCCCTCCGACCCCT-----GGACGCCGAGTGAGGTTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNGCGG---AGAGCTCACCCGCGGGCAGTGACCTCTACGAGGCTCCTTACAACCGCCGGGA-A-AGACCTCCTCCCTAAAAGGGCAGAGGAAAAGGTCGTCAGTAGGAGCTCCAGTACCCTGAG-ACGTCGACCCTGCCCT-GAACTCCACT-TACCCT------AGGGAACAACCGTAAGGAGAC-GGTGAACAATA--ATAAGGCCTTCCCTTTGAACTTAAGTATCCTCTACGGAAGCAAGTAAAACATGAGG-------------TAGTGACCGTAC-GCCATCCAAAGGCTCGTTGA------GGAGGCACCCCACGTCAGTAGTACTCCGCCTTACTACTGGTGGAAGTG--ACACTGCCTCCAAGAGTCCTCTAGCTAGTAGACGAGCACCGAGTCCCACCGGTGGCGGTGGAACGTCTAAAGGT-GGTCCG-CCGACGGGACGGAGTCGGGAC-GGACGGTCCGTATCGGCATCCCGAGGTGTAG-ACGG-TGCTACGGGAGCATCAACTGCACCTTCGGTTCGACAGGATTGGGTAAAGAAGCGTCTCAGTAGAGGTGTCTAA--GTGTCGGAAAGCT-TTGTTTGACGTTGGAGAACTGTTAGAGGAACTACAGGAGT---AACGATTACTTTGAGAGGTGGCTTTATGTAAGAGGAGCAAGCATCTTCTAAAAC-TCCGGCAACAGCCGACTTCAGGTCGGATAGT-GTAGAGACGTCCTTGTCAATTTTTGGTGCGAGAAACAGAGGTCAAGCTACTCGTGCTACTTGAGCAACCTCTAGGGATTCTCCGTCAGGTAGAGG-TACCGCCACAACA--TCAGAAAGCGGGCCTGACGGGTGTGCTACCGAGG-GTACGACGACTCGAGTCGAGGG-CCTATGCCGAA-TCTGAAAAGAAAGAAGAACCTCC-GGTCTCTCTACTTACCTTTGAACGGTCTCCCC-AGCTATCCCCACAACCACTGTAAGAAAAGACGGTCTAGAAAATCCATAAGGACCGCCCAGGCCCAACGGTACCGGGCCTTGAAAAGACTAAATACTCGTCG-TAAAAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGAC--TTGTGTAA--GCTGAATCCCTTCCACT--GGGGAAATCCCTACCCTC-CGGGTTTTCTTCCGTGTAGAGACCTGGACCACTGC--CGG-TGTGACATTGTCTGTGGCAGTCTCGTCCCTAATAC-GCGGGACTGTTACTGCTTTTGCACGACTTTCACTCTAGCC-TACAAAAAACCCGACTTGCCGACA----CGAGGTCCGAGAACCTTTTGTTAACAATGCTACAGTAAAGGTTATACGGCGA-AGGCATCCTC---------------------------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN----NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAAACCAC---AGTTCGACCCGTGCTATGAAGAGTTTAGGAAACTCGGCTTGAGCAAGAGGGTCCTCGACCTCCTTAAGGAGCTTCCGCCCCGGCCGGCGGAGCAACAACATGTAGAGGAGGAGCCACGAGACAGGCCGAGACGTCATGTACTACGGGTGAAAGACGACCATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGGACAAGTAGGTCAAAGTACTCA-------------------ACGAGAC---------------ACTGGAACCCCC--------------------------ACAAC-TTC--------------------------------------------------TGGTTGGCGTCCGTGA-CCT------GCAAC-------TCGAGCCCTTAGTGCACGAGGTCGTGAAGGAAGTCCCCGT---CAGGGGACCGACAGCTCACGAACC--GACATCTCCCTTAC-CGCAGGAGGTTCTGACTCGGAG--A-CTCTCAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAACATAC-CTC-TAGGTAAAAGAAAG-TAAAGAAGACCAAAGAGAGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAG-----AGTAGTCGGGGTTTTATCAAGGTCACAAGGAAGACCATCTTTTTGAAGGACATTATCATCCGCGGCTC-CAG-GTGCACA--AGGTGCTACATAGAGAAGTGAGACAGGTCTACGTTGTGGTTCAAGAGGAACCCGTGAAGTT-CATGACGGTGAGGACGCAAACACGTTACCCTCGATTCCCAGCTGAGTTTCCTCCG--TGA--AAGTGTAGGTG--ATTTCCTCGGACTTAACCT-----AAAGCTTTCCGACTAAAAGGAAAGGGGAAGAGGGGGCTAGAGTAAGGCTTTTATGCCTGTCGAGTAATGGTCAAGTAACGAAAACAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAACTCCAG---GAGGCACGACGGGT-AACCCCTCGATTTTGACCGACTTACTCTGTCTCC-AGGTCTCTGGTTCCTACGGTGCCGTCGACGGCTCCGGGGTCACCACCACAAAGAGGAGAATTATTG---CAACA----CGGACAAGTCTTAACGGAGCAAGTTTAGT--TTGTTATATGAG-ACCTGTA-GCAT-CACCCGTTTCGTGAAACCTGTGCAGGTCCCGGACTG-TCTC---AGAAGGCT-AGCCAATAGAAGTTCTAGACT-ACT---GAACGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--AGGCTTTTTACCTCCTCAGCGGTC-AAAGTCTGAACCTTGC-GGCGAACTCACCAAAAA--GGGAAAGGAATTTCTCGGACTCCAAGAGA--------GA-GACACTCTATCCTGGT------A-GATATTGTAGGTCTCTCTTTGGACTTG----GGAACGGGGGTGGCCGCTACTCCTTCAGTG----ACTCGAACCACCCCCAGTCTGGAG---CGAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTCTTCG-AA-AGAGG-CC---TACAGGGGACAAACAGTTAGCTCCTACACGATGGCATGAGGGACACGTCCCTGCCGGGTCAGAACATACACCCCCTCTATGTAAGTAGTTTCGATAGGTGCAGTTC-AAGTGAATACCACTATAGCGACAACGAGGCAACGGCATAAGACGCTTCTCGAAAGGACCCC-CTCGACCCCATCCGAAGTCCGTACTCCTTAGAGCTCAGGTACAAGAGGCTGCCAGAGCAGAACGAAAAGTC--GCACAAGACATACGAGAAGTCG--TACCG--ACTCAGTAACCTTAACTCCGATTG-TCTTG---A--AGGGACTGAGTGTCCTG-AT--GAAAACCTCCTGTCC-GAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATACAAAAGGAAAGCACCCCCGGTT-AGACGTTAGGA------ACGGGAGTGGG----GGTAAAACCCCTGTGGTCCCCAACTCCACTACACCACTCTTGGAAGCGAAGCGCCCCCAGGATGGCCTCGGTACTGACTCCTGGGATTCGACGGCAAGAC---CCTGGCCTTCGCAGCGTACATCTTCAGTAGCCACACCTGAAACC---CCCACGGCAG-ACACGGTTGTTGGCTCCGGGACAGTCACC--GTTGTCCGGAAAGACAAACACTGGCAAAGTTCGACCAGTAAA");
    // querys.push_back("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAA---GAAATCCAGAAGTGCGTCGTA-GAAGTC-GTACGGAAGGTTG-ACCAGGTGGAAGCTCCGGTTCCACCCTCTCCTCAGGAACGACNNNNNNN-NNNNNN-NNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNG------GAG-CATCCGAAACCGTCGTAGCTGGTTCAACCTGGTCA-GGTCTAG-CCGTCGTCACAGCCCGTCCCCGTAGTCCGGTCCCAGGCCGTATCCGG--CG---ACTCAGAGGACCTACAGCCACTCCCTC--GACAGGCTCCGTCTCTTGAGAGNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------GAATTTTAAGAAGACTCTCAGTCGGCAACC-CGAGAGGTCCAGCTAGT-CAAAGAAGTCCAGGAGGAGTTACCACCTGAACCT-CCAACCAAACCCCCTGCGCGCCTCCAGGGCCCACCTGCGGGCGTCGGACCGCCCTCGGTACAGCGACGGCAGTGAGACTACTGCGAAAGAACACA-GGTATTTAGGGTGCCTCA--ACTACCAGGGAAGTCCAAAAAAG--GTAGGTTAGTCT-ATGAACGACCACCGGGGGACCGCCCCCATCTTCCAGGAT-CCTGGTCTTCACCACCTCCTGCTACTACTC------TTCC-TACCTCGTTCGACTCCCCACCGCCGCCCCACGTCGAGA--AGGAACCCCCTCGAAGACCGTGACCTC--TTCTGCTGCCCCGACCGACTCACAGCGCTCGAGTACCATGACCACACCCTCCGCCGCAGTAAATACGA-CTTTCACGACAAGAGGGACCAGTCTCCCGACCTCGGTTGGACCGAGAACGACACCTACAGGTCCGAGTCCCTCCTAAGAC-GAAAGGCCAGGTACGACCCCGAGAGCCAG-AGCTCGCAGGC-TAGTCGGTACCT-CTAGTGGAGTGACCTTGACACGGTGTCCCCACACTGAAAGAGAAACGGGCCGGG---ACTCGCTCCTCTACT---CCGTCGAGGGTCCCTCTGCGACGAGGGCACCGATCCGGGTATCCTCCACCA-CAGCTACGGCGACAGTCGTCTGCACGGAACGTGTATCAACATCGACGACTCGAGCCTAAGACAAGCCAGGGGC---AGTATGGTGAAG-AGTAG---CGACAC-CCGCGAGCTACGCAACGACCTCTCC-CACAACGACGACCTCACCGAACTCATCACCGACAGTCTTACGGGAACAAACAAACCTCTTTCTACTTGGTGAC--CTAGGACGGCG-TAATTAGA-GGGTGACCCCGAGAGAGCCTTTCGACCGTC-CTGGAACCTTTG-TCGACCCTCC-ACCGTGCACCCAAGACATATAAGTACGACGA--GCT---GTTTTCTTGACTTTAGGACTCCAGGAGATATTGACATTCGTTAGACCAGCTGTTACGGCAGG-AGTGAAGAAACGGTCGA----CCCACACCTAAAGGGT-GACGGTGA--G---TAC---AGCCTGACCCACAAGACGGCCAGTT---ACCTTCGACCCCTA-GACGCCGACTGTGGTTTGGGAAGTCCGTACCTAAAGACAGACATTCACGGGGGCGGTCTCGGTCAGAGTGACCTACTTCTTTCCAACGACCCCCTCAATGCTGT---GACCCTAGAACGGTTACGGTACCTACTCTGCCAGTGTATCTACAGGCTTGGTCCTCTG-TCGGCGG---AGAGGTCACCCGCGGGCAGTGACCTCTACGAGGCTC--CTTACAACCGTCGAGAAAGACCTCCTCCGTAGAAAGGCAGAGGAAAGGGCCGTAAGTCGGCACTACAGTACCCCGA-ACCGTGGAC-GCGGCCCT-GGCCGCCCCGTACTCC------GGGGAACGACCGTGAGGCGACGG-TAAACAACAACAACGCCTTTCCTTTGAATTT--AAGCATCCTTTGAGGTAGTAAGTGAAACATGAGG-------------TAGTGTCCGTACGC-CTTCCAGAGCCTTATCGA------GGACGCACCCCACGTCAGTAGTACACCCCCCTACTACTGTTGGAAGTGGCAC--TGTCTCCAGGAGTCTTCTAGCTAGTAAACGAGTACTGAATCCCACCGGTGTCGCTGGAACGTCTAGAGGTGG-TCCGCCG-ACGGGACGGAGTCGGGACGG-ACGGTCCGTATCGGTATCCC-GAGA-TGTAGCCGGTGTTACGGGAGCATCAACTGTACCTTCGGCTCGACAGGTTCAGGTAAAGAAGCGTCTCAGTAGAGGTGGCTAAGTGTTGGAAAACTTTGTTTGACGTTGGCA---AACTGCTAGAGAAACTACAGGAGT---GGCAAGTACTTCGAGAGGTGGCTCTGTGTAAGAGGAGCGAGCATCTT-CTAAAACTCCGACCACAGTCAACTGCAGGTGGGATAGTGT-AGAGACGTCCTTGTCAACTTTTGGTGTGAGAAACAAAGAACAAGTTACTCCTGCTACTTGAGCAACCTGTGAGGGTCTTCCGTTAGATAGAGGTA-CCGACACAACATCAGAAAC--CGAGCGTGCTGGGTGTGTTACCGGGG-GTACGACGACTCGAGTCGAGGACCCATACCGAATCTAAAGA--GAGAAAAGAATCTCCGGTC-TCTCTACTTACCCTTGAACGGGCTCCCCAGCTATCCCCACAACCACTG-TAAGAAGAGACGGTCTAGAAAA-TTCATAAGGACCGC-TCAGGCTCACCGGTACCGGGCCTTGAAGAGGCTGAACACTCGTCGTAAAAGACGTAATTAGTGGAACCGGTTTTCCTTCAGGGAC-TTGTGTAA--ACTGAATCCCTTTCACT--GAGG---TCCCTACCCTCCAG-GTTTCCTTCCCTGTAGAGACCTAGACCAC--TG-CCGGTGCGACATTGTCTGCGACAGTCACGGTCCCAA-CACTCGGGACTGCTACTGCTTCTGCACGACTTTCACCCTAGCC-TACAAAAAGCCCGACTTACCAACACGAG---GTCCGAGAACCTTTT-GATAACAATGTTATAGTAAAGGTTACACGGCGAAGGCATC-CTC---------------------------GACAACGAACAACAACCCACACATGCCGTCGTATCATCTCTGTACCT-TGTACTAGAGTA--TCAGGAACATACAACA-CATGTCTCTTAC-CCAAGGTCACCTCAG--TCAAAACCACAGATCAACACGTGCCATAAAGAGCTTGGGAAACTCGGCTTGAGCAAGAGGGTCGTCGACTTCCTTAAGGAGCTTCCGGCCCGGCCGCCGAAGTAACAATATGTAGAGAAGAAGTCACGAGACGGGTCGAAACGTCATGTACTACGGGTGAAAGACGACTATCAAGTTG---------------------------------------------------------------------------------------------------------------------GGAACGAGTAGGTCAAAGTACTCG-------------------ACAAGAC---------------ACTGGAACCCCC--------------------------ACAACT-TC--------------------------------------------------CG---GTTC---GCGTCCGTG-ACCTGCGAC-------TCGAGTCCTTAGTGCACGAGGTCGTGAAGGAAGTCTCCGTCGGG---AGCCCGACAACTCACGAACCGACAACTTC-CTTACCGTAG--GAGGTCCTGGCTTGGGGA-ATCTCA-GTACTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAC----CTGGTTCTCAGAGT--AGTCG---CCTTAT---CAAGGTCACAAGGAAGACCATCTTTTTGA--AGGACATTATCATCCGCGGCTCTAGGTGCACGAGGTGCTAC-ATGGAAAAGTGAGCTA-GATTTACGTTGTGGTTCAAGAGGAACCCGTGAAGTTC-CTGACGGTGAGGTCGTAACTACGTCACC---GACTCCCATCTGAGTTTCCTCCG-TGAA---AGTGTAGGTG--ATTTCCTCGGTCTTAACC--TAA---AGCTGTCCGACTAGAAGGAAAGTGGGAGGGGTGGATAGAGTAAGGCTTTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTACGACTCC------GAAGGAAAC---AAGTCCAG---GAGGCACGACGGGTA-ACC---CGACTTCGACCGTCTTACCCTGTTTCC-CGGTCTCTGGT-TTCTGCG---TC-GCCGACGGCTCCGAGGACACCACCACAA---------GAATTATTGTAACA---CGGACAAATGTTAGCGGAGCAAGTCTAGTTTGTCATACGAG--ACCTGTA-GCATCA-CCCGTTTCGTAAACCCCGTGCAGGTCCCGGACTGC-GTT---AGAAGACTAGCCAATAGGAG-TTCTAGACTA-CT---AAAAGGTT---CAAG---T-AGTGGAAA-ACGTAAC--GCATTAA--ATGCTTTTTACCTACTCAGAGGTCA-GAGTCTGAACCTTGCG-GCGAACTCACCAAAGAGGG--AAAGGAATTTGTCGGACTCCACCAGAGA---------GACACTTTACCCCGGT---AGATAT----TGTAGGTTTCTCTTCGGTCTTGGG---AACGGGGGAGGCCA-CTACTCTTTTAGTGA----CTCGAACCAACCCCAACCTGGAG---CCAGTCTAAACAGTAGAAACATCGGGAACTTGTTCAGTTTCTTTGACAGAGGTCTACAGGGGACAAACAGTTAACTTCTACACGAGGGTATAAGGGACACGTCGCTACCGGGTGAGAACATACACCCACTCTATGTAAGCAGTTTCGATAGGTGTAGTTCGAGTGAATACCATTACAGTGATAACGAGGCGACGGCCTAAGACGCTTCGCGAAAGGATCCCCTCGACCCTATCCGAAGCCCGTAGTCTTTCGCCCTCAGGTACAAGAGGCGGCCGGAGCAGAACAAGAACTCACACAAGACATACGAGAAGTCGTACCGACTCAGTGACCTTAACTCCGATTGTCTTGAAGGGACTGAGTGTCCTGACGAAAACCTCCTG-------TTC-GAAAGTTCCGTTGAA--CTCCAAA-GAT--AAAGGACCGAACT-----TGCCG-A---A-CTAAGAAAATACAAAAGAAAGGCACC---CCCGGTT--AGGCGTTAAGA------ACGGGCTTGGGGGTAAAAACCCTGAGGCCCCCAGCTCCAATACACTACGCTCGGAAGC-AACGCCCCACCGGG--GTGTCCTCGGTACTGACTACTAGGATTCGATG-GTAAAACCC---TGGCCTTCGCGGCGTACATCTTCAGTAGCCACACCTGAAACC---CCCGCGGTAGACGCGGCTGGTACCTCCGGGACAGCCACCGGTGCCCAGAGAGCCAGACACTGGCAA-AGTCCGACCAG-TA---");
    // querys.push_back("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAA---GAAGTCCAGGAGCGCGTCGTA-GAAGTC-GTACGGAAGGTCG-ACCAGGTGGAACCTCCGGTCCCATCCGCTCCTCAGAAACGACGACGGCA-TGCTGA-AGTAGACCCGGTCTCGGACCAGT-AGCTCGACTAGCGGCAACCTCCGTCGTGAG-TCCCCTGCTACTAAAAGTAGTGATCGTCGTTTTTTGCT-CCGAGAGACCTG------GAG-TATCCGAAACCGCCGTAGATGGTTCAATCTGGTTA-GGTTTAG-ACGTCGTCACAGTCCGTCCCCGTAGTTCGGTCCCAGGCCGTACCCGA--CGGAAGCTCAGAGAACCTACAGCCACTCCCTC--GACAGACTCCGGCTCTTGAGAGGCACGTCCCTGAAGTAAGGGTT-CCACATCCTGCC-GGAGGCGTCGCCCTCGCTGAACCCC-TCCCTCCCCATGCACGA-CGACTTCTCCTGTAGCAACCCGTCCCGGACCAGGTCGTCGCTCCGGGAGCTGCACTTCTTCACGA-GGGCGACCGACGACATTTACGAGAGTAGGCTGTCGCAGGACACGTCGCGGGCGGACCTCCACCCCCACGACCCCTTCCCCCTGACGCCCCTCCTGTCGCGCACCTT------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAAGAACACA-GGTATTTAGGGTGCCTCA--AATACCAAGGAAGCCCGAAAAAG--GTAGGTTAGTTC-ATGAACGAGCACCGCGGGACCGCCCCCATCTTTCAGGAC-CCTGGCCTGCTCCTCCTCCTTCTACTCCTC------TTCC-TCCCGCGCTCGACCCCCCGCCACCGTCCTACGTCGAGA--AGGAACCCCCTTGAAGACCGTGACGAC--TTCTGCTGTCCCGACCGTCTCACGGCCCTCGAGTACCATGACCACACCCTCCGCCGTAGCAACTACGA-CTTTCACGACAAGAGGGACCAGTCTCCCGACCTCGGCCGGACCGAGAACGACACCTACAGGTCCGAGTCGCTCCTGAGGC-AAAAGGCCAGGTACGGTCCCGAGAGACAG-AGGTCGCAGGC-TAGACGGTACCT-TTAATGGAGTGACCTCGACACGGTGTCGCCACGGTGAAAGAGGAACGGACCGGG---ACGCGCTCCGCTGCTGCACCGTCGGGGGTCCCTCTGCCACGACGGTACCGACCCCGGCATCCTTCACCA-CAGCTACGGTGACAGCCGCCTCCACGGAACGTCCATCAACATCGACAACTCAAGTCTGAGGCACGCCAGGGGC---AGCATGGTGAAG-AGTAG---CGACAC-CCGTGAACTGCGCAACGACCTCTCC-CACAACGACGACCTCACCGAGCTCATCACCGACAGGCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGAC--CTAGGACGGCG-TAATTAGA-GGGTGAGCCCGAGAGAGCCTTTCGACCGTC-CTGGAACCTGTG-ACGGCCCTCT-ACCGTGCACCCGAGACACATAAGTACGACGA--GCTCCTGTTTTCTTAACTTTAGGACTCCAGGAGATATTGACATTCGGTAGACCAGCTGTTACGACAGG-AGTGAAGAAACGGTCGA----TCCACACCTAAAGGGC-GATGGTGA--C---TAC---AGGCTGACCCACAAGACGGCCAGTT---ACCTCCGACCCCTA-GACGCCGACTGCGGTTTGGGAAGTCCGTAGCTAAAGACAGACATTCACGGGGGTGGTCTCGGTCAGAGTGACCTACTTCTGTCCAACGACCCGCTCAATGCCGT---GACCCTAGAACGATTACGGTACCTACTCTGGCAGTGTATCTACAGGCTTGGTCCTCTG-TCGNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNN-NNNNNNNN-NNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------------NNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNN------NGAGGCGCCCCACGTCAGTAGCACCCCTCCCTACTACTGTTGGAAGTGGCAC--TGTCTCCAAGAGTCCTCTAGCTAGTAAACGAGTACCGAGTCTCACCGGTGCCGATGGAACGTCTAGAGATGG-TCCGCCG-ACGGGACGGAGTCGGGACGG-ACGGTCCGTATCGGCATCCC-GAGG-TGCAGACGGTGCTACGGGAGCATCAACTGCACCTTCGGTTCGACAGGATCGGGTAAAGAAGCGTCTCAGTAGAGGTGCCTAAGCGTTGGAAAACTTTGTTTGACGTTGGAG---AACTGTTAGAGAAACTATAGGAGC---GACAAGTACTTCGAGAGATGCCTCTACGTAAGAGGGGCGAGTATCTT-CTAAGACTCCGACCACAGGCGACTGCAGGTGGGATAGTGT-AGGGACGTCCTTGTCAACTTTTGGTGCGAGAACCAGAGGACGAGTTACTCGTGCTACTTGAGCAACCTCTAAGGATCGTCCGTCAGATAAAGGTA-CCGGCACAACATCAGAAAC--CGGGCCTGCCGGGTGTGCTACCGAGG-GTATGATGACTCGAGTCGAGGACCCATACCGAATCTGAAAA--AGAAAAAAAATCTCCGGTC-TCTCTACTTACCCTTGAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNN--NNNNNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNN--NN-NNNNNNNGACATTGTCTGCGACAGTCACGTCCCCAA-TACTCGGGACTGCTACTGCTTCTGCACGACTTTCACCCTAGCC-TACAAGAAGCCCGACTTACCGACACGAG---GGCCGAGAACCTTTT-GATAACAATGCTATAGTAAAGGTTACACGGCGAAGGAATT-CTC---------------------------GACAACGAACAACAACCCACACATGCCGTCGTACCATCTCTGTACCT-TGTATTAGAGTA--TCAGAAACATACAACA-CATGTCTCTTAC-CCAAGGTCACCTCAG--TCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNG---------------------------------------------------------------------------------------------------------------------GGAACGAGCAGGTCGAAGTACTCG-------------------ACGAGGC---------------ACTGAAACCCCC--------------------------ACAACT-TC--------------------------------------------------CG---GTCC---GCGTCCGTG-ACGTGCAAC-------TCGAGTCCTTAGTGCACGAGGTCGTGAAGGAAGTCCCCGTCGGG---GGACCGACAGCTCACGAAACGACAACTCC-CCTGCCGTAG--GAGGTCTTGTCTTGGAGA-GTCCCA-GTACTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN----NNNNNNNNNNNNNN--NNNNNNNNNNNNNN---NNNNNNNNNNNGGAAGACCATCTTTTTAA--AGGACATTATCATACGCGGCTCTAGGTGCACAAGGTGCTAC-ATGGAAAAGTGAGCTA-GGTCTACGTTGTCATCCAAGAGGAACCCGTGAAGTTC-GTGACGGTGAGGGCGCAAACACGTCACCCCCGACTCCCGTCTGAGCTT---CCG-TGAA---AGTGTGGGTG--ACTTCCTCGACCCTGACC--TGA---ATCTGTCCGACTAGAAGGAAAGGGGAAGTGGTGGATAGAGTAAGGCTTTCATTCCTGTCGAGTACTGTTCAAGTAACGAAAATAGTAGGGGGACCAGGTACGACTCC------GAAGGAAACCTCAAGTCTAG---GAGCCACGACGGGTA-ACCCCTCGATTTTGACCGACTTACCCTGTTTCC-GGGTCTCTGGT-TCCTACGGTGCC-GCCGACGGCTCCGAGGACACCACCACAAAGAGGA---GAATTACTGCAACA---CGGACAAGTGTTACCGGAGCAAGTT---TTTGTCATACGAG--ACCTGTA-GTATCA-CCCGTTTCGTGAAACCTGTGCAGGTCCCGGACTGC-CTC---AGAAGACTGGCCAACAGGAG-TTCCAGCCTC-CT---GAAAGGTT---CAAG---T-AGTGGAAA-CCGCAAC--GCATTGA--ATGCTTTTTACCTACTCAGAGGTCA-AAGTCTGAATCTTGCA-GCGAACTCACCAAAAAGGG--AAAGGAATTTTTCGGACTCCGAGAGAGA---------GACACTCTACCCAGGC---AGATAT----TGTAGGTTTCTTTTCGGTCTTGGA---AACGGGGGTGGCCG-CTACTCTTTCAGTGA----CTCGAACCAGCCCCAACCTGGAG---CCAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTTTTTGACAGAGGCCTACAGGGAACGAACAGTTAACTTCTACACGATGGTATAAGGGACACGTCACTACCGGGTGAGAACATACACCCACTCTATGTAAGTAGTTTCGACAGGTGCAGTTCGAGTGAATACCATTACAGCGATAACGAGGCGACGGCCTAAGACGCTTCTCGAAAGGATCCGCTCGACCCCATCCGAAGTCCGTAGTCCTTAGCCCTCAGGTACGAGAGGCGGCCGGAGCGGAACGAGAAGTCGCACAAGACATACGAGAAGTCGTACCGACTCAGTGACCTCAACTCCGAGTGTCTCGAAGGGACCGAGTGCCCTGACGAAAACCTCCTG-------TCC-GAGAGTTCCGTTGAC--CTCCAAA-GAT--AAAGGACCGATCT-----TGCCG-A---G-CTAAGAAAATNNNNNNNNNNNNNNNN---NNNNNNN--NNNNNNNNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNN--NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNGACAGTCACCGGTGTCCGGAAAGGCAGACACTGGCAC-AGTCCGACCAG-TA---");
    // querys.push_back("-------------------------------------------------------------------------ATGACCAGCTTGAAACGATC-GCAGACAGAAAGGCCCATTGCCACTGACAGGGCCTCTG--TTGTCGGCACA-GACAGCACCC-----------------CCAAGGTCCACACTGATGA-CTTC-TA-CATGCG-GCGCTTCCGGTCT--CAAAATGGCAGCTTAGGGTCATCCATCATGGCTCCTGTGGGACCCCCA-CG---A-AGTGAAGGC--TCTCACCATATC-ACTTCGACTCCTGGAGTCCCAAAGATGGGG-GTCAGGGCC------CG----GATTGCAGATTGGC-CA-CCGAGAAAGGAAAATATCAAA-GAATCT---A-GCCGTTCAAGCCAGGAAATAGACA-CCTCTAG-TTGCCTTGAGAGTATGTCCTCCAAGAGCAGCCCCGTG----AGCCAG---GGCAGCTCTGTGAGCCT-CA---ACTCCAGC--GACTCGGCCATGCTGAAGAG---CA--TACA-GA-GCACGC--T-GAACAGCAAGAC-GAGGCCGGCGGAGAGCATGGACTCC-CGGTTCCTCATGCCTGA-AGCCTA-TCCTGGCTCCCCGAGGAAGGCACTTCGCCG--GATCCGTCAGCGGAGCAACAGTGACA-TCACCATAAGTGAGCTTGATGTGGACAG-CTTTG--ATGA-ATGTATCTCCCCGACGTACAAGACGG-GATCATCGCTGCATAGGGAGTATGGTAGCACGT-CCTCCATCGACAAGCAGGGCACGTC-CGGAGAGAGCTTCTTCGATCTCTTAAAGGGTTACAAGGAGGACAAGGCGGACC---GAGGCCCAACA-CCGACCAAGCTTAGTGACTTTCTC-ATTGCCGG--TGGGGGCAAGGGCTCTGGCTTC-TCCTTGGATGTCATCGA---TGGGCCCATCTCACAGA-GGGAGAACCTCAGACTCTTTAAGGACAGGGAAAAGCCACTCAA--GC--GCCGGTCT-AAGTCGGAA-ACCGGGG-ACT-CCTCCATTTTCC-GTAA-GTTGCGCAATGCCAAAGGA-GAT---GAGC---TTGGGAA-G---TCGTCGGACCTTGAGGATAACCGGTCAGAGGA-----CTCCTCCAGGCCCTGGACTTGTCCTCGGTGCTTCTCACACTACGATGTGCAGAGCATACTGTTTGACTTGAATGAAGCCATCGTGAGTAGGC---A-CAATGTGATCAAGAGGAGGAACACTACCACAGGAGCCTCAGCGGCGGCTGTGGCATCCTTGGTCTCAGGGCCCCTGTCAC-ATTCTGCCAGCTTCAGCTCCCCCATGGGC-AGCACAGAGGA-CCTGAACTCCAA-AGGAAA------CCTGGGTGTGGACCCGGGAGATGACAAAAGCAACGAGCTGGTTTTGAGCTGCCCA-TACTTTCGGAACGAGATAGGTGGGGAAGGGGAGAGACAGATCAGCCTGTCCAAGGCCA---GCTCCGGCTCCTTCAGCA-GCTGTGATGGTGCCTCTTTTGAAGCCACCCTGACCTCCCACTGCACAAACGCAGGTGTG-GCTGTCC-TCGAG---GTACCC-AAGGAGAACTTGGTGCTGCATCTCGACAGAGTGAAAAGGT-ACGTCGTGGAGCACGTGGACCTTGGTGCATACTACTACAGGAAGTTCTTCTACCAGAA--GGAACACTG-G-AACTATTT-TGGGGCTGATGAGAATCTTGGTCCAGTGGC-TGTGAGCATTCG-AAGGGAAAAA-CCTGAAGAAATGAAAGAGAATGAGGCCTTGA--CC---------ATGGACACCAGGGCTGTCTGGCTCATGACGCTG-AGAGGTTCAGTCCTGGAGGATGCCATTCCCTCAACAGCCAAGCACTCAACAGCCAGAGGGCTGCCTCTGAAGGAAGTGCTGGAGCATGTGGTACCCGAGCT-------CAACGTGCAGTGCCTGCGACTGGC--------------------------------------------------CTTCAACA--------------------------CACCCAAGGTCA---------------CAGAGCA-------------------GCTCATGAAGTTGGACGAGCA-AGG-----------------------------------------------------------------------------------------------------------------------GC-TGAACTTCCAGCAGAAAGTGGGCATCATGTACTGCAGAGCTGGGCAGAGCACC-GAGGAGGAGA-TGTACAACAACGAGGCCGCTGGCCCGGCCTTCGAAGAGTTCCTCCAGCTCCTGGGAGAGCGGGTTCGGCTCAAGGGATTTGAGAAGTATCGTGCACAGCTTGACACCAAAACTGACTCCACTGGAACTCACTCTCTGTA-CACAACGTACAAGGACTATGAGATCATGTTCCATGTCTCCACCATGCTACC--CTACACACCTAACAACAAGCAGCAG---------------------------CTCCTTCGGAAGCGACACATTGGAAATGACATCGTCACCATCGTTTTCCAAGAGCCTGGAGCTCAGCCTTTCAGTCCAAAAAACATCCGCTCCCATTTCCAGCA-CGTCTTCGTCATCGTCCGGGCACACAACCCCT-GCTCTGACAGTGTCTGCTACAGCGTGGCTGTTACCAGGTCCAGAGACGTGCCTTCCTTTGGGCCTCCCATCCCAAAAGGGGTCACTTTCCCCAAGTCACATGTGTT-CAGGGACTTCCTTCTG-GTCAA----AG-TGAT-TAATGC-GGAAAACGC-TGCTCACAAATCGGAGA-AGTTCCGAGCCATG-GCGACACGTACCCGCCAGGAATAC---AAGGACCTGGCGGAGAAGAACGTCA-CCAATACCCCTATCGACCCGTCTGGCAAGTTCCCCTTTATCTCCCTG-GCCTCCAAAAAGAAGGAGAAGTCTAAGCCATACCCAGGTGC-TGAGCTCAGCAGCATGGGGG-CCATTGTGTGGGCT-GTCCGCGCCAAGGA-CTA-CAGCACAGCC-CTA-GAGATTGACTGCCTCTTGGGGA-TCTCTAATGAGGTCGTCGTCCTCATCG--AGCAGGAGACGAAGAGCGTGGTGTTCAACTGTTCCTGCAGGGACGTCATCGGGTGGACCTCAGCCAATG-C-TAGCCTCA--AGATCTTCTATGAGCGT--GGAGAGTGTGTCTCGGTGGAGAGCTTCATGGGCAG---CGAGG-ACATCAA--GG-AGATGGTCAGGAGGTTGCAGTTCGTCTCAAAAGGC--TGTGAGTCAGTGGAGATGACTCTGCGAAGAAACGGGCTGGGAC-AGCTTGGTTTCCATGTCAACTACGAGGG---TATCGTGGCTGATGTAGAGCCCTATGGTTATGCGTGGCAGGCAGGGCTGAGACAAGGCAGCCGCCTGGTAGAGATCTGCAAGGTGGCGGTGGCCACCCTGAGCCATGAGCAGATGATTGACCTCCTGAGGACATCTGTCACAGTGAAGGTCGTCATCATCCCCCCGCACGACGACTGCACCCCGCGCAGN--------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNGCTATCTCCTGGGTCGGACATTTATGTGACGGTTTCATCCTTGGCTCTGGCAATATCGCAG---AGTCCAAACTCCCCCAACAATGTGTCTTCTTT---TGAGACAGGCTCTGGGGGTGGGACATACAATTTCAAATCCATGCCCAAAGAGTTTGGCATGACCCGTAGGTCCCCTGCCTCCATGGACCGGCAGAACACCCAGTCGGACCT---CGGTGGCAGTGGGAAGTCCACCCCCAGCTGGCCAAGGAGTGAGGACAGCATTGCCGACCAGATGGCTTACGGTTATAGAGGACCTCAGGATTTCAATTCTTTTGTCCTCGAGCAGCATGAATACACAGAGCCCACATGCCACCTCCCAGCAGTGTCCAAGGTTCTGCCAGCTTTCCGAGAGAGCCCCAGTGGGAGATTAATGCGGCAGGATCCAGTGGTTCATTTGTCTCCAAACAGACAAGGGCATTCCGACAGCCACTATTCGAGCCACTCCAGCAGCAACACCCTGTCCAGCAACGCCTCCAGCGCCCATAGCGATGAGAAGTG-GTTCGA---TGCAGACCGCACAGAGGCCGAGCTCAACAGCTACAGCTATCTGCAGGGCACATCGGCTGACAGTGGCATTGACACCACATCCTATGGCCCCAGCCACGGCAGTACTGCCTCGCTGGGGGCCACTACGTCATCACCCCGCGCG---GGGCCGGGCAAGGAGAAGGTGGCGCCACCCTGGCACAGCTCCAGCGATGCCATCTCC-ATGGCGGATCGGCCCCTGGACACGGAGGGCCATGGCATGGAGCGAAAGGCAGAGTCCTCCCTGAGCCTGGACCTCCACAGCAAGGGACCTGGTGGCTCAAGCCCGCTGGCCAGGGAGAACAGC-ACGTTCAGCATC-GGCGACGTAGCCTCCCACACCAGCACTGAGAGCTTCCGACATTTTGCCAGTCCTGTGGTTTTCACCAGCGCC---AGTTCACCCAAAGAGGAGCTGCACCCGGCCGCCTCCTCCCAGCTCGCCCCAGCCTT------CTCCTCCTCCTCCTCCTCCTCTTCTGGCCC-CAGGAGCTTCTACCCACGCCAGGGTGCGACTAGCAAATACCTGATTGGATGGAAGAAACCTGAAGGAACCATCAACTCTGTGGGGTTTATGGACACGAGAAAGCGTCATCAAAGTGATGGCAATGAGATGGCCCACACCAGGCTGCGGGCCTCCACCAGGGACCTCCGGGCGTCCCCCAAGCCAGCCTCCAAGTCTACCATCGAGGAGGATCTGAAGAAACTGATTGATCTCGAAAGCCCAACTCCCGAATCACAGAAGAACTTTAAG------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN------NTCCAGAGAGCCTCATTTTTTGCTGCTAGTGATGAAAACCATCGACCCCTGAGTGCCGCATCTAACGGTGATCAGCTTGAGGATCTGCCTCTGGCCCAGATGAAACCGTACAGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNN--NNNNNN---NNN-------------GAGAAGGAGGACAAGGC---CCACCTGCAGGCTGAAGTTCAGCACCTGCGGGAGCACAACCTGCGGCTGCAGGAGGAGTCGCAGAATGCCTCTGACAAGCTGAAGAAGTTCACCGAGTGGGTCTTCAACACCATTGACATGAGCTA---");
    // querys.push_back("-ATGACCAGCTTGAAGAGGTCGCAGACCGAAAGGCCTGTCACTGCTGACAGAGCCTCTGTTGTCAGCACAGATGGCACC-C---CCAAAGTTCACACTGATGATTTCTACATG-CGTCGCTTCCGGTCCCAAAATGGCAGCTTAGGATCATCAGTC--------------ATGGCTGCAGTAGGGCCCCCTCGAAGTGAAGGTCCTCACCATATAACCTCAACCCCCGGAGTCCCCAAGATGGGGGTTAGGGCA------AGAATTGCAGATTGGCCT-CCAAGAAAGGAAAATGTAAAAGAATCT---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCCTGTCCTCCAAAGGCAGTCCTGTGAGTCAGGGGAGTTCTGTTAGCCTCAATTCCAATGACTCAGCCATGCTGAAGAGCATCCAGAACACCCTGAAGAACAAGACGGGGCCAGCGGAGAGCATGGACTCGAGATTCCTCATGCCTGAAG-CCTACCCCAGCTCCCCCAGGAAAGCCCTTCGCCGAATCCGTCAGCGTAGCAACAGTGATATCACCATAAGTGAGCTTGATGTGGATAGTTTTGATGAATGTATCTCCCCAACCTACAAGTCGGGGCCATCATTGCACAGGGAATATGGTAGCACATCTTCA-ATCGACAAGCAGGGGACATCCGGAGAAAGCTTCTTCGGTTTGTTAAAGGGC-TATAAAGACGACAGAGCTGACC---GAGGTCCAACTCCAACCAAACTCAGTGACTTCCTCATC-ACTGGTGGGGGCAAGGGTTCTGGTTTCTCCTTGGATGTGATCGA---TGGA-CCC-ATCTCACAGAGAGAGAACCTCAGGCTATTCAAGGAAAGGGAAAAACCACTCAAGCGACGCTCTAAGTCTGAAACTGGAGACTCATCCATTTTTCGTAAATTACGCAATGCCAAAGGTGAA---GAAC---TCGGGAAG---TCATCAGACCTTGAAGACAACA--GATCAGAGGA---TTCTGTGAGGCCCTGGACATGTCCAAAGTGCTTTGCCCACTATGATGTCCAGAGCATATTGTTTGACTTGAATGAAGC-CATTATGAACAGGC---ATAATGTGATTAAGAGGAGAAACACCACAACGGGAGCTTCCGCAGCCG--CTGTGGCGTCCTTAGTCTCGGGACCTTTG-TCTCACTCAGCCAGCTTCAGCTCTCCTATGGGCAGCACAGAGGACCTCAATTCCAAAGGAAG------CCTCGGCATGGACCAGGGAGA-TGACAAGAGCAATGAACTCGTCATGAGCTGTCCGTACTTTCGGAATGAGATTGGGGGAGAAGGTGAGAGGAAAATCAGCCTGTCCAAGTCCAATTCTGGCTCATTCAGTGGGTGTGAGAGCACGTCCTTTGAGTCTGCCCTCAGCTCTCACTGCACAAACGCAGGAGTGGCGGTTCTCGAAGTGCCCAAGGAGAGCTTGATGCTGCATCTAGACAGAGTGAGACGGTACACTGTGGAGCACGTGGATCTTGGAGCATACTATTACAGGAAGTGCTTCTACCAAAAGGAACACTGGAACTATTTTGGGGCTGATGAAAACCTCGGTCCAGTGGCTGTGAGCATTCGAAGGGAAAAACCAGAAGACATGAAG-GAAAATGGATCTCCCTACAACTACCGAATAATTTTTAGAACTAGTGAGCTCATGACACTG-AGAGGGTCTGTCCTGGAGGATGCCATTCCCTCCACAGCCAAGCACTCGACAGCCAGAGGATTGCCGCTGAAAGAGGTGCTGGAACATGTGATCCCAGAGCT-------CAACGTGCAGTGCCTGCGTTTGGC--------------------------------------------------CTTCAACA---------------------------CGCCCAAA-GTCA--------------------CGGAGCA----------------------GCTTATGAAACT-GGACGAGCAAGG---------------------------------------------------------------------------------------------------------------------GCTGAACTATCAGCAGA-AAGTAGGCATCATGTACTGCAAAGCAGGCCAGAGCACCGAGGAGGAGATGTACAATAATGAGTCCGCTGGCCCGGCCTTTGAGGAATTTCTTCAGC-TACTGGGAGAGCGAGTTCGCCTAAAAGGA-TTTGAAAAGTATCGTGCGCAGCTCGACACCAAAACTG-ACTCCACTGGAACCCACTCTCTGTACACAA-CCTACAAAGACTATGAGATAATGTT-CCACGTCTCCACCATGCTGCCCTACACACCTAACAACAAGCAACAG---------------------------CTCCTGAGG-AAGCGGCACATTGGGAATGACATTGTGACAATAGTTTTCCAGGAGCCTGGAGCACAACCA-TTTAGCCCGAAAAACATCCGATCTCACTTTCAGCATGTTTTTGTCATTGTCCGGGCTCACAACCCCTGCACT-GAGAGTGT-CTGTT-ACAGTG-TGGCAGTCACCAGGTCCAGAGAT-GTGCCTTCTTTTGGACCTCC--CATCCCTA-AAGGGGTCAC--CTTCCCCAAGTCAAACGT-GTTCAGGGACTTCCTTTTGGCCAAAGTGATTAATGCAGAAAACGCTG-CTCATAAATCAGAGAAGTTCAGGGCCA---TGGCAACAAGGACCCGCCAGGAATACCTGAAAGATCTGGCAGAAAAGAATGTCACCAACACACCTATTGACCCT-TCTGGCAAGTTTCCATTTATTTCTCTGGCTTCTAAGAAGAAGGAAAAGTCTAAGCCTTACCCAGGAGCTGAGCTCAGTAGCATGGGGGCCATTGTGTGGGCAGTCCGGGCCAAAG-ACTACAACAAGGCCATGGAATTCGACTGCCTCCTTGGAA--TCTCCAATGAGTTCATCGTCCTCATTG-AGCAGGAAACAAAGAGCGTGGTTTTCAATTGCTCCTGCAGAGATGTGATAGGGTGGACTTCCAGCGACAGCAGCCTTAAGATCTTCTATGAGCGTGGAGAATGCATTTCTGTGGAGAGCTTCATGAGCAG---TGAAGATATCAAAGAA-ATTGTCAAAAGGCTGCAGTTTGTTTCAAAAGGTTGTGAATCAGTGGAGATGACGCTGCGAAGAAATGGGCTAGGCCAGCTTGGCTTCCATGTCAACTATGAGGGCATCGTGGCGGATGTAGAACCCTATGGCTATGCATGGCAAGCAGGTTTGAGGCAGGGCAGCC--GCCTGGTGGAGATCTGCAAGGTGGCAGTGGCCACCCTGA-GCCATGAACAGATGATCGATCTCTTGAGAACGTCAGTC--ACA-GTGAAAGTTGTCATTATCCCTCCCCATG---ATGACTGCACCCCTC-G-----GAGG------AGTTGCTCAGAAACCTACCG-CATGCCAGTGAT-------------GGAGTACAAGATGAATGAAGGCGTTTCCTATGAGTACAAGTTTC-CCTTCCGGAGT-AACAACAAA-T--GGCAGCGTAATGCCGGCAAGGGT------G-CCCACTCGCCCCAGG-TTCCATCGCAGCTGCAGA-GCCCCATGACTTCACGGCTGAATGCTG-GGAAGGGAGATGGAAAGATGCCCCCTCCAGAAAGAGCCGCCAACATCCCTCGAAGCATCTCCAGTGATGGGCGCCCACTGGAAA---GGAG-GCTGTCTCCGGGTTCGGACATCTATGTGACAGTCTCATCCATGGCTTTGGCGAGATCCCAG---TGCCGTAACTCCCCCAGCAACTTGTCATCGTCCAGCGAGACTGGCTCTGGGGGTGGCACTTACAGACAAAAATCCATGCCTGAAGGGTTTGGAGTGAGCCGAAGATCCCCAGCTTCCATCGACAGGCAAAACACCCAGTCGGATAT---TGGTGGCAGTGGAAAATCCACACCCAGCTGGCAGAGAAGTGAGGACAGCCTTG-CCGACCAGATG-----------------------------------------------------------------------------------GAGCCTA-CATGCCA--TCTCCCAG-CAGTATC--TAAGG----TACTGCCTGCTTTC----CGAGAGAGCC-C---CAGTGG--GAGGT--TGATGCGGCAGGATC-CC-GTGGTTCACTTGTCTCCAAACAAGCAAGGC-CAT-TCTGACAGC--CACTACTCCAGCCACT-CCAGCAG--CAACACTCTGTCCAGCAACGCCTCGAGTGCCCACAGTGACGAGA--AGTGG--TAC-GA---TGGAGACCGCACCG-AGTCCGACCTCAACA--GCTATAACTACCTACAGGGCAC-ATCTG-CTGACAGCGG-CATCGACACCGCCTCATACGGCCTCAGCCATGGCAGCACGGCCTCCCTGGGAGCCTCCACATCCTCACCTCGTTCA---GGGCCAGGCAAAGAGAAGGTGGCACCCCTGTG--GCA--CAGCTCCAGTG--AGGTGCTC--TCCCTG--GC-AG---ATCGGACCTTAGAGACTG-AGGGCCATGGCATGGACAGGAAAACAGAGTCCTCCCTG-AGCCTGGACATCCA--CAGCAAGAGCCAAGGCGGTTCAAGCCCGCTGACAAGGGAGAACAGCACCTTCAGCATAA---ATGATGCCACATCT-CACACCAGTACCATGAGCTCCCGACACT-CTGCCA-GCCCAGTGGTGTTCTCCAGTGCCAGAAGTTCACCCAAAGAGGAGCTTCACCCCACCA---CGTCCTCCCAGCTCGCACCTTCC-TT------CTCTTCCTCTTCCTCATCTTCCTCTGG--ACCT-AG-GACT-TTCTATCCCCGCCAGGGCGCCACTAGCAAATACCTG-ATTGGATGGAAAA--AGCCAGAAGGAA-CCATTAACT--CCGTGGGATTTATGGATACAAGA-AAG-CGACA-TCAGAG--TGATGGTA--ATG-AGATAGCCCACACTAGGCTT-CGAGCCTCAACCAGGGACCTTCGGGCATCCCCAAAG--CCAACTTCCAAGTCCACCATTGAGGAGGA---------TCTAAAGAAACTCATCGACCTTGAG----AGCCCAACTC-CTGAATCACAGA---AAAACTTCAAG-------TTCCATGGACTG------TCCTCCCCACAGTCCCCGTTCCCCAGCACCCCTACCTCCCGGCGGGCCCTGCACAGGACT--CTGTCTGATGAGA-GCAT------TTACAGCAGCCAGAGGGAGCATTTCTTCACCTCCAGG-GCCTCACTTCTAGACCAAGCCCTGCCCAATGATGTC-CTCT-TCAGCAGCACCTACCCATCTCTCCCCAAGTCGCTTCCACTGAGGAGGCCATCTTACACATTGGGAATGAAGTCGCTGCATGGAGAGTTCTCTGCTTCCGAC-AGCTCCCTTACCGACATCCAGGAGACCCGGAGGCAGCCTATTCCTGACCCTGGCCTGATGCCCCTGCCTGACACTGCTTCG--GACTT--GGACTG-GTCCAACCTAG---TAGATGCCGCCAA--AGCCTATGAG------GTC-CAGAGAGCCTCATTTTTTGCTGCTAGTGATGAAAACCATCGCCCCTTGAGTGCGGCCTCCAACAGCGACCA-GCTAGAGGAGCAGGCC---CTGGTCC--AGATGAAGTCGTAC---AGCAGTAAGGACTCCTCTCCCACTCTGGCTTCTAAAGTGGACCAGTTGGAAGGTATGCTGAAAATGCTTCGAGAAGATTTGAAG---AAGGTAAATGTTCTCCGAAAAGGAAGACAA---GGCCCACCTGCAGGCAGAGGTTGA---------------------------------------------------------------------------------------------------------");
    // querys.push_back("-ATGACCAGTTTGAAGCGGTCGCAGACCGAAAGACCTGTCACCGCTGACAGAGCCTCTGTTGTCAGCACAGATGGCGCC-C---CCAAAGTCCACACCGATGACTTCTACATG-CGTCGCTTCCGCTCTCAGAATGGCAGCCTAGGATCATCAGTC--------------ATGGCTGCAGTGGGGCCCCCTCGAAGTGAAGGCCCTCACCATATCACCTCAACCCCCGGGGTCCCCAAGATGGGGGTTAGGGCA------AGAATAGCAGATTGGCCT-CCGAGAAAGGAAAATGTAAAAGAATCC---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCCTGTCCTCCAAAGGCAGTCCTGTGAGTCAGGGGAGTTCTGTTAGCCTCAATTCCAATGACTCAGCCATGCTGAAGAGCATACAGAACACCCTGAAGAACAAGACAGGGCCAGCGGAGAGCATGGACTCCAGATTCCTCATGCCTGAAG-CCTACCCCAGTTCCCCCAGGAAAGCCCTTCGCAGAATTCGGCAGCGCAGCAACAGTGATATCACCATAAGTGAGCTTGATGTGGATAGCTTCGATGAATGTATCTCCCCAACCTACAAGTCGGGGCCATCATTGCACAGGGAATATGGCAGCACATCTTCA-ATCGACAAGCAGGGAACATCCGGAGACAGCTTCTTCGATTTGTTAAAGGGC-TATAAAGATGACAGATCTGACC---GAGGTCCAACTCCAACCAAACTCAGTGACTTCCTCATC-ACTGGTGGGGGCAAGGGTTCTGGTTTCTCCTTGGATGTGATCGA---TGGC-CCC-ATCTCACAGAGAGAGAACCTCAGGCTTTTCAAGGAAAGGGAAAAACCACTCAAGCGACGCTCTAAGTCTGAGACTGGAGACTCGTCCATTTTTCGTAAATTGCGCAATGCCAAAGGTGAA---GAAC---TCGGGAAA---TCATCAGACCTTGAAGACAACA--GATCAGAAGA---TTCTGTGAGGCCCTGGACATGTCCAAAGTGCTTTGCCCACTATGATGTCCAGAGCATATTGTTTGACTTGAATGAAGC-CATTATGAACAGAC---ATAATGTGATTAAGAGGAGAAACACCACAACAGGAGCTTCGGCGGCTG--CGGTGGCATCCTTGGTCTCCGGACCTCTG-TCTCACTCAGCCAGCTTCAGCTCTCCCATGGGCAGCACAGAGGACCTCAACTCCAAAGGAAG------CCTTGGCATGGACCAGGGAGA-TGACAAGAGCAATGAACTCGTCATGAGCTGTCCGTATTTTCGGAATGAGATTGGGGGAGAAGGTGAGAGGAAGATCAGCCTGTCCAAGTCGAATTCTGGCTCATTTAGTGGGTGTGAGAGCACATCCTTTGAGTCTGCCCTCAGCTCTCACTGCACCAACGCGGGCGTGGCAGTTCTCGAAGTGCCCAAGGAAAGCTTGATGCTGCATCTGGACAGGGTGAAAAGGTACACCGTGGAACACGTGGATCTTGGCGCATACTATTACAGGAAGTTCTTCTACCAGAAGGAACACTGGAACTATTTTGGGGCTGATGAGAACCTCGGTCCAGTGGCTGTGAGCATTCGAAGGGAAAAACCAGAAGACATGAAG-GAAAACGGATCTCCATACAACTACCGAATAATATTCAGGACTAGTGAGCTCATGACGCTG-AGGGGGTCTGTCCTGGAGGATGCCATTCCCTCCACGGCCAAGCACTCGACAGCCAGGGGATTGCCTCTGAAAGAGGTGCTGGAACACGTGATCCCAGAGCT-------CAACGTGCAGTGCCTGCGCTTGGC--------------------------------------------------CTTCAACA---------------------------CACCCAAA-GTCA--------------------CAGAGCA----------------------GCTCATGAAACT-GGACGAGCAAGG---------------------------------------------------------------------------------------------------------------------GCTGAACTATCAGCAGA-AAGTAGGCATCATGTACTGCAAAGCAGGCCAGAGCACGGAGGAGGAGATGTACAACAACGAGTCTGCAGGCCCAGCCTTTGAGGAGTTCCTTCAGC-TGCTGGGGGAACGAGTCCGGCTAAAAGGA-TTCGAGAAGTATCGTGCGCAGCTTGACACCAAAACTG-ACTCCACTGGAACCCACTCTCTGTACACAA-CCTACAAAGACTATGAGATAATGTT-CCACGTCTCCACCATGCTGCCCTACACACCTAACAACAAGCAACAG---------------------------CTCCTGAGG-AAGCGGCACATTGGGAATGACATTGTGACAATAGTTTTCCAAGAGCCTGGAGCACAACCA-TTCAGCCCGAAAAACATCCGGTCTCACTTTCAGCATGTTTTTGTCATTGTCCGGGCTCACAACCCTTGCACT-GAGAGTGT-CTGTT-ACAGTG-TGGCAGTCACCAGGTCCAGAGAT-GTACCTTCTTTTGGACCTCC--CATCCCTA-AAGGGGTCAC--CTTCCCCAAGTCAAATGT-GTTCAGGGACTTCCTTTTGGCCAAAGTGATAAATGCAGAAAATGCTG-CTCATAAATCAGAGAAGTTCCGGGCCA---TGGCGACAAGGACCCGCCAGGAATACCTGAAAGATCTGGCAGAAAAGAATGTCACCAACACACCTATTGACCCT-TCTGGCAAGTTTCCATTTATTTCTCTGGCCTCCAAGAAGAAGGAAAAGTCTAAGCCTTATCCAGGAGCTGAGCTCAGTAGCATGGGGGCCATTGTGTGGGCTGTCCGGGCCAAAG-ACTACAACAAGGCCATGGAGTTCGACTGCCTCCTTGGGA--TCTCCAGCGAGTTCATCGTCCTCATTG-AGCAGGAGACAAAGAGTGTGGCTTTCAATTGCTCCTGCAGAGATGTCATAGGGTGGACTTCCAGCGACACCAGCCTCAAAATCTTCTATGAGCGGGGAGAATGTGTGTCGGTGGAGAGCTTCATTAGCGG---TGAAGATATCAAAGAA-ATTGTCAGAAGGCTGCAGTTTGTTTCAAAAGGTTGTGAATCTGTGGAAATGACTCTGCGAAGAAATGGGCTGGGGCAGCTTGGCTTCCATGTCAACTATGAGGGCATTGTGGCGGATGTAGAACCCTACGGCTACGCATGGCAAGCAGGGCTGAGGCAGGGCAGCC--GCCTGGTGGAGATCTGCAAGGTAGCAGTGGCCACCCTGA-GCCATGAACAGATGATCGATCTCTTGAGAACATCAGTC--ACA-GTGAAGGTTGTCATTATCCCTCCCCATG---ACGACTGCACCCCAC-G-----GAGG------AGTTGCTCAGAAACCTACCG-CATGCCAGTGAT-------------GGAGTACCAGATGAATGAAGGCATTTCCTACGAGTTCAAGTTTC-CCTTCCGGAAT-AATAACAAA-T--GGCAGCGGAATGCCAGCAAGGGT------G-CTCATTCGCCCCAGG-TTCCATCTCAGCTGCAGA-GTCCCATGACCTCACGACTGAATGCTG-GGAAGGGAGATGGGAAAATGCCCCCTCCAGAAAGAGCTGCCAACATCCCTCGAAGCATCTCCAGTGACGGGCGCCCACTGGAAA---GGAG-GCTGTCTCCTGGTTCGGACATCTATGTGACAGTCTCATCCATGGCTTTGGCGAGATCCCAG---TGCCGTAACTCTCCCAGCAACTTGTCTTCGTCCAGTGAGACTGGCTCTGGAGGTGGTACCTACAGACAAAAATCCATGCCCGAAGGGTTTGGGGTGAGCCGAAGATCCCCAGCTTCCATCGACAGGCAGAACACCCAGTCGGATAT---AAGTGGCAGTGGAAAATCCACTCCCAGCTGGCAGAGAAGTGAGGACAGCCTTG-CCGACCAGATG-----------------------------------------------------------------------------------GAGCCGA-CGTGCCA--TCTCCCAG-CAGTATC--GAAGG----TACTGCCTGCTTTC----CGAGAGAGCC-C---CAGTGG--GAGAT--TGATGCGGCAGGATC-CA-GTGGTTCACTTGTCTCCAAACAAACAAGGC-CAT-TCTGACAGC--CACTACTCCAGCCACT-CCAGCAG--CAACACGCTCTCCAGCAACGCCTCGAGTGCACACAGTGACGAGA--AGTGG--TAC-GA---TGGGGACCGCACGG-AGTCCGACCTCAACA--GCTACAACTACCTACAGGGCAC-GTCTG-CCGACAGCGG-GATTGACACCGCCTCCTACGGCCCCAGCCATGGCAGCACGGCCTCCCTGGGGGCCTCCACATCCTCACCTCGTTCA---GGGCCAGGCAAAGAAAAGGTGGCTCCCCTGTG--GCA--CAGCTCCAGTG--AAGTGCTC--TCCCTG--GC-AG---ATCGGACCTTAGAGACTG-AGGGCCACGGCATGGACAGGAAAGCAGAGTCCTCCCTG-AGCCTGGACATCCA--CAGCAAGAGCCAGGGCGGGTCAAGCCCGCTGAGCAGGGAGAACAGCACCTTCAGCATAA---ATGATGCTGCGTCC-CACACCAGTACCATGAGCTCCCGACACT-CTGCCA-GCCCAGTGGTATTCTCCAGTGCCAGAAGTTCCCCCAAAGAGGAGCTTCACCCCACCG---CATCCTCCCAGCTCGCACCGTCC-TT------TTCCTCTTCTTCCTCATCCTCCTCTGG--ACCT-AG-GACT-TTCTACCCTCGCCAGGGCGCCACTAGCAAATATCTG-ATTGGATGGAAAA--AGCCAGAAGGAA-CCATTAACT--CCGTGGGATTTATGGACACACGA-AAG-CGACA-TCAGAG--TGATGGCA--ATG-AGATAGCCCACACTAGGCTT-CGAGCCTCAACCAGGGACCTGCAGGCATCCCCAAAG--CCGACCTCCAAGTCTACCATTGAGGAAGA---------TCTAAAGAAACTCATCGACCTTGAG----AGCCCAACTC-CCGAATCCCAGA---AGAATTTCAAG-------TTCCATGCACTG------TCCTCCCCGCAGTCCCCGTTCCCCACTACCCCTACCTCCCGGCGGGCCCTGCACAGGACT--CTGTCAGATGAGA-GCAT------TTACAGCAGCCAGAGGGAGCATTTCTTCACCTCCAGG-GCTTCGCTTCTAGACCAAGCCCTGCCCAACGATGTC-CTCT-TCAGCAGTACCTACCCATCTCTCCCCAAGTCACTTCCACTGAGGAGGCCATCTTACACGTTGGGAATGAAGTCATTGCATGGAGAGTTCTCTGCCTCGGAC-AGCTCCCTCACCGACATCCAGGAGACCCGAAGGCAGCCTATCCCTGACCCTGGCCTGATGCCCCTGCCTGATGCAGCTTCA--GATTT--GGACTG-GTCCAACCTAG---TAGATGCCGCCAA--AGCCTATGAG------GTC-CAGAGAGCCTCATTTTTTGCTGCTAGTGATGAAAACCATCGCCCCCTGAGCGCGGCCTCCAACAGTGACCA-GCTGGAGGAGCAGGCC---CTGGTCC--AGATGAAGTCCTAC---AGCAGTAAGGACCCCTCTCCCACTCTGGCTTCTAAGGTGGACCAGCTGGAAGGTATGCTGAAAATGCTTCGAGAAGATTTGAAG---AAGGTAAATGCTCTCAGAAAAAGAAGACAA---GGCCCAGCTGCAGGCGGAAGTTGA---------------------------------------------------------------------------------------------------------");
    // querys.push_back("-GATCGAGTATAGATACCATAACTTCTGGGTGAGACACTTGAAGAAGTCGAACAGTCTCCGCAAGACACTGAGGAGGAC-GTTGGAGTCCAACAGGAGGGCGTTCACGACGTGAAGTCGGACTTCCACCCGAAACAGAAGAAAAAG---------------------------GAA---GAAGTCTAGAAGAGCTTCGTA-GAAGTC-GTATGGAAGGTTGA-CCAGGTGAAATCTTCGGTCTCACTCTCTCCTTAGGAATGACGACGACATACG-GA-AGTAGACGTGGTCTCGGACAAGTA-GTTCGACCAGTGACAACCTGCGGCGTGAGTCC-CCCGCTACCAAAAGTAGTGATCGTCGTTTTTTACTCC-GAGAGACCTG------GAG-TATCCGAAACCGTCGTAGATGGTTCAACCTGGTCAG-GTTTAG-ACGCCGTCACAGTCCGTCCCCGTAGTTCGGTCCCAGTCCGTACCCGA--CGGAAGCACAGAGGACCTACAGCCACTCCCTC--GACAGACTCCGGCTCTTGAGGGGTACGTTACTGAAGTAAGGATTGCA-CATACTGCCGGA-GGCGTCACCGTC---GAACCCCTCTCTAC-CCATGCACGAC-GACTTTTCCTGTAGTAACCCGTCGCGAACCAGGTCTTTACTTCGGGAACTCTACTTCTTCACGAGGGAG-ACCGACGACATTTACGAGAGCAGACTGTCACAGGACACGTTTCGAGCGGCGCTCTATCCGCAACCCCCTTTTCCTCTGACTCCCCTCCTGTCACGTACCTT------GAATTTTAAGAAGACACTAAGTCTTCAACC-CGA-AAGTTCCAGCTACTCAAAAAAATCTAGGAGAAGTTACCACCTGAACCTC-CAACCAAACCCCCTACGGGCCTCCAGGGCCCAACTCCGTGCGTCGGATCACACCCTATAGAGTAACGGTAGTGAGACTACTGCGAAAGAACACA-GGTATTTAGGGTGCCTCA--ATTATCAAGGAAGCCCAAAAAAG--GTAGGTTAGTCC-ATGAACGACCATCGAGGGACTGCCCCCATCTTACAAGAT-CCTGGTCTCCTTCTCCTCCTCCTCCTTCTC------TTCCTT-CCACGTTCGACACTCCTCCGCCGTCCTACTTCAAGA--AGAAATCCACTTGAAGACCGTGACCTCT--TTTGTTGACCCGACCGTCTCACGGCCCTCGAGTACCATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCGAA-ATGCCAGGTACGGCACCGAGAGACAG-AGATTTCAGGC-AAGACAGTACCTCT-CGTGAAGTGACCTCGATACGGTGTCCCCCCCGTGAAAGAGGAACGGACCGGG---ACTTGCCCCACTACTCCACCGTCGGGGATCCCTTCGACACGACGGTACCGACCCCGGTATCCTCCACCAC-AGCTACGGTGACAGTCGTCTCCACGGGACTTCTATCAACATCGATAACTCAAGTCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGTAGTGA-CACCC---GTGAGCTACGTAACGACCTCTCC-CATAACGATGACCTCACCGAGCTCATCACCGATAGTCGTACGGGAACAAACAAACCTCTGTTTACCTGGTGACCT--AGGACGGCA-TAATTAGA-GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTATG-ACGACCCTCT-ACCGTACATCCAANNNNN------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATAACCCAGCCGTTACGATAG-GAGTGAAGAAACGGTCGACCCACACCTAAA-AGGTGACGG---GGGT-----TCCAGACTGACCCACAAGA---CGGACAGGTACCTTCGACCACT-AGACGCTGAGTGAGGTTTGGGAAGTCCGTACCTAAAGACAGACATCCACGACGGGGGTCTCGGACAGAGCGACCTACTTCTGTCCAACGATCCCCTCAATGCTGT---GACTCTAGAACGATTTCGGTACCTACTCTGTCAGTGTATCTACAGGCTTGGGCCTCTG-TCCGCGG---AGAGATCACCCGCAGGCAGCGACCTCTACGAAGCTCCT--TACAACCGACGAGAAAGACCTCCTCCGTAAAAGGGTAGTGGAAAGGGCCATAAGTAGGCCCTCCAGTACCCTGAG-ACATGGACCCTACCCT-GAACTCCACT-TACTCC------TGGAAACAATCGTAAGGCGACGGT-AAACAATAAT--AAAGCCTTCCCTTTGAACTTGAGTATACTTTGCGGTAGTAAGTAAAATATGAGG-------------TAGTGACCGTACGC-CATCCAAAGGCTCGTTGA------GGAAGCACCCCACGTCAGTAGTACCCCCCCTTACTACTGTTGAAAGTGACACTGCCTACAGGA--GTCCTCTAGTTAGTAGACGAGTACCGAGTCCCACCGGTGTCGGTGGAACGTCTAGAGGT-GGTCCG-CCGACGGGACGGAGTCGGGACGG-ACGGTACGTATCGGTATGCC--AAGATGTAGACGGTGTTACGGGAGCATCAAGTGTACTTTCGGTTCGACAGGATTTGGTAAAGAAGCTTCTCAGTAAAGGTGACTAAGTGTTGGAAAACTTTGTTTGACGTCGGCA---AACTGTTAGAGAAACTATAGAAGT---AACGATTACTTTGAGAGATGGCTGTGCGTGAGAGGAGCGAGTATCTT-CTAAAACTGTGACCACAGTCAACTCCAGGTGGGATAGTGTA-GAGACGTCCTTGTTAACTTTTGGTGTGAGAAACAAAGGACAAGTTATTCCTGCTACTTGAGTAACCTCTAGGGCTCTTCTGTCAGATAAAGGTA-TCGGAACAATATCAGGAA--CCGGGCCTGACGGGTGTGTTACCGGGG-GTACGATGACTCGAGTCGAGGACCTATACCGAATCTGAAAAG--GAAGAAGAATCTCCGGTC-TCTCTACTTGCCTTTAAACGGTCTTCCCAGTTATCCCCACAACCACTGT-AAGAAAAGACGGTTTAG-AAAGTCCATAAGGACCGC-TCAGGCTCAACGGTATCGGGCCTTGAAAAGACTAAATACTCGTCGTAAGAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGACT-TGTGTAAAC--TGAACCCTTTACACTGG--GGAAAACCCTAT-CCTCCAGGTTTCCTTCCGTGTAGAGACCTGGACCAC--TGACGGT-GTGACATTGTCTGTGACAGGCACGTTCCCAA-CACTCGAGCCTGTTATTGCTTTTGTACGACTTTCACTCTGGCC-TACAAAAACCCCGATTTACCGACACGAG---GCCCGAGAACCTTTT-GTTAACAATGTTACAGTAACGGTTACACTGCGAAGG-CATCCTC---------------------------GACAACGAACAACAATCCACATATACCGTCGTATCATCTTTGTACCT-TGTATTAAAGT--ATCAGAAACAT-ACGGCACATGTCTCTTAC-CCAAGGTCACCTCAGTCAAAACCATAGTTCAACACGTGCTATGAAGAGTTTAGGAAAATCGGCTTGAGAAAGAGGGTTGTCAACTTCCTTAAGAAGTTTCCGAC-CCGGACGACTAAGTAACAACATGTAGAGAAGAAGTCATGAGACCGGTCGAAACGTCATGTACTACGGA-TGAAAGACGACTATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACGAGTAGGTCAAAGTACTCG-------------------ACGAGGC---------------ACTGAAACCCAC--------------------------ACAACT-TC--------------------------------------------------CGG------TCGGCGTCCGT-GACCTGTAAC-------TCGAGACCCTAGTGTACAAGGTCTTGGAGAAAGTCTCCGTCAGG---AGACCGACAACTTACGAAGCGACACCTCC-CTTACCGTAG--AAGGTCCTGTCTTGGGGA-ATCACA-GTACTCGAGTGATCAAGATTTTTAATAAGC-TAT-CAACATACC-TCTAGGTAAAAGAAAGT-AAAGAAGACCAAAAAGGGAAGCTTACGAATGTCGGTGGCCTGGTTCTAAGAGT--AGTCG--GGG-TTTTATCAAGGTCACAAGGAAGACCATCTTCTTAA--AGGATATTATCATTCGCGGCTCTAGGTGTACGAGGTGTT-ACATGGAAAAGTGAGATA-GATCTACGTTGTGGTTCAAAAGGAACCCCTGAAG-TTCATGACGGTGAGGGCGTAAACACGTTACCCTCGACTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTCGGTG--ATTTGCTCGGTCTTAA--ACTAA---AACTGTCCGACTAAAAGGAAAGTGGAAGAGGTGGATAGAGTAAGGCTTTTATGCCTGTTGAGTAATGTTCAAGTAACGAAAATAGTAGAGGGACCAGGTACGGTTCC------GAAGGGAACCTTAAGTCCAG---GAGACACGACGGGT-AACCACTCGACTTTGACCGGCTTACCCTGTTTCCA-GGTCTCTGGT-TCCTGCGGTGTC-GTCGACGGCTCCGAGGACACCACCACAAAGAGG---AGAATTATTGCAACA---CGGCCAAGTATTAGCGGAGTAAGTTTAGTTTATTATACGAG--ACCTGTA-GTAT-TACACGTTTCGTGAAACCTGTTCAGGTTCCGGACTG-TCTC---AGAAGACTAGCCAACAGAAG-TTCCAGGCT-GCT---AAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCATTAA--ATGCTTTTTACCTACTCAGAGGTC-AAAGTCTGAACCTGGCAGC-GAACTCACCAAAAAGG--GAAAGGAATTTTTCGGACTCTAAGAGAGA---------GACACTTTACCCGGGT---AGTTAT----TGTAGGTTTCTCTTTGGTCTTGG---GAACGGGGGTGGTCGT-TATTCTTTCAGTG----ACTCGAACCAACCTCAACCGGGTG---CTAGTCTAAACAGTAGAAACATCGGGAAATTGTTTAGTTTTTT---AAGAGGTCT------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAATCG--TACCG--ACTCAG-TAATCTTAACT---CCGATTGTCTTGA--AGG---TGAGTGTCCTGACGAAAACCTCCTGTACGAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCAA---TCTGAGAAAATACAAAAGGAAAGAACC---TCCGGTT--AGACGTTAGGA------ACGAGCTTGGGGGTAGAAACCCTGAGGTCCTCACCTCCAATATACTACTCTTGGAAG-TGAAGCACCCCCGGGATGCCCT--CGGTACTGACTACTAGGATTCG-ACGGTAAAACCC---TGGCCTTCGCCGCGTACATCTTTAGTAGGCACACCTGGAACC---CTCCCGGTAGGCACGGTTGTTGTCTCCAGGATAGTCACCGTTGTCCGGAAAGACAGACACTAGC-AAAGTTCGAACAGTA----");
    // querys.push_back("----CGAGTACAGTTACCACAACTTCTGGGTGAGGCACTTAAAGAAGTCGAACAGACTCCGCAAGACCCTGAGGAGAAC-GTCGGCGTCTAACAGAAGAGCTTCTACGACTTGGAGACGGACGTCCACCCGAAACAGAAGGAAAAG---------------------------GAA---GAAGTCTAGGAGTGCTTCGTAG-AAGTT-GTAGGGAAGGTTGACT-AGGTGAAATCTTCGGTCTCACCCACTTCTCAGGAATGACGACGACATACT-GA-AGTAGACTTGGTTTCGGACAAGGA-GTTCAACTAGTGACAACCTACGACGTGAGTTC-CCCGCTACCAAAAGTAGTGATCGTCGTTTTTTGCTCC-GAGCGACCTG------GAGT-ATCCGAAACCGTCGTAGCTGGTTCAATCTGGTTAG-GTTTAGA-CGTCGTCACAGTCCGTCCCCGTAGTCCGGTCCTAGCCCGTACCCGA--CGGAAGCTCACAGGACCTACAGTCAGTCCCTC--GACAGCCTCCGACTCTTGAGAGGTACTTCACTGAAGTAAGGATTCCA-CATACTTCCGGA-GGCGTCTCCCTCTCTGAACCCCTCCCTTC-CCATACACGATG-ACTTCTCCTGTAGTAACCCGTCCCGAACCAGGTCGTCACTTCGGGAGCTTCACTTCTTCACGAGAGAG-ACCGGTGACATCTATGAGAGTAGACTGTTGCAGGCCACGTTCCGGGCGGAGCTGCACCCCCACCATCCCTTTCCTCTGACTCCCCTTCTCTCGCCGACCTT------GAATTTTAAGAAGACACTAAGTCGTCAACC-CGA-AAGTTCTAGTTATTCGAAGAAATCTAGGAGGAGTTACCACCTAAACCTC-CAACCAAAACCCCTACGGGCCTCCAGGGATCAACTTCATGCATCAGACCACACCCTATAAAGTAACGGTAGTGAGACTACCGCGAAAGAACACA-GGTATTTAGGGTGTCTCA--ATTACCAAGGAAGTCCAAAAAAG--GTAGGTTAGTCCA-TGAACGATCATCGCGGGACCGCTCCCATCTTACAGGAT-CCTGGTCTCCTCCTGCTCCTCCTCCTC---------TTTCTT-CCACACTCGACACTCCTCCACCGTCCTACGTCGAGA--AGAAAACCACTTGAAGACCGTGACCTCT--TATGTTGACCCGACCGTCTCACGGCCCTCGAGTACCATGACCATACCCTTCGTCGTAGTAACTATGA-CTTTCACGACAAGAGGGACCAGTCTCCCGAGCTCGGTGGGACCGAGAACGACACCTACAGGTCCGAGTCGCTCCTGAGACGGAAG-GCTAGGTACGGTACCGAGAGACAG-AGATTCCAGGCTA-GCTGGTACCTGT-CGTGAAGTGACCTCGATACGGTGTCCCCATGGTGAAAGAGGAACGGACCAGG---GCTCGCTCCACTACTCCACCGTCGGGGGTCCCTCTGGCTCGACGATACCGACCCCGGTATTCTCCACCAT-AGCTAAGGTGACAGTCGACTCCAAGGGACGTCTATCAATATCGACAACTCAAGTCTAAGACACGCCAGGGGT---AGCATGGTGAAG-AGTAGTGA-CACCCG---TGAGCTACGTAACGACCTCTCC-CATAACGATGACCTCACCGAACTCATCACCGATAGTCTCACAGGAACAAACAAACCTCTGTTTACCTGGTGACCT--AGGACGGCGT-AATTAGA-GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTATG-ACGACCCTCT-ACCGTACATCCGAGACATATAAGTACGACGA--GCTCCTGTTTTCTTAACTTTAGGACTCCAGGAGATATTGACATTCGGTAGACCAGTCGT----TACGGTAG-GAGTGAAGAAACGGTCGACCCACACCTGAA-AGGTGACGG---TGGT-----TACAGACGGACCCACAAAA---CAGACAGTTATCTCCGACCTCT-GGAAGCCGAATGAGGTTTGGGAAGCCCGTACCTGAAAACGGACATTCATGGGGGTGGTCTGGGTCAGAGCGACCTACTTCTGTCCAACGATCCCCTCAATGCTGT---GACCCTAGAACGATTTCGGTACCTACTATGACAGTGTATCTACAGGCTTGGTCCTCTG-TCCGCGG---CAAGATCACCCGCGGGTAGTGACCTCTATGAAGCTCCT--TACAATCGCTGAGAAAGACCCCCTCCATAGAAGGGTAGATGAAAGGGTCGCAAGTTGGAACTCCATTACCCTGAG-ACGTGGACCCGACCTT-GAACTCCACG-TACTCC------GGGGAATGACCGTAAGGCGACGGT-AAACAACAAC--AAAGCCTTCCCTTTGAATTTAAGTATACTTTGTGGAAGTAAGTAAAACATAAGG-------------TAGTGACCGTACGC-CATCCAAAGTCTCGTTGA------GGAGGCACCCCACGTCAGTAGTACCCCCCCCTACTACTGTTGAAAGTGTCACTGTCTACAGGA--GTCCTCTAGCTAGTAGACGAGTACTGAGTCCCACCGGTGTCGGTGGAACGTCTAAAGGT-GGTCAG-CCGAAGGAACGGAGTCGGGACGG-ACGGTCCGTATCGGTATCCC--GAGGTGTAGACGCTGTTACGGGAGTATCAACTGTACCTTAGGATCGACAGGATCTGGTAAAGAAGCGTCTCAGTAGAGGTGACTAAGTGTTGGAAAACTTTGTTTGACGTTGGAA---AACTGTTAAAGAAACTATAGGAGT---AACGATTACTTTGAGAGGTGACTTTGTGTGAGGGGAGCGAGTATCTT-CTAAAACTCCGACCATAGTCAACTTCAGGTTGGATATTGTAG-AGACGTCCTTGTTAACTTTTGGTGCGAGAAACAAAGGACGAGTTACTCCTGTTACTTGAGTGACCTCTGTGGGTCTTCCGTCAGATAGAGGTAT-CGACACAACATCAGAAA--CCGGGCCTGCCGGGTGTGCTACCGGGG-GTACGATGATTCGAGTCGAGGACCTATACCGAATCTAAAAAG--GAAGAAGAACCTCCGGTC-TCTTTATTTACCCTTGAACGGACTTCCCAGTTATCCCCACAACCATTGT-AAGAAAAGACGATCCAG-AAAGTCCATAAGAACCGC-TCAAGCTCATCGGTACCGGGCCTTGAAAAGACTAAATACTCGTCGTAAAAGACGTAATTAGTGAAACCGGTTTTCTTTCAGGGACT-TATGCAAACTA--AACCCCTTTCACTGA--GGAAATCCCTAC-CCTCCAGGTTTCCTTCCGTGTAGAGATCTGGACCAC--TGTCGGTG-TGACATTGTCTGTGGTAGTCACGTCCCCAA-CACTCGGGACTGTTACTGCTTTTGTACGACTTTTACCCTAGAC-TACAAAAAACCCGACTTACCGACACGAG---GTCCGAGAACCTTTT-GTTAACAATGTTACAGTAAAGGTTACACGGCGAAGG-AATCCTC---------------------------GACAACGAACAACAATCCACACATGCCGTTGTACCATCTTTGTACCT-TGTACTAAAGT--ATCAGAAATAT-ACAACACATGTCTCTTAC-TCAAGGTCATCTCAGTCAAAACCACAGTTCGACACGTGCTATGAAGAGTTTAGGAAAGTCAGCTTGAGCAAGAGGGTTATCGACTTCCTTAAGGAGTTTCCGACCC-GGACGACTAAGTAACAATATGTAGAGAAGAAGTCACGAGACAGGTCGAAACGTCATGTACTACGG-TTGAAAAACGACTATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCAAAGTACTCG-------------------ACGAGCC---------------ACTGAAACCCTC--------------------------ACAACT-TT--------------------------------------------------CGG------TTGGCGTTCGT-GACCTGTAAC-------TCGAGTCCCTAGTGTACGAGGTCGTGAAGAAAGTCGCCGTTGGG---AGACCGTCAACTCACGAAACGACAACTTC-CTTACCGTAG--AAGATTCTGTCTGGGAGA-GTCACA-GTACTCGAGTGATCAAGATTTTTAATAAGC-TAT-CAACATACC-TCTAGGTAAAAGAAAGT-AAAGAAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGTTCTAAGAGT--AGTCG--GGG-TTTTATCAAGGTCACAAGAAAGACCATTTTTTTAA--AGGACATTATCATACGTGGTTCTAGGTGTACAAGGTGTT-ACATAGAAAAATGAGATA-GGTCCACGTTGTGATTCAAAAGGAACCCGTGAAG-TTCATGACGGTGAGGGCGTAAACACGTTACCCTCCATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTGGGTG--ATTTCCTCGGGCTTAA--ACTAA---AACTGTCCGACTAAAAGGAAAGTGGAAGAGGTGGATAGAGTAAGGCTTTTATACCTGTCGAGTAGTGTTCAAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTTAAGTCCAG---GAGGCACGACGGGT-AACCTCTCGACTTTGATCGACTTACTCTGTTTCCA-GGTCTCTGGT-TCCTACGGTGTC-GACGACGACTTCGAGGACATCACCACAAAGAGG---CGAATTACTGTAACA---CGGACAAGTATTAACGGAGTAAGTTTAGTTTATTATATGAG--ACTTGTA-GTAT-CACCCGTTTCGTGAAACCTGTACAGGTCCCAGACTG-TCTC---AGAAGACTAGCCAACAGAAG-TTCCAGTCT-ACT---AAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCATTAA--ATGCTTTTTACCTACTTAGAGGTC-AAAGCCTGAACCTTGCAGC-AAACTCACCAAAGAGG--GAAAGGAATTTTTCGGATTCCAAAAGAGA---------GACACTCTACCCGGGT---AGTTAT----TGTAGGTTTCTTTTTGGTCTTGG---GAACGAGGGTGGTCGT-TACTCTTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CCAGTCTAAATAGTAGAAATATCGAGAAATTGTTTAGTTTTTTCGAAAGAGGTCT------ACAAGGAACGAATAGTTAACTTCTACACGATGGTATAAGGGACACGTTACTACCCGGTCAGAACATACACCCACTCTATGTAAGTAGTTTCGATAGGTGCAGTTC-AAGTGAATACCACTATAGTGACAATGAGGCG-ACGGCATAAGACGCTTCTCGAAAGGATCCCCTCGACCCTATCCGAAGTCCGTATTCTTTAGACCTCAGGTACAAGAGGCTACCAGATCAAAACAAAAAGTCA--CACAAGACATATGAGAAGTCG--TACCG--ACTCAG-TAACCTTAACT---CCGATTGTCTTGA--AGGGACTGAGTGTCCTGATGAAAACCTCCTGTCCGAGAGTTCCGTTGGACTCCAAAGATAAAGGACCGAACTTGCCGA---CCTAAGGAAATACAAAAGGAAAGATCC---TCCGGTT--AGACGTTAGGA------ACGGGATTGGGGGTAAAAACCCTGTGGTCCTCATCTCCAATATACCACTCTTGGAAG-TGAAGCTCCCCCGGGATGTCCT--CGGTATTGACTACTAGGATTTGA-TGGTAAAACCC---TGGCCTTCGCTGCGTACATCTTTAGTAGTCACACCTGAAACC---CCCATGGTAGACACGGTTGTTGTCTCCGGGACAGTCACCGTTGTCCGGAAAGACAGACACTGGC-AAAGTTCGACCAGTA----");
    // querys.push_back("-GATCGAGTACAGATACCACAACTTCTGGGTGAGGCACTTGAAGAAGTCGAACAGACTCCGCAAGACCCTGAGGAGGAC-GTCGGAGTCCAACAGAAGGGCGTTCACGACGTGAAGACGGACGTCCACCCGAAACAGGAGGAAAAG---------------------------GAA---GAAGTCTAGGAGCGCTTCGTAG-AAGTC-GTATGGAAGGTTGACC-AGGTGAAATCTTCGGTCCCACCCTCTCCTTAGGAATGATGACGACATACT-AA-AGTAGACCCGGTCTCGGACCAGTA-GCTCGACCAGTGACAACCTACGACGTGAGTTC-CCGGCTACCAAAAGTAGTGATCGTCGTTTTTTACTCC-GAGAGACCTG------GAGT-ATCCGAAACCGTCGTAGATGGTTCAACCTGGTCAG-GTTTAGA-CGTCGTCACAGTCCGTCCCCGTATTCCGGCCCCAGACCGTATCCAA--CGGAAGCCCAGAGGACCTACAGTCACTCTCTC--GACAGACTCCGGCTCTTGAGAGGTACGTTGCTGAAGTAAGGATTCCA-CATCCTTCCGGA-GGCGTCGCCCTCTCTGAACCCCTCTCTTC-CCATGCACGACG-ACTTCTCCTGTAGCAACCCGTCCCGAACCAGGTCTTCACTTCGGGACCTCCACTTTTTCACGAGGGAG-ACCGGTGACATTTACGAGAGCAGTCTGTCTCAAGACACGTTTCGGGCGGAACTCCACCCTCACGACCCCTTTCCTCTGACTCCCCTCCTGTCGCGCACCTT------GAATTTTAAGAAGACCCTAAGTCCTCAACC-CGA-AAGTTCCAGCTACTCGAAGAAATCTAGGAGGAGTTACCACCTGAACCTC-CAACCGAACCCCCTGCGGGCCTCCAGGGACCAACTCCGGGCGTCGGACCACACCCGATAAAGTGACGGTAGTGAGACTACCGCAAAAGAACACA-GGTACTTAGGGTGTCTCA--ATTACCAAGGAAGTCCAAAAAAG--GTAGGTTAGTTCA-TGAACGAGCATCGCGGGACCGCTCCCATCTTGCAGGAT-CCTGGCCTCCTTCTCCTCCTCCTTCTCCTC------TTCCTT-CCACGCTCGACACTCCGCGACCGCCCCACTTCGAGG--AGGAAACCACTTGAAGCCCGTGACCTCT--TGTGTTGGCCCGACCGTCTCACGGCCCTCGAGTACCATGAACACACCCTTCGTCGCAGCAAATATGA-CTTCGACGACAAGAGGGAACAGTCTCCCGAGCTCGGTGGAACCGGGAACGACACCTATAGGTCCGAGTCCCTCCTGAGACGGAAG-GCCAGGTACGGCACCGAGAGACAG-AGATTGCAAGCTA-GACGGTACCTCT-CGTGAAGCGACCTCGATACGGTGTCCCCACGGTGAAAGAGGAACGGACCGGG---ACTCGCTCCTCTGCTCCACC---GGGGGTCTCTCCGACACGACGGCACCGACCCCGGCATCCTCCACCAC-AGCTACGGTGACAGTCGCCTCCACGGGACGTCTATCAATATCGACAACTCAAGCCTAAGGCACGCCAGGGGT---AGCATGGTGAAG-AGCAGTGA-CACCCG---TGACCTACGTAACGACCTCTCC-CATAACGATGACCTCACCGAGCTCATCACTGATAGTCTTACGGGAACAAACAAACCTCTGTTTACTTGGTGACCT--AGGACGGCGT-AATTAGA-GGGTGACCCCGAGAGAGCCTTTCGACCGTC-ATGGAAACTATG-ACGACCCTCT-ACTGTACATCCGAGACATATAAGTACGACGA--GCTCCTGTTTTCTTAACTTTAGGACTCCAGGAGATATTGACATTCNNNNNNNNNNNNNN----NNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNN---NNNN-----NNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNGGAAGCCCGTACCTAAAAACAGACATTCACGGGGGGGGTCTCGGTCCGAGCGACCTACTTCTGTCCAACGATCCCCTCAATGCTGT---GACCCTAGAACGATTTCGGTACCTACTCTGACAGTGTATCTACAGACTTGGTCCTCTG-TTGGCGG---AGAGATCGCCCGCGGGTAGTGACCTCTACGAAGCTCCT--TACAACCGCCGAGAAAGACCTCCTCCGTAGAAGGGCAGCGGAAAGGGTCGTAAATCGGCACTCCAGTACCCTGAG-ACATGGACCCTACCCT-GAACTCCCCT-TACTCT------GGGGAACAACCGTAAGGCGACGGT-AAACAATAAT--AACGCCTTCCCTTTAAACTTGAGTATACTTTGTGGAAGTAAGTAAAACATGAGG-------------TAGTGACCGTACGC-CATCCAAAGCCTCGTCGA------GGAGGCACCCCACGTCAGAAGTACCCCCCCTTACTACTGTTGAAAGTGACATTGCCTACAAGA--GTCCTCTAGTTAGTAGACGAGTACCGAGTCCCACCGGTGACGGTGGAACGTCTAAAGGT-GGTCCG-CCGACGGGACGGAGTCGGGACGG-ACGGTCCGTATCGGTATCCC--GAGGTGTAGACGGTGTTACGGGAGTATCAACTGTACCTTCGGTTCGACAGGATCGGGTAAAGAAGCGTCTCAGTAAAGGTGACTAAGTGTTGGAAAACTTTGTTTGACGTCGGAA---GACTGTTAAAGAAACTATAGGAGT---AACGATTATTTTGAAAGGTGACTTTGTGTGAGAGGAGCGAGTATCTT-CTAAAACTCCGAACACAGTCACCTTCAGGTGGGGTCGTGCAG-AGACGTCCTTGTCAACTTTTGGTGTGAGAAACAAAGGACAAGTTACTCCTGCTACTTAAGTAACCTCTAGGGATCCTCCGTCAGATCAAGGTAT-CGGAACCACATCAGAAA--CCGTGCCTGACGGGTGTGCTACCGGGG-GTACGATGACTCGAGTCGAGGACCTATACCGAAACTGAAAAG--GAAGAAGAACCTCCGGTC-TCTCTACTTGCCTTTGAACGGGCTTCCCAGCTATCCCCACAACCATTGT-AAGAAAAGCCGTTCTAG-AAAGTCCATAAGGACCGC-CCACGCTCAACGGTACCGGGCTTTGAAAAGACTAAATACTCGTCGTAAAAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGACT-TGTGTAAACTG--AACCCCTTTCACTGG--GGAAATCCCTAC-CCTCCAGGTTTCCTTCCGTGTAGAGACCTGGACCAC--TGTCGGTG-TGACATTGTCTGCAACAGCCACGTCCCCGA-CACTCGGGACTGTTACTGTTTCTGTACGACTTTTACCCTAGAC-TATAAAAAACCCGACTTACCGACACGAG---GTCCGAGAACCTTTT-GTTAACAATGTTATAGTAAAGGTTATACGGCGAAAG-CATCCTC---------------------------GACAACGAACAACAACCCACATATCCCGTCGTACCATCTGTGTACCT-TGTATTAGAGT--ATCAGAAACAT-ACAACACATGTCTCTTAC-CCAAGGTCACCTCAGTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-NNNNNNNNNNNNNNNNNNNNNG---------------------------------------------------------------------------------------------------------------------GGAACAAGTAGGTCGAAGTACTCG-------------------ACGAGAC---------------ACTGAAACCCTC--------------------------ACAACT-TC--------------------------------------------------CGG------TTGGCGTCCGT-AACCTGTAAC-------TCGAGCCCCTAGTGCACGAGGTCGTGAAGAAAGTTTCCTTCGGG---AGATCGGCAACTCACGAAACGACAACTTC-CTTATCGTAG--AAGGTTCTGTCTTGGAGA-GTCACA-GTACTCGAGTGATCAAGCTTTCTAGTAAGC-CAT-CAACATACC-TCTAGGTAAAAGGAAGT-AAAGAAGACCAAAAAGGGAAGCATACGAGTGTCGTTGACCTGGTTCTAAGAGT--AGCCG--GGG-TTTTATCAAGGTCACAAGGAAGACCATCTTTTTAA--AGGATATTATCATACGCGGTTCTAGGTGTACAAGGTGCT-ACATAGAAAAGTGAGACA-GATCTACGTTGTGGTTTAAAAGGAACCCGTGAAG-TTCATGACGGTGAGGGCGTAAACACGTTACCCTCCATTCCCATCTGAGTTTCCTCCG--TGA--AAGTGTGGGTG--ATTTCCTAGGTCTCAA--ACTAA---ACCTGTCCGACTAAAAGGAAAGTGGAAGAGGTGGATAGAGTAAGGCTTTTATGCCTGTCGAGTAATGTTCAAGTAACGAAAATAGTAGAGGGACCAGGTACGACTCC------GAAGGAAACCTCAAGTCCAG---GAGACACGACGGGT-AACCCCTCGACTTTGACCGACTTACCCTGTTTCCA-GGTCTCTGGT-TCCTACGGTGCC-GTCGACGGCTTCGAGGACACCACCATAAAGAGG---AAAATTATTGTAACA---CGGACAAGTATTAACGGAGTAAGTTTAGTTTATTATATGAG--ACCTGTA-GTAT-CACCCGTTTTGTGAAACCTGTACAGGTCCCGGACTG-TCTT---AGAAGACTAGCCAACAGAAG-ATCCAGACT-GCT---AAAGGGTT---CAAG----AAGTGGAAA-CCGTAA--CGCGTTAA--ACGCTTTTTACCTACTCAGAGGCC-AAAGTCTGAACCTTGCAGC-GAACTCACCAAAAAGG--GATAGGAATTTTTCGGACTCCAAGAGAGA---------GACACTCTACCCGGGT---AGCTAT----TGTAGGTTTCTCTTTGGCCTTGG---GAACGACGGTGGCCGT-TGCTATTTCAGTG----ACTCGAACCAACCTCAACCTGGAG---CTAGTCTAAACAGTAGAAACATCGGGAAATCGTCTAGTTTTTTCGAAAGAGGTCT------ACAAGGGACGAATAGTTAACTTCTACACGATGGTATAAGGGACACGTCTCTACCGGGTCAGAACATGCACCCACTCTATGTAAGTAGTTTCGATAGGTGTAGGTC-AAGTGAATACCACTATAGTGACAACGAGGCG-ACGGCATAAGACGCCTCTCGAAAGGATCCCCTCGACCCTATCCGAAGTCCGTATTCTTTAGACCTCAGGTACAAGAGGCTACCAGACCAGAACAAAAAGTCC--CACAAGACATACGAGAAGTCG--TACCG--ACTCAG-TAACCTTGACT---CCGATTGTCTTGA--AGGGACTGAGTGTCCTGACGAAAACCTCCTGTACGAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAGCTTGCCGA---TCTAAGAAAATACAAAAGGAAAGAACC---TCCGGTT--AGACGTTAGGA------ACGGGATTGGGGGTAAAACCCTTGAGGCCCCCAGCTCCAATATACCACTCTTGGAAG-TGAAGCTCCTCCGGGATGACCC--CGGTACTGACTACTAGGATTCGA-CGGTAAAACCC---TGGCCTTTGCCGCGTACATCTTCAGTAGTCATACCTGAAACC---CCCACGGTAGACACGGTTGTTGCCTCCAGGACAGTCACCGGTGTCCGGAAAGACAGACACTGGC-AAAGTTCGACCAGTA----");
    // querys.push_back("-ATGACCAGCTTGAAACGGTCACAGACAGAAAGGCCCGTTGCCACTGACAGGGCCTCCATTGTTGGCACAGATGGCACC-C---CCAAAGTCCACACCGATGACTTCTACATG-CGACGCTTCCGGTCCCAAAATGGCAGCTTAGGATCATCAGTA--------------ATGGCTCCAGTAGGGCCCCCTCGAAGTGAAGGTTCTCACCATATCACCTCAACCCCTGGAGTTCCAAAGATGGGGGTTAGGGCA------AGGATCGCAGATTGGCCT-CCAAGAAAGGAAAATATAAAAGAATCT---AGCCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCATGTCCTCCAAAAGCAGTCCTGTGAGTCAGGGAAGTTCTGTTAGCCTCAATTCCAATGACTCAGCCATGCTGAAAAGCATACAGAACACGCTGAAAAGCAAGACCAGACCGTCGGAGAACATGGACTCCAGATTTCTCATGCCTGAAG-CCTATCCCAGCTCCCCTCGGAAAGCTCTCCGCAGAATCCGGCAGCGGAGCAACAGCGACATCACCATAAGTGAACTTGATGTGGATAGCTTTGATGAATGTATCTCACCCACATACAAGAGTGGACCCTCACTGCACAGGGAATATGGCAGCACATCTTCCATT-GATAAGCAGGGAACGTCCGGAGACAGCTTTTTTGATTTGTTAAAGGGC-TACAAAGATGACAAATCTGATC---GAGGTCCAACTCCAACCAAGCTCAGTGACTTTCTCAT-TGCTGGTGGGGGCAAGGGTTCCGGTTTCTCTTTGGATGTTATAGA---CGGGCCCATCTC--ACAGAGAGAGAACCTCAGGCTTTTTAAGGAAAGGGAAAAACCACTCAAGCGGCGCTCGAAGTCTGAAACTGGAGACTCATCCATCTTTCGTAAATTACGCAATGCCAAAGGTGAA---GAAC---TCGGGAAA---TCGTCTGATCTAGAAGATAACC--GATCAGAAGA---CTCTGTCAGGCCCTGGACATGTCCAAAGTGCTTTGCCCACTACGATGTCCAGAGTATACTGTTTGATTTGAATGAGGC-GATTATGAACAGGC---ACAATGTTATTAAGAGGAGAAACACCACCACGGGAGCCTCGGCAGCTG--CCGTGGCATCCCTGGTCTCTGGACCTCTGTCC-CATTCAGCCAGCTTCAGCTCCCCCATGGGCAGCACAGAGGACCTGAATTCCAAAGGGAG------TCTCAGCATGGACCAGGGAGA-TGACAAAAGCAACGAACTGGTAATGAGCTGTCCATACTTTCGGAATGAGATAGGTGGAGAAGGTGAAAGAAAAATCAGCCTGTCAAAGTCTAACTCTGGCTCCTTTAGTGGATGTGAAAGTGCCTCCTTCGAATCTACTCTTAGCTCCCATTGCACAAACGCAGGAGTGGCAGTACTTGAAGTGCCCAAGGAAAACTTGGTGTTGCATTTAGATAGAGTGAAAAGATACATCGTGGAACACGTGGATTTGGGCGCATACTATTATAGAAAGTTTTTCTACCAGAAGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-CTTATGACACTG-AGAGGTTCTGTCCTGGAAGATGCCATTCCCTCCACAGCGAAGCACTCGACAGCCAGAGGGCTGCCTCTGAAGGAAGTGCTAGAGCACGTGATCCCTGAGCT-------CAATGTCCAGTGCCTGCGGTTGGC--------------------------------------------------CTTCAACA---------------------------CTCCCAAAG-TCA--------------------CAGAGCA---------------------AC-TCATGAAACTGG-ATGAACAAGG---------------------------------------------------------------------------------------------------------------------GCTGAACTATCAGCAGA-AAGTGGGCATCATGTACTGCAAAGCTGGACAAAGTACCGAAGAAGAGATGTACAACAATGAATCAGCCGGGCCAGCCTTTGAGGAATTCCTTCAGCTTT-TGGGAGAACGAGTTCGCCTCAAAGGA-TTTGAGAAGTATCGTGCACAGCTTGACACCAAAACTGA-CTCCACTGGAACCCATTCTCTGTATACAACA-TACAAAGACTATGAAATTATGTTTCAT-GTCTCTACCATGCTGCCATACACGCCTAACAACAAGCAACAG---------------------------CTCCTACGGAAA-CGGCACATTGGGAATGATATTGTGACAATTGTTTTCCAAGAGCCTGGCGCACAGCCA-TTCAGCCCAAAAAATATTCGATCCCACTTTCAGCACGTTTTTGTCGTTGTCAGGGCCCATAACCCGTGCACTG-ACAGTGT-CTGTT-ACAGTGT-GGCTGTCACCAGGTCCAGAGAT-GTGCCTTCCTTTGGGCCTCC--CATTCCTA-AAGGGGTCACT--TTCCCTAAGTCAAATGT-GTTCAGGGACTTCCTTTTGGCCAAAGTGATTAATGCAGAAAATGCTGCTC-ATAAATCAGAAAAGTTCCGGGCTA--TGGCAACTCGG-ACCCGCCAGGAATACCTGAAAGATCTAGCAGAAAAGAATGTCACCAACACACCTATCGACCCT-TCTGGCAAGTTTCCGTTCATCTCTCTGGCCTCTAAGAAGAAGGAAAAGTCTAAGCCCTATCCAGGAGCTGAGCTCAGTAGCATGGGGGCCATCGTGTGGGCAGTCCGGGCCAAAGAC-TACAACAAGGCTATGGAAATAGACTGCCTCTTAGGGA--TCTCCAATGAGTTCATCGTCCTCATCG-AGCAGGAAACAAAGAGCGTGGTTTTCAATTGTTCCTGCAGAGATGTGATAGGGTGGACTTCAACTGACACCAGTCTCAAGATTTTCTATGAACGAGGAGAATGTGTTTCAGTGGAGAGTTTCATTAGCAA---TGAGGATATCAAAGA-GATTGTCAAAAGGTTGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAGGGCAGCC--GCCTG---GAGATCTGCAAGGTGGCAGTGGCCACT---A-GTCACGAGCAGATGATCGATCTCCTGAGA---TCCGTC--ACC-GTGAAGGTCGTCATCATTCCCCCCCATG---ATGACTGCACCCCACGGAG------G------AGTTGCTCTGAAGCCTACCGCATGCCAG-TGAT-------------GGAGTACAAGATGAACGAAGGCGTGTCGTACGAATTCAAGTTTC-CCTTCCGAAAT-AACAACAAAT---GGCAGAGGAATGCCAACAAGGGG------C-CTCATTCACCTCAAG-TTCCGTCCCAGGTGCAGA-GTCCCATGACCTCACGGCTGAATGCTGGGAAAGGAGAAGGGAAGA-TGCCTCCTCCTGAAAGAGCCGCCAATATTCCTCGAAGCATTTCCAGTGACGGACGTCCACTAGAAA---GGCGG-CTGTCTCCTGGTTCGGACATCTATGTGACGGTCTCGTCCATGGCTTTAGCAAGATCCCAG---TGTCGTAACTCCCCTAGCAACCTGTCTTCATCCAGTGAGACTGGCTCTGGGGGGGGCACTTACAGACAAAAGTCCATGCCCGAAGGGTTTGGAGTGAGCCGCAGATCCCCAGCCTCCATCGACAGGCAGAACACTCAGTCAGATAT---CAGTGGAAGCGGAAAGTCGACACCCAGCTGGCAAAGAAGTGAGGATAGCATTGCCG-ACCAGATGGCTTACAGTT---ATAGAGGACCTCAGGATTT--CA--ATTCTTTTGTCCTCGAGCAGCATGAA--TATA-----------CAGAGCCAA-CATGCCA--TCTCCCAGCA-GTATCA--AAGG----TACTGCC----CACGTTCCGAGAGAGCC-CCAGTGG-----GAGATTA--ATGCGGCAGGATCCA--GTGGTTCATTTGTCTCCAAACAAACAAGGG-CAT-TCTGACAGC--CACTACTCGAGCCACT-CCAGTAG--CAATACCCTCTCCAGCAACGCATCGAGTGCCCACAGCGATGAGA--AGTGGT--AC-GA---TGGGGACCGCACAGAA-TCCGAACTCAACA--GCTACAACTATCTGCAGGGCACC-TCTGCG-GACAGCGGCA-TCGACACCACCTCCTATGGCCCCAGCCACGGCAGCACAGCCTCTCTGGGGGCTGCCACGTCATCGCCTCGCTCA---GGGCCAGGCAAGGAGAAAGTGGCGCCCCTGTG--GCAT--AGCTCTAGTG--AGGTGATCTCCATGG--CA--GA----TCGGACTTTAGAGACAGAGAG-CCACAGCATGGACCGGAAGGCAGAGTCCTCCCTG-AGTTTGGATATCCAC--AGCAAGAGCCAAGCCAGCTCGAACCCTCTGTCAAGAGAGAACAGCACTTTCAGTATAA---ATGACGCTGCTTCC-CACACAAGTACCATGAGCTCCCGACACTCTG-CCAG-CCCAGTCGTTTTCACCAGTGCCAGAAGTTCACCAAAAGAAGAGCTTCACCCTGCCACCCCCTCA---CAGCTCGCCCCTTCCTT-------CTCCTCGTCCTCCTCCTCCTCGTCTGG----TCCTAGGACT-TTCTACCCTCGCCAGGGCGCTACTAGCAAGTACCTG-ATTGGATGGAAAA--AACCTGAAGGAACC-ATAAACTC--CGTGGGATTTATGGATACAAGA-AAG-CGTCA-TCAGAG--TGATGGCA--ATGAA-ATAGCGCACACCAGGCTG-CGTGCCTCAACCAGGGACCTCCGTGCATCCCCCAAGCC--AACTTCCAAGTCCACCATCGAAGAGGA---------TCTAAAGAAACTGATCGACCTTGAGA----GCCCAACTC-CCGAATCGCAGAA--GAACTTTAA--G------TTCCACGCGCTG------TCCTCCCCGCAGTCTCCTTTCCCCCCGACCCCCACCTCCAGGCGGGCCCTGCACAGGAC--GCTGTCGGACGAGAG-CAT------TTACAGTGGCCAGAGGGAGCACTTCTTCACCTCCAGA-GCTTCACTTCTGGACCAAGCCCTGCCCAATGACGTC-CTCT-TCAGCAGCACATACCCCTCTCTCCCCAAGTCTCTCCCACTGCGGAGACCTTCCTACACCTTGGGGATGAAGTCACTGCATGGGGAGTTCTCAGCCTCAGAC-AGCTCCCTCACCGACATCCAGGAGACCCGGAGGCAGCCCATGCCCGACCCTGGCCTGATGCCTCTGCCTGACACCGCTGCGG--ATTT--GGATTG-GTCCAACTTGG---TAGATGCTGCCAA--AGCCTATGAG------GNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAGTAAGGATTCCTCTCCTACTCTGGCTTCTAAAGTGGACCAGTTGGAAGGTATGCTGAAGATGCTTCGGGAAGATTTGAAG---AAG-------------GAAAAGGAAGGCAAAGCCCCCCTTCAGGCGGAAGTCCAGCACTTGCGAGAGGACAACCTGAGGCTGCAGGAGGAGTCTCAGAATGCCTCAGACAAGCTGAAGAAGTTCACGGAATGGGTCTTCAACACCATAGACATGAGCTA---");

    // alignGrpToGrp(refs,querys, param, alignments);

    // std::string ref = "----------------------------------------------------------------------------------------------------------AGTTGGAGACGGACGTCCACCCGGAACAGAAGGAAAAGCCTCTTGTAAATGGAA---GAAGTTTAGAAGAGCTTCGTAAAAGTCGTATGGAAGGTTGACCAGGTGAAATCTTCGGTCTCACCCTCTCCTCAGGAATGACGA---CATGCTGAAGTAGACCTGGTCCCGGACGAGGAGATCGACCAGCGACAACCTCCGGCGTGAGTTCCCCGCTACCAAAAGTAGTGATCGTCGTTTTTTACTCCGAGAGACCTG------GAGTATCCGAAACCGCCGTAGATGATCCAACCTGGTCAGGTTCAGGCTTCGTCACAGTCCGTCCCCGTAGTCCGGTCCCAGTCCTTATCCGACGGAGGCCCAGAGGACCTACAGCCATTCCCTCGACAGCCTTCGTCTCTTGAGAGGTACGTCGCTGAAGTAAGGGTTACACATTCTACCGGAGGAGTCACCTTCGCTGAACCCCTCTCTACCCATCCACGACGACTTCTCCTGTAGTAACCCGTCCCGAACCAGATCTTCACTCCGGGACCTCCACTTCTTTACGAGGGAGACCGACGACATTTACGAGAGTAGTCTGTCTCAGGACACGTCCCGGGCGGCCCTCCATCCCCACGACCCCTTGCCCCTGACACCCCTCCTGTCAGGTACCTT------GAACTTCAAAAAGACACTAAGTCCTCAACCCGAGAGTTCCAGCTACTCAAAGAAATCTAGGAGGAGTTACCACCTGAACCTTCAACCGAAACCCCTACGGGCTTCCAGGGACCAACTCCGAGCTTCGGATCACACCCGATAGAGTAATGGTAGTGAGACTACAGCGAAAGAACATAGGTATTTAGGGTGCCTCAATTACCAAGGAAGACCGAAAAAGGTAGGTTAGTCCATAAACGATCACCGCGGGACCGCCCCTATCTTTCAGGATCCAGGTCTCCTTCTACTCCTTCTCCTTCTC------TTCCTTCCACGCTCGACCCTCCTGCACCACCCCACTTCGAGGAGAAACCCACTTGAAGACCGTGACCTCTTGTGGTGACCCGACCGTCTCACAGCCCTCGAGTACCATGACCACACTCTACACCGTAGTAAATACGACTTCCACGACAAGAGGGAACAGTCGCCCGAACTTGGCGGAACCGAGAACGACACCTACAGGTCCGAGTCCCTCCTGAGACAAAAGGACAGGTACGGTACCGGGAGTCAGAGATTCCAGGCTAGACGGTCCCTCTCGTGGAGTGACCTCGACACGGTGTCCCCACGGTGGAAGAGAAACGGACCGGG---ACTTGCTCCACTCCTACACCTCCGAGGGTCCCTCCGGCACGACGGTACCGACTCCGGCATACTCCGCCACAGCTACGGCGACAGTCGTCTACACGGGACATCCATCAATATCGACAACTCCAGCCTGAGCCACGCCAGAGGT---AGCATGGTGAAGAGCAGTGACACCCGTGAGCTCCGCAACGACCTGTCTCACAACGACGACCTCACCGACCTCATCACCGACAGTCTTACCGGAACGAACAAACCTCTGTTCACTTGGTGCCCTAGGACGGCGTAGTTGGAGGGTGACCCCGAGAGAGCCTTTCGTCCGTCATGGAATCTATGACGACCCTCTACCGTACATCCGAG---------------------------------------------------------------GTAGACCAGCCGTTCCGACAGGAGTGAAGAGACGGTCGACCCACACCTAAAAGGTGACGGTGGT---TATAGGCTGACCCACAAAACGGACAGCTACCTTCGACCCCTAGAAGCCGAGTGAGGTTTGGGAAGTCCGTACCTAAAAACAGACATTCACGGTGGGGGTCTCGGTCAGAGCGACCTGCTACTGTTCAACGACCCCCTCAATGCCGT---GACCCTAGAGCGGTTTCGGTACCTACTCTGACAGTGTATCTACAGGCTTGGGCCTCTGTCGGAGG---AAAGGTCACCCGCGGGTAGTGACCTCTACGAAGCTCCCTACAACCGCCGAGAAAGACCTCCCCCGTAGAAAGGTAGAGGGAAGGGTCGTAAGTCGGCACTTCAGTACCCCGAGACGTCGACGCTACCTT-GGACCCCGCTCACCCG------TGGGAACGGCCGTAATGCGACGGTAAACAACAATGAGGCCTTCCCTTTGAACATGAGTATCCTTTGCGGAAGTAAGTAGAACATGAGG-------------TAGTGACCGTACGCCATCCAAAGACTCGTTGA------GGAGGCTCCCCACGTCAGTAGTACCCCTCCCTATTACTGTTGAAAGTGACACTGACTGCAAGAGTTCTCTAGCTAGTAGACAAGTACCGAGTCCCACCGGTGACGGTGGAACGTCTAGAGGTGGTCCGCCGACGGGACGGAGTTTGGACGAACGGTACGTATCGGTATCCCAAGATGTAGGCGGTGCTACGGGAGTATCAACTGTACCTTCGGTTCGACCGGATCGGGTAAAGAAGCGTCGCAGTAGAGGTGACTAAGTGTTGGAAAACTTTGTTTGACGTCGGAAAACTGTTAAAGAAACTATAGAAGT---GACGAGTACTTCGAGAGGTGTCTTTACGTAAGAGGTGCGAGTATCTTCTAGAATTCCGACGACAGCGACCTTCAGGTGGGATAGTGTAGAGACGTCCTCGTTAACTTTTGGTGCGAGAAACAAAGGACGAGTTACTCCTGCTACTTGAGTAACCTCTAAGGTTCCTCCGTCAGCTTAAGGTACCGGAACAACATCAGAAACCGGGCCTGACGGGTGTGTTACCGGGGGTACGATGACTCGAGTCGAGGACCCATTCCGAATCTGAAAAGGAAGAAGAATCTTCGGTCTCTTTATTTACCTTTGAACGGTCTTCCCAGTTATCCACACAACCACTGTAAGAAAAGACGGTCTAGAAAGTCCATAAGGACCGCCCAGGAACAACGGTACCGGGACTTGAAGAGACTAAATACTCGTCGCAAAAGACGTAATTAGTGAAACCGGTTTTCCTTCAGGGACTTGTGCAAACTGAACCCCTTCCACTGGGGAAATCCCTACCCTCCAGGTTTTCTTCCGTGTAGAGACCTGGACCACTGACGGTGTGACATTGTCTGTGAGAGTCACGTCCCCAACACTCGGGCCTGTTACTGTTTTTGTACGACTTTCACTCTAGCCTACAAAAAGCCCGATTTACCAACACGAGGTCCGAGGACCTTTTGATAACAGTGTTACAGTAAGGGTTACACGGCGAAGGAGTCCTC---------------------------GACAACGAACAACAATCCACACATCCCGTCGTACCACCTCTGCACCTTGTAATAGAGTATCAGAAACATCCAACACATGTCTCTCACCCAAGGTCACCTCAGTCAAAACCACAGCTCGACGCGTGCTATGAAAAGTTTAGGAAAATCCGCTTGAGCGAGAGGGTCATCGACTTCTTTAAGGAGTTTCCGGCCCGGTCGCCTGAGTAATAACATGTAGAGGAGGAGCCACGAGACCGGACGAAACGTCATGTACTACGGATGAAAGACGACTATCAAGTCG---------------------------------------------------------------------------------------------------------------------GGAACGAGCAGGTCAAAGTATTCG-------------------ACGAGGC---------------ACTGAAACCCGC--------------------------ACAACTTC--------------------------------------------------CGGTTTGCGTCCGTGACGTGCAAC-------TCGAGACCCTAGTGTACAAGGTCGTGGAGAAAGTCGCCGTTAGGAGACCGACAGCTCACGAACCGACACCTCCCTTACCGTAGGAGGTCCTGTCTGGGAGA-GTCACAGTACTCGAGTGATCAAGATTTTTAATAAGCCATCAACATCCCTCTAGGTAAAAGGAAGTACAGAAGACCAAAAAGGGAAGCTTACGAGTGTCGGTGACCTGGCTCCAAAAGTAGTCGGGGTTTTATCAAGGTCACAAGGAAAACCATCTTCGTGAAGGACATTATCATACGAGGTTCTAGGTGCACGAGGTGTCACATGGCAGAGTGAGACAGATCTACGTCGTAGTTCGAGAGGAACCCGTGAAGCTCTTGGCGGTGAGGACGCAAACACGTCACTCTCGACTCCCGTCTGAGTTTCCTGCACGAGAGTGTGGGTGACTTACTCGGTCTTAACCTGAACCTGTCCGACTAAAAGGAGAGTGGAAGAGGGGGTTAGAGTAAGGCTTTCATGCCTGTCGAGTACTGCTCAAGTAACGAGAACAGTAGAGGGACCAGGTACGGCTCC------GAAGGAAACCTTAACTCCAGGAGACACGACGGGTATCCTCTCGACTTCGACCGACTCACTCTGTTTCCAGGGCTCTGATTCCTGCGGTGTCGCCGACGCCTTCGAGGGCAACACCACAAAGAGGAGAATTAGTGTAATA---CGGACAAGTATTACCGAAGTAAGTTCAGTTTGTTATACGAGACCTGTAGTATCACCCGTTTCGTGAAACCTGTACAGGTCCCGGAGTGTCTT---AGGAGACTAGACAACAGAAGTTCCAGACTACT---GAAGGGCT---CAAG---AAGTGGAAACCGTAACGCATTAAATGCTTTTTACCTACTCAGAGGTCAAAGTCTGAATCTCGCAGCGAACTCACCAAAAAGGGAAAGGAACTTATCGGACTCCAAGAGAGAGACACTCTACCCAGGT---AGCTAGTGTAGGTTCCTCTTTGGTCTTGGGAACGGGGGTGGTCACTACTCCTTCAGTGACTCAAACCAACCTCAACCTGGAG---CCAGTCGAGACAGCAGAAATATCGGGAAATTGTTTGGCTTCTTCGAAAGAGGCCTACAGGGGACGAACAGCTAACTTCTACACGATGGTATAAGGGACACGTTACTACCGGGGCTGAACATCCAACCCCTCTATGTAAGTAGTTTTGATAGGTGTAGTTCGAGTGAATACCACTATAGTGACAACGATGCGACTGCCTAAGCCGCTTCCCGAAAGGACCCCCTCGACCCCATCCGAAGTCCGTACTCCTTAGAGCTCAGGTACGAGAGGCGACCGGGGCAGAACAAGAAGTCCCACAAGACCTACGAGAAGTCGTACCGACTCAGTAACCTTAACTCCGATTGTCTTGAGGGGACTGAGTGTCCTGACGGAAACCTCCTGTCCGAGAGTTCCGTTGAACTCCAAAGATAAAGGACCGAACTTGCCGA---TCTAAGAAAATGTAAAAGGAAAGAACCTCCGGTTAGACGTTAAGA------ACGGGATTGGGGGTAGAACCCCTGAGGCCCCCAACTCCAATATACCACTCCTGGAAGTGAAGCTCCCCCGGGATGACGTCGGTACTGACTACTAGGATTCGACGGTAAAACCCTGGCCTTCGCTGCGTACATCTTTAGTAGTCACACTTGAAACC---CCCACGGTAGACACGACTGTTGTCTCCGAGACAGTCGTCACTGTCCGGAAAGCCAGACGCTGGAGAAGTTCGACCAGTA";
    // std::string query = "ATGACAAGCTTGAAACGATCACAGACAGAAAGGCCTGTTGCCACTGATAGGACCTCTGTTGTTGGCACGGATGGCCCTC---CCAAGGTCCACACGGATGATTTCTACATGCGCCGCTTCCGGTCCCAAAATGGCAGCTTAGGATCATCAGTCATGGCTCCCGTAGGGCCCCCACGAAGTGAAGGTTCTCATCATATAACCTCCACTCCTGGAGTCCCAAAGATGGGGGTTCGAGCA------AGGATTGCAGATTGGCCTCCAAGAAAGGAAAACATAAAAGAGTCT---AACCGTTCAAGCCAGGAAATAGAAACCTCAAGTTGCCTTGAGAGCATGTCCTCCAAAAGCAGTCCTGTGAGT---GGAAGTTCTGTTAGCCTCAATTCTAATGACTCAGCCATGCTAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCTGGAGAA---TTTTTTGATTTGTTAAAGGGCTACAAAGATGACAAATCTGATC---GTGGGCCAACTCCAACCAAGCTCAGTGACTTTCTTATTGCTGGTGGGGGCAAGGGTTCTGGTTTCTCTTTGGATGTTATTGA---TGGGCCCATTTCACAGAGAGAGAATCTCAGGCTTTTTAAGGAAAGGGAAAAACCACTCAAGCGACGGTCCAAGTCTGAAACTGGAGACTCATCCATTTTTCGTAAATTACGCAATGCCAAAGGTGAA---GAAC---TTGGGAAA---TCGTCGGACCTTGAAGACAACCGATCAGAAGA---CTCTGTCAGGCCTTGGACTTGTCCAAAGTGCTTTGCACATTATGATGTCCAGAGCATATTATTTGATTTGAATGAGGCGATTATGAACCGGC---ACAACGTTATTAAGAGGAGAAACACCACCACAGGAGCCTCGGCAGCTGCTGTGGCGTCCTTGGTCTCTGGACCTTTGTCCCATTCGGCCAGTTTCAGCTCACCAATGGGCAGCACAGAGGACCTGAATTCCAAGGGAAG------CCTTGGCATGGACCAGGGAGATGATAAAAGCAATGAACTTGTAATGAGTTGTCCGTATTTTCGGAATGAGATAGGTGGAGAAGGTGAAAGGAAAATCAGCCTGTCAAAATCAAATTCTGGCTCGTTTAGTGGCTGTGAAAGTGCCTCCTTTGAGTCTACCCTCAGCTCCCATTGCACAAATGCGGGAGTGGCAGTACTTGAAGTCCCCAAGGAAAACTTGGTGTTGCATCTAGATAGAGTGAAAAGGTACATTGTGGAGCATGTGGATCTCGGCGCTTACTATTATAGGAAATTCTTCTACCAGAAGGAACACTGGAACTATTTTGGGGCTGATGAGAATCTTGGTCCGGTGGCTGTAAGCATTCGAAGGGAAAAACCAGAAGAAATGAAAGAAAATGGATCTCCATACAACTATCGAATAATTTTTAGAACTAGTGAGCTCATGACACTA-AGGGGTTCTGTCCTGGAAGATGCCATTCCCTCCACAGCGAAGCATTCAACAGCCAGAGGACTGCCTCTGAAAGAGGTTCTGGAACATGTGATCCCAGAGCT-------CAATGTCCAGTGCCTGCGGCTGGC--------------------------------------------------CTTCAACA--------------------------CACCCAAAGTCA---------------CGGAGCA-------------------GCTCATGAAACTGGATGAGCAAGG---------------------------------------------------------------------------------------------------------------------GCTGAACTATCAGCAGAAAGTAGGCATCATGTACTGCAAAGCTGGCCAGAGTACTGAAGAAGAGATGTACAACAATGAATCAGCAGGCCCAGCCTTTGAAGAATTCCTTCAACTGTTGGGAGAAAGAGTTCGGCTAAAAGGATTTGAGAAGTATCGTGCACAACTTGATACCAAAACTGACTCCACTGGAACCCATTCTCTGTACACGGCATACAAAGACTATGAAATTATGTTCCATGTTTCTACTATGCTGCCATATACACCTAACAACAAGCAACAG---------------------------CTCCTACGGAAGCGTCACATTGGCAATGACATTGTAACAATTGTTTTCCAAGAGCCCGGAGCACAGCCATTTAGCCCCAAAAACATCCGGTCTCACTTTCAGCATGTTTTCGTTATTGTCCGAGCTCACAACCCTTGCACGGACAGTGTCTGTTACAGTGTGGCAGTCACCAGGTCCAGAGATGTGCCTTCCTTTGGACCTCCTATCCCAAAAGGGGTCACATTTCCCAAGTCAAATGTGTTCAGGGACTTCCTTTTGGCCAAAGTGATTAATGCAGAGAATGCTGCTCATAAATCAGAAAAGTTCCGGGCTATGGCAACTCGGACTCGCCAGGAATACCTGAAAGATTTGGCAGAAAAGAATGTCACCAACACCCCTATTGACCCTTCTGGCAAATTTCCGTTCATCTCTCTGGCCTCTAAGAAGAAGGAAAAGTCTAAGCCATATCCAGGAGCTGAGCTCAGTAGCATGGGGGCCATTGTGTGGGCAGTCCGGGCCAAGGACTATAACAAGGCTATGGAAATAGACTGTCTTCTCGGGATCTCCAATGAGTTCATCGTCCTTATTGAACAGGAAACAAAGAGTGTGGTTTTCAATTGTTCCTGCAGAGATGTGATAGGGTGGACCTCAACTGACACCAGTGTCAAAATCTTCTATGAGCGAGGAGAGTGCGTGTCGGTAGAGAGTTTCATTAGCAA---TGAAGATATCAAAGAGATTGTCAAACGGCTGCAGTTTGTTTCAAAAGGTTGTGAATCAGTGGAAATGACTCTTCGAAGAAATGGTTTAGGACAGCTTGGCTTTCATGTGAACTACGAGGGCATTGTGGCAGATGTAGAACCGTATGGCTATGCATGGCAGGCAGGGCTGAGGCAGGGCAGCCGCCTGGTGGAGATCTGCAAGGTGGCTGTGGCCACCCTGAGCCATGAGCAGATGATTGATCTCCTGAGGACATCCGTCACAGTGAAAGTTGTCATCATTCCCCCCCATGATGACTGCACCCCACGAAGG------AGTTGCTCGGAAACCTACCGCATGCCAGTGAT-------------GGAGTATAAAATGAATGATGGCGTTTCATATGAGTTCAAGTTTCCCTTCCGAAATAATAACAAATGGCAGCGGAATGCTAACAAAGGT------CCTCATTCACCTCAAG-TCCCATCCCAGGTACAGAGTCCCATGACCTCCCGGATGAATACCGGGAAAGGTGATGGGAAAATGCCTCCTCCAGAAAGAGCAGCCAACATTCCTCGAAGCATCTCCAGCGACGGACGCCCACTAGAGA---GGCGCCTGTCTCCGGGTTCGGACATCTATGTGACTGTCTCATCCATGGCTTTAGCAAGATCTCAG---TGTCGTAACTCCCCTAGCAACCTGTCTTCATCCAGCGAGACAGGCTCTGGGGGCAGCACCTACAGACAGAAATCCATGCCTGAAGGGTTTGGAGTGAGTCGCAGATCACCAGCTTCCATGGACAGGCAGAACACCCAGTCAGACCT---TGGGGGCAGTGGAAAATCCACACCCAGCTGGCAAAGAAGTGAGGATAGCATTGCCGACCCAATAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAACCTACATGCCATCTCCCAGCAGTATCAAAGGTACTGCCAGCTTTCCGAGAGAGCCCCAGTGGGAGATTAATACGGCAGGATCCAGTGGTCCATTTGTCTCCAAACAAACAAGGGCATGCTGATAGCCACTACTCGAGCCACTCCAGTAGCAATACCCTCTCCAGCAATGCATCGAGTGCCCACAGTGATGAGAAGTGGTACGA---TGGGGACCGCACAGAATCTGAACTCAATAGCTACAACTATCTTCAGGGCACCTCTGCTGACAGTGGCATCGACACCACCTCCTATGGCCCCAGCCATGGCAGCACAGCTTCCCTAGGGGCTGCCACCTCATCACCCCGTTCA---GGGCCAGGCAAGGAGAAAGTGCCCCCCCTGTGGCATAGCTCCAGTGAAGTGCTCTCCATGACAGAACGGACTTTAGAGACAGAGAGCCACGGCATGGACCGTAAAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTACCATGAGCTCCCGGCACTCTGCCAGCCCAGTTGTTTTCTCCAGTGCCAGAAGTTCACCTAAAGAAGAACTTCATCCTGCCGCCTCCTCACAGCTTGCACCTTCCTT------CTCTTCCTCCTCCTCCTCTTCCTCTGGTCCTAGAACATTCTACCCCCGTCAGGGAGCTACCAGCAAGTACCTGATTGGATGGAAAAAACCCGAAGGAACTATTAACTCCGTGGGATTTATGGACACAAGAAAGCGTCATCAGAGTGATGGCAATGAGATATCCCACACTAGGCTGCGTGCCTCAACCCGGGACCTCCGGGCATCCCCCAAACCAACCTCCAAGTCCACCATTGAAGAGGATCTAAAAAAACTCATCGACCTTGAAAGCCCAACTTCTGAATCACAGAAGAATTTTAAG------TTCCATGCACTGTCCTCCCCTCAGTCTCCTTTTCCCCCAACGCCTATCTCGCGGCGAGCTTTGCACAGGACACTGTCAGACGAGAGCATTTACAGCAGCCAGAGGGAGCACTTCTTCATCTCAAGGGCTTCATTTCTGGACCAAGCGCTGCCCAATGATGTCCTTTTCAGCAGCACGTACCCATCTCTCCCCAAG---CTGCCACTGCGGAGGCCGTCATACACGTTAGGAATGAAGTCATTGCATGGGGAGTTCTCGGCCTCAGACAGCTCCCTCACCGACATCCAGGAGACACGAAGGCAGCCCATGCCTGACCCTGGCTTGATGCCCCTGCCTGACACTGCCGCAGATTTGGACTGGTCCAACTTGGTAGATGCTGCCAAAGCCTATGAG------GTCCAGAGAGCCTCATTTTTTGCTGCTAGTGATGAAAACCATCGCCCCCTGAGTGCGGCGTCCAACAGTGACCAGCTTGATGAACAGGCTCTGGTGCAGATGAAGGCATACAGCAGCAGTAAGGATTCCTCTCTCACTCTGGCTTCTAAAGTGGACCAGTTGGAAGGTATGCTGAAGATGCTTCGAGAAGATCTGAAG---AAG-------------GAAAAAGAAGACAAAGCCCACCTTCAGGCTGAAGTGCAGCACTTGCGGGAGGACAACCTGAGGTTGCAGGAGGAGTCACAGAACGCCTCTGACAAGCTGAAGAAGTTCACAGAGTGGGTCTTCAATACCATAGATATGAGCTAG"; 
    // std::pair<std::string, std::string> alignment;
    // alignSeqToSeq(ref, query, param, alignment);
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " <<  mainTime.count() / 1000000000 << " s\n";
    return 0;
}


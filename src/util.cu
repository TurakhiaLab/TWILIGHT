#ifndef UTIL_HPP
#include "util.cuh"
#endif

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
        // std::vector<std::vector<uint16_t>> freq;
        std::vector<uint16_t*> freq;
        std::vector<std::pair<int32_t, int32_t>> seqIdx;
        std::vector<std::pair<int32_t, int32_t>> len;
        // store info to array 
        auto freqStart = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < alignSize; ++n) {
            int32_t nIdx = n + r*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t qryLen = tree->allNodes[nodes[nIdx].second->identifier]->msa[0].size();
            int32_t refLen = tree->allNodes[nodes[nIdx].first->identifier]->msa[0].size();
            refIdx = seqNum;
            
            // std::vector<uint16_t> temp;
            // for (int i = 0; i < 12*seqLen; ++i) temp.push_back(0);
            uint16_t *temp = new uint16_t[12*seqLen]; 
            for (int i = 0; i < 12*seqLen; ++i) temp[i]=0;
            // assert(temp.size() == 12*seqLen);
            // tbb::blocked_range<int> rangeRef(0, refLen);
            for (auto seq: tree->allNodes[nodes[nIdx].first->identifier]->msa) {
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (seq[s] == 'A' || seq[s] == 'a') temp[6*s+0]+=1;
                    else if (seq[s] == 'C' || seq[s] == 'c') temp[6*s+1]+=1;
                    else if (seq[s] == 'G' || seq[s] == 'g') temp[6*s+2]+=1;
                    else if (seq[s] == 'T' || seq[s] == 't') temp[6*s+3]+=1;
                    else if (seq[s] == 'N' || seq[s] == 'n') temp[6*s+4]+=1;
                    else                                     temp[6*s+5]+=1;
                }
                });
                ++seqNum;
                seqs.push_back(seq);
            }
            qryIdx = seqNum;
            for (auto seq: tree->allNodes[nodes[nIdx].second->identifier]->msa) {
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (seq[s] == 'A' || seq[s] == 'a') temp[6*(seqLen+s)+0]+=1;
                    else if (seq[s] == 'C' || seq[s] == 'c') temp[6*(seqLen+s)+1]+=1;
                    else if (seq[s] == 'G' || seq[s] == 'g') temp[6*(seqLen+s)+2]+=1;
                    else if (seq[s] == 'T' || seq[s] == 't') temp[6*(seqLen+s)+3]+=1;
                    else if (seq[s] == 'N' || seq[s] == 'n') temp[6*(seqLen+s)+4]+=1;
                    else                                     temp[6*(seqLen+s)+5]+=1;
                }
                });
                ++seqNum;
                seqs.push_back(seq);
            }
            // printf("len: (%d, %d), num: (%d, %d)\n", refLen, qryLen, refIdx, qryIdx);
            seqIdx.push_back(std::make_pair(refIdx, qryIdx));
            len.push_back(std::make_pair(refLen, qryLen));
            freq.push_back(temp);
        }
        auto freqEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds freqTime = freqEnd -freqStart;
        printf("Calculate frequency time : %d ms\n",  (freqTime.count() / 1000000));
        // Malloc
        uint16_t* hostFreq = (uint16_t*)malloc(12*seqLen * pairNum * sizeof(uint16_t));
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
        uint16_t* deviceFreq;
        int8_t* deviceAln;
        int32_t* deviceLen;
        int32_t* deviceAlnLen;
        int32_t* deviceSeqInfo;
        int16_t* deviceParam;
        auto kernelStart = std::chrono::high_resolution_clock::now();
        cudaMalloc((void**)&deviceFreq, 12*seqLen * pairNum * sizeof(uint16_t));
        cudaMalloc((void**)&deviceAln, 2*seqLen * pairNum * sizeof(int8_t));
        cudaMalloc((void**)&deviceLen, 2*pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceAlnLen, pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
        cudaMalloc((void**)&deviceParam, 6 * sizeof(int16_t));
        // Copy to device
        cudaMemcpy(deviceFreq, hostFreq, 12*seqLen * pairNum * sizeof(uint16_t), cudaMemcpyHostToDevice);
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
        auto reAlnStart = std::chrono::high_resolution_clock::now();
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
        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - reAlnStart;
        printf("Round. %d align %d pairs. reAlnTime: %d ms\n", r, alignSize, reAlnTime.count() / 1000000);
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


/*
void msaPostOrderTraversal_gpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{
    // assign msa to all nodes
    for (auto n: nodes) {
        if (n.first->children.size()==0) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else {
            if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
                Node* node = tree->allNodes[n.first->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                        if (node->children[childIndex]->msaIdx.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                        tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                        util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                        break;
                    }
                }
            }
        }
        if (n.second->children.size()==0) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else {
            if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
                Node* node = tree->allNodes[n.second->identifier];
                int grpID = node->grpID;
                for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                    if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                        if (node->children[childIndex]->msaIdx.size() == 0) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                        tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                        util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                        break;
                    }
                }
            }
        }
        for (auto k: tree->allNodes[n.first->identifier]->msaIdx) std::cout << k << ',';
        std::cout << '\n';
        for (auto k: tree->allNodes[n.second->identifier]->msaIdx) std::cout << k << ',';
        std::cout << '\n';
    }

    int numBlocks = 1024; 
    int blockSize = 512;

    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        // int32_t qryLen = tree->allNodes[n.second->identifier]->msa[0].size();
        // int32_t refLen = tree->allNodes[n.first->identifier]->msa[0].size();
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    std::cout << "Max Length:" << seqLen  << " NowStore: " << util->nowStore << '\n';
    int round = nodes.size() / numBlocks + 1;
    for (int r = 0; r < round; ++r) {
        int alignSize = (nodes.size() - r*numBlocks) < numBlocks ? (nodes.size() - r*numBlocks) : numBlocks;
        if (alignSize == 0) break;
        // store all sequences to array
        int32_t seqNum = 0;
        int32_t pairNum = alignSize;
        // std::vector<std::string> seqs;
        // std::vector<std::vector<uint16_t>> freq;
        std::vector<uint16_t*> freq;
        std::vector<std::pair<int32_t, int32_t>> seqIdx;
        std::vector<std::pair<int32_t, int32_t>> len;
        // store info to array 
        auto freqStart = std::chrono::high_resolution_clock::now();
        for (int n = 0; n < alignSize; ++n) {
            int32_t nIdx = n + r*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            refIdx = seqNum;
            
            // std::vector<uint16_t> temp;
            // for (int i = 0; i < 12*seqLen; ++i) temp.push_back(0);
            uint16_t *temp = new uint16_t[12*seqLen]; 
            for (int i = 0; i < 12*seqLen; ++i) temp[i]=0;
            // assert(temp.size() == 12*seqLen);
            // tbb::blocked_range<int> rangeRef(0, refLen);
            // for (auto seq: tree->allNodes[nodes[nIdx].first->identifier]->msa) {
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                int64_t start = util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][start+s] == 'A' || util->seqBuf[storage][start+s] == 'a') temp[6*s+0]+=1;
                    else if (util->seqBuf[storage][start+s] == 'C' || util->seqBuf[storage][start+s] == 'c') temp[6*s+1]+=1;
                    else if (util->seqBuf[storage][start+s] == 'G' || util->seqBuf[storage][start+s] == 'g') temp[6*s+2]+=1;
                    else if (util->seqBuf[storage][start+s] == 'T' || util->seqBuf[storage][start+s] == 't') temp[6*s+3]+=1;
                    else if (util->seqBuf[storage][start+s] == 'N' || util->seqBuf[storage][start+s] == 'n') temp[6*s+4]+=1;
                    else                                                                                     temp[6*s+5]+=1;
                }
                });
                ++seqNum;
            }
            qryIdx = seqNum;
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                int64_t start = sIdx*util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][start+s] == 'A' || util->seqBuf[storage][start+s] == 'a') temp[6*(seqLen+s)+0]+=1;
                    else if (util->seqBuf[storage][start+s] == 'C' || util->seqBuf[storage][start+s] == 'c') temp[6*(seqLen+s)+1]+=1;
                    else if (util->seqBuf[storage][start+s] == 'G' || util->seqBuf[storage][start+s] == 'g') temp[6*(seqLen+s)+2]+=1;
                    else if (util->seqBuf[storage][start+s] == 'T' || util->seqBuf[storage][start+s] == 't') temp[6*(seqLen+s)+3]+=1;
                    else if (util->seqBuf[storage][start+s] == 'N' || util->seqBuf[storage][start+s] == 'n') temp[6*(seqLen+s)+4]+=1;
                    else                                                                     temp[6*(seqLen+s)+5]+=1;
                }
                });
                ++seqNum;
            }

            seqIdx.push_back(std::make_pair(refIdx, qryIdx));
            len.push_back(std::make_pair(refLen, qryLen));
            freq.push_back(temp);
        }
        std::cout << 0 << '\n';
        for (int ppp = 0; ppp < 20; ++ppp){
            std::cout << ppp << " : ";
            for (int i = 0; i < 70; ++i) std::cout << util->seqBuf[0][ppp*util->memLen+i];
            std::cout << '\n';
        }
        std::cout << 1 << '\n';
        
        for (int ppp = 0; ppp < 20; ++ppp){
            std::cout << ppp << " : ";
            for (int i = 0; i < 70; ++i) std::cout << util->seqBuf[1][ppp*util->memLen+i];
            std::cout << '\n';
        }
        // for (int i = 0; i < 70; ++i) std::cout << freq[0][6*i];
        // std::cout << '\n';

        auto freqEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds freqTime = freqEnd -freqStart;
        printf("Calculate frequency time : %d ms\n",  (freqTime.count() / 1000000));
        // Malloc
        uint16_t* hostFreq = (uint16_t*)malloc(12*seqLen * pairNum * sizeof(uint16_t));
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
        uint16_t* deviceFreq;
        int8_t* deviceAln;
        int32_t* deviceLen;
        int32_t* deviceAlnLen;
        int32_t* deviceSeqInfo;
        int16_t* deviceParam;
        auto kernelStart = std::chrono::high_resolution_clock::now();
        cudaMalloc((void**)&deviceFreq, 12*seqLen * pairNum * sizeof(uint16_t));
        cudaMalloc((void**)&deviceAln, 2*seqLen * pairNum * sizeof(int8_t));
        cudaMalloc((void**)&deviceLen, 2*pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceAlnLen, pairNum * sizeof(int32_t));
        cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
        cudaMalloc((void**)&deviceParam, 6 * sizeof(int16_t));
        // Copy to device
        cudaMemcpy(deviceFreq, hostFreq, 12*seqLen * pairNum * sizeof(uint16_t), cudaMemcpyHostToDevice);
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
        auto reAlnStart = std::chrono::high_resolution_clock::now();
        for (int k = 0; k < pairNum; ++k) {
            // std::vector<std::string> alignment;
            int32_t refNum = seqIdx[k].second - seqIdx[k].first;
            int32_t qryNum = (k != pairNum-1) ? seqIdx[k+1].first - seqIdx[k].second : seqNum - seqIdx[k].second;
            int32_t refStart = seqIdx[k].first;
            int32_t qryStart = seqIdx[k].second;
            int32_t refIndex = 0;
            int32_t qryIndex = 0;
            // printf("k: %d, refNum: %d, qryNum: %d\n", k, refNum, qryNum);
            // printf("k: %d, length: %d\n", k, hostAlnLen[k]);
            // for (int j = 0; j < qryNum + refNum; ++j) alignment.push_back("");
            int nIdx = k + r*numBlocks;
            printf("k: %d, length: %d, %s\n", k, hostAlnLen[k], nodes[nIdx].first->identifier.c_str());
            
            if (hostAlnLen[k] <= 0) {
                std::vector<std::string> reference, query;
                std::vector<int8_t> aln;
                for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                    std::string s = "";
                    int64_t start = sIdx*util->memLen;
                    int storage = util->seqsStorage[sIdx];
                    for (int c = 0; c < util->seqsLen[nodes[nIdx].first->identifier]; ++c) {
                        s += util->seqBuf[sIdx][start+c];
                    }
                    reference.push_back(s);
                }
                for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                    std::string s = "";
                    int64_t start = sIdx*util->memLen;
                    int storage = util->seqsStorage[sIdx];
                    for (int c = 0; c < util->seqsLen[nodes[nIdx].second->identifier]; ++c) {
                        s += util->seqBuf[storage][start+c];
                    }
                    query.push_back(s);
                }
                
                Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, param.xdrop, param.marker);
                Talco_xdrop::Align (
                    talco_params,
                    reference,
                    query,
                    aln
                );
                for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                    int64_t start = sIdx*util->memLen;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - util->seqsStorage[sIdx];
                    int orgIdx = 0;
                    for (int j = 0; j < aln.size(); ++j) {
                        if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                            util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                            orgIdx++;
                        }
                        else {
                            util->seqBuf[storeTo][start+j] = '-';
                        }
                    }
                    util->seqsLen[nodes[nIdx].first->identifier] = aln.size();
                    util->changeStorage(sIdx);
                }
                for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                    int64_t start = sIdx*util->memLen;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - util->seqsStorage[sIdx];
                    int orgIdx = 0;
                    for (int j = 0; j < aln.size(); ++j) {
                        if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 1) {
                            util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                            orgIdx++;
                        }
                        else {
                            util->seqBuf[storeTo][start+j] = '-';
                        }
                    }
                    util->seqsLen[nodes[nIdx].second->identifier] = aln.size();
                    util->changeStorage(sIdx);
                }
                printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
            }
            else {
                for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                    int64_t start = sIdx*util->memLen;
                    int orgIdx = 0;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - util->seqsStorage[sIdx];
                    for (int j = 0; j < hostAlnLen[k]; ++j) {
                        if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[k*2*seqLen+j] & 0xFFFF) == 2) {
                            util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                            orgIdx++;
                        }
                        else {
                            util->seqBuf[storeTo][start+j] = '-';
                        }
                    }
                    util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[k];
                    util->changeStorage(sIdx);
                }
                for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                    int64_t start = sIdx*util->memLen;
                    int storeFrom = util->seqsStorage[sIdx];
                    int storeTo = 1 - util->seqsStorage[sIdx];
                    int orgIdx = 0;
                    for (int j = 0; j < hostAlnLen[k]; ++j) {
                        if ((hostAln[k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[k*2*seqLen+j] & 0xFFFF) == 1) {
                            util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                            orgIdx++;
                        }
                        else {
                            util->seqBuf[storeTo][start+j] = '-';
                        }
                    }
                    util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[k];
                    util->changeStorage(sIdx);
                }
            }
            tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
            for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) 
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
        }     
        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - reAlnStart;
        printf("Round. %d align %d pairs. reAlnTime: %d ms\n", r, alignSize, reAlnTime.count() / 1000000);
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
*/


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
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        std::cout << "Level:" << n.first->level << ',' << n.second->level << '\n';
        if (n.first->parent == nullptr) {
            tree->allNodes[n.first->identifier]->msa.clear();
            tree->allNodes[n.first->identifier]->msa = tree->allNodes[n.second->identifier]->msa;
            break;
        }
        bool sameLevel = (n.first->level == n.second->level); 
        std::vector<bool> allGapsR, allGapsQ;
        std::vector<std::string> reference, query;
        std::vector<std::string> alignment;
        for (auto s: tree->allNodes[n.first->identifier]->msa) reference.push_back(s);
        for (auto s: tree->allNodes[n.second->identifier]->msa) query.push_back(s);
        int32_t refLen = reference[0].size();
        int32_t qryLen = query[0].size();
        int32_t refNum = reference.size();
        int32_t qryNum = query.size();
        int32_t seqLen = max(refLen, qryLen);
        int32_t refCut = tree->allNodes[n.first->identifier]->refStartPos;
        int32_t qryCut = tree->allNodes[n.second->identifier]->refStartPos;
        // Overlapped sequence storage (up, down) = (self, parent)
        // Ref's level is always less or equal than Qry'level
        // Merge nodes with the same parent, overlapped sequences = down
        // Merge node with its parent, overlapped sequences = ref: up, qry: down
        int32_t overlapNumRef = (n.first->parent == nullptr) ? refNum : (sameLevel) ? refNum - refCut : refCut;
        int32_t overlapNumQry = qryNum - qryCut;
        std::cout << refNum << ',' << qryNum << ',' << refLen << ',' << qryLen << '\n';
        std::cout << refCut << ',' << qryCut << ',' << overlapNumRef << ',' << overlapNumQry << '\n';
        if ((overlapNumRef == qryNum || overlapNumQry == refNum) && tree->allNodes[n.first->identifier]->parent != nullptr) continue; 
        assert(overlapNumRef == overlapNumQry);

        for (int i = 0; i < refNum + qryNum - overlapNumRef; ++i) alignment.push_back("");

        auto calGapStart = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < seqLen; ++i) {
            if (i < refLen) {
                bool allGaps = true;
                int refStart = (sameLevel) ? refCut : 0;
                for (int j = refStart; j < refStart+overlapNumRef; ++j) {
                    if (reference[j][i] != '-') {
                        allGaps = false;
                        break;
                    }
                }
                allGapsR.push_back(allGaps);
            }
            if (i < qryLen) {
                bool allGaps = true;
                int qryStart = qryCut;
                for (int j = qryStart; j < qryStart+overlapNumQry; ++j) {
                    if (query[j][i] != '-') {
                        allGaps = false;
                        break;
                    }
                }
                allGapsQ.push_back(allGaps);
            }
        }
        auto calGapEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds calGapTime = calGapEnd - calGapStart;
        std::cout << "CalGapTime "<< calGapTime.count() / 1000000<< " ms\n";
        // std::cout << std::count (allGapsR.begin(), allGapsR.end(), true) << '\n';
        // std::cout << std::count (allGapsQ.begin(), allGapsQ.end(), true) << '\n';
        int32_t rIdx = 0, qIdx = 0;
        assert(allGapsR.size() == refLen);
        assert(allGapsQ.size() == qryLen);
        auto mergeStart = std::chrono::high_resolution_clock::now();
        while (rIdx < refLen && qIdx < qryLen) {
            if (allGapsR[rIdx] == false && allGapsQ[qIdx] == false) {
                for (size_t i=0; i<qryCut; i++)      alignment[i] += query[i][qIdx];
                for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += reference[i][rIdx];
                qIdx++;rIdx++;
            }
            else if (allGapsR[rIdx] == true && allGapsQ[qIdx] == false) {
                int consecGap = 0;
                int k = rIdx;
                while (allGapsR[k] && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                // std::cout << "R:" << rIdx << " consecGap: " << consecGap << '\n';
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<qryCut; i++)      alignment[i] += '-';
                    for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += reference[i][rIdx];            
                    // for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
                    // for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += '-';
                    // for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx]; ;
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
                
                // std::cout << "Q:" << qIdx << " consecGap: " << consecGap << '\n';
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<qryCut; i++)      alignment[i] += query[i][qIdx];
                    for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += '-';            
                    
                    // for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
                    // for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
                    // for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += '-';
                    qIdx += 1;
                }
                // if (k == qryLen - 1) break;
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (allGapsR[kr] && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                // std::cout << "dR:" << rIdx << " consecGap: " << consecGap << '\n';
                for (size_t g = 0; g < consecGap; ++g) {
                    // for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
                    // for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += '-';
                    // for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx];
                    for (size_t i=0; i<qryCut; i++)      alignment[i] += '-';
                    for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += reference[i][rIdx];            
                    rIdx += 1;
                }
                consecGap = 0;
                while (allGapsQ[kq] && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                // std::cout << "dQ:" << qIdx << " consecGap: " << consecGap << '\n';
                for (size_t g = 0; g < consecGap; ++g) {
                    // for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
                    // for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
                    // for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += '-';
                    for (size_t i=0; i<qryCut; i++)      alignment[i] += query[i][qIdx];
                    for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += '-';            
                    qIdx += 1;
                }
            }
        }
        
        // bool preGap = 0; // 0: ref, 1: qry
        // std::string rr = "", qr = ""; 
        // while (rIdx < refLen && qIdx < qryLen) {
        //     if (allGapsR[rIdx] == false && allGapsQ[qIdx] == false) {
        //         for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += reference[i][rIdx]; 
        //         for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
        //         for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx];
        //         qIdx++;rIdx++;
        //     }
        //     else if (allGapsR[rIdx] == true && allGapsQ[qIdx] == false) {
        //         for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //         for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += '-';
        //         for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx]; ;
        //         rIdx += 1;
        //         preGap = 0;
        //     }
        //     else if (allGapsR[rIdx] == false && allGapsQ[qIdx] == true) {
        //         for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //         for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
        //         for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += '-';
        //         qIdx += 1;
        //         preGap = 1;
        //     }
        //     else {
        //         int consecGapR = 0, consecGapQ = 0;
        //         int kr = rIdx, kq = qIdx;
        //         while (allGapsR[rIdx] && kr < refLen) {
        //             ++consecGapR;
        //             ++kr;
        //         }
        //         while (allGapsQ[qIdx] && kq < qryLen) {
        //             ++consecGapQ;
        //             ++kq;
        //         }
        //         if (!preGap) {
        //             for (size_t g = 0; g < consecGapR; ++g) {
        //                 for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //                 for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += '-';
        //                 for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx];
        //                 rIdx += 1;
        //             }
        //             for (size_t g = 0; g < consecGapQ; ++g) {
        //                 for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //                 for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
        //                 for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += '-';
        //                 qIdx += 1;
        //             }
        //         }
        //         else {
        //             for (size_t g = 0; g < consecGapQ; ++g) {
        //                 for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //                 for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += query[i][qIdx];
        //                 for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += '-';
        //                 qIdx += 1;
        //             }
        //             for (size_t g = 0; g < consecGapR; ++g) {
        //                 for (size_t i=refStart; i<refNum; i++) alignment[i-refStart]              += '-';
        //                 for (size_t i=0; i<qryStart; i++)      alignment[i+parentNumRef]          += '-';
        //                 for (size_t i=0; i<refStart; i++)      alignment[i+parentNumRef+qryStart] += reference[i][rIdx];
        //                 rIdx += 1;
        //             }
        //         }   
        //     }
        // }
        
        printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                for (size_t i=0; i<qryCut; i++)      alignment[i] += '-';
                for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += reference[i][g];            
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                for (size_t i=0; i<qryCut; i++)      alignment[i] += query[i][g];
                for (size_t i=0; i<refNum; i++)      alignment[i+qryCut] += '-';            
            }
        }
        
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        std::cout << "MergeTime "<< mergeTime.count() / 1000000<< " ms\n";
        tree->allNodes[n.first->identifier]->refStartPos = qryCut + refCut;
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

__global__ void calSPScore(char* seqs, int32_t* seqInfo, int64_t* result) {
    int32_t seqNum     = seqInfo[0]; 
    int32_t seqLen     = seqInfo[1]; 
    int32_t numBlocks  = seqInfo[2]; 
    int32_t numThreads = seqInfo[3]; 
    int32_t match      = seqInfo[4]; 
    int32_t mismatch   = seqInfo[5]; 
    int32_t gap        = seqInfo[6]; 

    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    int gs = gridDim.x;
    const int threadNum = 512;
    // int tidx = bx*bs+tx;
    if (bx < numBlocks) {
        __shared__ int64_t columnScore [threadNum];
        if (tx < numThreads) {
            columnScore[tx] = 0;
        }
        __syncthreads();
        for (int l = bx; l < seqLen; l=l+gs) {
            for (int i = 0; i < seqNum-1; ++i) {
                if (tx < numThreads) {
                    for (int j = i + 1 + tx; j < seqNum; j=j+bs) {
                        if (seqs[i*seqLen+l] == 'N' || seqs[j*seqLen+l] == 'N' ||
                           (seqs[i*seqLen+l] == '-' && seqs[j*seqLen+l] == '-' )) {
                            columnScore[tx] += 0;
                        }
                        else if (seqs[i*seqLen+l] == '-' || seqs[j*seqLen+l] == '-') {
                            columnScore[tx] += gap;
                        }
                        else if (seqs[i*seqLen+l] == seqs[j*seqLen+l]) {
                            columnScore[tx] += match;
                        }
                        else {
                            columnScore[tx] += mismatch;
                        }
                    }
                }
                __syncthreads();
            }
        }
        for (uint32_t r = threadNum/2; r > 0; r >>= 1) {
            if (tx < r) {
                columnScore[tx]   = columnScore[tx] + columnScore[tx+r];
            }
            __syncthreads();
        }
        if (tx == 0) result[bx] = columnScore[0];
    }

}

double getSPScore_gpu(std::vector<std::string>& alignment, Params& param) {
    auto scoreStart = std::chrono::high_resolution_clock::now();
    int numBlocks = 1024;
    int blockSize = 512;
    size_t seqNum = alignment.size();
    size_t seqLen = alignment[0].size();
    printf("(Num, Len) - (%d, %d)\n", seqNum, seqLen);
    char*    hostSeqs = (char*)malloc(seqLen * seqNum * sizeof(char));
    int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
    int64_t* hostResult = (int64_t*)malloc(numBlocks * sizeof(int64_t));
    int seqCount = 0;
    for (int j = 0; j < seqLen*seqNum; ++j) { 
        if (j%seqLen < alignment[seqCount].size()) {
            hostSeqs[j] = alignment[seqCount][j%seqLen];
        }
        else hostSeqs[j] = 0;
        if (j%seqLen == seqLen-1) ++seqCount;
    }
    

    hostSeqInfo[0] = seqNum;
    hostSeqInfo[1] = seqLen;
    hostSeqInfo[2] = numBlocks;
    hostSeqInfo[3] = blockSize;
    hostSeqInfo[4] = param.match;
    hostSeqInfo[5] = param.mismatch;
    hostSeqInfo[6] = param.gapExtend;
    for (int i = 0; i < numBlocks; ++i) hostResult[i] = 0;

    char*    deviceSeqs;
    int32_t* deviceSeqInfo;
    int64_t* deviceResult;

    cudaMalloc((void**)&deviceSeqs, seqLen * seqNum * sizeof(char));
    cudaMalloc((void**)&deviceSeqInfo, 7 * sizeof(int32_t));
    cudaMalloc((void**)&deviceResult, numBlocks * sizeof(int64_t));

    cudaMemcpy(deviceSeqs, hostSeqs, seqLen * seqNum * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceSeqInfo, hostSeqInfo, 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceResult, hostResult, numBlocks * sizeof(int64_t), cudaMemcpyHostToDevice);

    calSPScore<<<numBlocks, blockSize>>>(
        deviceSeqs, 
        deviceSeqInfo,
        deviceResult
    );
    cudaDeviceSynchronize();
    cudaMemcpy(hostResult, deviceResult, numBlocks * sizeof(int64_t), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    int64_t gpuScore = 0;
    for (int i = 0; i < numBlocks; ++i) {
        // std::cout << "Block: "<<i<<", score: " << hostResult[i] << '\n';
        gpuScore += hostResult[i];
    }
    double score = static_cast<double>(gpuScore)/static_cast<double>(seqNum*(seqNum-1)/2);
    // std::cout << "GPU score: " << score << '\n';
    auto scoreEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scoreTime = scoreEnd - scoreStart;
    
    printf("GPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
    
    cudaFree(deviceSeqs);
    cudaFree(deviceSeqInfo);
    cudaFree(deviceResult);
    free(hostSeqs);
    free(hostSeqInfo);
    free(hostResult);
    
    return score;
} 

double getSPScore_cpu(std::vector<std::string>& alignment, Params& param) {
    auto scoreStart = std::chrono::high_resolution_clock::now();
    size_t seqNum = alignment.size();
    size_t seqLen = alignment[0].size();
    printf("(Num, Len) - (%d, %d)\n", seqNum, seqLen);
    double score = 0;
    for (int l = 0; l < seqLen; l++) {
        double columnScore = 0;
        for (int i = 0; i < seqNum - 1; ++i) {
            for (int j = i + 1; j < seqNum; ++j) {
                if (alignment[i][l] == 'N' || alignment[j][l] == 'N' ||
                   (alignment[i][l] == '-' && alignment[j][l] == '-' )) {
                    columnScore += 0;
                }
                else if (alignment[i][l] == '-' || alignment[j][l] == '-') {
                    columnScore += static_cast<double>(param.gapExtend);
                }
                else if (alignment[i][l] == alignment[j][l]) {
                    columnScore +=  static_cast<double>(param.match);
                }
                else {
                    columnScore +=  static_cast<double>(param.mismatch);
                }
            }
        }
        int totalPairs = seqNum*(seqNum-1)/2;
        columnScore /= static_cast<double>(totalPairs);
        score += columnScore;
        if (l % 100 == 0) {
            auto scoreEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds scoreTime = scoreEnd - scoreStart;
            printf("CPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
        }
    }
    auto scoreEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds scoreTime = scoreEnd - scoreStart;
    printf("CPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
    return score;
} 




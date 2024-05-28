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

void getPostOrderList(Node* node, std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(node); 
    Node* current; 
  
    while (!s1.empty()) { 
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (auto ch: current->children) {
            if (ch->grpID == current->grpID) {
                s1.push(ch);
            }      
        }
    } 
    return;
}

void getPreOrderList(Node* node, std::vector<Node*>& preList) {
    if (node == NULL) return;
 
    std::stack<Node*> s1;
    s1.push(node);
    Node* current;
 
    while (s1.empty() == false) {
        current = s1.top();
        preList.push_back(current);
        s1.pop();
        for (auto ch: current->children) {
            if (ch->grpID == current->grpID) {
                s1.push(ch);
            }      
        }
    }
}

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

// Original version: using std::vector
void msaPostOrderTraversal_gpu_org(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
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
        printf("Preprocessing time : %d ms\n",  (freqTime.count() / 1000000));
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
            std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";
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
            for (auto qIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(qIdx);
            }
            // printf("(refNum, qryNum) = (%d, %d)\n", refNum, qryNum);
            // std::cout << "Finish: " << tree->allNodes[nodes[nIdx].first->identifier]->identifier << ',' << tree->allNodes[nodes[nIdx].second->identifier]->identifier << '\n';
        }     
        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - reAlnStart;
        printf("Alignment Time: %d ms\n", r, pairNum, reAlnTime.count() / 1000000);
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


void createOverlapMSA(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
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
        printf("Preprocessing time : %d ms\n",  (freqTime.count() / 1000000));
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
            std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";
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
                tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                printf("CPU fallback on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
            }
            else {
                std::vector<int8_t> aln;
                for (int j = 0; j < hostAlnLen[k]; ++j) {
                    aln.push_back(hostAln[k*2*seqLen+j]);
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
                    tree->allNodes[nodes[nIdx].first->identifier]->msaAln = aln;
                }
            }
            // tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
            // tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
            // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;
            tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
            tree->allNodes[nodes[nIdx].first->identifier]->msa.push_back(nodes[nIdx].first->identifier);
            tree->allNodes[nodes[nIdx].second->identifier]->msa.clear();
            tree->allNodes[nodes[nIdx].second->identifier]->msa.push_back(nodes[nIdx].second->identifier);
            
        }     
        
        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - reAlnStart;
        printf("Alignment Time: %d ms\n", reAlnTime.count() / 1000000);
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


void transitivityMerge_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util)
{
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        // std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        // std::cout << "Level:" << n.first->level << ',' << n.second->level << '\n';
        if (n.first->parent == nullptr) {
            tree->allNodes[n.first->identifier]->msa.clear();
            tree->allNodes[n.first->identifier]->msa = tree->allNodes[n.second->identifier]->msa;
            tree->allNodes[n.first->identifier]->msaIdx.clear();
            tree->allNodes[n.first->identifier]->msaIdx = tree->allNodes[n.second->identifier]->msaIdx;
            int seqLen = tree->allNodes[n.first->identifier]->msa[0].size();
            util->seqsLen[n.first->identifier] = seqLen;
            util->memCheck(seqLen);
            // std::cout << util->memNum << "  " << util->memLen << "  " << util->memLen*util->memNum << '\n';
            
            for (int i = 0; i < tree->allNodes[n.first->identifier]->msa.size(); ++i) {
                int sIdx = tree->allNodes[n.first->identifier]->msaIdx[i];
                int len = tree->allNodes[n.first->identifier]->msa[i].size();
                
                int64_t start = sIdx*util->memLen;
                // std::cout << i << "  " << sIdx << "  " << start << '\n';
                int storeTo = 1 - util->seqsStorage[sIdx];
                for (int j = 0; j < len; ++j) {
                    util->seqBuf[storeTo][start+j] = tree->allNodes[n.first->identifier]->msa[i][j];
                }
                util->changeStorage(sIdx);
            }
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
        // std::cout << refNum << ',' << qryNum << ',' << refLen << ',' << qryLen << '\n';
        // std::cout << refCut << ',' << qryCut << ',' << overlapNumRef << ',' << overlapNumQry << '\n';
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
        // std::cout << "CalGapTime "<< calGapTime.count() / 1000 << " us\n";
        // std::cout << std::count (allGapsR.begin(), allGapsR.end(), true) << '\n';
        // std::cout << std::count (allGapsQ.begin(), allGapsQ.end(), true) << '\n';
        int32_t rIdx = 0, qIdx = 0;
        assert(allGapsR.size() == refLen);
        assert(allGapsQ.size() == qryLen);

        // std::cout << allGapsR.size() << ':';
        // for (int i = 0; i < allGapsR.size(); ++i) {
        //     if (allGapsR[i])
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // std::cout << allGapsQ.size() << ':';
        // for (int i = 0; i < allGapsQ.size(); ++i) {
        //     if (allGapsQ[i])
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';

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
        // printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
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
        // std::cout << "Len: " << alignment[0].size() << '\n';

        std::vector<int> tempIdx;
        for (int q = 0; q < qryCut; ++q) tempIdx.push_back(tree->allNodes[n.second->identifier]->msaIdx[q]);
        for (int r = 0; r < refNum; ++r) tempIdx.push_back(tree->allNodes[n.first->identifier]->msaIdx[r]);
        tree->allNodes[n.first->identifier]->msaIdx.clear();
        tree->allNodes[n.first->identifier]->msaIdx = tempIdx;
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        // std::cout << "MergeTime "<< mergeTime.count() / 1000000<< " ms\n";
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



void transitivityMerge_cpu_mod(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util)
{
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        // std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        // std::cout << "Level:" << n.first->level << ',' << n.second->level << '\n';
        if (n.first->parent == nullptr) {
            // tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
            for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
            // tree->allNodes[n.first->identifier]->msaIdx.clear();
            // tree->allNodes[n.first->identifier]->msaIdx = tree->allNodes[n.second->identifier]->msaIdx;
            std::vector<int8_t> rootAln;
            for (int i = 0; i < tree->allNodes[n.second->identifier]->msaAln.size(); ++i) {
                if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 2) rootAln.push_back(1);
                else if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 1) rootAln.push_back(2);
                else rootAln.push_back(tree->allNodes[n.second->identifier]->msaAln[i]);
            }
            tree->allNodes[n.first->identifier]->msaAln = rootAln;
            
            int seqLen = tree->allNodes[n.first->identifier]->msaAln.size();
            util->seqsLen[n.first->identifier] = seqLen;
            util->memCheck(seqLen);
            break;
        }
        int8_t refGap = (newtree->allNodes[n.first->identifier]->level == newtree->allNodes[n.second->identifier]->level) ? 2 : 1; 
        std::vector<std::vector<int8_t>> refAln, qryAln;
        std::vector<std::vector<int8_t>> refNewAln, qryNewAln;
        // refAln.push_back(tree->allNodes[n.first->identifier]->msaAln);
        // qryAln.push_back(tree->allNodes[n.second->identifier]->msaAln);
        for (auto id: tree->allNodes[n.first->identifier]->msa)  refAln.push_back(tree->allNodes[id]->msaAln);
        for (auto id: tree->allNodes[n.second->identifier]->msa) qryAln.push_back(tree->allNodes[id]->msaAln);
        
        std::vector<int8_t> refOpAln = tree->allNodes[n.first->identifier]->msaAln;
        std::vector<int8_t> qryOpAln = tree->allNodes[n.second->identifier]->msaAln;
        int32_t refLen = refOpAln.size();
        int32_t qryLen = qryOpAln.size();
        int32_t refNum = refAln.size();
        int32_t qryNum = qryAln.size();
        int32_t seqLen = max(refLen, qryLen);
        for (int i = 0; i < refNum; ++i) {
            std::vector<int8_t> temp;
            refNewAln.push_back(temp);
        }
        for (int i = 0; i < qryNum; ++i) {
            std::vector<int8_t> temp;
            qryNewAln.push_back(temp);
        }
        // std::cout << refLen << ':';
        // for (int i = 0; i < refLen; ++i) {
        //     if ((refOpAln[i] & 0XFFFF) == refGap || (refOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // std::cout << qryLen << ':';
        // for (int i = 0; i < qryLen; ++i) {
        //     if ((qryOpAln[i] & 0XFFFF) == 2 || (qryOpAln[i] & 0XFFFF) == 3) 
        //         std::cout << i << ',';
        // }
        // std::cout << '\n';
        // aln = 0: match
        // aln = 1: ref gap
        // aln = 2: qry gap
        // aln = 3: permenant gap
        int32_t rIdx = 0, qIdx = 0;
        auto mergeStart = std::chrono::high_resolution_clock::now();
        
        while (rIdx < refLen && qIdx < qryLen) {
            if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                qIdx++;rIdx++;
            }
            else if (((refOpAln[rIdx] & 0xFFFF) == refGap || (refOpAln[rIdx] & 0xFFFF) == 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                int consecGap = 0;
                int k = rIdx;
                while (((refOpAln[k] & 0xFFFF) == refGap || (refOpAln[k] & 0xFFFF) == 3) && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
            }
            else if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) == 2 || (qryOpAln[qIdx] & 0xFFFF) == 3)) {
                int consecGap = 0;
                int k = qIdx;
                while (((qryOpAln[k] & 0xFFFF) == 2 || (qryOpAln[k] & 0xFFFF) == 3) && k < qryLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (((refOpAln[kr] & 0xFFFF) == refGap || (refOpAln[kr] & 0xFFFF) == 3) && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][rIdx]);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);
                    rIdx += 1;
                }
                consecGap = 0;
                while (((qryOpAln[kq] & 0xFFFF) == 2 || (qryOpAln[kq] & 0xFFFF) == 3) && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                    for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][qIdx]);
                    qIdx += 1;
                }
            }
        }
        // printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(refAln[i][g]);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(3);    
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                for (size_t i=0; i<refNum; i++)      refNewAln[i].push_back(3);
                for (size_t i=0; i<qryNum; i++)      qryNewAln[i].push_back(qryAln[i][g]);
            }
        }
        assert (refNewAln[0].size() == qryNewAln[0].size());
        for (int i = 0; i < tree->allNodes[n.first->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.first->identifier]->msa[i];
            tree->allNodes[id]->msaAln = refNewAln[i];
        } 
        for (int i = 0; i < tree->allNodes[n.second->identifier]->msa.size(); ++i) {
            std::string id = tree->allNodes[n.second->identifier]->msa[i];
            tree->allNodes[id]->msaAln = qryNewAln[i];
        } 
        for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
        // tree->allNodes[n.first->identifier]->msaAln = refNewAln[0];
    }
    
    // for (auto id: tree->root->msa) std::cout << id << ',' << tree->allNodes[id]->msaAln.size() << '\n';
    // std::cout << '\n';
    for (auto id: tree->root->msa) {
        std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
        // std::cout << id << '\n';
        // for (int a = 0; a < aln.size(); ++a) std::cout << (aln[a] & 0xFFFF);
        // std::cout << '\n';
        for (auto sIdx: tree->allNodes[id]->msaIdx) {
            int64_t start = sIdx*util->memLen;
            int orgIdx = 0;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            for (int j = 0; j < aln.size(); ++j) {
                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                    util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                    orgIdx++;
                }
                else {
                    util->seqBuf[storeTo][start+j] = '-';
                }
            }
            util->seqsLen[id] = aln.size();
            util->changeStorage(sIdx);
        }
        if (id != tree->root->identifier) {
            for (auto Idx: tree->allNodes[id]->msaIdx) {
                tree->root->msaIdx.push_back(Idx);
            }
        }
    }
    tree->root->msa.clear();
    


    return;
}


/*
void transitivityMerge_cpu_mod(Tree* tree, Tree* newtree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util)
{
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        // std::cout << "Level:" << n.first->level << ',' << n.second->level << '\n';
        if (n.first->parent == nullptr) {
            // tree->allNodes[n.first->identifier]->msa.push_back(n.first->identifier);
            for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
            // tree->allNodes[n.first->identifier]->msaIdx.clear();
            // tree->allNodes[n.first->identifier]->msaIdx = tree->allNodes[n.second->identifier]->msaIdx;
            std::vector<int8_t> rootAln;
            for (int i = 0; i < tree->allNodes[n.second->identifier]->msaAln.size(); ++i) {
                if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 2) rootAln.push_back(1);
                else if ((tree->allNodes[n.second->identifier]->msaAln[i] & 0XFFFF) == 1) rootAln.push_back(2);
                else rootAln.push_back(tree->allNodes[n.second->identifier]->msaAln[i]);
            }
            tree->allNodes[n.first->identifier]->msaAln = rootAln;
            
            int seqLen = tree->allNodes[n.first->identifier]->msaAln.size();
            util->seqsLen[n.first->identifier] = seqLen;
            util->memCheck(seqLen);
            break;
        }
        int8_t refGap = (newtree->allNodes[n.first->identifier]->level == newtree->allNodes[n.second->identifier]->level) ? 2 : 1; 
        std::vector<int8_t> refAln, qryAln;
        std::vector<int8_t> refNewAln, qryNewAln;
        refAln = tree->allNodes[n.first->identifier]->msaAln;
        qryAln = tree->allNodes[n.second->identifier]->msaAln;
        // for (auto id: tree->allNodes[n.first->identifier]->msa)  refAln.push_back(tree->allNodes[id]->msaAln);
        // for (auto id: tree->allNodes[n.second->identifier]->msa) qryAln.push_back(tree->allNodes[id]->msaAln);
        
        std::vector<int8_t> refOpAln = tree->allNodes[n.first->identifier]->msaAln;
        std::vector<int8_t> qryOpAln = tree->allNodes[n.second->identifier]->msaAln;
        int32_t refLen = refOpAln.size();
        int32_t qryLen = qryOpAln.size();
        // int32_t refNum = refAln.size();
        // int32_t qryNum = qryAln.size();
        int32_t seqLen = max(refLen, qryLen);
        std::cout << refLen << ':';
        for (int i = 0; i < refLen; ++i) {
            if ((refOpAln[i] & 0XFFFF) == refGap || (refOpAln[i] & 0XFFFF) == 3) 
                std::cout << i << ',';
        }
        std::cout << '\n';
        std::cout << qryLen << ':';
        for (int i = 0; i < qryLen; ++i) {
            if ((qryOpAln[i] & 0XFFFF) == 2 || (qryOpAln[i] & 0XFFFF) == 3) 
                std::cout << i << ',';
        }
        std::cout << '\n';
        // aln = 0: match
        // aln = 1: ref gap
        // aln = 2: qry gap
        // aln = 3: permenant gap
        int32_t rIdx = 0, qIdx = 0;
        auto mergeStart = std::chrono::high_resolution_clock::now();
        
        
        while (rIdx < refLen && qIdx < qryLen) {
            if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                refNewAln.push_back(refAln[rIdx]);
                qryNewAln.push_back(qryAln[qIdx]);
                qIdx++;rIdx++;
            }
            else if (((refOpAln[rIdx] & 0xFFFF) == refGap || (refOpAln[rIdx] & 0xFFFF) == 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) != 2 && (qryOpAln[qIdx] & 0xFFFF) != 3)) {
                int consecGap = 0;
                int k = rIdx;
                while (((refOpAln[k] & 0xFFFF) == refGap || (refOpAln[k] & 0xFFFF) == 3) && k < refLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    refNewAln.push_back(refAln[rIdx]);
                    qryNewAln.push_back(3);
                    tree->allNodes[n.second->identifier]->addGapPos.push_back(qIdx+g);
                    rIdx += 1;
                }
            }
            else if (((refOpAln[rIdx] & 0xFFFF) != refGap && (refOpAln[rIdx] & 0xFFFF) != 3)  &&
                     ((qryOpAln[qIdx] & 0xFFFF) == 2 || (qryOpAln[qIdx] & 0xFFFF) == 3)) {
                int consecGap = 0;
                int k = qIdx;
                while (((qryOpAln[k] & 0xFFFF) == 2 || (qryOpAln[k] & 0xFFFF) == 3) && k < qryLen) {
                    ++consecGap;
                    ++k;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    refNewAln.push_back(3);
                    qryNewAln.push_back(qryAln[qIdx]);
                    tree->allNodes[n.first->identifier]->addGapPos.push_back(rIdx+g);
                    qIdx += 1;
                }
            }
            else {
                int consecGap = 0;
                int kr = rIdx, kq = qIdx;
                while (((refOpAln[kr] & 0xFFFF) == refGap || (refOpAln[kr] & 0xFFFF) == 3) && kr < refLen) {
                    ++consecGap;
                    ++kr;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    refNewAln.push_back(refAln[rIdx]);
                    qryNewAln.push_back(3);
                    tree->allNodes[n.second->identifier]->addGapPos.push_back(qIdx+g);
                    rIdx += 1;
                }
                consecGap = 0;
                while (((qryOpAln[kq] & 0xFFFF) == 2 || (qryOpAln[kq] & 0xFFFF) == 3) && kq < qryLen) {
                    ++consecGap;
                    ++kq;
                }
                for (size_t g = 0; g < consecGap; ++g) {
                    refNewAln.push_back(3);
                    qryNewAln.push_back(qryAln[qIdx]);
                    tree->allNodes[n.first->identifier]->addGapPos.push_back(rIdx+g);
                    qIdx += 1;
                }
            }
        }
        // printf("rIdx:%d, qIdx:%d, refLen:%d, qryLen:%d, alnLen: %d\n", rIdx, qIdx, refLen, qryLen, alignment[0].size());
        if (rIdx < refLen) {
            for (size_t g = rIdx; g < refLen; ++g) {
                refNewAln.push_back(refAln[g]);
                qryNewAln.push_back(3);    
                tree->allNodes[n.second->identifier]->addGapPos.push_back(qIdx+g);
            }
        }
        if (qIdx < qryLen) {
            for (size_t g = qIdx; g < qryLen; ++g) {
                refNewAln.push_back(3);
                qryNewAln.push_back(qryAln[g]);
                tree->allNodes[n.first->identifier]->addGapPos.push_back(rIdx+g);
            }
        }
        assert (refNewAln.size() == qryNewAln.size());
        tree->allNodes[n.first->identifier]->msaAln = refNewAln;
        tree->allNodes[n.second->identifier]->msaAln = qryNewAln;
        
        for (auto id: tree->allNodes[n.second->identifier]->msa) tree->allNodes[n.first->identifier]->msa.push_back(id);
        // tree->allNodes[n.first->identifier]->msaAln = refNewAln[0];
    }
    
    // for (auto id: tree->root->msa) std::cout << id << ',' << tree->allNodes[id]->msaAln.size() << '\n';
    // std::cout << '\n';

    std::vector<Node*> preOrderList;
    getPreOrderList(newtree->root, preOrderList);

    for (auto idd: preOrderList) std::cout << idd->identifier << ',';
    std::cout << '\n';
    for (int i = 0; i < preOrderList.size(); ++i) {
        for (auto ch: preOrderList[i]->children) {
            std::cout << ch->identifier << '\n';
            std::vector<int8_t> newAln, oldAln = tree->allNodes[ch->identifier]->msaAln;
            int oldIdx = 0;
            std::cout << "Aln Len: " << oldAln.size() << ',';
            for (int j = 0; j < tree->allNodes[preOrderList[i]->identifier]->msaAln.size(); ++j) {
                if ((tree->allNodes[preOrderList[i]->identifier]->msaAln[j] & 0XFFFF) == 3 && ((oldAln[oldIdx] & 0XFFF) != 3 && (oldAln[oldIdx] & 0XFFF) != 1)) {
                    newAln.push_back(3);
                }
                else {
                    newAln.push_back(oldAln[oldIdx]);
                    ++oldIdx;
                }
            }
            std::cout << "oldIdx: " << oldIdx << '\n';
            tree->allNodes[ch->identifier]->msaAln = newAln;
        }
    }
    


    for (auto id: tree->root->msa) {
        std::vector<int8_t> aln = tree->allNodes[id]->msaAln;
        // std::cout << id << '\n';
        // for (int a = 0; a < aln.size(); ++a) std::cout << (aln[a] & 0xFFFF);
        // std::cout << '\n';
        for (auto sIdx: tree->allNodes[id]->msaIdx) {
            int64_t start = sIdx*util->memLen;
            int orgIdx = 0;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            for (int j = 0; j < aln.size(); ++j) {
                if ((aln[j] & 0xFFFF) == 0 || (aln[j] & 0xFFFF) == 2) {
                    util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                    orgIdx++;
                }
                else {
                    util->seqBuf[storeTo][start+j] = '-';
                }
            }
            util->seqsLen[id] = aln.size();
            util->changeStorage(sIdx);
        }
        if (id != tree->root->identifier) {
            for (auto Idx: tree->allNodes[id]->msaIdx) {
                tree->root->msaIdx.push_back(Idx);
            }
        }
    }
    tree->root->msa.clear();
    


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

// Modified version: malloc 2 memory

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
        // for (auto k: tree->allNodes[n.first->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
        // for (auto k: tree->allNodes[n.second->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
    }

    int numBlocks = 1024; 
    const int blockSize = 512;

    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
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
        std::vector<uint16_t*> freq;
        std::vector<std::pair<int32_t, int32_t>> seqIdx;
        std::vector<int32_t> seqIdxArr;
        std::vector<std::pair<int32_t, int32_t>> len;
        // for (int n = 0; n < alignSize; ++n) {
        //     freq.push_back(nullptr);
        //     // seqIdx.push_back(std::make_pair(0, 0));
        //     seqIdxArr.push_back(0);
        //     seqIdxArr.push_back(0);
        //     len.push_back(std::make_pair(0, 0));
        // }
        // store info to array 
        auto freqStart = std::chrono::high_resolution_clock::now();
        // tbb::parallel_for(tbb::blocked_range<int>(0, alignSize), [&](tbb::blocked_range<int> range) {
        // for (int n = range.begin(); n < range.end(); ++n) {
        for (int n = 0; n < alignSize; ++n) {
            int32_t nIdx = n + r*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t seqNum = 0;
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
                int64_t start = sIdx*util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                // for (int s = 0; s < refLen; ++s) {
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
            // refIdx = seqNum;
            qryIdx = seqNum;
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                int64_t start = sIdx*util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                // for (int s = 0; s < qryLen; ++s) {
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
            // qryIdx = seqNum;
            seqIdx.push_back(std::make_pair(refIdx, qryIdx));
            len.push_back(std::make_pair(refLen, qryLen));
            freq.push_back(temp);
            // seqIdx[n] = std::make_pair(refIdx, qryIdx);
            // seqIdxArr[n*2] = refIdx;
            // seqIdxArr[n*2+1] = qryIdx;
            // len[n] = std::make_pair(refLen, qryLen);
            // freq[n] = temp;
        }
        // });

        // int seqNum = 0;
        // seqIdx.push_back(std::make_pair(0, seqIdxArr[0]));
        // for (int s = 1; s < seqIdxArr.size(); s+=2) {
        //     seqNum += seqIdxArr[s];
        //     seqIdx.push_back(std::make_pair(seqIdxArr[s], seqIdxArr[s+1]));
        // }
        


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

        printf("Before kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
            deviceFreq,
            deviceAln, 
            deviceLen,
            deviceAlnLen,
            deviceSeqInfo, 
            deviceParam
        );

        cudaDeviceSynchronize();
        printf("After kernel %s\n", cudaGetErrorString(cudaGetLastError()));
        // Copy to host
        cudaMemcpy(hostAln, deviceAln, 2*seqLen * pairNum * sizeof(int8_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(hostAlnLen, deviceAlnLen, pairNum * sizeof(int32_t), cudaMemcpyDeviceToHost);
        auto kernelEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
        
        if (round > 1) {
            printf("Round. %d align %d pairs. GPU KernelTime: %d ms\n", r, alignSize, kernelTime.count() / 1000000);
        }
        else {
            std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";
        }
        int maxAlnLen = 0;
        for (int k = 0; k < pairNum; ++k) {
            if (hostAlnLen[k] > maxAlnLen) maxAlnLen = hostAlnLen[k];
        }
        util->memCheck(maxAlnLen);

        auto reAlnStart = std::chrono::high_resolution_clock::now();

        tbb::parallel_for(tbb::blocked_range<int>(0, pairNum), [&](tbb::blocked_range<int> range) {
        for (int k = range.begin(); k < range.end(); ++k) {
            // std::vector<std::string> alignment;
            int32_t refNum = seqIdx[k].second - seqIdx[k].first;
            int32_t qryNum = (k != pairNum-1) ? seqIdx[k+1].first - seqIdx[k].second : seqNum - seqIdx[k].second;
            int32_t refStart = seqIdx[k].first;
            int32_t qryStart = seqIdx[k].second;
            // int32_t refIndex = 0;
            // int32_t qryIndex = 0;
            // printf("k: %d, refNum: %d, qryNum: %d\n", k, refNum, qryNum);
            // printf("k: %d, length: %d\n", k, hostAlnLen[k]);
            // for (int j = 0; j < qryNum + refNum; ++j) alignment.push_back("");
            int nIdx = k + r*numBlocks;
            // printf("k: %d, length: %d, %s\n", k, hostAlnLen[k], nodes[nIdx].first->identifier.c_str());
            
            if (hostAlnLen[k] <= 0) {
                std::vector<std::string> reference, query;
                std::vector<int8_t> aln;
                for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                    std::string s = "";
                    int64_t start = sIdx*util->memLen;
                    int storage = util->seqsStorage[sIdx];
                    for (int c = 0; c < util->seqsLen[nodes[nIdx].first->identifier]; ++c) {
                        s += util->seqBuf[storage][start+c];
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
            // std::cout << "LenB : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
            //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
            for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) 
                tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
            // tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
            // std::cout << "LenA : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
            //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
            // tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
            // tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
            // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;
            
        }   
        });  

        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - reAlnStart;
        // std::chrono::nanoseconds reAlnTime = reAlnEnd - freqStart - kernelTime;
        printf("Alignment Time: %d ns\n", reAlnTime.count());
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
void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{
    // assign msa to all nodes
    // tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range){
    // for (int k = range.begin(); k < range.end(); ++k) { 
    //     auto n = nodes[k];
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
        // for (auto k: tree->allNodes[n.first->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
        // for (auto k: tree->allNodes[n.second->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
    }
    // });

    int numBlocks = 1024; 
    int blockSize = 512;

    int gpuNum;
    cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    int totalBlocks = numBlocks*gpuNum;
    int round = nodes.size() / totalBlocks + 1;


    
    for (int r = 0; r < round; ++r) {
        if (nodes.size() == totalBlocks*r) break;
        int* alignSize = new int[gpuNum];
        int32_t* seqNum = new int32_t[gpuNum];
        uint16_t** hostFreq = new uint16_t* [gpuNum];
        int8_t**   hostAln = new int8_t* [gpuNum];
        int32_t**  hostLen = new int32_t* [gpuNum];
        int32_t**  hostAlnLen = new int32_t* [gpuNum];
        int32_t**  hostSeqInfo = new int32_t* [gpuNum];
        int16_t**  hostParam =  new int16_t* [gpuNum];
        uint16_t** deviceFreq = new uint16_t* [gpuNum];
        int8_t**   deviceAln = new int8_t* [gpuNum];
        int32_t**  deviceLen = new int32_t* [gpuNum];
        int32_t**  deviceAlnLen = new int32_t* [gpuNum];
        int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
        int16_t**  deviceParam = new int16_t* [gpuNum];
        for (int gn = 0; gn < gpuNum; ++gn) {
            int pairsLeft = nodes.size() - r*totalBlocks - gn*numBlocks;
            if (pairsLeft < 0) alignSize[gn] = 0;
            else if (pairsLeft < numBlocks) alignSize[gn] = pairsLeft;
            else alignSize[gn] = numBlocks;
            seqNum[gn] = 0;
            hostFreq[gn] = (uint16_t*)malloc(12*seqLen * alignSize[gn] * sizeof(uint16_t));
            hostAln[gn] = (int8_t*)malloc(2*seqLen * alignSize[gn] * sizeof(int8_t));
            hostLen[gn] = (int32_t*)malloc(2*alignSize[gn] * sizeof(int32_t));
            hostAlnLen[gn] = (int32_t*)malloc(alignSize[gn] * sizeof(int32_t));
            hostSeqInfo[gn] = (int32_t*)malloc(7 * sizeof(int32_t));
            hostParam[gn] = (int16_t*)malloc(6 * sizeof(int16_t)); 
        }
        // store all sequences to array
        std::vector<std::vector<uint16_t*>> freq;
        std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
        std::vector<std::vector<std::pair<int32_t, int32_t>>> len;
        // store info to array 
        auto freqStart = std::chrono::high_resolution_clock::now();
        for (int gn = 0; gn < gpuNum; ++gn) {
            if (alignSize[gn] == 0) break;
            std::vector<uint16_t*> freqTemp;
            std::vector<std::pair<int32_t, int32_t>> seqIdxTemp;
            std::vector<std::pair<int32_t, int32_t>> lenTemp;
            
            for (int n = 0; n < alignSize[gn]; ++n) {
                int32_t nIdx = n + r*totalBlocks+gn*numBlocks;
                int32_t qryIdx = 0;
                int32_t refIdx = 0;
                int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                refIdx = seqNum[gn];
                uint16_t *temp = new uint16_t[12*seqLen]; 
                for (int i = 0; i < 12*seqLen; ++i) temp[i]=0;
                // assert(temp.size() == 12*seqLen);
                // tbb::blocked_range<int> rangeRef(0, refLen);
                for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                    int64_t start = sIdx*util->memLen;
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
                    seqNum[gn] += 1;
                }
                qryIdx = seqNum[gn];
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
                    seqNum[gn] += 1;
                }
                seqIdxTemp.push_back(std::make_pair(refIdx, qryIdx));
                lenTemp.push_back(std::make_pair(refLen, qryLen));
                freqTemp.push_back(temp);
            }
            freq.push_back(freqTemp);
            len.push_back(lenTemp);
            seqIdx.push_back(seqIdxTemp);
        }

        auto freqEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds freqTime = freqEnd -freqStart;
        printf("Calculate frequency time : %d ms\n",  (freqTime.count() / 1000000));
        
        auto kernelStart = std::chrono::high_resolution_clock::now();
        
        
        tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){
        for (int gn = range.begin(); gn < range.end(); ++gn) {
        // for (int gn = 0; gn < gpuNum; ++gn) {
            // Store Info to host mem
            for (int j = 0; j < 2*alignSize[gn]; ++j) { 
                if (j%2 == 0) hostLen[gn][j] = len[gn][j/2].first;
                else          hostLen[gn][j] = len[gn][j/2].second;
            }
            for (int j = 0; j < alignSize[gn]; ++j) {
                for (int l = 0; l < 12*seqLen; ++l) {
                    hostFreq[gn][12*seqLen*j+l] = freq[gn][j][l];
                }
            }
            for (int j = 0; j < 2*seqLen*alignSize[gn]; ++j) { 
                hostAln[gn][j] = 0;
            }
            for (int j = 0; j < alignSize[gn]; ++j) { 
                hostAlnLen[gn][j] = 0;
            }
            hostSeqInfo[gn][0] = seqLen;
            hostSeqInfo[gn][1] = seqNum[gn];
            hostSeqInfo[gn][2] = alignSize[gn];
            hostSeqInfo[gn][3] = numBlocks;
            hostSeqInfo[gn][4] = blockSize;
            hostSeqInfo[gn][5] = numBlocks;
            hostSeqInfo[gn][6] = blockSize;
            hostParam[gn][0] = param.match;
            hostParam[gn][1] = param.mismatch;
            hostParam[gn][2] = param.gapOpen;
            hostParam[gn][3] = param.gapExtend;
            hostParam[gn][4] = param.xdrop;
            hostParam[gn][5] = param.marker;

            for (int i = 0; i < alignSize[gn]; ++i) delete [] freq[gn][i];

            cudaSetDevice(gn);

            cudaMalloc((void**)&deviceFreq[gn], 12*seqLen * alignSize[gn] * sizeof(uint16_t));
            cudaMalloc((void**)&deviceAln[gn], 2*seqLen * alignSize[gn] * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn], 2*alignSize[gn] * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn], alignSize[gn] * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 7 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn], 6 * sizeof(int16_t));

            cudaMemcpy(deviceFreq[gn], hostFreq[gn], 12*seqLen * alignSize[gn] * sizeof(uint16_t), cudaMemcpyHostToDevice);
            cudaMemcpy(deviceAln[gn], hostAln[gn], 2*seqLen * alignSize[gn] * sizeof(int8_t), cudaMemcpyHostToDevice);
            cudaMemcpy(deviceLen[gn], hostLen[gn], 2*alignSize[gn] * sizeof(int32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(deviceAlnLen[gn], hostAlnLen[gn], alignSize[gn] * sizeof(int32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[gn], 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(deviceParam[gn], hostParam[gn], 6 * sizeof(int16_t), cudaMemcpyHostToDevice);

            alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                deviceFreq[gn],
                deviceAln[gn], 
                deviceLen[gn],
                deviceAlnLen[gn],
                deviceSeqInfo[gn], 
                deviceParam[gn]
            );
            cudaDeviceSynchronize();
            cudaMemcpy(hostAln[gn], deviceAln[gn], 2*seqLen * alignSize[gn] * sizeof(int8_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(hostAlnLen[gn], deviceAlnLen[gn], alignSize[gn] * sizeof(int32_t), cudaMemcpyDeviceToHost);
        }
        });
        
        for (int gn = 0; gn < gpuNum; ++gn) {
            cudaSetDevice(gn);
            cudaDeviceSynchronize();
        }
        auto kernelEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
        
        int totalPairs = 0;
        for (int gn = 0; gn < gpuNum; ++gn) totalPairs += alignSize[gn];
        if (round > 1) {
            printf("Round. %d align %d pairs. GPU KernelTime: %d ms\n", r, totalPairs, kernelTime.count() / 1000000);
        }
        else {
            std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";
        }
        auto reAlnStart = std::chrono::high_resolution_clock::now();
        int maxAlnLen = 0;
        for (int gn = 0; gn < gpuNum; ++gn) {
           for (int k = 0; k <  alignSize[gn]; ++k) {
                if (hostAlnLen[gn][k] > maxAlnLen) maxAlnLen = hostAlnLen[gn][k];
            }
        }
        util->memCheck(maxAlnLen);

        
        // for (int x = 0; x < alignSize[0]; ++x) printf("No. %d, len = %d\n", x, hostAlnLen[0][x]); 
        for (int gn = 0; gn < gpuNum; ++gn) {
            if (alignSize[gn] == 0) break;
            // std::cout << "alnSize: " << alignSize[gn] << '\n';
            tbb::parallel_for(tbb::blocked_range<int>(0, alignSize[gn]), [&](tbb::blocked_range<int> range) {
            // for (int k = 0; k < alignSize[gn]; ++k) {
            for (int k = range.begin(); k < range.end(); ++k) {
                // std::vector<std::string> alignment;
                int32_t refNum = seqIdx[gn][k].second - seqIdx[gn][k].first;
                int32_t qryNum = (k !=  alignSize[gn]-1) ? seqIdx[gn][k+1].first - seqIdx[gn][k].second : seqNum[gn] - seqIdx[gn][k].second;
                int32_t refStart = seqIdx[gn][k].first;
                int32_t qryStart = seqIdx[gn][k].second;
                int32_t refIndex = 0;
                int32_t qryIndex = 0;
                int32_t nIdx = k + r*totalBlocks+gn*numBlocks;
                if (hostAlnLen[gn][k] <= 0) {
                    std::vector<std::string> reference, query;
                    std::vector<int8_t> aln;
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                        std::string s = "";
                        int64_t start = sIdx*util->memLen;
                        int storage = util->seqsStorage[sIdx];
                        for (int c = 0; c < util->seqsLen[nodes[nIdx].first->identifier]; ++c) {
                            s += util->seqBuf[storage][start+c];
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
                    Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, 500, param.marker);
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
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 2) {
                                util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][start+j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                        int64_t start = sIdx*util->memLen;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        int orgIdx = 0;
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 1) {
                                util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][start+j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                }
                tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
                // std::cout << "LenB : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
                //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
                for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) 
                    tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                // tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                // std::cout << "LenA : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
                //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
                // tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
                // tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
                // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;

            }  
            });
        }   
        
        
        auto reAlnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds reAlnTime = reAlnEnd - kernelEnd;
        printf("Alignment Time: %d ns\n", r, alignSize, reAlnTime.count() / 1);
        for (int gn = 0; gn < gpuNum; ++gn) {
            // free memory
            cudaSetDevice(gn);
            cudaFree(deviceFreq[gn]);
            cudaFree(deviceAlnLen[gn]);
            cudaFree(deviceAln[gn]);
            cudaFree(deviceParam[gn]);
            cudaFree(deviceSeqInfo[gn]);
            cudaDeviceSynchronize();
            free(hostFreq[gn]);
            free(hostAlnLen[gn]);
            free(hostAln[gn]);
            free(hostParam[gn]);
            free(hostSeqInfo[gn]);
            // for (auto f: freq[gn]) delete [] f;   
        }
        delete [] alignSize;
        delete [] seqNum;
        delete [] deviceFreq;
        delete [] deviceAlnLen;
        delete [] deviceAln;
        delete [] deviceParam;
        delete [] deviceSeqInfo;
        delete [] hostFreq;
        delete [] hostAlnLen;
        delete [] hostAln;
        delete [] hostParam;
        delete [] hostSeqInfo;
    }
    return;
}
*/


void msaPostOrderTraversal_multigpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util, Params& param)
{

    auto freqStart = std::chrono::high_resolution_clock::now();
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
        // for (auto k: tree->allNodes[n.first->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
        // for (auto k: tree->allNodes[n.second->identifier]->msaIdx) std::cout << k << ',';
        // std::cout << '\n';
    }
    // });

    int numBlocks = 1024; 
    int blockSize = 1024;
    int gpuNum;
    cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    
    // get maximum sequence/profile length 
    int32_t seqLen = 0;
    for (auto n: nodes) {
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
        int32_t tempMax = max(qryLen, refLen);
        seqLen = max(seqLen, tempMax);
    }
    
    int roundGPU = nodes.size() / numBlocks + 1;
    if (nodes.size()%numBlocks == 0) roundGPU -= 1;
    int* alignSize = new int[roundGPU];
    int32_t* seqNum = new int32_t[roundGPU];
    uint16_t** hostFreq = new uint16_t* [roundGPU];
    int8_t**   hostAln = new int8_t* [roundGPU];
    int32_t**  hostLen = new int32_t* [roundGPU];
    int32_t**  hostAlnLen = new int32_t* [roundGPU];
    int32_t**  hostSeqInfo = new int32_t* [roundGPU];
    int16_t**  hostParam =  new int16_t* [roundGPU];
    std::vector<std::vector<uint16_t*>> freq;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> seqIdx;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> len;
    for (int gn = 0; gn < roundGPU; ++gn) {
        // std::cout << "Gn: " << gn << '\n';
        int pairsLeft = nodes.size() - gn*numBlocks;
        if (pairsLeft < numBlocks) alignSize[gn] = pairsLeft;
        else alignSize[gn] = numBlocks;
        seqNum[gn] = 0;
        hostFreq[gn] = (uint16_t*)malloc(12*seqLen * alignSize[gn] * sizeof(uint16_t));
        hostAln[gn] = (int8_t*)malloc(2*seqLen * alignSize[gn] * sizeof(int8_t));
        hostLen[gn] = (int32_t*)malloc(2*alignSize[gn] * sizeof(int32_t));
        hostAlnLen[gn] = (int32_t*)malloc(alignSize[gn] * sizeof(int32_t));
        hostSeqInfo[gn] = (int32_t*)malloc(7 * sizeof(int32_t));
        hostParam[gn] = (int16_t*)malloc(6 * sizeof(int16_t)); 
        // store all sequences to array
        std::vector<uint16_t*> freqTemp;
        std::vector<std::pair<int32_t, int32_t>> seqIdxTemp;
        std::vector<std::pair<int32_t, int32_t>> lenTemp;
        for (int n = 0; n < alignSize[gn]; ++n) {
            int32_t nIdx = n + gn*numBlocks;
            int32_t qryIdx = 0;
            int32_t refIdx = 0;
            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            // if (refLen <= 1) {
            //     std::cout << nodes[nIdx].first->identifier << ',' << nodes[nIdx].second->identifier  << '\n';
            //     printf("X: Gn: %d, n: %d, ref: %d, qry: %d, refNum: %d, qryNum: %d\n", gn, n ,refLen, qryLen, tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size(), tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size());
            //     for (auto idx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx)
            //         std::cout << idx << ',';
            //     std::cout << '\n';
            // }
            // if (qryLen <= 1) {
            //     std::cout << nodes[nIdx].first->identifier << ',' << nodes[nIdx].second->identifier  << '\n';
            //     printf("X: Gn: %d, n: %d, ref: %d, qry: %d, refNum: %d, qryNum: %d\n", gn, n ,refLen, qryLen, tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size(), tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size());
            //     for (auto idx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx)
            //         std::cout << idx << ',';
            //     std::cout << '\n';
            // }
            // std::cout << nodes[nIdx].first->identifier << ',' << nodes[nIdx].second->identifier  << '\n';
            // printf("Gn: %d, n: %d, ref: %d, qry: %d, refNum: %d, qryNum: %d\n", gn, n ,refLen, qryLen, tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size(), tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size());
            refIdx = seqNum[gn];
            uint16_t *temp = new uint16_t[12*seqLen]; 
            for (int i = 0; i < 12*seqLen; ++i) temp[i]=0;
            // assert(temp.size() == 12*seqLen);
            // tbb::blocked_range<int> rangeRef(0, refLen);
            for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) { 
                int64_t start = sIdx*util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, refLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][start+s] == 'A' || util->seqBuf[storage][start+s] == 'a') temp[6*s+0]+=1;
                    else if (util->seqBuf[storage][start+s] == 'C' || util->seqBuf[storage][start+s] == 'c') temp[6*s+1]+=1;
                    else if (util->seqBuf[storage][start+s] == 'G' || util->seqBuf[storage][start+s] == 'g') temp[6*s+2]+=1;
                    else if (util->seqBuf[storage][start+s] == 'T' || util->seqBuf[storage][start+s] == 't' ||
                             util->seqBuf[storage][start+s] == 'U' || util->seqBuf[storage][start+s] == 'u') temp[6*s+3]+=1;
                    else if (util->seqBuf[storage][start+s] == 'N' || util->seqBuf[storage][start+s] == 'n') temp[6*s+4]+=1;
                    else                                                                                     temp[6*s+5]+=1;
                }
                });
                seqNum[gn] += 1;
            }
            qryIdx = seqNum[gn];
            for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) { 
                int64_t start = sIdx*util->memLen;
                int storage = util->seqsStorage[sIdx];
                tbb::parallel_for(tbb::blocked_range<int>(0, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    if      (util->seqBuf[storage][start+s] == 'A' || util->seqBuf[storage][start+s] == 'a') temp[6*(seqLen+s)+0]+=1;
                    else if (util->seqBuf[storage][start+s] == 'C' || util->seqBuf[storage][start+s] == 'c') temp[6*(seqLen+s)+1]+=1;
                    else if (util->seqBuf[storage][start+s] == 'G' || util->seqBuf[storage][start+s] == 'g') temp[6*(seqLen+s)+2]+=1;
                    else if (util->seqBuf[storage][start+s] == 'T' || util->seqBuf[storage][start+s] == 't' ||
                             util->seqBuf[storage][start+s] == 'U' || util->seqBuf[storage][start+s] == 'u') temp[6*(seqLen+s)+3]+=1;
                    else if (util->seqBuf[storage][start+s] == 'N' || util->seqBuf[storage][start+s] == 'n') temp[6*(seqLen+s)+4]+=1;
                    else                                                                     temp[6*(seqLen+s)+5]+=1;
                }
                });
                seqNum[gn] += 1;
            }
            seqIdxTemp.push_back(std::make_pair(refIdx, qryIdx));
            lenTemp.push_back(std::make_pair(refLen, qryLen));
            freqTemp.push_back(temp);
        }
        freq.push_back(freqTemp);
        len.push_back(lenTemp);
        seqIdx.push_back(seqIdxTemp);
    }

    tbb::parallel_for(tbb::blocked_range<int>(0, roundGPU), [&](tbb::blocked_range<int> range){
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            for (int j = 0; j < 2*alignSize[gn]; ++j) { 
                if (j%2 == 0) hostLen[gn][j] = len[gn][j/2].first;
                else          hostLen[gn][j] = len[gn][j/2].second;
            }
            for (int j = 0; j < alignSize[gn]; ++j) {
                for (int l = 0; l < 12*seqLen; ++l) {
                    hostFreq[gn][12*seqLen*j+l] = freq[gn][j][l];
                }
            }
            for (int j = 0; j < 2*seqLen*alignSize[gn]; ++j) { 
                hostAln[gn][j] = 0;
            }
            for (int j = 0; j < alignSize[gn]; ++j) { 
                hostAlnLen[gn][j] = 0;
            }
            hostSeqInfo[gn][0] = seqLen;
            hostSeqInfo[gn][1] = seqNum[gn];
            hostSeqInfo[gn][2] = alignSize[gn];
            hostSeqInfo[gn][3] = numBlocks;
            hostSeqInfo[gn][4] = blockSize;
            hostSeqInfo[gn][5] = numBlocks;
            hostSeqInfo[gn][6] = blockSize;
            hostParam[gn][0] = param.match;
            hostParam[gn][1] = param.mismatch;
            hostParam[gn][2] = param.gapOpen;
            hostParam[gn][3] = param.gapExtend;
            hostParam[gn][4] = param.xdrop;
            hostParam[gn][5] = param.marker;
            
        }
    });
    auto freqEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds freqTime = freqEnd -freqStart;
    printf("Preprocessing time : %d ms\n",  (freqTime.count() / 1000000));        

    auto kernelStart = std::chrono::high_resolution_clock::now();
    uint16_t** deviceFreq = new uint16_t* [gpuNum];
    int8_t**   deviceAln = new int8_t* [gpuNum];
    int32_t**  deviceLen = new int32_t* [gpuNum];
    int32_t**  deviceAlnLen = new int32_t* [gpuNum];
    int32_t**  deviceSeqInfo = new int32_t* [gpuNum];
    int16_t**  deviceParam = new int16_t* [gpuNum];
    // int nowRound = 0;
    std::atomic<int> nowRound;
    nowRound.store(0);
    tbb::parallel_for(tbb::blocked_range<int>(0, gpuNum), [&](tbb::blocked_range<int> range){ 
        for (int gn = range.begin(); gn < range.end(); ++gn) {
            cudaSetDevice(gn);
            int nowMemSize = alignSize[gn];
            cudaMalloc((void**)&deviceFreq[gn], 12*seqLen * alignSize[gn] * sizeof(uint16_t));
            cudaMalloc((void**)&deviceAln[gn], 2*seqLen * alignSize[gn] * sizeof(int8_t));
            cudaMalloc((void**)&deviceLen[gn], 2*alignSize[gn] * sizeof(int32_t));
            cudaMalloc((void**)&deviceAlnLen[gn], alignSize[gn] * sizeof(int32_t));
            cudaMalloc((void**)&deviceSeqInfo[gn], 7 * sizeof(int32_t));
            cudaMalloc((void**)&deviceParam[gn], 6 * sizeof(int16_t));
            while (nowRound < roundGPU) {
                int rn = nowRound.fetch_add(1);
                // int rn = nowRound;
                // nowRound += 1;
                if (alignSize[rn] != nowMemSize) {
                    cudaSetDevice(gn);
                    cudaFree(deviceFreq[gn]);
                    cudaFree(deviceAlnLen[gn]);
                    cudaFree(deviceAln[gn]);
                    cudaFree(deviceParam[gn]);
                    cudaFree(deviceSeqInfo[gn]);
                    cudaDeviceSynchronize();
                    cudaMalloc((void**)&deviceFreq[gn], 12*seqLen * alignSize[rn] * sizeof(uint16_t));
                    cudaMalloc((void**)&deviceAln[gn], 2*seqLen * alignSize[rn] * sizeof(int8_t));
                    cudaMalloc((void**)&deviceLen[gn], 2*alignSize[rn] * sizeof(int32_t));
                    cudaMalloc((void**)&deviceAlnLen[gn], alignSize[rn] * sizeof(int32_t));
                    cudaMalloc((void**)&deviceSeqInfo[gn], 7 * sizeof(int32_t));
                    cudaMalloc((void**)&deviceParam[gn], 6 * sizeof(int16_t));
                }
                cudaMemcpy(deviceFreq[gn], hostFreq[rn], 12*seqLen * alignSize[rn] * sizeof(uint16_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAln[gn], hostAln[rn], 2*seqLen * alignSize[rn] * sizeof(int8_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceLen[gn], hostLen[rn], 2*alignSize[rn] * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceAlnLen[gn], hostAlnLen[rn], alignSize[rn] * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceSeqInfo[gn], hostSeqInfo[rn], 7 * sizeof(int32_t), cudaMemcpyHostToDevice);
                cudaMemcpy(deviceParam[gn], hostParam[rn], 6 * sizeof(int16_t), cudaMemcpyHostToDevice);

                alignGrpToGrp_talco<<<numBlocks, blockSize>>>(
                    deviceFreq[gn],
                    deviceAln[gn], 
                    deviceLen[gn],
                    deviceAlnLen[gn],
                    deviceSeqInfo[gn], 
                    deviceParam[gn]
                );
                cudaDeviceSynchronize();
                cudaMemcpy(hostAln[rn], deviceAln[gn], 2*seqLen * alignSize[rn] * sizeof(int8_t), cudaMemcpyDeviceToHost);
                cudaMemcpy(hostAlnLen[rn], deviceAlnLen[gn], alignSize[rn] * sizeof(int32_t), cudaMemcpyDeviceToHost);
            }
        }
    });

    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(gn);
        cudaDeviceSynchronize();
    }

    auto kernelEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    int totalPairs = 0;
    for (int gn = 0; gn < roundGPU; ++gn) totalPairs += alignSize[gn];
    std::cout << "GPU KernelTime "<< kernelTime.count() / 1000000<< " ms\n";

    auto reAlnStart = std::chrono::high_resolution_clock::now();
    
    int maxAlnLen = 0;
    for (int gn = 0; gn < roundGPU; ++gn) {
       for (int k = 0; k <  alignSize[gn]; ++k) {
            if (hostAlnLen[gn][k] > maxAlnLen) maxAlnLen = hostAlnLen[gn][k];
        }
    }
    util->memCheck(maxAlnLen);
        
    for (int gn = 0; gn < roundGPU; ++gn) {
        if (alignSize[gn] == 0) break;
        // std::cout << "Round: " << gn << " alnSize: " << alignSize[gn] << '\n';
        tbb::parallel_for(tbb::blocked_range<int>(0, alignSize[gn]), [&](tbb::blocked_range<int> range) {
            // for (int k = 0; k < alignSize[gn]; ++k) {
            for (int k = range.begin(); k < range.end(); ++k) {
                // std::vector<std::string> alignment;
                int32_t refNum = seqIdx[gn][k].second - seqIdx[gn][k].first;
                int32_t qryNum = (k !=  alignSize[gn]-1) ? seqIdx[gn][k+1].first - seqIdx[gn][k].second : seqNum[gn] - seqIdx[gn][k].second;
                int32_t refStart = seqIdx[gn][k].first;
                int32_t qryStart = seqIdx[gn][k].second;
                int32_t refIndex = 0;
                int32_t qryIndex = 0;
                int32_t nIdx = k + gn*numBlocks;
                if (hostAlnLen[gn][k] <= -10) {
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    std::vector<int8_t> aln;
                    alignGrpToGrp_traditional
                    (
                        freq[gn][k],
                        seqLen,
                        refLen,
                        qryLen,
                        param,
                        aln
                    );
                    int32_t alnLen = aln.size();
                    util->memCheck(alnLen);
                    std::reverse(aln.begin(), aln.end());
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
                    // for (auto a: aln) std::cout << (a&0xFFFF) << ',';
                    // std::cout << '\n';
                    printf("CPU fallback (traditional global alignment) on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                }
                else if (hostAlnLen[gn][k] <= 0) {
                    std::vector<int8_t> aln;
                    std::vector<std::vector<int>> freqRef;
                    std::vector<std::vector<int>> freqQry;
                    int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
                    int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
                    
                    for (int r = 0; r < refLen; r++) {
                        std::vector<int> temp;
                        for (int f = 0; f < 6; ++f) temp.push_back(freq[gn][k][6*r+f]);
                        freqRef.push_back(temp);
                    }
                    for (int q = 0; q < qryLen; q++) {
                        std::vector<int> temp;
                        for (int f = 0; f < 6; ++f) temp.push_back(freq[gn][k][6*(seqLen+q)+f]);
                        freqQry.push_back(temp);
                    }

                    Talco_xdrop::Params talco_params(param.match, param.mismatch, param.gapOpen, param.gapExtend, 1000, param.marker);
                    Talco_xdrop::Align_freq (
                        talco_params,
                        freqRef,
                        freqQry,
                        aln
                    );
                    util->memCheck(aln.size());
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
                    printf("CPU fallback (TALCO-Xdrop) on No. %d (%s), Alignment Length: %d\n", k, tree->allNodes[nodes[nIdx].first->identifier]->identifier.c_str(), aln.size());
                }
                else {
                    for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                        int64_t start = sIdx*util->memLen;
                        int orgIdx = 0;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 2) {
                                util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][start+j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].first->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                    for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                        int64_t start = sIdx*util->memLen;
                        int storeFrom = util->seqsStorage[sIdx];
                        int storeTo = 1 - util->seqsStorage[sIdx];
                        int orgIdx = 0;
                        for (int j = 0; j < hostAlnLen[gn][k]; ++j) {
                            if ((hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][k*2*seqLen+j] & 0xFFFF) == 1) {
                                util->seqBuf[storeTo][start+j] = util->seqBuf[storeFrom][start+orgIdx];
                                orgIdx++;
                            }
                            else {
                                util->seqBuf[storeTo][start+j] = '-';
                            }
                        }
                        util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][k];
                        util->changeStorage(sIdx);
                    }
                }
                tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
                // std::cout << "LenB : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
                //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
                for (auto q: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) 
                    tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.push_back(q);
                // tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.clear();
                // std::cout << "LenA : " << nodes[nIdx].first->identifier << '(' << tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size() << ')'
                //                       << nodes[nIdx].second->identifier << '(' << tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size() << ")\n";
                // tree->allNodes[nodes[nIdx].first->identifier]->refStartPos = refNum;
                // tree->allNodes[nodes[nIdx].first->identifier]->msa.clear();
                // tree->allNodes[nodes[nIdx].first->identifier]->msa = alignment;

            }  
        });
        for (int i = 0; i < alignSize[gn]; ++i) delete [] freq[gn][i];
    } 
    auto reAlnEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds reAlnTime = reAlnEnd - kernelEnd;
    printf("Alignment Time: %d us\n", reAlnTime.count() / 1000);

    // free memory  
    for (int gn = 0; gn < gpuNum; ++gn) {
        cudaSetDevice(gn);
        cudaFree(deviceFreq[gn]);
        cudaFree(deviceAlnLen[gn]);
        cudaFree(deviceAln[gn]);
        cudaFree(deviceParam[gn]);
        cudaFree(deviceSeqInfo[gn]);
        cudaDeviceSynchronize();  
    }

    for (int gn = 0; gn < roundGPU; ++gn) {
        free(hostFreq[gn]);
        free(hostAlnLen[gn]);
        free(hostAln[gn]);
        free(hostParam[gn]);
        free(hostSeqInfo[gn]);
    }  

    delete [] alignSize;
    delete [] seqNum;
    delete [] deviceFreq;
    delete [] deviceAlnLen;
    delete [] deviceAln;
    delete [] deviceParam;
    delete [] deviceSeqInfo;
    delete [] hostFreq;
    delete [] hostAlnLen;
    delete [] hostAln;
    delete [] hostParam;
    delete [] hostSeqInfo;
    return;
}


/*
void transitivityMerge_cpu(Tree* tree, std::vector<std::pair<Node*, Node*>> nodes, msa::utility* util)
{
    // auto kernelStart = std::chrono::high_resolution_clock::now();

    for (auto n: nodes) {
        std::cout << tree->allNodes[n.first->identifier]->identifier << ',' << tree->allNodes[n.second->identifier]->identifier << '\n';
        std::cout << "Level:" << n.first->level << ',' << n.second->level << '\n';
        if (n.first->parent == nullptr) {
            tree->allNodes[n.first->identifier]->msaIdx.clear();
            tree->allNodes[n.first->identifier]->msaIdx = tree->allNodes[n.second->identifier]->msaIdx;
            break;
        }
        bool sameLevel = (n.first->level == n.second->level); 
        std::vector<bool> allGapsR, allGapsQ;
        // std::vector<std::string> reference, query;
        // std::vector<std::string> alignment;
        // for (auto s: tree->allNodes[n.first->identifier]->msa) reference.push_back(s);
        // for (auto s: tree->allNodes[n.second->identifier]->msa) query.push_back(s);
        // int32_t refLen = reference[0].size();
        // int32_t qryLen = query[0].size();
        int32_t refLen = util->seqsLen[n.first->identifier];
        int32_t qryLen = util->seqsLen[n.second->identifier];
        int32_t refNum = tree->allNodes[n.first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[n.second->identifier]->msaIdx.size();
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

        // for (int i = 0; i < refNum + qryNum - overlapNumRef; ++i) alignment.push_back("");

        auto calGapStart = std::chrono::high_resolution_clock::now();
        
        for (int i = 0; i < seqLen; ++i) {
            if (i < refLen) {
                bool allGaps = true;
                int refStart = (sameLevel) ? refCut : 0;
                for (int j = refStart; j < refStart+overlapNumRef; ++j) {
                    int sIdx = tree->allNodes[n.first->identifier]->msaIdx[j];
                    int64_t start = sIdx*util->memLen;
                    int store = util->seqsStorage[sIdx];
                    if (util->seqBuf[store][i] != '-') {
                        allGaps = false;
                        break;
                    }
                    // if (reference[j][i] != '-') {
                    //     allGaps = false;
                    //     break;
                    // }
                }
                allGapsR.push_back(allGaps);
            }
            if (i < qryLen) {
                bool allGaps = true;
                int qryStart = qryCut;
                for (int j = qryStart; j < qryStart+overlapNumQry; ++j) {
                    int sIdx = tree->allNodes[n.second->identifier]->msaIdx[j];
                    int64_t start = sIdx*util->memLen;
                    int store = util->seqsStorage[sIdx];
                    if (util->seqBuf[store][i] != '-') {
                        allGaps = false;
                        break;
                    }
                    // if (query[j][i] != '-') {
                    //     allGaps = false;
                    //     break;
                    // }
                }
                allGapsQ.push_back(allGaps);
            }
        }
        auto calGapEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds calGapTime = calGapEnd - calGapStart;
        std::cout << "CalGapTime "<< calGapTime.count() / 1000000<< " ms\n";
        // std::cout << std::count (allGapsR.begin(), allGapsR.end(), true) << '\n';
        // std::cout << std::count (allGapsQ.begin(), allGapsQ.end(), true) << '\n';
        // int32_t rIdx = 0, qIdx = 0;
        assert(allGapsR.size() == refLen);
        assert(allGapsQ.size() == qryLen);
        auto mergeStart = std::chrono::high_resolution_clock::now();
        for (auto sIdx: tree->allNodes[n.first->identifier]->msaIdx) {
            int32_t rIdx = 0, qIdx = 0;
            int64_t start = sIdx*util->memLen;
            int newIdx = 0;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            while (rIdx < refLen && qIdx < qryLen) {
                if (allGapsR[rIdx] == false && allGapsQ[qIdx] == false) {
                    util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+rIdx];
                    qIdx++;rIdx++;newIdx++;
                }
                else if (allGapsR[rIdx] == true && allGapsQ[qIdx] == false) {
                    int consecGap = 0;
                    int k = rIdx;
                    while (allGapsR[k] && k < refLen) {
                        ++consecGap;
                        ++k;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+rIdx];
                        rIdx++;newIdx++;
                    }
                }
                else if (allGapsR[rIdx] == false && allGapsQ[qIdx] == true) {
                    int consecGap = 0;
                    int k = qIdx;
                    while (allGapsQ[k] && k < qryLen) {
                        ++consecGap;
                        ++k;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = '-';
                        qIdx++;newIdx++;
                    }
                }
                else {
                    int consecGap = 0;
                    int kr = rIdx, kq = qIdx;
                    while (allGapsR[kr] && kr < refLen) {
                        ++consecGap;
                        ++kr;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+rIdx];
                        rIdx++;newIdx++;
                    }
                    consecGap = 0;
                    while (allGapsQ[kq] && kq < qryLen) {
                        ++consecGap;
                        ++kq;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = '-';
                        qIdx++;newIdx++;
                    }
                }
            }
            if (rIdx < refLen) {
                for (size_t g = rIdx; g < refLen; ++g) {
                    util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+g];
                    newIdx++;          
                }
            }
            if (qIdx < qryLen) {
                for (size_t g = qIdx; g < qryLen; ++g) {
                    util->seqBuf[storeTo][start+newIdx] = '-';
                    newIdx++;   
                }
            }
            util->seqsLen[n.first->identifier] = newIdx;
            util->changeStorage(sIdx);
        }
        for (auto sIdx: tree->allNodes[n.second->identifier]->msaIdx) {
            int32_t rIdx = 0, qIdx = 0;
            int64_t start = sIdx*util->memLen;
            int newIdx = 0;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - util->seqsStorage[sIdx];
            while (rIdx < refLen && qIdx < qryLen) {
                if (allGapsR[rIdx] == false && allGapsQ[qIdx] == false) {
                    util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+qIdx];
                    qIdx++;rIdx++;newIdx++;
                }
                else if (allGapsR[rIdx] == true && allGapsQ[qIdx] == false) {
                    int consecGap = 0;
                    int k = rIdx;
                    while (allGapsR[k] && k < refLen) {
                        ++consecGap;
                        ++k;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = '-';
                        rIdx++;newIdx++;
                    }
                }
                else if (allGapsR[rIdx] == false && allGapsQ[qIdx] == true) {
                    int consecGap = 0;
                    int k = qIdx;
                    while (allGapsQ[k] && k < qryLen) {
                        ++consecGap;
                        ++k;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+qIdx];
                        qIdx++;newIdx++;
                    }
                }
                else {
                    int consecGap = 0;
                    int kr = rIdx, kq = qIdx;
                    while (allGapsR[kr] && kr < refLen) {
                        ++consecGap;
                        ++kr;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = '-';
                        rIdx++;newIdx++;
                    }
                    consecGap = 0;
                    while (allGapsQ[kq] && kq < qryLen) {
                        ++consecGap;
                        ++kq;
                    }
                    for (size_t g = 0; g < consecGap; ++g) {
                        util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+qIdx];
                        qIdx++;newIdx++;
                    }
                }
            }
            if (rIdx < refLen) {
                for (size_t g = rIdx; g < refLen; ++g) {
                    util->seqBuf[storeTo][start+newIdx] = '-';
                    newIdx++;          
                }
            }
            if (qIdx < qryLen) {
                for (size_t g = qIdx; g < qryLen; ++g) {
                    util->seqBuf[storeTo][start+newIdx] = util->seqBuf[storeFrom][start+g];
                    newIdx++;   
                }
            }
            util->seqsLen[n.first->identifier] = newIdx;
            util->changeStorage(sIdx);
        }

        printf("refLen:%d, qryLen:%d, alnLen: %d\n", refLen, qryLen, util->seqsLen[n.first->identifier]);
        
        
        tree->allNodes[n.first->identifier]->refStartPos = qryCut + refCut;
        std::vector<int> tempMSAIdx;
        for (int i = 0; i < qryCut; i++) tempMSAIdx.push_back(tree->allNodes[n.second->identifier]->msaIdx[i]);
        for (int i = 0; i < refNum; i++) tempMSAIdx.push_back(tree->allNodes[n.first->identifier]->msaIdx[i]);
        tree->allNodes[n.first->identifier]->msaIdx.clear();
        tree->allNodes[n.first->identifier]->msaIdx =tempMSAIdx;
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        std::cout << "MergeTime "<< mergeTime.count() / 1000000<< " ms\n";
        // std::cout << "Check (Length, SeqNum) = ("<< alignment[0].size() << ", " << alignment.size() << ')' << '\n';
    }
    
    // auto kernelEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds kernelTime = kernelEnd - kernelStart;
    // std::cout << "RunTime "<< kernelTime.count() / 1000<< " us\n";

    return;
}
*/

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& hier, std::stack<Node*> msaStack, int grpID, int mode) {
    int hierIdx = 0;
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
        std::vector<Node*> children;
        for (auto ch: node->children) {
            if (ch->grpID == grpID) children.push_back(ch);
        }
        if (children.empty()) {
            node->grpID = -2;
            msaStack.pop();
            continue;
        }
        else if (children.size() == 1 && node->parent != nullptr) {
            for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                if (node->parent->children[chIdx]->identifier == node->identifier) {
                    node->parent->children[chIdx] = children[0];
                }
            }
            msaStack.pop();
            continue;
        }
        size_t childIndex = 0;
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
    
    hierIdx = 0;
    std::stack<std::pair<Node*, int>> hierStack; 
    // Node* tempRoot = hier[0].first.first->parent;
    Node* preNode = hier[0].first.first;
    size_t prelevel = hier[0].first.first->level;
    hier[0].second = hierIdx;
    
    for (int k = 1; k < hier.size(); ++k) {
        // if (grpID == 1) std::cout << "Now: " << hier[k].first.first->identifier << '&' << hier[k].first.second->identifier << " at " << prelevel<<'\n';
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
                        // if (grpID == 1) {
                        //     std::cout << "TOP  "<< preNode->parent->identifier << ',' << hierIdx+1 << '\n';
                        // }
                        Node* tempNode = preNode;
                        while(true) {
                            Node* parent = tempNode->parent;
                            int childrenCount = 0;
                            for (auto ch: parent->children) {
                                if (ch->grpID == grpID) ++childrenCount;
                            }
                            if (childrenCount > 1) break;
                            tempNode = parent;
                        }
                        // hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                        hierStack.push(std::make_pair(tempNode->parent, (hierIdx+1)));
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
            // if (mode == 0) {
            //     if (hier[k].first.first->level <= prelevel) {
            //         hier[k].second = ++hierIdx;
            //         prelevel = hier[k].first.first->level;
            //     }
            //     else {
            //         hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
            //         hier[k].second = 0;
            //         hierIdx = 0;
            //         prelevel = hier[k].first.first->level;
            //     }
            // }
            // else {
                if (hier[k].first.first->level < prelevel || (hier[k].first.first->level == prelevel && hier[k].first.first->identifier == preNode->identifier)) {
                    hier[k].second = ++hierIdx;
                    prelevel = hier[k].first.first->level;
                }
                else {
                    // if (grpID == 1) {
                    //     std::cout << "TOP  "<< preNode->parent->identifier << ',' << hierIdx+1 << '\n';
                    // }
                    Node* tempNode = preNode;
                    while(true) {
                        Node* parent = tempNode->parent;
                        int childrenCount = 0;
                        for (auto ch: parent->children) {
                            if (ch->grpID == grpID) ++childrenCount;
                        }
                        if (childrenCount > 1) break;
                        tempNode = parent;
                    }
                    // hierStack.push(std::make_pair(preNode->parent, (hierIdx+1)));
                    hierStack.push(std::make_pair(tempNode->parent, (hierIdx+1)));
                    hier[k].second = 0;
                    hierIdx = 0;
                    prelevel = hier[k].first.first->level;
                }
            // }
        }
        preNode = hier[k].first.first;
    }
    // if (grpID == 1) {
    //     for (int h = 0; h < hier.size(); ++h) {
    //         std::cout << hier[h].first.first->identifier << ',' << hier[h].first.second->identifier << ',' << hier[h].second << '\n';
    //     } 
    // }
}


/*
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
}

*/



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

double getSPScore_gpu(std::vector<std::string>& alignment, msa::utility* util, Params& param) {
    auto scoreStart = std::chrono::high_resolution_clock::now();
    int numBlocks = 1024;
    int blockSize = 512;
    // size_t seqNum = util->memNum;
    // size_t seqLen = util->seqsLen["node_1"];
    size_t seqNum = alignment.size();
    size_t seqLen = alignment[0].size();
    
    printf("(Num, Len) - (%d, %d)\n", seqNum, seqLen);
    
    char*    hostSeqs = (char*)malloc(seqLen * seqNum * sizeof(char));
    int32_t* hostSeqInfo = (int32_t*)malloc(7 * sizeof(int32_t));
    int64_t* hostResult = (int64_t*)malloc(numBlocks * sizeof(int64_t));
    int seqCount = 0;
    // for (auto seq: util->seqs) {
    //     int sIdx = util->seqsIdx[seq.first]; 
    //     std::cout << sIdx << '\n';
    //     int storage = util->seqsStorage[sIdx];
    //     int start = sIdx*util->memLen;
    //     for (int j = 0; j < seqLen; ++j) {
    //         if (j > seqLen-100 && j < seqLen) std::cout << util->seqBuf[storage][start+j];
    //         hostSeqs[start+j] = util->seqBuf[storage][start+j];
    //         // if (j > seqLen-100 && j < seqLen) std::cout << hostSeqs[start+j];
    //     }
    //     std::cout << '\n';
    // }
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
    
    // printf("GPU score: %f, Runtime %d ms\n", score, scoreTime.count() / 1000000);
    
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




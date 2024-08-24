#ifndef UTIL_HPP
#include "../src/util.cuh"
#endif

po::options_description mainDesc("MSA Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,i", po::value<std::string>()->required(), "Input tip sequences - Fasta format (required)")
        ("cpu-num,c",  po::value<int>(), "Number of CPU threads")
        ("gpu-num,g",  po::value<int>(), "Number of GPUs")
        ("max-leaves,l",  po::value<int>(), "Maximum number of leaves per sub-subtree")
        ("max-subtree-size,m", po::value<int>()->default_value(1000000), "Maximum number of leaves per subtree")
        
        ("output,o", po::value<std::string>(), "Output file name")
        ("match",      po::value<paramType>()->default_value(2), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(0), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-8), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-1), "Gap extend penalty")
        ("trans",      po::value<paramType>(), "Transition score")
        ("xdrop",      po::value<paramType>(), "X-drop value")
        ("temp-dir", po::value<std::string>(), "Directory for storing temporary files")
        ("gpu-index", po::value<std::string>(), "Specify the GPU index, separated by commas. Ex. 0,2,3")
        ("read-batches", "Read sequences in batches and create temporary files")
        ("user-define-parameters", "Using user defined parameters. Please modify align.cuh")
        ("sum-of-pairs-score,s", "Calculate the sum-of-pairs score after the alignment")
        ("debug", "Enable debug on the final alignment")
        ("help,h", "Print help messages");

}


int main(int argc, char** argv) {

    auto mainStart = std::chrono::high_resolution_clock::now();
    
    parseArguments(argc, argv);
    po::variables_map vm;
    try{
        po::store(po::command_line_parser(argc, argv).options(mainDesc).run(), vm);
        po::notify(vm);
    }
    catch(std::exception &e){
        std::cerr << mainDesc << std::endl;
        if(vm.count("help"))
            return 0;
        else
            return 1;
    }

    msa::option* option = new msa::option;
    setOptions(vm, option);

    Params* param = setParameters(vm);

    // Partition tree into subtrees
    Tree* T = readNewick(vm);
    // printTree(T->root, -1);
    paritionInfo_t* P = new paritionInfo_t(option->maxSubtree, 0, 0, "centroid"); 
    partitionTree(T->root, P);
    Tree* newT = reconsturctTree(T->root, P->partitionsRoot);

    // int subtreeCount = 0;
    // for (auto p: P->partitionsRoot) {
    //     ++subtreeCount;
    //     std::cout << subtreeCount << '\t' << T->allNodes[p.first]->getNumLeaves() << '\n';
    // }
    // return;
    
    // Define MSA utility
    msa::utility* util = new msa::utility;

    
    std::unordered_map<std::string, std::string> beforeAln;
        
    // read sequences
    if (!option->readBatches) {
        readSequences(vm, util, T);
        if (T->m_numLeaves != util->seqsIdx.size()) {
            fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %lu)\n", T->m_numLeaves, util->seqsIdx.size()); 
            exit(1);
        }
        for (auto it = T->allNodes.begin(); it != T->allNodes.end(); ++it) {
            if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                printf("Missing Sequence %s.\n", it->first.c_str());
                exit(1);
            }
        }
        if (!option->debug) {
            for (auto s: T->allNodes) {
                if (s.second->is_leaf()) {
                    std::string seqName = s.first;
                    int sIdx = util->seqsIdx[s.first];
                    int storage = util->seqsStorage[sIdx];
                    std::string r = "";
                    int j = 0;
                    while (util->alnStorage[storage][sIdx][j] != 0) {
                        if (util->alnStorage[storage][sIdx][j] != '-') {
                            r += util->alnStorage[storage][sIdx][j];
                        }
                        ++j;
                    }
                    beforeAln[seqName] = r;
                }
            }   
        }
    }
    else {
        readSequencesNoutputTemp(vm, T, P);
    }
    // return;
    for (auto subRoot: P->partitionsRoot) {
        auto subtreeStart = std::chrono::high_resolution_clock::now();
        int subtree = T->allNodes[subRoot.first]->grpID;
        std::cout << "Start processing subtree No. " << subtree << '\n';
        
        Tree* subT;
        if (!option->readBatches) {
            subT = new Tree(subRoot.second.first);
            util->setSubtreeIdx(subtree);
        }
        else {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeTreeFile = tempDir + '/' + subtreeFileName + ".nwk";
            std::string subtreeSeqFile = tempDir + '/' + subtreeFileName + ".raw.fa";
            subT = readNewick(subtreeTreeFile);
            paritionInfo_t* tempP = new paritionInfo_t(INT_MAX, 0, 0, "centroid"); 
            partitionTree(subT->root, tempP);
            delete tempP;
            readSequences(subtreeSeqFile, util, subT);
            util->setSubtreeIdx(0);
            if (subT->m_numLeaves != util->seqsIdx.size()) {
                fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%lu != %lu)\n", subT->m_numLeaves, util->seqsIdx.size()); 
                exit(1);
            }
            for (auto it = subT->allNodes.begin(); it != subT->allNodes.end(); ++it) {
                if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                    printf("Missing Sequence %s.\n", it->first.c_str());
                    exit(1);
                }
            }
        }
        // printTree(subT->root, -1);
        std::cout << "Subtree No." << subtree << " contains "<< subT->m_numLeaves << " sequences.\n";
        // Partition subtree into sub-subtrees
        
        auto treeBuiltStart = std::chrono::high_resolution_clock::now();
        
        paritionInfo_t * subP = new paritionInfo_t(option->maxSubSubtree, 0, 0, "centroid");
        partitionTree(subT->root, subP);
       
        auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
        std::cout << "Partition the subtree into " << subP->partitionsRoot.size() << " trees in " <<  treeBuiltTime.count() / 1000000 << " ms\n";
        // Progressive alignment on each sub-subtree
        
        auto msaStart = std::chrono::high_resolution_clock::now();
        Tree* newSubT = reconsturctTree(subT->root, subP->partitionsRoot);
        // std::cout << newSubT->root->getNumLeaves() << '\n';
        msaOnSubtree(subT, util, option, subP, *param);
        auto msaEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds msaTime = msaEnd - msaStart;
        std::cout << "MSA on sub-subtree in " <<  msaTime.count() / 1000000000 << " s\n";
        
        if (subP->partitionsRoot.size() > 1) {
            // Align adjacent sub-subtrees to create overlap alignment
            auto alnStart = std::chrono::high_resolution_clock::now();
            alignSubtrees(subT, newSubT, util, option, *param);
            auto alnEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds alnTime = alnEnd - alnStart;
            std::cout << "Aligned adjacent sub-subtrees in " <<  alnTime.count() / 1000000 << " ms\n";

            auto mergeStart = std::chrono::high_resolution_clock::now();
            mergeSubtrees (subT, newSubT, util);
            auto mergeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
            int totalSeqs = subT->root->msaIdx.size();
            std::cout << "Merge " << newSubT->allNodes.size() << " sub-subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
            auto subtreeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds subtreeTime = subtreeEnd - subtreeStart;
            std::cout << "Finished the alignment on subtree No." << subtree << " in " << subtreeTime.count() / 1000000000 << " s\n";
        }
        
        for (auto sIdx: subT->root->msaIdx) T->allNodes[subT->root->identifier]->msaIdx.push_back(sIdx);
        // for (auto node: subT->allNodes) delete node.second;
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            std::string subtreeFileName = "subtree-" + std::to_string(subtree);
            std::string subtreeAlnFile = tempDir + '/' + subtreeFileName + ".temp.aln";
            std::string subtreeFreqFile = tempDir + '/' + subtreeFileName + ".freq.txt";
            outputFile(subtreeAlnFile, util, subT, 0);
            outputFreq(subtreeFreqFile, util, subT, subtree);
        }
        delete subT;
        delete subP;
        delete newSubT;
    }

    if (P->partitionsRoot.size() > 1) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        util->nowProcess = 2; // merge subtrees
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            readFreq(tempDir, T, P, util);
        }
        else {
            getFreq(T, P, util);
        }
        alignSubtrees(T, newT, util, option, *param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        std::cout << "Aligned adjacent subtrees in " <<  alnTime.count() / 1000000 << " ms\n";
        auto mergeStart = std::chrono::high_resolution_clock::now();
        mergeSubtrees (T, newT, util);
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        int totalSeqs = 0;
        if (option->readBatches) {
            std::string tempDir = (!vm.count("temp-dir")) ? "./temp" : vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            outputFinal (tempDir, T, P, util, totalSeqs);
        }
        else totalSeqs = T->root->msaIdx.size();
        std::cout << "Merge " << newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    }


    // post-alignment debugging
    if (option->debug) {
        auto dbgStart = std::chrono::high_resolution_clock::now();
        int alnLen = 0;
        // for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
        bool theFirst = true;
        for (auto s: beforeAln) {
            int sIdx = util->seqsIdx[s.first];
            int storage = util->seqsStorage[sIdx];
            std::string r = "";
            int offset = 0;
            while (util->alnStorage[storage][sIdx][offset] != 0) {
                if (util->alnStorage[storage][sIdx][offset] != '-') {
                    r += util->alnStorage[storage][sIdx][offset];
                }
                ++offset;
            }
            if (theFirst) {alnLen = offset; theFirst = false;}
            else {
                if (alnLen != offset) printf("seq: %s, the sequence length (%d) did not match (%d)\n", s.first.c_str(), offset, alnLen);
            }
            if (r != s.second) {
                printf("seq: %s, the sequence did not match\n", s.first.c_str());
                if (r.size() != s.second.size()) {
                    std::cout << "Wrong length. " << r.size() << '/' << s.second.size() << ".\n";
                }
                for (int i = 0; i < s.second.size(); ++i) {
                    if (r[i] != s.second[i]) {
                        std::cout << "Mismatch at position " << i << '\n';
                        break;
                    }
                }                
            }
            
        }
        auto dbgEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
        std::cout << "Completed checking " << T->m_numLeaves << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
    }
    
    // Calculate sum-of-pairs score
    if (vm.count("sum-of-pairs-score")) {
        auto spStart = std::chrono::high_resolution_clock::now();
        double score = getSPScore_gpu(util, *param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "Calculated Sum-of-Pairs-Score in " << spTime.count() / 1000000 << " ms. Score = " << score << ".\n";
    }
    
    // output MSA
    if (vm.count("output")) {
        std::string outFile = vm["output"].as<std::string>();
        if (outFile == "") outFile = "output.aln";
        auto outStart = std::chrono::high_resolution_clock::now();
        std::string subtreeFreqFile = outFile + ".freq.txt";
        outputFile(outFile, util, T, -1);
        outputFreq(subtreeFreqFile, util, T, -1);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    }
    // release memory
    util->seqsFree();
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return;
}

    // paritionInfo_t * partition = new paritionInfo_t(maxSubtreeSize, 0, 0, "centroid");
    
    // partitionTree(T->root, partition);
    
   
    
    // // int totalLeaves = 0;
    // // for (auto p: partition->partitionsRoot) {
    // //     printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     totalLeaves += p.second.second;
    // // }
    // // std::cout << "Total Leaves: " << totalLeaves << '\n';

    // Tree* newT = reconsturctTree(T->root, partition->partitionsRoot);    
    
    // // Start MSA on subtrees
    // auto msaStart = std::chrono::high_resolution_clock::now();
    // std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    // for (auto &p: partition->partitionsRoot) {
    //     std::stack<Node*> msaStack;
    //     getPostOrderList(p.second.first, msaStack);
        
    //     std::vector<std::pair<std::pair<Node*, Node*>, int>> subhier;
    //     int grpID = p.second.first->grpID;
    //     getMsaHierachy(subhier, msaStack, grpID, 0);
    //     for (auto h: subhier) {
    //         while (hier.size() < h.second+1) {
    //             std::vector<std::pair<Node*, Node*>> temp;
    //             hier.push_back(temp);
    //         }
    //         hier[h.second].push_back(h.first);
    //     }
    // }
    // int level = 0;
    
    // for (auto m: hier) {
    //     // std::cout << "Aln level: " << level << '\n';
    //     auto alnStart = std::chrono::high_resolution_clock::now();
    //     msaPostOrderTraversal_multigpu(T, m, util, param);
    //     auto alnEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    //     if (m.size() > 1) std::cout << "Level "<< level << ", " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
    //     else              std::cout << "Level "<< level << ", " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
    //     ++level;
    // }
    // auto msaEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds msaTime = msaEnd - msaStart;
    // std::cout << "MSA in " <<  msaTime.count() / 1000000000 << " s\n";

    // // Push MSA to partition roots
    // for (auto p: partition->partitionsRoot) {
    //     std::stack<Node*> msaStack;
    //     getPostOrderList(p.second.first, msaStack);
    //     std::vector<Node*> msaArray;
    //     while (!msaStack.empty()) {
    //         msaArray.push_back(msaStack.top());
    //         msaStack.pop();
    //     }
    //     if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
    //         if (msaArray.size() == 2) {
    //             T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
    //             util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
    //             break;
    //         }
    //         for (int m = msaArray.size()-2; m >=0; --m) {
    //             if (msaArray[m]->msaIdx.size()>0) {
    //                 T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
    //                 util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
    //                 break;
    //             }
    //         }
    //     }
    // }
    

    // // type 1 alignment
    // // int totalLeaves = 0;
    // // for (auto n: newT->allNodes) {
    // //     std::cout << n.first << ':' << T->allNodes[n.first]->msaIdx.size() << '\n';
    // //     totalLeaves += T->allNodes[n.first]->msaIdx.size();
    // //     // std::cout << '\n';
    // // }
    // // std::cout << "totalLeaves: " << totalLeaves << '\n';
    
    // // int totalLeaves = 0;
    // // for (auto p: partition->partitionsRoot) {
    // //     // std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\n';
    // //     // printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     if (T->allNodes[p.first]->msaIdx.size() != p.second.second) {
    // //         std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\t';
    // //         printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    // //     }
    // //     totalLeaves += T->allNodes[p.first]->msaIdx.size();
    // // }
    // // std::cout << "Total Leaves: " << totalLeaves << '\n';
    

    // // Create overlapped alignment (MSA between parent and children)
    // if (partition->partitionsRoot.size() > 1) {
    //     auto alnStart = std::chrono::high_resolution_clock::now();
    //     std::vector<std::pair<Node*, Node*>> type1Aln;
    //     // Copy from util->seqBuf to Node->msa
    //     for (auto n: newT->allNodes) {
    //         for (auto m: n.second->children) {
    //             if (newT->allNodes[m->identifier]->grpID == newT->allNodes[n.second->identifier]->grpID) {
    //                 type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
    //             }
    //         }
    //         T->allNodes[n.second->identifier]->msa.clear();
    //         for (int sIdx: T->allNodes[n.second->identifier]->msaIdx) {
    //             std::string r = "";
    //             int storage = util->seqsStorage[sIdx];
    //             int offset = 0;
    //             while (util->seqBuf[storage][sIdx][offset] != 0) {
    //                 r += util->seqBuf[storage][sIdx][offset];
    //                 ++offset;
    //             }
    //             T->allNodes[n.second->identifier]->msa.push_back(r);
    //         }
            
    //     }
    //     // for (auto n: type1Aln) {
    //     //     std::cout << n.first->identifier << '(' << T->allNodes[n.first->identifier]->msa.size() << ')'
    //     //               << n.second->identifier << '(' << T->allNodes[n.second->identifier]->msa.size() << ")\n";
    //     // }
    //     createOverlapMSA(T, type1Aln, util, param);
    //     auto alnEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds alnTime = alnEnd - alnStart;
    //     if (type1Aln.size() > 1) std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pairs in " <<  alnTime.count() / 1000000 << " ms\n";
    //     else                     std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pair in "  <<  alnTime.count() / 1000000 << " ms\n";
    //     hier.clear();

    //     // Create a new tree with all roots of subtrees
    //     auto mergeStart = std::chrono::high_resolution_clock::now();
    //     paritionInfo_t * newPartition = new paritionInfo_t(std::numeric_limits<size_t>::max(), 0, 0, "centroid");
    //     partitionTree(newT->root, newPartition); 
    //     // printTree(newT->root);
    //     // newPartition->partitionsRoot[midNode->identifier] = std::make_pair(newT->allNodes[midNode->identifier], newT->allNodes[midNode->identifier]->getNumNodes());
    //     // newPartition->partitionsRoot[newT->root->identifier] = std::make_pair(newT->root, newT->root->getNumNodes());
    //     std::vector<std::pair<Node*, Node*>> mergePairs;
    //     std::vector<std::pair<std::pair<Node*, Node*>, float>> sortedMergePairs;
    //     for (auto n: newT->allNodes) {
    //         if (n.second->children.size() > 1) {
    //             mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
    //             for (int i = 1; i < n.second->children.size(); ++i) {
    //                 mergePairs.push_back(std::make_pair(n.second->children[0], n.second->children[i]));
    //             }
    //         }
    //         else if (n.second->children.size() == 1) {
    //             mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
    //         }
    //     }
        
    //     std::vector<std::pair<Node*, Node*>> singleLevel;
    //     std::map<std::string, char> addedNodes;
    //     // std::cout << "\nDEBUG:MERGE PAIRS\n";
    //     // for (auto n: mergePairs) {
    //     //     std::cout << n.first->identifier << ',' << n.second->identifier << ':' << n.second->branchLength << '\n';
    //     // }
    //     while (true) {
    //         auto roundStart = std::chrono::high_resolution_clock::now();
    //         addedNodes.clear();
    //         singleLevel.clear();
    //         for (auto it = mergePairs.begin(); it != mergePairs.end();) {
    //             Node* a = it->first;
    //             Node* b = it->second;
    //             if ((a->parent != nullptr && b->parent != nullptr) &&  
    //                 (addedNodes.find(a->identifier) == addedNodes.end() && addedNodes.find(b->identifier) == addedNodes.end())) {
    //                 singleLevel.push_back(std::make_pair(a, b));
    //                 for (auto id: T->allNodes[a->identifier]->msa) addedNodes[id] = 0;
    //                 for (auto id: T->allNodes[b->identifier]->msa) addedNodes[id] = 0;
    //                 mergePairs.erase(it);
    //             }
    //             else {
    //                 ++it;
    //             }
    //         }
    //         bool breakLoop = false;
    //         if (singleLevel.empty()) {
    //             for (auto mp: mergePairs) {
    //                 singleLevel.push_back(mp);
    //                 // std::cout << mp.first->identifier << ':' << mp.second->identifier << '\n';
    //             }
    //             breakLoop = true;
    //         }
    //         transitivityMerge_cpu_mod(T, newT, singleLevel, util);
    //         // mergeLevels.push_back(singleLevel);
    //         auto roundEnd = std::chrono::high_resolution_clock::now();
    //         std::chrono::nanoseconds roundTime = roundEnd - roundStart;
    //         std::cout << "Merge "<< singleLevel.size() << " edges in " << roundTime.count() / 1000000 << " ms\n";
    //         if (breakLoop) break;
    //     }
    //     // transitivityMerge_cpu(T, mergePairs, util);
    //     // transitivityMerge_cpu_mod(T, newT, mergePairs, util);
    //     auto mergeEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
    //     int totalSeqs = 0;
    //     for (auto &p: newPartition->partitionsRoot) totalSeqs += T->allNodes[p.first]->msaIdx.size();
    //     // std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << T->root->msaIdx.size() << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
    //     std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";

    //     T->root->msa.clear();
    //     for (int sIdx: T->root->msaIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }
    // else {
    //     for (int sIdx: T->root->msaIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }   

    // }
    
    // bool regressive = false;
    // if (regressive) {
    //     int clusterSize = vm["max-leaves"].as<int>();
    //     if (clusterSize == 0) clusterSize = 1000;
    //     std::cout << "Cluster Size: " << clusterSize << '\n';
    //     getLongestDescendent(T, util);
    //     std::cout << "Finish longest descendent\n";
    //     std::map<int, std::vector<std::vector<Node*>>> subMSAs;
    //     getSubMSAs(subMSAs, T, clusterSize);
    //     std::cout << "Finish subMSAs, subMSA size: " << subMSAs.size() << "\n";
    //     for (auto level = subMSAs.begin(); level != subMSAs.end(); ++level) {
    //         std::cout << "Level: " << level->first << " Size: " << level->second.size() << '\n';
    //         std::vector<std::vector<std::pair<Node*, Node*>>> alnPairs;
    //         getAlnPairs(alnPairs, level->second);
    //         for (auto pairs: alnPairs) {
    //             msaPostOrderTraversal_multigpu_regressive(T, pairs, util, param);
    //         }
    //         std::vector<Node*> nodes;
    //         for (auto msa: level->second) {
    //             nodes.push_back(msa[0]);
    //         }
    //         storeMSA(T, nodes, util, level->first);
    //         resetSeqMem(nodes, util);
    //         for (auto n: T->allNodes) {
    //             n.second->msaIdx.clear();
    //             n.second->msa.clear();
    //         }
    //     }
    //     std::cout << "Finish MSA\n";
    //     transitivityMerge_regressive(T, subMSAs, util);
    //     std::cout << "Finish Merger\n";
    //     T->root->msa.clear();
    //     for (int sIdx = 0; sIdx < util->memNum; ++sIdx) {
    //         std::string r = "";
    //         int storage = util->seqsStorage[sIdx];
    //         int offset = 0;
    //         while (util->seqBuf[storage][sIdx][offset] != 0) {
    //             r += util->seqBuf[storage][sIdx][offset];
    //             ++offset;
    //         }
    //         T->root->msa.push_back(r);
    //     }
    // }

    


    // Calculate Sum-of-pairs score
    
    // auto spStart = std::chrono::high_resolution_clock::now();
    // double score;
    // if (calSP == "True" || calSP == "TRUE" || calSP == "true" || calSP == "T" || calSP == "t") {
    //     score = getSPScore_gpu(T->root->msa, util, param);
    //     // getSPScore_cpu(T->root->msa, param);
    //     auto spEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds spTime = spEnd - spStart;
    //     std::cout << "SP-score: " << score << ". Runtime " << spTime.count() / 1000000 << " ms\n";
    // }


    // // output file
    // std::string outFile = vm["output"].as<std::string>();
    // // if (outFile == "") outFile = "output.aln"; // default output ifle name
    // if (outFile != "") {
    //     auto outStart = std::chrono::high_resolution_clock::now();
    //     outputFile(outFile, util, T, -1);
    //     auto outEnd = std::chrono::high_resolution_clock::now();
    //     std::chrono::nanoseconds outTime = outEnd - outStart;
    //     std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    // }
    
    
    // util->seqFree();
    // auto mainEnd = std::chrono::high_resolution_clock::now();
    // std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    // std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    // return 0;
// }
#ifndef UTIL_HPP
#include "../src/util.cuh"
#endif

#include "../src/treePartition.cpp"
#include <tbb/task_scheduler_init.h>
#include <tbb/task_arena.h>


namespace po = boost::program_options;

KSEQ_INIT2(, gzFile, gzread)

po::options_description mainDesc("MSA Command Line Arguments");

void parseArguments(int argc, char** argv)
{
    // Setup boost::program_options
    mainDesc.add_options()
        ("tree,t", po::value<std::string>()->required(), "Initial Tree - Newick format (required)")
        ("sequences,i", po::value<std::string>()->required(), "Input tip sequences - Fasta format (required)")
        ("num-threads,c",  po::value<int>()->default_value(0), "Number of CPU threads")
        ("max-leaves,l",  po::value<int>()->default_value(0), "Maximum number of leaves per subtree")
        ("sofp-score,s",  po::value<std::string>()->default_value("false"), "Calculate the sum-of-pairs score (True, False, Only)")
        ("output,o", po::value<std::string>()->default_value(""), "Output file name")
        ("match",      po::value<paramType>()->default_value(2.5), "Match score")
        ("mismatch",   po::value<paramType>()->default_value(0), "Mismatch penalty")
        ("gap-open",   po::value<paramType>()->default_value(-4), "Gap open penalty")
        ("gap-extend", po::value<paramType>()->default_value(-1), "Gap extend penalty")
        ("trans",      po::value<paramType>()->default_value(0), "Transition score")
        ("xdrop",      po::value<paramType>()->default_value(0), "X-drop value")
        ("help,h", "Print help messages");

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

    int seqNum = 0, maxLen = 0;
    uint64_t totalLen = 0;
    
    std::map<std::string, std::string> seqs;
    
    while (kseq_read(kseq_rd) >= 0) {
        // ++seqNum;
        size_t seqLen = kseq_rd->seq.l;
        // seqs.push_back(std::make_pair(kseq_rd->name.s, std::string(kseq_rd->seq.s, seqLen)));
        seqs[kseq_rd->name.s] = std::string(kseq_rd->seq.s, seqLen);
        if (seqLen > maxLen) maxLen = seqLen;
        totalLen += seqLen;
    }
    
    seqNum = seqs.size();
    uint32_t avgLen = totalLen/seqNum;
    // printf("(Num, MaxLen, AvgLen) = (%d, %d, %f)\n", seqNum, maxLen, avgLen);
    std::cout << "(Num, MaxLen, AvgLen) = (" << seqNum << ", " << maxLen << ", " << avgLen << "\n";
    util->seqMalloc(seqNum, maxLen);

    
    // kseq_t* kseq_rd2 = kseq_init(f_rd);
    // uint32_t s = 0;
    // while (kseq_read(kseq_rd2) >= 0) {
    //     size_t seqLen = kseq_rd2->seq.l;
    //     for (int i = 0; i < util->memLen; ++i) {
    //         if (i < kseq_rd2->seq.l) util->seqBuf[0][s][i] = std::string(kseq_rd2->seq.s, seqLen)[i];
    //         else                     util->seqBuf[0][s][i] = 0; 
    //         util->seqBuf[1][s][i] = 0;  
    //     }
    //     // std::cout << seq.first << ':' << s << ':' << seq.second.size() << '\n';
    //     util->seqsIdx[kseq_rd2->name.s] = s;
    //     util->seqsLen[kseq_rd2->name.s] = seqLen;
    //     ++s;
    // }

    
    uint32_t s = 0;
    for (auto seq: seqs) {
        for (int i = 0; i < util->memLen; ++i) {
            if (i < seq.second.size()) util->seqBuf[0][s][i] = seq.second[i];
            else                       util->seqBuf[0][s][i] = 0; 
            util->seqBuf[1][s][i] = 0;  
        }
        // std::cout << seq.first << ':' << s << ':' << seq.second.size() << '\n';
        util->seqsIdx[seq.first] = s;
        util->seqsLen[seq.first] = seq.second.size();
        ++s;
        // if (s%1000 == 999) std::cout << s+1 << '\n';
    }

    
    
    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;
    std::cout << "Sequences read in: " <<  seqReadTime.count() / 1000000 << " ms\n";
}

bool cmp(std::string a, std::string b) {
    if (a.size() != b.size()) return a.size() < b.size();
    return a < b;
}

void outputFile(std::string fileName, msa::utility* util, Tree* T, int grpID) {
    std::ofstream outFile(fileName);
    std::vector<std::string> seqs;
    for (auto seq: util->seqsIdx) {
        if (grpID == -1 || T->allNodes[seq.first]->grpID == grpID) seqs.push_back(seq.first);
    }
    std::sort(seqs.begin(), seqs.end(), cmp);
    for (int s = 0; s < seqs.size(); ++s) {
        // outFile << '>' << seq.first << "\n";
        // int sIdx = util->seqsIdx[seq.first]; 
        outFile << '>' << seqs[s] << "\n";
        int sIdx = util->seqsIdx[seqs[s]]; 
        int storage = util->seqsStorage[sIdx];
        int i = 0;
        while (util->seqBuf[storage][sIdx][i] != 0) {
            outFile << util->seqBuf[storage][sIdx][i];
            ++i;
        }
        // std::cout << seq.first << ':' << i << '\n'; 
        outFile << '\n';
    }
    outFile.close();
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

    int readThreads = vm["num-threads"].as<int>();
    printf("Maximun available threads: %d\n", tbb::this_task_arena::max_concurrency());
    int numThreads = tbb::this_task_arena::max_concurrency();
    if (readThreads > 0 && readThreads < tbb::this_task_arena::max_concurrency()) {
        numThreads = readThreads;
    }
    printf("Using %d threads.\n", numThreads);
    tbb::task_scheduler_init init(numThreads);


    // Define MSA utility
    msa::utility* util = new msa::utility;

    // Read Tree (Newick String)
    // std::cout << "Start read tree...\n";
    Tree* T = readNewick(vm);
    
    // Read Input Sequences (Fasta format)
    readSequences(vm, util);

    if (T->m_numLeaves != util->seqsIdx.size()) {
        fprintf(stderr, "Error: Mismatch between the number of leaves and the number of sequences, (%d != %d)\n", T->m_numLeaves, util->seqsIdx.size()); 
        exit(1);
    }

    std::cout << "Node Num :" << T->allNodes.size() << " Size: " << sizeof(Node) <<  '\n';


    bool debug = false;
    // For debugging
    std::vector<std::string> beforeAln, afterAln;
    if (debug) {
        for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
            std::string r = "";
            int storage = util->seqsStorage[sIdx];
            int offset = 0;
            while (util->seqBuf[storage][sIdx][offset] != 0) {
                if (util->seqBuf[storage][sIdx][offset] != '-') {
                    r += util->seqBuf[storage][sIdx][offset];
                }
                ++offset;
            }
            beforeAln.push_back(r);
        }
        for (auto it = T->allNodes.begin(); it != T->allNodes.end(); ++it) {
            if (util->seqsIdx.find(it->first) == util->seqsIdx.end() && it->second->is_leaf()) {
                printf("Missing Sequence %s.\n", it->first.c_str());
                exit(1);
            }
        }
    }

    
    // printTree(T->root);
    paramType mat = vm["match"].as<paramType>();
    paramType mis = vm["mismatch"].as<paramType>();
    paramType gapOp = vm["gap-open"].as<paramType>();
    paramType gapEx = vm["gap-extend"].as<paramType>();
    paramType trans = vm["trans"].as<paramType>();
    paramType xdrop = vm["xdrop"].as<paramType>();
    // check parameters
    float gapOp_int =  round(gapOp);
    float gapEx_int =  round(gapEx);
    float xdrop_int =  round(xdrop);
    float gapOp_diff = gapOp-gapOp_int;
    float gapEx_diff = gapEx-gapEx_int;
    float xdrop_diff = xdrop-xdrop_int;
    if (gapOp_diff != 0) printf("WARNING: Floating point gap open penalty is not allowed, %.0f is used intead of %f.\n", gapOp_int, gapOp);
    if (gapEx_diff != 0) printf("WARNING: Floating point gap extend penalty is not allowed, %.0f is used intead of %f.\n", gapEx_int, gapEx);
    if (xdrop_diff != 0) printf("WARNING: Floating point xdrop is not allowed, %.0f is used intead of %f.\n", xdrop_int, xdrop);
    if (xdrop_int == 0) xdrop_int = round((FRONT_WAVE_LEN/3)*(-gapEx_int));
    if (trans == 0) trans = mis + (mat-mis)/2;
    int scoreMode = 0;
    Params param(mat,mis,trans,gapOp_int,gapEx_int,xdrop_int,scoreMode);

    std::string calSP = vm["sofp-score"].as<std::string>();
    
    // Calculate Sum-of-Pairs score and exit
    if (calSP == "Only" || calSP == "only") {
        auto spStart = std::chrono::high_resolution_clock::now();
        double score = 0;
        std::vector<std::string> alignment;
        for (int id = 0; id <  util->memNum; ++id) {
            std::string seq = "";
            int j = 0;
            while (util->seqBuf[0][id][j] != 0) {
                seq += util->seqBuf[0][id][j];
                ++j;
            }
            alignment.push_back(seq);
        }
        // for (auto seq: util->seqs) alignment.push_back(seq.second);
        bool same_length = true;
        int length = alignment[0].size(); 
        for (auto seq: alignment) {
            if (seq.size() != length) {
                same_length = false;
                break;
            }
        }
        if (!same_length) {
            fprintf(stderr, "Error: Non-uniform sequence lengths\n"); 
            exit(1);
        }
        score = getSPScore_gpu(alignment, util, param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "SP-score: " << score << ". Runtime " << spTime.count() / 1000000 << " ms\n";
        return;
    }


    bool progressive = true;
    if (progressive) {
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    int maxSubtreeSize = vm["max-leaves"].as<int>();
    if (maxSubtreeSize == 0) maxSubtreeSize = INT_MAX;
    paritionInfo_t * partition = new paritionInfo_t(maxSubtreeSize, 0, 0, "centroid");
    
    partitionTree(T->root, partition);
    
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cout << "Partition the tree in: " <<  treeBuiltTime.count() << " ns\n";
    
    // int totalLeaves = 0;
    // for (auto p: partition->partitionsRoot) {
    //     printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    //     totalLeaves += p.second.second;
    // }
    // std::cout << "Total Leaves: " << totalLeaves << '\n';

    Tree* newT = reconsturctTree(T->root, partition->partitionsRoot);    
    
    // Start MSA on subtrees
    auto msaStart = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::pair<Node*, Node*>>> hier;
    for (auto &p: partition->partitionsRoot) {
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
    int level = 0;
    int gpuNum;
    cudaGetDeviceCount(&gpuNum); // number of CUDA devices
    printf("Print GPU information\n");
    for(int i=0;i<gpuNum;i++) {
        // Query the device properties. 
        cudaDeviceProp prop; 
        cudaGetDeviceProperties(&prop, i);
        std::cout << "Device ID: "  << i << std::endl;
        std::cout << "Device Name: " << prop.name << std::endl; 
    }
    for (auto m: hier) {
        // std::cout << "Aln level: " << level << '\n';
        auto alnStart = std::chrono::high_resolution_clock::now();
        msaPostOrderTraversal_multigpu(T, m, util, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (m.size() > 1) std::cout << "Level "<< level << ", " << m.size() << " pairs in " <<  alnTime.count() / 1000000 << " ms\n";
        else              std::cout << "Level "<< level << ", " << m.size() << " pair in " <<  alnTime.count() / 1000000 << " ms\n";
        ++level;
    }
    auto msaEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds msaTime = msaEnd - msaStart;
    std::cout << "MSA in " <<  msaTime.count() / 1000000000 << " s\n";

    // Push MSA to partition roots
    for (auto p: partition->partitionsRoot) {
        std::stack<Node*> msaStack;
        getPostOrderList(p.second.first, msaStack);
        std::vector<Node*> msaArray;
        while (!msaStack.empty()) {
            msaArray.push_back(msaStack.top());
            msaStack.pop();
        }
        if (msaArray.back()->msaIdx.size() == 0 && msaArray.size() > 1) {
            if (msaArray.size() == 2) {
                T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[0]->msaIdx;
                util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[0]->identifier];
                break;
            }
            for (int m = msaArray.size()-2; m >=0; --m) {
                if (msaArray[m]->msaIdx.size()>0) {
                    T->allNodes[msaArray.back()->identifier]->msaIdx = msaArray[m]->msaIdx;
                    util->seqsLen[msaArray.back()->identifier] = util->seqsLen[msaArray[m]->identifier];
                    break;
                }
            }
        }
    }
    

    // type 1 alignment
    // int totalLeaves = 0;
    // for (auto n: newT->allNodes) {
    //     std::cout << n.first << ':' << T->allNodes[n.first]->msaIdx.size() << '\n';
    //     totalLeaves += T->allNodes[n.first]->msaIdx.size();
    //     // std::cout << '\n';
    // }
    // std::cout << "totalLeaves: " << totalLeaves << '\n';
    
    // int totalLeaves = 0;
    // for (auto p: partition->partitionsRoot) {
    //     // std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\n';
    //     // printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    //     if (T->allNodes[p.first]->msaIdx.size() != p.second.second) {
    //         std::cout << p.first << ':' << T->allNodes[p.first]->msaIdx.size() << '\t';
    //         printf("%s, leaves: %d, grpID: %d\n", p.first.c_str(), p.second.second, p.second.first->grpID);
    //     }
    //     totalLeaves += T->allNodes[p.first]->msaIdx.size();
    // }
    // std::cout << "Total Leaves: " << totalLeaves << '\n';
    

    // Create overlapped alignment (MSA between parent and children)
    if (partition->partitionsRoot.size() > 1) {
        auto alnStart = std::chrono::high_resolution_clock::now();
        std::vector<std::pair<Node*, Node*>> type1Aln;
        // Copy from util->seqBuf to Node->msa
        for (auto n: newT->allNodes) {
            for (auto m: n.second->children) {
                if (newT->allNodes[m->identifier]->grpID == newT->allNodes[n.second->identifier]->grpID) {
                    type1Aln.push_back(std::make_pair(T->allNodes[m->identifier], T->allNodes[n.second->identifier]));
                }
            }
            T->allNodes[n.second->identifier]->msa.clear();
            for (int sIdx: T->allNodes[n.second->identifier]->msaIdx) {
                std::string r = "";
                int storage = util->seqsStorage[sIdx];
                int offset = 0;
                while (util->seqBuf[storage][sIdx][offset] != 0) {
                    r += util->seqBuf[storage][sIdx][offset];
                    ++offset;
                }
                T->allNodes[n.second->identifier]->msa.push_back(r);
            }
            
        }
        // for (auto n: type1Aln) {
        //     std::cout << n.first->identifier << '(' << T->allNodes[n.first->identifier]->msa.size() << ')'
        //               << n.second->identifier << '(' << T->allNodes[n.second->identifier]->msa.size() << ")\n";
        // }
        createOverlapMSA(T, type1Aln, util, param);
        auto alnEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds alnTime = alnEnd - alnStart;
        if (type1Aln.size() > 1) std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pairs in " <<  alnTime.count() / 1000000 << " ms\n";
        else                     std::cout << "Create type 2 alignment "<<  type1Aln.size() <<" pair in "  <<  alnTime.count() / 1000000 << " ms\n";
        hier.clear();

        // Create a new tree with all roots of subtrees
        auto mergeStart = std::chrono::high_resolution_clock::now();
        paritionInfo_t * newPartition = new paritionInfo_t(std::numeric_limits<size_t>::max(), 0, 0, "centroid");
        partitionTree(newT->root, newPartition); 
        // printTree(newT->root);
        // newPartition->partitionsRoot[midNode->identifier] = std::make_pair(newT->allNodes[midNode->identifier], newT->allNodes[midNode->identifier]->getNumNodes());
        // newPartition->partitionsRoot[newT->root->identifier] = std::make_pair(newT->root, newT->root->getNumNodes());
        std::vector<std::pair<Node*, Node*>> mergePairs;
        std::vector<std::pair<std::pair<Node*, Node*>, float>> sortedMergePairs;
        for (auto n: newT->allNodes) {
            if (n.second->children.size() > 1) {
                mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
                for (int i = 1; i < n.second->children.size(); ++i) {
                    mergePairs.push_back(std::make_pair(n.second->children[0], n.second->children[i]));
                }
            }
            else if (n.second->children.size() == 1) {
                mergePairs.push_back(std::make_pair(n.second, n.second->children[0]));
            }
        }
        
        std::vector<std::pair<Node*, Node*>> singleLevel;
        std::map<std::string, char> addedNodes;
        // std::cout << "\nDEBUG:MERGE PAIRS\n";
        // for (auto n: mergePairs) {
        //     std::cout << n.first->identifier << ',' << n.second->identifier << ':' << n.second->branchLength << '\n';
        // }
        while (true) {
            auto roundStart = std::chrono::high_resolution_clock::now();
            addedNodes.clear();
            singleLevel.clear();
            for (auto it = mergePairs.begin(); it != mergePairs.end();) {
                Node* a = it->first;
                Node* b = it->second;
                if ((a->parent != nullptr && b->parent != nullptr) &&  
                    (addedNodes.find(a->identifier) == addedNodes.end() && addedNodes.find(b->identifier) == addedNodes.end())) {
                    singleLevel.push_back(std::make_pair(a, b));
                    for (auto id: T->allNodes[a->identifier]->msa) addedNodes[id] = 0;
                    for (auto id: T->allNodes[b->identifier]->msa) addedNodes[id] = 0;
                    mergePairs.erase(it);
                }
                else {
                    ++it;
                }
            }
            bool breakLoop = false;
            if (singleLevel.empty()) {
                for (auto mp: mergePairs) {
                    singleLevel.push_back(mp);
                    // std::cout << mp.first->identifier << ':' << mp.second->identifier << '\n';
                }
                breakLoop = true;
            }
            transitivityMerge_cpu_mod(T, newT, singleLevel, util);
            // mergeLevels.push_back(singleLevel);
            auto roundEnd = std::chrono::high_resolution_clock::now();
            std::chrono::nanoseconds roundTime = roundEnd - roundStart;
            std::cout << "Merge "<< singleLevel.size() << " edges in " << roundTime.count() / 1000000 << " ms\n";
            if (breakLoop) break;
        }
        // transitivityMerge_cpu(T, mergePairs, util);
        // transitivityMerge_cpu_mod(T, newT, mergePairs, util);
        auto mergeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds mergeTime = mergeEnd - mergeStart;
        int totalSeqs = 0;
        for (auto &p: newPartition->partitionsRoot) totalSeqs += T->allNodes[p.first]->msaIdx.size();
        // std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << T->root->msaIdx.size() << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";
        std::cout << "Merge "<< newT->allNodes.size() << " subtrees (total " << totalSeqs << " sequences) in " << mergeTime.count() / 1000000 << " ms\n";

        T->root->msa.clear();
        for (int sIdx: T->root->msaIdx) {
            std::string r = "";
            int storage = util->seqsStorage[sIdx];
            int offset = 0;
            while (util->seqBuf[storage][sIdx][offset] != 0) {
                r += util->seqBuf[storage][sIdx][offset];
                ++offset;
            }
            T->root->msa.push_back(r);
        }
    }
    else {
        for (int sIdx: T->root->msaIdx) {
            std::string r = "";
            int storage = util->seqsStorage[sIdx];
            int offset = 0;
            while (util->seqBuf[storage][sIdx][offset] != 0) {
                r += util->seqBuf[storage][sIdx][offset];
                ++offset;
            }
            T->root->msa.push_back(r);
        }
    }   

    }
    
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

    // debugging
    if (debug) {
        auto dbgStart = std::chrono::high_resolution_clock::now();
        int alnLen = 0;
        for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
            std::string r = "";
            int storage = util->seqsStorage[sIdx];
            int offset = 0;
            while (util->seqBuf[storage][sIdx][offset] != 0) {
                if (util->seqBuf[storage][sIdx][offset] != '-') {
                    r += util->seqBuf[storage][sIdx][offset];
                }
                ++offset;
            }
            if (sIdx == 0) alnLen = offset;
            else {
                if (alnLen != offset) printf("No: %d, the sequence length (%d) did not match (%d)\n", sIdx, offset, alnLen);
            }
            afterAln.push_back(r);
        }
        for (int sIdx = 0; sIdx < T->m_numLeaves; sIdx++) {
            if (beforeAln[sIdx] != afterAln[sIdx]) {
                printf("No: %d, the sequence did not match\n", sIdx);
            }
        }
        auto dbgEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
        std::cout << "Complete checking " << T->m_numLeaves << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
    }


    // Calculate Sum-of-pairs score
    
    auto spStart = std::chrono::high_resolution_clock::now();
    double score;
    if (calSP == "True" || calSP == "TRUE" || calSP == "true" || calSP == "T" || calSP == "t") {
        score = getSPScore_gpu(T->root->msa, util, param);
        // getSPScore_cpu(T->root->msa, param);
        auto spEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds spTime = spEnd - spStart;
        std::cout << "SP-score: " << score << ". Runtime " << spTime.count() / 1000000 << " ms\n";
    }


    // output file
    std::string outFile = vm["output"].as<std::string>();
    // if (outFile == "") outFile = "output.aln"; // default output ifle name
    if (outFile != "") {
        auto outStart = std::chrono::high_resolution_clock::now();
        outputFile(outFile, util, T, -1);
        auto outEnd = std::chrono::high_resolution_clock::now();
        std::chrono::nanoseconds outTime = outEnd - outStart;
        std::cout << "Output file in " <<  outTime.count() / 1000000 << " ms\n";
    }
    
    
    util->seqFree();
    auto mainEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds mainTime = mainEnd - mainStart;
    std::cout << "Total Execution in " << std::fixed << std::setprecision(6) << static_cast<float>(mainTime.count()) / 1000000000.0 << " s\n";
    return 0;
}
#ifndef MSA_HPP
#define MSA_HPP

#include "phylogeny.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>
#include <zlib.h>

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);

const float BLOSUM62[20][20] = {
    //    A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    {4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2},
    {0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2},
    {-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3},
    {-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2},
    {-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3},
    {0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3},
    {-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2},
    {-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1},
    {-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2},
    {-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1},
    {-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1},
    {-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2},
    {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3},
    {-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1},
    {-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2},
    {1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2},
    {0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2},
    {0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1},
    {-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2},
    {-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7}};

namespace msa
{
    using Node = phylogeny::Node;
    using Tree = phylogeny::Tree;
    using PartitionInfo = phylogeny::PartitionInfo;

    using NodePair = std::pair<phylogeny::Node *, phylogeny::Node *>;
    using NodePairVec = std::vector<NodePair>;

    // using gappyColumnQueue = std::pair<std::queue<std::tuple<int, int, float*>>, std::queue<std::tuple<int, int, float*>>>;
    // using gappyColumnList = std::vector<gappyColumnQueue>;
    using stringPair = std::pair<std::string, std::string>;
    using stringPairVec = std::vector<stringPair>;
    using IntPair = std::pair<int, int>;
    using IntPairVec = std::vector<IntPair>;
    using FloatPair = std::pair<float, float>;
    using alnPath = std::vector<int8_t>;
    using alnPathVec = std::vector<alnPath>;
    using Profile = std::vector<std::vector<float>>;

    struct Option
    {
        // Mode Options
        int alnMode; // 0: MSA from raw sequences, 1: merge multiple MSA files, 2: add new sequences to existing MSA
        // Hardware Options
        int gpuNum;
        int cpuNum;
        std::vector<int> gpuIdx;
        std::vector<size_t> gpuMem;
        bool cpuOnly;
        // Alignment Options
        int maxSubtree;
        float gappyVertical;
        float lenDev;
        float maxAmbig;
        bool debug;
        bool noFilter;
        bool reroot;
        bool compressed;
        char type; // 'n' for dna/rna, 'p' for protein
        bool alignGappy;
        // File Names
        std::string treeFile;
        std::string seqFile;
        std::string outFile;
        std::string backboneAlnFile;
        std::string msaDir;
        std::string tempDir;
        std::string subtreeDir;
        std::string merger;
        std::string outType;
        // Runtime Options
        bool printDetail;
        bool deleteTemp;
        // Constructor
        Option(po::variables_map &vm);
        // Only used in GPU version
        void getGpuInfo(po::variables_map &vm);
    };

    struct Params
    {
        float gapOpen;
        float gapExtend;   // for gap-affine
        float gapBoundary; // gap penalty at ends
        float xdrop;       // optional for now
        float scaleFactor;
        float **scoringMatrix;
        int matrixSize;
        Params(po::variables_map &vm, char type);
        ~Params();
    };

    struct SequenceDB
    {
        struct SequenceInfo
        {
            int id;
            std::string name;
            std::string unalignedSeq;
            int len;
            bool lowQuality;
            int subtreeIdx;
            float weight;
            // Storage
            bool storage; // where the sequence is stored, buffer 0 or 1
            bool backbone; // Used in PLACE_W_TREE
            int memLen;
            char *alnStorage[2];
            const float timesBigger = 2.0;
            // Functions
            void changeStorage();
            void memCheck(int len);
            void allocate_and_store(const std::string &seq);
            SequenceInfo(int id_, const std::string &name_, std::string &seq, int subtreeIdx_, float _weight, bool debug, int alnMode);
            ~SequenceInfo();
        };
        int currentTask = {0};
        std::vector<SequenceInfo*> sequences;
        std::vector<Node*> fallback_nodes;
        std::unordered_map<int, SequenceInfo *> id_map;
        std::unordered_map<std::string, SequenceInfo *> name_map;

        void addSequence(int id, const std::string &name, std::string &seq, int subtreeIdx, float weight, bool debug, int alnMode);
        void debug();
        void storeSubtreeProfile(Tree* subT, char type, int subtreeIdx);
        Tree* getPlacementTree(Tree* T);
        void cleanSubtreeDB();
        

        // For merging step
        std::unordered_map<int, alnPath> subtreeAln;
        std::vector<std::pair<std::string, int>> subAlnFiles; // subalignment file name, subtreeIdx
        
        SequenceDB() {};
        ~SequenceDB() {};

        

        // For debugging
        // std::unordered_map<std::string, std::string> rawSeqs;

        // std::unordered_map<std::string, int> seqsIdx;
        // std::unordered_map<int, std::string> seqsName;
        // std::unordered_map<std::string, int> seqsLen;
        // std::unordered_map<int, bool> seqsStorage;
        // std::map<int, std::vector<std::vector<float>>> profileFreq;
        // std::unordered_map<std::string, std::string> seqsCIGAR;
        // std::unordered_map<std::string, std::string> backboneAln;
        // std::unordered_map<std::string, std::string> placedSeqs;
        // std::unordered_map<int, bool> lowQuality;

        // char **alnStorage[2] = {nullptr, nullptr};
        // int nowProcess = 0;
        // int subtreeIdx = 0;
        // int memLen = 0;
        // int *seqMemLen = 0;
        // int memNum = 0;
        // int seqLen = 0;
        // int maxRawSeqLen = 0;
        // const float timesBigger = 2.0;

        // void setSubtreeIdx(int idx);
        // void seqsFree();
        // void seqFree(int i);
        // void seqMalloc(int seqNum, int seqLen, option* option);
        // void seqMalloc(int seqLen, int idx);
        // void seqsMallocNStore(size_t seqLen, std::map<std::string, std::pair<std::string, int>>& seqsMap, option* option);
        // void storeCIGAR();
        // void memCheck(int seqLen, option* option);
        // void memCheck(int seqLen, int idx);
        // void clearAll();
        // void debug(int& debugNum);
    };

    namespace io
    {
        void readSequenceNames(std::string seqFile, std::unordered_set<std::string> &seqNames);
        void readSequences(std::string fileName, SequenceDB *database, Option *option, Tree *&T);
        void readAlignment(std::string alnFileName, SequenceDB* database, Option* option, Node*& node);
        Tree* readAlignments_and_buildTree(SequenceDB *database, Option *option);
        void readBackboneAlignment(Tree* T, SequenceDB *database, Option *option);

        void writePrunedTree(Tree *T, msa::Option *option);
        void writeSubtrees(Tree *tree, PartitionInfo *partition, Option *option);
        void writeSubAlignments(SequenceDB* database, Option* option, int subtreeIdx, int alnLen);
        void writeFinalAlignments(SequenceDB* database, Option* option, int& totalSeqs);
        void writeAlignment(std::string fileName, SequenceDB* database, int alnLen, bool compressed);
        void writeAlignment(std::string fileName, stringPairVec& seqs, bool compressed, bool append);
        void update_and_writeAlignment(SequenceDB* database, Option* option, std::string fileName, int subtreeIdx);
        void writeWholeAlignment(SequenceDB* database, Option* option, int alnLen);
    }

    using alnFunction = std::function<void(Tree *, NodePairVec &, SequenceDB *, Option *, Params &)>;

    namespace alignment_helper
    {
        constexpr int _CAL_PROFILE_TH = 1000;

        void calculateProfile(float *profile, NodePair &nodes, SequenceDB *database, Option *option, int32_t memLen);
        void removeGappyColumns(float *hostFreq, NodePair &nodes, Option *option, std::pair<IntPairVec, IntPairVec> &gappyColumns, int32_t memLen, IntPair &lens, int currentTask);
        void calculatePSGP(float *hostFreq, float *hostGapOp, float *hostGapEx, NodePair &nodes, SequenceDB* database, Option *option, int memLen, IntPair offset, IntPair lens, Params &param);
        void getConsensus(Option *option, float *profile, std::string &consensus, int len);
        void pairwiseGlobal(const std::string &seq1, const std::string &seq2, alnPath &alnPath, Params &param);
        void addGappyColumnsBack(alnPath &aln_before, alnPath &aln_after, std::pair<IntPairVec, IntPairVec> &gappyColumns, Params &param, IntPair rgcLens, stringPair orgSeqs);
        void updateAlignment(NodePair &nodes, SequenceDB *database, alnPath &aln);
        void updateFrequency(NodePair &nodes, SequenceDB *database, alnPath &aln, FloatPair weights);
        void fallback2cpu(std::vector<int>& fallbackPairs,  NodePairVec& nodes, SequenceDB* database, Option* option);
    }

    namespace progressive
    {

        void collectPostOrder(Node *node, std::stack<Node *> &postStack);
        void getProgressivePairs(std::vector<std::pair<NodePair, int>> &alnOrder, std::stack<Node *> postStack, int grpID, int currentTask);
        void scheduling(Node *root, std::vector<NodePairVec> &alnPairsPerLevel, int currentTask);
        void updateNode(Tree *tree, NodePairVec &nodes, SequenceDB *database);
        void progressiveAlignment(Tree *T, SequenceDB *database, Option *option, std::vector<NodePairVec> &alnPairPerLevel, Params &param, alnFunction alignmentKernel);
        void msaOnSubtree(Tree *T, SequenceDB *database, Option *option, Params &param, alnFunction alignmentKernel);

        namespace cpu 
        {
            void allocateMemory_and_Initialize(float*& freq, float*& gapOp, float*& gapEx, int memLen, int profileSize);
            void freeMemory(float*& freq, float*& gapOp, float*& gapEx);
            void parallelAlignmentCPU(Tree *T, NodePairVec &alnPairs, SequenceDB *database, Option *option, Params &param);
            void alignmentKernel_CPU(Tree *T, NodePairVec &alnPairs, SequenceDB *database, Option *option, Params &param);
        }
        namespace gpu
        {
            // Constants
            constexpr int _MAX_LENGTH_N = (1 << 16);
            constexpr int _MAX_LENGTH_P = (1 << 14);
            // Pointers to GPU memory
            struct GPU_pointers {
                float *deviceFreq [8];
                float **devicePointers = nullptr;
                float *deviceGapOp = nullptr;
                float *deviceGapEx = nullptr;
                float *deviceParam = nullptr; // parameters
                int8_t  *deviceAln = nullptr;
                int32_t *deviceLen = nullptr;
                int32_t *deviceNum = nullptr;
                int32_t *deviceAlnLen = nullptr;
                int32_t *deviceSeqInfo = nullptr;
                // Pointers to Host memory
                float *hostFreq [8];
                float **hostPointers = nullptr;
                float *hostGapOp = nullptr;
                float *hostGapEx = nullptr; // gap extend
                float *hostParam = nullptr;
                int8_t  *hostAln = nullptr;
                int32_t *hostLen = nullptr;
                int32_t *hostNum = nullptr;
                int32_t *hostAlnLen = nullptr;
                int32_t *hostSeqInfo = nullptr;
                int gpu_index;
                int memBlock = {1};
                void allocateMemory(Option *option, int len, int numBlocks);
                void initializeHostMemory(Option *option, int len, int numBlocks, Params &param);
                void freeMemory();
                void memcpyHost2Device(Option *option, int len, int alnPairs, int numBlocks);
                void memcpyDevice2Host(Option *option, int len, int alnPairs);
                GPU_pointers(int index, int block){
                    this->gpu_index = index; 
                    this->memBlock = block;
                    for (int i = 0; i < 8; ++i) {
                        deviceFreq[i] = nullptr;
                        hostFreq[i] = nullptr;
                    }
                };
                ~GPU_pointers(){};
            };
            void alignmentKernel_GPU(Tree *T, NodePairVec &alnPairs, SequenceDB *database, Option *option, Params &param);
            void parallelAlignmentGPU(Tree *T, NodePairVec &alnPairs, SequenceDB *database, Option *option, Params &param);
        }
    }
}

#endif
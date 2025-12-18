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

namespace msa
{
    enum Type {
        DEFAULT_ALN = 0,
        MERGE_MSA   = 1,
        PLACE_WO_TREE = 2,
        PLACE_W_TREE  = 3
    };
    
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
        int maxLen;
        int minLen;
        bool writeFiltered;
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
            bool storage;  // where the sequence is stored, buffer 0 or 1
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
    };

    namespace io
    {
        void readSequenceNames(std::string seqFile, std::unordered_set<std::string> &seqNames);
        void readSequences(std::string fileName, SequenceDB *database, Option *option, Tree *&T, int subtree=-1);
        void readAlignment(std::string alnFileName, SequenceDB* database, Option* option, Node*& node);
        Tree* readAlignments_and_buildTree(SequenceDB *database, Option *option);
        void readBackboneAlignment(Tree* T, SequenceDB *database, Option *option);

        void writePrunedTree(Tree *T, msa::Option *option);
        void writeSubtrees(Tree *tree, PartitionInfo *partition, Option *option);
        void writeSubAlignments(SequenceDB* database, Option* option, int subtreeIdx, int alnLen);
        void writeAlignment(std::string fileName, SequenceDB* database, int alnLen, bool compressed);
        void writeAlignment(std::string fileName, stringPairVec& seqs, bool compressed, bool append);
        int update_and_writeAlignment(SequenceDB* database, Option* option, std::string fileName, int subtreeIdx);
        void update_and_writeAlignments(SequenceDB* database, Option* option, int& totalSeqs);
        void writeFinalMSA(SequenceDB* database, Option* option, int alnLen);
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
        void updateAlignment(NodePair &nodes, SequenceDB *database, Option *option, alnPath &aln);
        void updateFrequency(NodePair &nodes, SequenceDB *database, alnPath &aln, FloatPair weights);
        void fallback2cpu(std::vector<int>& fallbackPairs,  NodePairVec& nodes, SequenceDB* database, Option* option);
        void mergeInsertions(SequenceDB *database, Node* root);
        // void getStartingAndEndingPoints(std::string& ref, std::string& qry, IntPair& starting, IntPair& ending);
        // void alignEnds (NodePair &nodes, SequenceDB *database,  char type, IntPair& starting, IntPair& ending, std::map<int, std::string>& startAln, std::map<int, std::string>& endAln);
        // IntPair smith_waterman(const std::string &seq1, const std::string &seq2);
        // alnPath smith_waterman_tb(const std::string &seq1, const std::string &seq2);
    }

    namespace progressive
    {
        void getProgressivePairs(std::vector<std::pair<NodePair, int>> &alnOrder, std::stack<Node *> postStack, int grpID, int currentTask);
        void scheduling(Node *root, std::vector<NodePairVec> &alnPairsPerLevel, int currentTask);
        void updateNode(Tree *tree, NodePairVec &nodes, SequenceDB *database);
        void progressiveAlignment(Tree *T, SequenceDB *database, Option *option, std::vector<NodePairVec> &alnPairPerLevel, Params &param, alnFunction alignmentKernel);
        void msaOnSubtree(Tree *T, SequenceDB *database, Option *option, Params &param, alnFunction alignmentKernel, int subtree=-1);

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
                float  *deviceFreq [8];
                float  *deviceGapOp [8];
                float  *deviceGapEx [8];
                int8_t *deviceAln [8];

                float  **deviceFreqPointers = nullptr;
                float  **deviceGapOpPointers = nullptr; 
                float  **deviceGapExPointers = nullptr;
                int8_t **deviceAlnPointers = nullptr;
                
                float   *deviceParam = nullptr; // parameters
                int32_t *deviceLen = nullptr;
                int32_t *deviceNum = nullptr;
                int32_t *deviceAlnLen = nullptr;
                int32_t *deviceSeqInfo = nullptr;
                // Pointers to Host memory
                float  *hostFreq [8];
                float  *hostGapOp [8];
                float  *hostGapEx [8];
                int8_t *hostAln [8];

                float  **hostFreqPointers = nullptr;
                float  **hostGapOpPointers = nullptr; 
                float  **hostGapExPointers = nullptr;
                int8_t **hostAlnPointers = nullptr;
                
                float *hostParam = nullptr;
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
                void memcpyDevice2Host(Option *option, int len, int alnPairs, int numBlocks);
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
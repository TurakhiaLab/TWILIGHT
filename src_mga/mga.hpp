#ifndef MGA_HPP
#define MGA_HPP

#include <iostream>
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

#include "phylogeny.hpp"
#include "option.hpp"
#include "block.hpp"

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);



namespace mga
{
    using IntVec = std::vector<int>;
    using IntPair = std::pair<int, int>;
    using CigarOp = std::pair<int, char>;
    using Cigar = std::vector<CigarOp>;

    using Tree = phylogeny::Tree;
    using Node = phylogeny::Node;
    using NodePair = std::pair<Node*, Node*>;
    using NodePairVec = std::vector<NodePair>;
    
    struct Alignment {
        int identifier;
        IntPair refIdx;  // first: start, second: end (not included)
        IntPair qryIdx;
        bool inverse;
        Cigar CIGAR;
        float alnScore;

        // For chaining (Identifying primary alignments)
        bool used;
        float chainScore;
        

        // Primary or Secondary Alignment
        bool primary;  // True: primary, False: Secondary -> default: False
        std::unordered_set<int> duplications;  // Store secondary alignments that are duplications
        // For spliting
        bool valid;
        // Function
        void show();
        int mapCoordinate(int pos, bool inputIsRef);

        // Constructor
        Alignment(): used(false), chainScore(0), primary(false), valid(true) {};
        
    };

    struct AlnChain {
        int identifier;
        IntVec chainedAln;
        float score;
        AlnChain(int id, float score, IntVec& chain);
    };

    // --- Helper: Interval Management ---
    struct IntervalIndex {
        // Stores sorted, non-overlapping intervals [start, end)
        std::vector<IntPair> intervals;
        // Adds an interval and merges overlaps (assuming the input doesn't overlap existing significantly, 
        // strictly used for tracking occupied space after we decided to keep an alignment)
        void add(int start, int end);
        // Check overlap with existing intervals
        // Returns: Total overlap amount
        // Outputs: requiredTrimFront (amount to cut from start), requiredTrimBack (amount to cut from end)
        int checkOverlap(int start, int end, int& trimFront, int& trimBack);
    };

    using alnVec = std::vector<Alignment>;
    using chainVec = std::vector<AlnChain>;


    namespace parser {
        Cigar parseCigar(const std::string& cigar);
        // alnVec parseSAM (const std::string& filename);
        // alnVec parsePAF (const std::string& filename);
        alnVec parseMinimap2(const std::string& filename, bool isSAM);
    }

    struct SplitResult {
        Alignment head;   
        Alignment middle; 
        Alignment tail;
    };

    struct OverlapInfo {
        int primaryIdx;
        bool onRef; // true = ref, false = qry
        int overlapLen;
        int start;  // absolute coordinate of overlap start
        int end;    // absolute coordinate of overlap end
    };

    // Helper struct for cross-axis processing
    struct CrossOverlap {
        int pIdx;         // The index of the 'Other' Primary (PA, PB)
        int sStart, sEnd; // The overlap range on S's CURRENT axis (The "Other" axis)
        int pStart, pEnd; // The overlap range mapped to P_Main's axis (The "Best" axis)
        bool onRef;       // The axis of pStart/pEnd (same as BestOverlap)
    };

    chainVec getAlignmentChains(alnVec& alignments);
    void identifyPrimaryAlignments(alnVec& alignments, chainVec& chains);
    void trimHead(Alignment& aln, int targetRefTrim, int targetQryTrim);
    void trimTail(Alignment& aln, int targetRefTrim, int targetQryTrim);
    std::pair<Alignment, Alignment> singleSplit(const Alignment& parent, int splitPos, bool onRef);
    SplitResult splitAlignment(const Alignment& aln, int start, int end, bool onRef);
    void detectDuplications(alnVec& alignments);


    namespace io
    {
        BlockManager& readSequences(std::string& fileName, Option& option, Tree& T);
    }

    namespace progressive {
        void getProgressivePairs(std::vector<std::pair<NodePair, int>> &alnOrder, std::stack<Node *> postStack, int grpID, int currentTask);
        void scheduling(Node* root, std::vector<NodePairVec>& levels, int mode);
        void updateNode(NodePairVec& nodes, BlockManager& blockManager);
        void progressiveAlignment(Tree& T, Option& option, std::vector<NodePairVec>& alnPairsPerLevel, BlockManager& blockManager);
        void msaOnSubtree(Tree& T, Option& option, BlockManager& blockManager, int subtree);
        void alignmentKernel(NodePairVec& alnPairs, BlockManager& blockManager, Option& option);
    }
}

#endif
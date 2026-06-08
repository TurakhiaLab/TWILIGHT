#ifndef MGA_HPP
#define MGA_HPP


#include <boost/program_options.hpp>

#include "type.hpp"
#include "block.hpp"
#include "phylogeny.hpp"
#include "option.hpp"

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);

namespace mga {

    /*
    void identifyPrimaryAlignments(alnVec& alignments, BlockSet* ref_blockset, BlockSet* qry_blockset, stringMap& ref_seqs, stringMap& qry_seqs);
    void recoverSuperCoordinates(alnVec& alignments, BlockSet* ref_blockset, BlockSet* qry_blockset);

    void resolveGlobalOverlaps(alnVec& alignments);
    std::pair<Alignment, Alignment> singleSplit(const Alignment& parent, int splitPos, bool onRef);
    //vSplitResult splitAlignment(const Alignment& aln, int start, int end, bool onRef);
    // void detectDuplications(alnVec& alignments);
    void fillUnalignedRegions(alnVec& alignments, int refTotalLen, int qryTotalLen);
    bool validateCoverage(const alnVec& alignments, int refTotalLen, int qryTotalLen);
    void collectCutPoints(const std::vector<mga::Alignment>& alignments, std::set<int>& refCuts, std::set<int>& qryCuts);

    std::vector<Alignment> splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts, const std::set<int>& qryCuts);
    void snapAlignment(mga::Alignment& aln, int r_pad_left, int r_pad_right, int q_pad_left, int q_pad_right);
    std::vector<Alignment> splitAlignmentsByCuts(const std::vector<Alignment>& alignments, const std::set<int>& refCuts, const std::set<int>& qryCuts);
    */

    namespace io
    {
        std::unique_ptr<BlockManager> readSequences(std::string& fileName, Option& option, Tree& tree);
        void writeAlignment(std::string fileName, StringPairs& seqs, bool compressed, bool append);
        void writeMAF(BlockSet* blockSet, const std::string& outputFileName);
    }

    namespace progressive {
        void getProgressivePairs(std::vector<std::pair<NodePair, int>> &alnOrder, std::stack<Node *>& postStack, int grpID, int currentTask);
        void scheduling(Node* root, std::vector<NodePairs>& levels, int mode);
        void updateNode(NodePairs& nodes, BlockManager* blockManager);
        void progressiveAlignment(Tree& T, Option& option, std::vector<NodePairs>& alnPairsPerLevel, BlockManager* blockManager);
        void msaOnSubtree(Tree& T, Option& option, BlockManager* blockManager, int subtree);
        void alignmentKernel(NodePairs& alnPairs, BlockManager* blockManager, Option& option, Tree& tree);
    }

    // std::vector<std::pair<int, char>> miniGlobalAlignment(const std::string& ref, const std::string& qry);
    // int getNearestBoundary(int pos, const std::set<int>& bnds, int threshold);
    // void snapAlignmentsToBlockBoundaries(alnVec& alignments, BlockSet* bs1, BlockSet* bs2, const std::string& refSeq, const std::string& qrySeq, int threshold = 100);
}


// mga::Cigar adjustCigarWithVariations(mga::Cigar& origCigar, Segment& refSeg, Segment& qrySeg, bool qryInverse, int qryConsLen);
// mga::Cigar extractSubCigar(const mga::Cigar& origCigar, int refOffset, int refLen);

#endif // MGA_HPP

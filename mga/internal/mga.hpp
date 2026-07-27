#ifndef MGA_HPP
#define MGA_HPP


#include <boost/program_options.hpp>

#include "type.hpp"
#include "block_manager.hpp"
#include "phylogeny.hpp"
#include "option.hpp"

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);

namespace mga {

    namespace io
    {
        std::unique_ptr<BlockManager> readSequences(std::string& fileName, Option& option, Tree& tree);
        void writeAlignment(std::string fileName, StringPairs& seqs, bool compressed, bool append);
        void writeMAF(BlockSet* blockSet, const std::string& outputFileName, bool grouped=false);
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

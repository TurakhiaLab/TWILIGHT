#ifndef PROGRESSIVE_HPP
#define PROGRESSIVE_HPP

#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <stack>
#include <vector>

namespace mga {
namespace progressive {

    void getProgressivePairs(std::vector<std::pair<NodePair, int>>& alnOrder, std::stack<Node*>& postStack, int grpID, int mode);
    
    void scheduling(Node* root, std::vector<NodePairVec>& levels, int mode);
    
    void updateNode(NodePairVec& nodes, BlockManager& blockManager);
    
    void progressiveAlignment(Tree& T, Option& option, std::vector<NodePairVec>& alnPairsPerLevel, BlockManager& blockManager);
    
    void msaOnSubtree(Tree& T, Option& option, BlockManager& blockManager, int subtree = -1);

} // namespace progressive
} // namespace mga

#endif // PROGRESSIVE_HPP
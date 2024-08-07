#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <stack> 
#include <utility> 
#include <limits> 

#ifndef TREE_HPP
#include "tree.hpp"
#endif




class paritionInfo_t
{
public:
    size_t maxPartitionSize;
    size_t minPartitionSize;
    size_t numPartitions;
    std::stack<Node*> partitionsStack;
    std::unordered_map<std::string, std::pair<Node*, size_t>> partitionsRoot;
    std::string partitionOption;

    paritionInfo_t (size_t t_maxPartitionSize, size_t t_minPartitionSize, size_t t_numPartitions, std::string option):
        maxPartitionSize(t_maxPartitionSize), minPartitionSize(t_minPartitionSize), numPartitions(t_numPartitions) {
            if (option == "centroid" || option == "longest") partitionOption = option;
            else {fprintf(stderr, "Error: Wrong partition option: %s\n", option.c_str()); exit(1);} 
        }

    void createPartition(Tree* T, paritionInfo_t* partition);
    
};

void findLongestBranchLength(Node* node, int grpID, float& maxBranchLength, Node ** nodeWithLongestBranchLength);

size_t getNumLeaves(Node* node, int grpID);
bool nodeWithNoLeafTraversal(Node* node);
void internalNodeResolution(Node* node, paritionInfo_t* partition);
void updateCentroidEdge (Node* node, Node* root, size_t halfTaxa, size_t& Imbalance, Node*& breakEdge);
Node* getCentroidEdge(Node* node, int rootID);
void updateLongestEdge (Node* node, Node* root, size_t totalLeaves, size_t minSize, float& longestLen, Node*& breakEdge);
Node* getLowestChild(Node* node, std::string identifier);
Node* getLongestEdge(Node* node, int rootID, int minSize);
Node* getBreakingEdge(Node* root, std::string option, int minSize);
void setChildrenGrpID(Node*& node, int ID_org, int ID);
void bipartition(Node* root, Node* edge, Node*& tree1Root, Node*& tree2Root, paritionInfo_t*& partition);
void partitionTree(Node* root, paritionInfo_t* partition); 
void preOrderTraversal(Node* parent, Node* node, Tree*& T, std::unordered_map<std::string, std::pair<Node*, size_t>>& nodes);
Tree* reconsturctTree(Node* root, std::unordered_map<std::string, std::pair<Node*, size_t>>& partitionRoots); 


void preOrderTraversal(Node* parent, Node* node, Tree*& T, std::map<std::string, std::string>& nodes);
Tree* reconsturctTree(Node* root, std::map<std::string, std::string>& partitionRoots); 

#endif

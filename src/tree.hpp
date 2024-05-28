#ifndef TREE_HPP
#define TREE_HPP

#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <iostream>
#include <chrono>


class Node 
{
public:
    Node(std::string id, float len);
    Node(std::string id, Node* par, float len);
    size_t getNumLeaves();
    bool is_leaf() {return !(identifier.substr(0,4) == "node");}
    void setNumleaves() {numLeaves = getNumLeaves();};

    float branchLength;
    size_t level;


    std::string identifier;
    Node* parent;
    std::vector< Node* > children;
    std::vector<std::string> msa; //use this to store identifier
    std::vector<int> msaIdx;
    std::vector<int8_t> msaAln;
    std::vector<int> addGapPos;
    // std::vector<std::vector<uint16_t>> freq;
    size_t refStartPos = {0};
    size_t numLeaves = {0};

    /*Partition*/
    int grpID;
    Node * partitionParent;
};

class Tree
{
public:

    size_t m_currInternalNode{ 0 };
    size_t m_maxDepth{ 0 };
    size_t m_numLeaves{ 0 };
    float m_meanDepth{ 0 };
    std::string newInternalNodeId() { return "node_" + std::to_string(++m_currInternalNode);}

    Node* root;
    std::unordered_map< std::string, Node* > allNodes;
    Tree(std::string newick);
    Tree() {root = nullptr;}
};

#endif
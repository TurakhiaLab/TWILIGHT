#ifndef TREE_HPP
#define TREE_HPP

#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <map>


class Node 
{
public:
    Node(std::string id, float len);
    Node(std::string id, Node* par, float len);
    size_t getNumLeaves();
    size_t getNumNodes();
    bool is_leaf() {return !(identifier.substr(0,4) == "node");}
    void setNumleaves() {numLeaves = getNumLeaves();};
    void setLongestDescendant(Node* n) {longestDescendant = n;};

    float branchLength;
    size_t level;
    std::string identifier;
    Node* parent;
    Node* longestDescendant;
    std::vector< Node* > children;
    std::vector<std::string> msa; //use this to store identifier
    std::vector<int> msaIdx;
    std::vector<int8_t> msaAln;
    std::vector<std::vector<float>> msaFreq;
    size_t numLeaves = {0};
    float weight = {0};

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
    void calLeafNum();
    void calSeqWeight();
    Tree(std::string newick);
    Tree(Node* node);
    Tree() {root = nullptr;}
    ~Tree();
};

#endif
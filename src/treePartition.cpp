#ifndef TREE_HPP
#include "tree.hpp"
#endif

#include <stack> 
#include <utility> 
#include <limits> 


class paritionInfo_t
{
public:
    size_t maxPartitionSize;
    size_t numPartitions;
    std::stack<Node*> partitionsStack;
    std::unordered_map<std::string, std::pair<Node*, size_t>> partitionsRoot;

    paritionInfo_t (size_t t_maxPartitionSize, size_t t_numPartitions):
        maxPartitionSize(t_maxPartitionSize), numPartitions(t_numPartitions) {}

    void createPartition(Tree* T, paritionInfo_t* partition);
};

void findLongestBranchLength(Node* node, int grpID, float& maxBranchLength, Node ** nodeWithLongestBranchLength)
{
    if (node->children.size() == 0) return;
    if (!(node->grpID == -1 || node->grpID == grpID)) return;
    
    for (auto& child: node->children)
    {
        findLongestBranchLength(child, grpID, maxBranchLength, nodeWithLongestBranchLength);
    }
    if (node->branchLength > maxBranchLength)
    {
        *nodeWithLongestBranchLength = node;
        maxBranchLength = node->branchLength;
    }

}

size_t getNumLeaves(Node* node, int grpID)
{
    size_t children = 0;
    // std::cout << node->identifier << "\t" << node->grpID << "\t" << grpID << "\n";

    if (!(node->grpID==grpID || node->grpID ==-1)) return 0;

    if (node->children.size()==0) return 1;

    for (auto &child: node->children)
        children += getNumLeaves(child, grpID);

    // std::cout << node->identifier << "\t" << children << "\n"; 
    return children;
}

bool nodeWithNoLeafTraversal(Node* node)
{
    if (node->children.size() == 0) return false;

    if (node->grpID != -1) return true;

    bool nodeWithNoLeaf = true;
    for (auto &child: node->children)
    {
        if (!nodeWithNoLeafTraversal(child)) 
        {
            nodeWithNoLeaf = false;
            break;
        }
    }

    return nodeWithNoLeaf;
}

void internalNodeResolution(Node* node, paritionInfo_t* partition)
{
    if (node->grpID != -1) return;

    if (!nodeWithNoLeafTraversal(node)) return;

    Node* newNode = nullptr; 
    float bl = std::numeric_limits<float>::infinity();
    for (auto &child: node->children)
    {
        if (child->grpID != 1)
        {
            if (child->branchLength < bl)
            {
                bl = child->branchLength;
                newNode = child;
            }
        }
    }

    if (newNode != nullptr)
    {
        partition->partitionsRoot[node->identifier] = std::make_pair(node, partition->partitionsRoot[newNode->identifier].second);
        partition->partitionsRoot.erase(partition->partitionsRoot.find(newNode->identifier));
    }

    return;

}

void paritionInfo_t::createPartition(Tree* T, paritionInfo_t* partition)
{
    Node* root = T->root;
    partition->partitionsRoot[root->identifier] = std::make_pair(root, T->m_numLeaves);
    root->grpID = partition->numPartitions++;
    partitionsStack.push(root);


    while (!partitionsStack.empty())
    {
        Node* node = partitionsStack.top();
        std::cout << node->identifier << "\t" << node->grpID << "\t" <<  partition->partitionsRoot[node->identifier].second << "\t";
        if (partition->partitionsRoot[node->identifier].second < partition->maxPartitionSize) 
        {
            partition->partitionsStack.pop();
            std::cout << "\n";
            continue;
        }
        
        int grpID = node->grpID;
        float maxBranchLength = 0.0;
        Node* nodeWithLongestBranchLength = nullptr;
        findLongestBranchLength(node, grpID, maxBranchLength, &nodeWithLongestBranchLength);
        
        nodeWithLongestBranchLength->grpID = partition->numPartitions++;
        size_t numChildNodes = getNumLeaves(nodeWithLongestBranchLength, nodeWithLongestBranchLength->grpID);
        
        std::cout << nodeWithLongestBranchLength->identifier << "\t" << nodeWithLongestBranchLength->grpID << "\t" << numChildNodes << "\t";
        
        partition->partitionsRoot[nodeWithLongestBranchLength->identifier] = std::make_pair(nodeWithLongestBranchLength ,numChildNodes);
        partition->partitionsRoot[node->identifier].second -= numChildNodes;
        if (numChildNodes > partition->maxPartitionSize)
            partition->partitionsStack.push(nodeWithLongestBranchLength);
        if (partition->partitionsRoot[node->identifier].second < partition->maxPartitionSize) partitionsStack.pop();

        std::cout << node->identifier << "\t" << node->grpID << "\t" << partition->partitionsRoot[node->identifier].second << "\n";

    }

    /* To resolve internal node partition */
    for (auto &p: partition->partitionsRoot)
    {
        internalNodeResolution(p.second.first, partition);
    }
    
    /* Assigning partition parent node */
    for (auto &p: partition->partitionsRoot)
    {
        if (p.second.first->identifier == root->identifier) {p.second.first->partitionParent = nullptr; continue;}
        
        Node* parentNode = p.second.first->parent;
        while (true)
        {
            if (parentNode->identifier == root->identifier) {p.second.first->partitionParent = root; break;}
            else if (parentNode->grpID != -1) {p.second.first->partitionParent = parentNode;break;}
            else parentNode = parentNode->parent;
        }
    }
    

    
    for (auto &p: partition->partitionsRoot) std::cout << p.first << " " << p.second.second << "\t" << (p.second.first->partitionParent==nullptr? "NullPtr":p.second.first->partitionParent->identifier) << "\n";
    
    return;
}



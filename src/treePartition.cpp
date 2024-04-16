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
    std::string partitionOption;

    paritionInfo_t (size_t t_maxPartitionSize, size_t t_numPartitions, std::string option):
        maxPartitionSize(t_maxPartitionSize), numPartitions(t_numPartitions) {
            if (option == "centroid" || option == "longest") partitionOption = option;
            else {fprintf(stderr, "Error: Wrong partition option: %s\n", option.c_str()); exit(1);} 
        }

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

    // if (!(node->grpID==grpID || node->grpID ==-1)) return 0; /* For Sumit's partition */
    if (node->grpID != grpID) return 0;

    if (node->children.size() == 0) return 1;

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

void updateCentroidEdge (Node* node, Node* root, size_t halfTaxa, size_t& Imbalance, Node*& breakEdge){
    if (node->grpID != root->grpID || node->children.size() == 0) return;
    
    for (auto ch: node->children){
        updateCentroidEdge(ch, root, halfTaxa, Imbalance, breakEdge);
    }
    
    // if (node->parent->identifier == root->identifier) return;
    size_t numDescendants = getNumLeaves(node, root->grpID);
    size_t tempImbalance = (halfTaxa > numDescendants) ? (halfTaxa - numDescendants) : (numDescendants - halfTaxa);
    // std::cout << "DEBUG: updateCentroid " << tempImbalance << '\t' << Imbalance << '\n';
    // std::cout << "DEBUG: " << node->identifier << '\t' << breakEdge->identifier << '\n' ; 
    if (tempImbalance < Imbalance) {
        breakEdge = node;
        Imbalance = tempImbalance;
    }

    return;
}

Node* getCentroidEdge(Node* node, int rootID){
    Node* centroidEdge = node;
    size_t numLeaves = getNumLeaves(node, node->grpID);
    // std::cout << "DEBUG: getCentroid " << node->identifier << '\t' << numLeaves << '\n';
    size_t centroidImbalance = numLeaves;
    size_t halfTaxa = numLeaves/2;
    halfTaxa = (halfTaxa == 0) ? 1 : halfTaxa;
    // size_t numDescendants = 0;
    // for (auto ch: node->children) {
        updateCentroidEdge (node, node, halfTaxa, centroidImbalance, centroidEdge);
    // }
    return centroidEdge;
}

void updateLongestEdge (Node* node, Node* root, size_t totalLeaves, float& longestLen, Node*& breakEdge){
    // std::cout << "DEBUG: updateLongestEdge " << node->grpID << '\t' << rootID << '\n';
    if (node->grpID != root->grpID ) return;
    size_t numLeaves = getNumLeaves(node, node->grpID);
    if ((node->children).size() != 0 && numLeaves == 0) return;
    if ((node->children).size() != 0) { 
        for (auto ch: node->children){
            updateLongestEdge(ch, root, totalLeaves, longestLen, breakEdge);
        }
    }
    if (node->identifier == root->identifier) return;
    if (numLeaves == totalLeaves) return;
    if (node->branchLength > longestLen) {
        // std::cout << "DEBUG: branchLength\t" << node->branchLength << "\tLongest\t" << longestLen << '\n';
        longestLen = node->branchLength;
        breakEdge = node;
    }
    return;
}

Node* getLowestChild(Node* node, std::string identifier) {
    size_t remainNumChildren = 0;
    Node* remainChild = nullptr;
    // if (node->children.size() == 0) return node;
    // std::cout << "Lowest: " << node->identifier << '\n';
    for (auto ch: node->children) {
        // std::cout << "Lowest-ch: " << ch->identifier << '\t' << ch->grpID << '\t' << getNumLeaves(ch, node->grpID) << '\t' << identifier<< '\n';
        if (ch->grpID == node->grpID && ch->identifier != identifier && getNumLeaves(ch, node->grpID) > 0) {
            remainNumChildren++;
            // std::cout << "Lowest-ch: " << ch->identifier << '\t' << ch->grpID << '\t' << getNumLeaves(ch, node->grpID)  << '\n';
            remainChild = ch;
        }
        // if (remainNumChildren > 1) return node;
    }
    
    if (remainNumChildren == 1) {
        remainChild = getLowestChild(remainChild, identifier);
        
        return remainChild;
    }
    else {
        return node;
    }
    
}

Node* getLongestEdge(Node* node, int rootID){
    Node* longestEdge = node;
    float longestLen = 0;
    size_t totalLeaves = getNumLeaves(node, node->grpID);
    updateLongestEdge (node, node, totalLeaves, longestLen, longestEdge);
    // std::cout << "DEBUG: getLongestEdge " << longestEdge->identifier << '\t' << longestEdge->branchLength << '\n'; 
    // Update branch Length
    Node* traceback = longestEdge;
    bool downTraversal = false;
    bool doAdjustBranch = false;
    Node* lowestChild1 = nullptr;
    Node* lowestChild2 = nullptr;
    if (longestEdge->parent->parent != nullptr) {
        while (getNumLeaves(traceback->parent, traceback->parent->grpID) == getNumLeaves(longestEdge, longestEdge->grpID) && traceback->parent->grpID == longestEdge->grpID) {
            if (traceback->parent->parent == nullptr) {
                downTraversal = true;
                break;
            }
            traceback = traceback->parent;
        }
    }
    else {
        downTraversal = true;
    }
    if (downTraversal == false && doAdjustBranch == false) {
        if (traceback->identifier == longestEdge->identifier) {
            // check whether parents are dead route, i.e. no other leaves 
            Node* tb_traceback = traceback->parent;
            bool deadparent = true;
            
            while (tb_traceback->parent->grpID == traceback->grpID) {
                
                if (getNumLeaves(tb_traceback->parent, tb_traceback->parent->grpID) != getNumLeaves(tb_traceback, tb_traceback->grpID)) {
                    deadparent = false;
                    break;
                }
                if (tb_traceback == nullptr) break;
                tb_traceback = tb_traceback->parent;
            }
            if (deadparent == true) {
                traceback = getLowestChild(longestEdge->parent, longestEdge->identifier);
                size_t numChildren = 0;
                lowestChild1 = nullptr;
                lowestChild2 = nullptr;

                if (traceback->children.size() > 0) {
                    for (auto ch: traceback->children) {
                        if (ch->grpID == traceback->grpID) {
                            ++numChildren;
                            if (numChildren >= 3) break;
                            if (lowestChild1 == nullptr) lowestChild1 = ch;
                            else lowestChild2 = ch;
                        }

                    }
                    if (numChildren == 2) {
                        // Always add shorter lenght to longer one
                        if (lowestChild1->branchLength <= lowestChild2->branchLength) {
                            Node* temp = lowestChild1;
                            lowestChild1 = lowestChild2;
                            lowestChild2 = temp;
                        }
                        doAdjustBranch = true;
                    }
                }
            }
            else {
                lowestChild1 = getLowestChild(traceback->parent, traceback->identifier);
                lowestChild2 = traceback->parent;
                doAdjustBranch = true;
            }
            
        }
        else {
            lowestChild1 = getLowestChild(traceback->parent, traceback->identifier);
            lowestChild2 = traceback->parent;
            doAdjustBranch = true;
        }
        if (doAdjustBranch) {
            for (auto ch: (lowestChild1->parent)->children) {
                if (ch->identifier == lowestChild1->identifier) {
                    // std::cout << "DEBUG: " << ch->identifier << '+' << lowestChild2->identifier << '\n';
                    // std::cout << "DEBUG: " << ch->branchLength << '+' << lowestChild2->branchLength  << '=' << ch->branchLength + lowestChild2->branchLength<< '\n'; 
                    ch->branchLength += lowestChild2->branchLength;
                }
            }
        }
    }
    else {
        size_t remainNumChildren = 0;
        for (auto ch: longestEdge->parent->children) {
            if (ch->grpID == node->grpID && ch->identifier != longestEdge->identifier) {
                remainNumChildren++;
                if (lowestChild1 == nullptr) {
                    lowestChild1 = getLowestChild(ch, ch->identifier);
                }
                else {
                    lowestChild2 = getLowestChild(ch, ch->identifier);
                }
            }
        }
        if (remainNumChildren == 2) {
            if (lowestChild1->branchLength < lowestChild2->branchLength) {
                Node* temp = lowestChild1;
                lowestChild1 = lowestChild2;
                lowestChild2 = temp;
            }
            for (auto ch: (lowestChild1->parent)->children) {
                if (ch->identifier == lowestChild1->identifier) {
                    // std::cout << "DEBUG: " << ch->branchLength << '+' << lowestChild2->branchLength << '=' << ch->branchLength + lowestChild2->branchLength<< '\n'; 
                    ch->branchLength += lowestChild2->branchLength;
                }
            }
        }
        else if (remainNumChildren == 1) {
            traceback = getLowestChild(lowestChild1, lowestChild1->identifier);
            size_t numChildren = 0;
            lowestChild1 = nullptr;
            lowestChild2 = nullptr;
            if (traceback->children.size() > 0) {
                for (auto ch: traceback->children) {
                    if (ch->grpID == traceback->grpID) {
                        ++numChildren;
                        if (numChildren >= 3) break;
                        if (lowestChild1 == nullptr) lowestChild1 = ch;
                        else lowestChild2 = ch;
                    }
                }
                if (numChildren == 2) {
                    // Always add shorter lenght to longer one
                    if (lowestChild1->branchLength <= lowestChild2->branchLength) {
                        Node* temp = lowestChild1;
                        lowestChild1 = lowestChild2;
                        lowestChild2 = temp;
                    }
                    doAdjustBranch = true;
                }
            }
            for (auto ch: (lowestChild1->parent)->children) {
                if (ch->identifier == lowestChild1->identifier) {
                    // std::cout << "DEBUG: " << ch->identifier << '+' << lowestChild2->identifier << '\n';
                    // std::cout << "DEBUG: " << ch->branchLength << '+' << lowestChild2->branchLength  << '=' << ch->branchLength + lowestChild2->branchLength<< '\n'; 
                    ch->branchLength += lowestChild2->branchLength;
                }
            }
        }
        
    }
    // Check another side
    
    lowestChild1 = nullptr;
    lowestChild2 = nullptr;
    size_t numChildren = 0;
    if (longestEdge->children.size() > 1) {
        for (auto ch: longestEdge->children) {
            if (ch->grpID == longestEdge->grpID) {
                ++numChildren;
                if (numChildren >= 3) break;
                if (lowestChild1 == nullptr) {
                    lowestChild1 = ch;  
                }
                else {
                    lowestChild2 = ch;
                }
            }
            
        }
        if (numChildren == 2) {
            if (lowestChild1->branchLength > lowestChild2->branchLength) {
                Node* temp = lowestChild1;
                lowestChild1 = lowestChild2;
                lowestChild2 = temp;
            }
            for (auto ch: (lowestChild2->parent)->children) {
                if (ch->identifier == lowestChild2->identifier) {
                    // std::cout << "DEBUG: " << ch->branchLength << '+' << lowestChild1->branchLength  << '=' << ch->branchLength + lowestChild1->branchLength<< '\n'; 
                    ch->branchLength += lowestChild1->branchLength;

                }
            }
            

        }
    }
    
    return longestEdge;
}

Node* getBreakingEdge(Node* root, std::string option){
    if (option == "centroid") return getCentroidEdge(root, root->grpID);
    else if (option == "longest") return getLongestEdge(root, root->grpID);
    else {fprintf(stderr, "Error: Wrong partition option: %s\n", option.c_str()); exit(1);} 
}

void setChildrenGrpID(Node*& node, size_t ID_org, int ID) {
    if (node->grpID != ID_org) return;
    node->grpID = ID;
    if (node->children.size() == 0) return;
    for (auto ch: node->children) {
        setChildrenGrpID(ch, ID_org, ID);
    }
    return;
}

void bipartition(Node* root, Node* edge, Node*& tree1Root, Node*& tree2Root, paritionInfo_t*& partition) {
    
    size_t tree1ID = (root->grpID == -1) ? 0 : root->grpID;
    size_t tree2ID = (root->grpID == -1) ? 1 : partition->numPartitions + 1;
    partition->numPartitions += 1;
    // std::cout << "DEBUG: ID " << tree1ID << '\t' << tree2ID << '\n';
    Node* head = edge->parent;
    int headID = edge->parent->grpID;
    while (true){
        if (head->parent == nullptr) break;
        if (head->parent->grpID != headID) break;
        head = head->parent;
    }

    tree1Root = head;
    tree2Root = edge;
    size_t tree2ID_org = tree2Root->grpID;
    size_t tree1ID_org = tree1Root->grpID;
    setChildrenGrpID(tree2Root, tree2ID_org, tree2ID);
    if (tree1Root->grpID == -1) {
        setChildrenGrpID(tree1Root, tree1ID_org, tree1ID);
    }
    // std::cout << "DEBUG: ID1 " << tree2Root->identifier << '\t' << tree2Root->grpID << '\n';
    // std::cout << "DEBUG: ID2 " << tree1Root->identifier << '\t' << tree1Root->grpID << '\n';
    
    
}

void partitionTree(Node* root, paritionInfo_t* partition) {
    size_t totalLeaves = getNumLeaves(root, root->grpID);
    if (totalLeaves <= partition->maxPartitionSize) {
        if (partition->partitionsRoot.empty()) {
            setChildrenGrpID(root, root->grpID, 0);
            size_t numLeaves = getNumLeaves(root, root->grpID);
            partition->partitionsRoot[root->identifier] = std::make_pair(root, numLeaves);
        }
        return;
    }
    // std::cout << "DEBUG: Now partiotion at " << root->identifier << '\n';
    Node* breakEdge = getBreakingEdge(root, partition->partitionOption);
    
    // std::cout << "DEBUG: get breaking edge at " << breakEdge->identifier << '\n';
    Node* tree1 = nullptr;
    Node* tree2 = nullptr;
    bipartition(root, breakEdge, tree1, tree2, partition);
    
    size_t numTree1Leaves = getNumLeaves(tree1, tree1->grpID);
    size_t numTree2Leaves = getNumLeaves(tree2, tree2->grpID);
    // std::cout << "DEBUG: Tree1 " << tree1->identifier << '\t' << numTree1Leaves 
    //           << "\tTree2 "      << tree2->identifier << '\t' << numTree2Leaves << '\n';
    // printf("=============\nInput tree has %ld leaf nodes.\n", totalLeaves);
    // printf("Subtree 1 tree has %ld leaf nodes.\n", numTree2Leaves);
    // printf("Subtree 2 tree has %ld leaf nodes.\n", numTree1Leaves);
    // printf("==============\n");
    if (root->parent == nullptr) {
        partition->partitionsRoot[tree1->identifier] = std::make_pair(tree1, numTree2Leaves);
    }
    partition->partitionsRoot[tree2->identifier] = std::make_pair(tree2, numTree2Leaves);
    partition->partitionsRoot[tree1->identifier].second = numTree1Leaves;
    if (numTree2Leaves > partition->maxPartitionSize) {
        partitionTree(tree2, partition);
    }
    if (numTree1Leaves > partition->maxPartitionSize) {
        partitionTree(tree1, partition);
    }
    
    
    return;
}

void preOrderTraversal(Node* parent, Node* node, Tree*& T, std::unordered_map<std::string, std::pair<Node*, size_t>>& nodes) {
    if (nodes.find(node->identifier) != nodes.end()) {
        Node* NodeCopy;
        if (T->allNodes.size() == 0) {
            NodeCopy = new Node(node->identifier, 0);
            NodeCopy->grpID = -1;
            T->root = NodeCopy;
        }
        else {
            NodeCopy = new Node(node->identifier, parent, 0);
            NodeCopy->grpID = -1;
            // T->allNodes[parent->identifier]->children.push_back(NodeCopy);
        }
        parent = NodeCopy;
        T->allNodes[NodeCopy->identifier] = NodeCopy;
    }
    for (auto ch: node->children) {
        preOrderTraversal(parent, ch, T, nodes);
    }
    return;
}

Tree* reconsturctTree(Node* root, std::unordered_map<std::string, std::pair<Node*, size_t>>& partitionRoots) {
    Tree* T = new Tree();
    Node* nullNode = nullptr;
    preOrderTraversal(nullNode, root, T, partitionRoots);
    return T;
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
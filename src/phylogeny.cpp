#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

void phylogeny::pruneTree(Tree*& T, std::unordered_set<std::string>& seqs) {
    Tree* prunedT = T->prune(seqs);
    delete T;
    T = prunedT;
    return;
}


void preOrderTraversal(phylogeny::Node* parent, phylogeny::Node* node, phylogeny::Tree*& T, std::unordered_map<std::string, std::pair<phylogeny::Node*, size_t>>& nodes) {
    if (nodes.find(node->identifier) != nodes.end()) {
        phylogeny::Node* NodeCopy;
        if (T->allNodes.size() == 0) {
            NodeCopy = new phylogeny::Node(node->identifier, node->branchLength);
            NodeCopy->grpID = -1;
            T->root = NodeCopy;
        }
        else {
            NodeCopy = new phylogeny::Node(node->identifier, parent, node->branchLength);
            NodeCopy->grpID = -1;
        }
        parent = NodeCopy;
        T->allNodes[NodeCopy->identifier] = NodeCopy;
    }
    for (auto ch: node->children) {
        preOrderTraversal(parent, ch, T, nodes);
    }
    return;
}

phylogeny::Tree* phylogeny::constructTreeFromPartitions(Node* root, PartitionInfo* P) {
    Tree* T = new Tree();
    Node* nullNode = nullptr;
    preOrderTraversal(nullNode, root, T, P->partitionsRoot);
    return T;
}

phylogeny::Tree* phylogeny::getPlacementTree(Tree* T) {
    for (auto node: T->allNodes) {
        if (node.second->is_leaf() && node.second->placed) {
            Node* current = node.second;
            while (current->parent != nullptr) {
                if (current->parent->placed) break;
                current->parent->placed = true;
                current = current->parent;
            }
        }
    }
    Tree* placement_T = new Tree();
    Node* rootNode = new Node(T->root);
    placement_T->root = rootNode;
    placement_T->allNodes[rootNode->identifier] = rootNode;
    for (auto node: T->allNodes) {
        if (node.second->placed) {
            Node* copyNode = new Node(node.second);
            for (auto c: node.second->children) 
                if (c->placed) copyNode->children.push_back(c);
            placement_T->allNodes[copyNode->identifier] = copyNode;
        }
    }
    return placement_T;
}


phylogeny::Node* buildInducedTree(phylogeny::Node* origNode, 
                                  std::unordered_map<std::string, std::pair<phylogeny::Node*, size_t>>& targetNodes, 
                                  phylogeny::Tree* T) {
    if (origNode == nullptr) return nullptr;

    std::vector<phylogeny::Node*> validChildrenCopies;
    
    for (auto ch : origNode->children) {
        phylogeny::Node* childCopy = buildInducedTree(ch, targetNodes, T);
        if (childCopy != nullptr) {
            validChildrenCopies.push_back(childCopy);
        }
    }

    bool isTarget = targetNodes.find(origNode->identifier) != targetNodes.end();

    if (isTarget) {
         phylogeny::Node* targetTip = new phylogeny::Node(origNode->identifier, 0.0); 
        targetTip->grpID = -1;
        T->allNodes[targetTip->identifier] = targetTip;

        if (validChildrenCopies.empty()) {
            targetTip->branchLength = origNode->branchLength; 
            return targetTip;
        } else {
            std::string internalID = origNode->identifier + "_LCA"; 
            phylogeny::Node* internalNode = new phylogeny::Node(internalID, origNode->branchLength);
            internalNode->grpID = -1;
            
            targetTip->parent = internalNode;
            internalNode->children.push_back(targetTip);

            for (auto chCopy : validChildrenCopies) {
                chCopy->parent = internalNode;
                internalNode->children.push_back(chCopy);
            }
            
            T->allNodes[internalNode->identifier] = internalNode;
            return internalNode; 
        }
    } 
    else if (validChildrenCopies.size() > 1) {
        phylogeny::Node* nodeCopy = new phylogeny::Node(origNode->identifier, origNode->branchLength);
        nodeCopy->grpID = -1;
        for (auto chCopy : validChildrenCopies) {
            chCopy->parent = nodeCopy;
            nodeCopy->children.push_back(chCopy);
        }
        
        T->allNodes[nodeCopy->identifier] = nodeCopy;
        return nodeCopy;
    } 
    else if (validChildrenCopies.size() == 1) {
        phylogeny::Node* singleChild = validChildrenCopies[0];
        singleChild->branchLength += origNode->branchLength; 
        return singleChild;
    }

    return nullptr;
}

phylogeny::Tree* phylogeny::constructTreeFromPartitions_new(Node* root, PartitionInfo* P) {
    Tree* T = new Tree();
    
    Node* newRoot = buildInducedTree(root, P->partitionsRoot, T);
    
    if (newRoot != nullptr) {
        newRoot->parent = nullptr;
        T->root = newRoot;
    }
    
    return T;
}
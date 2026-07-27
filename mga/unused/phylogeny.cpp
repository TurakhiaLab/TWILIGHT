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
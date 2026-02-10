#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

#include <functional>

std::unique_ptr<phylogeny::Node> phylogeny::Node::deepCopy(std::unordered_map<std::string, Node*>& nodeTracker, int targetGrp) const {
    // Create the new node (Independent memory)
    auto newNode = std::make_unique<Node>(this->identifier, this->branchLength);
    
    // Copy Properties
    newNode->grpID = -1; // Reset as per your logic
    newNode->weight = this->weight;
    
    // Register in the map
    nodeTracker[newNode->identifier] = newNode.get();

    // Recursively Copy Children
    for (const auto& child : this->children) {
        if (child->grpID == targetGrp) {
            // Generate the deep copy of the child
            auto newChild = child->deepCopy(nodeTracker, targetGrp);
            // Linkage: Set Parent/Level
            newChild->parent = newNode.get();
            newChild->level = newNode->level + 1;
            // Ownership: Move child into this node's vector
            newNode->children.push_back(std::move(newChild));
        }
    }
    return newNode;
}


void phylogeny::Node::collectPostOrder(std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(this); 
    Node* current; 
    while (!s1.empty()) {
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == current->grpID) s1.push(current->children[i].get());     
        }
    } 
    return;
}
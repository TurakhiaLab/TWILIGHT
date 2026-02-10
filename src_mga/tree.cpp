#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

#include <fstream>
#include <algorithm>
#include <queue>
#include <limits>
#include <chrono>
#include <iomanip>
#include <functional>
#include <cassert>


void stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0, temp_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        std::string sub;
        if (temp_pos == 0) {
            sub = s.substr(start_pos, end_pos-start_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 1) {
                temp_pos = start_pos;
            }
            else {
                words.emplace_back(sub);
            }
        }
        else {
            sub = s.substr(temp_pos, end_pos-temp_pos);
            if (std::count(sub.begin(), sub.end(), '\'') % 2 == 0) {
                temp_pos = 0;
                words.emplace_back(sub);
            }
        }
        // words.emplace_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
}

std::string stripString(std::string s){
    while(s.length() && s[s.length() - 1] == ' '){
        s.pop_back();
    }
    for(size_t i = 0; i < s.length(); i++){
        if(s[i] != ' '){
            return s.substr(i);
        }
    }
    return s;
}

void phylogeny::Tree::parseNewick(std::string& newickString) {

    newickString = stripString(newickString);
    
    // Reset Root
    this->root = nullptr; 
    this->allNodes.clear();

    std::vector<std::string> segments;
    stringSplit(newickString, ',', segments);

    std::vector<size_t> numOpen; numOpen.reserve(segments.size());
    std::vector<size_t> numClose; numClose.reserve(segments.size());
    std::vector<std::string> leaves; leaves.reserve(segments.size());
    std::vector<std::queue<float>> branchLen(128);

    size_t level = 0;

    for (const auto& s : segments) {
        size_t no = 0;
        size_t nc = 0;
        size_t leafDepth = 0;
        bool stop = false;
        bool branchStart = false;
        bool nameZone = false;
        bool hasApo = false;
        std::string leaf = "", branch = "";
        
        for (auto c: s) {
            if (nameZone) {
                leaf += c;
                if (c == '\'') nameZone = false;
            } else if (c == '\'' && !nameZone) {
                nameZone = true;
                hasApo = true;
                leaf += c;
            } else if (c == ':') {
                stop = true;
                branch = "";
                branchStart = true;
            } else if (c == '(') {
                no++;
                level++;
                if (branchLen.size() <= level) {
                    branchLen.resize(level*2);
                }
            } else if (c == ')') {
                stop = true;
                nc++;
                // float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
                float len = (branch.size() > 0) ? std::stof(branch) : 0.0;
                if (len == 0) len = 1.0;
                branchLen[level].push(len);
                level--;
                branchStart = false;
            } else if (!stop) {
                leaf += c;
                branchStart = false;
                leafDepth = level;

            } else if (branchStart) {
                if (isdigit(c)  || c == '.' || c == 'e' || c == 'E' || c == '-' || c == '+') {
                    branch += c;
                }
            }
        }
        if (hasApo && leaf[0] == '\'' && leaf[leaf.length()-1] == '\'') leaf = leaf.substr(1, leaf.length()-2);

        leaves.push_back(std::move(leaf));
        numOpen.push_back(no);
        numClose.push_back(nc);
        
        float len = (branch.empty()) ? 0.0f : std::stof(branch);
        branchLen[level].push(len);

        m_maxDepth = std::max(m_maxDepth, leafDepth);
        m_meanDepth += leafDepth;
    }

    if (!leaves.empty()) m_meanDepth /= leaves.size();
    if (level != 0) throw std::runtime_error("ERROR: Incorrect Newick format!");

    m_numLeaves = leaves.size();

    // Build Tree with unique_ptr
    
    // Stack holds RAW pointers. These are "Observers". 
    // They point to memory OWNED by the tree structure.
    std::stack<Node*> parentStack;
    for (size_t i = 0; i < leaves.size(); ++i) {
        auto leafName = leaves[i];
        const auto no = numOpen[i];
        const auto nc = numClose[i];

        // 1. Create Internal Nodes
        for (size_t j = 0; j < no; ++j) {
            std::string nid = newInternalNodeId();
            float currentLen = branchLen[level].front();
            
            // Create the node (it is not attached to anything yet)
            auto newNode = std::make_unique<Node>(nid, currentLen);
            Node* newNodePtr = newNode.get(); // Keep a raw handle

            // Link Parent/Child
            if (parentStack.empty()) {
                // This is the root
                newNode->level = 1;
                // Move ownership to Tree::root
                this->root = std::move(newNode); 
            } else {
                // This is a child
                Node* parent = parentStack.top();
                newNode->parent = parent;
                newNode->level = parent->level + 1;
                // Move ownership into Parent's children vector
                parent->children.push_back(std::move(newNode));
            }

            // Store raw pointer for lookups
            newNodePtr->grpID = -1;
            allNodes[nid] = newNodePtr;
            parentStack.push(newNodePtr);
            
            branchLen[level].pop();
            level++;
        }

        // 2. Handle Leaf Duplicates
        if (allNodes.find(leafName) != allNodes.end()) {
            std::cerr << "WARNING: Duplicate leaf: " << leafName << "\n";
            leafName += "_dup_" + std::to_string(allNodes.size());
        }

        // 3. Create Leaf Node
        float leafLen = branchLen[level].front();
        auto leafNode = std::make_unique<Node>(leafName, leafLen);
        Node* leafNodePtr = leafNode.get();

        if (!parentStack.empty()) {
            Node* parent = parentStack.top();
            leafNodePtr->parent = parent;
            leafNodePtr->level = parent->level + 1;
            parent->children.push_back(std::move(leafNode));
        } else {
            // Edge case: Tree with single node
            leafNodePtr->level = 1;
            this->root = std::move(leafNode);
        }

        leafNodePtr->grpID = -1;
        allNodes[leafName] = leafNodePtr;
        branchLen[level].pop();

        // 4. Backtrack
        for (size_t j = 0; j < nc; ++j) {
            parentStack.pop();
            level--;
        }
    }

    if (!this->root) throw std::runtime_error("Tree found empty!");

    this->root->branchLength = 0;

    // --- PHASE 3: Post Processing ---
    float minBrLen = std::numeric_limits<float>::max();
    bool anyPositive = false;

    for (const auto& [id, node] : allNodes) {
        if (node->branchLength > 0) {
            minBrLen = std::min(minBrLen, node->branchLength);
            anyPositive = true;
        }
    }

    for (auto& [id, node] : allNodes) {
        if (node->identifier == this->root->identifier) continue;
        if (!anyPositive) node->branchLength = 1.0f;
        else if (node->branchLength <= 0.0f) node->branchLength = minBrLen;
    }

    this->calLeafNum();
    this->calSeqWeight();

    // Recalculate levels for the whole tree
    std::function<void(Node*, size_t)> fixLevels = [&](Node* n, size_t lvl) {
        n->level = lvl;
        for (const auto& child : n->children) fixLevels(child.get(), lvl + 1);
    };
    fixLevels(this->root.get(), 1);
}

phylogeny::Tree::Tree(const std::string& treeFileName) {
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::ifstream inputStream(treeFileName);
    if (!inputStream.is_open()) { 
        throw std::runtime_error("Error: Failed to open tree file: " + treeFileName);
    }
    std::string newick; 
    std::getline(inputStream, newick);
    
    this->parseNewick(newick);

    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(treeBuiltEnd - treeBuiltStart);
    std::cerr << "Newick string read in: " << duration.count() << " ms\n";
}
  
phylogeny::Tree::Tree(Node* sourceNode, bool reroot) {
    if (!sourceNode) throw std::runtime_error("Source node is null");
    // 1. Prepare
    this->allNodes.clear();
    int targetGrp = sourceNode->grpID;

    // 2. Deep Copy the Subtree
    this->root = sourceNode->deepCopy(this->allNodes, targetGrp);

    // 3. Final Root Setup
    this->root->level = 1;
    this->root->parent = nullptr; // Root has no parent
    
    // 4. Finalize
    if (reroot) {
        this->reroot();
    } else {
        this->calLeafNum();
        this->calSeqWeight();
    }
}

void phylogeny::Tree::calLeafNum() {
    if (!this->root) return;
    // postorder traversal
    std::stack<Node*> postorder;
    this->root->collectPostOrder(postorder);
    Node* current;
    // get numbers of leaves
    while (!postorder.empty()) {
        current = postorder.top(); 
        postorder.pop();
        if (current->is_leaf()) this->allNodes[current->identifier]->numLeaves = 1;
        else {
            int leaves = 0;
            for (const auto& ch: current->children) leaves += ch->numLeaves;
            this->allNodes[current->identifier]->numLeaves = leaves;
        }
    }
    this->m_numLeaves = this->root->numLeaves;
    return;
};

void phylogeny::Tree::calSeqWeight() {
    float maxWeight = 0;
    for (auto node: this->allNodes) {
        if (!node.second->is_leaf()) continue;
        float w = 0;
        Node* current = node.second;
        while (true) {
            w += current->branchLength / current->numLeaves;
            // w += 1.0 / current->numLeaves;
            current = current->parent;
            if (current == nullptr) break; 
        }
        // w = 1.0;
        this->allNodes[node.second->identifier]->weight = w;
        if (w > maxWeight) maxWeight = w;
    }
    float normFactor = maxWeight / 1.0;
    for (auto node: this->allNodes) {
        if (!node.second->is_leaf()) continue;
        this->allNodes[node.second->identifier]->weight /= normFactor;
    }
    return;
};

void phylogeny::Tree::showTree() {
    std::function<void(const std::unique_ptr<Node>&)> showTreeImpl = [&](const std::unique_ptr<Node>& node) {
        std::cerr << std::left << std::setw(12) << node->identifier     // Identifier
                  << std::right << std::setw(10) << std::fixed << std::setprecision(4) << node->branchLength  // Branch length
                  << std::setw(10) << (node->parent ? node->parent->identifier : "ROOT")  // Parent or ROOT
                  << std::setw(8)  << node->level
                  << std::setw(8)  << node->grpID
                  << std::setw(10) << std::fixed << std::setprecision(3) << node->children.size()
                  << '\n';
        for (auto& c : node->children) {
            showTreeImpl(c);
        }
    };
    std::cerr << std::left << std::setw(12) << "Identifier"
              << std::right << std::setw(10) << "Length"
              << std::setw(10) << "Parent"
              << std::setw(8)  << "Level"
              << std::setw(8)  << "Group"
              << std::setw(10) << "Weight" << '\n';
    std::cerr << std::string(50, '-') << '\n';
    if (this->root) showTreeImpl(this->root);
}

phylogeny::Tree* phylogeny::Tree::prune(std::unordered_set<std::string>& seqs) {
    
    auto treePruneStart = std::chrono::high_resolution_clock::now();

    std::unordered_map<std::string, bool> keepMap;
    
    std::function<bool(const Node*)> markNodes = [&](const Node* node) -> bool {
        if (!node) return false;
        bool keep = false;
        if (node->is_leaf()) {
            keep = (seqs.find(node->identifier) != seqs.end());
        } else {
            for (const auto& child : node->children) {
                if (markNodes(child.get())) {
                    keep = true;
                }
            }
        }
        keepMap[node->identifier] = keep;
        return keep;
    };
    if (this->root) markNodes(this->root.get());
    
    Tree* pT = new Tree();
    std::function<std::unique_ptr<Node>(const Node*, bool)> rebuild = [&](const Node* curr, bool isRoot) -> std::unique_ptr<Node> {
        // If not marked, this whole branch is pruned
        if (!keepMap[curr->identifier]) return nullptr;

        // Recursively process children
        std::vector<std::unique_ptr<Node>> validChildren;
        for (const auto& child : curr->children) {
            auto keptChild = rebuild(child.get(), false);
            if (keptChild) {
                validChildren.push_back(std::move(keptChild));
            }
        }

        // Case A: The Root Node
        if (isRoot) {
            auto newNode = std::make_unique<Node>(curr->identifier, curr->branchLength);
            newNode->grpID = -1; 
            newNode->level = 1;
            
            // Attach children
            for (auto& child : validChildren) {
                child->parent = newNode.get();
                child->level = newNode->level + 1;
                newNode->children.push_back(std::move(child));
            }
            
            // Register in map
            pT->allNodes[newNode->identifier] = newNode.get();
            return newNode;
        }

        // Case B: Leaf Node (Must be in 'seqs')
        if (validChildren.empty()) {
            auto newNode = std::make_unique<Node>(curr->identifier, curr->branchLength);
            newNode->grpID = -1;
            pT->allNodes[newNode->identifier] = newNode.get();
            return newNode;
        }

        // Case C: Path Compression (Internal Node with exactly 1 child)
        // If this node has 1 child, we remove this node and extend the child's branch.
        // This effectively "pulls" the child up to the parent.
        if (validChildren.size() == 1) {
            auto child = std::move(validChildren[0]);

            // Absorb current branch length into the child
            child->branchLength += curr->branchLength;
            
            // We do NOT add 'curr' to pT->allNodes because 'curr' effectively ceases to exist.
            // We return the child directly to be attached to 'curr's parent.
            return child;
        }

        // Case D: Junction Node (Internal Node with > 1 children)
        // We keep this node as a structural branch point.
        auto newNode = std::make_unique<Node>(curr->identifier, curr->branchLength);
        newNode->grpID = -1;
        
        for (auto& child : validChildren) {
            child->parent = newNode.get();
            child->level = newNode->level + 1;
            newNode->children.push_back(std::move(child));
        }
        
        pT->allNodes[newNode->identifier] = newNode.get();
        return newNode;
    };
    
    pT->root = rebuild(this->root.get(), true);

    if (!pT->root) {
        std::cerr << "ERROR: No sequences from the input sequence file are found in the tree file.\n";
        exit(1);
    }

    pT->calLeafNum(); 
    pT->calSeqWeight();

    if (pT->m_numLeaves != seqs.size()) {
        std::cerr << "WARNING: " << (seqs.size() - pT->m_numLeaves) 
                  << " sequences are missing from the tree and will be ignored.\n";
    }
    
    std::cerr << "Number of Leaves: " << this->m_numLeaves 
              << " (before) -> " << pT->m_numLeaves << " (after)\n";
    auto treePruneEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treePruneTime = treePruneEnd - treePruneStart;
    std::cerr << "Tree pruned in: " <<  treePruneTime.count() / 1000000 << " ms\n";

    return pT;
}

std::string phylogeny::Tree::getNewickString() {
    std::string newickString = "";
    std::function<void(std::unique_ptr<Node>&)> buildNewick = [&](std::unique_ptr<Node>& node) {
        if(node->children.size() != 0) {
	    	newickString += "(";
	    	for(int n = 0; n < node->children.size(); ++n) {
	    		if(n != 0) newickString += ",";
	    		buildNewick(node->children[n]);
	    	}
	    	if (node->parent != nullptr) newickString += ("):" + std::to_string(node->branchLength));
            else newickString += ")";
        }
	    else {
            auto name = node->identifier;
            bool specialChar = (name.find(',') != std::string::npos)  || \
                               (name.find(':') != std::string::npos)  || \
                               (name.find('(') != std::string::npos)  || \
                               (name.find(')') != std::string::npos);
            if (specialChar) name = '\'' + name + '\'';
            newickString += (name + ':' + std::to_string(node->branchLength));
        }
    };
    buildNewick(this->root);
    newickString += ";";
    return newickString;
}


void phylogeny::Tree::convert2binaryTree() {
    std::stack<Node*> postOrder;
    this->root->collectPostOrder(postOrder);

    while (!postOrder.empty()) {
        Node* node = postOrder.top();
        postOrder.pop(); 

        // ---------------------------------------------------------
        // CASE 1: Binarization (Node has > 2 children)
        // ---------------------------------------------------------
        if (node->children.size() > 2) {
            int grp = node->grpID;
            // 1. Extract all children
            std::vector<std::unique_ptr<Node>> currentLayer = std::move(node->children);
            node->children.clear(); // Vector is now empty
            // 2. Iteratively group them in pairs
            while (currentLayer.size() > 2) {
                std::vector<std::unique_ptr<Node>> nextLayer;
                for (size_t i = 0; i < currentLayer.size() - 1; i += 2) {
                    // Create intermediate node
                    auto name = this->newInternalNodeId();
                    auto newNode = std::make_unique<Node>(name, 0.0);
                    newNode->grpID = grp;
                    Node* newNodePtr = newNode.get();

                    // Register in lookup map
                    this->allNodes[name] = newNodePtr;

                    // Link children to this new intermediate node
                    // Child 1
                    currentLayer[i]->parent = newNodePtr;
                    newNode->children.push_back(std::move(currentLayer[i]));
                    
                    // Child 2
                    currentLayer[i+1]->parent = newNodePtr;
                    newNode->children.push_back(std::move(currentLayer[i+1]));

                    // Add new node to next layer
                    nextLayer.push_back(std::move(newNode));
                }

                // Handle odd one out (just carry it over to the next layer)
                if (currentLayer.size() % 2 == 1) {
                    nextLayer.push_back(std::move(currentLayer.back()));
                }

                currentLayer = std::move(nextLayer);
            }

            assert(currentLayer.size() == 2);

            // 3. Attach the final 2 nodes back to the original parent
            for (auto& child : currentLayer) {
                child->parent = node;
                node->children.push_back(std::move(child));
            }
        }
        
        // ---------------------------------------------------------
        // CASE 2: Remove Useless Intermediate Nodes (1 child)
        // ---------------------------------------------------------
        else if (node->children.size() == 1 && node->parent != nullptr) { 
            Node* parent = node->parent;
            
            // Find 'node' in the parent's children vector
            auto it = std::find_if(parent->children.begin(), parent->children.end(),
                [&](const std::unique_ptr<Node>& ptr) { return ptr.get() == node; });

            if (it != parent->children.end()) {
                // 1. Steal the grandchild (ownership transfer)
                auto grandchild = std::move(node->children[0]);
                
                // 2. Update grandchild properties (absorb branch length)
                grandchild->branchLength += node->branchLength;
                grandchild->parent = parent;

                // 3. Remove 'node' from map (pointer becomes invalid soon)
                this->allNodes.erase(node->identifier);

                // 4. Overwrite 'node' with 'grandchild' in the parent's vector
                // This implicitly destroys 'node'
                *it = std::move(grandchild);
            }
        }
        
        // ---------------------------------------------------------
        // CASE 3: Clean Empty Internal Nodes
        // ---------------------------------------------------------
        else if (node->children.empty() && !node->is_leaf()) {
             if (node->parent != nullptr) {
                Node* parent = node->parent;

                // Find and remove 'node' from parent
                auto it = std::find_if(parent->children.begin(), parent->children.end(),
                    [&](const std::unique_ptr<Node>& ptr) { return ptr.get() == node; });

                if (it != parent->children.end()) {
                    this->allNodes.erase(node->identifier);
                    parent->children.erase(it); // Destroys node
                }
            }
        }
    }

    // Recalculate levels for the whole tree
    std::function<void(Node*, size_t)> fixLevels = [&](Node* n, size_t lvl) {
        n->level = lvl;
        for (const auto& child : n->children) fixLevels(child.get(), lvl + 1);
    };
    fixLevels(this->root.get(), 1);
}

void phylogeny::Tree::reroot(bool placement)
{
    // -----------------------------------------------------
    // 1. Pre-calculation & Binarization
    // -----------------------------------------------------
    int beforeConvert = 0, beforeReroot = 0, afterReroot = 0;
    for (const auto& n : this->allNodes) beforeConvert = std::max(beforeConvert, (int)n.second->level);

    this->convert2binaryTree();
    
    for (const auto& n : this->allNodes) beforeReroot = std::max(beforeReroot, (int)n.second->level);

    // -----------------------------------------------------
    // 2. Select Start Node for BFS
    // -----------------------------------------------------
    Node* start = nullptr;
    if (placement) {
        for (const auto& n : this->allNodes) {
            if (n.second->is_leaf()) { start = n.second; break; }
        }
    } else {
        // Fallback or default
        for (const auto& n : this->allNodes) {
            if (n.second->is_leaf()) { start = n.second; break; }
        }
    }

    if (!start) {
        std::cerr << "WARNING: Could not find a valid start node for rerooting.\n";
        return;
    }

    // -----------------------------------------------------
    // 3. BFS Function (Unchanged logic, just cleanup)
    // -----------------------------------------------------
    auto bfs = [&](Node* startNode, std::unordered_map<Node*, Node*>& parentOut) -> Node* {
        std::queue<Node*> q;
        std::unordered_map<Node*, int> dist;
        
        parentOut.clear();
        q.push(startNode);
        dist[startNode] = 0;
        parentOut[startNode] = nullptr;
        
        Node* farthest = startNode;

        while (!q.empty()) {
            Node* u = q.front(); q.pop();

            // Neighbors = Children + Parent (Treat as undirected graph)
            // Note: We use raw pointers here, which is safe for traversal
            std::vector<Node*> neighbors;
            for (const auto& child : u->children) neighbors.push_back(child.get());
            if (u->parent) neighbors.push_back(u->parent);

            for (Node* v : neighbors) {
                if (dist.find(v) == dist.end()) {
                    dist[v] = dist[u] + 1;
                    parentOut[v] = u;
                    q.push(v);

                    int d_f = dist[farthest];
                    int d_v = dist[v];
                    
                    if (d_v > d_f) {
                        farthest = v;
                    }
                }
            }
        }
        return farthest;
    };

    // -----------------------------------------------------
    // 4. Find Diameter (A -> B) and Midpoint
    // -----------------------------------------------------
    std::unordered_map<Node*, Node*> parentA, parentB;
    Node* A = bfs(start, parentA); // Farthest from start
    Node* B = bfs(A, parentB);     // Farthest from A (Diameter end)

    // Reconstruct path from B back to A
    std::vector<Node*> path;
    for (Node* cur = B; cur != nullptr; cur = parentB[cur]) {
        path.push_back(cur);
    }
    // Now path is [B, ..., A]. The center is simply size/2
    Node* newRootNode = path[path.size() / 2];

    if (newRootNode == this->root.get()) return; // Already rooted correctly

    // -----------------------------------------------------
    // 5. Invert Ownership Path (The Hard Part)
    // -----------------------------------------------------
    
    // Helper to steal a child from its parent
    auto detachNode = [&](Node* child) -> std::unique_ptr<Node> {
        Node* parent = child->parent;
        
        // If parent is null, this child is the current Tree Root
        if (!parent) {
            return std::move(this->root);
        }

        // Find child in parent's vector
        auto& siblings = parent->children;
        auto it = std::find_if(siblings.begin(), siblings.end(), 
            [&](const std::unique_ptr<Node>& ptr){ return ptr.get() == child; });
        
        if (it == siblings.end()) throw std::runtime_error("Topology error: child not found in parent");

        // Move out and erase
        std::unique_ptr<Node> ptr = std::move(*it);
        siblings.erase(it);
        return ptr;
    };

    // We walk UP from newRoot to oldRoot.
    // At each step, we detach the node, and make it the PARENT of its old parent.
    
    // 1. Detach the New Root from the tree. It is now floating.
    float nextBranchLength = newRootNode->branchLength; // Save length to apply to the edge overlap
    Node* curr = newRootNode;
    Node* oldParent = curr->parent;

    // We must hold the new root safely. 
    // Note: We can't assign to this->root yet because we need to traverse the old root structure.
    std::unique_ptr<Node> newRootPtr = detachNode(curr);

    // 2. Traverse up
    while (oldParent != nullptr) {
        // Save the *next* parent up the chain
        Node* nextParent = oldParent->parent;
        float tempLen = oldParent->branchLength;

        // Detach oldParent from ITS parent (or root)
        std::unique_ptr<Node> oldParentPtr = detachNode(oldParent);

        // INVERT RELATIONSHIP:
        // oldParent becomes a child of curr
        oldParent->parent = curr;
        oldParent->branchLength = nextBranchLength; // Give it the length of the edge we just traversed
        curr->children.push_back(std::move(oldParentPtr));

        // Move pointers up
        curr = oldParent;
        oldParent = nextParent;
        nextBranchLength = tempLen;
    }

    // 3. Finalize New Root
    newRootNode->parent = nullptr;
    newRootNode->branchLength = 0.0;
    
    // Assign to Tree
    this->root = std::move(newRootPtr);

    // -----------------------------------------------------
    // 6. Post-Processing
    // -----------------------------------------------------
    
    // Recalculate levels for the whole tree
    std::function<void(Node*, size_t)> fixLevels = [&](Node* n, size_t lvl) {
        n->level = lvl;
        for (const auto& child : n->children) fixLevels(child.get(), lvl + 1);
    };
    fixLevels(this->root.get(), 1);

    // Re-binarize if necessary (rerooting can create a ternary root)
    this->convert2binaryTree();
    this->calLeafNum();
    this->calSeqWeight();

    for (const auto& n : this->allNodes) afterReroot = std::max(afterReroot, (int)n.second->level);

    std::cerr << "======== Tree Depth ========\n";
    std::cerr << "Original: " << beforeConvert << '\n';
    std::cerr << "Binary:   " << beforeReroot << '\n';
    std::cerr << "Reroot:   " << afterReroot << '\n';
}
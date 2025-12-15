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
    Node* treeRoot = nullptr;

    std::vector<std::string> leaves;
    std::vector<size_t> numOpen;
    std::vector<size_t> numClose;
    std::vector<std::queue<float>> branchLen (128);  // will be resized later if needed
    size_t level = 0;

    std::vector<std::string> s1;
    stringSplit(newickString, ',', s1);

    numOpen.reserve(s1.size());
    numClose.reserve(s1.size());
    

    for (auto s: s1) {
        size_t no = 0;
        size_t nc = 0;
        size_t leafDepth = 0;

        bool stop = false;
        bool branchStart = false;
        bool nameZone = false;
        bool hasApo = false;
        std::string leaf = "";
        std::string branch = "";

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
        // float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
        float len = (branch.size() > 0) ? std::stof(branch) : 0.0;
        // if (len == 0) len = 1.0;
        branchLen[level].push(len);

        // Adjusting max and mean depths
        m_maxDepth = std::max(m_maxDepth, leafDepth);
        m_meanDepth += leafDepth;

    }

    m_meanDepth /= leaves.size();

    if (level != 0) {
        fprintf(stderr, "ERROR: incorrect Newick format!\n");
        exit(1);
    }

    m_numLeaves = leaves.size();

    std::stack<Node*> parentStack;

    for (size_t i=0; i<leaves.size(); i++) {
        auto leaf = leaves[i];
        auto no = numOpen[i];
        auto nc = numClose[i];
        for (size_t j=0; j<no; j++) {
            std::string nid = newInternalNodeId();
            Node* newNode = nullptr;
            if (parentStack.size() == 0) {
                newNode = new Node(nid, branchLen[level].front());
                treeRoot = newNode;
            } else {
                newNode = new Node(nid, parentStack.top(), branchLen[level].front());
                // if (branchLen[level].front() < 0.001) std::cout << nid << '\t' << branchLen[level].front() << '\n';
        
            }
            branchLen[level].pop();
            level++;

            /* Group Id */
            newNode->grpID = -1;

            allNodes[nid] = newNode;
            parentStack.push(newNode);
        }
        if (allNodes.find(leaf) != allNodes.end()) {
            printf("WARNING: duplicate leaf names found in the tree! Leaf name: %s. All duplicate leaves will be removed during processing.\n", leaf.c_str());
            leaf += "_dup_" + std::to_string(allNodes.size());
        }
        
        Node* leafNode = new Node(leaf, parentStack.top(), branchLen[level].front());
        // if (branchLen[level].front() < 0.001) std::cout << leaf << '\t' << branchLen[level].front() << '\n';
        /* Group Id */
        leafNode->grpID = -1;
        

        allNodes[leaf] = leafNode;

        branchLen[level].pop();
        for (size_t j=0; j<nc; j++) {
            parentStack.pop();
            level--;
        }
    }

    if (treeRoot == nullptr) {
        fprintf(stderr, "WARNING: Tree found empty!\n");
    }

    treeRoot->branchLength = 0;
    root = treeRoot;

    float minBrLen = std::numeric_limits<float>::max();
    bool allZero = true;
    for (auto node: this->allNodes) {
        if (node.second->branchLength > 0 && node.second->branchLength < minBrLen) minBrLen = node.second->branchLength;
        if (allZero) if(node.second->branchLength > 0) allZero = false;
    }
    if (allZero) {
        for (auto node: this->allNodes) {
            if (node.second->identifier != this->root->identifier) node.second->branchLength = 1.0;
        }
    }
    else {
        for (auto node: this->allNodes) {
            if (node.second->identifier != this->root->identifier && node.second->branchLength == 0) node.second->branchLength = minBrLen;
        }
    }

    this->calLeafNum();
    this->calSeqWeight();
}

phylogeny::Tree::Tree (std::string treeFileName)
{
    auto treeBuiltStart = std::chrono::high_resolution_clock::now();
    std::ifstream inputStream(treeFileName);
    if (!inputStream) { fprintf(stderr, "Error: Failed to open file: %s\n", treeFileName.c_str()); exit(1); }
    std::string newick; 
    // inputStream >> std::noskipws >> newick;
    std::getline(inputStream, newick);
    this->parseNewick(newick);
    auto treeBuiltEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treeBuiltTime = treeBuiltEnd - treeBuiltStart;
    std::cerr << "Newick string read in: " <<  treeBuiltTime.count() / 1000000 << " ms\n";
}
  
phylogeny::Tree::Tree(Node* node, bool reroot) {
    Node* root = new Node(node->identifier, node->branchLength);
    int grp = node->grpID;
    root->grpID = -1;
    this->allNodes[node->identifier] = root;
    this->root = root;
    std::stack<Node*> s1;
    s1.push(node); 
    Node* current; 
    while (!s1.empty()) { 
        current = s1.top();
        if (current->identifier != this->root->identifier) {
            Node* copyNode = new Node(current->identifier, this->allNodes[current->parent->identifier], current->branchLength);
            copyNode->grpID = -1; copyNode->level = current->level - (node->level-1);
            copyNode->weight = current->weight; // copyNode->numLeaves = current->numLeaves;
            this->allNodes[current->identifier] = copyNode;
        }
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == grp) s1.push(current->children[i]);     
        }
    }
    int max_interNode = 0;
    for (auto a: this->allNodes) {
        if (!a.second->is_leaf()) max_interNode = std::max(std::stoi(a.first.substr(5)), max_interNode);
    }
    this->m_currInternalNode = max_interNode;
    if (reroot) this->reroot();
    else {
        this->calLeafNum();
        this->calSeqWeight();
    }
    return;
}

phylogeny::Tree::Tree(std::unordered_set<std::string>& seqNames) {
    Node* rootNode = new Node("node_1", 0);
    rootNode->grpID = 0;
    for (auto name: seqNames) {
        Node* node = new Node(name, rootNode, 1.0);
        node->weight = 1.0;
        node->grpID = rootNode->grpID;
        this->allNodes[name] = node;
    }
    this->root = rootNode;
    this->allNodes[rootNode->identifier] = rootNode;
    
}

phylogeny::Tree::~Tree() {
    for (auto n: this->allNodes) {
        delete n.second;
    }
    this->allNodes.clear();
}

void phylogeny::Tree::calLeafNum() {
    
    // postorder traversal
    std::stack<Node*> s;
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
            for (auto ch: current->children) leaves += ch->numLeaves;
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
    std::function<void(Node*)> showTreeImpl = [&](Node* node) {
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

    Tree* pT = new Tree();
    pT->root = new Node(this->root->identifier, this->root->branchLength);
    pT->root->grpID = -1;
    pT->allNodes[pT->root->identifier] = pT->root;
    
    std::unordered_map<std::string, bool> keepNode;
    
    for (const auto& n : this->allNodes) {
        if (n.second->is_leaf()) {
            keepNode[n.second->identifier] = (seqs.find(n.second->identifier) != seqs.end());
        }
    }

    std::function<bool(Node*)> hasKeepDescendant = [&](Node* node) -> bool {
        if (node->is_leaf()) return keepNode[node->identifier];
        bool keep = false;
        for (Node* child : node->children) {
            if (hasKeepDescendant(child)) {
                keep = true;
            }
        }
        keepNode[node->identifier] = keep;
        return keep;
    };
    hasKeepDescendant(this->root);
    
    std::function<void(Node*, Node*)> buildPrunedTree = [&](Node* origNode, Node* newParent) {
        if (!keepNode[origNode->identifier]) return;
        
        if (origNode->identifier == this->root->identifier) {
            for (Node* child : this->root->children) {
                buildPrunedTree(child, this->root);
            }
            return;
        } 
        else {
            // Process children
            std::vector<Node*> keepChildren;
            for (Node* child : this->allNodes[origNode->identifier]->children) {
                if (keepNode[child->identifier]) {
                    keepChildren.push_back(child);
                }
            }

            // Handle nodes with only one child (merge them)
            if (keepChildren.size() == 0) {
                if (origNode->is_leaf()) {
                    Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                    newNode->grpID = -1;
                    pT->allNodes[newNode->identifier] = newNode;
                }
                else return;
            }
            else if (keepChildren.size() == 1) {
                Node* onlyChild = keepChildren[0];
                float combinedLength = origNode->branchLength;
                while (true) {
                    std::vector<Node*> temp;
                    combinedLength += onlyChild->branchLength;
                    for (Node* child : onlyChild->children) {
                        if (keepNode[child->identifier]) {
                            temp.push_back(child);
                        }
                    }
                    if (temp.size() > 1) {
                        Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                        newChild->grpID = -1;
                        pT->allNodes[newChild->identifier] = newChild;
                        break;
                    }
                    else if (temp.empty()) {
                        if (onlyChild->is_leaf()) {
                            Node* newChild = new Node(onlyChild->identifier, pT->allNodes[newParent->identifier], combinedLength);
                            newChild->grpID = -1;
                            pT->allNodes[newChild->identifier] = newChild;
                            break;
                        }
                        else return;
                    }
                    onlyChild = temp[0];
                }
                // Process the child's children
                for (Node* grandchild : this->allNodes[onlyChild->identifier]->children) {
                    buildPrunedTree(grandchild, this->allNodes[onlyChild->identifier]);
                }
            } 
            else {
                Node* newNode = new Node(origNode->identifier, pT->allNodes[newParent->identifier], origNode->branchLength);
                newNode->grpID = -1;
                pT->allNodes[newNode->identifier] = newNode;
                for (Node* child : this->allNodes[origNode->identifier]->children) {
                    buildPrunedTree(child, this->allNodes[origNode->identifier]);
                }
            }
        }
    };
    
    buildPrunedTree(pT->root, nullptr);
    
    // Update tree properties
    pT->m_numLeaves = 0;
    for (const auto& nodePair : pT->allNodes) {
        if (nodePair.second->is_leaf()) {
            pT->m_numLeaves++;
        }
    }
    
    pT->calLeafNum();
    pT->calSeqWeight();

    std::cerr << "Number of Leaves: " << this->m_numLeaves << " (before pruning) -> " << pT->m_numLeaves << " (after pruning)\n";
    if (pT->m_numLeaves == 0) {
        std::cerr << "ERROR: No sequences from the input sequence file are found in the tree file.\n";
        exit(1);
    }
    if (pT->m_numLeaves != seqs.size()) std::cerr << "WARNING: " << (seqs.size() - pT->m_numLeaves) << " sequences are missing from the tree and will be ignored.\n";
    auto treePruneEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds treePruneTime = treePruneEnd - treePruneStart;
    std::cerr << "Tree pruned in: " <<  treePruneTime.count() / 1000000 << " ms\n";
    return pT;
}

std::string phylogeny::Tree::getNewickString() {
    std::string newickString = "";
    std::function<void(Node*)> buildNewick = [&](Node* node) {
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

void phylogeny::updateSubrootInfo(Node*& subroot, Tree* subT, int subtreeIdx) {
    subroot->seqsIncluded.push_back(subtreeIdx);
    subroot->alnLen = subT->root->alnLen;
    subroot->alnNum = subT->root->seqsIncluded.size();
    subroot->msaFreq = subT->root->msaFreq;
    subroot->alnWeight = subT->root->alnWeight;
    return;
}

void phylogeny::Tree::convert2binaryTree() {
    std::stack<Node*> postOrder;
    this->root->collectPostOrder(postOrder);
    int maxDepth = this->m_maxDepth;
    while(!postOrder.empty()) {
        Node* node = postOrder.top();
        if (node->children.size() > 2) {
            int grp = node->grpID;
            auto temp = node->children;
            while (temp.size() > 2) {
                std::vector<Node*> nodeLeft;
                for (int i = 0; i < temp.size()-1; i+=2) {
                    auto name = this->newInternalNodeId();
                    Node* newNode = new Node (name, 0.0);
                    newNode->children.push_back(temp[i]);
                    newNode->children.push_back(temp[i+1]);
                    newNode->grpID = grp;
                    this->allNodes[name] = newNode;
                    temp[i]->parent = newNode;
                    temp[i+1]->parent = newNode;
                    nodeLeft.push_back(newNode);
                }
                if (temp.size()%2 == 1) nodeLeft.push_back(temp.back());
                temp = nodeLeft;
            }
            assert(temp.size() == 2);
            node->children.clear();
            node->children.push_back(temp[0]);
            node->children.push_back(temp[1]);
            temp[0]->parent = node;
            temp[1]->parent = node;
        }
        // Remove useless intermediate nodes
        else if (node->children.size() == 1 && node->parent != nullptr) { 
            for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                if (node->parent->children[chIdx]->identifier == node->identifier) {
                    node->parent->children[chIdx] = node->children[0];
                    node->children[0]->branchLength += node->branchLength;
                    node->children[0]->parent = node->parent;
                    break;
                }
            }
        }
        else if (node->children.empty() && !node->is_leaf() && !node->seqsIncluded.empty()) {
            std::vector<Node*> children_temp;
            for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                if (node->parent->children[chIdx]->identifier != node->identifier) {
                    children_temp.push_back(node->parent->children[chIdx]);
                }
            }
            node->parent->children = children_temp;
        }
        postOrder.pop();
    }

    std::queue<Node*> q;
    updateLevels(this->root, 1);
    return;
}

void phylogeny::Tree::reroot(bool placement)
{
    int beforeConvert = 0, beforeReroot = 0, afterReroot = 0;
    for (auto n: this->allNodes) beforeConvert = (n.second->level > beforeConvert) ? n.second->level : beforeConvert;
    this->convert2binaryTree();
    for (auto n: this->allNodes) beforeReroot = (n.second->level > beforeReroot) ? n.second->level : beforeReroot;
    // pick any leaf to start BFS
    Node* start = nullptr;
    if (placement) {
        for (auto n : this->allNodes) {
            if (n.second->is_leaf() && n.second->placed) { start = n.second; break; }
        }
    }
    else {
        for (auto n : this->allNodes) {
            if (n.second->is_leaf()) { start = n.second; break; }
        }
    }
    // -----------------------------------------------------
    // BFS function
    // -----------------------------------------------------
    auto bfs = [&](Node* start, std::unordered_map<Node*, Node*>& parentOut) {
        std::queue<Node*> q;
        std::unordered_map<Node*, int> dist;
        parentOut.clear();
        q.push(start);
        dist[start] = 0;
        parentOut[start] = nullptr;
        Node* farthest = start;
        while (!q.empty()) {
            Node* u = q.front(); q.pop();
            // neighbors = parent + children
            std::vector<Node*> neigh = u->children;
            if (u->parent != nullptr) neigh.push_back(u->parent);
            for (Node* v : neigh) {
                if (!dist.count(v)) {
                    dist[v] = dist[u] + 1;
                    parentOut[v] = u;
                    q.push(v);
                    if (placement) {
                        if (dist[v] > dist[farthest] && v->placed) farthest = v;
                    }
                    else {
                        if (dist[v] > dist[farthest]) farthest = v;
                    }

                }
            }
        }
        return farthest;
    };

    // -----------------------------------------------------
    // Find diameter endpoints (A and B)
    // -----------------------------------------------------
    std::unordered_map<Node*,Node*> parentA, parentB;
    Node* A = bfs(start, parentA);
    Node* B = bfs(A, parentB);

    // -----------------------------------------------------
    // Extract the A → B path
    // -----------------------------------------------------
    std::vector<Node*> path;
    for (Node* cur = B; cur != nullptr; cur = parentB[cur])
        path.push_back(cur);
    std::reverse(path.begin(), path.end());
    // -----------------------------------------------------
    // Choose the center node as new root
    // -----------------------------------------------------
    Node* newRoot = path[path.size() / 2];
    if (newRoot->identifier == this->root->identifier) return;

    std::vector<Node*> newRoot_to_oldRoot;
    newRoot_to_oldRoot.push_back(newRoot);
    Node* current = newRoot;
    while (current->parent != nullptr) {
        newRoot_to_oldRoot.push_back(current->parent);
        current = current->parent;
    }
    std::reverse(newRoot_to_oldRoot.begin(), newRoot_to_oldRoot.end());
    for (int i = 0; i < newRoot_to_oldRoot.size()-1; ++i) {
        Node* node = newRoot_to_oldRoot[i]; 
        node->parent = newRoot_to_oldRoot[i+1];
        node->children.erase(std::remove(node->children.begin(), node->children.end(), newRoot_to_oldRoot[i+1]), node->children.end());
        node->branchLength = node->parent->branchLength;
        if (i > 0) node->children.push_back(newRoot_to_oldRoot[i-1]);
    }
    newRoot->children.push_back(newRoot->parent);
    newRoot->parent = nullptr;
    newRoot->branchLength = 0.0;
    updateLevels(newRoot, 1);
    Node* oldRoot = this->root;
    std::string rootName = oldRoot->identifier;
    oldRoot->identifier = newRoot->identifier;
    newRoot->identifier = rootName;
    this->allNodes.erase(rootName);
    this->allNodes.erase(newRoot->identifier);
    this->allNodes[oldRoot->identifier] = oldRoot;
    this->allNodes[newRoot->identifier] = newRoot;
    this->root = newRoot;
    this->convert2binaryTree();
    this->calLeafNum();
    this->calSeqWeight();
    for (auto n: this->allNodes) afterReroot = (n.second->level > afterReroot) ? n.second->level : afterReroot;
    std::cerr << "======== Tree Depth ========\n";
    std::cerr << "Original: " << beforeConvert << '\n';
    std::cerr << "Binary: " << beforeReroot << '\n';
    std::cerr << "Reroot: " << afterReroot << '\n';
}

void phylogeny::Tree::extractResult(Tree* placementT) {
    this->root->seqsIncluded = placementT->root->seqsIncluded;
    if (!placementT->root->msaFreq.empty()) this->root->msaFreq = placementT->root->msaFreq;
    this->root->alnLen = placementT->root->alnLen;
    this->root->alnNum = placementT->root->alnNum;
    this->root->alnWeight = placementT->root->alnWeight;
}

void phylogeny::updateLevels(Node* node, size_t currentLevel) {
    if (!node) return;
    node->level = currentLevel;
    for (Node* child : node->children) {
        updateLevels(child, currentLevel + 1);
    }
}
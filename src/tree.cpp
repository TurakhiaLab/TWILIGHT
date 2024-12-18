#ifndef TREE_HPP
#include "tree.hpp"
#endif

Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
    longestDescendant = nullptr;
}


Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
    longestDescendant = nullptr;
}

size_t Node::getNumLeaves(){
    size_t num_leaves = 0;
    if (children.size() == 0) return num_leaves;
    for (auto ch: children){
        if (ch->is_leaf()) num_leaves += 1;
        else num_leaves += ch->getNumLeaves();
    }
    return num_leaves;
}

size_t Node::getNumNodes(){
    size_t num_nodes = 1;
    if (children.size() == 0) return num_nodes;
    for (auto ch: children){
        num_nodes += ch->getNumNodes();
    }
    return num_nodes;
}

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

Tree::Tree(std::string newickString) {
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
                float len = (branch.size() > 0) ? std::stof(branch) : 1.0;
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
        float len = (branch.size() > 0) ? std::stof(branch) : 1.0;
        if (len == 0) len = 1.0;
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

    treeRoot->branchLength = -1;
    root = treeRoot;

    this->calLeafNum();
    this->calSeqWeight();
}

Tree::Tree(Node* node) {
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
            copyNode->grpID = -1; copyNode->level = current->level - node->level;
            // copyNode->weight = current->weight; copyNode->numLeaves = current->numLeaves;
            this->allNodes[current->identifier] = copyNode;
        }
        if (current->is_leaf()) this->m_numLeaves += 1;
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == grp) s1.push(current->children[i]);     
        }
        // for (auto ch: current->children) {
        //     if (ch->grpID == grp) s1.push(ch);      
        // }
    } 
    this->calLeafNum();
    this->calSeqWeight();
}

Tree::Tree(std::unordered_map<std::string,int>& seqsLen, std::unordered_map<std::string,int>& seqsIdx) {
    Node* treeRoot = new Node("node_1", 1.0);
    this->root = treeRoot;
    this->root->grpID = -1;
    this->allNodes["node_1"] = treeRoot;
    for (auto seq: seqsIdx) {
        // int grpID = seq.second;
        int grpID = this->root->grpID;
        std::string nodeName = seq.first;
        int msaLen = seqsLen[nodeName];
        Node* newNode = new Node(nodeName, this->root, 1.0);
        newNode->grpID = grpID;
        this->allNodes[nodeName] = newNode;
        this->allNodes[nodeName]->msaAln = std::vector<int8_t> (msaLen, 0);
    }
}

Tree::~Tree() {
    for (auto n: this->allNodes) {
        delete n.second;
    }
    this->allNodes.clear();
}


void Tree::calLeafNum() {
    
    // postorder traversal
    std::stack<Node*> s;
    std::stack<Node*> postorder;
    s.push(this->root); 
    Node* current; 
    while (!s.empty()) { 
        current = s.top(); 
        postorder.push(current);
        s.pop(); 
        for (auto ch: current->children) {
            if (ch->grpID == current->grpID) {
                s.push(ch);
            }      
        }
    }
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
    return;
};

void Tree::calSeqWeight() {
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

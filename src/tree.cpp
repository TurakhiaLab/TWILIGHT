#ifndef TREE_HPP
#include "tree.hpp"
#endif

Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
}

void stringSplit (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        if (end_pos >= s.length()) {
            break;
        }
        words.emplace_back(s.substr(start_pos, end_pos-start_pos));
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
        std::string leaf = "";
        std::string branch = "";

        for (auto c: s) {
            if (c == ':') {
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
                float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
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
        leaves.push_back(std::move(leaf));
        numOpen.push_back(no);
        numClose.push_back(nc);
        float len = (branch.size() > 0) ? std::stof(branch) : -1.0;
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
            }
            branchLen[level].pop();
            level++;

            /* Group Id */
            newNode->grpID = -1;

            allNodes[nid] = newNode;
            parentStack.push(newNode);
        }
        Node* leafNode = new Node(leaf, parentStack.top(), branchLen[level].front());
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
}

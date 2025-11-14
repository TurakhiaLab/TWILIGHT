#ifndef PHYLO_HPP
#define PHYLO_HPP

#include <iostream>
#include <vector>
#include <stack>
#include <map>
#include <unordered_map>
#include <unordered_set>


namespace phylogeny {

    using Profile = std::vector<std::vector<float>>;

    struct Node 
    {
        Node(std::string id, float len);
        Node(std::string id, Node* par, float len);
        size_t getNumLeaves();
        size_t getNumNodes();
        bool is_leaf() {return !(identifier.substr(0,4) == "node");}
        void setNumLeaves() {numLeaves = getNumLeaves();};
        
        std::string identifier;
        Node* parent;
        float branchLength;
        size_t level;
        std::vector<Node*> children;
        size_t numLeaves = {0};
        float weight = {0};

        /*Partition*/
        int grpID;
        Node * partitionParent;

        std::vector<int> seqsIncluded; // Indices of sequences included in this node
        Profile msaFreq;
        int alnLen = {0};
        int alnNum = {0};
        float alnWeight = {0};
        int getAlnNum(int currentTask) {
            return (currentTask == 2) ? alnNum : seqsIncluded.size();
        };
        int getAlnLen(int currentTask) {
            return (currentTask == 2 && !msaFreq.empty()) ? msaFreq.size() : alnLen;
        };

        // Unused
        // Node* longestDescendant;
        // void setLongestDescendant(Node* n) {longestDescendant = n;};
        // std::vector<std::string> msa; //use this to store identifier
        // std::vector<int> msaIdx;
        // std::vector<int8_t> msaAln;
        // std::vector<std::vector<float>> msaFreq;
    };

    struct PartitionInfo
    {
        size_t maxPartitionSize;
        size_t minPartitionSize;
        size_t numPartitions;
        std::stack<Node*> partitionsStack;
        std::unordered_map<std::string, std::pair<Node*, size_t>> partitionsRoot;

        PartitionInfo (size_t t_maxPartitionSize, size_t t_minPartitionSize, size_t t_numPartitions):
            maxPartitionSize(t_maxPartitionSize), minPartitionSize(t_minPartitionSize), numPartitions(t_numPartitions) {}
        ~PartitionInfo();

        void partitionTree(Node* root);
        void bipartition(Node* root, Node* edge, Node*& tree1Root, Node*& tree2Root);
    };

    struct Tree
    {
        size_t m_currInternalNode{ 0 };
        size_t m_maxDepth{ 0 };
        size_t m_numLeaves{ 0 };
        float m_meanDepth{ 0 };
        std::string newInternalNodeId() { return "node_" + std::to_string(++m_currInternalNode);}
        Node* root;
        std::unordered_map<std::string, Node*> allNodes;

        void calLeafNum();
        void calSeqWeight();
        void parseNewick(std::string& newickString, bool reroot);
        void showTree();
        void reroot();
        std::string getNewickString();
        Tree* prune(std::unordered_set<std::string>& seqs);
        // Empty tree, used for align alignments
        Tree() {root = nullptr;}
        // Used for whole Tree
        Tree(std::string newickFileName, bool reroot);
        // Used for subtree
        Tree(Node* node);
        // Used for placement without tree
        Tree(std::unordered_set<std::string>& seqNames);

        // Tree(std::unordered_map<std::string,int>& seqsLen, std::unordered_map<std::string,int>& seqsIdx);
        // Tree(std::vector<std::string>& refSeq, std::vector<std::string>& qrySeq);
        ~Tree();
    };

    void pruneTree(Tree*& T, std::unordered_set<std::string>& seqs);
    Tree* constructTreeFromPartitions(Node* root, PartitionInfo* P);
    void updateSubrootInfo(Node*& subroot, Tree* subT, int subtreeIdx);
}


#endif
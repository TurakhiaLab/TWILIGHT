#ifndef PHYLO_HPP
#define PHYLO_HPP

#include <memory> // Required for unique_ptr
#include <string>
#include <iostream>
#include <vector>
#include <stack>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>


namespace phylogeny {

    using Profile = std::vector<std::vector<float>>;

    struct Node {
        // Constructors
        Node(std::string id, float len) 
            : identifier(id), branchLength(len), parent(nullptr), 
              level(0), numLeaves(0), weight(0),
              grpID(-1), partitionParent(nullptr), blockId(0) {}

        // Destructor
        // unique_ptr automatically cleans up children.
        ~Node() = default; 

        // Methods
        size_t getNumLeaves();
        size_t getNumNodes();
        bool is_leaf() const { return !(identifier.substr(0, 4) == "node"); } // Marked const
        void setNumLeaves() { numLeaves = getNumLeaves(); }

        std::unique_ptr<Node> deepCopy( std::unordered_map<std::string, Node*>& nodeTracker, int targetGrp) const;

        // Data Members
        std::string identifier;
        Node* parent; // Raw pointer is correct here (Observer, prevents ownership cycles)
        // THE BIG CHANGE: Owning pointers to children
        std::vector<std::unique_ptr<Node>> children; 

        float branchLength;
        size_t level;
        size_t numLeaves;
        float weight;

        /* Partition */
        int grpID;
        Node* partitionParent; // Raw pointer (Observer)

        uint64_t blockId;

        void collectPostOrder(std::stack<Node*>& postStack);
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
        std::unique_ptr<Node> root;
        std::unordered_map<std::string, Node*> allNodes;

        void calLeafNum();
        void calSeqWeight();
        void parseNewick(std::string& newickString);
        void showTree();
        void reroot(bool placement=false);
        void convert2binaryTree();
        void extractResult(Tree* placementT);
        std::string getNewickString();
        Tree* prune(std::unordered_set<std::string>& seqs);
        // Empty tree, used for align alignments
        Tree() {root = nullptr;}
        // Used for whole Tree
        Tree(const std::string& newickFileName);
        // Used for subtree
        Tree(Node* node, bool reroot);
        // Used for placement without tree
        Tree(std::unordered_set<std::string>& seqNames);

        ~Tree() = default;
    };

    void pruneTree(Tree*& T, std::unordered_set<std::string>& seqs);
    Tree* constructTreeFromPartitions(Node* root, PartitionInfo* P);
    void updateSubrootInfo(Node*& subroot, Tree* subT, int subtreeIdx);
    Tree* getPlacementTree(Tree*);
    void updateLevels(Node* node, size_t currentLevel);
}


#endif
#ifndef MGA_HPP
#define MGA_HPP

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>
#include <zlib.h>

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <list>

#include <boost/program_options.hpp>

#include "phylogeny.hpp"
#include "option.hpp"

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);

// Forward declarations
class Block;
class BlockSet;
class BlockManager;
class SequenceInfo;
class Segment;
class Variation;

using variationList = std::vector<Variation>;


// Helper to check if a CIGAR op consumes Reference
bool consumesRef(char op);

// Helper to check if a CIGAR op consumes Query
bool consumesQry(char op);

inline char complement(char c) {
    switch (c) {
        case 'A': return 'T'; case 'a': return 't';
        case 'C': return 'G'; case 'c': return 'g';
        case 'G': return 'C'; case 'g': return 'c';
        case 'T': return 'A'; case 't': return 'a';
        case 'N': return 'N'; case 'n': return 'n';
        case '-': return '-'; // GAP
        default: return c;
    }
}

std::shared_ptr<Block> stitchBlocks(const std::vector<std::shared_ptr<Block>>& blocks);
void mergeAdjacentBlocks(std::list<std::shared_ptr<Block>>& list, const std::set<int>& cuts);


namespace mga {
    enum ALN_TYPE {
        INIT = 0,
        PRIMARY = 1,
        SECONDARY = 2,
        UNALIGNED = 3,
        REMAINING_ALN = 4
    };

    using IntVec = std::vector<int>;
    using IntPair = std::pair<int, int>;
    using CigarOp = std::pair<int, char>;
    using Cigar = std::vector<CigarOp>;

    using Tree = phylogeny::Tree;
    using Node = phylogeny::Node;
    using NodePair = std::pair<Node*, Node*>;
    using NodePairVec = std::vector<NodePair>;
    using stringPairVec = std::vector<std::pair<std::string, std::string>>;
    
    struct Alignment {

        std::string refName; 
        std::string qryName;

        int identifier;
        IntPair refIdx;  // first: start, second: end (not included)
        IntPair qryIdx;
        bool inverse;
        Cigar CIGAR;
        int alnScore;
        int mis;
        int ins;
        int del;

        // For chaining (Identifying primary alignments)
        bool used;
        float chainScore;
        
        int type;

        struct DuplicationInfo {
            int targetAlnID;      // 重疊的目標 Primary ID
            bool onRef;           // true: 在 Ref 軸重疊; false: 在 Query 軸重疊
            std::pair<int, int> range; // 重疊的區間
        };
        std::vector<DuplicationInfo> dupLinks;


        // For spliting
        bool valid;
        // Function
        void show();

        // Constructor
        Alignment(): used(false), chainScore(0), type(INIT), valid(true) {};
        void setValid2False() { valid = false; }


        std::vector<int> duplications;
        std::vector<int> paralogs;
        
    };

    struct AlnChain {
        int identifier;
        IntVec chainedAln;
        float score;
        AlnChain(int id, float score, IntVec& chain);
    };

    using alnVec = std::vector<Alignment>;
    using chainVec = std::vector<AlnChain>;


    namespace parser {
        Cigar parseCigar(const std::string& cigar);
        alnVec parseMinimap2PAF(const std::string& filename);
    }

    struct SplitResult {
        Alignment head;   
        Alignment middle; 
        Alignment tail;
    };

    struct OverlapInfo {
        int primaryIdx;
        bool onRef; // true = ref, false = qry
        int overlapLen;
        int start;  // absolute coordinate of overlap start
        int end;    // absolute coordinate of overlap end
    };

    // Helper struct for cross-axis processing
    struct CrossOverlap {
        int pIdx;         // The index of the 'Other' Primary (PA, PB)
        int sStart, sEnd; // The overlap range on S's CURRENT axis (The "Other" axis)
        int pStart, pEnd; // The overlap range mapped to P_Main's axis (The "Best" axis)
        bool onRef;       // The axis of pStart/pEnd (same as BestOverlap)
    };

    chainVec getAlignmentChains(alnVec& alignments);
    void identifyPrimaryAlignments(alnVec& alignments, chainVec& chains);
    void resolveGlobalOverlaps(alnVec& alignments);
    std::pair<Alignment, Alignment> singleSplit(const Alignment& parent, int splitPos, bool onRef);
    //vSplitResult splitAlignment(const Alignment& aln, int start, int end, bool onRef);
    // void detectDuplications(alnVec& alignments);
    void fillUnalignedRegions(alnVec& alignments, int refTotalLen, int qryTotalLen);
    bool validateCoverage(const alnVec& alignments, int refTotalLen, int qryTotalLen);
    void collectCutPoints(const std::vector<mga::Alignment>& alignments, std::set<int>& refCuts, std::set<int>& qryCuts);

    std::vector<Alignment> splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts, const std::set<int>& qryCuts);
    std::vector<Alignment> splitAlignmentsByCuts(const std::vector<Alignment>& alignments, const std::set<int>& refCuts, const std::set<int>& qryCuts);

    namespace io
    {
        std::unique_ptr<BlockManager> readSequences(std::string& fileName, Option& option, Tree& tree);
        void writeAlignment(std::string fileName, stringPairVec& seqs, bool compressed, bool append);
    }

    namespace progressive {
        void getProgressivePairs(std::vector<std::pair<NodePair, int>> &alnOrder, std::stack<Node *>& postStack, int grpID, int currentTask);
        void scheduling(Node* root, std::vector<NodePairVec>& levels, int mode);
        void updateNode(NodePairVec& nodes, BlockManager* blockManager);
        void progressiveAlignment(Tree& T, Option& option, std::vector<NodePairVec>& alnPairsPerLevel, BlockManager* blockManager);
        void msaOnSubtree(Tree& T, Option& option, BlockManager* blockManager, int subtree);
        void alignmentKernel(NodePairVec& alnPairs, BlockManager* blockManager, Option& option);
    }
}


class SequenceInfo {

    using ID = std::string;

    private:
        ID sequence_id;
        std::map<int, Segment> info;
        
    public:
        ID getID() const;
        Segment& getSegment(int start_coordinate);
        std::map<int, Segment>& getSegments();

        SequenceInfo(ID id): sequence_id(id) {};
        SequenceInfo() = default;
        ~SequenceInfo() = default;
        bool addSegments(int start, int end, variationList& var);
        bool addSegments(int start, int end);
        bool addSegments(Segment& seg);
};

class Segment {
    private:
        int start_coordinate;
        int end_coordinate;
        bool is_reverse;
        std::vector<Variation> variations;
        std::weak_ptr<Block> prev_blocks_;
        std::weak_ptr<Block> next_blocks_;
    public:
        Segment(int start, int end): start_coordinate(start), end_coordinate(end), is_reverse(false), variations(){};
        Segment(int start, int end, variationList& var): start_coordinate(start), end_coordinate(end), is_reverse(false), variations(var) {};
        Segment() = default;
        ~Segment() = default;
        Segment(const Segment&) = default;
        Segment& operator=(const Segment&) = default;
        std::pair<Segment, Segment> split(int localCut);
        std::vector<Variation>& getVariations() { return variations; }
        int getStart() { return start_coordinate; }
        int getEnd() { return end_coordinate; }
        void setStart(int st) { start_coordinate = st; }
        void setEnd(int en) { end_coordinate = en; }
        bool isReverse() { return is_reverse; }
        std::pair<int,int> getRange() { return {start_coordinate, end_coordinate}; }
        void resetReverse() { is_reverse = false; }
        void reverseComplement(int consLen);
        void setPrevBlock(std::shared_ptr<Block> block) { prev_blocks_ = block; }
        void setNextBlock(std::shared_ptr<Block> block) { next_blocks_ = block; }
        std::weak_ptr<Block> getPrevBlock() { return prev_blocks_; }
        std::weak_ptr<Block> getNextBlock() { return next_blocks_; }
};

class Variation {
    public:
        enum Type { SNV = 1, GAP = 2 };
    private:
        Type type;
        int start;
        int end;
        char alt; 
    public:
        Variation(int pos, char alt): type(SNV), start(pos), end(pos+1), alt(alt) {};
        Variation(int start, int end): type(GAP), start(start), end(end), alt('-') {};
        ~Variation() = default;
        
        Type getType() { return type; }
        int getStart() { return start; }
        int getEnd() { return end; }
        char getAlt() { return alt; }
        std::pair<int,int> getRange() { return {start, end}; }
        
        static Variation createGap(int start, int end) {
            Variation v(start, end);
            return v;
        }
        void reverseComplement(int consLen) {
            int oldStart = start;
            start = consLen - end;
            end = consLen - oldStart;
            if (type == SNV) alt = complement(alt);
        }
        void shift(int offset) {
            start += offset;
            end += offset;
        }
};

class Block: public std::enable_shared_from_this<Block> {
    public:
        using ID = uint64_t;

        Block(ID id, std::string consensus);
        ~Block() = default;
        Block(const Block&) = delete;
        Block& operator=(const Block&) = delete;

        // std::shared_ptr<Block> clone(ID new_id) const {
        //     auto new_block = std::make_shared<Block>(new_id, consensus_sequence_);
        //     for (const auto& seq : sequences_) {
        //         new_block->addSequence(seq);
        //     }
        //     return new_block;
        // }

        // --- Getters ---
        ID getId() const;
        std::string getConsensus() const;
        std::string getConsensusAsFasta(const std::string& name) const;
        std::unordered_map<std::string, SequenceInfo>& getSequences();
        std::vector<std::weak_ptr<Block>> getPrevBlocks() const;
        std::vector<std::weak_ptr<Block>> getNextBlocks() const;

        // --- Modifiers ---
        void addSequence(const SequenceInfo& seq_info);
        void addPrevBlock(std::shared_ptr<Block> block);
        void addNextBlock(std::shared_ptr<Block> block);
        void removePrevBlock(const std::shared_ptr<Block>& block);
        void removeNextBlock(const std::shared_ptr<Block>& block);
        void clearLinkages();

        void print(std::ostream& os = std::cout) const;

        std::pair<std::shared_ptr<Block>, std::shared_ptr<Block>> split(int offset, ID new_id_1, ID new_id_2);


        void addDuplication(ID id) { duplications_.insert(id); }
        void addParalog(ID id) { paralogs_.insert(id); }

        // Check source type (用於判斷 Case 2 vs Case 3)
        bool isFromPrimary() const { return createdFromPrimary_; }
        void setFromPrimary(bool v) { createdFromPrimary_ = v; }
        void reverseComplement() {
            int len = consensus_sequence_.length();
            std::string rc_cons = "";
            for (auto it = consensus_sequence_.rbegin(); it != consensus_sequence_.rend(); ++it) {
                rc_cons += complement(*it);
            }
            consensus_sequence_ = rc_cons;
            for (auto& seq : sequences_) {
                for (auto& seg : seq.second.getSegments()) {
                    seg.second.reverseComplement(len);
                }
            }
        }

    private:
        const ID id_;
        std::string consensus_sequence_;
        std::unordered_map<std::string, SequenceInfo> sequences_;

        // Pangenome graph connections
        std::vector<std::weak_ptr<Block>> prev_blocks_;
        std::vector<std::weak_ptr<Block>> next_blocks_;

        bool createdFromPrimary_ = false; // True if this block was created via Primary Merge
        std::set<ID> duplications_;       // IDs of blocks that are copies of this one
        std::set<ID> paralogs_;           // IDs of blocks that are paralogous to this one
};

class BlockSet {
    public:
        using SetId = std::string;
        // 【新增結構】：字典的節點，記錄原始座標區間與對應的 Block
        struct SegNode {
            int start;
            int end;
            Block::ID blkId;
        };
    
        BlockSet(SetId id);
        ~BlockSet() = default;
    
        BlockSet(const BlockSet&) = delete;
        BlockSet& operator=(const BlockSet&) = delete;
    
        SetId getId() const;

        void selfMapping(Option& option);
    
        std::shared_ptr<Block> createBlock(const std::string& consensus);
        std::shared_ptr<Block> getBlock(Block::ID id);
        std::vector<std::shared_ptr<Block>> getAllBlocks();
        std::vector<Block::ID> getRepresentativeBlocks();
        std::vector<std::shared_ptr<Block>> getRemainingBlocks();
        void clearBlocks();
        std::shared_ptr<Block> addBlock(std::shared_ptr<Block> block);
        bool deleteBlock(Block::ID id);
        size_t getSequenceCount() {return seqs.size(); }
    
        void getRepresentativeAndRemaining(std::vector<std::pair<std::string, std::string>>& representative, std::vector<std::pair<std::string, std::string>>& remaining);
        
        void print(std::ostream& os = std::cout) const;
        std::map<int, uint64_t> splitBlocksByCuts(const std::set<int>& cuts);
        std::pair<uint64_t, uint64_t> splitSingleBlock(int parentID, int localCut);
        std::shared_ptr<Block> mergeTwoBlocks(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, const mga::Cigar& cigar, bool inverse); 
        void updateSegmentLinks(std::shared_ptr<Block> oldBlk, std::shared_ptr<Block> newBlk);
        void addSequence(std::string seqName) {seqs.push_back(seqName); }
        void updateLongestSequence(std::unordered_map<std::string, int>& sequence_lengths);
        void debugValidateSegments(bool verbose);
        void debugValidateLinkages(bool verbose);
        void refine();
    private:
        friend class BlockManager; // Allow BlockManager to modify id_
    
        SetId id_; // No longer const
        std::string longest_sequence_ = "";
        std::vector<std::string> seqs;
        std::atomic<Block::ID> next_block_id_{1};
        std::unordered_map<Block::ID, std::shared_ptr<Block>> blocks_;
        void rebuildDictionary(std::map<int, SegNode>& dict, const std::string& targetSeqName);

        // std::vector<Block::ID> consensus_path_;
};

class BlockManager {
    public:
        BlockManager() = default;
        ~BlockManager() = default;

        BlockManager(const BlockManager&) = delete;
        BlockManager& operator=(const BlockManager&) = delete;

        BlockSet* createBlockSet(BlockSet::SetId id);
        BlockSet* getBlockSet(BlockSet::SetId id) const;
        std::vector<BlockSet*> getAllBlockSets() const;
        bool changeBlockSetId(BlockSet::SetId old_id, BlockSet::SetId new_id);

        void print(std::ostream& os = std::cout) const;

        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, mga::alnVec& alignments);

        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, const std::vector<mga::Alignment>& alignments);
        BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, std::vector<mga::Alignment>& alignments);
        void integrateRemainingBlocks(BlockSet* refSet,  BlockSet* qrySet,  BlockSet* mergedSet,  const std::vector<mga::Alignment>& remainingAlns, int mergedConsensusLen);
        void updateLongestSequences();
        void addSequenceLength(std::string seqName, int seqLen) {
            sequence_lengths[seqName] = seqLen;
        }

    private:
        std::unordered_map<BlockSet::SetId, std::unique_ptr<BlockSet>> block_sets_;
        std::unordered_map<std::string, int> sequence_lengths;

        void applySplits(std::list<std::shared_ptr<Block>>& list, const std::set<int>& cuts);

        std::shared_ptr<Block> mergeBlockPair(BlockSet* resultSet, std::shared_ptr<Block> rBlk, std::shared_ptr<Block> qBlk, const mga::Alignment& subAln);

        std::list<std::shared_ptr<Block>> cloneConsensusPath(BlockSet* set);
        std::list<std::shared_ptr<Block>> prepareBlocks(BlockSet* set, const std::vector<mga::Alignment>& alns, bool isRef);

        // 建立索引 (Pos -> Block)
        struct BlockEntry { int start; int end; std::shared_ptr<Block> blk; };
        std::vector<BlockEntry> indexBlocks(const std::list<std::shared_ptr<Block>>& list) {
            std::vector<BlockEntry> index;
            int pos = 0;
            for (auto& b : list) {
                int len = b->getConsensus().length();
                index.push_back({pos, pos + len, b});
                pos += len;
            }
            return index;
        }

        // 範圍查找
        std::vector<std::shared_ptr<Block>> getBlocksInRange(const std::vector<BlockEntry>& index, std::pair<int, int> range) {
            std::vector<std::shared_ptr<Block>> result;
            int start = std::min(range.first, range.second);
            int end = std::max(range.first, range.second);

            // 簡單線性搜尋 (因為預先切過，應該會精準命中)
            for (const auto& entry : index) {
                if (entry.end > start && entry.start < end) {
                    result.push_back(entry.blk);
                }
            }
            return result;
        }

        int getBlockGlobalStart(const std::vector<BlockEntry>& index, std::shared_ptr<Block> targetBlock) {
            if (!targetBlock) return -1;

            // 簡單的線性搜尋。因為 index 裡的數量通常不多，直接掃描是最快的
            for (const auto& entry : index) {
                // 透過 ID 比對確認是不是同一個 Block
                if (entry.blk->getId() == targetBlock->getId()) {
                    return entry.start;
                }
            }

            // 如果找不到 (理論上不應該發生，除非 targetBlock 不在這個 index 裡)
            std::cerr << "[Warning] Block ID " << targetBlock->getId() << " not found in index.\n";
            return -1; 
        }
};

mga::Cigar adjustCigarWithVariations(mga::Cigar& origCigar, Segment& refSeg, Segment& qrySeg, bool qryInverse, int qryConsLen);
mga::Cigar extractSubCigar(const mga::Cigar& origCigar, int refOffset, int refLen);

#endif

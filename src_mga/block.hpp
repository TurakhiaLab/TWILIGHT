#ifndef BLOCK_HPP
#define BLOCK_HPP

#include <iostream>
#include <memory>
#include <map>
#include <algorithm>
#include <atomic>

#include "type.hpp"
#include "alignment.hpp"

// Forward declaration


std::shared_ptr<Block> extractBlockFromSuper(BlockSet* bSet, std::shared_ptr<Block> superBlock, int start, int end);


class Variant {
    private:
        VariantType type;
        Range range;
        char alt; 
    public:
        Variant(int pos_, char alt_): type(VariantType::SNV), range({pos_, pos_+1}), alt(alt_) {};
        Variant(int start_, int end_): type(VariantType::GAP), range({start_, end_}), alt('-') {};
        ~Variant() = default;
        
        // --- Getter ---
        VariantType getType() { return type; }
        Range getRange() { return range; }
        int getStart() { return range.first; }
        int getEnd() { return range.second; }
        char getAlt() { return alt; }
        
        // --- Setter ---
        void setStart(int st_) { range.first = st_; }
        void setEnd(int en_) { range.second = en_; }
        void setRange(int start_, int end_) { range = {start_, end_}; }
        
        // --- Helper ---
        static Variant createGap(int start, int end) {
            Variant v(start, end);
            return v;
        }
        void reverseComplement(int consLen) {
            int oldStart = range.first;
            range.first = consLen - range.second;
            range.second = consLen - oldStart;
            if (type == VariantType::SNV) alt = complement(alt);
        }
        void shift(int offset) {
            range.first += offset;
            range.second += offset;
        }
};

class Segment {
    private:
        Range range;
        bool is_reverse;
        Variants variations;
        BlockWeakPtr prev_block;
        BlockWeakPtr next_block;

    public:
        Segment(int start, int end): range({start, end}), is_reverse(false), variations(){};
        Segment(int start, int end, Variants& var): range({start, end}), is_reverse(false), variations(var) {};
        Segment(int start, int end, bool reverse): range({start, end}), is_reverse(reverse), variations(){};
        Segment(int start, int end, Variants& var, bool reverse): range({start, end}), is_reverse(reverse), variations(var) {};
        Segment() = default;
        ~Segment() = default;
        Segment(const Segment&) = default;
        Segment& operator=(const Segment&) = default;

        // --- Getter ---
        Variants& getVariants() { return variations; }
        Range getRange() { return range; }
        int getStart() { return range.first; }
        int getEnd() { return range.second; }
        BlockWeakPtr getPrevBlock() { return prev_block; }
        BlockWeakPtr getNextBlock() { return next_block; }
        bool isReverse() { return is_reverse; }

        // --- Setter ---
        void setStart(int st) { range.first = st; }
        void setEnd(int en) { range.second = en; }
        void setReverse(bool reverse) { is_reverse = reverse; }
        void resetReverse() { is_reverse = false; }
        void setPrevBlock(BlockPtr block) { prev_block = block; }
        void setNextBlock(BlockPtr block) { next_block = block; }
        void setPrevBlock(BlockWeakPtr block) { prev_block = block; }
        void setNextBlock(BlockWeakPtr block) { next_block = block; }
        
        // --- Helper ---
        void reverseVariants(int consLen) {
            is_reverse = !is_reverse;
            for (auto& var : variations) {
                var.reverseComplement(consLen);
            }
            std::reverse(variations.begin(), variations.end());
        }
        std::pair<Segment, Segment> split(int localCut);      
};

class Sequence {
    private:
        SequenceID ID;
        Segments segments;

    public:
        Sequence(SequenceID id_): ID(id_) {};
        Sequence() = default;
        ~Sequence() = default;

        // --- Getter ---
        SequenceID getId() const { return ID; };
        Segment& getSegment(int start_coordinate) { return segments[start_coordinate]; };
        std::map<int, Segment>& getSegments() { return segments; };

        // --- Helper ---
        bool addSegment(int start, int end, Variants& var) {
            if (segments.find(start) != segments.end()) return false;
            Segment seg(start, end, var);
            segments[start] = seg;
            return true;
        }
        bool addSegment(int start, int end) {
            if (segments.find(start) != segments.end()) return false;
            Segment seg(start, end);
            segments[start] = seg;
            return true;
        }
        bool addSegment(Segment& seg) {
            if (segments.find(seg.getStart()) != segments.end()) return false;
            segments[seg.getStart()] = seg;
            return true;
        }
        void reverseSegments(int len) {
            for (auto& seg : segments) {
                seg.second.reverseVariants(len);
            }
        }
};



class Block {
    private:
        BlockID ID;
        FamilyID family_ID;
        std::string consensus;
        Sequences sequences;
        BlockWeakPtr prev_block;
        BlockWeakPtr next_block;
        bool distant;

    public:
        Block(BlockID id, std::string consensus): ID(id), family_ID(0), consensus(std::move(consensus)), distant(false) {};
        Block(BlockID id, FamilyID f_id, std::string consensus): ID(id), family_ID(f_id), consensus(std::move(consensus)), distant(false) {};
        ~Block() = default;
        Block(const Block&) = delete;
        Block& operator=(const Block&) = delete;

        // --- Getter ---
        BlockID getId() const { return ID; }
        FamilyID getFamilyId() const { return family_ID; }
        std::string getConsensus() const { return consensus; }
        std::string getConsensusAsFasta(const std::string& name) const { return (">" + name + "\n" + consensus + "\n"); }
        Sequences& getSequences() { return sequences;}
        BlockWeakPtr getPrevBlock() const { return prev_block; };
        BlockWeakPtr getNextBlock() const { return next_block; };
        bool isDistant() const { return distant; }



        // --- Setter ---
        void setId(BlockID id_) { ID = id_; }
        void setFamilyId(FamilyID id_) { family_ID = id_; }
        void setConsensus(std::string& consensus_, bool copy = false) { consensus = (copy) ? consensus_ : std::move(consensus_); }
        void setPrevBlock(BlockPtr block_ptr) { prev_block = block_ptr; }
        void setNextBlock(BlockPtr block_ptr) { next_block = block_ptr; }
        void setSequences(Sequences& sequences_) { sequences.clear(); sequences = sequences_; }
        void setDistant(bool dis) { distant = dis; }

        // --- Helper ---
        void addSequence(const Sequence& seq) {
            sequences[seq.getId()] = seq;
        }
        void reverse() {
            int len = consensus.length();
            consensus = getReverseComplement(consensus);
            for (auto& seq : sequences) {
                seq.second.reverseSegments(len);
            }
        }

        struct ColumnSplitScore {
            double perfect_bonus;
            double id_left;      
            double id_right;     
        };

        std::vector<ColumnSplitScore> calculateSplittingScores(int local_start, int local_end);

        void print(std::ostream& os = std::cout) const;
        bool normalizeStrand();
        void refine();

        // std::pair<std::shared_ptr<Block>, std::shared_ptr<Block>> split(int offset, ID new_id_1, ID new_id_2);
    
};

struct BlockBoundary {
    BlockID leftBlockId;
    BlockID rightBlockId;
    
    int leftConsensusEndPos; 
    
    BoundaryType type = BoundaryType::FLEXIBLE;
    
    int leftFamilyId = 0;  
    int rightFamilyId = 0; 
    
    std::string reason; // For debugging

    bool isHomoLeft() const { return leftFamilyId != 0; }
    bool isHomoRight() const { return rightFamilyId != 0; }
    bool isFamilySeam() const { return leftFamilyId != 0 || rightFamilyId != 0; }

    void print() const {
        std::cout << "[Boundary] " << leftBlockId;
        if (leftFamilyId != 0) std::cout << "(Fam:" << leftFamilyId << ")";
        
        std::cout << " -> " << rightBlockId;
        if (rightFamilyId != 0) std::cout << "(Fam:" << rightFamilyId << ")";
        
        std::cout << " | AbsPos: " << leftConsensusEndPos << " | Type: ";
        switch(type) {
            case BoundaryType::FLEXIBLE: std::cout << "FLEXIBLE"; break;
            case BoundaryType::STRICT: std::cout << "STRICT"; break;
            case BoundaryType::PUSH_RIGHT_ONLY: std::cout << "PUSH_RIGHT_ONLY"; break;
            case BoundaryType::PUSH_LEFT_ONLY: std::cout << "PUSH_LEFT_ONLY"; break;
        }
        std::cout << " | Reason: " << reason << "\n";
    }
};

class BlockSet {

    public:
        struct SegNode; // Forward declaration

    private:
        friend class BlockManager; // Allow BlockManager to modify ID

        BlockSetID ID;
        Blocks blocks;

        std::vector<std::string> sequence_names;
        std::unordered_map<FamilyID, BlockIDs> family_index;
    
        bool is_cached;
        BlockIDs linear_block_cache;
        
        std::atomic<BlockID> next_block_id_{1};
        
        // void rebuildDictionary(std::map<int, SegNode>& dict, const std::string& targetSeqName);
        

        // --- Merger Helper ---
        struct MergeMappingData {
            std::string mergedConsensus;
            std::vector<int> refOldToNew;
            std::vector<int> qryOldToNew;
            Variants newRefGaps;
            Variants newQryGaps;
            CigarString refConsensusChanges;
            CigarString qryConsensusChanges;
        };
        // Helper 1: 走訪 CIGAR、進行 TBB 投票並建構映射表
        MergeMappingData calculateMergedConsensusAndMappings(
            std::shared_ptr<Block> refBlock, 
            std::shared_ptr<Block> qryBlock, 
            const CigarString& cigar, 
            bool inverse, 
            const std::vector<Segment*>& refSegsFlat, 
            const std::vector<Segment*>& qrySegsFlat
        );
        // Helper 2: 平行更新所有 Segment 的 Variations 與座標
        void updateSegmentVariations(
            std::vector<Segment*>& refSegsFlat,
            std::vector<Segment*>& qrySegsFlat,
            const MergeMappingData& mapData,
            bool inverse,
            int refConsLen,
            int qryConsLen
        );
        
        // (選用) 保持版面乾淨的 Debug 列印函數
        void debugValidateBlockMerge(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock); 
    
    public:
        struct SegNode {
            int start;
            int end;
            BlockID blkId;
        };
    
        BlockSet(BlockSetID id) : ID(id), is_cached(false) {};
        ~BlockSet() = default;

        // BlockSet(const BlockSet&) = delete;
        // BlockSet& operator=(const BlockSet&) = delete;
    
        // --- Getter ---
        BlockSetID getId() const { return ID; }
        BlockPtr getBlock(BlockID id) {return (blocks.find(id) == blocks.end()) ? nullptr : blocks[id]; };
        std::vector<std::string> getSequences() {return sequence_names; }
        BlockIDs getLinearizeBlocks();
        BlockIDs getAncestralBlocks();
        StringPairs getRepresentativeConsensus();
        StringPairs getRemainingBlockConsensus();
        
        BlockWeakPtrs getAllBlocks() {
            BlockWeakPtrs all_blocks;
            all_blocks.reserve(blocks.size());
            for (const auto& pair : blocks) all_blocks.push_back(pair.second);
            return all_blocks;
        };
        size_t getSequenceCount() {return sequence_names.size(); }
        
        


        // --- Setter ---
        void invalidateRepCache() { is_cached = false; }
        void setDistantBlocks(Tree& tree, int lookdownDepth = 2);


        // --- Helper ---
        // void selfMapping(Option& option);
        void print(std::ostream& os) const;
        void rebuildAllPointers();
        BlockPtr createBlock(const std::string& consensus) {
            BlockID new_id = next_block_id_++;
            auto new_block = std::make_shared<Block>(new_id, consensus);
            blocks[new_id] = new_block;
            invalidateRepCache();
            return new_block;
        }
        bool deleteBlock(BlockID id) {
            invalidateRepCache();
            return blocks.erase(id) > 0;
        }
        void clearBlocks() {
            blocks.clear();
            invalidateRepCache();
            next_block_id_ = 1; 
        }
        void addSequenceName(std::string seqName) { sequence_names.push_back(seqName); }
        BlockPtr addBlock(BlockPtr oldBlock, FamilyID familyId = 0);
        BlockBoundaries extractBlockBoundaries();
        void normalizeFamilyIDs();
        


        // std::vector<std::shared_ptr<Block>> getRemainingBlocks();
        std::string reconstructSequence(const std::string& seqName);
        BlockPtr concatenateBlocks(BlockID superId);
        std::pair<BlockID, BlockID> splitSingleBlock(int parentID, int localCut);
        
        

        // --- Debug/Validate ---
        void debugValidateSegments(bool verbose = false);
        void debugValidateLinkages(bool verbose = false);
        void debugValidateQuality(bool verbose = false);
        void debugValidateSequences(BlockManager* manager, bool verbose = false);
        void debugValidateLinearizedBlocks(bool verbose = false);
        void debugValidateBlocks(bool verbose = false);


        // --- Merger ---
        BlockPtr mergeTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse);
        void linkTwoBlocks(BlockPtr refBlock, BlockPtr qryBlock, const CigarString& cigar, bool inverse);

        // --- refine ---
        void refineGraph();
        void splitBlocksByLongGaps();
        void extractMicroSegments();
        BlockIDs absorbMicroBlocksIndividual();
        void cleanGappedColumns();
        BlockIDs splitMultiBlocks(BlockID parentID, const std::vector<int>& cuts);
        void reconnectBlocks();

        void realignBlocks(std::string tempDir);
        bool realignBlock(BlockID blkId, std::string tempDir, int iterations=1);

        
        
        // std::shared_ptr<Block> addBlock(std::shared_ptr<Block> oldBlock, BlockID familyId = 0);
        // void updateFamilyIndex(Block::ID blkId, Block::ID famId);
        // std::vector<std::shared_ptr<Block>> getBlocksByFamily(Block::ID famId);


        // const std::vector<int>& getChunkStarts() const { return chunk_starts_; }

        // void splitBlocksByLongGaps();
        // void extractMicroSegments();

        /*
        bool realign(Block::ID blkId, int iterations, std::string tempDir);
        void realignBlocks(std::string tempDir);

        // void refine();
        // void refine_new(BlockSet* refSet, BlockSet* qrySet);
        // void refineFast();
        // void refineBlocks();
        

        
    
        void getRepresentativeAndRemaining(std::vector<std::pair<std::string, std::string>>& representative, std::vector<std::pair<std::string, std::string>>& remaining);
        std::vector<Block::ID> getLinearizeBlock();
        
        void printBlocks(std::ostream& os = std::cout);
        void printBlock(Block::ID blockId, std::ostream& os = std::cout);

        void print(std::ostream& os = std::cout) const;
        std::map<int, uint64_t> splitBlocksByCuts(const std::set<int>& cuts);
        std::vector<Block::ID> splitMultiBlocks(Block::ID parentID, const std::vector<int>& cuts);

        std::shared_ptr<Block> mergeTwoBlocks(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, const mga::Cigar& cigar, bool inverse); 
        void linkTwoBlocks(std::shared_ptr<Block> refBlock, std::shared_ptr<Block> qryBlock, const mga::Cigar& cigar, bool inverse);
        
        void updateSegmentLinks(std::shared_ptr<Block> oldBlk, std::shared_ptr<Block> newBlk);
        void addSequence(std::string seqName) {seqs.push_back(seqName); }
        void updateLongestSequence(std::unordered_map<std::string, int>& sequence_lengths);
        void debugValidateQuality(bool verbose);
        void debugValidateBubble(bool verbose);
        


        
        
        std::vector<Block::ID> absorbMicroBlocks();
        void absorbSingletons();
        

        void realignShallowBlocks(std::string tempDir);
        void realignAllToAll(std::string tempDir);

        void writeGFA(const std::string& outputFileName);
        void writeMAF(const std::string& outputFileName);
        */
        

    
};

class BlockManager {
    private:
        BlockSets blocksets;
        std::unordered_map<std::string, int> sequence_lengths;
        std::unordered_map<std::string, std::string> sequences;

    public:
        BlockManager() = default;
        ~BlockManager() = default;

        BlockManager(const BlockManager&) = delete;
        BlockManager& operator=(const BlockManager&) = delete;

        // --- Getter ---
        BlockSet* getBlockSet(BlockSetID id) const { return (blocksets.find(id) == blocksets.end()) ? nullptr : blocksets.at(id).get(); };
        BlockSets getBlockSets() const { return blocksets; };
        BlockSetPtrs getAllBlockSetPtrs() const {    
            std::vector<BlockSet*> all_sets;
            all_sets.reserve(blocksets.size());
            for (const auto& pair : blocksets) all_sets.push_back(pair.second.get());
            return all_sets;
        }
        

        BlockSet* createBlockSet(BlockSetID id) {
            auto new_set = std::make_unique<BlockSet>(id);
            auto ptr = new_set.get();
            blocksets[id] = std::move(new_set);
            return ptr;
        }
        bool changeBlockSetId(BlockSetID old_id, BlockSetID new_id) {
            if (blocksets.count(new_id) || old_id == new_id) return false;
            auto node_handle = blocksets.extract(old_id);
            if (node_handle.empty()) return false;
            node_handle.key() = new_id;
            node_handle.mapped()->ID = new_id;
            blocksets.insert(std::move(node_handle));
            return true;
        }
        bool removeBlockSet(BlockSetID id) {
            auto it = blocksets.find(id);
            if (it != blocksets.end()) {
                blocksets.erase(it);
                return true;
            }
            return false;
        }

        void print(std::ostream& os = std::cout) const;
        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, AlignmentCollection& alnCollection, int L_min=100);
        BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, BlockBoundaries& refBounds, BlockBoundaries& qryBounds, AlignmentCollection& alnCollection, int L_min=100);
        void orientCircularGenomes(Option& option, std::string refSequenceName = "");

        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, mga::alnVec& alignments);

        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, mga::AlignmentCollection& alnCollection, int L_min);
        // BlockSet* merge(BlockSet* refSet, BlockSet* qrySet, std::vector<mga::Alignment>& alignments);
        // void integrateRemainingBlocks(BlockSet* refSet,  BlockSet* qrySet,  BlockSet* mergedSet,  const std::vector<mga::Alignment>& remainingAlns, int mergedConsensusLen);
        void updateLongestSequences();
        void addSequenceLength(std::string seqName, int seqLen) { sequence_lengths[seqName] = seqLen;}
        int getSequenceLength(std::string seqName) { return sequence_lengths[seqName];}
        void addSequence(std::string& seqName, std::string& seq) {
            if (sequences.find(seqName) != sequences.end()) {
                std::cerr << "ERROR: Sequence " << seqName << " already exists.\n";
                return;
            }
            sequences[seqName] = seq;
        };
        std::string& getSequence(std::string seqName) { return sequences[seqName];}    
};



#endif
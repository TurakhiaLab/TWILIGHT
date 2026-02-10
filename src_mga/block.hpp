#ifndef BLOCK_H
#define BLOCK_H

#include <string>
#include <vector>
#include <set>
#include <memory>
#include <atomic>
#include <iostream>

// TBB Headers
#include <tbb/spin_rw_mutex.h>
#include <tbb/concurrent_hash_map.h>

// Variations
struct Variation {
    int position;
    char reference_base;
    char alternate_base;
};

struct SequenceInfo {
    std::string sequence_id;
    int start_coordinate;
    int end_coordinate;
    std::vector<Variation> variations;
};

class Block {
    public:
        using ID = uint64_t;

        Block(ID id, std::string consensus);
        ~Block() = default;

        Block(const Block&) = delete;
        Block& operator=(const Block&) = delete;

        // --- Getters ---
        ID getId() const;
        std::string getConsensus() const;
        std::string getConsensusAsFasta(const std::string& name) const;
        std::vector<SequenceInfo> getSequences() const;
        std::set<ID> getDuplicationIds() const;

        // --- Modifiers ---
        void addSequence(const SequenceInfo& seq);
        void addDuplicationId(ID dup_id);


        void absorb(const Block& other);
        void print(std::ostream& os = std::cout) const;

    private:
        const ID id_;
        std::string consensus_sequence_;

        std::vector<SequenceInfo> sequences_;
        std::set<ID> duplication_block_ids_;

        mutable tbb::spin_rw_mutex internal_mutex_;
};

// -----------------------------------------------------------
// Wrapper / Manager Class (TBB Version)
// -----------------------------------------------------------

class BlockManager {
    public:
        static BlockManager& instance();

        // 創建 Block
        std::shared_ptr<Block> createBlock(const std::string& consensus);

        // 獲取 Block (Thread-Safe)
        std::shared_ptr<Block> getBlock(Block::ID id);

        // 刪除 Block
        bool deleteBlock(Block::ID id);

        // 合併 Blocks
        bool mergeBlocks(Block::ID to_id, Block::ID from_id);

        // 用於 Debug
        std::vector<Block::ID> getAllBlockIds();

        void printStatus(std::ostream& os = std::cout);

    private:
        BlockManager() = default;
        ~BlockManager() = default;

        // 全局 ID 計數器
        std::atomic<Block::ID> next_id_{1};

        // TBB Concurrent Hash Map
        // Key: Block::ID, Value: shared_ptr<Block>
        // 這個容器本身就是 Thread-safe 的，不需要額外的 manager_mutex
        using BlockMap = tbb::concurrent_hash_map<Block::ID, std::shared_ptr<Block>>;
        BlockMap blocks_;
};

#endif // BLOCK_H
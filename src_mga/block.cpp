
#include "block.hpp"
#include <iomanip>

std::string truncate(const std::string& str, size_t width) {
    if (str.length() > width) {
        return str.substr(0, width - 3) + "...";
    }
    return str;
}

// =======================
// Block Implementation
// =======================

Block::Block(ID id, std::string consensus)
    : id_(id), consensus_sequence_(std::move(consensus)) {}

Block::ID Block::getId() const { return id_; }

std::string Block::getConsensus() const { return consensus_sequence_; }

std::string Block::getConsensusAsFasta(const std::string& name) const {
    return ">" + name + "\n" + consensus_sequence_ + "\n";
}

std::vector<SequenceInfo> Block::getSequences() const {
    // false = Reader Lock (Shared)
    tbb::spin_rw_mutex::scoped_lock lock(internal_mutex_, false); 
    return sequences_;
}

std::set<Block::ID> Block::getDuplicationIds() const {
    // false = Reader Lock (Shared)
    tbb::spin_rw_mutex::scoped_lock lock(internal_mutex_, false);
    return duplication_block_ids_;
}

void Block::addSequence(const SequenceInfo& seq) {
    // true = Writer Lock (Exclusive)
    tbb::spin_rw_mutex::scoped_lock lock(internal_mutex_, true);
    sequences_.push_back(seq);
}

void Block::addDuplicationId(ID dup_id) {
    // true = Writer Lock (Exclusive)
    tbb::spin_rw_mutex::scoped_lock lock(internal_mutex_, true);
    duplication_block_ids_.insert(dup_id);
}

void Block::absorb(const Block& other) {
    // 這裡我們鎖定 "this" (寫鎖)
    tbb::spin_rw_mutex::scoped_lock my_lock(internal_mutex_, true);

    // 讀取 other 的資料 (這會觸發 other 的讀鎖)
    // 為了避免死鎖，最好的方式是不要在持有自己的鎖時去鎖別人
    // 所以我們透過 public method 拿一份 copy 出來
    // (如果 copy 開銷太大，可以使用 std::defer_lock 與 scoped_lock 進行雙重鎖定，但比較複雜)
    
    // 這裡為了邏輯清晰，我們先暫存 other 的數據 (這不會鎖住 this)
    my_lock.release(); // 先釋放自己的鎖
    
    auto other_seqs = other.getSequences();       // 獲取 other 的讀鎖並拷貝
    auto other_dups = other.getDuplicationIds();  // 獲取 other 的讀鎖並拷貝

    my_lock.acquire(internal_mutex_, true); // 重新獲取自己的寫鎖
    
    sequences_.insert(sequences_.end(), other_seqs.begin(), other_seqs.end());
    duplication_block_ids_.insert(other_dups.begin(), other_dups.end());
}

void Block::print(std::ostream& os) const {
    // 1. Acquire Read Lock (shared lock)
    // using "false" means we are NOT writing
    tbb::spin_rw_mutex::scoped_lock lock(internal_mutex_, false);

    const std::string border = "==================================================";
    const std::string subBorder = "--------------------------------------------------";

    os << "\n" << border << "\n";
    os << " BLOCK ID      : " << id_ << "\n";
    os << " CONSENSUS     : " << consensus_sequence_ << "\n";
    os << " LINKED DUPS   : ";
    
    if (duplication_block_ids_.empty()) {
        os << "[None]\n";
    } else {
        os << "{ ";
        for (auto dup : duplication_block_ids_) os << dup << " ";
        os << "}\n";
    }

    os << subBorder << "\n";
    os << " CONTAINED SEQUENCES (" << sequences_.size() << ")\n";
    os << subBorder << "\n";

    if (sequences_.empty()) {
        os << "  (No sequences mapped)\n";
    }

    for (const auto& seq : sequences_) {
        os << "  > ID: " << std::left << std::setw(15) << seq.sequence_id 
           << " Range: [" << seq.start_coordinate << "-" << seq.end_coordinate << "]\n";
        
        if (seq.variations.empty()) {
            os << "    (Perfect Match - No Variations)\n";
        } else {
            os << "    Variations:\n";
            os << "    +----------+-----+-----+\n";
            os << "    | Position | Ref | Alt |\n";
            os << "    +----------+-----+-----+\n";
            for (const auto& var : seq.variations) {
                os << "    | " << std::right << std::setw(8) << var.position << " | " 
                   << std::setw(3) << var.reference_base << " | " 
                   << std::setw(3) << var.alternate_base << " |\n";
            }
            os << "    +----------+-----+-----+\n";
        }
        os << "\n";
    }
    os << border << "\n";
}

// =======================
// BlockManager Implementation
// =======================

BlockManager& BlockManager::instance() {
    static BlockManager instance;
    return instance;
}

std::shared_ptr<Block> BlockManager::createBlock(const std::string& consensus) {
    Block::ID new_id = next_id_++;
    auto new_block = std::make_shared<Block>(new_id, consensus);

    // concurrent_hash_map 的插入操作是 Thread-safe 的
    typename BlockMap::accessor a; // accessor 用於寫入，但在這裡我們主要用 insert
    if (blocks_.insert(a, new_id)) {
        a->second = new_block;
    }
    // accessor 在這裡析構，釋放該 bucket 的鎖
    return new_block;
}

std::shared_ptr<Block> BlockManager::getBlock(Block::ID id) {
    typename BlockMap::const_accessor a; // const_accessor 用於讀取
    if (blocks_.find(a, id)) {
        return a->second;
    }
    return nullptr;
}

bool BlockManager::deleteBlock(Block::ID id) {
    // erase 是 atomic 的
    return blocks_.erase(id);
}

bool BlockManager::mergeBlocks(Block::ID to_id, Block::ID from_id) {
    if (to_id == from_id) return false;

    std::shared_ptr<Block> target_ptr;
    std::shared_ptr<Block> source_ptr;

    // 步驟 1: 獲取 Target 指標
    {
        typename BlockMap::const_accessor a;
        if (blocks_.find(a, to_id)) {
            target_ptr = a->second;
        } else {
            return false; // Target 不存在
        }
    }

    // 步驟 2: 獲取 Source 指標並從 Map 中移除 (原子操作)
    // 我們直接嘗試 erase。如果我們拿到指標後再 erase，中間可能會被別人刪掉。
    // 但 concurrent_hash_map 如果要原子地 "取出並刪除"，需要一點技巧。
    // 比較簡單的方式是先 Find 再 Erase。
    
    {
        typename BlockMap::accessor a;
        if (blocks_.find(a, from_id)) {
            source_ptr = a->second;
            blocks_.erase(a); // 透過 accessor 刪除，保證原子性
        } else {
            return false; // Source 不存在
        }
    }
    
    // 此時 source_ptr 已經從全域 Map 移除，新的線程找不到它了。
    // 但 source_ptr 本身還活著，因為我們持有 shared_ptr。
    
    // 步驟 3: 執行合併
    if (target_ptr && source_ptr) {
        target_ptr->absorb(*source_ptr);
        return true;
    }

    return false;
}

std::vector<Block::ID> BlockManager::getAllBlockIds() {
    std::vector<Block::ID> ids;
    for (auto it = blocks_.begin(); it != blocks_.end(); ++it) {
        ids.push_back(it->first);
    }
    return ids;
}

void BlockManager::printStatus(std::ostream& os) {
    // Note: iterating concurrent_hash_map is thread-safe regarding the container structure,
    // but the data inside blocks might change. 
    
    os << "\n";
    os << "############################################################\n";
    os << "               BLOCK MANAGER SYSTEM STATUS                  \n";
    os << "############################################################\n\n";

    // Table Header
    os << std::left 
       << std::setw(10) << "Block ID" 
       << "| " << std::setw(20) << "Consensus (Preview)" 
       << "| " << std::setw(10) << "Seqs" 
       << "| " << std::setw(10) << "Dups" << "\n";
    
    os << std::string(10, '-') << "+-" 
       << std::string(20, '-') << "+-" 
       << std::string(10, '-') << "+-" 
       << std::string(10, '-') << "\n";

    size_t total_blocks = 0;

    for (auto it = blocks_.begin(); it != blocks_.end(); ++it) {
        // Accessor is not strictly needed just to read pointers if we trust they exist,
        // but here we get the shared_ptr.
        std::shared_ptr<Block> block = it->second;

        // We invoke getters which handle their own locking
        std::string cons = block->getConsensus();
        size_t seq_count = block->getSequences().size();
        size_t dup_count = block->getDuplicationIds().size();

        os << std::left 
           << std::setw(10) << block->getId() 
           << "| " << std::setw(20) << truncate(cons, 18) 
           << "| " << std::setw(10) << seq_count 
           << "| " << std::setw(10) << dup_count << "\n";

        total_blocks++;
    }

    os << "\n> Total Blocks Active: " << total_blocks << "\n";
    os << "############################################################\n";
}
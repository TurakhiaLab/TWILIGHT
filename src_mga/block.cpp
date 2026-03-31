#include "mga.hpp"
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <functional>
#include <unordered_set>


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

std::unordered_map<std::string, SequenceInfo>& Block::getSequences() {
    return sequences_;
}

void Block::addSequence(const SequenceInfo& seq_info) {
    sequences_[seq_info.getID()] = seq_info;
}


void Block::addNextBlock(std::shared_ptr<Block> block) {
    for (const auto& weak_b : next_blocks_) {
        if (auto b = weak_b.lock()) {
            if (b == block) return;
        }
    }
    next_blocks_.push_back(block);
}

void Block::addPrevBlock(std::shared_ptr<Block> block) {
    for (const auto& weak_b : prev_blocks_) {
        if (auto b = weak_b.lock()) {
            if (b == block) return;
        }
    }
    prev_blocks_.push_back(block);
}

void Block::clearLinkages() {
    prev_blocks_.clear();
    next_blocks_.clear();
}

void Block::removePrevBlock(const std::shared_ptr<Block>& block) {
    if (!block) return;
    auto id_to_remove = block->getId();
    prev_blocks_.erase(
        std::remove_if(prev_blocks_.begin(), prev_blocks_.end(),
            [id_to_remove](const std::weak_ptr<Block>& weak_b) {
                if (auto b = weak_b.lock()) {
                    return b->getId() == id_to_remove;
                }
                return true; // Remove expired pointers
            }),
        prev_blocks_.end()
    );
}

void Block::removeNextBlock(const std::shared_ptr<Block>& block) {
    if (!block) return;
    auto id_to_remove = block->getId();
    next_blocks_.erase(
        std::remove_if(next_blocks_.begin(), next_blocks_.end(),
            [id_to_remove](const std::weak_ptr<Block>& weak_b) {
                if (auto b = weak_b.lock()) {
                    return b->getId() == id_to_remove;
                }
                return true; // Remove expired pointers
            }),
        next_blocks_.end()
    );
}


void Block::print(std::ostream& os) const {
    if (sequences_.size() <= 1) return;

    const std::string border = "==================================================";
    const std::string subBorder = "--------------------------------------------------";

    os << "\n" << border << "\n";
    os << " BLOCK ID      : " << id_ << "\n";
    os << " CONSENSUS     : " << consensus_sequence_.size() << "\n";

    auto print_links = [&](const std::string& label, const auto& links) {
        os << label;
        if (links.empty()) {
            os << "[None]\n";
        } else {
            os << "{ ";
            for (const auto& weak_link : links) {
                if (auto link = weak_link.lock()) {
                    os << link->getId() << " ";
                } else {
                    os << "[exp] ";
                }
            }
            os << "}\n";
        }
    };
    
    print_links(" PREV BLOCKS   : ", prev_blocks_);
    print_links(" NEXT BLOCKS   : ", next_blocks_);

    os << subBorder << "\n";
    os << " CONTAINED SEQUENCES (" << sequences_.size() << ")\n";
    os << subBorder << "\n";

    if (sequences_.empty()) {
        os << "  (No sequences mapped)\n";
    }

    // for (const auto& seq : sequences_) {
    //     os << "  > ID: " << std::left << std::setw(15) << seq.sequence_id 
    //        << " Range: [" << seq.start_coordinate << "-" << seq.end_coordinate << "]\n";
    // }
    os << border << "\n";
}


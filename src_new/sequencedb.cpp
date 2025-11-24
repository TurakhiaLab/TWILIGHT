#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include <chrono>
#include <tbb/parallel_for.h>

msa::SequenceDB::SequenceInfo::SequenceInfo(int id_, const std::string &name_, std::string &seq, int subtreeIdx_, float weight_, bool debug, int alnMode)
{
    this->id = id_;
    this->name = name_;
    this->len = seq.length();
    this->subtreeIdx = subtreeIdx_;
    this->storage = 0;
    this->lowQuality = false;
    this->allocate_and_store(seq);
    this->weight = weight_;
    if (debug) {
        if (alnMode == 3) {
            size_t k = 0;
            for (size_t i = 0; i < seq.size(); ++i) {
                if (seq[i] != '-') seq[k++] = seq[i];
            }
            seq.resize(k);
        }
        this->unalignedSeq = seq;
    }
}

msa::SequenceDB::SequenceInfo::~SequenceInfo()
{
    delete[] this->alnStorage[0];
    delete[] this->alnStorage[1];
}


void msa::SequenceDB::SequenceInfo::allocate_and_store(const std::string &seq)
{
    int len = seq.length();
    this->memLen = static_cast<int>(len * this->timesBigger);
    this->alnStorage[0] = new char[this->memLen];
    this->alnStorage[1] = new char[this->memLen];
    for (int i = 0; i < this->memLen; ++i)
    {
        this->alnStorage[0][i] = (i < len) ? seq[i] : '\0';
        this->alnStorage[1][i] = '\0';
    }
    return;
}

void msa::SequenceDB::SequenceInfo::changeStorage()
{
    this->storage = (this->storage == 1) ? 0 : 1;
    return;
}

void msa::SequenceDB::SequenceInfo::memCheck(int len)
{
    if (this->memLen < len)
    {
        int adjustLen = static_cast<int>(len * this->timesBigger);
        char *temp_0 = new char[adjustLen];
        char *temp_1 = new char[adjustLen];
        for (int j = 0; j < adjustLen; ++j)
        {
            temp_0[j] = (j < this->memLen) ? alnStorage[0][j] : 0;
            temp_1[j] = (j < this->memLen) ? alnStorage[1][j] : 0;
        }
        delete[] this->alnStorage[0];
        delete[] this->alnStorage[1];
        this->alnStorage[0] = temp_0;
        this->alnStorage[1] = temp_1;
        this->memLen = adjustLen;
    }
    return;
}

void msa::SequenceDB::addSequence(int id, const std::string &name, std::string &seq, int subtreeIdx, float weight, bool debug, int alnMode)
{
    SequenceInfo *newSeq = new SequenceInfo(id, name, seq, subtreeIdx, weight, debug, alnMode);
    sequences.push_back(newSeq);
    id_map[id] = newSeq;
    name_map[name] = newSeq;
    return;
}

void msa::SequenceDB::debug() {
    auto dbgStart = std::chrono::high_resolution_clock::now();
    int debugNum = 0;
    int alnLen = 0, offset = 0;
    bool theFirst = true;
    for (auto seqInfo: this->sequences) {
        std::string r = "";
        if (!seqInfo->lowQuality) {
            int storage = seqInfo->storage;
            offset = 0;
            while ((seqInfo->alnStorage[storage][offset] >= 'A' && seqInfo->alnStorage[storage][offset] <= 'Z') || 
                   (seqInfo->alnStorage[storage][offset] >= 'a' && seqInfo->alnStorage[storage][offset] <= 'z') || 
                   (seqInfo->alnStorage[storage][offset] == '-')) {
                if (seqInfo->alnStorage[storage][offset] != '-') {
                    r += seqInfo->alnStorage[storage][offset];
                }
                ++offset;
            }
            if (theFirst) {alnLen = offset; theFirst = false;}
            else {
                // if (alnLen != offset) printf("%s: the sequence length (%d) did not match the MSA length(%d)\n", seqInfo->name.c_str(), offset, alnLen);
            }
            if (r != seqInfo->unalignedSeq) {
                printf("%s: after removing the gaps, the alignment did not match the original sequence.\n", seqInfo->name.c_str());    
                // std::cout << r << '\n' << seqInfo->unalignedSeq << '\n';  
            }
            ++debugNum;
        }
    }
    auto dbgEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds dbgTime = dbgEnd - dbgStart;
    std::cerr << "Completed checking " << debugNum << " sequences in " << dbgTime.count() / 1000000 << " ms.\n";
    return;
}

void msa::SequenceDB::storeSubtreeProfile(Tree* subT, char type, int subtreeIdx) {
    this->subtreeAln[subtreeIdx] = alnPath(subT->root->alnLen, 0);
    if (!subT->root->msaFreq.empty()) return;
    int profileSize = (type == 'n') ? 6 : 22;
    subT->root->msaFreq = Profile(subT->root->alnLen, std::vector<float>(profileSize, 0.0));
    for (int sIdx: subT->root->seqsIncluded) {
        int storage = this->id_map[sIdx]->storage;
        float w = this->id_map[sIdx]->weight;
        tbb::parallel_for(tbb::blocked_range<int>(0, subT->root->alnLen), [&](tbb::blocked_range<int> r) {
        for (int j = r.begin(); j < r.end(); ++j) {
            int letterIndex = letterIdx(type, toupper(this->id_map[sIdx]->alnStorage[storage][j]));
            subT->root->msaFreq[j][letterIndex] += 1.0 * w;
        }
        });
    }
    return;
}

void msa::SequenceDB::cleanSubtreeDB() {
    for (auto seq: this->sequences) delete seq;
    this->sequences.clear();
    this->fallback_nodes.clear();
    this->id_map.clear();
    this->name_map.clear();
};

phylogeny::Tree* msa::SequenceDB::getPlacementTree(Tree* T) {
    // Identify nodes on the path from placed sequence to the root
    for (auto node: T->allNodes) {
        if (node.second->is_leaf() && node.second->placed) {
            Node* current = node.second;
            while (current->parent != nullptr) {
                if (current->parent->placed) break;
                current->parent->placed = true;
                current = current->parent;
            }
        }
    }
    // Append backbone aligned sequences to nodes
    std::function<void(Node*, Node*)> addBackboneSeqs = [&](Node* root, Node* node) {
        if (node->is_leaf() && !node->placed) root->seqsIncluded.push_back(this->name_map[node->identifier]->id);
        for (auto& c : node->children)
            if (!c->placed) addBackboneSeqs(root, c);
    };
    for (auto node: T->allNodes) {
        if (node.second->placed) {
            addBackboneSeqs(node.second, node.second);
        }
    }
    // Remove all gap columns
    for (auto node: T->allNodes) {
        if (node.second->placed && !node.second->is_leaf() && !node.second->seqsIncluded.empty()) {
            int length = this->id_map[node.second->seqsIncluded.front()]->len;
            std::vector<uint8_t> allGaps (length, 1);
            tbb::parallel_for(tbb::blocked_range<size_t>(0, length), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t j = range.begin(); j != range.end(); ++j) {
                for (auto sIdx: node.second->seqsIncluded) {
                    if (this->id_map[sIdx]->alnStorage[0][j] != '-') {
                        allGaps[j] = 0;
                        break;
                    }
                }
            }
            });

            tbb::parallel_for(tbb::blocked_range<size_t>(0, node.second->seqsIncluded.size()), [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                int sIdx = node.second->seqsIncluded[i];
                int newIdx = 0;
                for (int j = 0; j < length; ++j) {
                    if (allGaps[j] == 0) {
                        this->id_map[sIdx]->alnStorage[0][newIdx] = this->id_map[sIdx]->alnStorage[0][j];
                        newIdx++;
                    }
                }
                for (int j = newIdx; j < length; ++j) {
                    this->id_map[sIdx]->alnStorage[0][j] = 0;
                }    
                this->id_map[sIdx]->len = newIdx;
            }
            });

            node.second->alnLen = this->id_map[node.second->seqsIncluded.front()]->len;
            node.second->alnNum = node.second->seqsIncluded.size();
            node.second->alnWeight = 0.0;
            for (auto sIdx: node.second->seqsIncluded) node.second->alnWeight += this->id_map[sIdx]->weight;
            // std::cout << node.first << '\t' << node.second->seqsIncluded.size() << '\t' << node.second->alnLen << '/' << length << '\t' << this->id_map[node.second->seqsIncluded.front()]->unalignedSeq.size() << '\t' << node.second->alnNum << '\t' << node.second->alnWeight << '\n';
        }
    }
    // Assign to a new tree
    Tree* placement_T = new Tree();
    for (auto node: T->allNodes) {
        if (node.second->placed) {
            Node* copyNode = new Node (node.second);
            placement_T->allNodes[copyNode->identifier] = copyNode;
        }
    }
    for (auto& node: placement_T->allNodes) {
        for (auto c: T->allNodes[node.first]->children) {
            if (c->placed) node.second->children.push_back(placement_T->allNodes[c->identifier]);
        }
        if (T->allNodes[node.first]->parent != nullptr) {
            auto parentName = T->allNodes[node.first]->parent->identifier;
            node.second->parent = placement_T->allNodes[parentName];
        }
        else {
            node.second->parent = nullptr;
            placement_T->root = node.second;
        }
    }
    
    placement_T->reroot();
    size_t depth = 0, numLeaves = 0;
    for (auto node: placement_T->allNodes) {
        depth = std::max(depth, node.second->level);
        if (node.second->is_leaf()) numLeaves++;
    }
    placement_T->m_numLeaves = numLeaves;
    placement_T->m_maxDepth = depth;
    // placement_T->showTree();

    for (auto seq: this->sequences) if (seq->len <  seq->unalignedSeq.size()) std::cout << seq->name << ':' << seq->len << '\t' << seq->unalignedSeq.size() << '\n';
    
    return placement_T;
}
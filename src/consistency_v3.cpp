#include "msa.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <atomic>
#include <mutex>
#include <numeric>

#include <tbb/parallel_for.h>

std::string msa::getCurrentSequence(const SequenceDB::SequenceInfo* sequence)
{
    return std::string(sequence->alnStorage[sequence->storage], sequence->len);
}

std::vector<int> msa::accurate::gatherClustersFromTree(Node* node, SequenceDB* database, std::size_t targetSize, std::vector<std::vector<int>>& clusters) {
    if (!node) return {};
    
    std::vector<int> local_seqs;
    
    if (node->is_leaf()) {
        local_seqs.push_back(database->name_map[node->identifier]->id);
    } else {
        for (auto ch : node->children) {
            auto child_seqs = gatherClustersFromTree(ch, database, targetSize, clusters);
            local_seqs.insert(local_seqs.end(), child_seqs.begin(), child_seqs.end());
        }
    }
    
    if (local_seqs.size() >= targetSize) {
        clusters.push_back(local_seqs);
        return {};
    }
    return local_seqs;
}

int msa::accurate::DirectPairLibrary::createRecord(int refID, int qryID, int refLen, int qryLen) {
    int pairId = idx(refID, qryID);
    int recordIdx = records.size();
    SparseAlignmentRecord rec;
    rec.forward.assign(refLen, -1);
    rec.backward.assign(qryLen, -1);
    records.push_back(std::move(rec));
    pair_to_record_idx[pairId] = recordIdx;
    return recordIdx;
}

// --- Constructor ---
void msa::accurate::DirectPairLibrary::init (int sequenceCount, std::vector<int>& activeSeqIdx, std::vector<int>& seqLengths) {
    this->totalSequence = sequenceCount;
    totalPairs = sequenceCount * (sequenceCount - 1) / 2;
    pair_to_record_idx.assign(totalPairs, -1);
    
    if (sequenceCount <= DENSE_LIMIT) {
        records.reserve(totalPairs);
    } else {
        size_t rep_pairs = (TARGET_CLUSTERS * (TARGET_CLUSTERS - 1)) / 2;
        size_t avg_cluster_size = std::max(1, sequenceCount / TARGET_CLUSTERS) + 1; 
        size_t internal_pairs_per_cluster = (avg_cluster_size * (avg_cluster_size - 1)) / 2;
        size_t total_cluster_pairs = TARGET_CLUSTERS * internal_pairs_per_cluster;
        size_t expected_total_pairs = rep_pairs + total_cluster_pairs;
        records.reserve(std::min(static_cast<size_t>(totalPairs), expected_total_pairs));
    }

    residueSupport.resize(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i) {
        residueSupport[i].assign(seqLengths[i], 0.0f);
    }

    active_seqs = activeSeqIdx;
    global_to_local.assign(sequenceCount, -1);
    for (int local_id = 0; local_id < totalSequence; ++local_id) {
        global_to_local[activeSeqIdx[local_id]] = local_id;
    }

    seq_to_cluster.assign(totalSequence, -1);
    is_rep.assign(totalSequence, false);
    cluster_members.clear();
    cluster_reps.clear();
    all_reps.clear();
}
        

// --- Getter ---
int msa::accurate::DirectPairLibrary::getSequence() {return totalSequence;}

int msa::accurate::DirectPairLibrary::getPairs() {return totalPairs;}

inline float msa::accurate::DirectPairLibrary::getWeight(int refID, int qryID) const {
    int pairId = idx(refID, qryID);
    if (pairId == -1) return 0.0f;
    int recordIdx = pair_to_record_idx[pairId];
    if (recordIdx == -1) return 0.0f;
    return records[recordIdx].weight;
}

inline int msa::accurate::DirectPairLibrary::getAlignedPos(int refID, int qryID, int refPos) const {
    int pairId = idx(refID, qryID);
    if (pairId == -1) return -1;
    int recordIdx = pair_to_record_idx[pairId];
    if (recordIdx == -1) return -1;
    const auto& rec = records[recordIdx];
    if (refID < qryID) {
        if (refPos >= rec.forward.size()) return -1;
        return rec.forward[refPos];
    } else {
        if (refPos >= rec.backward.size()) return -1;
        return rec.backward[refPos];
    }
}

// --- Modifier ---
void msa::accurate::DirectPairLibrary::addLocalAlignmentResult(int refID, int qryID, AlignmentResult& result) {
    int pairId = idx(refID, qryID);
    int recordIdx = pair_to_record_idx[pairId];
    
    assert(recordIdx != -1 && "Fatal: Alignment record was not pre-allocated!");

    auto& rec = records[recordIdx];
    rec.weight = result.identity;

    if (refID < qryID) { 
        for (const auto& alignedPair : result.alignedPairs) {
            rec.forward[alignedPair.refIndex] = alignedPair.qryIndex;
            rec.backward[alignedPair.qryIndex] = alignedPair.refIndex;
        }
    } else { 
        for (const auto& alignedPair : result.alignedPairs) {
            rec.forward[alignedPair.qryIndex] = alignedPair.refIndex;
            rec.backward[alignedPair.refIndex] = alignedPair.qryIndex;
        }
    }
}



std::shared_ptr<msa::accurate::SubtreeAccurateState> msa::accurate::buildSubtreeAccurateState(SequenceDB* database, Option* option, Tree* tree, int subtreeIdx, Params& params)
{
    auto time0 = std::chrono::high_resolution_clock::now();
    Aligner aligner;

    const auto& sequences = database->sequences;

    std::vector<int> activeSeqIdx;
    activeSeqIdx.reserve(sequences.size());
    for (std::size_t seqIdx = 0; seqIdx < sequences.size(); ++seqIdx) {
        if (!sequences[seqIdx]->lowQuality || option->noFilter) activeSeqIdx.push_back(seqIdx);
    }

    // activeSeqIdx: vectorIdx -> seqIdx
    // currentSequences: seqIdx -> raw sequence

    std::vector<std::string> currentSequences(sequences.size());
    int maxLength = 0;
    int maxSeqID = 0;
    
    std::vector<int> seqLengths(activeSeqIdx.size());
    for (std::size_t idx = 0; idx < activeSeqIdx.size(); ++idx) {
        const std::size_t seqIdx = activeSeqIdx[idx];
        auto seq = getCurrentSequence(sequences[seqIdx]);
        currentSequences[seqIdx] = seq;
        seqLengths[seqIdx] = seq.size();
        maxSeqID = std::max(maxSeqID, static_cast<int>(seqIdx));
    }

    int maxSeqID_add_1 = maxSeqID + 1;
    auto accurateState = std::make_shared<SubtreeAccurateState>(subtreeIdx, maxSeqID_add_1, activeSeqIdx, seqLengths);
    
    auto& ConsistencyLibrary = accurateState->directLib;

    std::vector<std::pair<std::size_t, std::size_t>> pairJobs;
    
    auto addPair = [&](std::size_t a, std::size_t b) {
        if (ConsistencyLibrary.idx(a, b) == -1) return;
        if (a > b) std::swap(a, b);
        if (a != b) pairJobs.push_back({a, b});
    };


    if (activeSeqIdx.size() <= ConsistencyLibrary.DENSE_LIMIT) {
        // Dense Mode
        ConsistencyLibrary.cluster_members.push_back({});
        for (std::size_t seqIdx : activeSeqIdx) {
            ConsistencyLibrary.cluster_members[0].push_back(seqIdx);
            int local_id = ConsistencyLibrary.global_to_local[seqIdx];
            ConsistencyLibrary.all_reps.push_back(local_id);
            ConsistencyLibrary.is_rep[local_id] = true;
        }
        ConsistencyLibrary.cluster_reps.push_back(ConsistencyLibrary.cluster_members[0]);
        
        pairJobs.reserve((activeSeqIdx.size() * (activeSeqIdx.size() - 1)) / 2);
        for (std::size_t i = 0; i < activeSeqIdx.size(); ++i) {
            for (std::size_t j = i + 1; j < activeSeqIdx.size(); ++j) {
                addPair(activeSeqIdx[i], activeSeqIdx[j]);
            }
        }
    } else {
        // Sparse Mode
        std::size_t targetClusters = ConsistencyLibrary.TARGET_CLUSTERS;
        std::size_t targetSize = std::max(static_cast<std::size_t>(1), activeSeqIdx.size() / targetClusters);
        
        std::vector<std::vector<int>> clusters;
        auto leftover = gatherClustersFromTree(tree->root, database, targetSize, clusters);
        if (!leftover.empty()) {
            if (clusters.empty()) clusters.push_back(leftover);
            else clusters.back().insert(clusters.back().end(), leftover.begin(), leftover.end());
        }

        int max_reps_per_cluster = ConsistencyLibrary.REPS_PER_CLUSTER;
        std::vector<int> representatives;
        representatives.reserve(clusters.size() * max_reps_per_cluster);
        
        int cluster_id = 0;
        for (const auto& cl : clusters) {
            std::vector<int> current_members;
            std::vector<int> current_reps;
            
            for (std::size_t s : cl) {
                int local_s = ConsistencyLibrary.global_to_local[s];
                ConsistencyLibrary.seq_to_cluster[local_s] = cluster_id;
                current_members.push_back(s);
            }
            ConsistencyLibrary.cluster_members.push_back(current_members);

            if (cl.size() <= static_cast<std::size_t>(max_reps_per_cluster)) {
                for (std::size_t s : cl) current_reps.push_back(s);
            } else {
                std::size_t chunk_size = cl.size() / max_reps_per_cluster;
                for (int i = 0; i < max_reps_per_cluster; ++i) {
                    std::size_t start = i * chunk_size;
                    std::size_t end = (i == max_reps_per_cluster - 1) ? cl.size() : (start + chunk_size);
                    std::size_t best_seq = cl[start];
                    int max_len = -1;
                    for (std::size_t j = start; j < end; ++j) {
                        int l = currentSequences[cl[j]].size();
                        if (l > max_len) { max_len = l; best_seq = cl[j]; }
                    }
                    current_reps.push_back(best_seq);
                }
            }
            
            for (int r : current_reps) {
                representatives.push_back(r);
                int local_r = ConsistencyLibrary.global_to_local[r];
                ConsistencyLibrary.all_reps.push_back(local_r);
                ConsistencyLibrary.is_rep[local_r] = true;
            }
            ConsistencyLibrary.cluster_reps.push_back(current_reps);
            cluster_id++;
            
            // Cluster: All-to-all
            for (std::size_t i = 0; i < cl.size(); ++i) {
                for (std::size_t j = i + 1; j < cl.size(); ++j) {
                    addPair(cl[i], cl[j]);
                }
            }
        }
        
        // Rep: All-to-all
        for (std::size_t i = 0; i < representatives.size(); ++i) {
            for (std::size_t j = i + 1; j < representatives.size(); ++j) {
                addPair(representatives[i], representatives[j]);
            }
        }
        std::sort(pairJobs.begin(), pairJobs.end());
        pairJobs.erase(std::unique(pairJobs.begin(), pairJobs.end()), pairJobs.end());
    }

    std::atomic<size_t> progress{0};
    std::mutex cout_mutex;
    size_t total = pairJobs.size();

    // Preallocate memory before TBB parallel_for
    ConsistencyLibrary.records.resize(pairJobs.size());
    for (size_t i = 0; i < pairJobs.size(); ++i) {
        
        int refID = pairJobs[i].first;
        int qryID = pairJobs[i].second;
        int pairId = ConsistencyLibrary.idx(refID, qryID);
        // std::cout << refID << " " << qryID << " " << pairId << " " << ConsistencyLibrary.pair_to_record_idx.size() << std::endl;

        ConsistencyLibrary.pair_to_record_idx[pairId] = i;
        ConsistencyLibrary.records[i].forward.assign(currentSequences[refID].size(), -1);
        ConsistencyLibrary.records[i].backward.assign(currentSequences[qryID].size(), -1);
    }

    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, total), [&](const tbb::blocked_range<std::size_t>& range) {
        for (std::size_t pairIdx = range.begin(); pairIdx < range.end(); ++pairIdx) {
            const auto [refIdx, qryIdx] = pairJobs[pairIdx];
            
            auto alignment = aligner.align_affine_local(
                currentSequences[refIdx],
                currentSequences[qryIdx],
                option->type,
                params
            );
            
            ConsistencyLibrary.addLocalAlignmentResult(refIdx, qryIdx, alignment);

            // print progress
            size_t current = ++progress;
            if ((current % 10 == 0 || current == total) && pairIdx == range.begin()) {
                double percent = 100.0 * current / total;
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cout << "1. All-to-all Pairwise Alignment: ["
                          << current << "/" << total
                          << " pairs aligned] ("
                          << std::fixed << std::setprecision(1)
                          << percent << "%)\r"
                          << std::flush;
            }
        }
    });

    std::cout << std::endl;

    auto time1 = std::chrono::high_resolution_clock::now();
    accurateState.get()->directLib.computeResidueSupport();
    auto time2 = std::chrono::high_resolution_clock::now();
    accurateState.get()->directLib.computePairWeights();
    auto time3 = std::chrono::high_resolution_clock::now();

    if (option->printDetail) {
        auto ms = [](auto a, auto b) {
            return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count() / 1000000;
        };
        std::cerr << "Built accurate-mode direct library for subtree " << subtreeIdx
                  << " with " << ConsistencyLibrary.records.size()
                  << " pairs in " << ms(time0, time3) << " ms.\n"
                  << "Time breakdown: \n"
                  << "  1. Build pairwise library: " << ms(time0, time1) << " ms\n"
                  << "  2. Compute residue support: " << ms(time1, time2) << " ms\n"
                  << "  3. Compute pair weights: " << ms(time2, time3) << " ms\n";
    }

    return accurateState;
}

// -------
void msa::accurate::DirectPairLibrary::computeResidueSupport() {
    std::cerr << "2. Compute Residue Support: ";

    for (int u = 0; u < totalSequence; ++u) {
        for (int v = u + 1; v < totalSequence; ++v) {
            int pairId = u * totalSequence - u * (u + 1) / 2 + (v - u - 1);
            int recordIdx = pair_to_record_idx[pairId];
            
            if (recordIdx == -1) continue; 

            const auto& rec = records[recordIdx];
            float w = rec.weight;
            if (w <= 0.0f) continue;

            for (size_t posI = 0; posI < rec.forward.size(); ++posI) {
                int posJ = rec.forward[posI];
                if (posJ != -1) {
                    residueSupport[u][posI] += w; 
                    residueSupport[v][posJ] += w;
                }
            }
        }
    }
    std::cerr << "Done.\n";
}

void msa::accurate::removeColumns(ColumnProvenance& provenance, const IntPairVec& removedColumns)
{
    if (removedColumns.empty()) return;
    ColumnProvenance filtered;
    filtered.reserve(provenance.size());
    int removeIdx = 0;
    int nextRemoveStart = removedColumns[removeIdx].first;
    int nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
    for (int col = 0; col < static_cast<int>(provenance.size()); ++col) {
        while (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveEnd) {
            ++removeIdx;
            if (removeIdx < static_cast<int>(removedColumns.size())) {
                nextRemoveStart = removedColumns[removeIdx].first;
                nextRemoveEnd = nextRemoveStart + removedColumns[removeIdx].second;
            }
        }
        const bool removed = (removeIdx < static_cast<int>(removedColumns.size()) && col >= nextRemoveStart && col < nextRemoveEnd);
        if (!removed) filtered.push_back(std::move(provenance[col]));
    }
    provenance = std::move(filtered);
}

std::vector<std::vector<float>> msa::accurate::buildConsistencyTable(
    const ColumnProvenance& refProvenance,
    const ColumnProvenance& qryProvenance,
    SubtreeAccurateState& accurateState)
{
    auto& directLib = accurateState.directLib;
    int totalSeq = directLib.totalSequence;
    int refLen = refProvenance.size();
    int qryLen = qryProvenance.size();
    
    bool is_dense = (totalSeq <= directLib.DENSE_LIMIT);

    std::vector<std::vector<float>> consistencyTable(refLen, std::vector<float>(qryLen, 0.0f));

    struct QryResInfo {
        int qryCol; int seqB; int localB; int posB; float weight; float S_B;
    };
    std::vector<QryResInfo> qryResList;
    std::vector<std::vector<int>> qryMap(totalSeq);
    std::vector<float> qryColWeightSum(qryLen, 0.0f);
    std::vector<float> refColWeightSum(refLen, 0.0f);
    std::vector<int> qrySeqsList; 

    for (int i = 0; i < totalSeq; ++i) {
        qryMap[i].assign(directLib.residueSupport[i].size(), -1);
    }

    std::vector<bool> inQry(totalSeq, false);
    for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
        for (const auto& qryRes : qryProvenance[qryCol]) {
            int seqB = qryRes.seqId; 
            int localB = directLib.global_to_local[seqB]; 
            if (localB == -1) continue;

            int posB = qryRes.residueIndex;
            int qID = qryResList.size();
            float S_B = directLib.residueSupport[localB][posB]; 
            
            qryResList.push_back({qryCol, seqB, localB, posB, qryRes.weight, S_B});
            qryMap[localB][posB] = qID;
            qryColWeightSum[qryCol] += qryRes.weight;
            
            if (!inQry[localB]) {
                inQry[localB] = true;
                qrySeqsList.push_back(seqB);
            }
        }
    }

    for (int refCol = 0; refCol < refLen; ++refCol) {
        for (const auto& refRes : refProvenance[refCol]) {
            refColWeightSum[refCol] += refRes.weight;
        }
    }

    int qryResCount = qryResList.size();

    // 進入 TBB 平行運算
    tbb::parallel_for( tbb::blocked_range<std::size_t>(0, refLen), [&](const tbb::blocked_range<std::size_t>& range) {
        std::vector<float> num_B(qryResCount, 0.0f);
        std::vector<int> active_qIDs;
        active_qIDs.reserve(totalSeq);
        std::vector<float> col_scores(qryLen, 0.0f);

        for (std::size_t refCol = range.begin(); refCol < range.end(); ++refCol) {
            float denom_ref = refColWeightSum[refCol];
            if (refProvenance[refCol].empty() || denom_ref <= 0.0f) continue;
            
            std::fill(col_scores.begin(), col_scores.end(), 0.0f);

            for (const auto& refRes : refProvenance[refCol]) {
                int seqA = refRes.seqId; 
                int localA = directLib.global_to_local[seqA]; 
                if (localA == -1) continue; 

                int posA = refRes.residueIndex;
                float wA = refRes.weight;
                float S_A = directLib.residueSupport[localA][posA];

                if (S_A <= 0.0f) continue;

                // 提取 A 的身分特徵
                bool isRepA = directLib.is_rep[localA];
                int clusterA = directLib.seq_to_cluster[localA];

                for (int seqB : qrySeqsList) {
                    if (seqA == seqB) continue;
                    int localB = directLib.global_to_local[seqB];
                    if (localB == -1) continue;

                    if (is_dense) {
                        // =========================================================
                        // 🚀 【Dense 模式】: 全部查表
                        // =========================================================
                        int pId = directLib.idx(seqA, seqB);
                        int rId = directLib.pair_to_record_idx[pId];
                        
                        if (rId != -1) {
                            const auto& rec = directLib.records[rId];
                            if (localA < localB) {
                                if (posA < rec.extendedForward.size()) {
                                    for (const auto& ext : rec.extendedForward[posA]) {
                                        if (qryMap[localB][ext.targetPos] != -1) {
                                            int qID = qryMap[localB][ext.targetPos];
                                            if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                            num_B[qID] += ext.weight; 
                                        }
                                    }
                                }
                            } else {
                                if (posA < rec.extendedBackward.size()) {
                                    for (const auto& ext : rec.extendedBackward[posA]) {
                                        if (qryMap[localB][ext.targetPos] != -1) {
                                            int qID = qryMap[localB][ext.targetPos];
                                            if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                            num_B[qID] += ext.weight;
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        // =========================================================
                        // 🛣️ 【Sparse 模式】: 完美落實你的 5 條規則！
                        // =========================================================
                        bool isRepB = directLib.is_rep[localB];
                        int clusterB = directLib.seq_to_cluster[localB];

                        // 【規則 1 & 5】：兩個 seq 都是 Rep
                        if (isRepA && isRepB) {
                            // (1) 從 Dense Mode Table 讀取 Rep-Rep 的預計算分數
                            int pId = directLib.idx(seqA, seqB);
                            int rId = directLib.pair_to_record_idx[pId];
                            if (rId != -1) {
                                const auto& rec = directLib.records[rId];
                                if (localA < localB) {
                                    if (posA < rec.extendedForward.size()) {
                                        for (const auto& ext : rec.extendedForward[posA]) {
                                            if (qryMap[localB][ext.targetPos] != -1) {
                                                int qID = qryMap[localB][ext.targetPos];
                                                if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                                num_B[qID] += ext.weight;
                                            }
                                        }
                                    }
                                } else {
                                    if (posA < rec.extendedBackward.size()) {
                                        for (const auto& ext : rec.extendedBackward[posA]) {
                                            if (qryMap[localB][ext.targetPos] != -1) {
                                                int qID = qryMap[localB][ext.targetPos];
                                                if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                                num_B[qID] += ext.weight;
                                            }
                                        }
                                    }
                                }
                            }

                            // 【規則 5 延伸】：如果在同個 Cluster，補上 Cluster 內部非 Rep Member 的貢獻
                            if (clusterA == clusterB) {
                                for (int seqC : directLib.cluster_members[clusterA]) {
                                    if (seqC == seqA || seqC == seqB) continue;
                                    if (directLib.is_rep[directLib.global_to_local[seqC]]) continue; // Rep 已算過
                                    
                                    int posC = directLib.getAlignedPos(seqA, seqC, posA);
                                    if (posC == -1) continue;
                                    int posB = directLib.getAlignedPos(seqC, seqB, posC);
                                    if (posB != -1 && qryMap[localB][posB] != -1) {
                                        int qID = qryMap[localB][posB];
                                        if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                        num_B[qID] += std::min(directLib.getWeight(seqA, seqC), directLib.getWeight(seqC, seqB));
                                    }
                                }
                            }
                        } 
                        // 【規則 2, 3, 4】：不全為 Rep 的情況
                        else {
                            // (A) 先補上 Direct Weight (只有同群組才有)
                            if (clusterA == clusterB) {
                                int pos_direct = directLib.getAlignedPos(seqA, seqB, posA);
                                if (pos_direct != -1 && qryMap[localB][pos_direct] != -1) {
                                    int qID = qryMap[localB][pos_direct];
                                    if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                    num_B[qID] += directLib.getWeight(seqA, seqB);
                                }
                            }

                            // (B) 找出合法的第三方橋樑 C
                            const std::vector<int>* bridge_candidates = nullptr;

                            if (clusterA == clusterB) {
                                // 【規則 4】：同群組的 Member，Search Space 限定在 Cluster 內
                                bridge_candidates = &directLib.cluster_members[clusterA];
                            } else {
                                // 不同群組
                                if (isRepA && !isRepB) {
                                    // 【規則 3】：A是Rep, B是Member，只看 B 的 Cluster Reps
                                    bridge_candidates = &directLib.cluster_reps[clusterB];
                                } else if (!isRepA && isRepB) {
                                    // 【規則 3】：A是Member, B是Rep，只看 A 的 Cluster Reps
                                    bridge_candidates = &directLib.cluster_reps[clusterA];
                                } else {
                                    // 【規則 2】：兩個都是 Member，連橋樑都沒有 -> 直接 skip 看都不用看！
                                    bridge_candidates = nullptr; 
                                }
                            }

                            // (C) 執行 On-the-fly 橋樑計算
                            if (bridge_candidates) {
                                for (int seqC : *bridge_candidates) {
                                    if (seqC == seqA || seqC == seqB) continue;
                                    int posC = directLib.getAlignedPos(seqA, seqC, posA);
                                    if (posC == -1) continue;
                                    int posB = directLib.getAlignedPos(seqC, seqB, posC);
                                    
                                    if (posB != -1 && qryMap[localB][posB] != -1) {
                                        int qID = qryMap[localB][posB];
                                        if (num_B[qID] == 0.0f) active_qIDs.push_back(qID);
                                        num_B[qID] += std::min(directLib.getWeight(seqA, seqC), directLib.getWeight(seqC, seqB));
                                    }
                                }
                            }
                        }
                    }
                }

                for (int qID : active_qIDs) {
                    float num = num_B[qID];
                    num_B[qID] = 0.0f; 

                    const auto& qInfo = qryResList[qID];
                    float S_B = qInfo.S_B;

                    float denom = std::min(S_A, S_B);
                    float tcs = (denom > 0.0f) ? std::clamp(num / denom, 0.0f, 1.0f) : 0.0f;
                    col_scores[qInfo.qryCol] += wA * qInfo.weight * tcs;
                }
                active_qIDs.clear();
            }
            
            for (int qryCol = 0; qryCol < qryLen; ++qryCol) {
                float denom = denom_ref * qryColWeightSum[qryCol];
                if (denom > 0.0f) {
                    consistencyTable[refCol][qryCol] = col_scores[qryCol] / denom;
                }
            }
        }
    });

    return consistencyTable;
}

void msa::accurate::DirectPairLibrary::computePairWeights() {
    bool is_dense = (totalSequence <= DENSE_LIMIT);
    
    std::vector<int> target_locals;
    if (is_dense) {
        std::cerr << "3. Compute Pair Weights (Dense Mode - All Pairs):\n";
        target_locals.resize(totalSequence);
        std::iota(target_locals.begin(), target_locals.end(), 0);
    } else {
        std::cerr << "3. Compute Pair Weights (Sparse Mode - Rep Highway):\n";
        target_locals = all_reps;
        
        std::sort(target_locals.begin(), target_locals.end());
        target_locals.erase(std::unique(target_locals.begin(), target_locals.end()), target_locals.end());
    }

    std::atomic<size_t> progress{0};
    std::mutex cout_mutex;
    size_t total_targets = target_locals.size();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, total_targets), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t i = range.begin(); i < range.end(); ++i) {
            int u = target_locals[i];
            for (size_t j = i + 1; j < total_targets; ++j) {
                int v = target_locals[j];
                
                // 現在有了修復 1，u 絕對小於 v，這個數學公式 100% 安全
                int pairId = u * totalSequence - u * (u + 1) / 2 + (v - u - 1);
                int recIdx = pair_to_record_idx[pairId];
                if (recIdx == -1) continue;

                auto& rec = records[recIdx];
                int lenU = rec.forward.size();
                int lenV = rec.backward.size();

                rec.extendedForward.assign(lenU, std::vector<ExtendedWeight>());
                rec.extendedBackward.assign(lenV, std::vector<ExtendedWeight>());

                std::vector<float> num_V(lenV, 0.0f);
                std::vector<int> active_V;
                active_V.reserve(lenV);

                float w_uv = rec.weight;

                for (size_t posU = 0; posU < lenU; ++posU) {
                    int directV = rec.forward[posU];
                    if (directV != -1) {
                        num_V[directV] = w_uv;
                        active_V.push_back(directV);
                    }

                    for (int c : target_locals) {
                        if (c == u || c == v) continue;

                        auto getPos = [&](int a, int b, int posA) {
                            if (a > b) {
                                int pId = b * totalSequence - b * (b + 1) / 2 + (a - b - 1);
                                int rId = pair_to_record_idx[pId];
                                if (rId == -1) return -1;
                                // 🔥 【修復 2】：加上強制的陣列邊界防護，遇到異常直接略過
                                if (posA < 0 || posA >= records[rId].backward.size()) return -1;
                                return records[rId].backward[posA];
                            } else {
                                int pId = a * totalSequence - a * (a + 1) / 2 + (b - a - 1);
                                int rId = pair_to_record_idx[pId];
                                if (rId == -1) return -1;
                                // 🔥 【修復 2】：加上強制的陣列邊界防護
                                if (posA < 0 || posA >= records[rId].forward.size()) return -1;
                                return records[rId].forward[posA];
                            }
                        };

                        int posC = getPos(u, c, posU);
                        if (posC == -1) continue;
                        
                        int posV = getPos(c, v, posC);
                        if (posV == -1) continue;
                        
                        // 雙重保險：確保算出來的 posV 不會超過 num_V 邊界
                        if (posV < 0 || posV >= lenV) continue;

                        auto getW = [&](int a, int b) {
                            int _a = std::min(a, b); int _b = std::max(a, b);
                            int pId = _a * totalSequence - _a * (_a + 1) / 2 + (_b - _a - 1);
                            int rId = pair_to_record_idx[pId];
                            return (rId == -1) ? 0.0f : records[rId].weight;
                        };
                        
                        float w_uc = getW(u, c);
                        float w_cv = getW(c, v);
                        if (w_uc > 0.0f && w_cv > 0.0f) {
                            if (num_V[posV] == 0.0f && directV != posV) active_V.push_back(posV);
                            num_V[posV] += std::min(w_uc, w_cv);
                        }
                    }

                    for (int posV : active_V) {
                        float w = num_V[posV];
                        if (w > 0.0f) {
                            rec.extendedForward[posU].push_back({posV, w});
                            rec.extendedBackward[posV].push_back({static_cast<int32_t>(posU), w});
                        }
                        num_V[posV] = 0.0f;
                    }
                    active_V.clear();
                }
            }
            
            size_t current = ++progress;
            if ((current % 10 == 0 || current == total_targets) && i == range.begin()) {
                double percent = 100.0 * current / total_targets;
                std::lock_guard<std::mutex> lock(cout_mutex);
                std::cerr << "  [" << current << "/" << total_targets << " seqs computed] ("
                          << std::fixed << std::setprecision(1) << percent << "%)\r"
                          << std::flush;
            }
        }
    });
    std::cerr << "\n";
}
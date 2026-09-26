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
    
    records.clear(); 

    residueSupport.resize(sequenceCount);
    residueSupportRep.resize(sequenceCount);
    for (int i = 0; i < sequenceCount; ++i) {
        residueSupport[i].assign(seqLengths[i], 0.0f);
        residueSupportRep[i].assign(seqLengths[i], 0.0f);
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

    std::vector<std::string> currentSequences(sequences.size());
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

    // DENSE Mde or Not
    bool is_dense = (activeSeqIdx.size() <= ConsistencyLibrary.DENSE_LIMIT);

    // =========================================================================
    // All-to-all Pairwise Alignment
    // =========================================================================
    auto alignJobs = [&](const std::vector<std::pair<std::size_t, std::size_t>>& jobs, const std::string& phaseName) {
        if (jobs.empty()) return;
        
        size_t start_idx = ConsistencyLibrary.records.size();
        ConsistencyLibrary.records.resize(start_idx + jobs.size());
        
        for (size_t i = 0; i < jobs.size(); ++i) {
            int refID = jobs[i].first;
            int qryID = jobs[i].second;
            int pairId = ConsistencyLibrary.idx(refID, qryID);

            ConsistencyLibrary.pair_to_record_idx[pairId] = start_idx + i;
            ConsistencyLibrary.records[start_idx + i].forward.assign(currentSequences[refID].size(), -1);
            ConsistencyLibrary.records[start_idx + i].backward.assign(currentSequences[qryID].size(), -1);
        }

        std::atomic<size_t> progress{0};
        std::mutex cout_mutex;
        size_t total = jobs.size();

        tbb::parallel_for( tbb::blocked_range<std::size_t>(0, total), [&](const tbb::blocked_range<std::size_t>& range) {
            Aligner aligner;
            for (std::size_t pairIdx = range.begin(); pairIdx < range.end(); ++pairIdx) {
                const auto [refIdx, qryIdx] = jobs[pairIdx];
                
                auto alignment = aligner.align_affine_local(
                    currentSequences[refIdx], currentSequences[qryIdx], option->type, params
                );

                // auto alignment = aligner.align_linear_local(
                //     currentSequences[refIdx], currentSequences[qryIdx], option->type, params
                // );
                
                
                ConsistencyLibrary.addLocalAlignmentResult(refIdx, qryIdx, alignment);

                size_t current = ++progress;
                if ((current % 10 == 0 || current == total) && pairIdx == range.begin()) {
                    double percent = 100.0 * current / total;
                    std::lock_guard<std::mutex> lock(cout_mutex);
                    std::cout << phaseName << " [" << current << "/" << total << "] ("
                              << std::fixed << std::setprecision(1) << percent << "%)\r" << std::flush;
                }
            }
        });
        std::cout << std::endl;
    };

    if (is_dense) {
        // [Dense Mode]: All-to-all
        size_t exact_pairs = (activeSeqIdx.size() * (activeSeqIdx.size() - 1)) / 2;
        ConsistencyLibrary.records.reserve(exact_pairs);

        ConsistencyLibrary.cluster_members.push_back({});
        for (std::size_t seqIdx : activeSeqIdx) {
            ConsistencyLibrary.cluster_members[0].push_back(seqIdx);
            int local_id = ConsistencyLibrary.global_to_local[seqIdx];
            ConsistencyLibrary.all_reps.push_back(local_id);
            ConsistencyLibrary.is_rep[local_id] = true;
        }
        ConsistencyLibrary.cluster_reps.push_back(ConsistencyLibrary.cluster_members[0]);
        
        std::vector<std::pair<std::size_t, std::size_t>> denseJobs;
        denseJobs.reserve((activeSeqIdx.size() * (activeSeqIdx.size() - 1)) / 2);
        for (std::size_t i = 0; i < activeSeqIdx.size(); ++i) {
            for (std::size_t j = i + 1; j < activeSeqIdx.size(); ++j) {
                std::size_t u = activeSeqIdx[i], v = activeSeqIdx[j];
                if (u > v) std::swap(u, v);
                denseJobs.push_back({u, v});
            }
        }
        alignJobs(denseJobs, "1. [Dense  Mode] All-to-all Alignment:");
    } 
    else {
        // [Sparse Mode]: 2-Phase
        size_t est_reps_pairs = (ConsistencyLibrary.TARGET_TOTAL_REPS * ConsistencyLibrary.TARGET_TOTAL_REPS) / 2;
        size_t est_num_clusters = std::max((size_t)1, activeSeqIdx.size() / ConsistencyLibrary.MAX_CLUSTER_SIZE);
        size_t est_intra_pairs_per_cluster = (ConsistencyLibrary.MAX_CLUSTER_SIZE * ConsistencyLibrary.MAX_CLUSTER_SIZE) / 2;
        size_t total_est_pairs = est_reps_pairs + (est_num_clusters * est_intra_pairs_per_cluster);
        
        ConsistencyLibrary.records.reserve(total_est_pairs * 1.1);

        std::vector<std::vector<int>> clusters;
        auto leftover = gatherClustersFromTree(tree->root, database, ConsistencyLibrary.MAX_CLUSTER_SIZE, clusters);
        if (!leftover.empty()) {
            if (clusters.empty()) clusters.push_back(leftover);
            else clusters.back().insert(clusters.back().end(), leftover.begin(), leftover.end());
        }

        // =========================================================================
        // 🔥 動態配額分配 (Dynamic Quota Allocation)
        // =========================================================================
        int num_clusters = clusters.size();
        std::vector<int> allocated_reps(num_clusters, 0);

        if (num_clusters >= ConsistencyLibrary.TARGET_TOTAL_REPS) {
            // 情境 1：Cluster 數量爆炸 (>= 500)
            // 為了保證多樣性，每個 Cluster 強制給 1 個代表，允許總數超過 500
            for (int i = 0; i < num_clusters; ++i) {
                allocated_reps[i] = 1;
            }
        } else {
            int base_quota = ConsistencyLibrary.TARGET_TOTAL_REPS / num_clusters;
            int remainder_quota = ConsistencyLibrary.TARGET_TOTAL_REPS % num_clusters;
            std::vector<std::pair<int, int>> size_priority;
            for (int i = 0; i < num_clusters; ++i) {
                size_priority.push_back({clusters[i].size(), i});
            }
            std::sort(size_priority.rbegin(), size_priority.rend());

            for (int i = 0; i < num_clusters; ++i) {
                int c_idx = size_priority[i].second;
                int c_size = clusters[c_idx].size();
                int expected_quota = base_quota + (i < remainder_quota ? 1 : 0);
                allocated_reps[c_idx] = std::min(expected_quota, c_size);
            }
        }
        // =========================================================================
        // --- Phase 1: Intra-cluster Alignment ---
        std::vector<std::pair<std::size_t, std::size_t>> intraJobs;
        for (int cluster_id = 0; cluster_id < num_clusters; ++cluster_id) {
            const auto& cl = clusters[cluster_id];
            ConsistencyLibrary.cluster_members.push_back(cl);
            
            for (std::size_t s : cl) {
                ConsistencyLibrary.seq_to_cluster[ConsistencyLibrary.global_to_local[s]] = cluster_id;
            }

            // Cluster Size <= reps_per_cluster
            if (cl.size() <= allocated_reps[cluster_id]) {
                continue; 
            }

            // Intra-cluster All-to-all
            for (size_t i = 0; i < cl.size(); ++i) {
                for (size_t j = i + 1; j < cl.size(); ++j) {
                    size_t u = cl[i], v = cl[j];
                    if (u > v) std::swap(u, v);
                    intraJobs.push_back({u, v});
                }
            }
        }
        std::sort(intraJobs.begin(), intraJobs.end());
        intraJobs.erase(std::unique(intraJobs.begin(), intraJobs.end()), intraJobs.end());

        alignJobs(intraJobs, "1-1. Intra-cluster Alignment:");

        // --- Rep Selection ---
        for (int cluster_id = 0; cluster_id < num_clusters; ++cluster_id) {
            const auto& cl = ConsistencyLibrary.cluster_members[cluster_id];
            std::vector<std::pair<float, int>> scores; 
            
            for (int seqA : cl) {
                float total_score = 0.0f;
                for (int seqB : cl) {
                    if (seqA == seqB) continue;
                    int pId = ConsistencyLibrary.idx(seqA, seqB);
                    int rId = ConsistencyLibrary.pair_to_record_idx[pId];
                    if (rId != -1) {
                        int aln_len = 0;
                        const auto& rec = ConsistencyLibrary.records[rId];
                        const auto& map_array = (seqA < seqB) ? rec.forward : rec.backward;
                        
                        for (int pos : map_array) {
                            if (pos != -1) aln_len++;
                        }
                        // Criteria: local alignment length * identity
                        total_score += aln_len * rec.weight; 
                    }
                }
                float avg_score = (cl.size() > 1) ? total_score / (cl.size() - 1) : 0.0f;
                scores.push_back({avg_score, seqA});
            }
            
            std::sort(scores.rbegin(), scores.rend());
            
            std::vector<int> current_reps;
            int reps_to_pick = allocated_reps[cluster_id];
            for (int i = 0; i < reps_to_pick; ++i) {
                int rep_id = scores[i].second;
                current_reps.push_back(rep_id);
                ConsistencyLibrary.all_reps.push_back(ConsistencyLibrary.global_to_local[rep_id]);
                ConsistencyLibrary.is_rep[ConsistencyLibrary.global_to_local[rep_id]] = true;
            }
            ConsistencyLibrary.cluster_reps.push_back(current_reps);
        }

        // =========================================================================
        // 🔥 Sparse Mode Debug Message
        // =========================================================================
        if (option->printDetail) {
            int total_seqs_in_sparse = 0;
            int total_reps_picked = 0;

            std::cerr << "\n--- [Sparse Mode Allocation Info] ---\n";
            std::cerr << "Total Clusters      : " << num_clusters << "\n";
            
            for (int i = 0; i < num_clusters; ++i) {
                int c_size = ConsistencyLibrary.cluster_members[i].size();
                int r_picked = ConsistencyLibrary.cluster_reps[i].size();
                
                total_seqs_in_sparse += c_size;
                total_reps_picked += r_picked;

                std::cerr << "  - Cluster " << std::setw(3) << i 
                          << " | Size: " << std::setw(4) << c_size 
                          << " | Reps picked: " << std::setw(4) << r_picked << "\n";
            }
            
            std::cerr << "-------------------------------------\n";
            std::cerr << "Total Sequences       : " << total_seqs_in_sparse << "\n";
            std::cerr << "Total Target Reps (T) : " << total_reps_picked 
                      << " (Target Limit: " << ConsistencyLibrary.TARGET_TOTAL_REPS << ")\n\n";
        }
        // =========================================================================

        // --- Phase 2: Inter-Rep Alignment ---
        std::vector<std::pair<std::size_t, std::size_t>> repJobs;
        for (size_t i = 0; i < ConsistencyLibrary.all_reps.size(); ++i) {
            for (size_t j = i + 1; j < ConsistencyLibrary.all_reps.size(); ++j) {
                int u = activeSeqIdx[ConsistencyLibrary.all_reps[i]];
                int v = activeSeqIdx[ConsistencyLibrary.all_reps[j]];
                if (u > v) std::swap(u, v);
                
                // If already compared in a cluster, then skip
                if (ConsistencyLibrary.pair_to_record_idx[ConsistencyLibrary.idx(u, v)] == -1) {
                    repJobs.push_back({u, v});
                }
            }
        }
        std::sort(repJobs.begin(), repJobs.end());
        repJobs.erase(std::unique(repJobs.begin(), repJobs.end()), repJobs.end());

        alignJobs(repJobs, "1-2. Inter-Rep Alignment:");
    }

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
        bool isRepU = is_rep[u]; // 提取身份
        for (int v = u + 1; v < totalSequence; ++v) {
            int pairId = u * totalSequence - u * (u + 1) / 2 + (v - u - 1);
            int recordIdx = pair_to_record_idx[pairId];
            
            if (recordIdx == -1) continue; 

            const auto& rec = records[recordIdx];
            float w = rec.weight;
            if (w <= 0.0f) continue;

            bool isRepV = is_rep[v]; // 提取身份
            bool bothReps = (isRepU && isRepV); // 判斷是否皆為 Rep

            for (size_t posI = 0; posI < rec.forward.size(); ++posI) {
                int posJ = rec.forward[posI];
                if (posJ != -1) {
                    // 全域 Support (維持原樣)
                    residueSupport[u][posI] += w; 
                    residueSupport[v][posJ] += w;

                    // 🔥 純 Rep Support (只有兩人都是 Rep 時才記錄)
                    if (bothReps) {
                        residueSupportRep[u][posI] += w;
                        residueSupportRep[v][posJ] += w;
                    }
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
        int qryCol; int seqB; int localB; int posB; float weight; float S_B; float S_B_rep;
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
            float S_B_rep = directLib.residueSupportRep[localB][posB];
            
            qryResList.push_back({qryCol, seqB, localB, posB, qryRes.weight, S_B, S_B_rep});
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
                    int localB = qInfo.localB;
                    float S_B = qInfo.S_B;
                    float S_B_rep = qInfo.S_B_rep; // 🔥 取出 S_B_rep

                    // 預設使用全域的 min
                    float denom = std::min(S_A, S_B);

                    if (!is_dense) {
                        bool isRepB = directLib.is_rep[localB];
                        int clusterB = directLib.seq_to_cluster[localB];

                        // Cross-Cluster
                        if (clusterA != clusterB) {
                            if (isRepA && isRepB) {
                                float S_A_rep = directLib.residueSupportRep[localA][posA];
                                // denom = std::min(S_A_rep, S_B_rep);
                                denom = (S_A_rep + S_B_rep) / 2.0f;
                            }
                            else if (isRepA && !isRepB) {
                                denom = static_cast<float>(directLib.cluster_reps[clusterB].size());
                            } 
                            else if (!isRepA && isRepB) {
                                denom = static_cast<float>(directLib.cluster_reps[clusterA].size());
                            }
                        }
                        else {
                            if ((isRepA && isRepB) || (!isRepA && !isRepB)) {
                                float S_A = directLib.residueSupport[localA][posA];
                                denom = (S_A + S_B) / 2.0f;
                            }
                        }
                    }
                    else {
                        denom = (S_A + S_B) / 2.0f;
                    }
                    // =========================================================
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

    // 輔助函式：快速計算 Pair ID
    auto get_pair_id = [this](int a, int b) {
        int _a = std::min(a, b);
        int _b = std::max(a, b);
        return _a * totalSequence - _a * (_a + 1) / 2 + (_b - _a - 1);
    };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, total_targets), [&](const tbb::blocked_range<size_t>& range) {
        
        // 🔥 優化 2：Thread-Local 緩衝區，避免頻繁 Heap Allocation
        std::vector<float> num_V;
        std::vector<int> active_V;
        
        // 定義快取結構，使用 Raw Pointer 達到極致速度
        struct BridgeCache {
            const int* uc_map;
            int uc_size;
            const int* cv_map;
            int cv_size;
            float weight;
        };
        std::vector<BridgeCache> valid_bridges;

        for (size_t i = range.begin(); i < range.end(); ++i) {
            int u = target_locals[i];
            for (size_t j = i + 1; j < total_targets; ++j) {
                int v = target_locals[j];
                
                int pId_uv = get_pair_id(u, v);
                int rIdx_uv = pair_to_record_idx[pId_uv];
                if (rIdx_uv == -1) continue;

                auto& rec_uv = records[rIdx_uv];
                int lenU = rec_uv.forward.size();
                int lenV = rec_uv.backward.size();

                rec_uv.extendedForward.assign(lenU, std::vector<ExtendedWeight>());
                rec_uv.extendedBackward.assign(lenV, std::vector<ExtendedWeight>());

                // 複用 Thread-Local 緩衝區
                // 只需要確保陣列夠大就好。
                if (num_V.size() < lenV) {
                    num_V.resize(lenV, 0.0f); 
                }
                active_V.clear();
                valid_bridges.clear();

                // 🔥 優化 1：將 Bridge 的查表與權重計算「提取」到 posU 迴圈外部
                for (int c : target_locals) {
                    if (c == u || c == v) continue;

                    int rIdx_uc = pair_to_record_idx[get_pair_id(u, c)];
                    if (rIdx_uc == -1) continue;
                    
                    int rIdx_cv = pair_to_record_idx[get_pair_id(c, v)];
                    if (rIdx_cv == -1) continue;

                    float w_uc = records[rIdx_uc].weight;
                    float w_cv = records[rIdx_cv].weight;
                    
                    if (w_uc > 0.0f && w_cv > 0.0f) {
                        BridgeCache cache;
                        cache.weight = std::min(w_uc, w_cv);

                        // 判斷方向並獲取 Raw Pointer 與 Size
                        const auto& rec_uc = records[rIdx_uc];
                        if (u < c) {
                            cache.uc_map = rec_uc.forward.data();
                            cache.uc_size = rec_uc.forward.size();
                        } else {
                            cache.uc_map = rec_uc.backward.data();
                            cache.uc_size = rec_uc.backward.size();
                        }

                        const auto& rec_cv = records[rIdx_cv];
                        if (c < v) {
                            cache.cv_map = rec_cv.forward.data();
                            cache.cv_size = rec_cv.forward.size();
                        } else {
                            cache.cv_map = rec_cv.backward.data();
                            cache.cv_size = rec_cv.backward.size();
                        }

                        valid_bridges.push_back(cache); // 只把有用的橋樑存起來
                    }
                }

                float w_uv = rec_uv.weight;

                // 🔥 優化 3：極速的最內層迴圈
                for (int posU = 0; posU < lenU; ++posU) {
                    int directV = rec_uv.forward[posU];
                    if (directV != -1) {
                        num_V[directV] = w_uv;
                        active_V.push_back(directV);
                    }

                    // 現在只需遍歷「確定有連接」的橋樑，且只剩下純粹的陣列索引操作
                    for (const auto& bc : valid_bridges) {
                        if (posU >= bc.uc_size) continue;
                        
                        int posC = bc.uc_map[posU];
                        if (posC == -1 || posC >= bc.cv_size) continue;
                        
                        int posV = bc.cv_map[posC];
                        if (posV == -1 || posV < 0 || posV >= lenV) continue;

                        if (num_V[posV] == 0.0f && directV != posV) {
                            active_V.push_back(posV);
                        }
                        num_V[posV] += bc.weight;
                    }

                    for (int posV : active_V) {
                        float w = num_V[posV];
                        if (w > 0.0f) {
                            rec_uv.extendedForward[posU].push_back({posV, w});
                            rec_uv.extendedBackward[posV].push_back({static_cast<int32_t>(posU), w});
                        }
                        num_V[posV] = 0.0f; // 寫完立刻歸零，省去 memset 開銷
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
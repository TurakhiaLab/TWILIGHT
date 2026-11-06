#include "block_manager.hpp"
#include "mga.hpp"
#include "type.hpp"
#include "message_manager.hpp"
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cctype>
#include <chrono>
#include <unordered_set>

// Helper to truncate strings for printing
std::string truncate(const std::string& str, size_t width) {
    if (str.length() > width) {
        return str.substr(0, width - 3) + "...";
    }
    return str;
}





void BlockManager::orientCircularGenomes(Option& option, std::string refSequenceName) {
    auto func_start_time = std::chrono::high_resolution_clock::now();
    bool filter = false;
    bool write_rotate = option.writeOriented;
    auto start_time = func_start_time;
    auto log_step_time = [&](const char* step_name) {
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = current_time - start_time;
        mga::msg::debug("%-30s | Elapsed: %.2f ms\n", step_name, elapsed.count());
    };

    mga::msg::debug("\n--- Start orientCircularGenomes ---\n");

    if (blocksets.empty() || sequence_lengths.empty()) {
        mga::msg::warn("No blocksets or sequences to orient.\n");
        return;
    }

    // 1. Determine Reference Sequence Name
    std::string bestRefName = refSequenceName;
    if (!bestRefName.empty() && blocksets.find(bestRefName) == blocksets.end()) {
        mga::msg::warn("Specified reference '%s' not found in blocksets. Falling back to the longest genome.\n",
                       bestRefName.c_str());
        bestRefName.clear();
    }

    if (bestRefName.empty()) { 
        int maxLen = -1;
        for (const auto& pair : sequence_lengths) {
            if (blocksets.find(pair.first) != blocksets.end() && pair.second > maxLen) {
                maxLen = pair.second;
                bestRefName = pair.first;
            }
        }
    }

    if (bestRefName.empty() || blocksets.find(bestRefName) == blocksets.end()) {
        mga::msg::error("No valid reference genome found in blocksets for circular orientation.\n");
        return;
    }

    mga::msg::debug("Selected Reference: %s (Length: %d)\n", bestRefName.c_str(), sequence_lengths[bestRefName]);
    log_step_time("Reference Selected");

    // 2. Extract Reference and Query sequences directly from this->sequences (zero copy using SequenceRefs, avoids costly rebuildAllPointers)
    SequenceRefs reference_genome;
    reference_genome.push_back({bestRefName, this->sequences[bestRefName]});
    
    SequenceRefs query_genomes;
    query_genomes.reserve(this->sequences.size() - 1);
    for (const auto& pair : this->sequences) {
        if (pair.first != bestRefName) {
            query_genomes.push_back({pair.first, pair.second});
        }
    }
    log_step_time("Prepared Queries for minimap2");

    // 3. Run Minimap2 (orientation only requires chaining length, CIGAR is disabled and max_gap is reduced for speed)
    mga::msg::debug("Running chaining-only minimap2 on %zu queries...\n", query_genomes.size());
    
    mm_config_orientation fastConfig(option);
    
    auto allAlignments = runMinimap2(reference_genome, query_genomes, "reference_genome", "query_genomes", option, fastConfig);
    log_step_time("Minimap2 Completed");

    // 4. Strand Voting & Finding Hybrid Anchor
    struct StrandStats {
        long long forwardLength = 0;
        long long reverseLength = 0;
    };
    std::map<std::string, StrandStats> qryStatsMap;
    for (const auto& aln : allAlignments) {
        if (!aln.valid) continue; 
        if (aln.inverse) {
            qryStatsMap[aln.qryName].reverseLength += aln.alnLength;
        } else {
            qryStatsMap[aln.qryName].forwardLength += aln.alnLength;
        }
    }
    
    std::map<std::string, bool> rcDecisionMap;
    for (const auto& pair : qryStatsMap) {
        rcDecisionMap[pair.first] = (pair.second.reverseLength > pair.second.forwardLength);
    }

    // Find an alignment that is sufficiently long and closest to the reference start as the shift baseline
    std::map<std::string, const Alignment*> bestAlnMap; 
    std::map<std::string, int> minRefStartMap; // Track closest distance to position 0
    
    for (const auto& aln : allAlignments) {
        if (!aln.valid) continue;
        bool needsRC = rcDecisionMap[aln.qryName];
        if (aln.inverse != needsRC) continue; // Ignore noise alignments in the wrong orientation

        // Threshold: minimum alignment length requirement
        int lengthThreshold = 10000;
        if (aln.alnLength < lengthThreshold) continue;

        // Among candidates meeting length criteria, pick the one with ref.start closest to 0
        auto it = minRefStartMap.find(aln.qryName);
        if (it == minRefStartMap.end() || aln.refIdx.first < it->second) {
            minRefStartMap[aln.qryName] = aln.refIdx.first;
            bestAlnMap[aln.qryName] = &aln;
        }
    }
    log_step_time("Strand Voting & Hybrid Anchor Found");
    
    // Backup original sequences before reverse complement and rotation
    std::unordered_map<std::string, std::string> original_sequences = this->sequences;

    // 5. Apply Reverse Complement and Circular Shift (Rotation)
    int rcCount = 0;
    int shiftCount = 0;

    for (const auto& pair : this->blocksets) {
        std::string qryName = pair.first;
        if (qryName == bestRefName) continue; // Skip reference itself

        bool needsRC = rcDecisionMap[qryName];
        std::string seq = this->sequences[qryName];
        int L = sequence_lengths[qryName];
        int raw_shift = 0;

        // Compute shift offset if an anchor alignment was found
        if (bestAlnMap.count(qryName)) {
            const Alignment* bestAln = bestAlnMap[qryName];
            if (needsRC) {
                // Coordinate shift conversion under reverse-complement state
                raw_shift = (L - bestAln->qryIdx.second) - bestAln->refIdx.first;
            } else {
                raw_shift = bestAln->qryIdx.first - bestAln->refIdx.first;
            }
        }

        // Normalize shift offset to a non-negative integer (circular modulo)
        int true_shift = ((raw_shift % L) + L) % L;

        // --- Apply Transformations ---
        if (needsRC) {
            seq = getReverseComplement(seq);
            rcCount++;
        }
        
        if (true_shift != 0) {
            // Circular sequence shift (rotate the first true_shift bases to the end)
            seq = seq.substr(true_shift) + seq.substr(0, true_shift);
            shiftCount++;
        }

        mga::msg::debug("%-15s | RC: %s (Forward: %lld bp) vs (Reverse: %lld bp) | Shifted: %d bp\n",
                        qryName.c_str(), (needsRC ? "YES" : "NO "),
                        qryStatsMap[qryName].forwardLength, qryStatsMap[qryName].reverseLength,
                        true_shift);

        // Write back to original data structures
        BlockSet* qryBlockSet = pair.second.get();
        auto first_block = qryBlockSet->getAllBlocks()[0];
        if (auto sp = first_block.lock()) {
            sp->setConsensus(Consensus(seq)); 
        }
        qryBlockSet->rebuildAllPointers();
        this->sequences[qryName] = seq;
    }

    mga::msg::debug("Inverted: %d | Shifted: %d genomes.\n", rcCount, shiftCount);
    log_step_time("Transformations Applied");

    // ==========================================================
    // 6. Sequence Filtering (retains sequences with clear strand and >90% coverage when filter == true)
    // ==========================================================

    long long refLen = sequence_lengths[bestRefName];

    // Collect qualified sequence names (Reference itself is the first)
    std::vector<std::string> qualified_names;
    qualified_names.push_back(bestRefName);

    for (const auto& pair : this->sequences) {
        std::string qryName = pair.first;
        if (qryName == bestRefName) continue;

        if (filter) {
            long long f_len = qryStatsMap[qryName].forwardLength;
            long long r_len = qryStatsMap[qryName].reverseLength;

            bool pass = (refLen > 0) && (((double)f_len / refLen > 0.9) || ((double)r_len / refLen > 0.9));
            if (pass) {
                qualified_names.push_back(qryName);
            } else {
                mga::msg::debug("[FILTERED OUT] %s (Forward: %lld bp, Reverse: %lld bp, RefLen: %lld bp)\n",
                                qryName.c_str(), f_len, r_len, refLen);
            }
        } else {
            qualified_names.push_back(qryName);
        }
    }

    if (filter || mga::msg::isVerbose()) {
        mga::msg::debug("\n================ Sequence Filtering Summary ================\n"
                        "  Filter Status      : %s\n"
                        "  Total input genomes: %zu\n"
                        "  Filtered out       : %zu\n"
                        "  Remaining genomes  : %zu\n"
                        "============================================================\n\n",
                        filter ? "ENABLED" : "DISABLED",
                        this->sequences.size(),
                        this->sequences.size() - qualified_names.size(),
                        qualified_names.size());
    }

    // Write all qualified sequences to FASTA files (output all remaining sequences without target count restriction)
    if (filter) {
        std::string outDir = option.tempDir.empty() ? "." : option.tempDir;
        size_t count = qualified_names.size();
        StringPairs orig_seqs;
        StringPairs orient_seqs;
        orig_seqs.reserve(count);
        orient_seqs.reserve(count);

        for (size_t i = 0; i < count; ++i) {
            const std::string& name = qualified_names[i];
            orig_seqs.push_back({name, original_sequences[name]});
            orient_seqs.push_back({name, this->sequences[name]});
        }

        std::string origFileName = outDir + "/filtered.original.fa";
        std::string orientFileName = outDir + "/filtered.orientated.fa";

        bool isCompressed = false;
        bool appendMode = false;

        mga::io::writeAlignment(origFileName, orig_seqs, isCompressed, appendMode);
        mga::io::writeAlignment(orientFileName, orient_seqs, isCompressed, appendMode);

        mga::msg::debug("Wrote all %zu remaining sequences to '%s' and '%s'\n",
                        count, origFileName.c_str(), orientFileName.c_str());
        log_step_time("Output FASTA Files Written");

        // Filter out unqualified sequences from internal structures
        std::unordered_set<std::string> qualified_set(qualified_names.begin(), qualified_names.end());
        for (auto it = this->sequences.begin(); it != this->sequences.end(); ) {
            if (qualified_set.find(it->first) == qualified_set.end()) {
                it = this->sequences.erase(it);
            } else {
                ++it;
            }
        }
        for (auto it = this->blocksets.begin(); it != this->blocksets.end(); ) {
            if (qualified_set.find(it->first) == qualified_set.end()) {
                it = this->blocksets.erase(it);
            } else {
                ++it;
            }
        }
        for (auto it = sequence_lengths.begin(); it != sequence_lengths.end(); ) {
            if (qualified_set.find(it->first) == qualified_set.end()) {
                it = sequence_lengths.erase(it);
            } else {
                ++it;
            }
        }
    }

    mga::msg::debug("--- End orientCircularGenomes ---\n\n");
    
    /*
    // ==========================================================
    // 6. Output converted full sequences (Reference + Oriented Queries)
    // ==========================================================
    const size_t SAMPLE_SIZE = 50; // Maximum output sample size (including Reference)
    StringPairs output_seqs;
    
    // 6-1. Place Reference as first entry to maintain alignment coordinate baseline
    output_seqs.push_back({bestRefName, this->sequences[bestRefName]});
    
    // 6-2. Collect all other Query sequences (apply quality pass and sample size limit)
    int passed_quality_count = 0;
    for (const auto& pair : this->sequences) {
        if (pair.first != bestRefName) {
            std::string qryName = pair.first;
            
            // Get forward and reverse alignment lengths for this query
            long long f_len = qryStatsMap[qryName].forwardLength;
            long long r_len = qryStatsMap[qryName].reverseLength;
            
            // Quality Pass: length must overwhelmingly favor one strand (> 3x)
            bool pass_quality = (f_len > 3 * r_len) || (r_len > 3 * f_len);

            if (pass_quality) {
                passed_quality_count++;
                // Check if sample count limit is reached
                if (output_seqs.size() < SAMPLE_SIZE) {
                    output_seqs.push_back({qryName, pair.second});
                }
            } else {
                // Failed quality pass (e.g. 50-50 inversion or no alignment)
                if (write_rotate) {
                    mga::msg::debug("[FILTERED] %-15s excluded from output (F: %lld vs R: %lld). Strand ambiguous.\n",
                                    qryName.c_str(), f_len, r_len);
                }
            }
        }
    }

    // 6-3. Invoke IO to write file
    std::string outFileName = option.tempDir + "/oriented_genomes.fasta"; 
    bool isCompressed = false; 
    bool appendMode = false;   
    
    if (write_rotate) mga::io::writeAlignment(outFileName, output_seqs, isCompressed, appendMode);
    
    if (write_rotate) {
        mga::msg::debug("Total queries passed quality: %d\n", passed_quality_count);
        mga::msg::debug("Wrote %zu oriented sequences (capped at %zu) to '%s'\n",
                        output_seqs.size(), SAMPLE_SIZE, outFileName.c_str());
        log_step_time("Output Sequences Written");
    }
    */
    // ==========================================================

    auto func_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> func_elapsed = func_end_time - func_start_time;
    mga::msg::info("orientCircularGenomes execution time: %.2f ms (%.3f s)\n",
                   func_elapsed.count(), func_elapsed.count() / 1000.0);
}

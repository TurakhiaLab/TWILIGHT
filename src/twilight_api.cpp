#include "twilight.h"

#include "msa.hpp"
#include "score.hpp"
#include "version.hpp"

#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

// RAII guard to restore MessageManager state when leaving C API calls
struct LogStateGuard {
    LogLevel prevLevel;
    bool prevVerbose;
    bool prevQuiet;
    bool prevExitOnError;

    LogStateGuard(bool quiet, bool verbose) {
        auto& mgr = MessageManager::getInstance();
        prevLevel = mgr.getLogLevel();
        prevVerbose = mgr.isVerbose();
        prevQuiet = mgr.isQuiet();
        prevExitOnError = mgr.isExitOnError();

        mgr.setExitOnError(false);
        if (quiet) {
            mgr.setLogLevel(LogLevel::NONE);
            mgr.setQuiet(true);
            mgr.setVerbose(false);
        } else {
            mgr.setLogLevel(LogLevel::INFO);
            mgr.setQuiet(false);
            mgr.setVerbose(verbose);
        }
    }

    ~LogStateGuard() {
        auto& mgr = MessageManager::getInstance();
        mgr.waitEmpty();
        mgr.setLogLevel(prevLevel);
        mgr.setVerbose(prevVerbose);
        mgr.setQuiet(prevQuiet);
        mgr.setExitOnError(prevExitOnError);
    }
};

tw_status_t mapErrorCode(msa::ErrorCode code) {
    switch (code) {
        case msa::ErrorCode::OK:          return TW_OK;
        case msa::ErrorCode::INVALID_ARG: return TW_ERR_INVALID_ARG;
        case msa::ErrorCode::EMPTY_INPUT: return TW_ERR_EMPTY_INPUT;
        case msa::ErrorCode::TREE_PARSE:  return TW_ERR_TREE_PARSE;
        case msa::ErrorCode::IO_ERROR:    return TW_ERR_IO;
        case msa::ErrorCode::OOM:         return TW_ERR_OOM;
        case msa::ErrorCode::UNSUPPORTED: return TW_ERR_UNSUPPORTED;
        case msa::ErrorCode::INTERNAL:
        default:                          return TW_ERR_INTERNAL;
    }
}

char detectSeqTypeFromMemory(const tw_seq_t* seqs, size_t n_seqs) {
    size_t scanned = 0;
    for (size_t i = 0; i < n_seqs && scanned < 100; ++i) {
        if (!seqs[i].seq) continue;
        size_t len = (seqs[i].seq_len > 0) ? seqs[i].seq_len : std::strlen(seqs[i].seq);
        if (len == 0) continue;
        for (size_t j = 0; j < len; ++j) {
            char c = static_cast<char>(std::toupper(static_cast<unsigned char>(seqs[i].seq[j])));
            char detected = checkOnly(c);
            if (detected != 'x') {
                return detected;
            }
        }
        ++scanned;
    }
    return 'n';
}

msa::CliOptions buildCliOptions(const tw_opt_t& opt, char resolvedType) {
    msa::CliOptions cli;
    cli.inMemory = true;
    cli.verbose = (opt.verbose != 0);
    cli.quiet = (opt.quiet != 0);
    cli.cpuOnly = (opt.cpu_only != 0);

    if (opt.num_threads > 0) {
        cli.cpuNum = opt.num_threads;
        cli.hasCpu = true;
    }

    cli.typeStr = std::string(1, resolvedType);
    cli.maxSubtree = INT32_MAX;
    cli.hasMaxSubtree = false;

    cli.removeGappy = opt.remove_gappy;
    cli.wildcard = (opt.wildcard != 0);
    cli.rooted = (opt.reroot == 0);
    cli.prune = (opt.prune_tree != 0);

    if (opt.len_dev > 0.0f) {
        cli.lenDev = opt.len_dev;
        cli.hasLenDev = true;
    }
    cli.maxAmbig = opt.max_ambig;
    cli.hasMaxAmbig = (opt.max_ambig != 0.1f && opt.max_ambig != 1.0f);
    if (opt.max_len > 0) {
        cli.maxLen = opt.max_len;
        cli.hasMaxLen = true;
    }
    if (opt.min_len > 0) {
        cli.minLen = opt.min_len;
        cli.hasMinLen = true;
    }
    cli.filter = (opt.filter_low_qual != 0);

    cli.match = opt.match;
    cli.mismatch = opt.mismatch;
    cli.transition = opt.transition;
    if (opt.gap_open <= 0.0f) {
        cli.gapOpen = opt.gap_open;
        cli.hasGapOpen = true;
    }
    if (opt.gap_extend <= 0.0f) {
        cli.gapExtend = opt.gap_extend;
        cli.hasGapExtend = true;
    }
    if (opt.gap_ends <= 0.0f) {
        cli.gapEnds = opt.gap_ends;
        cli.hasGapEnds = true;
    }
    if (opt.scale > 0.0f) {
        cli.scale = opt.scale;
        cli.hasScale = true;
    }
    cli.xdrop = opt.xdrop;
    cli.blosum = opt.blosum;

    cli.accuracyLevelStr = std::to_string(opt.accuracy_level);
    cli.hasAccuracyLevel = true;

    cli.sampleMode = opt.sample_mode;
    cli.hasSampleMode = true;
    if (opt.sample_max_group > 0) {
        cli.sampleMaxGroup = opt.sample_max_group;
        cli.hasSampleMaxGroup = true;
    }
    if (opt.sample_rate > 0.0f) {
        cli.sampleRate = opt.sample_rate;
        cli.hasSampleRate = true;
    }
    if (opt.sample_top_centrality > 0.0f) {
        cli.sampleTopCentrality = opt.sample_top_centrality;
        cli.hasSampleTopCentrality = true;
    }
    cli.consistencyWeight = opt.consistency_weight;

    if (opt.refine_subtree == 0) {
        cli.noRefineSubtree = true;
    } else {
        cli.noRefineSubtree = false;
    }
    cli.refineMethod = (opt.refine_method == TW_REFINE_NJ) ? "nj" : "upgma";

    if (opt.iterative_refine == 0) {
        cli.noIterativeRefine = true;
        cli.iterativeRefine = false;
    } else if (opt.iterative_refine == 1) {
        cli.noIterativeRefine = false;
        cli.iterativeRefine = true;
    } else {
        cli.noIterativeRefine = false;
        cli.iterativeRefine = false;
    }
    cli.refineIters = opt.refine_iters;
    cli.refineBatchSize = opt.refine_batch_size;

    if (opt.compute_score > 0) {
        cli.score = opt.compute_score;
        cli.hasScore = true;
    }

    return cli;
}

void loadSequencesFromMemory(
    const tw_seq_t* seqs,
    size_t n_seqs,
    msa::SequenceDB* database,
    msa::Option* option,
    phylogeny::Tree*& tree,
    std::vector<int>& dbIdxToInputIdx)
{
    auto seqReadStart = std::chrono::high_resolution_clock::now();

    int seqNum = 0;
    int maxLen = 0;
    int minLen = INT32_MAX;
    uint64_t totalLen = 0;
    std::vector<int> seqsLens;
    seqsLens.reserve(n_seqs);
    dbIdxToInputIdx.clear();
    dbIdxToInputIdx.reserve(n_seqs);

    for (size_t i = 0; i < n_seqs; ++i) {
        if (!seqs[i].name || !seqs[i].seq) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Input sequence entry has NULL name or seq pointer at index " + std::to_string(i) + ".");
        }
        size_t nameLen = (seqs[i].name_len > 0) ? seqs[i].name_len : std::strlen(seqs[i].name);
        size_t seqLen = (seqs[i].seq_len > 0) ? seqs[i].seq_len : std::strlen(seqs[i].seq);
        if (nameLen == 0) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Input sequence entry has empty name at index " + std::to_string(i) + ".");
        }

        std::string seqNameFull(seqs[i].name, nameLen);
        size_t spacePos = seqNameFull.find(' ');
        std::string seqNameNoBlank = (spacePos == std::string::npos) ? seqNameFull : seqNameFull.substr(0, spacePos);

        std::string matchedName;
        if (tree->allNodes.find(seqNameFull) != tree->allNodes.end()) {
            matchedName = seqNameFull;
        } else if (tree->allNodes.find(seqNameNoBlank) != tree->allNodes.end()) {
            matchedName = seqNameNoBlank;
        }

        if (!matchedName.empty()) {
            if (database->name_map.find(matchedName) != database->name_map.end()) {
                msg::warn("Duplicate leaf names found in input sequences! Leaf name: %s. Only the first occurrence will be kept.\n",
                          matchedName.c_str());
            } else {
                int subtreeIdx = tree->allNodes[matchedName]->grpID;
                int sLen = static_cast<int>(seqLen);
                if (sLen > maxLen) maxLen = sLen;
                if (sLen < minLen) minLen = sLen;
                if (sLen == 0) msg::warn("Null sequence: %s\n", matchedName.c_str());
                totalLen += seqLen;
                database->addSequenceFromBuffer(
                    seqNum,
                    matchedName.data(),
                    matchedName.size(),
                    seqs[i].seq,
                    seqLen,
                    subtreeIdx,
                    tree->allNodes[matchedName]->weight,
                    option->debug,
                    option->alnMode);
                tree->allNodes[matchedName]->placed = false;
                dbIdxToInputIdx.push_back(static_cast<int>(i));
                ++seqNum;
                seqsLens.push_back(sLen);
            }
        }
    }

    // Prune tree if some tree leaves were not provided in input sequences
    if (tree->m_numLeaves != static_cast<size_t>(seqNum)) {
        msg::warn("Mismatch between the number of leaves and the number of sequences, (%lu != %d)\n",
                  tree->m_numLeaves, seqNum);
        std::unordered_set<std::string> seqNames;
        seqNames.reserve(database->name_map.size());
        for (const auto& kv : database->name_map) {
            seqNames.insert(kv.first);
        }
        phylogeny::pruneTree(tree, seqNames);
    }

    if (seqNum == 0) {
        throw msa::TwilightError(msa::ErrorCode::EMPTY_INPUT,
            "No input sequences matched any leaf in the guide tree.");
    }

    std::sort(seqsLens.begin(), seqsLens.end());
    uint32_t avgLen = static_cast<uint32_t>(totalLen / seqNum);
    uint32_t medLen = static_cast<uint32_t>(seqsLens[seqNum / 2]);

    bool min_max = (option->minLen > 0 || option->maxLen < INT32_MAX);
    int minLenTh = min_max ? option->minLen : ((option->lenDev > 0) ? static_cast<int>(medLen * (1.0f - option->lenDev)) : option->minLen);
    int maxLenTh = min_max ? option->maxLen : ((option->lenDev > 0) ? static_cast<int>(medLen * (1.0f + option->lenDev)) : option->maxLen);

    std::atomic<int> numLowQ(0);
    tbb::parallel_for(tbb::blocked_range<int>(0, seqNum), [&](tbb::blocked_range<int> range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            auto* seq = database->sequences[i];
            int ambig = (option->type == 'n') ? 4 : 20;
            seq->lowQuality = (seq->len > maxLenTh || seq->len < minLenTh);
            if (!seq->lowQuality) {
                int ambigCount = 0;
                for (int j = 0; j < seq->len; ++j) {
                    if (letterIdx(option->type, static_cast<char>(std::toupper(static_cast<unsigned char>(seq->alnStorage[0][j])))) == ambig) {
                        ambigCount++;
                    }
                }
                seq->lowQuality = (ambigCount > (seq->len * option->maxAmbig));
            }
            if (seq->lowQuality) {
                numLowQ.fetch_add(1);
                if (!option->noFilter) {
                    seq->len = 0;
                }
            }
        }
    });

    auto seqReadEnd = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds seqReadTime = seqReadEnd - seqReadStart;

    if (option->accurate) {
        option->resolveAccurateParams(seqNum);
    }

    msg::info("===== Sequence Summary =====\n");
    msg::info("Number : %d\n", seqNum);
    msg::info("Max. Length: %d\n", maxLen);
    msg::info("Min. Length: %d\n", minLen);
    msg::info("Avg. Length: %u\n", avgLen);
    msg::info("Med. Length: %u\n", medLen);
    if (option->noFilter)
        msg::info("Deferred sequences: %d\n", numLowQ.load());
    else
        msg::info("Excluded sequences: %d\n", numLowQ.load());
    msg::info("Sequences loaded in %lld ms\n", static_cast<long long>(seqReadTime.count() / 1000000));
}

tw_msa_t* packMsaOutput(
    const msa::SequenceDB* database,
    const std::vector<int>& dbIdxToInputIdx,
    size_t alnLen,
    double score)
{
    // Count valid (non-excluded) sequences and total name bytes
    size_t outCount = 0;
    size_t totalNameBytes = 0;
    for (size_t i = 0; i < database->sequences.size(); ++i) {
        const auto* seq = database->sequences[i];
        if (seq && !seq->lowQuality) {
            ++outCount;
            totalNameBytes += seq->name.size() + 1;
        }
    }

    tw_msa_t* msa = static_cast<tw_msa_t*>(std::calloc(1, sizeof(tw_msa_t)));
    if (!msa) return nullptr;

    msa->n_seqs = outCount;
    msa->aln_len = alnLen;
    msa->score = score;

    if (outCount == 0) {
        return msa;
    }

    size_t seqStride = alnLen + 1;
    msa->names = static_cast<char**>(std::malloc(outCount * sizeof(char*)));
    msa->seqs = static_cast<char**>(std::malloc(outCount * sizeof(char*)));
    msa->input_indices = static_cast<int*>(std::malloc(outCount * sizeof(int)));
    msa->_seq_buf = static_cast<char*>(std::malloc(outCount * seqStride));
    msa->_name_buf = static_cast<char*>(std::malloc(totalNameBytes));

    if (!msa->names || !msa->seqs || !msa->input_indices || !msa->_seq_buf || !msa->_name_buf) {
        tw_msa_free(msa);
        return nullptr;
    }

    size_t row = 0;
    size_t nameOffset = 0;
    for (size_t i = 0; i < database->sequences.size(); ++i) {
        const auto* seq = database->sequences[i];
        if (!seq || seq->lowQuality) continue;

        char* nameDst = msa->_name_buf + nameOffset;
        std::memcpy(nameDst, seq->name.data(), seq->name.size());
        nameDst[seq->name.size()] = '\0';
        msa->names[row] = nameDst;
        nameOffset += seq->name.size() + 1;

        char* seqDst = msa->_seq_buf + row * seqStride;
        if (alnLen > 0) {
            std::memcpy(seqDst, seq->alnStorage[seq->storage], alnLen);
        }
        seqDst[alnLen] = '\0';
        msa->seqs[row] = seqDst;

        msa->input_indices[row] = dbIdxToInputIdx[i];
        ++row;
    }

    return msa;
}

} // anonymous namespace

extern "C" {

void tw_opt_init(tw_opt_t *opt) {
    if (!opt) return;
    std::memset(opt, 0, sizeof(tw_opt_t));

    opt->num_threads = 0;
    opt->cpu_only = 0;
    opt->quiet = 1;
    opt->verbose = 0;

    opt->seq_type = TW_SEQ_AUTO;
    opt->reroot = 1;
    opt->prune_tree = 1;

    opt->remove_gappy = 0.95f;
    opt->max_ambig = 0.1f;
    opt->len_dev = 0.0f;
    opt->min_len = 0;
    opt->max_len = 0;
    opt->filter_low_qual = 0;
    opt->wildcard = 0;

    opt->match = 18.0f;
    opt->mismatch = -8.0f;
    opt->transition = -4.0f;
    opt->gap_open = 1.0f;    /* > 0 means use type-specific default */
    opt->gap_extend = 1.0f;  /* > 0 means use type-specific default */
    opt->gap_ends = 1.0f;    /* > 0 means same as gap_extend */
    opt->scale = 0.0f;       /* <= 0 means use type-specific default */
    opt->xdrop = 600.0f;
    opt->blosum = 62;

    opt->accuracy_level = TW_ACC_STANDARD;
    opt->sample_mode = 0;
    opt->sample_max_group = 0;
    opt->sample_rate = 0.0f;
    opt->sample_top_centrality = 0.0f;
    opt->consistency_weight = 5.0f;
    opt->refine_subtree = -1;
    opt->refine_method = TW_REFINE_UPGMA;
    opt->iterative_refine = -1;
    opt->refine_iters = 1;
    opt->refine_batch_size = 4;

    opt->compute_score = 0;
}

tw_status_t tw_opt_set_preset(tw_opt_t *opt, const char *preset) {
    if (!opt || !preset) return TW_ERR_INVALID_ARG;
    std::string p(preset);
    std::transform(p.begin(), p.end(), p.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (p == "fast" || p == "standard" || p == "default" || p == "1") {
        opt->accuracy_level = TW_ACC_STANDARD;
        return TW_OK;
    }
    if (p == "accurate" || p == "mafft" || p == "pairwise" || p == "2") {
        opt->accuracy_level = TW_ACC_PAIRWISE;
        return TW_OK;
    }
    if (p == "tcs" || p == "transitive" || p == "3") {
        opt->accuracy_level = TW_ACC_TRANSITIVE;
        return TW_OK;
    }
    return TW_ERR_INVALID_ARG;
}

tw_status_t tw_align(
    const tw_opt_t *opt,
    const char *newick_tree,
    const tw_seq_t *seqs,
    size_t n_seqs,
    tw_msa_t **out_msa)
{
    if (!out_msa) return TW_ERR_INVALID_ARG;
    *out_msa = nullptr;

    if (!newick_tree || newick_tree[0] == '\0' || !seqs) {
        return TW_ERR_INVALID_ARG;
    }
    if (n_seqs < 2) {
        return TW_ERR_EMPTY_INPUT;
    }

    tw_opt_t localOpt;
    if (opt) {
        localOpt = *opt;
    } else {
        tw_opt_init(&localOpt);
    }

    if (localOpt.seq_type != TW_SEQ_AUTO && localOpt.seq_type != TW_SEQ_NUCL && localOpt.seq_type != TW_SEQ_PROTEIN) {
        return TW_ERR_INVALID_ARG;
    }
    if (localOpt.refine_method != TW_REFINE_UPGMA && localOpt.refine_method != TW_REFINE_NJ) {
        return TW_ERR_INVALID_ARG;
    }

    LogStateGuard logGuard(localOpt.quiet != 0, localOpt.verbose != 0);

    try {
        char resolvedType = 'n';
        if (localOpt.seq_type == TW_SEQ_NUCL) {
            resolvedType = 'n';
        } else if (localOpt.seq_type == TW_SEQ_PROTEIN) {
            resolvedType = 'p';
        } else {
            resolvedType = detectSeqTypeFromMemory(seqs, n_seqs);
        }

        msa::CliOptions cliOpt = buildCliOptions(localOpt, resolvedType);
        std::unique_ptr<msa::Option> option(new msa::Option(cliOpt));
        std::unique_ptr<msa::SequenceDB> database(new msa::SequenceDB());

#if defined(CUDA) || defined(__HIP_PLATFORM_AMD__)
        option->getGpuInfo(cliOpt);
#endif
        std::unique_ptr<msa::Params> param(new msa::Params(cliOpt, option->type));

        size_t finalAlnLen = 0;
        double finalScore = 0.0;
        std::vector<int> dbIdxToInputIdx;

        tbb::task_arena arena(option->cpuNum);
        arena.execute([&]() {
            std::unique_ptr<phylogeny::Tree> T(phylogeny::Tree::fromNewickString(std::string(newick_tree)));
            if (cliOpt.prune) {
                std::unordered_set<std::string> seqNames;
                seqNames.reserve(n_seqs);
                for (size_t i = 0; i < n_seqs; ++i) {
                    if (!seqs[i].name) {
                        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                            "Input sequence entry has NULL name at index " + std::to_string(i) + ".");
                    }
                    size_t nameLen = (seqs[i].name_len > 0) ? seqs[i].name_len : std::strlen(seqs[i].name);
                    std::string fullName(seqs[i].name, nameLen);
                    seqNames.insert(fullName);
                    size_t spacePos = fullName.find(' ');
                    if (spacePos != std::string::npos) {
                        seqNames.insert(fullName.substr(0, spacePos));
                    }
                }
                phylogeny::Tree* rawT = T.release();
                phylogeny::pruneTree(rawT, seqNames);
                T.reset(rawT);
            }

            T->convert2binaryTree();
            T->calLeafNum();

            std::unique_ptr<phylogeny::Tree> subT(new phylogeny::Tree(T->root, option->reroot));
            phylogeny::Tree* subTRaw = subT.release();
            try {
                loadSequencesFromMemory(seqs, n_seqs, database.get(), option.get(), subTRaw, dbIdxToInputIdx);
                if (!option->accurate) {
                    msa::progressive::msaOnSubtree(subTRaw, database.get(), option.get(), *param, nullptr, -1);
                } else {
                    msa::accurate::msaOnSubtree_accurate(subTRaw, database.get(), option.get(), *param, nullptr, -1);
                }
                finalAlnLen = static_cast<size_t>(subTRaw->root->getAlnLen(database->currentTask));
                if (option->scoreMode > 0) {
                    auto scoreRes = msa::score::computeAlignmentScore(
                        database.get(), option.get(), param.get(), static_cast<int>(finalAlnLen), option->scoreMode);
                    finalScore = scoreRes.score;
                }
            } catch (...) {
                delete subTRaw;
                throw;
            }
            delete subTRaw;
        });

        tw_msa_t* result = packMsaOutput(database.get(), dbIdxToInputIdx, finalAlnLen, finalScore);
        if (!result) {
            return TW_ERR_OOM;
        }
        *out_msa = result;
        return TW_OK;

    } catch (const msa::TwilightError& e) {
        return mapErrorCode(e.code());
    } catch (const std::bad_alloc&) {
        return TW_ERR_OOM;
    } catch (const std::exception&) {
        return TW_ERR_INTERNAL;
    } catch (...) {
        return TW_ERR_INTERNAL;
    }
}

void tw_msa_free(tw_msa_t *msa) {
    if (!msa) return;
    std::free(msa->names);
    std::free(msa->seqs);
    std::free(msa->input_indices);
    std::free(msa->_seq_buf);
    std::free(msa->_name_buf);
    std::free(msa);
}

const char *tw_strerror(tw_status_t status) {
    switch (status) {
        case TW_OK:              return "Success";
        case TW_ERR_INVALID_ARG: return "Invalid argument";
        case TW_ERR_EMPTY_INPUT: return "Empty input or no sequences matched the guide tree";
        case TW_ERR_TREE_PARSE:  return "Failed to parse Newick guide tree";
        case TW_ERR_IO:          return "I/O error";
        case TW_ERR_OOM:         return "Out of memory";
        case TW_ERR_UNSUPPORTED: return "Unsupported mode or option";
        case TW_ERR_INTERNAL:    return "Internal error";
        default:                 return "Unknown error";
    }
}

const char *tw_version(void) {
    return PROJECT_VERSION;
}

} /* extern "C" */

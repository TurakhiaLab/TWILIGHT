#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include <zlib.h>
#include <sys/stat.h>
#include <filesystem>
#include <tbb/parallel_for.h>

namespace fs = std::filesystem;

msa::Option::Option(const CliOptions& opt) {

    // Detect Alignment Mode
    if (opt.inMemory) {
        this->alnMode = DEFAULT_ALN;
    } else {
        bool hasTree = !opt.treeFile.empty(), hasSeq = !opt.seqFile.empty(), hasFile = !opt.msaDir.empty(), hasAln = !opt.backboneAlnFile.empty();
        int8_t mask = (hasFile << 3) | (hasTree << 2) | (hasSeq << 1) | (hasAln << 0);
        switch (mask) {
            case 0b0110: this->alnMode = 0; break;
            case 0b1000: this->alnMode = 1; break;
            case 0b0011: this->alnMode = 2; break;
            case 0b0111: this->alnMode = 3; break;
            default:
                throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                    "Unrecognized alignment mode based on the provided options.\n"
                    "Valid combinations are:\n"
                    "  [1] --tree and --sequences (for building MSA from unaligned sequences)\n"
                    "  [2] --files (for merging MSAs)\n"
                    "  [3] --sequences and --alignment (for adding new sequences to an existing MSA without guide tree)\n"
                    "  [4] --sequences, --alignment and --tree (for adding new sequences to an existing MSA using a given guide tree)\n"
                    "Please check the input arguments or run with --help for usage.");
        }
    }
    // Read Input File Names
    this->treeFile        = opt.treeFile;
    this->seqFile         = opt.seqFile;
    this->msaDir          = opt.msaDir;
    this->backboneAlnFile = opt.backboneAlnFile;

    // Read Options
    this->gpuNum = 0;
    this->gpuIdx = std::vector<int> (0);
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    this->cpuNum = (opt.hasCpu) ? opt.cpuNum : maxCpuThreads;
    if (this->cpuNum <= 0 || cpuNum > maxCpuThreads) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid number of CPU cores. Please request between 1 and " + std::to_string(maxCpuThreads) + " (got " + std::to_string(this->cpuNum) + ").");
    }

    this->maxSubtree = (opt.hasMaxSubtree) ? opt.maxSubtree : INT32_MAX;
    if (this->maxSubtree <= 0) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --max-subtree. The value of --max-subtree should be a positive integer (got " + std::to_string(this->maxSubtree) + ").");
    }
    this->gappyVertical = opt.removeGappy;
    if (this->gappyVertical > 1 || this->gappyVertical <= 0) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --remove-gappy. The value of --remove-gappy should be in (0,1]");
    }
    this->lenDev = opt.hasLenDev ? opt.lenDev : 0;
    if (this->lenDev < 0) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --length-deviation. The value of --length-deviation should be non-negative");
    }
    this->maxAmbig = opt.maxAmbig;
    if (this->maxAmbig > 1 || this->maxAmbig <= 0) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --max-ambig. The value of --max-ambig should be in (0,1]");
    }
    this->maxLen = opt.hasMaxLen ? opt.maxLen : INT32_MAX;
    if (this->maxLen <= 0 && opt.hasMaxLen) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --max-len. The value of --max-len should be larger than 0 (got " + std::to_string(this->maxLen) + ").");
    }
    this->minLen =  opt.hasMinLen ? opt.minLen : 0;
    if (this->minLen <= 0 && opt.hasMinLen) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --min-len. The value of --min-len should be larger than 0 (got " + std::to_string(this->minLen) + ").");
    }
    if (this->minLen >= this->maxLen) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --min-len and --max-len. --max-len should be larger than --min-len (got --max-len: " + std::to_string(this->maxLen) + ", --min-len: " + std::to_string(this->minLen) + ").");
    }

    if ((opt.hasMinLen || opt.hasMaxLen) && opt.hasLenDev) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid arguments. --length-deviation cannot be used together with --min-len or --max-len.");
    }
    
    this->reroot = !opt.rooted;
    this->debug = opt.check;
    this->cpuOnly = opt.cpuOnly;
    this->printDetail = opt.verbose;
    this->quiet = opt.quiet;
    this->deleteTemp = !opt.keepTemp;
    this->alignGappy = !opt.noAlignGappy;
    this->compressed = opt.compress;
    this->noFilter = !opt.filter;
    this->writeFiltered = opt.writeFiltered;

    // Parse accuracy-level
    this->accuracyLevel = 1;
    this->accurate = false;
    this->accurateMethod = "pairwise";

    if (opt.hasAccuracyLevel) {
        std::string lvlStr = opt.accuracyLevelStr;
        std::string lvlLower = lvlStr;
        std::transform(lvlLower.begin(), lvlLower.end(), lvlLower.begin(), ::tolower);
        if (lvlLower == "1" || lvlLower == "standard" || lvlLower == "fast" || lvlLower == "default") {
            this->accuracyLevel = 1;
            this->accurate = false;
            this->accurateMethod = "pairwise";
        } else if (lvlLower == "2" || lvlLower == "mafft" || lvlLower == "pairwise" || lvlLower == "consistent") {
            this->accuracyLevel = 2;
            this->accurate = true;
            this->accurateMethod = "pairwise";
        } else if (lvlLower == "3" || lvlLower == "tcs" || lvlLower == "transitive") {
            this->accuracyLevel = 3;
            this->accurate = true;
            this->accurateMethod = "tcs";
        } else {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Invalid value for --accuracy-level: \"" + lvlStr +
                "\". Usage: --accuracy-level <1|2|3> (1: traditional sum-of-pairs scoring, 2: + pairwise consistency, 3: + transitive consistency).");
        }
    }

    this->consistencyWeight = opt.consistencyWeight;

    // Parse sample-mode
    this->sampleMode = opt.hasSampleMode ? opt.sampleMode : 0;
    this->accMode = this->sampleMode;
    if (this->sampleMode < 0 || this->sampleMode > 3) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --sample-mode. Must be 0 (auto), 1 (all-to-all), 2 (high sample rate), or 3 (low sample rate) (got " + std::to_string(this->sampleMode) + ").");
    }

    this->hasUserAccMaxGroup = opt.hasSampleMaxGroup;
    this->hasUserAccSampleRate = opt.hasSampleRate;
    this->hasUserAccSampleTop = opt.hasSampleTopCentrality;

    if (this->hasUserAccMaxGroup) {
        this->userAccMaxGroup = opt.sampleMaxGroup;
        if (this->userAccMaxGroup <= 1) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Invalid value for --sample-max-group. The value of --sample-max-group should be > 1 (got " + std::to_string(this->userAccMaxGroup) + ").");
        }
    }
    if (this->hasUserAccSampleRate) {
        this->userAccSampleRate = opt.sampleRate;
        if (this->userAccSampleRate <= 0.0f || this->userAccSampleRate >= 1.0f) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "--sample-rate must be greater than 0 and less than 1 (got " + std::to_string(this->userAccSampleRate) + ").");
        }
    }
    if (this->hasUserAccSampleTop) {
        this->userAccSampleTop = opt.sampleTopCentrality;
        if (this->userAccSampleTop <= 0.0f || this->userAccSampleTop > 1.0f) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "--sample-top-centrality must be greater than 0 and less than or equal to 1 (got " + std::to_string(this->userAccSampleTop) + ").");
        }
    }

    // Parse refine-subtree options: default ON in accurate mode (-A 2, -A 3), OFF otherwise.
    if (opt.noRefineSubtree) {
        this->refineSubtree = false;
    } else if (this->accurate) {
        this->refineSubtree = true;
    } else {
        this->refineSubtree = false;
    }
    this->refineMethod = opt.refineMethod;
    std::transform(this->refineMethod.begin(), this->refineMethod.end(), this->refineMethod.begin(), ::tolower);
    if (this->refineMethod != "upgma" && this->refineMethod != "nj") {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Invalid value for --refine-method: \"" + this->refineMethod + "\". Valid options are: \"upgma\" or \"nj\".");
    }

    // Parse alignment iterative refinement options: default ON in accurate mode (-A 2, -A 3), OFF otherwise.
    if (opt.noIterativeRefine) {
        this->iterativeRefine = false;
    } else if (this->accurate) {
        this->iterativeRefine = true;
    } else if (opt.iterativeRefine) {
        this->iterativeRefine = true;
    } else {
        this->iterativeRefine = false;
    }
    this->refineIters = std::max(1, opt.refineIters);
    this->refineBatchSize = (opt.refineBatchSize > 0) ? opt.refineBatchSize : 4;

    // Initialize with 0 sequences first (will be dynamically resolved per dataset / subtree)
    this->resolveAccurateParams(0);

    if (this->accurate) {
        if (opt.hasMaxAmbig || opt.hasLenDev || opt.hasMinLen || opt.hasMaxLen || opt.filter || opt.writeFiltered) {
            throw msa::TwilightError(msa::ErrorCode::UNSUPPORTED,
                "Sequence filtering/deferring options (--max-ambig, --length-deviation, --min-len, --max-len, --filter, --write-filtered) are currently not supported in accurate mode (-A 2 / -A 3).");
        }
        // Disable ambiguity-based deferring in accurate mode by default
        this->maxAmbig = 1.0f;
        if (this->alnMode != DEFAULT_ALN && this->alnMode != ACCURATE_ALN) {
            throw msa::TwilightError(msa::ErrorCode::UNSUPPORTED,
                "Accuracy levels 2 and 3 (consistency mode) are currently only supported for guide-tree alignment (--tree and --sequences).");
        }
        if (this->alnMode == DEFAULT_ALN) {
            this->alnMode = ACCURATE_ALN;
        }
    }

    // Detect Data Type
    if (!opt.typeStr.empty()) {
        if (opt.typeStr != "n" && opt.typeStr != "p") {
            if (!tempDir.empty()) fs::remove_all(tempDir);
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Unrecognized data type \"" + opt.typeStr + "\".");
        }
        else this->type = opt.typeStr[0];
    }
    else {
        std::string seqFile = "";
        if (this->msaDir == "") seqFile = this->seqFile;
        else {
            for (const auto& file : fs::directory_iterator(this->msaDir)) {
                seqFile = file.path().string();
                break;
            }
        }
        if (seqFile.size() >= 3 && seqFile.substr(seqFile.size()-3, 3) == ".gz") {
            gzFile file = gzopen(seqFile.c_str(), "r");
            if (!file) {
                throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                    "Failed to open file " + seqFile + ".");
            }
            std::string line;
            int lineCount = 0;
            bool fin = false;
            this->type = 'n';
            char buffer[4096];
            while (gzgets(file, buffer, 4096)) {
                line = buffer;
                if (line.empty() || line[0] == '>') continue;
                for (char c : line) {
                    c = toupper(c);
                    char seqType = checkOnly(c);
                    if (seqType != 'x') {
                        this->type = seqType;
                        fin = true;
                        break;
                    }
                }
                ++lineCount;
                if (fin || lineCount == 100) break;
            }
            gzclose(file);
        }
        else {
            std::ifstream file(seqFile);
            if (!file.is_open()) {
                if (!tempDir.empty()) fs::remove_all(tempDir);
                throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                    "Cannot open " + seqFile);
            }
            std::string line;
            int lineCount = 0;
            bool fin = false;
            this->type = 'n';
            while (getline(file, line)) {
                if (line.empty() || line[0] == '>') continue;
                for (char c : line) {
                    c = toupper(c);
                    char seqType = checkOnly(c);
                    if (seqType != 'x') {
                        this->type = seqType;
                        fin = true;
                        break;
                    }
                }
                ++lineCount;
                if (fin || lineCount == 100) break;
            }
        }
    }

    // Check Output File (skip filesystem checks in in-memory API mode)
    if (opt.inMemory) {
        this->outFile = "";
        this->scoreMode = opt.hasScore ? opt.score : 0;
        if (this->scoreMode != 0 && this->scoreMode != 1 && this->scoreMode != 2) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Invalid value for score mode: " + std::to_string(this->scoreMode) +
                ". Valid options are 0 (disabled), 1 (Consistency objective score), or 2 (WSP).");
        }
        this->scoreFile = "";
    } else {
    if (opt.outFile.empty()) {
        throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
            "Output file name is required.");
    }
    this->outFile = opt.outFile;

    // Parse score options
    if (opt.hasScore) {
        this->scoreMode = opt.score;
        if (this->scoreMode != 1 && this->scoreMode != 2) {
            throw msa::TwilightError(msa::ErrorCode::INVALID_ARG,
                "Invalid value for --score: " + std::to_string(this->scoreMode) +
                ". Valid options are 1 (Consistency objective score) or 2 (WSP).");
        }
    } else if (!opt.scoreFile.empty()) {
        this->scoreMode = 1; // Default to consistency score if score-file is given but not score
    } else {
        this->scoreMode = 0; // Disabled by default
    }

    if (!opt.scoreFile.empty()) {
        this->scoreFile = opt.scoreFile;
    } else if (this->scoreMode > 0) {
        this->scoreFile = this->outFile + ".score";
    } else {
        this->scoreFile = "";
    }
    if (!opt.overwrite) {
        std::string outFileStr = (this->compressed) ? this->outFile + ".gz" : this->outFile;
        if (fs::exists(outFileStr)) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                outFileStr + " already exists. Please use another file name or add --overwrite to overwrite the existing file.");
        }
    }
    if (this->compressed) {
        gzFile outFile = gzopen(this->outFile.c_str(), "wb"); // "wb" = write binary
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                "Failed to open file: " + this->outFile);
        }
        gzclose(outFile);
        fs::remove(this->outFile);
    }
    else {
        std::ofstream outFile(this->outFile);
        if (!outFile) {
            throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                "Failed to open file: " + this->outFile);
        }
        outFile.close();
        fs::remove(this->outFile);
    }
    

    // Create Temporary Directory if Needed
    if (this->maxSubtree < INT32_MAX || !opt.msaDir.empty() || this->alnMode == 2) {
        std::string tempDir;
        if (opt.tempDir.empty()) {
            int idx = 1;
            std::string tempDir_org = "./twilight_temp";
            tempDir = tempDir_org;
            while (true) {
                if (mkdir(tempDir.c_str(), 0777) == -1) {
                    if( errno == EEXIST ) {
                        tempDir = tempDir_org + '_' + std::to_string(idx);
                        ++idx;
                    }
                    else {
                        throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                            "Can't create directory: " + tempDir);
                    }
                }
                else break;
            }
        }
        else {
            tempDir = opt.tempDir;
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            if (fs::exists(tempDir)) {
                if (!opt.overwrite) {
                    throw msa::TwilightError(msa::ErrorCode::IO_ERROR,
                        tempDir + " already exists. In order to prevent your file from being overwritten, please delete this folder or use another folder name.");
                }
            }
            fs::create_directories(tempDir);
        }
        msg::info("%s created for storing temporary alignments\n", tempDir.c_str());
        this->tempDir = tempDir;
    }
    }

    
    msg::info("====== Configuration =======\n");
    if (this->accurate) {
        std::string methodDisplay = (this->accuracyLevel == 3) ? "3 (+ transitive consistency)" : "2 (+ pairwise consistency)";
        std::string modeStr = (this->sampleMode == 0) ? "0 (auto)" : (this->sampleMode == 1) ? "1 (all-to-all)" : (this->sampleMode == 2) ? "2 (high sample rate)" : "3 (low sample rate)";
        {
            auto s = msg::stream(LogLevel::INFO);
            s << "Accuracy level: " << methodDisplay << " (sample-mode=" << modeStr;
            if (this->hasUserAccMaxGroup) s << ", sample-max-group=" << this->userAccMaxGroup << " [user]";
            if (this->hasUserAccSampleRate) s << ", sample-rate=" << this->userAccSampleRate << " [user]";
            if (this->hasUserAccSampleTop) s << ", sample-top-centrality=" << this->userAccSampleTop << " [user]";
            s << ", consistency weight=" << this->consistencyWeight << ")\n";
        }
        msg::info("Layer 1 subtree refinement: %s\n", (this->refineSubtree ? ("enabled (" + this->refineMethod + ")") : std::string("disabled")).c_str());
        msg::info("Layer 1 iterative refinement: %s\n", (this->iterativeRefine ? ("enabled (max " + std::to_string(this->refineIters) + " iters, batch-size=" + std::to_string(this->refineBatchSize) + ")") : std::string("disabled")).c_str());
    } else {
        msg::info("Accuracy level: 1 (traditional sum-of-pairs scoring)\n");
    }
    if (this->maxSubtree != INT32_MAX)
        msg::info("Max-subtree: %d\n", this->maxSubtree);
    if (this->gappyVertical == 1)
        msg::info("Disable removing gappy columns.\n");
    else
        msg::stream(LogLevel::INFO) << "Threshold for removing gappy columns: " << this->gappyVertical << '\n';
    if (this->lenDev > 0)
        msg::stream(LogLevel::INFO) << "Allowed deviation from the median length: " << (this->lenDev * 100) << "%\n";
    else if (this->minLen > 0 || this->maxLen < INT32_MAX)
        msg::info("Allowed sequence length range: [%d, %d]\n", this->minLen, this->maxLen);
    if (this->maxAmbig < 1)
        msg::stream(LogLevel::INFO) << "Allowed proportion of ambiguous characters: " << (this->maxAmbig * 100) << "%\n";
    msg::info("Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);
}

void msa::Option::resolveAccurateParams(int numSeqs) {
    int effectiveMode = this->sampleMode;
    if (effectiveMode == 0) {
        if (numSeqs <= 500) {
            effectiveMode = 1;
        } else if (numSeqs <= 100000) {
            effectiveMode = 2;
        } else {
            effectiveMode = 3;
        }
    }
    this->resolvedAccMode = effectiveMode;

    int baseMaxGroup = 200;
    float baseSampleRate = 0.5f;
    float baseSampleTop = 1.0f;

    if (effectiveMode == 1) {
        baseMaxGroup = (numSeqs > 0) ? (numSeqs + 1) : INT32_MAX;
        baseSampleRate = 0.5f;
        baseSampleTop = 1.0f;
    } else if (effectiveMode == 2) {
        baseMaxGroup = 200;
        baseSampleRate = 0.5f;
        baseSampleTop = 1.0f;
    } else if (effectiveMode == 3) {
        baseMaxGroup = 200;
        baseSampleRate = 0.2f;
        baseSampleTop = 0.3f;
    }

    this->accMaxGroup = this->hasUserAccMaxGroup ? this->userAccMaxGroup : baseMaxGroup;
    this->accSampleRate = this->hasUserAccSampleRate ? this->userAccSampleRate : baseSampleRate;
    this->accSampleTop = this->hasUserAccSampleTop ? this->userAccSampleTop : baseSampleTop;
    this->sampleMaxGroup = this->accMaxGroup;
    this->sampleRate = this->accSampleRate;
    this->sampleTopCentrality = this->accSampleTop;
}
#ifndef MSA_HPP
#include "msa.hpp"
#endif

#include <zlib.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>

namespace fs = boost::filesystem;

msa::Option::Option(po::variables_map& vm) {

    // Detect Alignment Mode
    bool hasTree = vm.count("tree"), hasSeq = vm.count("sequences"), hasFile = vm.count("files"), hasAln = vm.count("alignment");
    int8_t mask = (hasFile << 3) | (hasTree << 2) | (hasSeq << 1) | (hasAln << 0);
    switch (mask) {
        case 0b0110: this->alnMode = 0; break;
        case 0b1000: this->alnMode = 1; break;
        case 0b0011: this->alnMode = 2; break;
        case 0b0111: this->alnMode = 3; break;
        default:
            std::cerr << "ERROR: Unrecognized alignment mode based on the provided options.\n"
                      << "Valid combinations are:\n"
                      << "  [1] --tree and --sequences (for building MSA from unaligned sequences)\n"
                      << "  [2] --files (for merging MSAs)\n"
                      << "  [3] --sequences and --alignment (for adding new sequences to an existing MSA without guide tree)\n"
                      << "  [4] --sequences, --alignment and --tree (for adding new sequences to an existing MSA using a given guide tree)\n"
                      << "Please check the input arguments or run with --help for usage.\n";
            exit(1);
            break;
    }
    // Read Input File Names
    this->treeFile        = (vm.count("tree"))      ? vm["tree"].as<std::string>()      : "";
    this->seqFile         = (vm.count("sequences")) ? vm["sequences"].as<std::string>() : "";
    this->msaDir          = (vm.count("files"))     ? vm["files"].as<std::string>()     : "";
    this->backboneAlnFile = (vm.count("alignment")) ? vm["alignment"].as<std::string>() : "";

    // Read Options
    this->gpuNum = 0;
    this->gpuIdx = std::vector<int> (0);
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    this->cpuNum = (vm.count("cpu")) ? vm["cpu"].as<int>() : maxCpuThreads;
    if (this->cpuNum <= 0 || cpuNum > maxCpuThreads) {
        std::cerr << "ERROR: Invalid number of CPU cores. Please request between 1 and " << maxCpuThreads << " (got " << this->cpuNum << ").\n";
        exit(1);
    }

    this->maxSubtree = (vm.count("max-subtree")) ? vm["max-subtree"].as<int>() : INT32_MAX;
    if (this->maxSubtree <= 0) {
        std::cerr << "ERROR: Invalid value for --max-subtree. The value of --max-subtree should be a positive integer (got " << this->maxSubtree << ").\n";
        exit(1);
    }
    this->gappyVertical = vm["remove-gappy"].as<float>();
    if (this->gappyVertical > 1 || this->gappyVertical <= 0) {
        std::cerr << "ERROR: Invalid value for --remove-gappy. The value of --remove-gappy should be in (0,1]\n";
        exit(1);
    }
    this->lenDev = vm.count("length-deviation") ? vm["length-deviation"].as<float>() : 0;
    if (this->lenDev < 0) {
        std::cerr << "ERROR: Invalid value for --length-deviation. The value of --length-deviation should be non-negative\n";
        exit(1);
    }
    this->maxAmbig = vm["max-ambig"].as<float>();
    if (this->maxAmbig > 1 || this->maxAmbig <= 0) {
        std::cerr << "ERROR: Invalid value for --max-ambig. The value of --max-ambig should be in (0,1]\n";
        exit(1);
    }
    this->maxLen = vm.count("max-len") ? vm["max-len"].as<int>() : INT32_MAX;
    if (this->maxLen <= 0 && vm.count("max-len")) {
        std::cerr << "ERROR: Invalid value for --max-len. The value of --max-len should be larger than 0 (got " << this->maxLen << ").\n";
        exit(1);
    }
    this->minLen =  vm.count("min-len") ? vm["min-len"].as<int>() : 0;
    if (this->minLen <= 0 && vm.count("min-len")) {
        std::cerr << "ERROR: Invalid value for --min-len. The value of --min-len should be larger than 0 (got " << this->minLen << ").\n";
        exit(1);
    }
    if (this->minLen >= this->maxLen) {
        std::cerr << "ERROR: Invalid value for --min-len and --max-len. --max-len should be larger than --min-len (got --max-len: " << this->maxLen << ", --min-len: " << this->minLen << ").\n";
        exit(1);
    }

    if ((vm.count("min-len") || vm.count("max-len")) && vm.count("length-deviation")) {
        std::cerr << "ERROR: Invalid arguments. --length-deviation cannot be used together with --min-len or --max-len.\n";
        exit(1);
    }
    
    this->reroot = !vm.count("rooted");
    this->debug = vm.count("check");
    this->cpuOnly = vm.count("cpu-only");
    this->printDetail = vm.count("verbose");
    this->deleteTemp = !vm.count("keep-temp");
    this->alignGappy = !vm.count("no-align-gappy");
    this->compressed = vm.count("compress");
    this->noFilter = !vm.count("filter");
    this->writeFiltered = vm.count("write-filtered");

    // Parse accuracy-level
    this->accuracyLevel = 1;
    this->accurate = false;
    this->accurateMethod = "pairwise";

    if (vm.count("accuracy-level")) {
        std::string lvlStr = vm["accuracy-level"].as<std::string>();
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
            std::cerr << "ERROR: Invalid value for --accuracy-level: \"" << lvlStr 
                      << "\". Usage: --accuracy-level <1|2|3> (1: traditional sum-of-pairs scoring, 2: + pairwise consistency, 3: + transitive consistency).\n";
            exit(1);
        }
    }

    this->consistencyWeight = (vm.count("consistency-weight")) ? vm["consistency-weight"].as<float>() : 5.0f;

    // Parse sample-mode
    this->sampleMode = (vm.count("sample-mode")) ? vm["sample-mode"].as<int>() : 0;
    this->accMode = this->sampleMode;
    if (this->sampleMode < 0 || this->sampleMode > 3) {
        std::cerr << "ERROR: Invalid value for --sample-mode. Must be 0 (auto), 1 (all-to-all), 2 (high sample rate), or 3 (low sample rate) (got " << this->sampleMode << ").\n";
        exit(1);
    }

    this->hasUserAccMaxGroup = vm.count("sample-max-group");
    this->hasUserAccSampleRate = vm.count("sample-rate");
    this->hasUserAccSampleTop = vm.count("sample-top-centrality");

    if (this->hasUserAccMaxGroup) {
        this->userAccMaxGroup = vm["sample-max-group"].as<int>();
        if (this->userAccMaxGroup <= 1) {
            std::cerr << "ERROR: Invalid value for --sample-max-group. The value of --sample-max-group should be > 1 (got " << this->userAccMaxGroup << ").\n";
            exit(1);
        }
    }
    if (this->hasUserAccSampleRate) {
        this->userAccSampleRate = vm["sample-rate"].as<float>();
        if (this->userAccSampleRate <= 0.0f || this->userAccSampleRate >= 1.0f) {
            std::cerr << "ERROR: --sample-rate must be greater than 0 and less than 1 (got " << this->userAccSampleRate << ").\n";
            exit(1);
        }
    }
    if (this->hasUserAccSampleTop) {
        this->userAccSampleTop = vm["sample-top-centrality"].as<float>();
        if (this->userAccSampleTop <= 0.0f || this->userAccSampleTop > 1.0f) {
            std::cerr << "ERROR: --sample-top-centrality must be greater than 0 and less than or equal to 1 (got " << this->userAccSampleTop << ").\n";
            exit(1);
        }
    }

    // Initialize with 0 sequences first (will be dynamically resolved per dataset / subtree)
    this->resolveAccurateParams(0);

    if (this->accurate) {
        if (this->alnMode != DEFAULT_ALN && this->alnMode != ACCURATE_ALN && this->alnMode != PLACE_WO_TREE) {
            std::cerr << "ERROR: Accuracy levels 2 and 3 (consistency mode) are currently supported for guide-tree alignment and sequence placement without tree.\n";
            exit(1);
        }
        if (this->alnMode == DEFAULT_ALN) {
            this->alnMode = ACCURATE_ALN;
        }
    }

    // Detect Data Type
    if (vm.count("type")) {
        if (vm["type"].as<std::string>() != "n" && vm["type"].as<std::string>() != "p") {
            std::cerr << "ERROR: Unrecognized data type \"" << vm["type"].as<std::string>() << "\".\n";
            fs::remove_all(tempDir);
            exit(1);
        }
        else this->type = vm["type"].as<std::string>()[0];
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
        if (seqFile.substr(seqFile.size()-3, 3) == ".gz") {
            gzFile file = gzopen(seqFile.c_str(), "r");
            if (!file) {
                std::cerr << "ERROR: Failed to open file " << seqFile << ".\n";
                exit(1);
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
                std::cerr << "Error: cannot open " << seqFile << std::endl;
                fs::remove_all(tempDir);
                exit(1);
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

    // Check Output File
    if (!vm.count("output")) {
        std::cerr << "ERROR: Output file name is required.\n";
        exit(1);
    }
    this->outFile = vm["output"].as<std::string>();
    if (!vm.count("overwrite")) {
        std::string outFileStr = (this->compressed) ? this->outFile + ".gz" : this->outFile;
        if (fs::exists(outFileStr)) {
            std::cerr << "ERROR: " << outFileStr << " already exists. Please use another file name or add --overwrite to overwrite the existing file.\n";
            exit(1);
        }
    }
    if (this->compressed) {
        gzFile outFile = gzopen(this->outFile.c_str(), "wb"); // "wb" = write binary
        if (!outFile) {
            fprintf(stderr, "ERROR: failed to open file: %s\n", this->outFile.c_str());
            exit(1);
        }
        gzclose(outFile);
        fs::remove(this->outFile);
    }
    else {
        std::ofstream outFile(this->outFile);
        if (!outFile) {
            fprintf(stderr, "ERROR: failed to open file: %s\n", this->outFile.c_str());
            exit(1);
        }
        outFile.close();
        fs::remove(this->outFile);
    }
    

    // Create Temporary Directory if Needed
    if (this->maxSubtree < INT32_MAX || vm.count("files") || this->alnMode == 2) {
        std::string tempDir;
        if (!vm.count("temp-dir")) {
            int idx = 1;
            std::string tempDir_org = "./twilight_temp";
            tempDir = tempDir_org;
            while (true) {
                if (mkdir(tempDir.c_str(), 0777) == -1) {
                    if( errno == EEXIST ) {
                        tempDir = tempDir_org + '_' + std::to_string(idx);
                        ++idx;
                    }
                    else { fprintf(stderr, "ERROR: Can't create directory: %s\n", tempDir.c_str()); exit(1); }
                }
                else break;
            }
        }
        else {
            tempDir = vm["temp-dir"].as<std::string>();
            if (tempDir[tempDir.size()-1] == '/') tempDir = tempDir.substr(0, tempDir.size()-1);
            if (fs::exists(tempDir)) {
                if (!vm.count("overwrite")) {
                    std::cerr << "ERROR: " << tempDir << " already exists. In order to prevent your file from being overwritten, please delete this folder or use another folder name.\n";
                    exit(1);
                }
            }
            fs::create_directories(tempDir);
        }
        std::cout << tempDir << " created for storing temporary alignments\n";
        this->tempDir = tempDir;
    }

    
    std::cerr << "====== Configuration =======\n";
    if (this->accurate) {
        std::string methodDisplay = (this->accuracyLevel == 3) ? "3 (+ transitive consistency)" : "2 (+ pairwise consistency)";
        std::string modeStr = (this->sampleMode == 0) ? "0 (auto)" : (this->sampleMode == 1) ? "1 (all-to-all)" : (this->sampleMode == 2) ? "2 (high sample rate)" : "3 (low sample rate)";
        std::cerr << "Accuracy level: " << methodDisplay << " (sample-mode=" << modeStr;
        if (this->hasUserAccMaxGroup) std::cerr << ", sample-max-group=" << this->userAccMaxGroup << " [user]";
        if (this->hasUserAccSampleRate) std::cerr << ", sample-rate=" << this->userAccSampleRate << " [user]";
        if (this->hasUserAccSampleTop) std::cerr << ", sample-top-centrality=" << this->userAccSampleTop << " [user]";
        std::cerr << ", consistency weight=" << this->consistencyWeight << ")\n";
    } else {
        std::cerr << "Accuracy level: 1 (traditional sum-of-pairs scoring)\n";
    }
    if (this->maxSubtree != INT32_MAX) 
    std::cerr << "Max-subtree: " << this->maxSubtree << '\n';
    if (this->gappyVertical == 1) 
    std::cerr << "Disable removing gappy columns.\n";
    else
    std::cerr << "Threshold for removing gappy columns: " << this->gappyVertical << '\n';
    if (this->lenDev > 0)   std::cerr << "Allowed deviation from the median length: " << (this->lenDev * 100) << "%\n";
    else if (this->minLen > 0 || this->maxLen < INT32_MAX) std::cerr << "Allowed sequence length range: [" << this->minLen << ", " << this->maxLen << "]\n";
    if (this->maxAmbig < 1) std::cerr << "Allowed proportion of ambiguous characters: " << (this->maxAmbig * 100) << "%\n";
    fprintf(stderr, "Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);
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

    int baseMaxGroup = 100;
    float baseSampleRate = 0.5f;
    float baseSampleTop = 1.0f;

    if (effectiveMode == 1) {
        baseMaxGroup = (numSeqs > 0) ? (numSeqs + 1) : INT32_MAX;
        baseSampleRate = 0.5f;
        baseSampleTop = 1.0f;
    } else if (effectiveMode == 2) {
        baseMaxGroup = 100;
        baseSampleRate = 0.5f;
        baseSampleTop = 1.0f;
    } else if (effectiveMode == 3) {
        baseMaxGroup = 100;
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
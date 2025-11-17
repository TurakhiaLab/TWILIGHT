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
        std::cerr << "ERROR: Invalid number of CPU cores. Please request between 1 and " << maxCpuThreads << ".\n";
        exit(1);
    }

    this->maxSubtree = (vm.count("max-subtree")) ? vm["max-subtree"].as<int>() : INT32_MAX;
    if (this->maxSubtree <= 0) {
        std::cerr << "ERROR: Invalid value for --max-subtree. The value of --max-subtree should be a positive integer.\n";
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
    
    this->reroot = !vm.count("rooted");
    this->debug = vm.count("check");
    this->cpuOnly = vm.count("cpu-only");
    this->printDetail = vm.count("verbose");
    this->deleteTemp = !vm.count("keep-temp");
    this->alignGappy = !vm.count("no-align-gappy");
    this->compressed = vm.count("compress");
    this->noFilter = !vm.count("filter");

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
    if (this->maxSubtree != INT32_MAX) 
    std::cerr << "Max-subtree: " << this->maxSubtree << '\n';
    if (this->gappyVertical == 1) 
    std::cerr << "Disable removing gappy columns.\n";
    else
    std::cerr << "Threshold for removing gappy columns: " << this->gappyVertical << '\n';
    if (this->lenDev > 0)   std::cerr << "Allowed deviation from the average length: " << (this->lenDev * 100) << "%\n";
    if (this->maxAmbig < 1) std::cerr << "Allowed proportion of ambiguous characters: " << (this->maxAmbig * 100) << "%\n";
    fprintf(stderr, "Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);
}
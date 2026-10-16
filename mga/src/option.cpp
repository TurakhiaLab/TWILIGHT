#ifndef OPTION_HPP
#include "option.hpp"
#endif

#include "message_manager.hpp"
#include <zlib.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>
#include <algorithm>
#include <cctype>

namespace fs = boost::filesystem;

Option::Option(po::variables_map& vm) {

    // Read Input File Names
    this->treeFile        = (vm.count("tree"))      ? vm["tree"].as<std::string>()      : "";
    this->seqFile         = (vm.count("sequences")) ? vm["sequences"].as<std::string>() : "";

    // Read Logging & Verbosity Options
    this->verbose         = (vm.count("verbose"));
    this->quite           = (vm.count("quiet"));
    this->logFile         = (vm.count("log"))       ? vm["log"].as<std::string>()       : "";

    // Initialize mga::msg logging system based on command line options
    mga::msg::init(this->verbose, this->quite, this->logFile);

    // Read Hardware & Alignment Options
    this->gpuNum = 0;
    this->gpuIdx = std::vector<int> (0);
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    this->cpuNum = (vm.count("cpu")) ? vm["cpu"].as<int>() : maxCpuThreads;
    if (this->cpuNum <= 0 || cpuNum > maxCpuThreads) {
        mga::msg::error("Invalid number of CPU cores. Please request between 1 and %d (got %d).\n", maxCpuThreads, this->cpuNum);
    }

    if ((vm.count("min-len") || vm.count("max-len")) && vm.count("length-deviation")) {
        mga::msg::error("Invalid arguments. --length-deviation cannot be used together with --min-len or --max-len.\n");
    }
    mga::msg::info("Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);

    this->cpuOnly = (vm.count("cpu-only"));
    this->circular = (vm.count("circular"));
    this->writeOriented = (vm.count("write-oriented"));

    // Read Minimap2 ASM Preset Divergence Level (5, 10, 20)
    this->mmAsm = (vm.count("mm-asm")) ? vm["mm-asm"].as<int>() : 5;
    if (this->mmAsm != 5 && this->mmAsm != 10 && this->mmAsm != 20) {
        mga::msg::error("Invalid --mm-asm value: %d. Supported values are 5, 10, or 20 (corresponding to asm5, asm10, asm20).\n", this->mmAsm);
    }
    this->mmAsmPreset = "asm" + std::to_string(this->mmAsm);

    // Check Prefix (Required)
    if (!vm.count("prefix") || vm["prefix"].as<std::string>().empty()) {
        mga::msg::error("Prefix for output files is required. Please specify using -p or --prefix.\n");
    }
    this->prefix = vm["prefix"].as<std::string>();

    // Check Output Format (maf or panman, default: output both)
    if (vm.count("output")) {
        std::string out = vm["output"].as<std::string>();
        std::string outLower = out;
        std::transform(outLower.begin(), outLower.end(), outLower.begin(), [](unsigned char c) {
            return std::tolower(c);
        });
        if (outLower == "maf") {
            this->outputMAF = true;
            this->outputPanman = false;
            this->outputFormat = "maf";
        } else if (outLower == "panman") {
            this->outputMAF = false;
            this->outputPanman = true;
            this->outputFormat = "panman";
        } else {
            mga::msg::error("Invalid output format: '%s'. Supported formats are 'maf' or 'panman'.\n", out.c_str());
        }
    } else {
        this->outputMAF = true;
        this->outputPanman = true;
        this->outputFormat = "both";
    }

    // Check and create working directory
    std::string workDir;
    if (!vm.count("dir")) {
        int idx = 1;
        std::string baseDir = "twilight_working_dir_";
        while (true) {
            std::string candidate = baseDir + std::to_string(idx);
            if (mkdir(candidate.c_str(), 0777) == 0) {
                workDir = candidate;
                break;
            } else if (errno == EEXIST) {
                ++idx;
            } else {
                mga::msg::error("Can't create directory: %s\n", candidate.c_str());
            }
        }
    } else {
        workDir = vm["dir"].as<std::string>();
        while (workDir.size() > 1 && workDir.back() == '/') {
            workDir.pop_back();
        }
        if (fs::exists(workDir)) {
            if (!vm.count("overwrite")) {
                mga::msg::error("Directory '%s' already exists. In order to prevent your file from being overwritten, please specify another folder name or use --overwrite to remove the current directory.\n", workDir.c_str());
            }
            boost::system::error_code ec;
            fs::remove_all(workDir, ec);
            if (ec) {
                mga::msg::error("Failed to remove directory '%s': %s\n", workDir.c_str(), ec.message().c_str());
            }
        }
        boost::system::error_code ec;
        if (!fs::create_directories(workDir, ec)) {
            mga::msg::error("Can't create directory: %s\n", workDir.c_str());
        }
    }
    mga::msg::info("%s created for storing working and output files.\n", workDir.c_str());
    this->dir = workDir;
    this->tempDir = workDir;
    
    this->minimap2_time = 0;
    this->refineBlockSet_time = 0;
    this->refineBlock_time = 0;
    this->merge_time = 0;
    this->debug_time = 0;
    
}
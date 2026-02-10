#ifndef OPTION_HPP
#include "option.hpp"
#endif

#include <zlib.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <tbb/parallel_for.h>

namespace fs = boost::filesystem;

Option::Option(po::variables_map& vm) {

    // Read Input File Names
    this->treeFile        = (vm.count("tree"))      ? vm["tree"].as<std::string>()      : "";
    this->seqFile         = (vm.count("sequences")) ? vm["sequences"].as<std::string>() : "";

    // Read Options
    this->gpuNum = 0;
    this->gpuIdx = std::vector<int> (0);
    int maxCpuThreads = tbb::this_task_arena::max_concurrency();
    this->cpuNum = (vm.count("cpu")) ? vm["cpu"].as<int>() : maxCpuThreads;
    if (this->cpuNum <= 0 || cpuNum > maxCpuThreads) {
        std::cerr << "ERROR: Invalid number of CPU cores. Please request between 1 and " << maxCpuThreads << " (got " << this->cpuNum << ").\n";
        exit(1);
    }


    if ((vm.count("min-len") || vm.count("max-len")) && vm.count("length-deviation")) {
        std::cerr << "ERROR: Invalid arguments. --length-deviation cannot be used together with --min-len or --max-len.\n";
        exit(1);
    }
    fprintf(stderr, "Maximum available CPU cores: %d. Using %d CPU cores.\n", maxCpuThreads, cpuNum);

    this->verbose = (vm.count("verbose"));
    this->quite = (vm.count("quiet"));
    this->cpuOnly = (vm.count("cpu-only"));
    
}
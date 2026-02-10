#ifndef OPTION_HPP
#define OPTION_HPP



#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <sys/stat.h>
#include <zlib.h>

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

char checkOnly(char inChar);
int letterIdx(char type, char inChar);

struct Option {
    // Hardware Options
    int gpuNum;
    int cpuNum;
    std::vector<int> gpuIdx;
    std::vector<size_t> gpuMem;
    bool cpuOnly;
    bool verbose;
    bool quite;
    // File Names
    std::string treeFile;
    std::string seqFile;
    std::string outFile;
    // Constructor
    Option(po::variables_map &vm);
    // Only used in GPU version
    void getGpuInfo(po::variables_map &vm);
};


#endif
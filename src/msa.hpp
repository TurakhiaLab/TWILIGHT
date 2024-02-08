#ifndef MSA_HPP
#define MSA_HPP

#include <string>
#include <unordered_map>

namespace msa
{
    struct utility
    {
        std::unordered_map<std::string, std::string> seqs;
    };
} 

#endif
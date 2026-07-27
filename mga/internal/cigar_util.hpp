#ifndef CIGAR_UTIL_HPP
#define CIGAR_UTIL_HPP

#include <string>
#include <algorithm>
#include <iostream>
#include <cctype>

#include "type.hpp"

// =========================
// CIGAR & Sequence Utilities
// =========================

namespace parser {
inline CigarString parseCigar(const std::string& cigarString) {
    CigarString ops;
    size_t i = 0;
    while (i < cigarString.length()) {
        size_t j = i;
        while (j < cigarString.length() && std::isdigit(cigarString[j])) {
            j++;
        }
        if (i == j) break;

        int length = std::stoi(cigarString.substr(i, j - i));
        char op = cigarString[j];
        ops.push_back({length, op});
        i = j + 1;
    }
    return ops;
}
} // namespace parser


inline bool consumesRef(char op) {
    return (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X');
}

inline bool consumesQry(char op) {
    return (op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X');
}

inline char complement(char base) {
    switch (base) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a';
        case 'c': return 'g'; case 'g': return 'c';
        default: return base;
    }
}

inline std::string getReverseComplement(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        c = complement(c);
    }
    return seq;
}

inline CigarString compressCigar(const CigarString& cigar) {
    CigarString compressed;
    for (auto& op : cigar) {
        if (!compressed.empty() && compressed.back().second == op.second) {
            compressed.back().first += op.first;
        } else if (op.first > 0) {
            compressed.push_back(op);
        }
    }
    return compressed;
}

inline std::ostream& operator<<(std::ostream& os, const CigarString& cigar) {
    for (const auto& op : cigar) {
        os << op.first << op.second;
    }
    return os;
}

inline void printCIGAR(const CigarString& cigar, bool changeLine = true) {
    std::cout << cigar;
    if (changeLine) std::cout << std::endl;
}

inline String cigarToStr(const CigarString& c) {
    std::string s = "";
    for (auto& op : c) s += std::to_string(op.first) + op.second;
    return s;
}

inline int scoreCIGAR(const CigarString& c) {
    const int penalty = -2;
    const int bonus = 1;
    int score = 0;
    for (auto& op : c) score += ((op.second == 'M') ? op.first * bonus : op.first * penalty);
    return score;
}

inline double identityCIGAR(const CigarString& c) {
    int match_count = 0;
    int total_aln_len = 0;
    for (const auto& op : c) {
        if (op.second == 'M' || op.second == '=') {
            match_count += op.first;
            total_aln_len += op.first;
        } else if (op.second == 'X' || op.second == 'I' || op.second == 'D') {
            total_aln_len += op.first;
        }
    }
    double identity = total_aln_len > 0 ? (double)match_count / total_aln_len : 0.0;
    return identity;
}

inline int mapCoordinate(const CigarString& cigar, int q_start, int r_start, int target, bool is_target_ref) {
    int source_start = is_target_ref ? r_start : q_start;
    int dest_start   = is_target_ref ? q_start : r_start;

    if (target <= source_start) return dest_start;

    int curr_q = q_start;
    int curr_r = r_start;

    for (const auto& op_pair : cigar) {
        int len = op_pair.first;
        char op = op_pair.second;

        bool cons_q = consumesQry(op);
        bool cons_r = consumesRef(op);

        int next_q = curr_q + (cons_q ? len : 0);
        int next_r = curr_r + (cons_r ? len : 0);

        bool cons_source = is_target_ref ? cons_r : cons_q;
        bool cons_dest   = is_target_ref ? cons_q : cons_r;
        int curr_source  = is_target_ref ? curr_r : curr_q;
        int curr_dest    = is_target_ref ? curr_q : curr_r;
        int next_source  = is_target_ref ? next_r : next_q;

        if (cons_source && target <= next_source) {
            int offset = target - curr_source;
            return cons_dest ? (curr_dest + offset) : curr_dest;
        }

        curr_q = next_q;
        curr_r = next_r;
    }
    
    return is_target_ref ? curr_q : curr_r;
}

inline int getCigarLength(const CigarString& cigar) {
    int len = 0;
    for (const auto& op : cigar) {
        len += op.first;
    }
    return len;
}

#endif // CIGAR_UTIL_HPP

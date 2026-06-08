#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "type.hpp"

const int INF = -1e9;

const int MATCH_SCORE = 1;
const int MISMATCH_PENALTY = -4;
const int GAP_OPEN = -10;
const int GAP_EXT = -2;

// ==========================================
// 1. Global Alignment with Affine Gap
// ==========================================
CigarString runGlobalAlignment(const std::string& ref, const std::string& qry); 

// ==========================================
// 2. Semi-Global Alignment with Affine Gap
// ==========================================
CigarString runSemiGlobalAlignment(const std::string& ref, const std::string& qry); 


#endif // GLOBAL_ALIGNMENT_H
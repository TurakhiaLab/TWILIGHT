#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "type.hpp"
#include "consensus.hpp"

const int INF = -1e9;

const int MATCH_SCORE = 1;
const int MISMATCH_PENALTY = -4;
const int GAP_OPEN = -10;
const int GAP_EXT = -2;
const int GAP_PENALTY = -4;

// ==========================================
// 1. Global Alignment with Affine Gap
// ==========================================
CigarString runGlobalAlignment(const std::string& ref, const std::string& qry); 
CigarString runGlobalAlignment(const Consensus& refCons, const Consensus& qryCons);

// ==========================================
// 1b. Global Alignment with Linear Gap (No Affine)
// ==========================================
CigarString runGlobalAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty = GAP_PENALTY); 
CigarString runGlobalAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty = GAP_PENALTY);

// ==========================================
// 2. Semi-Global Alignment with Affine Gap
// ==========================================
CigarString runSemiGlobalAlignment(const std::string& ref, const std::string& qry); 
CigarString runSemiGlobalAlignment(const Consensus& refCons, const Consensus& qryCons);

// ==========================================
// 2b. Semi-Global Alignment with Linear Gap (No Affine)
// ==========================================
CigarString runSemiGlobalAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty = GAP_PENALTY); 
CigarString runSemiGlobalAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty = GAP_PENALTY);

// ==========================================
// 3. Tiling Alignment (GACT)
// ==========================================
CigarString runTilingAlignment(const std::string& ref, const std::string& qry);
CigarString runTilingAlignment(const Consensus& refCons, const Consensus& qryCons);

// ==========================================
// 3b. Tiling Alignment with Linear Gap (No Affine)
// ==========================================
CigarString runTilingAlignmentLinearGap(const std::string& ref, const std::string& qry, int gapPenalty = GAP_PENALTY);
CigarString runTilingAlignmentLinearGap(const Consensus& refCons, const Consensus& qryCons, int gapPenalty = GAP_PENALTY);


#endif // GLOBAL_ALIGNMENT_H
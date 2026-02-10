#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

// Helper to parse CIGAR string into vector of ops
mga::Cigar mga::parser::parseCigar(const std::string& cigarString) {
    Cigar ops;
    size_t i = 0;
    while (i < cigarString.length()) {
        size_t j = i;
        // Get number
        while (j < cigarString.length() && std::isdigit(cigarString[j])) {
            j++;
        }
        if (i == j) break; // Safety break

        int length = std::stoi(cigarString.substr(i, j - i));
        char op = cigarString[j];
        ops.push_back({length, op});
        i = j + 1;
    }
    return ops;
}

mga::alnVec mga::parser::parseMinimap2(const std::string& filename, bool isSAM) {
    std::vector<Alignment> alignments;
    std::ifstream file(filename);

    // Scoring constants
    const float alpha = 1.0f;
    const float beta = 1.0f;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return alignments;
    }

    std::string line;
    int alignmentCounter = 0;

    while (std::getline(file, line)) {
        // Skip header lines (SAM usually starts with @)
        if (line.empty() || line[0] == '@') continue;

        std::stringstream ss(line);
        std::string segment;
        std::vector<std::string> fields;

        // Split line by tabs
        while (std::getline(ss, segment, '\t')) fields.push_back(segment);

        // Basic validation
        size_t minFields = isSAM ? 11 : 12;
        if (fields.size() < minFields) continue;

        Alignment aln;
        aln.identifier = alignmentCounter++;
        aln.used = false;
        aln.chainScore = 0;
        aln.primary = false;

        std::string cigarString;
        int nm_val = 0;
        int pafRefStart = 0, pafRefEnd = 0, pafQryStart = 0, pafQryEnd = 0;

        // ==========================================
        // 1. EXTRACT RAW FIELDS BASED ON FORMAT
        // ==========================================
        if (isSAM) {
            // --- SAM FORMAT ---
            // Flag (Field 1)
            int flag = std::stoi(fields[1]);
            aln.inverse = (flag & 16) != 0; // 0x10 = Reverse Strand

            // Ref Start (Field 3) - SAM is 1-based, convert to 0-based
            aln.refIdx.first = std::stoi(fields[3]) - 1;

            // CIGAR (Field 5)
            cigarString = fields[5];

            // Parse Tags (starting at Field 11) for NM
            for (size_t k = 11; k < fields.size(); ++k) {
                if (fields[k].rfind("NM:i:", 0) == 0) {
                    try { nm_val = std::stoi(fields[k].substr(5)); } catch (...) { nm_val = 0; }
                    break; 
                }
            }

        } 
        else {
            // --- PAF FORMAT ---
            // Query Coords (Fields 2, 3) - 0-based
            pafQryStart = std::stoi(fields[2]);
            pafQryEnd = std::stoi(fields[3]);

            // Orientation (Field 4)
            aln.inverse = (fields[4] == "-");

            // Ref Coords (Fields 7, 8) - 0-based
            pafRefStart = std::stoi(fields[7]);
            pafRefEnd   = std::stoi(fields[8]);

            // Set explicit PAF coords immediately
            aln.refIdx.first  = pafRefStart;
            aln.refIdx.second = pafRefEnd;
            aln.qryIdx.first  = pafQryStart;
            aln.qryIdx.second = pafQryEnd;

            // Parse Tags (starting at Field 12) for cg:Z (CIGAR) and NM:i
            for (size_t k = 12; k < fields.size(); ++k) {
                if (fields[k].rfind("cg:Z:", 0) == 0) {
                    cigarString = fields[k].substr(5);
                } else if (fields[k].rfind("NM:i:", 0) == 0) {
                    try { nm_val = std::stoi(fields[k].substr(5)); } catch (...) { nm_val = 0; }
                }
            }
        }

        aln.CIGAR = parseCigar(cigarString);
        if (aln.CIGAR.empty()) continue;

        // ==========================================
        // 2. UNIFIED CIGAR PROCESSING & SCORING
        // ==========================================
        // We use the SAM logic here to calculate lengths and offsets.
        // For SAM, this is required to determine End coords and Qry Start.
        // For PAF, this is used for Scoring (totLen, totIndel).

        int currentRefLen = 0;
        int currentQryLen = 0;
        int samQryOffset = 0; // Calculated from initial S/H clips
        bool firstOp = true;

        int totLen = 0;    // Alignment length (M + I + D + X + =)
        int totIndel = 0;  // Total Indel length (I + D)

        for (auto op : aln.CIGAR) {
            int length = op.first;
            char opType = op.second;

            switch (opType) {
                case 'M': 
                case 'X':
                case '=':
                    currentRefLen += length;
                    currentQryLen += length;
                    totLen += length; 
                    firstOp = false;
                    break;
                case 'D': // Deletion from ref
                case 'N': // Skipped region
                    currentRefLen += length;
                    totLen += length;
                    totIndel += length;
                    firstOp = false;
                    break;
                case 'I': // Insertion to ref
                    currentQryLen += length;
                    totLen += length;
                    totIndel += length;
                    firstOp = false;
                    break;
                case 'S': // Soft clipping
                    if (firstOp) samQryOffset += length;
                    break;
                case 'H': // Hard clipping
                    if (firstOp) samQryOffset += length;
                    break;
                default:
                    break;
            }
            
            // If we hit a 'real' alignment op, stop counting initial clips
            if (opType != 'S' && opType != 'H') {
                firstOp = false;
            }
        }

        // ==========================================
        // 3. FINALIZE COORDINATES & CLEANUP
        // ==========================================
        
        // Remove starting and ending H/S from vector
        if (!aln.CIGAR.empty() && (aln.CIGAR.front().second == 'H' || aln.CIGAR.front().second == 'S')) 
            aln.CIGAR.erase(aln.CIGAR.begin());
        if (!aln.CIGAR.empty() && (aln.CIGAR.back().second == 'H'  || aln.CIGAR.back().second == 'S'))  
            aln.CIGAR.pop_back();
        
        if (aln.CIGAR.empty()) continue;

        // For SAM, we must derive the missing coordinates now
        if (isSAM) {
            aln.refIdx.second = aln.refIdx.first + currentRefLen;
            aln.qryIdx.first  = samQryOffset;
            aln.qryIdx.second = samQryOffset + currentQryLen;
            
            // Debug output for inverse SAM (from your original code)
            if (aln.inverse) {
                // std::cout << aln.qryIdx.first << ',' << aln.qryIdx.second << '\n';
            }
        }

        // Calculate Final Score
        // Score = totLen - alpha * NM - beta * totIndel
        aln.alnScore = (float)totLen - (alpha * nm_val) - (beta * totIndel);

        alignments.push_back(aln);

        if (aln.inverse) std::cerr << "Inverse alignment: " <<aln.qryIdx.first << ',' << aln.qryIdx.second << '\n';
    }

    // Sort alignments: Group by RefStart, then by QryStart
    std::sort(alignments.begin(), alignments.end(), [&](const auto &a, const auto &b) {
        if (a.refIdx.first == b.refIdx.first) return a.qryIdx.first < b.qryIdx.first;
        return (a.refIdx.first < b.refIdx.first); 
    });

    return alignments;
}
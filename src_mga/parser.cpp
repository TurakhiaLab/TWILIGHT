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

mga::alnVec mga::parser::parseMinimap2PAF(const std::string& filename) {
    std::vector<Alignment> alignments;
    std::ifstream file(filename);

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

        Alignment aln;
        aln.identifier = alignmentCounter++;
        aln.used = false;
        aln.chainScore = 0;

        std::string cigarString;
        int numMismatch = 0, numIndel = 0, numEdits = 0, score = 0;
        int pafRefStart = 0, pafRefEnd = 0, pafQryStart = 0, pafQryEnd = 0;

        aln.qryName = fields[0]; // [NEW] Query sequence name
        aln.refName = fields[5]; // [NEW] Target sequence name
        // Query Coords (Fields 2, 3) - 0-based
        pafQryStart = std::stoi(fields[2]);
        pafQryEnd = std::stoi(fields[3]);

        // Orientation (Field 4)
        aln.inverse = (fields[4] == "-");

        // Ref Coords (Fields 7, 8) - 0-based
        pafRefStart = std::stoi(fields[7]);
        pafRefEnd   = std::stoi(fields[8]);

        // Set coords 
        aln.refIdx = {pafRefStart, pafRefEnd};
        aln.qryIdx = {pafQryStart, pafQryEnd};
    
        // Parse Tags (starting at Field 12) for cg:Z (CIGAR) and NM:i
        for (size_t k = 12; k < fields.size(); ++k) {
            if (fields[k].rfind("cg:Z:", 0) == 0) {
                cigarString = fields[k].substr(5);
            } 
            else if (fields[k].rfind("AS:i:", 0) == 0) {
                score = std::stoi(fields[k].substr(5));
            }
            else if (fields[k].rfind("NM:i:", 0) == 0) {
                numEdits = std::stoi(fields[k].substr(5));
            }
        }


        aln.CIGAR = parseCigar(cigarString);
        if (aln.CIGAR.empty()) continue;

        int totLen = 0, ins = 0, del = 0;

        for (auto op : aln.CIGAR) {
            int length = op.first;
            char opType = op.second;

            switch (opType) {
                case 'M': 
                case 'X':
                case '=':
                    totLen += length; 
                    break;
                case 'D': // Deletion from ref
                case 'N': // Skipped region
                    totLen += length;
                    del += length;
                    break;
                case 'I': // Insertion to ref
                    totLen += length;
                    ins += length;
                    break;
                default:
                    break;
            }
        }

        aln.alnScore = score;
        aln.ins = ins;
        aln.del = del;
        aln.mis = (numEdits - ins - del);

        alignments.push_back(aln);

        // if (aln.inverse) std::cerr << "Inverse alignment: " <<aln.qryIdx.first << ',' << aln.qryIdx.second << '\n';
    }

    // Sort alignments: Group by RefStart, then by QryStart
    std::sort(alignments.begin(), alignments.end(), [&](const auto &a, const auto &b) {
        if (a.refIdx.first == b.refIdx.first) return a.qryIdx.first < b.qryIdx.first;
        return (a.refIdx.first < b.refIdx.first); 
    });

    return alignments;
}
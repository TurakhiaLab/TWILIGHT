#ifndef MGA_HPP
#include "mga.hpp"
#endif

#include <list>
#include <fstream>

void printCoordinate(int st, int en) {
    std::cerr << "(" << st << "," << en << "]";
}

void printCoordinate(std::pair<int,int> r) {
    std::cerr << "(" << r.first << "," << r.second << "]";
}

// Helper to check if a CIGAR op consumes Reference
bool consumesRef(char op) {
    return (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X');
}

// Helper to check if a CIGAR op consumes Query
bool consumesQry(char op) {
    return (op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X');
}

mga::AlnChain::AlnChain(int id, float score, IntVec& chain) {
    this->identifier = id;
    this->score = score;
    this->chainedAln = chain;
}

void mga::Alignment::show() {
    std::string cigarString = "";
    for (auto op: this->CIGAR) cigarString += (std::to_string(op.first) + op.second);
    std::cerr << "Alignment{id=" << this->identifier
              << ", ref=[" << this->refIdx.first << "," << this->refIdx.second << ")"
              << ", qry=[" << this->qryIdx.first << "," << this->qryIdx.second << ")"
              << ", inverse=" << this->inverse
              << ", CIGAR=" << cigarString
              << ", alnScore=" << this->alnScore << "\n";
    return;
}

int mga::Alignment::mapCoordinate(int pos, bool inputIsRef) {
        
    // 1. Bounds Check & Target Offset Calculation
    int targetOffset = 0;
    if (inputIsRef) {
        // Check Ref Bounds
        if (pos < refIdx.first || pos >= refIdx.second) return -1;
        targetOffset = pos - refIdx.first;
    } 
    else {
        // Check Query Bounds
        if (pos < qryIdx.first || pos >= qryIdx.second) return -1;
        
        // Calculate Offset based on Strand
        if (!inverse) {
            // Forward: Simple offset from start
            targetOffset = pos - qryIdx.first;
        } else {
            // Inverse: Offset is calculated from the END moving backwards
            // (CIGAR start aligns with Query End on reverse strand)
            targetOffset = (qryIdx.second - 1) - pos; 
        }
    }
    // 2. Iterate CIGAR to find the operation containing 'targetOffset'
    int currR = 0;
    int currQ = 0;
    for (const auto& op : CIGAR) {
        int len = op.first;
        char type = op.second;
        // Does this op consume Reference / Query?
        bool cRef = consumesRef(type);
        bool cQry = consumesQry(type);
        // Select the tracking variable relevant to the Input Axis
        int currInput = inputIsRef ? currR : currQ;
        int stepInput = inputIsRef ? (cRef ? len : 0) : (cQry ? len : 0);
        // Check if our target lies within this operation
        if (currInput + stepInput > targetOffset) {
            // Found the segment! Now check if it actually maps to the other axis.
            // Case A: Input is Ref, outputting Query
            // We need the op to consume Query (e.g., M, =, X). If it's D (deletion), return -1.
            if (inputIsRef && !cQry) return -1; 
            // Case B: Input is Query, outputting Ref
            // We need the op to consume Ref. If it's I (insertion), return -1.
            if (!inputIsRef && !cRef) return -1;
            // 3. Calculate Final Result
            int internalOffset = targetOffset - currInput;
            if (inputIsRef) {
                // Outputting Query Coordinate
                if (!inverse) return qryIdx.first + currQ + internalOffset;
                else          return (qryIdx.second - 1) - (currQ + internalOffset);
            } else {
                // Outputting Reference Coordinate (Always Forward)
                return refIdx.first + currR + internalOffset;
            }
        }
        // Advance trackers
        currR += (cRef ? len : 0);
        currQ += (cQry ? len : 0);
    }
    return -1; // Fallback (should typically be caught by bounds check)
}

mga::chainVec mga::getAlignmentChains(alnVec& alignments) {
    chainVec chains;

    if (alignments.empty()) return chains;

    // Sort alignments by Strand and Reference Start 
    std::sort(alignments.begin(), alignments.end(), [](const Alignment& a, const Alignment& b) {
        if (a.inverse == b.inverse) return a.refIdx.second < b.refIdx.second;
        return (b.inverse);
    });

    int n = alignments.size();
    std::vector<int> parents(n, -1);
    float w_avg = 0;
    // Constants for Chaining
    const int MAX_DIST = 10000; // Max distance to look for a chain (heuristic)
    const int MAX_SKIP = 50;     // Max number of predecessors to check (heuristic)
    
    // Gap Penalties (Minimap2)
    std::function<float(const Alignment, const Alignment)> alpha = [&](const Alignment i, const Alignment j) -> float {
        float a = std::min((i.qryIdx.second-j.qryIdx.second),(i.refIdx.second-j.refIdx.second));
        float b = i.alnScore;
        return std::min(a, b);
    };
    std::function<float(const Alignment, const Alignment)> beta = [&](const Alignment i, const Alignment j) -> float {
        std::function<float(const int l)> gamma_c = [&](const int l) -> float {
            if (l != 0) return (0.01 * w_avg * std::abs(l) + 0.5 * std::log2(abs(l)));
            return 0;
        };
        int l = (i.qryIdx.second-j.qryIdx.second) - (i.refIdx.second-j.refIdx.second);
        return  gamma_c(l);
    };
    
    // Dynamic Programming
    for (int i = 0; i < n; ++i) {
        // Initialize
        alignments[i].chainScore = alignments[i].alnScore;
        alignments[i].primary = false; // Reset
        alignments[i].used = false;

        int skip_count = 0;

        // Look backwards for the best predecessor i
        for (int j = i - 1; j >= 0; --j) {
            if (skip_count > MAX_SKIP) break;
            // If reference distance is too large
            if (alignments[i].refIdx.second - alignments[j].refIdx.second > MAX_DIST) break;
            // If not on the same strand
            if (alignments[j].inverse != alignments[i].inverse) continue;
            // If overlapping
            int d_ref = alignments[i].refIdx.first - alignments[j].refIdx.second;
            int d_qry = (!alignments[i].inverse) ? alignments[i].qryIdx.first - alignments[j].qryIdx.second : 
                                                   alignments[j].qryIdx.first - alignments[i].qryIdx.second;
            if (d_ref < 0 || d_qry < 0) continue; 
            // Calculate Chain Score
            int gap_len = std::min(d_ref, d_qry);
            if (gap_len > MAX_DIST) continue; // Gap too big

            float chain_score = alignments[j].chainScore + alpha(alignments[i],alignments[j]) - beta(alignments[i],alignments[j]);
            float new_score = std::max(alignments[i].alnScore, chain_score);

            if (new_score > alignments[i].chainScore) {
                alignments[i].chainScore = new_score;
                parents[i] = j;
            }
            skip_count++;
        }
    }

    

    while (true) {
        // Find Global Maximum Score
        int id = chains.size();
        int best_idx = -1;
        float max_score = -1e9;
        for (int k = 0; k < n; ++k) {
            if (!alignments[k].used) {
                if (alignments[k].chainScore > max_score) {
                    max_score = alignments[k].chainScore;
                    best_idx = k;
                }
            }
        }
        if (best_idx == -1) break;
        IntVec aln;
        int curr = best_idx;
        while (curr != -1) {
            alignments[curr].used = true;
            aln.push_back(curr);
            curr = parents[curr];
            if (alignments[curr].used) break;
        }
        AlnChain chain (id, max_score, aln);
        chains.push_back(chain);
    }

    return chains;
}

void mga::identifyPrimaryAlignments(alnVec& alignments, chainVec& chains) {
    
    const int MAX_OVERLAP = 200; // Threshold: Discard if overlap > 200, Trim otherwise

    IntVec primaryIdx;

    // Iterate Chains (Highest Score First)
    for (const auto& chain : chains) {
        
        for (int alnId : chain.chainedAln) {
            Alignment& aln = alignments[alnId];
            
            bool discard = false;
            
            // We track the maximum overlap found on the Reference axis
            int maxHeadRefOv = 0; 
            int maxTailRefOv = 0;

            for (auto pIdx: primaryIdx) {
                const auto& pAln = alignments[pIdx];

                // 1. Calculate Reference Overlap
                int rStart = std::max(aln.refIdx.first, pAln.refIdx.first);
                int rEnd = std::min(aln.refIdx.second, pAln.refIdx.second);
                int rOv = std::max(0, rEnd - rStart);

                // 2. Calculate Query Overlap
                // (Used primarily for the Discard decision)
                int qStart = std::max(aln.qryIdx.first, pAln.qryIdx.first);
                int qEnd = std::min(aln.qryIdx.second, pAln.qryIdx.second);
                int qOv = std::max(0, qEnd - qStart);

                // --- Discard Logic ---
                if (rOv > MAX_OVERLAP || qOv > MAX_OVERLAP) {
                    discard = true;
                    break;
                }

                // --- Accumulate Trim Requirements ---
                // We use Reference coordinates as the "Anchor" for trimming.
                if (rOv > 0) {
                    // If P starts before Aln, it overlaps Aln's Head
                    if (pAln.refIdx.first <= aln.refIdx.first) {
                        maxHeadRefOv = std::max(maxHeadRefOv, rOv);
                    } else {
                        // P overlaps Aln's Tail
                        maxTailRefOv = std::max(maxTailRefOv, rOv);
                    }
                }
            }

            if (!discard) {
                // Determine the "Keep Range" (Absolute Reference Coordinates)
                int keepStart = aln.refIdx.first + maxHeadRefOv;
                int keepEnd   = aln.refIdx.second - maxTailRefOv;

                // Check if we actually need to trim
                if (maxHeadRefOv > 0 || maxTailRefOv > 0) {
                    
                    // Validate range before splitting
                    if (keepEnd > keepStart) {
                        // REUSE SPLIT FUNCTION:
                        // We want the segment between keepStart and keepEnd.
                        // splitAlignment returns {Head, Middle, Tail}. We want Middle.
                        SplitResult parts = splitAlignment(aln, keepStart, keepEnd, true); // true = onRef

                        if (parts.middle.valid) {
                            std::cerr << "Trimmed: " << alnId << '\n';
                            printCoordinate(aln.refIdx.first, aln.refIdx.second);
                            std::cerr << "\t"; 
                            printCoordinate(parts.middle.refIdx.first, parts.middle.refIdx.second);
                            std::cerr << "\n"; 
                            // Update the alignment with the trimmed version
                            aln = parts.middle;
                            aln.primary = true;
                            primaryIdx.push_back(alnId);
                        }
                    } 
                    // If keepEnd <= keepStart, the overlap consumed the whole alignment -> Do nothing (implicitly discarded)
                } 
                else {
                    // No overlap found, keep original
                    aln.primary = true;
                    primaryIdx.push_back(alnId);
                }
            }
        }
    }

    // Debug Output
    int pCount = 0;
    for (const auto& a: alignments) { if (a.primary) pCount++; }
    std::cout << pCount << '\\' << (alignments.size()-pCount) << '\n';
}


// Trims the alignment from the front (Head)
// targetRefTrim: minimum Ref bases to remove
// targetQryTrim: minimum Query bases to remove
void mga::trimHead(Alignment& aln, int targetRefTrim, int targetQryTrim) {
    int rRemoved = 0;
    int qRemoved = 0;
    int opIdx = 0;
    // Iterate CIGAR to find cut point
    while (opIdx < aln.CIGAR.size() && (rRemoved < targetRefTrim || qRemoved < targetQryTrim)) {
        auto& op = aln.CIGAR[opIdx];
        int len = op.first;
        char type = op.second;

        bool refC = consumesRef(type);
        bool qryC = consumesQry(type);

        int rStep = refC ? len : 0;
        int qStep = qryC ? len : 0;

        // Check if this full operation is within the trim range
        // We continue if we haven't met EITHER the ref requirement OR the qry requirement
        if ((rRemoved + rStep <= targetRefTrim) && (qRemoved + qStep <= targetQryTrim)) {
            // Remove whole op
            rRemoved += rStep;
            qRemoved += qStep;
            opIdx++; 
        } else {
            // Partial trim required on this op
            // Calculate how much we NEED to take to satisfy the deepest trim requirement
            int remainingRefNeeded = std::max(0, targetRefTrim - rRemoved);
            int remainingQryNeeded = std::max(0, targetQryTrim - qRemoved);

            int cutLen = 0;
            
            // Logic: We must cut enough to satisfy the larger requirement
            if (refC && qryC) {
                cutLen = std::max(remainingRefNeeded, remainingQryNeeded);
            } else if (refC) {
                cutLen = remainingRefNeeded;
            } else if (qryC) {
                cutLen = remainingQryNeeded;
            }

            // Apply partial cut
            op.first -= cutLen;
            rRemoved += (refC ? cutLen : 0);
            qRemoved += (qryC ? cutLen : 0);
            break; // Done
        }
    }

    // Remove fully consumed ops
    if (opIdx > 0) {
        aln.CIGAR.erase(aln.CIGAR.begin(), aln.CIGAR.begin() + opIdx);
    }

    // Update Coordinates
    aln.refIdx.first += rRemoved;
    
    // Update Query Coordinates based on strand
    if (!aln.inverse) {
        // Forward: Head trim increases start index
        aln.qryIdx.first += qRemoved;
    } else {
        // Reverse: Head trim (in CIGAR/Ref space) corresponds to the END of the query segment physically
        // Assuming qryIdx is [start, end) on the forward strand.
        // If mapped to reverse, CIGAR start aligns with qryIdx.second (high coord) moving down.
        aln.qryIdx.second -= qRemoved;
    }
}

// Trims the alignment from the back (Tail)
void mga::trimTail(Alignment& aln, int targetRefTrim, int targetQryTrim) {
    int rRemoved = 0;
    int qRemoved = 0;
    
    // Iterate CIGAR backwards
    while (!aln.CIGAR.empty() && (rRemoved < targetRefTrim || qRemoved < targetQryTrim)) {
        auto& op = aln.CIGAR.back();
        int len = op.first;
        char type = op.second;

        bool refC = consumesRef(type);
        bool qryC = consumesQry(type);

        int rStep = refC ? len : 0;
        int qStep = qryC ? len : 0;

        if ((rRemoved + rStep <= targetRefTrim) && (qRemoved + qStep <= targetQryTrim)) {
            // Remove whole op
            rRemoved += rStep;
            qRemoved += qStep;
            aln.CIGAR.pop_back();
        } else {
            // Partial trim
            int remainingRefNeeded = std::max(0, targetRefTrim - rRemoved);
            int remainingQryNeeded = std::max(0, targetQryTrim - qRemoved);

            int cutLen = 0;
            if (refC && qryC) cutLen = std::max(remainingRefNeeded, remainingQryNeeded);
            else if (refC) cutLen = remainingRefNeeded;
            else if (qryC) cutLen = remainingQryNeeded;

            op.first -= cutLen;
            rRemoved += (refC ? cutLen : 0);
            qRemoved += (qryC ? cutLen : 0);
            break;
        }
    }

    // Update Coordinates
    aln.refIdx.second -= rRemoved;

    if (!aln.inverse) {
        aln.qryIdx.second -= qRemoved;
    } else {
        // Reverse: Tail trim (CIGAR end) corresponds to qryIdx.first (low coord)
        aln.qryIdx.first += qRemoved;
    }
}


// Internal helper to perform a single split at a specific absolute coordinate
// Returns {Left Part, Right Part}
std::pair<mga::Alignment, mga::Alignment> mga::singleSplit(const Alignment& parent, int splitPos, bool onRef) {
    Alignment left = parent; 
    Alignment right = parent;
    
    // Bounds check to determine if split is even possible
    int pStart = onRef ? parent.refIdx.first : parent.qryIdx.first;
    int pEnd   = onRef ? parent.refIdx.second : parent.qryIdx.second;

    if (splitPos <= pStart) {
        left.valid = false; left.CIGAR.clear(); left.refIdx = {0,0}; left.qryIdx = {0,0};
        return {left, right}; // Split is at/before start -> Left is empty, Right is full Parent
    } 
    if (splitPos >= pEnd) {
        right.valid = false; right.CIGAR.clear(); right.refIdx = {0,0}; right.qryIdx = {0,0};
        return {left, right}; // Split is at/after end -> Left is full Parent, Right is empty
    }

    // Calculate Target Relative Offset
    // We convert the absolute splitPos into a "consumed bases needed" value
    int targetDist = 0;
    if (onRef) {
        targetDist = splitPos - parent.refIdx.first;
    } 
    else {
        // For Query, direction depends on strand
        if (!parent.inverse) targetDist = splitPos - parent.qryIdx.first;
        else targetDist = parent.qryIdx.second - splitPos; // Inverse: moves high -> low
    }

    // Iterate CIGAR to find cut point
    int currentDist = 0; // Tracks consumption of the *Split Axis*
    // Trackers for updating the *Other Axis*
    int otherAxisConsumed = 0;
    int cutOpIdx = -1;
    int cutOpOffset = -1; // How much of the OP goes to Left/Head

    for (int i = 0; i < parent.CIGAR.size(); ++i) {
        const auto& op = parent.CIGAR[i];
        int len = op.first;
        char type = op.second;
        
        int step = onRef ? (consumesRef(type) ? len : 0) : (consumesQry(type) ? len : 0);

        // Check if target falls inside this op
        if (currentDist + step > targetDist) {
            cutOpIdx = i;
            cutOpOffset = targetDist - currentDist; 
            
            // Add proportional consumption for the *Other Axis* for the partial part
            // (e.g. if we cut a Match, both axes advance. If we cut a Deletion, only Ref advances)
            int otherStepC = 0;
            if (onRef) otherStepC = consumesQry(type) ? cutOpOffset : 0; // If splitting Ref, how much Qry used?
            else       otherStepC = consumesRef(type) ? cutOpOffset : 0; // If splitting Qry, how much Ref used?
            
            otherAxisConsumed += otherStepC;
            currentDist += cutOpOffset;
            break;
        }

        currentDist += step;
        
        // Add full consumption for Other Axis
        if (onRef) otherAxisConsumed += consumesQry(type) ? len : 0;
        else       otherAxisConsumed += consumesRef(type) ? len : 0;
    }

    // Construct Left (Head) 
    left.CIGAR.clear();
    for(int i=0; i<cutOpIdx; ++i) left.CIGAR.push_back(parent.CIGAR[i]);
    if (cutOpOffset > 0) left.CIGAR.push_back({cutOpOffset, parent.CIGAR[cutOpIdx].second});

    // Update Indices Left
    // The "Split Axis" coordinate is simply set to splitPos
    // The "Other Axis" coordinate is calculated via otherAxisConsumed
    if (onRef) {
        left.refIdx.second = splitPos;
        if (!parent.inverse) left.qryIdx.second = parent.qryIdx.first + otherAxisConsumed;
        else                 left.qryIdx.first  = parent.qryIdx.second - otherAxisConsumed;
    } 
    else {
        if (!parent.inverse) left.qryIdx.second = splitPos;
        else                 left.qryIdx.first  = splitPos;
        left.refIdx.second = parent.refIdx.first + otherAxisConsumed;
    }
    left.identifier = parent.identifier * 10 + 1; // ID Logic

    // Construct Right (Tail)
    right.CIGAR.clear();
    int remLen = parent.CIGAR[cutOpIdx].first - cutOpOffset;
    if (remLen > 0) right.CIGAR.push_back({remLen, parent.CIGAR[cutOpIdx].second});
    for(int i=cutOpIdx+1; i<parent.CIGAR.size(); ++i) right.CIGAR.push_back(parent.CIGAR[i]);

    // Update Indices Right
    if (onRef) {
        right.refIdx.first = splitPos;
        if (!parent.inverse) right.qryIdx.first = parent.qryIdx.first + otherAxisConsumed;
        else                 right.qryIdx.second = parent.qryIdx.second - otherAxisConsumed;
    } 
    else {
        if (!parent.inverse) right.qryIdx.first = splitPos;
        else                 right.qryIdx.second = splitPos;
        right.refIdx.first = parent.refIdx.first + otherAxisConsumed;
    }
    right.identifier = parent.identifier * 10 + 2; // ID Logic
    return {left, right};
}

mga::SplitResult mga::splitAlignment(const Alignment& aln, int start, int end, bool onRef) {
    SplitResult result;

    // First Split: Cut at 'end' -> separates {Head+Middle} from {Tail}
    // We assign this to temporary variables
    auto split1 = singleSplit(aln, end, onRef);
    Alignment tempHead = split1.first;
    result.tail = split1.second; // This is the final Tail

    // Second Split: Cut 'tempHead' at 'start' -> separates {Head} from {Middle}
    if (!tempHead.valid) {
        // The cut was way before the alignment started or invalid range
        result.head.valid = false;
        result.middle.valid = false;
        return result; 
    }

    auto split2 = singleSplit(tempHead, start, onRef);
    result.head = split2.first;   // Final Head
    result.middle = split2.second; // Final Middle

    // ID Cleanup (Optional)
    // The helpers generate recursive IDs (e.g., 101, 102). You might want to reset them here.
    if (result.head.valid)   result.head.identifier   = aln.identifier * 10 + 1;
    if (result.middle.valid) result.middle.identifier = aln.identifier * 10 + 2;
    if (result.tail.valid)   result.tail.identifier   = aln.identifier * 10 + 3;
    return result;
}



void mga::detectDuplications(alnVec& alignments) {
    const int FUZZY_THRESHOLD = 100; 
    const int MIN_OVERLAP = 200; 

    // Separate Indices
    std::list<int> primaryIndices;
    std::vector<int> secondaryIndices;
    for (int i = 0; i < alignments.size(); ++i) {
        if (!alignments[i].valid) continue;
        if (alignments[i].primary) primaryIndices.push_back(i);
        else secondaryIndices.push_back(i);
    }

    
    primaryIndices.sort( [&] (int a, int b) {
        return alignments[a].refIdx.first < alignments[b].refIdx.first;
    });
    // std::ofstream outFile ("chain.txt");
    // outFile << 1 << ':' << 0 << "\n";
    for (auto& p: primaryIndices) {
        // outFile << '(' << alignments[p].refIdx.first << ',' << alignments[p].refIdx.second << "),(" 
        //                << alignments[p].qryIdx.first << ',' << alignments[p].qryIdx.second << ")\n";
        printCoordinate(alignments[p].refIdx.first, alignments[p].refIdx.second);
        std::cerr << '\t';
        printCoordinate(alignments[p].qryIdx.first, alignments[p].qryIdx.second);
        std::cerr << '\n';
    }
    // outFile.close();

    // Sort Secondaries (High -> Low Score)
    std::sort(secondaryIndices.begin(), secondaryIndices.end(), [&](int a, int b) {
        return alignments[a].alnScore > alignments[b].alnScore;
    });

    std::list<int> processQueue(secondaryIndices.begin(), secondaryIndices.end());

    int sCount = 0;
    while (!processQueue.empty()) {
        int sIdx = processQueue.front();
        processQueue.pop_front();
       
        
        Alignment* currS = &alignments[sIdx];
        if (!currS->valid) continue;
        sCount++;
        std::cerr << "[" << sCount << "]: " << currS->identifier << '\t'; printCoordinate(currS->refIdx.first, currS->refIdx.second);
        std::cerr << '\t'; printCoordinate(currS->qryIdx.first, currS->qryIdx.second);
        std::cerr << '\n';
        // --- Step 1: Find overlaps to Primaries ---
        std::vector<OverlapInfo> overlaps;
        for (int pIdx : primaryIndices) {
            Alignment& P = alignments[pIdx];
            if (!P.primary || !P.valid) continue;

            // Ref
            int rStart = std::max(currS->refIdx.first, P.refIdx.first);
            int rEnd = std::min(currS->refIdx.second, P.refIdx.second);
            if (rEnd > rStart) overlaps.push_back({pIdx, true, rEnd - rStart, rStart, rEnd});
            // Qry
            int qStart = std::max(currS->qryIdx.first, P.qryIdx.first);
            int qEnd = std::min(currS->qryIdx.second, P.qryIdx.second);
            if (qEnd > qStart) overlaps.push_back({pIdx, false, qEnd - qStart, qStart, qEnd});
        }

        if (overlaps.empty()) {
            // No overlaps, S is an orphan secondary
            alignments[sIdx].primary = true;
            primaryIndices.push_back(sIdx);
            std::cerr << "Added to Primary\n";
            continue; 
        }

        // --- Step 2: Select Best Overlap (Longest) ---
        auto bestOv = *std::max_element(overlaps.begin(), overlaps.end(), 
            [](const OverlapInfo& a, const OverlapInfo& b){ return a.overlapLen < b.overlapLen; });
        
        std::cerr << "Total overlaps: " << overlaps.size() << '\n';
        for (auto op: overlaps) {
            std::cout << op.overlapLen << '\t' << op.primaryIdx << '\t' << op.onRef << '\t'; printCoordinate(op.start, op.end);
            std::cerr << '\n';
        }
        if (bestOv.overlapLen < MIN_OVERLAP) continue;
        std::cout << "BEST on " << bestOv.primaryIdx << ": " << bestOv.overlapLen << ',' << bestOv.onRef << '\t'; printCoordinate(bestOv.start, bestOv.end);
        std::cerr << '\n';
        

        // --- Step 3: Split Secondary (S) to create the "Bridge" ---
        SplitResult resS = splitAlignment(*currS, bestOv.start, bestOv.end, bestOv.onRef);
        
        if (!resS.middle.valid) continue; 
        currS->valid = false; // Invalidate original S

        // Push Head/Tail back to queue
        if (resS.head.valid) {
            alignments.push_back(resS.head);
            processQueue.push_front(alignments.size() - 1);
            std::cerr << "Add secondary head.\n";
        }
        if (resS.tail.valid) {
            alignments.push_back(resS.tail);
            processQueue.push_front(alignments.size() - 1);
            std::cerr << "Add secondary tail.\n";
        }

        // Add Middle (Bridge)
        alignments.push_back(resS.middle);
        int sMidIdx = alignments.size() - 1;
        Alignment* sMid = &alignments[sMidIdx];

        // --- Step 4: Find Cross-Overlaps on the OTHER axis ---
        // If bestOv was on Ref, we look at Query. If BestOv on Query, look at Ref.
        bool scanRef = !bestOv.onRef;
        std::vector<CrossOverlap> crossOverlaps;

        for (int pIdx : primaryIndices) {
            Alignment& P = alignments[pIdx];
            if (!P.primary || !P.valid) continue;

            // Calculate overlap range on the 'Other' Axis
            int ovStart, ovEnd, len;
            if (scanRef) {
                ovStart = std::max(sMid->refIdx.first, P.refIdx.first);
                ovEnd = std::min(sMid->refIdx.second, P.refIdx.second);
            } else {
                ovStart = std::max(sMid->qryIdx.first, P.qryIdx.first);
                ovEnd = std::min(sMid->qryIdx.second, P.qryIdx.second);
            }
            len = ovEnd - ovStart;

            if (len > MIN_OVERLAP) {
                // IMPORTANT: We found an overlap on the 'Other' axis.
                // We use mapCoordinate to verify collision on P_Main's axis 
                // and to determine where to cut P_Main later.

                // mapCoordinate takes (pos, inputIsRef). 
                // If scanRef is true, we are scanning Reference coords, so inputIsRef=true.
                int mapS = sMid->mapCoordinate(ovStart, scanRef);
                int mapE = sMid->mapCoordinate(ovEnd, scanRef);

                // Check mapping validity (might return -1 if inside a deletion/insertion)
                if (mapS == -1 || mapE == -1) continue;

                // Handle Inversion: mapS might be > mapE if sMid is inverse
                int pS = std::min(mapS, mapE);
                int pE = std::max(mapS, mapE);

                // --- Constraint Check: Overlapping regions on same primary ---
                // If this P is the SAME as the BestOverlap P, we must ensure the 
                // 'Other' overlap doesn't collide with the 'Best' overlap region on the Main Axis.
                if (pIdx == bestOv.primaryIdx) {
                    // Check intersection: [pS, pE] vs [bestOv.start, bestOv.end]
                    // Standard intersection logic: max(Starts) < min(Ends)
                    if (std::max(pS, bestOv.start) < std::min(pE, bestOv.end)) {
                        continue; // They collide, skip this cross-overlap
                    }
                }

                // If valid, store it
                crossOverlaps.push_back({pIdx, ovStart, ovEnd, pS, pE, !scanRef});
            }
        }

        for (auto ov: crossOverlaps) {
            std::cout << "Cross on " << ov.pIdx << ": " << ov.onRef << '\t'; 
            printCoordinate(ov.sStart, ov.sEnd);
            printCoordinate(ov.pStart, ov.pEnd);
            std::cerr << '\n';
        }

        // --- Step 5: Process Cross Overlaps ---
        
        if (crossOverlaps.empty()) {
            std::cerr << "No crossoverlap. Mark duplications.\n";
            // No cross overlaps -> Simple Duplication (P_Main <-> S_Mid)
            
            int currP_Idx = bestOv.primaryIdx;
            Alignment* currP = &alignments[currP_Idx];
            
            // Fuzzy Logic
            int cutS = bestOv.start;
            int cutE = bestOv.end;
            int pStart = bestOv.onRef ? currP->refIdx.first : currP->qryIdx.first;
            int pEnd   = bestOv.onRef ? currP->refIdx.second : currP->qryIdx.second;
            int offsetStart = cutS - pStart;
            int offsetEnd = pEnd - cutE;
            if (std::abs(offsetStart) < FUZZY_THRESHOLD) cutS = pStart;
            if (std::abs(offsetEnd) < FUZZY_THRESHOLD)   cutE = pEnd;

            // Split P
            if (cutS > pStart || cutE < pEnd) {
                 SplitResult resP = splitAlignment(*currP, cutS, cutE, bestOv.onRef);
                 if (resP.middle.valid) {
                     currP->valid = false;
                     // Add Head/Tail
                     if (resP.head.valid) { resP.head.primary = true; alignments.push_back(resP.head); primaryIndices.push_back(alignments.size()-1); std::cerr << "Add primary head.\n";}
                     if (resP.tail.valid) { resP.tail.primary = true; alignments.push_back(resP.tail); primaryIndices.push_back(alignments.size()-1); std::cerr << "Add primary tail.\n";}
                     // Add Middle
                     resP.middle.primary = true; alignments.push_back(resP.middle);
                     currP_Idx = alignments.size() - 1;
                 }
            }
            // Mark
            alignments[currP_Idx].duplications.insert(sMid->identifier);
            // sMid->valid = false; // Consumed
        } 
        else {
            // Complex Case: S_Mid bridges P_Main to multiple P_Cross
            // 1. Sort Cross Overlaps by P_Main's Axis Coordinates (Low -> High)
            std::sort(crossOverlaps.begin(), crossOverlaps.end(), 
                [](const CrossOverlap& a, const CrossOverlap& b){ return a.pStart < b.pStart; });

            // 2. Iterate and Cut P_Main
            int currP_Idx = bestOv.primaryIdx; 

            for (const auto& co : crossOverlaps) {
                Alignment* currP = &alignments[currP_Idx];
                if (!currP->valid) break; 

                // Cut Coords (Mapped from Cross Overlap using mapCoordinate earlier)
                int cutS = co.pStart;
                int cutE = co.pEnd;
                
                // Fuzzy Logic
                int pStart = bestOv.onRef ? currP->refIdx.first : currP->qryIdx.first;
                int pEnd   = bestOv.onRef ? currP->refIdx.second : currP->qryIdx.second;
                if (std::abs(cutS - pStart) < FUZZY_THRESHOLD) cutS = pStart;
                if (std::abs(pEnd - cutE) < FUZZY_THRESHOLD)   cutE = pEnd;

                // Split P_Main to match CrossOverlap
                int matchP_Idx = currP_Idx; 
                
                if (cutS > pStart || cutE < pEnd) {
                    SplitResult resP = splitAlignment(*currP, cutS, cutE, bestOv.onRef);
                    if (resP.middle.valid) {
                        currP->valid = false;
                        
                        // Head (Keep as primary)
                        if (resP.head.valid) { 
                            resP.head.primary = true; 
                            alignments.push_back(resP.head); 
                            primaryIndices.push_back(alignments.size()-1); 
                            std::cerr << "Cross (Add head): "; printCoordinate(resP.head.refIdx.first, resP.head.refIdx.second);
                            std::cerr << '\n';
                        }
                        
                        // Middle (The Match)
                        resP.middle.primary = true;
                        alignments.push_back(resP.middle);
                        matchP_Idx = alignments.size() - 1;
                        std::cerr << "Cross (Add middle): "; printCoordinate(resP.middle.refIdx.first, resP.middle.refIdx.second);
                        std::cerr << '\n';

                        // Tail (Remainder for next loop)
                        if (resP.tail.valid) {
                            resP.tail.primary = true;
                            alignments.push_back(resP.tail);
                            currP_Idx = alignments.size() - 1; 
                            primaryIndices.push_back(currP_Idx);
                            std::cerr << "Cross (Add head): "; printCoordinate(resP.tail.refIdx.first, resP.tail.refIdx.second);
                            std::cerr << '\n';
                        } else {
                            currP_Idx = -1; 
                        }
                    }
                }

                // 3. Split the OTHER Primary (P_Cross) to match S
                // We split it at sStart/sEnd (the coords on the Other Axis)
                int crossP_Idx = co.pIdx;
                Alignment* crossP = &alignments[crossP_Idx];
                
                // Fuzzy for CrossP
                int cStart = scanRef ? crossP->refIdx.first : crossP->qryIdx.first;
                int cEnd   = scanRef ? crossP->refIdx.second : crossP->qryIdx.second;
                int cCutS = co.sStart;
                int cCutE = co.sEnd;
                if (std::abs(cCutS - cStart) < FUZZY_THRESHOLD) cCutS = cStart;
                if (std::abs(cEnd - cCutE) < FUZZY_THRESHOLD)   cCutE = cEnd;

                if (cCutS > cStart || cCutE < cEnd) {
                     SplitResult resC = splitAlignment(*crossP, cCutS, cCutE, scanRef);
                     if (resC.middle.valid) {
                         crossP->valid = false;
                         if (resC.head.valid) { resC.head.primary=true; alignments.push_back(resC.head); primaryIndices.push_back(alignments.size()-1); }
                         if (resC.tail.valid) { resC.tail.primary=true; alignments.push_back(resC.tail); primaryIndices.push_back(alignments.size()-1); }
                         resC.middle.primary = true; alignments.push_back(resC.middle);
                         crossP_Idx = alignments.size() - 1;
                     }
                }

                // 4. Mark Duplication
                alignments[matchP_Idx].duplications.insert(alignments[crossP_Idx].identifier);
                alignments[crossP_Idx].duplications.insert(alignments[matchP_Idx].identifier);
            }
            
            // Consumed
             sMid->valid = false; 
        }
    }
}
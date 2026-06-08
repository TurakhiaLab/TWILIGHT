
#include "alignment.hpp"
#include "mga.hpp"

#include <fstream>
#include <sstream>
#include <chrono>

// Helper to check if a CIGAR op consumes Reference
bool consumesRef(char op) {
    return (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X');
}

// Helper to check if a CIGAR op consumes Query
bool consumesQry(char op) {
    return (op == 'M' || op == 'I' || op == 'S' || op == '=' || op == 'X');
}

int mapCoordinate(const CigarString& cigar, int q_start, int r_start, int target, bool is_target_ref) {
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

std::string getReverseComplement(std::string seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break; case 'T': c = 'A'; break;
            case 'C': c = 'G'; break; case 'G': c = 'C'; break;
            case 'a': c = 't'; break; case 't': c = 'a'; break;
            case 'c': c = 'g'; break; case 'g': c = 'c'; break;
        }
    }
    return seq;
}

CigarString compressCigar(const CigarString& cigar) {
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

void printCIGAR(const CigarString& cigar, bool changeLine) {
    std::string cigar_str;
    for (const auto& op : cigar) {
        cigar_str += std::to_string(op.first) + op.second;
    }
    if (changeLine) std::cout << cigar_str << std::endl;
    else std::cout << cigar_str;
}

// Helper to parse CIGAR string into vector of ops
CigarString parser::parseCigar(const std::string& cigarString) {
    CigarString ops;
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

Alignments parser::parseMinimap2PAF(const std::string& filename) {
    Alignments alignments;
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
        aln.ID = alignmentCounter++;

        std::string cigarString;
        int numMismatch = 0, numIndel = 0, numEdits = 0, score = 0;
        int pafRefStart = 0, pafRefEnd = 0, pafQryStart = 0, pafQryEnd = 0;

        aln.qryName = fields[0]; // Query sequence name
        aln.refName = fields[5]; // Reference sequence name
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

        aln.alnLength = totLen;
        aln.alnScore = score;
        aln.ins = ins;
        aln.del = del;
        aln.mis = (numEdits - ins - del);

        alignments.push_back(aln);
    }

    // Sort alignments: Group by RefStart, then by QryStart
    std::sort(alignments.begin(), alignments.end(), [&](const auto &a, const auto &b) {
        if (a.refIdx.first == b.refIdx.first) return a.qryIdx.first < b.qryIdx.first;
        return (a.refIdx.first < b.refIdx.first); 
    });

    return alignments;
}

CigarString adjustCigarWithVariations(CigarString& origCigar, Segment& refSeg, Segment& qrySeg, bool qryInverse, int qryConsLen) {
    CigarString newCigar;
    
    std::map<int, int> refGaps;
    for (auto& v : refSeg.getVariants()) {
        if (v.getType() == VariantType::GAP) {
            refGaps[v.getStart()] = v.getEnd() - v.getStart();
        }
    }

    std::map<int, int> qryGaps;
    for (auto& v : qrySeg.getVariants()) {
        if (v.getType() == VariantType::GAP) {
            int s = v.getStart();
            int e = v.getEnd();
            if (qryInverse) {
                int old_s = s;
                s = qryConsLen - e;
                e = qryConsLen - old_s;
            }
            qryGaps[s] = e - s;
        }
    }

    int rConsPos = 0; 
    int qConsPos = 0; 

    auto appendOp = [&](int len, char op) {
        if (len > 0) {
            if (!newCigar.empty() && newCigar.back().second == op) {
                newCigar.back().first += len;
            } else {
                newCigar.push_back({len, op});
            }
        }
    };

    for (const auto& op : origCigar) {
        int remain = op.first;
        char type = op.second;

        while (remain > 0) {
            int nextRefGapDist = remain + 1;
            int nextQryGapDist = remain + 1;
            
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            if (consumesRef) {
                auto it = refGaps.lower_bound(rConsPos);
                if (it != refGaps.end() && it->first < rConsPos + remain) {
                    nextRefGapDist = it->first - rConsPos;
                }
            }
            if (consumesQry) {
                auto it = qryGaps.lower_bound(qConsPos);
                if (it != qryGaps.end() && it->first < qConsPos + remain) {
                    nextQryGapDist = it->first - qConsPos;
                }
            }

            int step = std::min({remain, nextRefGapDist, nextQryGapDist});

            if (step > 0) {
                appendOp(step, type);
                if (consumesRef) rConsPos += step;
                if (consumesQry) qConsPos += step;
                remain -= step;
            }

            if (consumesRef && refGaps.count(rConsPos)) {
                int gLen = refGaps[rConsPos];
                appendOp(gLen, 'D'); 
                rConsPos += gLen;
            }
            if (consumesQry && qryGaps.count(qConsPos)) {
                int gLen = qryGaps[qConsPos];
                appendOp(gLen, 'I'); 
                qConsPos += gLen;
            }
        }
    }

    while (refGaps.count(rConsPos)) {
        int gLen = refGaps[rConsPos];
        appendOp(gLen, 'D');
        rConsPos += gLen;
    }
    while (qryGaps.count(qConsPos)) {
        int gLen = qryGaps[qConsPos];
        appendOp(gLen, 'I');
        qConsPos += gLen;
    }

    return newCigar;
}

CigarString extractSubCigar(const CigarString& origCigar, int refOffset, int refLen) {
    CigarString subCigar;
    int currentRef = 0;
    int extractedRef = 0;

    auto appendOp = [&](int len, char type) {
        if (len <= 0) return;
        if (!subCigar.empty() && subCigar.back().second == type) {
            subCigar.back().first += len;
        } else {
            subCigar.push_back({len, type});
        }
    };

    for (const auto& op : origCigar) {
        if (extractedRef >= refLen) break; 

        int opLen = op.first;
        char type = op.second;

        bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
        int refOpLen = consumesRef ? opLen : 0;
        
        if (consumesRef) {
            if (currentRef + refOpLen <= refOffset) {
                currentRef += refOpLen;
                continue;
            }
        } else {
            if (currentRef < refOffset) {
                continue;
            }
        }

        int useLen = opLen;

        if (consumesRef && currentRef < refOffset) {
            int trim = refOffset - currentRef;
            useLen -= trim;
            currentRef += trim;
        }

        if (consumesRef && (extractedRef + useLen > refLen)) {
            useLen = refLen - extractedRef;
        }

        appendOp(useLen, type);

        if (consumesRef) {
            currentRef += useLen;
            extractedRef += useLen;
        }
    }
    
    return subCigar;
}

Alignments runMinimap2(StringPairs& ref, StringPairs& qry, std::string refName, std::string qryName, Option& option) {
    // 1. Prepare temporary files for minimap2
    std::string temp_dir = option.tempDir;
    std::string refFile_path = temp_dir + "/" + refName + ".fa";
    std::string qryFile_path = temp_dir + "/" + qryName + ".fa";
        
    mga::io::writeAlignment(refFile_path, ref, false, false);
    mga::io::writeAlignment(qryFile_path, qry, false, false);
        
    // whole-genome alignment between 1 and 2
    std::string PAF_path = temp_dir + "/output_" + refName + "_" + qryName + ".paf";

    // 2. Run minimap2
    const char* home_dir = getenv("HOME");
    if (home_dir == nullptr) {
        std::cerr << "Error: Could not get HOME directory." << std::endl;
        exit(1);
    }
        
    // std::string minimap2_path = std::string(home_dir) + "/bin/minimap2";
    std::string minimap2_path = "/home/y3tseng@AD.UCSD.EDU/minimap2/minimap2";
    std::string command;
    int system_ret; 

    // minimap2 -cx asm5 -g 500 -r 500 -n 5 -m 50 -N 20 -p 0.8 asm1.fa asm2.fa > aln.paf
    // whole-genome alignment between 1 and 2
    command = minimap2_path + " -cx asm5 " + refFile_path + " " + qryFile_path + " > " + PAF_path + " 2> /dev/null";
        
    auto minimap2_start = std::chrono::high_resolution_clock::now();
    system_ret = system(command.c_str());
    auto minimap2_end = std::chrono::high_resolution_clock::now();
    option.minimap2_time += std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count();
        
    if (system_ret != 0) {
        std::cerr << "Error: minimap2 execution failed for command: " << command << std::endl;
        exit(1);
    }
        
    // Clear sequence file
    // std::remove(refFile_path.c_str());
    // std::remove(qryFile_path.c_str());
    // 3. Parse PAF output
    auto minimap2_alignments = parser::parseMinimap2PAF(PAF_path);
    return minimap2_alignments;
}

Alignments splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts,const std::set<int>& qryCuts)  {
    Alignments frags; // 用來裝切碎的子片段

    // 1. 篩選並排序切點 (WGA Edge-Cut Fix)
    std::vector<int> rCuts;
    for (int c : refCuts) {
        // Ref 永遠是正向掃描，所以允許在終點 (second) 切割
        if (c > aln.refIdx.first && c <= aln.refIdx.second) rCuts.push_back(c);
    }

    std::vector<int> qCuts;
    for (int c : qryCuts) {
        if (aln.inverse) {
            // 反向掃描 (從 second 往 down 走到 first)，允許在終點 (first) 切割
            if (c >= aln.qryIdx.first && c < aln.qryIdx.second) qCuts.push_back(c);
        } else {
            // 正向掃描，允許在終點 (second) 切割
            if (c > aln.qryIdx.first && c <= aln.qryIdx.second) qCuts.push_back(c);
        }
    }

    if (aln.inverse) {
        std::sort(qCuts.rbegin(), qCuts.rend());
    } else {
        std::sort(qCuts.begin(), qCuts.end());
    }

    if (rCuts.empty() && qCuts.empty()) {
        frags.push_back(aln);
        return frags;
    }

    // 2. 準備走訪 CIGAR 進行動態切割
    int rCutIdx = 0;
    int qCutIdx = 0;

    int rPos = aln.refIdx.first;
    int qPos = aln.inverse ? aln.qryIdx.second : aln.qryIdx.first; 
    int qDir = aln.inverse ? -1 : 1;

    int currRStart = rPos;
    int currQStart = qPos;

    Alignment currAln = aln;
    currAln.CIGAR.clear(); 

    // 3. 逐一消耗 CIGAR Operations
    for (const auto& op : aln.CIGAR) {
        int len = op.first;
        char type = op.second;

        while (len > 0) {
            bool consumesRef = (type == 'M' || type == '=' || type == 'X' || type == 'D');
            bool consumesQry = (type == 'M' || type == '=' || type == 'X' || type == 'I');

            int step = len;

            if (consumesRef && rCutIdx < rCuts.size()) {
                int distR = rCuts[rCutIdx] - rPos;
                if (distR > 0 && distR < step) step = distR;
            }

            if (consumesQry && qCutIdx < qCuts.size()) {
                int distQ = std::abs(qCuts[qCutIdx] - qPos);
                if (distQ > 0 && distQ < step) step = distQ;
            }

            if (!currAln.CIGAR.empty() && currAln.CIGAR.back().second == type) {
                currAln.CIGAR.back().first += step; 
            } else {
                currAln.CIGAR.push_back({step, type});
            }

            if (consumesRef) rPos += step;
            if (consumesQry) qPos += step * qDir;
            len -= step;

            // 4. 檢查是否精準踩到切點
            bool hitRef = (consumesRef && rCutIdx < rCuts.size() && rPos == rCuts[rCutIdx]);
            bool hitQry = (consumesQry && qCutIdx < qCuts.size() && qPos == qCuts[qCutIdx]);

            if (hitRef || hitQry) {
                currAln.refIdx.first = currRStart;
                currAln.refIdx.second = rPos;

                if (aln.inverse) {
                    currAln.qryIdx.first = qPos;        
                    currAln.qryIdx.second = currQStart; 
                } else {
                    currAln.qryIdx.first = currQStart;  
                    currAln.qryIdx.second = qPos;       
                }

                if (!currAln.CIGAR.empty()) {
                    frags.push_back(currAln); // 改存進 frags
                }

                currRStart = rPos;
                currQStart = qPos;
                currAln = aln; 
                currAln.CIGAR.clear();

                if (hitRef) rCutIdx++;
                if (hitQry) qCutIdx++;
            }
        }
    }

    // 5. 收尾
    if (!currAln.CIGAR.empty()) {
        currAln.refIdx.first = currRStart;
        currAln.refIdx.second = rPos;

        if (aln.inverse) {
            currAln.qryIdx.first = qPos;
            currAln.qryIdx.second = currQStart;
        } else {
            currAln.qryIdx.first = currQStart;
            currAln.qryIdx.second = qPos;
        }
        frags.push_back(currAln); // 改存進 frags
    }

    // 6. Update alignment length
    for (auto& frag : frags) {
        frag.updateAlnLength();
    }

    return frags;
}

void snapAlignment(Alignment& aln, int r_pad_left, int r_pad_right, int q_pad_left, int q_pad_right) {
    std::vector<std::pair<int, char>> prepend_ops;
    std::vector<std::pair<int, char>> append_ops;

    // 1. Ref 的補丁 (Ref 永遠是正向)
    if (r_pad_left > 0)  prepend_ops.push_back({r_pad_left, 'D'});
    if (r_pad_right > 0) append_ops.push_back({r_pad_right, 'D'});

    // 2. Qry 的補丁 (根據 Inverse 決定加在頭還是尾)
    if (aln.inverse) {
        if (q_pad_right > 0) prepend_ops.push_back({q_pad_right, 'I'});
        if (q_pad_left > 0)  append_ops.push_back({q_pad_left, 'I'});
    } else {
        if (q_pad_left > 0)  prepend_ops.push_back({q_pad_left, 'I'});
        if (q_pad_right > 0) append_ops.push_back({q_pad_right, 'I'});
    }

    // 3. 將 D/I 補丁貼上 CIGAR
    if (!prepend_ops.empty()) aln.CIGAR.insert(aln.CIGAR.begin(), prepend_ops.begin(), prepend_ops.end());
    if (!append_ops.empty())  aln.CIGAR.insert(aln.CIGAR.end(), append_ops.begin(), append_ops.end());

    // 4. 更新座標 (因為你保證了 first 永遠小於 second，所以直接加減即可)
    aln.refIdx.first  -= r_pad_left;
    aln.refIdx.second += r_pad_right;

    aln.qryIdx.first  -= q_pad_left;
    aln.qryIdx.second += q_pad_right;

    // 5. 更新對齊總長度
    aln.updateAlnLength(); 
}
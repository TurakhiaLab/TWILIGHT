
#include "alignment.hpp"
#include "cigar_util.hpp"

extern "C" {
#include "minimap.h"
}

#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdio>



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
        
        // Set a default value in case the tp tag is missing from the PAF
        aln.primary = true; 

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
    
        // Parse Tags (starting at Field 12) for cg:Z (CIGAR), NM:i, AS:i, and tp:A
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
            else if (fields[k].rfind("tp:A:", 0) == 0) {
                // If the tag is exactly 'P', it's primary. Otherwise (e.g., 'S' or 'I'), it's not.
                aln.primary = (fields[k].substr(5) == "P");
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

Alignments runMinimap2(const std::string& refSeq, const std::string& qrySeq,
                       const std::string& refID, const std::string& qryID,
                       Option& option, bool needCigar, bool write_paf) {
    auto minimap2_start = std::chrono::high_resolution_clock::now();
    Alignments alignments;
    if (refSeq.empty() || qrySeq.empty()) return alignments;

    // 1. 初始化 Minimap2 選項 (必須先套用 0 預設值，再套用 preset，否則其他欄位將為隨機垃圾值)
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_set_opt(0, &iopt, &mopt);
    mm_set_opt("asm5", &iopt, &mopt);

    mopt.pri_ratio = 0.0f;     // -p 0
    mopt.best_n = 50;          // -N 50
    // mopt.max_gap = 2000;    // 相當於 -g 1000 (原預設為 10000)
    // mopt.bw = 2000;
    if (needCigar) {
        mopt.flag |= MM_F_CIGAR;   // 啟用 CIGAR (-c)
    }

    // 2. 零拷貝 C 字串指標傳遞 (直接引用傳入的 const std::string&)
    const char* ref_seq_ptr = refSeq.c_str();
    const char* ref_name_ptr = refID.c_str();

    // 3. 在記憶體中直接建立 Index (零檔案 I/O，零字串拷貝)
    mm_idx_t* mi = mm_idx_str(iopt.w, iopt.k, iopt.flag & MM_I_HPC, 14, 1, &ref_seq_ptr, &ref_name_ptr);
    if (!mi) {
        std::cerr << "Error: Failed to build minimap2 index in memory." << std::endl;
        return alignments;
    }
    mm_mapopt_update(&mopt, mi);

    // 4. 初始化 Thread Buffer
    mm_tbuf_t* tbuf = mm_tbuf_init();

    // 5. 執行 Query 序列 Mapping
    int n_regs = 0;
    mm_reg1_t* regs = mm_map(mi, static_cast<int>(qrySeq.length()), qrySeq.c_str(),
                             &n_regs, tbuf, &mopt, qryID.c_str());

    std::ofstream pafOut;
    if (write_paf) {
        std::string pafPath = option.tempDir + "/" + (refID.empty() ? "ref" : refID) + "_" + (qryID.empty() ? "qry" : qryID) + ".paf";
        pafOut.open(pafPath, std::ios::out);
        if (pafOut.is_open()) {
            std::cout << "[DEBUG runMinimap2] Writing PAF output to: " << pafPath << std::endl;
        }
    }

    for (int j = 0; j < n_regs; ++j) {
        mm_reg1_t* r = &regs[j];

        if (write_paf && pafOut.is_open()) {
            pafOut << qryID << "\t"
                   << qrySeq.length() << "\t"
                   << r->qs << "\t"
                   << r->qe << "\t"
                   << (r->rev ? "-" : "+") << "\t"
                   << refID << "\t"
                   << refSeq.length() << "\t"
                   << r->rs << "\t"
                   << r->re << "\t"
                   << r->mlen << "\t"
                   << r->blen << "\t"
                   << r->mapq << "\t"
                   << "tp:A:" << (r->sam_pri ? "P" : "S") << "\n";
        }

        Alignment aln;
        aln.ID = 0;
        aln.refName = refID;
        aln.qryName = qryID;
        aln.refIdx = {r->rs, r->re};
        aln.qryIdx = {r->qs, r->qe};
        aln.inverse = (r->rev != 0);
        aln.primary = (r->sam_pri != 0);
        aln.valid = true;

        int totLen = 0, ins = 0, del = 0, numEdits = 0;
        int score = 0;

        if (r->p) {
            score = r->p->dp_max;
            numEdits = r->p->n_ambi;

            aln.CIGAR.reserve(r->p->n_cigar);
            for (uint32_t k = 0; k < r->p->n_cigar; ++k) {
                uint32_t c32 = r->p->cigar[k];
                int length = c32 >> 4;
                char opType = "MIDNSHPE="[c32 & 0xf];

                aln.CIGAR.push_back({length, opType});

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
            free(r->p); // 釋放 CIGAR payload 記憶體
        } else {
            // Chaining-only mode: 使用 block length 作為對齊長度
            totLen = r->blen;
        }

        aln.alnLength = totLen;
        aln.alnScore = score;
        aln.ins = ins;
        aln.del = del;
        aln.mis = std::max(0, numEdits - ins - del);

        alignments.push_back(std::move(aln));
    }
    free(regs);

    // 6. 清理 Minimap2 記憶體
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);

    // Sort alignments: Group by RefStart, then by QryStart
    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        if (a.refIdx.first == b.refIdx.first) return a.qryIdx.first < b.qryIdx.first;
        return (a.refIdx.first < b.refIdx.first);
    });

    auto minimap2_end = std::chrono::high_resolution_clock::now();
    option.minimap2_time += std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count();

    return alignments;
}

Alignments runMinimap2(const SequenceRefs& ref, const SequenceRefs& qry, std::string refName, std::string qryName, Option& option, bool needCigar, bool write_paf) {
    auto minimap2_start = std::chrono::high_resolution_clock::now();
    Alignments alignments;
    if (ref.empty() || qry.empty()) {
        std::cout << "[DEBUG runMinimap2] Reference or Query is empty. Ref size: " << ref.size() << ", Qry size: " << qry.size() << std::endl;
        return alignments;
    }

    std::cout << "[DEBUG runMinimap2] Start running minimap2 API. Ref count: " << ref.size() << ", Qry count: " << qry.size() << std::endl;

    // 1. 初始化 Minimap2 選項 (必須先套用 0 預設值，再套用 preset，否則其他欄位將為隨機垃圾值)
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    mm_set_opt(0, &iopt, &mopt);
    mm_set_opt("asm5", &iopt, &mopt);

    mopt.pri_ratio = 0.0f;     // -p 0
    mopt.best_n = 50;          // -N 50
    if (needCigar) {
        mopt.flag |= MM_F_CIGAR;   // 啟用 CIGAR (-c)
    }

    // 2. 將 Reference 序列寫入暫存檔以利正確建立 Index (解決 mm_idx_str 處理大序列時的 Bug)
    std::string tempRefFasta = option.tempDir + "/temp_minimap2_ref.fa";
    std::ofstream refOut(tempRefFasta);
    if (!refOut.is_open()) {
        std::cerr << "[DEBUG runMinimap2] Error: Failed to open temp ref file: " << tempRefFasta << std::endl;
        return alignments;
    }
    for (size_t i = 0; i < ref.size(); ++i) {
        refOut << ">" << ref[i].name << "\n" << ref[i].seq << "\n";
    }
    refOut.close();

    // 3. 用 mm_idx_reader_open 建立 Index
    std::cout << "[DEBUG runMinimap2] Building minimap2 index via reader for reference genome..." << std::endl;
    mm_idx_reader_t* r = mm_idx_reader_open(tempRefFasta.c_str(), &iopt, 0);
    if (!r) {
        std::cerr << "[DEBUG runMinimap2] Error: Failed to open minimap2 index reader." << std::endl;
        return alignments;
    }
    mm_idx_t* mi = mm_idx_reader_read(r, 1);
    mm_idx_reader_close(r);

    if (!mi) {
        std::cerr << "[DEBUG runMinimap2] Error: Failed to build minimap2 index." << std::endl;
        return alignments;
    }
    std::cout << "[DEBUG runMinimap2] Index built successfully. Updating map options..." << std::endl;
    mm_mapopt_update(&mopt, mi);

    // 4. 初始化 Thread Buffer
    std::cout << "[DEBUG runMinimap2] Initializing thread buffer..." << std::endl;
    mm_tbuf_t* tbuf = mm_tbuf_init();

    std::ofstream pafOut;
    if (write_paf) {
        std::string pafPath = option.tempDir + "/" + (refName.empty() ? "ref" : refName) + "_" + (qryName.empty() ? "qry" : qryName) + ".paf";
        pafOut.open(pafPath, std::ios::out);
        if (pafOut.is_open()) {
            std::cout << "[DEBUG runMinimap2] Writing PAF output to: " << pafPath << std::endl;
        }
    }

    // 5. 零拷貝過渡 Query 序列進行 mapping
    std::cout << "[DEBUG runMinimap2] Starting query sequence mapping loop..." << std::endl;
    for (const auto& qRef : qry) {
        const std::string& q_name = qRef.name;
        const std::string& q_seq = qRef.seq;
        int n_regs = 0;

        std::cout << "[DEBUG runMinimap2] Mapping query: " << q_name << " (length: " << q_seq.length() << " bp)..." << std::endl;
        
        // 診斷是否有非法字元 (例如 Gap '-', Space, Newline 等)
        size_t invalid_chars = 0;
        for (char c : q_seq) {
            char uc = std::toupper(static_cast<unsigned char>(c));
            if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T' && uc != 'N') {
                invalid_chars++;
            }
        }
        if (invalid_chars > 0) {
            std::cout << "[DEBUG runMinimap2] WARNING: Query " << q_name << " contains " 
                      << invalid_chars << " non-ACGTN characters (e.g. gaps or formatting symbols)!" << std::endl;
        }
        if (!q_seq.empty()) {
            std::cout << "[DEBUG runMinimap2] Query Prefix (100bp): " 
                      << q_seq.substr(0, std::min<size_t>(100, q_seq.length())) << std::endl;
        }

        mm_reg1_t* regs = mm_map(mi, static_cast<int>(q_seq.length()), q_seq.c_str(),
                                 &n_regs, tbuf, &mopt, q_name.c_str());
        std::cout << "[DEBUG runMinimap2] Query: " << q_name << " mapping complete. Found " << n_regs << " alignment regions." << std::endl;

        for (int j = 0; j < n_regs; ++j) {
            mm_reg1_t* aln_reg = &regs[j];

            if (write_paf && pafOut.is_open()) {
                std::string r_name = (aln_reg->rid >= 0 && aln_reg->rid < mi->n_seq) ? mi->seq[aln_reg->rid].name : "unknown_ref";
                uint32_t r_len     = (aln_reg->rid >= 0 && aln_reg->rid < mi->n_seq) ? mi->seq[aln_reg->rid].len : 0;

                pafOut << q_name << "\t"
                       << q_seq.length() << "\t"
                       << aln_reg->qs << "\t"
                       << aln_reg->qe << "\t"
                       << (aln_reg->rev ? "-" : "+") << "\t"
                       << r_name << "\t"
                       << r_len << "\t"
                       << aln_reg->rs << "\t"
                       << aln_reg->re << "\t"
                       << aln_reg->mlen << "\t"
                       << aln_reg->blen << "\t"
                       << aln_reg->mapq << "\t"
                       << "tp:A:" << (aln_reg->sam_pri ? "P" : "S") << "\n";
            }

            Alignment aln;
            aln.ID = 0;
            // 安全性檢查：防範 rid 越界讀取
            if (aln_reg->rid >= 0 && aln_reg->rid < mi->n_seq) {
                aln.refName = mi->seq[aln_reg->rid].name;
            } else {
                aln.refName = "unknown_ref";
            }
            aln.qryName = q_name;
            aln.refIdx = {aln_reg->rs, aln_reg->re};
            aln.qryIdx = {aln_reg->qs, aln_reg->qe};
            aln.inverse = (aln_reg->rev != 0);
            aln.primary = (aln_reg->sam_pri != 0);
            aln.valid = true;

            int totLen = 0, ins = 0, del = 0, numEdits = 0;
            int score = 0;

            if (aln_reg->p) {
                score = aln_reg->p->dp_max;
                numEdits = aln_reg->p->n_ambi;

                aln.CIGAR.reserve(aln_reg->p->n_cigar);
                for (uint32_t k = 0; k < aln_reg->p->n_cigar; ++k) {
                    uint32_t c32 = aln_reg->p->cigar[k];
                    int length = c32 >> 4;
                    char opType = "MIDNSHPE="[c32 & 0xf];

                    aln.CIGAR.push_back({length, opType});

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
                free(aln_reg->p); // 釋放 CIGAR payload 記憶體
            } else {
                // Chaining-only mode: 使用 block length 作為對齊長度
                totLen = aln_reg->blen;
            }

            aln.alnLength = totLen;
            aln.alnScore = score;
            aln.ins = ins;
            aln.del = del;
            aln.mis = std::max(0, numEdits - ins - del);

            alignments.push_back(std::move(aln));
        }
        free(regs);
    }

    // 6. 清理 Minimap2 記憶體
    std::cout << "[DEBUG runMinimap2] Cleaning up minimap2 buffers and index..." << std::endl;
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
    // ⚠️ 註解掉 mm_idx_reader_close(r) 以避免 zlib / fclose 對 EOF 流關閉時觸發 ASan 的 SEGV。
    // 這是一次性比對，少量的記憶體洩漏不會對程式碼有任何影響，且暫存檔依然會在下面被 std::remove 移除。
    // mm_idx_reader_close(r);

    // 🌟 Hybrid Fallback 註解 (CLI Debug 關閉)
    /*
    if (alignments.empty() && pafLines > 0) {
        std::cout << "[DEBUG runMinimap2] C-API returned 0 alignments, but CLI generated " 
                  << pafLines << " lines. Falling back to parsing CLI PAF..." << std::endl;
        alignments = parser::parseMinimap2PAF(option.tempDir + "/minimap2_cli_test.paf");
    }
    */

    // 清理產生的暫存 Reference 檔案
    std::remove((option.tempDir + "/temp_minimap2_ref.fa").c_str());

    // Sort alignments: Group by RefStart, then by QryStart
    std::cout << "[DEBUG runMinimap2] Sorting " << alignments.size() << " alignment results..." << std::endl;
    std::sort(alignments.begin(), alignments.end(), [](const auto &a, const auto &b) {
        if (a.refIdx.first == b.refIdx.first) return a.qryIdx.first < b.qryIdx.first;
        return (a.refIdx.first < b.refIdx.first);
    });

    auto minimap2_end = std::chrono::high_resolution_clock::now();
    option.minimap2_time += std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count();

    std::cout << "[DEBUG runMinimap2] Minimap2 API finished. Elapsed time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(minimap2_end - minimap2_start).count() << " ms" << std::endl;

    return alignments;
}

Alignments runMinimap2(StringPairs& ref, StringPairs& qry, std::string refName, std::string qryName, Option& option, bool needCigar, bool write_paf) {
    SequenceRefs refRefs;
    refRefs.reserve(ref.size());
    for (const auto& pair : ref) {
        refRefs.push_back({pair.first, pair.second});
    }
    SequenceRefs qryRefs;
    qryRefs.reserve(qry.size());
    for (const auto& pair : qry) {
        qryRefs.push_back({pair.first, pair.second});
    }
    return runMinimap2(refRefs, qryRefs, refName, qryName, option, needCigar, write_paf);
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
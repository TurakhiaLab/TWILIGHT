#include "medoid_realigner.hpp"
#include "cigar_util.hpp"
#include "global_alignment.hpp"
#include "message_manager.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>
#include <vector>
#include <string>

// Helper to append operations and merge adjacent identical ops
static inline void appendCigarOp(CigarString &cigar, int len, char op) {
    if (len <= 0) return;
    if (!cigar.empty() && cigar.back().second == op) {
        cigar.back().first += len;
    } else {
        cigar.push_back({len, op});
    }
}

// Count actual base matches (identity) given an alignment
static inline int countActualMatches(const std::string &ref, const std::string &qry, const CigarString &cigar) {
    int matches = 0;
    int r = 0, q = 0;
    for (const auto &op : cigar) {
        if (op.second == 'M' || op.second == '=' || op.second == 'X') {
            for (int k = 0; k < op.first; ++k) {
                if (r + k < static_cast<int>(ref.length()) && q + k < static_cast<int>(qry.length())) {
                    if (ref[r + k] == qry[q + k]) matches++;
                }
            }
            r += op.first;
            q += op.first;
        } else if (consumesRef(op.second)) {
            r += op.first;
        } else if (consumesQry(op.second)) {
            q += op.first;
        }
    }
    return matches;
}

// Perform Head/Tail bounded test for extreme length asymmetry
static CigarString realignAsymmetricWindow(
    const std::string &ref_sub,
    const std::string &qry_sub,
    int F_L, int F_R,
    int L_R, int L_Q)
{
    // If Qry is short (<= 5) and Ref is large (>= 30)
    if (L_Q <= 5 && L_R >= 30) {
        int head_ref_take = std::min(static_cast<int>(ref_sub.length()), F_L + 20);
        int head_qry_take = std::min(static_cast<int>(qry_sub.length()), F_L + L_Q);
        std::string ref_head = ref_sub.substr(0, head_ref_take);
        std::string qry_head = qry_sub.substr(0, head_qry_take);

        CigarString head_cigar = runGlobalAlignment(ref_head, qry_head);

        int head_matches = countActualMatches(ref_head, qry_head, head_cigar);
        int head_ref_consumed = 0;
        int head_qry_consumed = 0;
        for (const auto &op : head_cigar) {
            if (consumesRef(op.second)) head_ref_consumed += op.first;
            if (consumesQry(op.second)) head_qry_consumed += op.first;
        }

        // If head matched extra bases beyond the left flank
        if (head_matches >= F_L + L_Q - 1 && head_qry_consumed == F_L + L_Q) {
            CigarString res = head_cigar;
            int rem_ref = static_cast<int>(ref_sub.length()) - head_ref_consumed - F_R;
            if (rem_ref > 0) appendCigarOp(res, rem_ref, 'D');
            if (F_R > 0) appendCigarOp(res, F_R, 'M');
            return compressCigar(res);
        }

        // Tail alignment test
        int tail_ref_take = std::min(static_cast<int>(ref_sub.length()), F_R + 20);
        int tail_qry_take = std::min(static_cast<int>(qry_sub.length()), F_R + L_Q);
        std::string ref_tail = ref_sub.substr(ref_sub.length() - tail_ref_take);
        std::string qry_tail = qry_sub.substr(qry_sub.length() - tail_qry_take);

        CigarString tail_cigar = runGlobalAlignment(ref_tail, qry_tail);
        int tail_matches = countActualMatches(ref_tail, qry_tail, tail_cigar);
        int tail_ref_consumed = 0;
        int tail_qry_consumed = 0;
        for (const auto &op : tail_cigar) {
            if (consumesRef(op.second)) tail_ref_consumed += op.first;
            if (consumesQry(op.second)) tail_qry_consumed += op.first;
        }

        if (tail_matches >= F_R + L_Q - 1 && tail_qry_consumed == F_R + L_Q) {
            CigarString res;
            if (F_L > 0) appendCigarOp(res, F_L, 'M');
            int rem_ref = static_cast<int>(ref_sub.length()) - tail_ref_consumed - F_L;
            if (rem_ref > 0) appendCigarOp(res, rem_ref, 'D');
            for (const auto &op : tail_cigar) {
                appendCigarOp(res, op.first, op.second);
            }
            return compressCigar(res);
        }

        // Neither head nor tail matches cleanly: keep pure indel to avoid shattering
        CigarString res;
        if (F_L > 0) appendCigarOp(res, F_L, 'M');
        if (L_Q > 0) appendCigarOp(res, L_Q, 'I');
        if (L_R > 0) appendCigarOp(res, L_R, 'D');
        if (F_R > 0) appendCigarOp(res, F_R, 'M');
        return compressCigar(res);
    }

    // Symmetric case: Ref is short (<= 5) and Qry is large (>= 30)
    if (L_R <= 5 && L_Q >= 30) {
        int head_ref_take = std::min(static_cast<int>(ref_sub.length()), F_L + L_R);
        int head_qry_take = std::min(static_cast<int>(qry_sub.length()), F_L + 20);
        std::string ref_head = ref_sub.substr(0, head_ref_take);
        std::string qry_head = qry_sub.substr(0, head_qry_take);

        CigarString head_cigar = runGlobalAlignment(ref_head, qry_head);
        int head_matches = countActualMatches(ref_head, qry_head, head_cigar);
        int head_ref_consumed = 0;
        int head_qry_consumed = 0;
        for (const auto &op : head_cigar) {
            if (consumesRef(op.second)) head_ref_consumed += op.first;
            if (consumesQry(op.second)) head_qry_consumed += op.first;
        }

        if (head_matches >= F_L + L_R - 1 && head_ref_consumed == F_L + L_R) {
            CigarString res = head_cigar;
            int rem_qry = static_cast<int>(qry_sub.length()) - head_qry_consumed - F_R;
            if (rem_qry > 0) appendCigarOp(res, rem_qry, 'I');
            if (F_R > 0) appendCigarOp(res, F_R, 'M');
            return compressCigar(res);
        }

        int tail_ref_take = std::min(static_cast<int>(ref_sub.length()), F_R + L_R);
        int tail_qry_take = std::min(static_cast<int>(qry_sub.length()), F_R + 20);
        std::string ref_tail = ref_sub.substr(ref_sub.length() - tail_ref_take);
        std::string qry_tail = qry_sub.substr(qry_sub.length() - tail_qry_take);

        CigarString tail_cigar = runGlobalAlignment(ref_tail, qry_tail);
        int tail_matches = countActualMatches(ref_tail, qry_tail, tail_cigar);
        int tail_ref_consumed = 0;
        int tail_qry_consumed = 0;
        for (const auto &op : tail_cigar) {
            if (consumesRef(op.second)) tail_ref_consumed += op.first;
            if (consumesQry(op.second)) tail_qry_consumed += op.first;
        }

        if (tail_matches >= F_R + L_R - 1 && tail_ref_consumed == F_R + L_R) {
            CigarString res;
            if (F_L > 0) appendCigarOp(res, F_L, 'M');
            int rem_qry = static_cast<int>(qry_sub.length()) - tail_qry_consumed - F_L;
            if (rem_qry > 0) appendCigarOp(res, rem_qry, 'I');
            for (const auto &op : tail_cigar) {
                appendCigarOp(res, op.first, op.second);
            }
            return compressCigar(res);
        }

        CigarString res;
        if (F_L > 0) appendCigarOp(res, F_L, 'M');
        if (L_R > 0) appendCigarOp(res, L_R, 'D');
        if (L_Q > 0) appendCigarOp(res, L_Q, 'I');
        if (F_R > 0) appendCigarOp(res, F_R, 'M');
        return compressCigar(res);
    }

    // Default fallback to standard global alignment
    return runGlobalAlignment(ref_sub, qry_sub);
}

Alignments realignMedoidAlignmentsToConsensus(
    const Alignments &medAlignments,
    BlockSet *refBlockSet,
    BlockSet *qryBlockSet,
    const Option &option)
{
    Alignments consensusAlignments;
    consensusAlignments.reserve(medAlignments.size());

    const std::string &ref_cons_all = refBlockSet->getAncestralSequence();
    const std::string &qry_cons_all = qryBlockSet->getAncestralSequence();

    for (const auto &med_aln : medAlignments) {
        if (!med_aln.valid || med_aln.CIGAR.empty()) continue;

        int m_r_start = med_aln.refIdx.first;
        int m_r_end   = med_aln.refIdx.second;
        int m_q_start = med_aln.qryIdx.first;
        int m_q_end   = med_aln.qryIdx.second;
        bool is_rev   = med_aln.inverse;

        if (m_r_start >= m_r_end || m_q_start >= m_q_end) continue;

        // 1. Liftover coordinates from Medoid to Consensus space
        int c_r_start = refBlockSet->medoidToConsensusOffset(m_r_start);
        int c_r_last  = refBlockSet->medoidToConsensusOffset(m_r_end - 1);
        int c_r_end   = c_r_last + 1;

        int c_q_start = qryBlockSet->medoidToConsensusOffset(m_q_start);
        int c_q_last  = qryBlockSet->medoidToConsensusOffset(m_q_end - 1);
        int c_q_end   = c_q_last + 1;

        int expected_ref_bases = c_r_end - c_r_start;
        int expected_qry_bases = c_q_end - c_q_start;

        if (expected_ref_bases <= 0 || expected_qry_bases <= 0) continue;

        // 2. Initial Liftover of CIGAR base-by-base
        CigarString raw_consensus_cigar;
        int curr_m_r = m_r_start;
        int curr_m_q = is_rev ? (m_q_end - 1) : m_q_start;

        int last_c_r = c_r_start - 1;
        int last_c_q = is_rev ? c_q_end : (c_q_start - 1);

        for (const auto &op : med_aln.CIGAR) {
            int op_len = op.first;
            char op_type = op.second;

            for (int step = 0; step < op_len; ++step) {
                if (op_type == 'M' || op_type == '=' || op_type == 'X') {
                    int pos_c_r = refBlockSet->medoidToConsensusOffset(curr_m_r);
                    int pos_c_q = qryBlockSet->medoidToConsensusOffset(curr_m_q);

                    int r_gap = pos_c_r - (last_c_r + 1);
                    int q_gap = is_rev ? (last_c_q - (pos_c_q + 1)) : (pos_c_q - (last_c_q + 1));

                    if (r_gap > 0) appendCigarOp(raw_consensus_cigar, r_gap, 'D');
                    if (q_gap > 0) appendCigarOp(raw_consensus_cigar, q_gap, 'I');

                    appendCigarOp(raw_consensus_cigar, 1, op_type);

                    last_c_r = pos_c_r;
                    last_c_q = pos_c_q;
                    curr_m_r++;
                    if (is_rev) curr_m_q--; else curr_m_q++;
                }
                else if (op_type == 'D' || op_type == 'N') {
                    int pos_c_r = refBlockSet->medoidToConsensusOffset(curr_m_r);
                    int r_gap = pos_c_r - (last_c_r + 1);

                    if (r_gap > 0) appendCigarOp(raw_consensus_cigar, r_gap, 'D');
                    appendCigarOp(raw_consensus_cigar, 1, op_type);

                    last_c_r = pos_c_r;
                    curr_m_r++;
                }
                else if (op_type == 'I' || op_type == 'S') {
                    int pos_c_q = qryBlockSet->medoidToConsensusOffset(curr_m_q);
                    int q_gap = is_rev ? (last_c_q - (pos_c_q + 1)) : (pos_c_q - (last_c_q + 1));

                    if (q_gap > 0) appendCigarOp(raw_consensus_cigar, q_gap, 'I');
                    appendCigarOp(raw_consensus_cigar, 1, op_type);

                    last_c_q = pos_c_q;
                    if (is_rev) curr_m_q--; else curr_m_q++;
                }
            }
        }

        // Account for any trailing deleted bases up to c_r_end / c_q_end
        int trailing_r = (c_r_end - 1) - last_c_r;
        int trailing_q = is_rev ? (last_c_q - c_q_start) : ((c_q_end - 1) - last_c_q);
        if (trailing_r > 0) appendCigarOp(raw_consensus_cigar, trailing_r, 'D');
        if (trailing_q > 0) appendCigarOp(raw_consensus_cigar, trailing_q, 'I');

        raw_consensus_cigar = compressCigar(raw_consensus_cigar);

        // 3. Scan for dirty regions (adjacent or clustered indels) for Flank-Guarded Realignment
        CigarString realigned_cigar;
        int curr_r_pos = c_r_start;
        int curr_q_pos = is_rev ? c_q_end : c_q_start;

        size_t i = 0;
        while (i < raw_consensus_cigar.size()) {
            // Check if current op is an Indel
            if (raw_consensus_cigar[i].second != 'D' && raw_consensus_cigar[i].second != 'I') {
                appendCigarOp(realigned_cigar, raw_consensus_cigar[i].first, raw_consensus_cigar[i].second);
                if (consumesRef(raw_consensus_cigar[i].second)) curr_r_pos += raw_consensus_cigar[i].first;
                if (consumesQry(raw_consensus_cigar[i].second)) {
                    if (is_rev) curr_q_pos -= raw_consensus_cigar[i].first;
                    else curr_q_pos += raw_consensus_cigar[i].first;
                }
                i++;
                continue;
            }

            // Found an Indel: expand dirty window to include any neighboring indels or micro-matches (<= 10bp)
            size_t win_start_op = i;
            size_t win_end_op = i;
            int win_del = 0;
            int win_ins = 0;

            while (win_end_op < raw_consensus_cigar.size()) {
                char op_c = raw_consensus_cigar[win_end_op].second;
                int op_l = raw_consensus_cigar[win_end_op].first;

                if (op_c == 'D') {
                    win_del += op_l;
                    win_end_op++;
                } else if (op_c == 'I') {
                    win_ins += op_l;
                    win_end_op++;
                } else if (op_l <= 10 && (win_end_op + 1 < raw_consensus_cigar.size()) &&
                           (raw_consensus_cigar[win_end_op + 1].second == 'D' || raw_consensus_cigar[win_end_op + 1].second == 'I')) {
                    // Small intervening match between indels: include in dirty window
                    win_end_op++;
                } else {
                    break;
                }
            }

            // If only pure deletion or pure insertion, no mutual realignment needed
            if (win_del == 0 || win_ins == 0) {
                for (size_t k = win_start_op; k < win_end_op; ++k) {
                    appendCigarOp(realigned_cigar, raw_consensus_cigar[k].first, raw_consensus_cigar[k].second);
                    if (consumesRef(raw_consensus_cigar[k].second)) curr_r_pos += raw_consensus_cigar[k].first;
                    if (consumesQry(raw_consensus_cigar[k].second)) {
                        if (is_rev) curr_q_pos -= raw_consensus_cigar[k].first;
                        else curr_q_pos += raw_consensus_cigar[k].first;
                    }
                }
                i = win_end_op;
                continue;
            }

            // Both Del and Ins exist in this window: determine flanking Match bounds
            int F_L = 0;
            if (!realigned_cigar.empty() && realigned_cigar.back().second == 'M') {
                F_L = std::min(15, realigned_cigar.back().first);
                realigned_cigar.back().first -= F_L;
                if (realigned_cigar.back().first == 0) realigned_cigar.pop_back();
            }

            int F_R = 0;
            if (win_end_op < raw_consensus_cigar.size() && raw_consensus_cigar[win_end_op].second == 'M') {
                F_R = std::min(15, raw_consensus_cigar[win_end_op].first);
                raw_consensus_cigar[win_end_op].first -= F_R;
            }

            // Calculate exact slice spans in Consensus Ref and Qry
            int orig_win_r_len = 0;
            int orig_win_q_len = 0;
            for (size_t k = win_start_op; k < win_end_op; ++k) {
                if (consumesRef(raw_consensus_cigar[k].second)) orig_win_r_len += raw_consensus_cigar[k].first;
                if (consumesQry(raw_consensus_cigar[k].second)) orig_win_q_len += raw_consensus_cigar[k].first;
            }

            int win_r_start = curr_r_pos - F_L;
            int total_r_len = orig_win_r_len + F_L + F_R;

            int win_q_start, total_q_len;
            if (!is_rev) {
                win_q_start = curr_q_pos - F_L;
                total_q_len = orig_win_q_len + F_L + F_R;
            } else {
                win_q_start = (curr_q_pos - orig_win_q_len) - F_R;
                total_q_len = orig_win_q_len + F_L + F_R;
            }

            // Extract consensus substrings
            bool valid_bounds = (win_r_start >= 0 && win_r_start + total_r_len <= static_cast<int>(ref_cons_all.length())) &&
                                (win_q_start >= 0 && win_q_start + total_q_len <= static_cast<int>(qry_cons_all.length()));

            if (valid_bounds && (total_r_len + total_q_len <= 500)) {
                std::string ref_sub = ref_cons_all.substr(win_r_start, total_r_len);
                std::string qry_sub;
                if (!is_rev) {
                    qry_sub = qry_cons_all.substr(win_q_start, total_q_len);
                } else {
                    qry_sub = getReverseComplement(qry_cons_all.substr(win_q_start, total_q_len));
                }

                CigarString local_res;
                if (std::min(orig_win_r_len, orig_win_q_len) <= 5 && std::max(orig_win_r_len, orig_win_q_len) >= 30) {
                    local_res = realignAsymmetricWindow(ref_sub, qry_sub, F_L, F_R, orig_win_r_len, orig_win_q_len);
                } else {
                    local_res = runGlobalAlignment(ref_sub, qry_sub);
                }

                // Verify base counts of local realignment
                int local_r_count = 0, local_q_count = 0;
                for (const auto &op : local_res) {
                    if (consumesRef(op.second)) local_r_count += op.first;
                    if (consumesQry(op.second)) local_q_count += op.first;
                }

                if (local_r_count == total_r_len && local_q_count == total_q_len) {
                    for (const auto &op : local_res) {
                        appendCigarOp(realigned_cigar, op.first, op.second);
                    }
                } else {
                    // Fallback to initial ops if count mismatch
                    if (F_L > 0) appendCigarOp(realigned_cigar, F_L, 'M');
                    for (size_t k = win_start_op; k < win_end_op; ++k) {
                        appendCigarOp(realigned_cigar, raw_consensus_cigar[k].first, raw_consensus_cigar[k].second);
                    }
                    if (F_R > 0) appendCigarOp(realigned_cigar, F_R, 'M');
                }
            } else {
                // Fallback without realignment for large/out-of-bounds windows
                if (F_L > 0) appendCigarOp(realigned_cigar, F_L, 'M');
                for (size_t k = win_start_op; k < win_end_op; ++k) {
                    appendCigarOp(realigned_cigar, raw_consensus_cigar[k].first, raw_consensus_cigar[k].second);
                }
                if (F_R > 0) appendCigarOp(realigned_cigar, F_R, 'M');
            }

            curr_r_pos += orig_win_r_len + F_R;
            if (!is_rev) curr_q_pos += orig_win_q_len + F_R;
            else curr_q_pos -= (orig_win_q_len + F_R);

            i = win_end_op;
            if (i < raw_consensus_cigar.size() && raw_consensus_cigar[i].first == 0) i++;
        }

        realigned_cigar = compressCigar(realigned_cigar);

        // 4. Invariant Verification: Total bases consumed must exactly match interval lengths
        int final_ref_count = 0;
        int final_qry_count = 0;
        int ins_count = 0, del_count = 0;
        for (const auto &op : realigned_cigar) {
            if (consumesRef(op.second)) final_ref_count += op.first;
            if (consumesQry(op.second)) final_qry_count += op.first;
            if (op.second == 'I') ins_count += op.first;
            if (op.second == 'D') del_count += op.first;
        }

        Alignment final_aln = med_aln;
        final_aln.refIdx = {c_r_start, c_r_end};
        final_aln.qryIdx = {c_q_start, c_q_end};

        if (final_ref_count == expected_ref_bases && final_qry_count == expected_qry_bases) {
            final_aln.CIGAR = std::move(realigned_cigar);
            final_aln.ins = ins_count;
            final_aln.del = del_count;
            final_aln.updateAlnLength();
        } else {
            // Safety fallback to raw liftover if any discrepancy is detected
            mga::msg::warn("Realignment verification failed for alignment %d: ref (%d vs %d), qry (%d vs %d). Falling back to initial liftover.\n",
                           med_aln.ID, final_ref_count, expected_ref_bases, final_qry_count, expected_qry_bases);
            final_aln.CIGAR = std::move(raw_consensus_cigar);
            final_aln.updateAlnLength();
        }

        consensusAlignments.push_back(std::move(final_aln));
    }

    return consensusAlignments;
}

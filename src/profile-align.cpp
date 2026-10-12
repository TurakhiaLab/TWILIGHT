/*
 * profile-align.cpp - Profile-Profile Alignment Engine
 *
 * This module implements dynamic programming for progressive profile alignment
 * with terminal gap costs and consistency score guidance.
 *
 * Acknowledgments & Attribution:
 * The gap penalty model (terminal vs. internal gaps) and compressed DP matrix
 * representation implemented here are adapted from the FAMSA algorithm:
 *   Deorowicz, S., Debudaj-Grabysz, A. & Gudyś, A.
 *   "FAMSA: Fast and accurate multiple sequence alignment of huge protein families."
 *   Sci Rep 6, 33964 (2016). https://doi.org/10.1038/srep33964
 *
 * Extensions:
 * Integrated with TWILIGHT's multi-level progressive alignment architecture and
 * positional consistency weighting library (MAFFT / TCS consistency bonuses).
 */

#include "align-method.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {
    using namespace msa;

    struct ProfileDPMatrix {
        size_t n_rows;
        size_t n_cols;
        std::vector<uint8_t> raw_data;

        ProfileDPMatrix(size_t rows, size_t cols) : n_rows(rows), n_cols(cols), raw_data(rows * cols, 0) {}

        inline int get_dir_D(size_t r, size_t c) const {
            return raw_data[r * n_cols + c] & 0x03;
        }
        inline int get_dir_H(size_t r, size_t c) const {
            return (raw_data[r * n_cols + c] >> 2) & 0x03;
        }
        inline int get_dir_V(size_t r, size_t c) const {
            return (raw_data[r * n_cols + c] >> 4) & 0x03;
        }

        inline void set_dir_D(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xFC) | (dir & 0x03);
        }
        inline void set_dir_H(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xF3) | ((dir & 0x03) << 2);
        }
        inline void set_dir_V(size_t r, size_t c, int dir) {
            uint8_t& p = raw_data[r * n_cols + c];
            p = (p & 0xCF) | ((dir & 0x03) << 4);
        }
        inline void set_dir_all(size_t r, size_t c, int dir) {
            uint8_t x = dir & 0x03;
            raw_data[r * n_cols + c] = x | (x << 2) | (x << 4);
        }
    };

    inline void solve_gaps_start(
        size_t source_col_id, size_t prof_width, float prof_size,
        const msa::AlnProfile& profile,
        float& n_gap_open, float& n_gap_ext, float& n_gap_term_open, float& n_gap_term_ext)
    {
        n_gap_open = n_gap_ext = n_gap_term_open = n_gap_term_ext = 0.0f;
        if (source_col_id >= prof_width) {
            const float* source_col = profile.get_counters(source_col_id);
            float cnt = source_col[GAP_TERM_OPEN_IDX] + source_col[GAP_TERM_EXT_IDX];
            n_gap_term_ext = cnt;
            n_gap_term_open += (prof_size - cnt);
        } else {
            const float* next_col = profile.get_counters(source_col_id + 1);
            n_gap_term_open += next_col[GAP_TERM_OPEN_IDX];

            const float* source_col = profile.get_counters(source_col_id);
            n_gap_term_ext += source_col[GAP_TERM_OPEN_IDX];
            n_gap_term_ext += source_col[GAP_TERM_EXT_IDX];

            n_gap_ext = source_col[GAP_OPEN_IDX];
            n_gap_ext += source_col[GAP_EXT_IDX];

            n_gap_open = prof_size - n_gap_ext - n_gap_term_open - n_gap_term_ext;
        }
    }

    inline void solve_gaps_cont(
        size_t source_col_id, size_t prof_width, float prof_size,
        const msa::AlnProfile& profile,
        float& n_gap_ext, float& n_gap_term_ext)
    {
        n_gap_ext = n_gap_term_ext = 0.0f;
        if (source_col_id == prof_width) {
            n_gap_term_ext = prof_size;
            n_gap_ext = 0.0f;
        } else {
            const float* next_col = profile.get_counters(source_col_id + 1);
            n_gap_term_ext = next_col[GAP_TERM_OPEN_IDX];

            const float* source_col = profile.get_counters(source_col_id);
            n_gap_term_ext += source_col[GAP_TERM_OPEN_IDX];
            n_gap_term_ext += source_col[GAP_TERM_EXT_IDX];

            n_gap_ext = prof_size - n_gap_term_ext;
        }
    }
} // anonymous namespace

std::vector<int8_t> msa::alignProfile(
    const AlnProfile& profile1,
    const AlnProfile& profile2,
    const Params& param,
    const std::vector<std::vector<float>>* consistencyTable,
    float consistencyWeight,
    int bandWidth)
{
    size_t prof1_width = profile1.width;
    size_t prof2_width = profile2.width;

    if (prof1_width == 0 && prof2_width == 0) return {};
    if (prof1_width == 0) return std::vector<int8_t>(prof2_width, 1);
    if (prof2_width == 0) return std::vector<int8_t>(prof1_width, 2);

    float prof1_card = (profile1.totalWeight > 0.0f) ? profile1.totalWeight : (profile1.numSeqs > 0 ? static_cast<float>(profile1.numSeqs) : 1.0f);
    float prof2_card = (profile2.totalWeight > 0.0f) ? profile2.totalWeight : (profile2.numSeqs > 0 ? static_cast<float>(profile2.numSeqs) : 1.0f);

    float gap_open = param.gapOpen;
    float gap_ext = param.gapExtend;
    float gap_term_open = (param.gapTerminal != 0.0f) ? param.gapTerminal : (gap_open * (0.66f / 14.85f));
    float gap_term_ext = (param.gapTerminal != 0.0f) ? param.gapTerminal : (gap_ext * (0.66f / 1.25f));

    // Dynamic Banded Alignment Determination:
    // - bandWidth == 0: Force full-matrix
    // - bandWidth > 0: Explicit user-specified band width
    // - bandWidth < 0: Auto adaptive mode (full-matrix if both widths <= 1500, dynamic adaptive band if larger)
    constexpr size_t FULL_MATRIX_THRESHOLD = 1500;
    bool use_banded = false;
    size_t half_band = prof2_width;

    if (bandWidth == 0) {
        use_banded = false;
    } else if (bandWidth > 0) {
        use_banded = true;
        half_band = static_cast<size_t>(std::max(1, bandWidth / 2));
    } else {
        if (prof1_width > FULL_MATRIX_THRESHOLD || prof2_width > FULL_MATRIX_THRESHOLD) {
            size_t min_len = std::min(prof1_width, prof2_width);
            size_t diff_len = (prof1_width > prof2_width) ? (prof1_width - prof2_width) : (prof2_width - prof1_width);
            // Adaptive half-band: diff_len + 4 * sqrt(min_len) + 10% min_len buffer + 64
            half_band = diff_len + static_cast<size_t>(4.0f * std::sqrt(static_cast<float>(min_len))) + 64;
            half_band = std::max(half_band, diff_len + static_cast<size_t>(0.10f * min_len));
            if (half_band < prof2_width) {
                use_banded = true;
            }
        }
    }

    ProfileDPMatrix matrix(prof1_width + 1, prof2_width + 1);

    const float infty = 1e20f;
    struct dp_row_elem_t {
        float D, H, V;
    };
    std::vector<dp_row_elem_t> prev_row(prof2_width + 2, {-infty, -infty, -infty});
    std::vector<dp_row_elem_t> curr_row(prof2_width + 2, {-infty, -infty, -infty});

    struct dp_gap_costs {
        float open, ext, term_open, term_ext;
    };
    struct dp_gap_corrections {
        float n_gap_start_open, n_gap_start_ext, n_gap_start_term_open, n_gap_start_term_ext;
        float n_gap_cont_ext, n_gap_cont_term_ext;
    };

    std::vector<dp_gap_costs> prof2_gaps(prof2_width + 1);
    for (size_t j = 0; j <= prof2_width; ++j) {
        const float* s2 = profile2.get_scores(j);
        prof2_gaps[j].open      = s2[GAP_OPEN_IDX];
        prof2_gaps[j].ext       = s2[GAP_EXT_IDX];
        prof2_gaps[j].term_open = s2[GAP_TERM_OPEN_IDX];
        prof2_gaps[j].term_ext  = s2[GAP_TERM_EXT_IDX];
    }

    std::vector<dp_gap_corrections> gap_corrections(prof2_width + 1);
    std::vector<float> n_gaps_prof2_to_change(prof2_width + 1, 0.0f);
    std::vector<float> n_gaps_prof2_term_to_change(prof2_width + 1, 0.0f);
    std::vector<float> gaps_prof2_change(prof2_width + 1, 0.0f);

    for (size_t j = 1; j <= prof2_width; ++j) {
        solve_gaps_start(j, prof2_width, prof2_card, profile2,
            gap_corrections[j].n_gap_start_open, gap_corrections[j].n_gap_start_ext,
            gap_corrections[j].n_gap_start_term_open, gap_corrections[j].n_gap_start_term_ext);

        solve_gaps_cont(j, prof2_width, prof2_card, profile2,
            gap_corrections[j].n_gap_cont_ext, gap_corrections[j].n_gap_cont_term_ext);

        const float* c2 = profile2.get_counters(j);
        n_gaps_prof2_to_change[j]      = c2[GAP_OPEN_IDX];
        n_gaps_prof2_term_to_change[j] = c2[GAP_TERM_OPEN_IDX];

        gaps_prof2_change[j] = n_gaps_prof2_to_change[j] * (gap_ext - gap_open) +
                               n_gaps_prof2_term_to_change[j] * (gap_term_ext - gap_term_open);
    }

    // Boundary conditions for row 0
    prev_row[0].D = 0.0f;
    prev_row[0].H = -infty;
    prev_row[0].V = -infty;

    size_t row0_end = use_banded ? std::min(prof2_width, half_band) : prof2_width;
    if (prof2_width >= 1 && row0_end >= 1) {
        prev_row[1].D = -infty;
        prev_row[1].H = prev_row[0].D + prof2_gaps[1].term_open * prof1_card;
        prev_row[1].V = -infty;
        matrix.set_dir_all(0, 1, 1);
    }

    for (size_t j = 2; j <= row0_end; ++j) {
        prev_row[j].D = -infty;
        prev_row[j].H = prev_row[j - 1].H + prof2_gaps[j].term_ext * prof1_card;
        prev_row[j].V = -infty;
        matrix.set_dir_all(0, j, 1);
    }
    for (size_t j = row0_end + 1; j <= prof2_width; ++j) {
        prev_row[j].D = -infty;
        prev_row[j].H = -infty;
        prev_row[j].V = -infty;
    }
    prev_row[prof2_width].H = -infty;

    // Calculate matrix interior
    for (size_t i = 1; i <= prof1_width; ++i) {
        const float* s1_i = profile1.get_scores(i);
        float prof1_gap_open_curr      = s1_i[GAP_OPEN_IDX];
        float prof1_gap_term_open_curr = s1_i[GAP_TERM_OPEN_IDX];
        float prof1_gap_ext_curr       = s1_i[GAP_EXT_IDX];
        float prof1_gap_term_ext_curr  = s1_i[GAP_TERM_EXT_IDX];

        size_t j_start = 1;
        size_t j_end = prof2_width;
        if (use_banded) {
            size_t j_center = (i * 2 * prof2_width + prof1_width) / (2 * prof1_width);
            j_start = (j_center > half_band) ? std::max((size_t)1, j_center - half_band) : (size_t)1;
            j_end = std::min(prof2_width, j_center + half_band);
        }

        if (j_start == 1) {
            curr_row[0].D = -infty;
            curr_row[0].H = -infty;
            matrix.set_dir_all(i, 0, 2);

            if (i < prof1_width) {
                if (i == 1)
                    curr_row[0].V = std::max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_open_curr * prof2_card;
                else
                    curr_row[0].V = std::max(prev_row[0].D, prev_row[0].V) + prof1_gap_term_ext_curr * prof2_card;
            } else {
                curr_row[0].V = -infty;
            }
        } else {
            curr_row[j_start - 1].D = -infty;
            curr_row[j_start - 1].H = -infty;
            curr_row[j_start - 1].V = -infty;
        }

        if (j_end < prof2_width) {
            curr_row[j_end + 1].D = -infty;
            curr_row[j_end + 1].H = -infty;
            curr_row[j_end + 1].V = -infty;
        }

        const float* c1_i = profile1.get_counters(i);
        std::array<std::pair<int, float>, 32> col1;
        size_t col1_size = 0;
        float col1_n_non_gaps = 0.0f;
        for (int k = 0; k < NO_SYMBOLS; ++k) {
            float cnt = c1_i[k];
            if (cnt > 0.0f) {
                col1[col1_size++] = {k, cnt};
                if (k < param.matrixSize)
                    col1_n_non_gaps += cnt;
            }
        }

        float n_gap_prof1_start_open = 0.0f, n_gap_prof1_start_ext = 0.0f;
        float n_gap_prof1_start_term_open = 0.0f, n_gap_prof1_start_term_ext = 0.0f;
        solve_gaps_start(i, prof1_width, prof1_card, profile1,
            n_gap_prof1_start_open, n_gap_prof1_start_ext, n_gap_prof1_start_term_open, n_gap_prof1_start_term_ext);

        float n_gap_prof1_cont_ext = 0.0f, n_gap_prof1_cont_term_ext = 0.0f;
        solve_gaps_cont(i, prof1_width, prof1_card, profile1,
            n_gap_prof1_cont_ext, n_gap_prof1_cont_term_ext);

        float n_gaps_prof1_to_change      = c1_i[GAP_OPEN_IDX];
        float n_gaps_prof1_term_to_change = c1_i[GAP_TERM_OPEN_IDX];

        for (size_t j = j_start; j <= j_end; ++j) {
            const float* scores2_column = profile2.get_scores(j);

            float t = 0.0f;
            for (size_t k = 0; k < col1_size; ++k) {
                t += col1[k].second * scores2_column[col1[k].first];
            }

            // Consistency bonus if available (cVal is per-pair average, so multiply by prof1_card * prof2_card to match FAMSA cardinality)
            if (consistencyTable && i - 1 < consistencyTable->size() && j - 1 < (*consistencyTable)[i - 1].size()) {
                float cVal = (*consistencyTable)[i - 1][j - 1];
                if (cVal > 0.0f) {
                    t += cVal * consistencyWeight * prof1_card * prof2_card;
                }
            }

            // State D
            float t_D = prev_row[j - 1].D + t;

            float t_H = prev_row[j - 1].H;
            if (n_gaps_prof1_to_change > 0.0f || n_gaps_prof1_term_to_change > 0.0f) {
                float delta_t = n_gaps_prof1_to_change * (scores2_column[GAP_EXT_IDX] - scores2_column[GAP_OPEN_IDX]) +
                                n_gaps_prof1_term_to_change * (scores2_column[GAP_TERM_EXT_IDX] - scores2_column[GAP_TERM_OPEN_IDX]);
                t_H = t_H + t + delta_t;
            } else {
                t_H = t_H + t;
            }

            float t_V = prev_row[j - 1].V + t + gaps_prof2_change[j] * col1_n_non_gaps;

            if (t_D > t_H && t_D > t_V) {
                curr_row[j].D = t_D;
                matrix.set_dir_D(i, j, 0);
            } else if (t_H > t_V) {
                curr_row[j].D = t_H;
                matrix.set_dir_D(i, j, 1);
            } else {
                curr_row[j].D = t_V;
                matrix.set_dir_D(i, j, 2);
            }

            // State H
            float gap_corr_H = prof2_gaps[j].open * n_gap_prof1_start_open +
                               prof2_gaps[j].ext * n_gap_prof1_start_ext +
                               prof2_gaps[j].term_open * n_gap_prof1_start_term_open +
                               prof2_gaps[j].term_ext * n_gap_prof1_start_term_ext;

            float t_D_H = curr_row[j - 1].D + gap_corr_H;
            float t_H_H = curr_row[j - 1].H + prof2_gaps[j].ext * n_gap_prof1_cont_ext +
                                             prof2_gaps[j].term_ext * n_gap_prof1_cont_term_ext;

            if (i > 1 && j > 1) {
                float t_V_H = curr_row[j - 1].V + gap_corr_H;
                if (t_D_H > t_H_H && t_D_H > t_V_H) {
                    curr_row[j].H = t_D_H;
                    matrix.set_dir_H(i, j, 0);
                } else if (t_V_H > t_H_H) {
                    curr_row[j].H = t_V_H;
                    matrix.set_dir_H(i, j, 2);
                } else {
                    curr_row[j].H = t_H_H;
                    matrix.set_dir_H(i, j, 1);
                }
            } else {
                if (t_D_H > t_H_H) {
                    curr_row[j].H = t_D_H;
                    matrix.set_dir_H(i, j, 0);
                } else {
                    curr_row[j].H = t_H_H;
                    matrix.set_dir_H(i, j, 1);
                }
            }

            // State V
            float gap_corr_V = prof1_gap_open_curr * gap_corrections[j].n_gap_start_open +
                               prof1_gap_ext_curr * gap_corrections[j].n_gap_start_ext +
                               prof1_gap_term_open_curr * gap_corrections[j].n_gap_start_term_open +
                               prof1_gap_term_ext_curr * gap_corrections[j].n_gap_start_term_ext;

            float t_D_V = prev_row[j].D + gap_corr_V;
            float t_V_V = prev_row[j].V + prof1_gap_ext_curr * gap_corrections[j].n_gap_cont_ext +
                                         prof1_gap_term_ext_curr * gap_corrections[j].n_gap_cont_term_ext;

            if (i > 1 && j > 1) {
                float t_H_V = prev_row[j].H + gap_corr_V;
                if (t_D_V > t_H_V && t_D_V > t_V_V) {
                    curr_row[j].V = t_D_V;
                    matrix.set_dir_V(i, j, 0);
                } else if (t_H_V > t_V_V) {
                    curr_row[j].V = t_H_V;
                    matrix.set_dir_V(i, j, 1);
                } else {
                    curr_row[j].V = t_V_V;
                    matrix.set_dir_V(i, j, 2);
                }
            } else {
                if (t_D_V > t_V_V) {
                    curr_row[j].V = t_D_V;
                    matrix.set_dir_V(i, j, 0);
                } else {
                    curr_row[j].V = t_V_V;
                    matrix.set_dir_V(i, j, 2);
                }
            }
        }
        curr_row.swap(prev_row);
    }

    // Traceback
    std::vector<int8_t> path;
    path.reserve(prof1_width + prof2_width);

    int cur_dir = 0;
    if (prev_row[prof2_width].D >= prev_row[prof2_width].H && prev_row[prof2_width].D >= prev_row[prof2_width].V)
        cur_dir = 0;
    else if (prev_row[prof2_width].H > prev_row[prof2_width].V)
        cur_dir = 1;
    else
        cur_dir = 2;

    size_t i = prof1_width, j = prof2_width;
    while (i > 0 || j > 0) {
        if (cur_dir == 0) {
            path.push_back(0); // Match
            cur_dir = matrix.get_dir_D(i, j);
            --i; --j;
        } else if (cur_dir == 1) {
            path.push_back(1); // Gap in Profile 1 (Query residue)
            cur_dir = matrix.get_dir_H(i, j);
            --j;
        } else {
            path.push_back(2); // Gap in Profile 2 (Ref residue)
            cur_dir = matrix.get_dir_V(i, j);
            --i;
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

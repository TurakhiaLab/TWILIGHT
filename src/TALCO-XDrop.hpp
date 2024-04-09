#include <stdint.h>
#include "timer.hpp"
#include <vector>
#include <string>

namespace Talco_xdrop {
    struct Params 
    {
        int16_t match;
        int16_t mismatch;
        int16_t gapOpen;
        int16_t gapExtend;

        int16_t xdrop;
        int16_t marker;

        Params(int16_t t_match, int16_t t_mismatch, int16_t t_gapOpen, int16_t t_gapExtend, int16_t t_xdrop, int16_t marker) : 
            match(t_match), mismatch(t_mismatch), gapOpen(t_gapOpen), gapExtend(t_gapExtend), xdrop(t_xdrop), marker(marker) {}
    };

    
    void Align (
        Params params,

        const std::vector<std::string>& reference,
        const std::vector<std::string>& query,
        std::vector<int8_t>& aln
        // size_t num_alignments
    );

    void Tile (
        // const std::vector<std::string>& reference,
        // const std::vector<std::string>& query,
        const std::vector<std::vector<int>>& reference, 
        const std::vector<std::vector<int>>& query, 
        Params params,
        int32_t &reference_idx,
        int32_t &query_idx,
        std::vector<int8_t> &aln,
        int8_t &state,
        bool &last_tile,
        const int &tile
    );

    int32_t Reduction_tree (
        const int32_t *C,
        const int32_t L,
        const int32_t U
    );
    void Traceback(
        const std::vector<int32_t> &ftr_length, 
        const std::vector<int32_t> &ftr_lower_limit, 
        const int32_t tb_start_addr, 
        const int16_t tb_start_ftr,
        const int8_t tb_state,
        const int16_t tb_start_idx,
        const int16_t ref_start_idx,
        const std::vector<int8_t> &tb,
        std::vector<int8_t> &aln
    );

    int Score (
        Params params,
        const std::vector<int8_t> &aln,
        const std::string &reference,
        const std::string &query,
        const int ref_idx,
        const int query_idx
    );

    
}


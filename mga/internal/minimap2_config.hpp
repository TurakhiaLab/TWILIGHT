#ifndef MINIMAP2_CONFIG_HPP
#define MINIMAP2_CONFIG_HPP

#include <string>
#include "option.hpp"

extern "C" {
#include "minimap.h"
}

/**
 * @brief Wrapper for minimap2 indexing and mapping configurations.
 */
class Minimap2Config {
public:
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    // Constructors
    explicit Minimap2Config(const std::string& preset = "asm5", bool needCigar = true) {
        applyPreset(preset);
        setCigar(needCigar);
    }

    explicit Minimap2Config(int asm_level, bool needCigar = true) {
        applyPreset("asm" + std::to_string(asm_level));
        setCigar(needCigar);
    }

    explicit Minimap2Config(const Option& opt, bool needCigar = true) {
        applyPreset(opt.mmAsmPreset);
        setCigar(needCigar);
    }

    // Apply minimap2 preset (e.g., "asm5", "asm10", "asm20", "sr", "map-ont")
    Minimap2Config& applyPreset(const std::string& preset = "asm5") {
        mm_set_opt(0, &iopt, &mopt);
        if (!preset.empty()) {
            mm_set_opt(preset.c_str(), &iopt, &mopt);
        }
        // Project defaults
        mopt.pri_ratio = 0.0f; // -p 0
        mopt.best_n = 10;      // default best_n
        return *this;
    }

    // Fluent Setter Methods
    Minimap2Config& setCigar(bool enable = true) {
        if (enable) {
            mopt.flag |= MM_F_CIGAR;
        } else {
            mopt.flag &= ~MM_F_CIGAR;
        }
        return *this;
    }

    Minimap2Config& setMaxGap(int max_gap) {
        mopt.max_gap = max_gap;
        return *this;
    }

    Minimap2Config& setBandwidth(int bw) {
        mopt.bw = bw;
        return *this;
    }

    Minimap2Config& setPriRatio(float pri_ratio) {
        mopt.pri_ratio = pri_ratio;
        return *this;
    }

    Minimap2Config& setBestN(int best_n) {
        mopt.best_n = best_n;
        return *this;
    }

    Minimap2Config& setMinChainScore(int min_chain_score) {
        mopt.min_chain_score = min_chain_score;
        return *this;
    }

    Minimap2Config& setKmerSize(int k) {
        iopt.k = k;
        return *this;
    }

    Minimap2Config& setWindowSize(int w) {
        iopt.w = w;
        return *this;
    }

    // Static Factory Methods for ASM presets
    static Minimap2Config ASM5(bool needCigar = true) {
        return Minimap2Config("asm5", needCigar);
    }

    static Minimap2Config ASM10(bool needCigar = true) {
        return Minimap2Config("asm10", needCigar);
    }

    static Minimap2Config ASM20(bool needCigar = true) {
        return Minimap2Config("asm20", needCigar);
    }

    // Purpose-based Static Factory Methods
    static Minimap2Config orientation(int asm_level = 5, int max_gap = 500) {
        return Minimap2Config(asm_level, false).setMaxGap(max_gap);
    }

    static Minimap2Config orientation(const Option& opt, int max_gap = 500) {
        return Minimap2Config(opt, false).setMaxGap(max_gap);
    }

    static Minimap2Config whole_genome(int asm_level = 5, int max_gap = 1000, int best_n = 10, int min_chain_score = 100, bool needCigar = true) {
        return Minimap2Config(asm_level, needCigar).setMaxGap(max_gap).setBestN(best_n).setMinChainScore(min_chain_score);
    }

    static Minimap2Config whole_genome(const Option& opt, int max_gap = 1000, int best_n = 10, int min_chain_score = 100, bool needCigar = true) {
        return Minimap2Config(opt, needCigar).setMaxGap(max_gap).setBestN(best_n).setMinChainScore(min_chain_score);
    }

    static Minimap2Config refine_block(int asm_level = 5, int max_gap = 2000, int best_n = 50, bool needCigar = true) {
        return Minimap2Config(asm_level, needCigar).setMaxGap(max_gap).setBestN(best_n);
    }

    static Minimap2Config refine_block(const Option& opt, int max_gap = 2000, int best_n = 50, bool needCigar = true) {
        return Minimap2Config(opt, needCigar).setMaxGap(max_gap).setBestN(best_n);
    }
};

/**
 * @brief Purpose-based Preset Classes
 */

// 1. Orientation configuration: CIGAR calculation is disabled for fast chaining & orientation detection.
class mm_config_orientation : public Minimap2Config {
public:
    explicit mm_config_orientation(int asm_level = 5, int max_gap = 500)
        : Minimap2Config(asm_level, false) {
        setMaxGap(max_gap);
    }

    explicit mm_config_orientation(const std::string& preset, int max_gap = 500)
        : Minimap2Config(preset, false) {
        setMaxGap(max_gap);
    }

    explicit mm_config_orientation(const Option& opt, int max_gap = 500)
        : Minimap2Config(opt, false) {
        setMaxGap(max_gap);
    }
};

// 2. Whole-genome alignment configuration: full alignment with CIGAR for main whole-genome comparisons.
class mm_config_whole_genome : public Minimap2Config {
public:
    explicit mm_config_whole_genome(int asm_level = 5, int max_gap = 1000, int best_n = 10, int min_chain_score = 100, bool needCigar = true)
        : Minimap2Config(asm_level, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
        setMinChainScore(min_chain_score);
    }

    explicit mm_config_whole_genome(const std::string& preset, int max_gap = 1000, int best_n = 10, int min_chain_score = 100, bool needCigar = true)
        : Minimap2Config(preset, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
        setMinChainScore(min_chain_score);
    }

    explicit mm_config_whole_genome(const Option& opt, int max_gap = 1000, int best_n = 10, int min_chain_score = 100, bool needCigar = true)
        : Minimap2Config(opt, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
        setMinChainScore(min_chain_score);
    }
};

// 3. Block refinement configuration: tuned for refining blocks and local alignments.
class mm_config_refine_block : public Minimap2Config {
public:
    explicit mm_config_refine_block(int asm_level = 5, int max_gap = 2000, int best_n = 50, bool needCigar = true)
        : Minimap2Config(asm_level, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }

    explicit mm_config_refine_block(const std::string& preset, int max_gap = 2000, int best_n = 50, bool needCigar = true)
        : Minimap2Config(preset, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }

    explicit mm_config_refine_block(const Option& opt, int max_gap = 2000, int best_n = 50, bool needCigar = true)
        : Minimap2Config(opt, needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

// Backward-compatibility aliases
using mm_config_1 = mm_config_whole_genome;
using mm_config_2 = mm_config_whole_genome;
using mm_config_3 = mm_config_refine_block;
using mm_config_4 = mm_config_whole_genome;
using mm_config_5 = mm_config_refine_block;

#endif // MINIMAP2_CONFIG_HPP

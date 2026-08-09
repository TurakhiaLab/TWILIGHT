#ifndef MINIMAP2_CONFIG_HPP
#define MINIMAP2_CONFIG_HPP

#include <string>

extern "C" {
#include "minimap.h"
}

class Minimap2Config {
public:
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    // 構造函式：初始化 mm_idxopt_t 與 mm_mapopt_t 並套用指定 preset
    explicit Minimap2Config(const std::string& preset = "asm5", bool needCigar = true) {
        applyPreset(preset);
        setCigar(needCigar);
    }

    // 套用 preset (如 "asm5", "asm10", "sr", "map-ont" 等)
    Minimap2Config& applyPreset(const std::string& preset = "asm5") {
        mm_set_opt(0, &iopt, &mopt);
        if (!preset.empty()) {
            mm_set_opt(preset.c_str(), &iopt, &mopt);
        }
        // 專案標準預設調整
        mopt.pri_ratio = 0.0f; // -p 0
        mopt.best_n = 50;     // -N 50
        return *this;
    }

    // 鏈式 Setter 方法 (Fluent API)
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

    Minimap2Config& setKmerSize(int k) {
        iopt.k = k;
        return *this;
    }

    Minimap2Config& setWindowSize(int w) {
        iopt.w = w;
        return *this;
    }

    // 靜態 Factory 方法
    static Minimap2Config ASM5(bool needCigar = true) {
        return Minimap2Config("asm5", needCigar);
    }

    static Minimap2Config ASM10(bool needCigar = true) {
        return Minimap2Config("asm10", needCigar);
    }

    static Minimap2Config ASM20(bool needCigar = true) {
        return Minimap2Config("asm20", needCigar);
    }

    // 自訂測試 Preset 靜態 Factory 方法 (預設基於 asm5，可調整 max_gap 與 best_n)
    static Minimap2Config mm_config_1(int max_gap = 1000, int best_n = 20, bool needCigar = true) {
        return Minimap2Config("asm5", needCigar).setMaxGap(max_gap).setBestN(best_n);
    }

    static Minimap2Config mm_config_2(int max_gap = 1000, int best_n = 50, bool needCigar = true) {
        return Minimap2Config("asm5", needCigar).setMaxGap(max_gap).setBestN(best_n);
    }

    static Minimap2Config mm_config_3(int max_gap = 2000, int best_n = 50, bool needCigar = true) {
        return Minimap2Config("asm5", needCigar).setMaxGap(max_gap).setBestN(best_n);
    }

    static Minimap2Config mm_config_4(int max_gap = 1000, int best_n = 20, bool needCigar = true) {
        return Minimap2Config("asm5", needCigar).setMaxGap(max_gap).setBestN(best_n);
    }

    static Minimap2Config mm_config_5(int max_gap = 2000, int best_n = 100, bool needCigar = true) {
        return Minimap2Config("asm5", needCigar).setMaxGap(max_gap).setBestN(best_n);
    }
};

// 預定義之標籤類別 (方便呼叫端直接作為物件傳入)
class Minimap2ConfigASM5 : public Minimap2Config {
public:
    explicit Minimap2ConfigASM5(bool needCigar = true) : Minimap2Config("asm5", needCigar) {}
};

class Minimap2ConfigASM10 : public Minimap2Config {
public:
    explicit Minimap2ConfigASM10(bool needCigar = true) : Minimap2Config("asm10", needCigar) {}
};

// 5 個自訂 Preset 類別 (方便作為物件傳入或自行調整參數)
class mm_config_1 : public Minimap2Config {
public:
    explicit mm_config_1(int max_gap = 500, int best_n = 50, bool needCigar = true)
        : Minimap2Config("asm5", needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

class mm_config_2 : public Minimap2Config {
public:
    explicit mm_config_2(int max_gap = 1000, int best_n = 50, bool needCigar = true)
        : Minimap2Config("asm5", needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

class mm_config_3 : public Minimap2Config {
public:
    explicit mm_config_3(int max_gap = 2000, int best_n = 50, bool needCigar = true)
        : Minimap2Config("asm5", needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

class mm_config_4 : public Minimap2Config {
public:
    explicit mm_config_4(int max_gap = 1000, int best_n = 20, bool needCigar = true)
        : Minimap2Config("asm5", needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

class mm_config_5 : public Minimap2Config {
public:
    explicit mm_config_5(int max_gap = 2000, int best_n = 100, bool needCigar = true)
        : Minimap2Config("asm5", needCigar) {
        setMaxGap(max_gap);
        setBestN(best_n);
    }
};

#endif // MINIMAP2_CONFIG_HPP

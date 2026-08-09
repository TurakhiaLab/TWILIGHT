#ifndef MINIMAP2_UTIL_HPP
#define MINIMAP2_UTIL_HPP

#include <string>
#include <vector>
#include "type.hpp"
#include "option.hpp"
#include "minimap2_config.hpp"

struct Alignment;
using Alignments = std::vector<Alignment>;

// 使用 Minimap2Config 的 runMinimap2 介面
Alignments runMinimap2(StringPairs& ref, StringPairs& qry, std::string refName, std::string qryName, Option& option, const Minimap2Config& config, bool write_paf = false);
Alignments runMinimap2(const SequenceRefs& ref, const SequenceRefs& qry, std::string refName, std::string qryName, Option& option, const Minimap2Config& config, bool write_paf = false);
Alignments runMinimap2(const std::string& refSeq, const std::string& qrySeq, const std::string& refID, const std::string& qryID, Option& option, const Minimap2Config& config, bool write_paf = false);

// 保持相容性的原 runMinimap2 介面 (內部自動帶入 Minimap2ConfigASM5(needCigar))
Alignments runMinimap2(StringPairs& ref, StringPairs& qry, std::string refName, std::string qryName, Option& option, bool needCigar = true, bool write_paf = false);
Alignments runMinimap2(const SequenceRefs& ref, const SequenceRefs& qry, std::string refName, std::string qryName, Option& option, bool needCigar = true, bool write_paf = false);
Alignments runMinimap2(const std::string& refSeq, const std::string& qrySeq, const std::string& refID, const std::string& qryID, Option& option, bool needCigar = true, bool write_paf = false);

namespace parser {
    Alignments parseMinimap2PAF(const std::string& filename);
}

#endif // MINIMAP2_UTIL_HPP

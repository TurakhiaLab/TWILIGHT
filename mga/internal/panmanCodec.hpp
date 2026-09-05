#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace mga {
namespace panman {

class PanmanCodec {
public:
    // IUPAC DNA 4-bit symbol mapping
    static uint8_t symbolToCode(char c);
    static char codeToSymbol(uint8_t code);

    // Pack a consensus DNA string into 32-bit packed words (8 symbols per word)
    static std::vector<uint32_t> packConsensusSequence(const std::string& gapFreeSeq);

    // Compute compacted payload bits for wire format
    static uint32_t compactMutationPayloadForWire(uint32_t packedMutation, uint8_t len);

    // Combine symbol mask and mutation type/len into 32-bit mutInfo for Cap'n Proto
    static uint32_t buildMutInfoWire(uint32_t symbolMask, uint8_t mutType, uint8_t len);
};

} // namespace panman
} // namespace mga

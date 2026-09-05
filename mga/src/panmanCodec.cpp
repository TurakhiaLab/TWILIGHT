#include "panmanCodec.hpp"
#include <cctype>

namespace mga {
namespace panman {

uint8_t PanmanCodec::symbolToCode(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 1;  // 0001
        case 'C': return 2;  // 0010
        case 'G': return 4;  // 0100
        case 'T':
        case 'U': return 8;  // 1000
        case 'R': return 5;  // A or G
        case 'Y': return 10; // C or T
        case 'S': return 6;  // C or G
        case 'W': return 9;  // A or T
        case 'K': return 12; // G or T
        case 'M': return 3;  // A or C
        case 'B': return 14; // C, G, T
        case 'D': return 13; // A, G, T
        case 'H': return 11; // A, C, T
        case 'V': return 7;  // A, C, G
        case 'N': return 15; // Any
        case '-': return 0;  // Gap / empty
        default: return 15;  // Default N
    }
}

char PanmanCodec::codeToSymbol(uint8_t code) {
    static const char symbols[] = {
        '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    };
    return symbols[code & 0x0F];
}

std::vector<uint32_t> PanmanCodec::packConsensusSequence(const std::string& gapFreeSeq) {
    std::vector<uint32_t> packed;
    size_t len = gapFreeSeq.length();
    size_t numWords = (len + 7) / 8;
    packed.reserve(numWords);

    for (size_t i = 0; i < numWords; i++) {
        uint32_t word = 0;
        for (size_t j = 0; j < 8; j++) {
            size_t idx = i * 8 + j;
            uint8_t code = (idx < len) ? symbolToCode(gapFreeSeq[idx]) : 0;
            uint8_t shift = static_cast<uint8_t>(4 * (7 - j));
            word |= (static_cast<uint32_t>(code) << shift);
        }
        packed.push_back(word);
    }
    return packed;
}

uint32_t PanmanCodec::compactMutationPayloadForWire(uint32_t packedMutation, uint8_t len) {
    uint8_t usedBits = len * 4; // 4-bit per symbol for DNA
    uint8_t storageBits = 24;   // 24 bits wire slot capacity
    if (usedBits >= storageBits) {
        return packedMutation;
    }
    return (packedMutation >> (storageBits - usedBits));
}

uint32_t PanmanCodec::buildMutInfoWire(uint32_t symbolMask, uint8_t mutType, uint8_t len) {
    uint32_t baseInfo = (static_cast<uint32_t>(len) << 4) | (mutType & 0x07);
    uint32_t wirePayload = compactMutationPayloadForWire(symbolMask, len);
    return (wirePayload << 8) | baseInfo;
}

} // namespace panman
} // namespace mga

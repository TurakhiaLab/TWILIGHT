#pragma once

#include <string>
#include <vector>
#include "block_set.hpp"
#include "phylogeny.hpp"

namespace mga {
namespace panman {

class PanmanWriter {
public:
    /**
     * Exports TWILIGHT BlockSet and Tree topology into a valid .panman file
     * (Cap'n Proto binary with LZMA compression) compatible with panmanUtils.
     *
     * @param rootBlockSet Pointer to the root BlockSet containing blocks & consensus.
     * @param tree Pointer to the TWILIGHT Tree structure (can be nullptr if newickString given).
     * @param newickString Newick representation of the phylogenetic tree.
     * @param outputPath Output path for the .panman file.
     * @return true on success, false on error.
     */
    static bool writePanMAN(BlockSet* rootBlockSet, 
                            Tree* tree, 
                            const std::string& newickString, 
                            const std::string& outputPath);
};

} // namespace panman
} // namespace mga

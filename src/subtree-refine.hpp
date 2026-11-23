#ifndef SUBTREE_REFINE_HPP
#define SUBTREE_REFINE_HPP

#include "msa.hpp"
#include "phylogeny.hpp"

namespace msa {
namespace accurate {

/**
 * Reconstruct a new guide subtree for a Layer 1 group using pairwise
 * local alignment scores / distances from LocalHomTable.
 * 
 * Supports UPGMA (default) and Neighbor-Joining (NJ).
 * Uses MAFFT score2dist normalization to prevent short-fragment high-identity artifacts.
 */
phylogeny::Tree* refineGroupSubtreeFromHomTable(
    HierarchyGroup& group,
    const LocalHomTable& localHomTable,
    phylogeny::Tree* oldSubT,
    SequenceDB* database,
    Option* option,
    const Params& param
);

/**
 * Perform iterative alignment refinement on a Layer 1 subtree.
 * 
 * Uses mini-batch rolling coordinate descent across candidate branch cuts.
 * Evaluates candidate bipartitions in batches of size teamThreads, applying
 * improvements immediately to cascade beneficial realignments into subsequent cuts.
 * Pre-allocates flat consistency memory pools to eliminate dynamic heap allocations.
 */
void iterativeRefineGroupSubtree(
    HierarchyGroup& group,
    phylogeny::Tree* subT,
    std::shared_ptr<LocalHomTable> localHomTable,
    SequenceDB* database,
    Option* option,
    const Params& param,
    int teamThreads = 8
);

} // namespace accurate
} // namespace msa

#endif // SUBTREE_REFINE_HPP


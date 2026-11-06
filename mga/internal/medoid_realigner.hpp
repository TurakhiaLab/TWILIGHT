#pragma once

#include "alignment.hpp"
#include "block_set.hpp"
#include "option.hpp"

// ============================================================================
// Medoid to Consensus Liftover and Flank-Guarded Realignment
// ============================================================================
// Transforms Minimap2 alignments generated on Medoid sequences into Consensus
// space alignments for downstream consumption by AlignmentCollection and BlockMerger.
//
// Invariants guaranteed:
//   1. Sum of Ref-consuming bases in CIGAR == (aln.refIdx.second - aln.refIdx.first)
//   2. Sum of Qry-consuming bases in CIGAR == std::abs(aln.qryIdx.second - aln.qryIdx.first)
//   3. Boundary-flanked Needleman-Wunsch realignment normalizes micro-indels.
//   4. Head/Tail boundary testing for extreme length asymmetry (e.g. 3bp vs 100bp)
//      prevents short fragments from shattering large contiguous indels.
// ============================================================================

Alignments realignMedoidAlignmentsToConsensus(
    const Alignments &medAlignments,
    BlockSet *refBlockSet,
    BlockSet *qryBlockSet,
    const Option &option);

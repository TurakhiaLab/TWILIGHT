#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string>
#include <vector>
#include <set>
#include <queue>

#include "type.hpp"
#include "option.hpp"


// =========================
// Alignment Utils
// =========================
CigarString adjustCigarWithVariations(CigarString& origCigar, Segment& refSeg, Segment& qrySeg, bool qryInverse, int qryConsLen);
CigarString extractSubCigar(const CigarString& origCigar, int refOffset, int refLen);
Alignments runMinimap2(StringPairs& ref, StringPairs& qry, std::string refName, std::string qryName, Option& option, bool needCigar = true, bool write_paf = false);
Alignments runMinimap2(const SequenceRefs& ref, const SequenceRefs& qry, std::string refName, std::string qryName, Option& option, bool needCigar = true, bool write_paf = false);
Alignments runMinimap2(const std::string& refSeq, const std::string& qrySeq, const std::string& refID, const std::string& qryID, Option& option, bool needCigar = true, bool write_paf = false);
Alignments splitSingleAlignment(const Alignment& aln,const std::set<int>& refCuts,const std::set<int>& qryCuts);
void snapAlignment(Alignment& aln, int r_pad_left, int r_pad_right, int q_pad_left, int q_pad_right);

// =========================
// End of Alignment Utils
// =========================


namespace parser {
    CigarString parseCigar(const std::string& cigar);
    Alignments parseMinimap2PAF(const std::string& filename);
}

struct Alignment {

    // ========= Variables ========= 
    AlignmentId ID;
    std::string refName; 
    std::string qryName;
    Range refIdx;
    Range qryIdx;
    CigarString CIGAR;
    bool inverse;
    bool primary;
        
    int alnLength;
    int alnScore;
    int mis;
    int ins;
    int del;

    bool valid;
    double energy;
    
    // ========= Function ========= 
    Alignment(): valid(true), energy(0) {};
    void setValid2False() { valid = false; }     
        
    int countVariationsInRange(BlockSet* blockSet, int aln_start, int aln_end);
    void updateEnergy(BlockSet* refBlockSet, BlockSet* qryBlockSet, double beta=10.0);
    void updateAlnLength();
    int mapCoordinate(int targetPos, bool targetIsRef);
    void finalizeAlignment(CoverageTracker& ref_coverageTracker, CoverageTracker& qry_coverageTracker, const std::map<int, std::string>& r_cut_reasons = {}, const std::map<int, std::string>& q_cut_reasons = {}, int snapLength=30);

};

struct CoordTracker {
    uint64_t blkId = 0;
    int localPos = -1;
};

struct CoverageTracker {
    struct SegmentInfo {
        int end;
        uint64_t blockId; // Which Block this region belongs to
    };

    // key: start_coordinate, value: {end_coordinate, alnId}
    std::map<int, SegmentInfo> intervals; 
    void add(int start, int end, uint64_t alnId);
    void overwrite(int start, int end, uint64_t newBlockId);
    int distLeft(int pos) const;
    int distRight(int pos) const;
    int getLeftOverlap(int start, int end) const;
    int getRightOverlap(int start, int end) const;
    std::set<uint64_t> getOverlappingIds(int qStart, int qEnd) const;
    void getCuts(int start, int end, std::set<int>& cuts) const;
    int getOverlapLength(int start, int end) const;
    bool isCovered(int start, int end) const;
    // void syncFromMap(const std::vector<CoordTracker>& globalMap);
    void syncFromMap(const CoordinateManager& coordMgr, bool isRefAxis);
};

struct CompareEnergy {
    bool operator()(const Alignment& a, const Alignment& b) {
        return a.energy > b.energy; // C++ priority_queue 預設是 Max-Heap，用 > 反轉為 Min-Heap
    }
};

struct AlignmentCollection {
    std::priority_queue<Alignment, std::vector<Alignment>, CompareEnergy> queue_;
    BlockSet* ref_blockset;
    BlockSet* qry_blockset;
    CoverageTracker ref_coverageTracker;
    CoverageTracker qry_coverageTracker;
    AlignmentCollection(Alignments& alignments, BlockSet* refBlockSet, BlockSet* qryBlockSet);
    // Alignment getBestAlignment(BlockSet* refBlockSet, BlockSet* qryBlockSet, int L_min=100);
    // Alignment getBestAlignment(BlockSet* refBlockSet, BlockSet* qryBlockSet, const BlockBoundaries& refBounds, const BlockBoundaries& qryBounds, int L_min=100);
    // Alignments getBestAlignments(BlockSet* refBlockSet, BlockSet* qryBlockSet, const BlockBoundaries& refBounds, const BlockBoundaries& qryBounds, int L_min=100);
    Alignments getBestAlignments(BlockSet* refBlockSet, BlockSet* qryBlockSet, CoordinateManager& coordMgr, int L_min=100);
    void addAlignment(Alignment aln);
};

struct SplitResult {
    Alignment head;   
    Alignment middle; 
    Alignment tail;
};

struct OverlapInfo {
    int primaryIdx;
    bool onRef; // true = ref, false = qry
    int overlapLen;
    int start;  // absolute coordinate of overlap start
    int end;    // absolute coordinate of overlap end
};

// Helper struct for cross-axis processing
struct CrossOverlap {
    int pIdx;         // The index of the 'Other' Primary (PA, PB)
    int sStart, sEnd; // The overlap range on S's CURRENT axis (The "Other" axis)
    int pStart, pEnd; // The overlap range mapped to P_Main's axis (The "Best" axis)
    bool onRef;       // The axis of pStart/pEnd (same as BestOverlap)
};

/*
    void identifyPrimaryAlignments(alnVec& alignments, BlockSet* ref_blockset, BlockSet* qry_blockset, stringMap& ref_seqs, stringMap& qry_seqs);
    void recoverSuperCoordinates(alnVec& alignments, BlockSet* ref_blockset, BlockSet* qry_blockset);

    void resolveGlobalOverlaps(alnVec& alignments);
    std::pair<Alignment, Alignment> singleSplit(const Alignment& parent, int splitPos, bool onRef);
    //vSplitResult splitAlignment(const Alignment& aln, int start, int end, bool onRef);
    // void detectDuplications(alnVec& alignments);
    void fillUnalignedRegions(alnVec& alignments, int refTotalLen, int qryTotalLen);
    bool validateCoverage(const alnVec& alignments, int refTotalLen, int qryTotalLen);
    void collectCutPoints(const std::vector<mga::Alignment>& alignments, std::set<int>& refCuts, std::set<int>& qryCuts);

    std::vector<Alignment> splitSingleAlignment(const Alignment& aln, const std::set<int>& refCuts, const std::set<int>& qryCuts);
    std::vector<Alignment> splitAlignmentsByCuts(const std::vector<Alignment>& alignments, const std::set<int>& refCuts, const std::set<int>& qryCuts);
*/


#endif
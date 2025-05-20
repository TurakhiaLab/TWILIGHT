#ifndef ALNUTIL_HPP
#define ALNUTIL_HPP

#include <fstream>
#include <complex>
#include "pocketfft_hdronly.h"

#ifndef TREE_HPP
#include "tree.hpp"
#endif

#ifndef SETTING_HPP
#include "setting.hpp"
#endif

#ifndef TALCO_HPP
#include "TALCO-XDrop.hpp"
#endif

#ifndef PARTITION_HPP
#include "treePartition.hpp"
#endif

class Seed
{
    public:
        int rs;
        int qs;
        int len;
        float score;
        Seed(int _rs, int _qs, int _len, float _score): rs(_rs), qs(_qs), len(_len), score(_score) {};
        int rend() {return this->rs + this->len;}
        int qend() {return this->qs + this->len;}
};

class Segment
{
    public:
        int rs;
        int re;
        int qs;
        int qe;
        Segment(int _rs, int _re, int _qs, int _qe): rs(_rs), re(_re), qs(_qs), qe(_qe) {};
        void update (int _rs, int _re, int _qs, int _qe) {
            this->rs = _rs; this->re = _re; this->qs = _qs; this->qe = _qe;
        }
};

bool compSeeds(const Seed* a, const Seed* b);

const int32_t PROFILE_LEN_TH = 1000;
const int32_t SEQNUM_TH = 100;

void updateNode(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util);
void calculateProfileFreq(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, char type, int32_t profileLen, int32_t profileSize, std::pair<int, int> startPos);
void calculatePSGOP(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t seqLen, std::pair<int32_t,int32_t> offset, std::pair<int32_t,int32_t>lens, Params& param);

void removeGappyColumns(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<int,int>& rawLens);
void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx);

void mergeInsertions (Tree* tree, Node* nodeRef, std::vector<Node*>& nodes, msa::utility* util, std::vector<std::vector<int8_t>>& alnBad);
void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln);
void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx);

void fallback2cpu(std::vector<int>& fallbackPairs, Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option);
void getPostOrderList(Node* node, std::stack<Node*>& postStack);
void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID);
double calColumnSimilarity(Tree* tree, Node* node, msa::utility* util, Params& param);

// void getReAlignRegion(std::vector<int8_t>& aln, std::vector<bool>& conserved, std::vector<reAlign*>& realign,  std::vector<std::pair<int,int>>& regions);
// void conservedColumns(std::vector<std::vector<float>>& freqs, std::vector<bool>& conserved, float weight);

void compute_vec (std::vector<std::vector<float>>& freq, std::vector<double>& polarity, std::vector<double>& volume, size_t fft_size);
void r2cFFT(const std::vector<double>& input, std::vector<std::complex<double>>& output);
void conjugate_and_multiply(std::vector<std::complex<double>>& fft1, std::vector<std::complex<double>>& fft2, std::vector<std::complex<double>>& result);
void c2rFFT(const std::vector<std::complex<double>>& input, std::vector<double>& output);
void getPeaks(const std::vector<double>& corr, std::vector<int>& peaks);
void getSeeds(std::vector<std::vector<float>>& freq1, std::vector<std::vector<float>>& freq2, std::vector<int>& peaks, std::vector<Seed*>& seeds);
void chainSeeds(const std::vector<Seed*>& seeds, std::vector<Seed*>& chainedSeeds);
void getAlnSegments(const std::vector<Seed*>& chainedSeeds, std::vector<Segment*>& segments, std::pair<int, int> lens);
#endif
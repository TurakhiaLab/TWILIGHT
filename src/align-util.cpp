#ifndef ALNUTIL_HPP
#include "align-util.hpp"
#endif

bool compSeeds(const Seed* a, const Seed* b) {
    return (a->rs < b->rs) || (a->rs == b->rs && a->qs < b->qs);
}

bool compSegments(const Segment* a, const Segment* b) {
    return (a->rs < b->rs) || (a->rs == b->rs && a->qs < b->qs);
}

constexpr std::array<float, 22> make_polarity_lookup() {
    std::array<float, 22> table{};
    table[0]  = (8.1  - 8.325) / 2.6905; // A (Ala)
    table[1]  = (5.5  - 8.325) / 2.6905; // C (Cys)
    table[2]  = (13.0 - 8.325) / 2.6905; // D (Asp)
    table[3]  = (12.3 - 8.325) / 2.6905; // E (Glu)
    table[4]  = (5.2  - 8.325) / 2.6905; // F (Phe)
    table[5]  = (9.0  - 8.325) / 2.6905; // G (Gly)    
    table[6]  = (10.4 - 8.325) / 2.6905; // H (Hys)
    table[7]  = (5.2  - 8.325) / 2.6905; // I (Ile)
    table[8]  = (11.3 - 8.325) / 2.6905; // K (Lys)
    table[9]  = (4.9  - 8.325) / 2.6905; // L (Leu)
    table[10] = (5.7  - 8.325) / 2.6905; // M (Met)
    table[11] = (11.6 - 8.325) / 2.6905; // N (Asn)
    table[12] = (8.0  - 8.325) / 2.6905; // P (Pro)
    table[13] = (10.5 - 8.325) / 2.6905; // Q (Gln)
    table[14] = (10.5 - 8.325) / 2.6905; // R (Arg)
    table[15] = (9.2  - 8.325) / 2.6905; // S (Ser)
    table[16] = (8.6  - 8.325) / 2.6905; // T (Thr)
    table[17] = (5.9  - 8.325) / 2.6905; // V (Val) 
    table[18] = (5.4  - 8.325) / 2.6905; // W (Trp)
    table[19] = (6.2  - 8.325) / 2.6905; // Y (Tyr)
    table[20] = 0.0;
    table[21] = 0.0;
    return table;
}

constexpr std::array<float, 22> make_volume_lookup() {
    std::array<float, 22> table{};
    table[0]  = (31   - 84.025) / 42.953; // A (Ala)
    table[1]  = (55   - 84.025) / 42.953; // C (Cys)
    table[2]  = (54   - 84.025) / 42.953; // D (Asp)
    table[3]  = (83   - 84.025) / 42.953; // E (Glu)
    table[4]  = (132  - 84.025) / 42.953; // F (Phe)
    table[5]  = (3    - 84.025) / 42.953; // G (Gly)    
    table[6]  = (96   - 84.025) / 42.953; // H (Hys)
    table[7]  = (111  - 84.025) / 42.953; // I (Ile)
    table[8]  = (119  - 84.025) / 42.953; // K (Lys)
    table[9]  = (111  - 84.025) / 42.953; // L (Leu)
    table[10] = (105  - 84.025) / 42.953; // M (Met)
    table[11] = (56   - 84.025) / 42.953; // N (Asn)
    table[12] = (32.5 - 84.025) / 42.953; // P (Pro)
    table[13] = (85   - 84.025) / 42.953; // Q (Gln)
    table[14] = (124  - 84.025) / 42.953; // R (Arg)
    table[15] = (32   - 84.025) / 42.953; // S (Ser)
    table[16] = (61   - 84.025) / 42.953; // T (Thr)
    table[17] = (84   - 84.025) / 42.953; // V (Val) 
    table[18] = (170  - 84.025) / 42.953; // W (Trp)
    table[19] = (136  - 84.025) / 42.953; // Y (Tyr)
    table[20] = 0.0;
    table[21] = 0.0;
    return table;
}

constexpr auto polarity_lookup = make_polarity_lookup();
constexpr auto volume_lookup = make_volume_lookup();

void compute_vec (std::vector<std::vector<float>>& freq, std::vector<double>& polarity, std::vector<double>& volume, size_t fft_size) {
    polarity = std::vector<double>(fft_size, 0.0);
    volume = std::vector<double>(fft_size, 0.0);
    float weight = 0.0;
    for (auto w: freq[0]) weight += w;
    for (int i = 0; i < freq.size(); ++i) {
        for (int j = 0; j < 20; ++j) {
            polarity[i] += polarity_lookup[j] * freq[i][j] / weight;
            volume[i] += volume_lookup[j] * freq[i][j] / weight;
        }
    }
    return;
}

void r2cFFT(const std::vector<double>& input, std::vector<std::complex<double>>& output) {
    size_t N = input.size();
    output = std::vector<std::complex<double>>(N/2+1);
    pocketfft::detail::shape_t shape = {N};
    pocketfft::detail::stride_t stride_in = {sizeof(double)};
    pocketfft::detail::stride_t stride_out = {sizeof(std::complex<double>)};
    pocketfft::detail::shape_t axes = {0};
    pocketfft::detail::r2c<double>(
        shape, 
        stride_in, 
        stride_out, 
        axes, 
        true,          
        input.data(), 
        output.data(),
        1.0
    );
    // for (int i =  N/2+1; i < fft_size; ++i) output.push_back({0,0});
    // for (auto i: output) std::cout << i << ',';
    // std::cout << '\n';
    return;
}

void conjugate_and_multiply(std::vector<std::complex<double>>& fft1, std::vector<std::complex<double>>& fft2, std::vector<std::complex<double>>& result) {
    if (fft1.size() < fft2.size()) for (int i = 0; i < (fft2.size()-fft1.size()); ++i) fft1.push_back({0,0});
    if (fft1.size() > fft2.size()) for (int i = 0; i < (fft1.size()-fft2.size()); ++i) fft2.push_back({0,0});
    result = std::vector<std::complex<double>>(fft1.size());
    for (size_t i = 0; i < fft1.size(); ++i) {
        result[i] = fft1[i] * std::conj(fft2[i]);  // cross-correlation
    }
    return;
}

void c2rFFT(const std::vector<std::complex<double>>& input, std::vector<double>& output) {
    size_t N = (input.size()-1) * 2;
    output = std::vector<double>(N);
    pocketfft::detail::shape_t shape = {N};
    pocketfft::detail::stride_t stride_in = {sizeof(double)};
    pocketfft::detail::stride_t stride_out = {sizeof(std::complex<double>)};
    pocketfft::detail::shape_t axes = {0};
    pocketfft::detail::c2r<double>(
        shape,
        stride_out,
        stride_in,
        axes,
        false,
        input.data(),
        output.data(),
        1.0f / static_cast<float>(N)
    );
    return;
}

void getPeaks(const std::vector<double>& corr, std::vector<int>& peaks) {
    const size_t K = 5;
    double max_peak = 0.0;
    const float min_height = 0.5; // only store peaks that height > max_peak * min_height
    auto comp = [](const std::pair<double, size_t>& a, const std::pair<double, size_t>& b) {
        return a.first > b.first;
    };
    std::priority_queue<std::pair<double, size_t>, std::vector<std::pair<double, size_t>>, decltype(comp)> min_heap(comp);
    for (size_t i = 0; i < corr.size(); ++i) {
        if (corr[i] > max_peak) max_peak = corr[i];
        if (min_heap.size() < K) {
            min_heap.emplace(corr[i], i);
        } else if (corr[i] > min_heap.top().first) {
            min_heap.pop();
            min_heap.emplace(corr[i], i);
        }
    }
    while (!min_heap.empty()) {
        // std::cout << min_heap.top().second << ':' << min_heap.top().first << '\n';
        if (min_heap.top().first >= min_height * max_peak) peaks.push_back(min_heap.top().second);
        min_heap.pop();
    }
    return;
}

void getSeeds(std::vector<std::vector<float>>& freq1, std::vector<std::vector<float>>& freq2, std::vector<int>& peaks, std::vector<Seed*>& seeds) {
    // float entropy = 0.0;
    float max_entropy = 0.7; // larger than this means not very conserved
    int count = 0;
    int L1 = freq1.size(), L2 = freq2.size();
    int maxL = std::max(L1, L2), minL = std::min(L1, L2);
    std::vector<float> entropy;
    const int window = 30;
    const float avg_entropy = 1.0;
    // Amino Acids             |  Group
    // A, F, I, L, M, V, Y, W  |   0
    // N, Q, S, T              |   1 
    // C, G, P                 |   2
    // H, K, R                 |   3
    // D, E                    |   4
    //                                  A      C      D      E      F      G      H      I      K     L       M
    std::map<int, int> amino_groups = {{0,0}, {1,2}, {2,4}, {3,4}, {4,0}, {5,2}, {6,3}, {7,0}, {8,3}, {9,0}, {10,0}, 
    //                                  N       P       Q       R       S       T       V       W       Y
                                       {11,1}, {12,2}, {13,1}, {14,3}, {15,1}, {16,1}, {17,0}, {18,0}, {19,0}};
    float group_ratio [5] = {0.0};
    float weight = 0.0;
    for (auto f: freq1[0]) weight += f;
    for (auto f: freq2[0]) weight += f;
    for (auto it = peaks.rbegin(); it != peaks.rend(); ++it) {
        auto p = *it;
        entropy.clear();
        int shift = (p <= maxL / 2) ? p : -(maxL - p);
        int st1 = std::max(0, -shift);
        int st2 = std::max(0, shift);
        int validL1 = L1 - st1;
        int validL2 = L2 - st2;
        for (int k = 0; k < std::min(validL1, validL2); ++k) {
            for (int i = 0; i < 5; ++i) group_ratio[i] = 0;
            float entropy_temp = 0.0;
            for (int i = 0; i < 20; ++i) group_ratio[amino_groups[i]] += (freq1[st1+k][i] / (weight)) + (freq2[st2+k][i] / weight);
            for (int i = 0; i < 5; ++i) if (group_ratio[i] > 0) entropy_temp -= (group_ratio[i] * std::log2(group_ratio[i]));
            entropy.push_back(entropy_temp);
            // std::cout << k << ": " << entropy_temp << '\n';
        }
        // Sliding window
        int start = 0;
        // std::cout << "Shift: " << shift << ", Seeds:\n";
        while (start < entropy.size()) {
            float sum = 0.0, sum_window = 0.0;
            int end = start;
            float avg = 0.0, avg_window = 0.0;
            while (end < entropy.size()) {
                int len = end - start + 1;
                sum += entropy[end];
                if (len <= window) sum_window += entropy[end];
                else sum_window += (entropy[end] - entropy[end-window]);
                avg = sum / len;
                avg_window = sum_window / window;
                if (avg_window < avg_entropy) ++end;
                else break;
            }
            int len = end - start;
            if (len >= window) {
                Seed* newSeed = new Seed(st1+start, st2+start, len, avg); // use avg_entropy-avg, larger is better for chaining
                // printf("(%d,%d),%d, %f, %f\n", newSeed->rs, newSeed->qs, newSeed->len, newSeed->score, avg_window);
                seeds.push_back(newSeed); // end - 1 is last valid index
                start = end; // move start beyond current seed
            } else {
                start += 1; // shift start by 1 if no valid seed
            }
        }
    }
    sort(seeds.begin(), seeds.end(), compSeeds);
    return;
}

void chainSeeds(const std::vector<Seed*>& seeds, std::vector<Seed*>& chainedSeeds) {
    
    int n = seeds.size();
    if (n == 0) return;
    std::vector<float> DP (n);        // Best score ending at i
    std::vector<int> prev (n, -1);     // To trace back the path

    float avgLen = 0.0;
    for (int i = 0; i < n; ++i) avgLen += seeds[i]->len;
    avgLen /= n;

    for (int i = 0; i < n; ++i) DP[i] = static_cast<float>(seeds[i]->len);

    // TODO: optimized this part
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            if (seeds[j]->rs < seeds[i]->rs && seeds[j]->qs < seeds[i]->qs) {
                float newScore = DP[j] + std::min(std::min(seeds[i]->rs-seeds[j]->rs, seeds[i]->qs-seeds[j]->qs), seeds[i]->len);
                int l = abs((seeds[i]->rs-seeds[j]->rs) - (seeds[i]->qs-seeds[j]->qs));
                if (l > 0) newScore -= (0.01 * avgLen * l + 0.5 * std::log2(l));
                if (newScore > DP[i]) {
                    DP[i] = newScore;
                    prev[i] = j;
                }
            }
        }
    }

    // Find lowest entropy
    int bestIdx = std::max_element(DP.begin(), DP.end()) - DP.begin();
    // Traceback
    chainedSeeds.clear();
    int nowIdx = bestIdx;
    while (true) {
        chainedSeeds.push_back(seeds[nowIdx]);
        nowIdx = prev[nowIdx];
        if (nowIdx == -1) break;
    }
    std::reverse(chainedSeeds.begin(), chainedSeeds.end());
    std::cout << "Chained Seeds: \n";
    for (auto s: chainedSeeds) printf("(%d,%d),%d,%f\n", s->rs, s->qs, s->len, s->score);
    return;
}


void getAlnSegments(const std::vector<Seed*>& chainedSeeds, std::vector<Segment*>& segments, std::pair<int, int> lens) {
    if (chainedSeeds.empty()) {
        Segment* current = new Segment(0, lens.first, 0, lens.second);
        segments.push_back(current);
        return;
    }
    const int max_dist = 10;
    const int overlapped = 10;
        
    if (chainedSeeds.size() == 1) {
        int new_rs = std::max(chainedSeeds[0]->rs - max_dist, 0);
        int new_qs = std::max(chainedSeeds[0]->qs - max_dist, 0);
        int new_re = std::min(chainedSeeds[0]->rend() + max_dist, lens.first);
        int new_qe = std::min(chainedSeeds[0]->qend() + max_dist, lens.second);
        Segment* current = new Segment(new_rs, new_re, new_qs, new_qe);
        segments.push_back(current);
    }
    else {
        Segment* current = new Segment(chainedSeeds[0]->rs, chainedSeeds[0]->rend(), chainedSeeds[0]->qs, chainedSeeds[0]->qend());
        for (size_t i = 1; i < chainedSeeds.size(); ++i) {
            Seed* next = chainedSeeds[i];
            bool overlap1 = current->re >= next->rs - max_dist;
            bool overlap2 = current->qe >= next->qs - max_dist;
            if (overlap1 && overlap2) {
                // Merge: take min start, max end, average score
                int new_rs = std::min(current->rs, next->rs);
                int new_qs = std::min(current->qs, next->qs);
                int new_re = std::max(current->re, next->rend());
                int new_qe = std::max(current->qe, next->qend());
                current->update(new_rs, new_re, new_qs, new_qe);
                
            } 
            if (!(overlap1 && overlap2) || i == chainedSeeds.size()-1) {
                int new_rs = std::max(current->rs - max_dist, 0);
                int new_qs = std::max(current->qs - max_dist, 0);
                int new_re = std::min(current->re + max_dist, lens.first);
                int new_qe = std::min(current->qe + max_dist, lens.second);
                current->update(new_rs, new_re, new_qs, new_qe);
                segments.push_back(current);
                current = new Segment(next->rs, next->rend(), next->qs, next->qend());
            }
        }
    }
    std::vector<Segment*> temp;
    for (int i = 0; i < segments.size(); ++i) {
        auto s = segments[i];
        if (i == 0) {
            if (s->rs > 0 || s->qs > 0) {
                Segment* head = new Segment(0, s->rs + overlapped, 0, s->qs + overlapped);
                temp.push_back(head);
            }
        }
        else {
            auto p = segments[i-1];
            Segment* middle = new Segment(std::max(0, p->re - overlapped), s->rs + overlapped, std::max(0, p->qe - overlapped), s->qs + overlapped);
            temp.push_back(middle);
        }
        if (i == segments.size() - 1) {
            if (s->re < lens.first || s->qe < lens.second) {
                Segment* tail = new Segment(std::max(0, s->re - overlapped), lens.first, std::max(0, s->qe - overlapped), lens.second);
                temp.push_back(tail);
            }
        }
    }
    for (auto t: temp) segments.push_back(t);
    std::sort(segments.begin(), segments.end(), compSegments);
    std::cout << "Alignment Segments: \n";
    for (auto s: segments) printf("(%d,%d),(%d,%d)\n", s->rs, s->re, s->qs, s->qe);
    return;    
}


void updateNode(Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util) {
    for (auto n_pair: nodes) {
        auto n = std::make_pair(tree->allNodes[n_pair.first->identifier], tree->allNodes[n_pair.second->identifier]);
        if (n.first->is_leaf() && n.first->msaIdx.empty()) {
            tree->allNodes[n.first->identifier]->msaIdx.push_back(util->seqsIdx[n.first->identifier]);
        }
        else if (tree->allNodes[n.first->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.first->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.second->identifier)) {
                    // if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.first->identifier]->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    tree->allNodes[n.first->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.first->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
        if (n.second->is_leaf() && n.second->msaIdx.empty()) {
            tree->allNodes[n.second->identifier]->msaIdx.push_back(util->seqsIdx[n.second->identifier]);
        }
        else if (tree->allNodes[n.second->identifier]->msaIdx.size() == 0) {
            Node* node = tree->allNodes[n.second->identifier];
            int grpID = node->grpID;
            for (int childIndex=0; childIndex<node->children.size(); childIndex++) {
                if ((node->children[childIndex]->grpID == -1 || node->children[childIndex]->grpID == grpID) && (node->children[childIndex]->identifier != n.first->identifier)) {
                    // if (node->children[childIndex]->is_leaf()) tree->allNodes[node->children[childIndex]->identifier]->msaIdx.push_back(util->seqsIdx[node->children[childIndex]->identifier]);
                    tree->allNodes[n.second->identifier]->msaFreq = node->children[childIndex]->msaFreq;
                    node->children[childIndex]->msaFreq.clear();
                    tree->allNodes[n.second->identifier]->msaIdx = node->children[childIndex]->msaIdx;
                    util->seqsLen[n.second->identifier] = util->seqsLen[node->children[childIndex]->identifier];
                    break;
                }
            }
        }
    }
}

void fallback2cpu(std::vector<int>& fallbackPairs, Tree* tree, std::vector<std::pair<Node*, Node*>>& nodes, msa::utility* util, msa::option* option) {
    int totalSeqs = 0;
    std::sort(fallbackPairs.begin(), fallbackPairs.end());
    for (int i = 0; i < fallbackPairs.size(); ++i) {
        int nIdx = fallbackPairs[i];
        int grpID = tree->allNodes[nodes[nIdx].first->identifier]->grpID;
        int32_t refNum = tree->allNodes[nodes[nIdx].first->identifier]->msaIdx.size();
        int32_t qryNum = tree->allNodes[nodes[nIdx].second->identifier]->msaIdx.size();
        // if (refNum < qryNum) {
        if ((refNum < qryNum) || (refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]])) {
            bool fallback = !(refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]]) || option->noFilter;
            if (fallback) {
                if (util->badSequences.find(grpID) == util->badSequences.end()) {
                    util->badSequences[grpID] = std::vector<std::string> (0);
                }
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                if (refNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]]) {
                    util->lowQuality[tree->allNodes[nodes[nIdx].first->identifier]->msaIdx[0]] = false;
                }
            }
            
            int32_t refLen = util->seqsLen[nodes[nIdx].first->identifier];
            int32_t qryLen = util->seqsLen[nodes[nIdx].second->identifier];
            util->seqsLen[nodes[nIdx].second->identifier] = refLen;
            util->seqsLen[nodes[nIdx].first->identifier] = qryLen;
            auto temp = nodes[nIdx].second->msaIdx;
            tree->allNodes[nodes[nIdx].second->identifier]->msaIdx = nodes[nIdx].first->msaIdx;
            tree->allNodes[nodes[nIdx].first->identifier]->msaIdx = temp;
            auto tempFreq = nodes[nIdx].second->msaFreq;
            tree->allNodes[nodes[nIdx].second->identifier]->msaFreq = nodes[nIdx].first->msaFreq;
            tree->allNodes[nodes[nIdx].first->identifier]->msaFreq = tempFreq;
            totalSeqs += refNum;
        }
        else {
            // std::cout << "Deferring " << nodes[nIdx].second->identifier << '\n';
            bool fallback = !(qryNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]]) || option->noFilter;
            if (fallback) {
                if (util->badSequences.find(grpID) == util->badSequences.end()) {
                    util->badSequences[grpID] = std::vector<std::string> (0);
                }
                util->badSequences[grpID].push_back(nodes[nIdx].second->identifier);
                if (qryNum == 1 && util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]]) {
                    util->lowQuality[tree->allNodes[nodes[nIdx].second->identifier]->msaIdx[0]] = false;
                }
            }
            totalSeqs += qryNum;
        }
    }
    if (option->printDetail) {
        if (fallbackPairs.size() == 1 && totalSeqs == 1) std::cout << "Deferring/excluding 1 pair (1 sequence).\n";
        else if (fallbackPairs.size() == 1 && totalSeqs > 1) printf("Deferring/excluding 1 pair (%d sequences).\n", totalSeqs); 
        else printf("Deferring/excluding %lu pair (%d sequences).\n", fallbackPairs.size(), totalSeqs); 
    }
    return;
}

void updateAlignment(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln) {
    if (util->nowProcess == 0) {
        if (aln.size() != util->seqsLen[nodes.first->identifier]) {
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.first->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
            for (int idx = range.begin(); idx < range.end(); ++idx) {
                int sIdx = tree->allNodes[nodes.first->identifier]->msaIdx[idx];  
                util->memCheck(aln.size(), sIdx);
            // for (auto sIdx: tree->allNodes[nodes[nIdx].first->identifier]->msaIdx) {
                int orgIdx = 0;
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - storeFrom;
                for (int k = 0; k < aln.size(); ++k) {
                    if (aln[k] == 0 || aln[k] == 2) {
                        util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][k] = '-';
                    }
                }
                util->changeStorage(sIdx);
            }
            });
            });
        }
        if (aln.size() != util->seqsLen[nodes.second->identifier]) {
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodes.second->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
            for (int idx = range.begin(); idx < range.end(); ++idx) {
                int sIdx = tree->allNodes[nodes.second->identifier]->msaIdx[idx];
                util->memCheck(aln.size(), sIdx);
            // for (auto sIdx: tree->allNodes[nodes[nIdx].second->identifier]->msaIdx) {
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - util->seqsStorage[sIdx];
                int orgIdx = 0;
                for (int k = 0; k < aln.size(); ++k) {
                    // if ((hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 0 || (hostAln[gn][n*2*seqLen+j] & 0xFFFF) == 1) {
                    if ((aln[k] & 0xFFFF) == 0 || (aln[k] & 0xFFFF) == 1) {
                        util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else {
                        util->alnStorage[storeTo][sIdx][k] = '-';
                    }
                }
                // util->seqsLen[nodes[nIdx].second->identifier] = hostAlnLen[gn][n];
                util->changeStorage(sIdx);
            }
            });
            });
        }
        util->seqsLen[nodes.first->identifier] = aln.size();
        util->seqsLen[nodes.second->identifier] = aln.size();  
        for (auto q: tree->allNodes[nodes.second->identifier]->msaIdx) {
            tree->allNodes[nodes.first->identifier]->msaIdx.push_back(q);
        }
        tree->allNodes[nodes.second->identifier]->msaIdx.clear();
    }
    else if (util->nowProcess == 2) {
        // int seqLenR = util->seqsLen[nodes.first->identifier];
        // int seqLenQ = util->seqsLen[nodes.second->identifier];
        if (tree->allNodes[nodes.first->identifier]->msa.empty())  tree->allNodes[nodes.first->identifier]->msa.push_back(nodes.first->identifier);
        if (tree->allNodes[nodes.second->identifier]->msa.empty()) tree->allNodes[nodes.second->identifier]->msa.push_back(nodes.second->identifier);
        for (auto seqName: tree->allNodes[nodes.first->identifier]->msa) {
            std::vector<int8_t> seqAln;
            int orgIdx = 0;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 2) {
                    seqAln.push_back(tree->allNodes[seqName]->msaAln[orgIdx]);
                    orgIdx++;
                }
                else {
                    seqAln.push_back(3);
                }
            }
            assert(orgIdx == tree->allNodes[seqName]->msaAln.size());
            tree->allNodes[seqName]->msaAln = seqAln;
            util->seqsLen[seqName] = aln.size();
        }
        for (auto seqName: tree->allNodes[nodes.second->identifier]->msa) {
            std::vector<int8_t> seqAln;
            int orgIdx = 0;
            for (int k = 0; k < aln.size(); ++k) {
                if (aln[k] == 0 || aln[k] == 1) {
                    seqAln.push_back(tree->allNodes[seqName]->msaAln[orgIdx]);
                    orgIdx++;
                }
                else {
                    seqAln.push_back(3);
                }
            }
            assert(orgIdx == tree->allNodes[seqName]->msaAln.size());
            tree->allNodes[seqName]->msaAln = seqAln;
            util->seqsLen[seqName] = aln.size();
        }
        for (auto q: tree->allNodes[nodes.second->identifier]->msa) {
            tree->allNodes[nodes.first->identifier]->msa.push_back(q);
        }
        for (auto q: tree->allNodes[nodes.first->identifier]->msa) {
            if (q != nodes.first->identifier) tree->allNodes[q]->msa = tree->allNodes[nodes.first->identifier]->msa;
        }
    }
    return;
}

void updateFrequency(Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, std::vector<int8_t>& aln, float refWeight, float qryWeight, std::pair<int, int>& debugIdx) {
    int rIdx = 0, qIdx = 0;
    if (tree->allNodes[nodes.first->identifier]->msaFreq.empty() || tree->allNodes[nodes.second->identifier]->msaFreq.empty()) return;
    int m_size = tree->allNodes[nodes.first->identifier]->msaFreq[0].size();
    std::vector<std::vector<float>> mergeFreq(aln.size(), std::vector<float>(m_size,0.0));
    if (util->nowProcess < 2) {
        for (int j = 0; j < aln.size(); ++j) {
            if (aln[j] == 0) {
                for (int k = 0; k < m_size; ++k) mergeFreq[j][k] = 
                (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k] + 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k]);
                ++rIdx; ++qIdx;
            }
            else if (aln[j] == 1) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k];
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][m_size-1] + 1.0 * refWeight);
                ++qIdx;
            }
            else if (aln[j] == 2) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][m_size-1] + 1.0 * qryWeight);
                ++rIdx;
            }
        }  
        debugIdx.first = rIdx; debugIdx.second = qIdx;
        tree->allNodes[nodes.first->identifier]->msaFreq.clear();
        tree->allNodes[nodes.second->identifier]->msaFreq.clear();
        tree->allNodes[nodes.first->identifier]->msaFreq = mergeFreq;
    }
    else {
        float refNum_f = 0, qryNum_f = 0;
        for (int t = 0; t < m_size; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < m_size; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        int32_t refNum = static_cast<int32_t>(round(refNum_f));
        int32_t qryNum = static_cast<int32_t>(round(qryNum_f));
        for (int j = 0; j < aln.size(); ++j) {
            if (aln[j] == 0) {
                for (int k = 0; k < m_size; ++k) mergeFreq[j][k] = 
                (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k] + 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k]);
                ++rIdx; ++qIdx;
            }
            else if (aln[j] == 1) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][k];
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.second->identifier]->msaFreq[qIdx][m_size-1] + 1.0 * refNum);
                ++qIdx;
            }
            else if (aln[j] == 2) {
                for (int k = 0; k < m_size-1; ++k) mergeFreq[j][k] = 
                tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][k]; 
                mergeFreq[j][m_size-1] = (tree->allNodes[nodes.first->identifier]->msaFreq[rIdx][m_size-1] + 1.0 * qryNum);
                ++rIdx;
            }
        }  
        debugIdx.first = rIdx; debugIdx.second = qIdx;
        for (auto node: tree->allNodes[nodes.first->identifier]->msa) {
            tree->allNodes[node]->msaFreq = mergeFreq;
            assert(tree->allNodes[node]->msaFreq.size() == tree->allNodes[node]->msaAln.size());
        }
    }
    return;
}

void getMsaHierachy(std::vector<std::pair<std::pair<Node*, Node*>, int>>& alnOrder, std::stack<Node*> postOrder, int grpID) {
    
    std::map<std::string, int> NodeAlnOrder;

    while(!postOrder.empty()) {
        Node* node = postOrder.top();
        // Leaf or not belongs to the subtree
        if (!(node->grpID==-1 || node->grpID==grpID) || node->is_leaf()) {
            postOrder.pop();
            continue;
        }
        // Collect children that belong to the subtree
        std::vector<Node*> children;
        for (int i = 0; i < node->children.size(); ++i) {
            if (node->children[i]->grpID == grpID) children.push_back(node->children[i]);
        }
        // Useless node, remove from the subtree
        if (children.empty()) {
            node->grpID = -2;
            postOrder.pop();
            for (auto it = node->parent->children.begin(); it != node->parent->children.end(); ++it) {
                if ((*it)->identifier == node->identifier) {
                    node->parent->children.erase(it);
                    break;
                }
            }
            continue;
        }
        // Only one child, merge the child and the parent 
        if (children.size() == 1 && node->parent != nullptr) {
            if (node->parent->grpID == grpID) {
                for (int chIdx = 0; chIdx < node->parent->children.size(); ++chIdx) {
                    if (node->parent->children[chIdx]->identifier == node->identifier) {
                        node->parent->children[chIdx] = children[0];
                        children[0]->branchLength += node->branchLength;
                        children[0]->parent = node->parent;
                        break;
                    }
                }
                postOrder.pop();
                continue;
            }
        }
        // Pair children
        while (children.size() > 1) {
            std::vector<Node*> nodeLeft;
            for (int i = 0; i < children.size()-1; i+=2) {
                int firstIdx  = (NodeAlnOrder.find(children[i]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i]->identifier]+1 : 0;
                int secondIdx = (NodeAlnOrder.find(children[i+1]->identifier) != NodeAlnOrder.end()) ? NodeAlnOrder[children[i+1]->identifier]+1 : 0;
                int maxIdx = std::max(firstIdx, secondIdx);
                NodeAlnOrder[children[i]->identifier] = maxIdx;
                NodeAlnOrder[children[i+1]->identifier] = maxIdx;
                alnOrder.push_back(std::make_pair(std::make_pair(children[i], children[i+1]), maxIdx));
                nodeLeft.push_back(children[i]);
            }
            if (children.size()%2 == 1) nodeLeft.push_back(children.back());
            children = nodeLeft;
        }
        NodeAlnOrder[node->identifier] = NodeAlnOrder[children[0]->identifier];
        postOrder.pop();
    }
}

void getPostOrderList(Node* node, std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(node); 
    Node* current; 
  
    while (!s1.empty()) {
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == current->grpID) s1.push(current->children[i]);     
        }
    } 
    return;
}

double calColumnSimilarity(Tree* tree, Node* node, msa::utility* util, Params& param) {
    int len = util->seqsLen[node->identifier];
    std::vector<std::vector<float>> freq;
    int seqNum = tree->allNodes[node->identifier]->msaIdx.size();
    float totalWeight = 0.0;
    for (auto sIdx: tree->allNodes[node->identifier]->msaIdx)  totalWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
    if (tree->allNodes[node->identifier]->msaFreq.empty()) {
        freq = std::vector<std::vector<float>> (len, std::vector<float>(6,0)); 
        for (auto sIdx: tree->allNodes[node->identifier]->msaIdx) { 
            int storage = util->seqsStorage[sIdx];
            std::string name = util->seqsName[sIdx];
            float w = tree->allNodes[name]->weight / totalWeight * seqNum;
            tbb::this_task_arena::isolate( [&]{
            tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
            for (int s = r.begin(); s < r.end(); ++s) {
                if      (util->alnStorage[storage][sIdx][s] == 'A' || util->alnStorage[storage][sIdx][s] == 'a') freq[s][0]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'C' || util->alnStorage[storage][sIdx][s] == 'c') freq[s][1]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'G' || util->alnStorage[storage][sIdx][s] == 'g') freq[s][2]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == 'T' || util->alnStorage[storage][sIdx][s] == 't' ||
                         util->alnStorage[storage][sIdx][s] == 'U' || util->alnStorage[storage][sIdx][s] == 'u') freq[s][3]+=1.0*w;
                else if (util->alnStorage[storage][sIdx][s] == '-')                                              freq[s][5]+=1.0*w;
                else                                                                                             freq[s][4]+=1.0*w;
            }
            });
            });
        }
    }
    else {
        freq = tree->allNodes[node->identifier]->msaFreq;
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
        for (int s = r.begin(); s < r.end(); ++s) {
            for (int t = 0; t < 6; ++t) freq[s][t] = freq[s][t] / totalWeight * seqNum;
        }
        });
        });
    }

    std::vector<double> spscore (len, 0);
    tbb::parallel_for(tbb::blocked_range<int>(0, len), [&](tbb::blocked_range<int> r) {
    for (int s = r.begin(); s < r.end(); ++s) {
        for (int l = 0; l < 5; l++) spscore[s] += param.scoringMatrix[l][l] * (freq[s][l] * (freq[s][l] - 1) / 2);
        for (int l = 0; l < 5; l++) spscore[s] += param.gapExtend * (freq[s][l] * (freq[s][5]));
        for (int l = 0; l < 5; l++) {
            for (int m = 0; m < 5; m++) {
                if (m != l) spscore[s] += param.scoringMatrix[l][m] * (freq[s][m] * freq[s][l]);
            }
        }  
    }
    });
    double normScore = 0;
    for (int i = 0; i < len; i++) normScore += spscore[i];
    float matchAvg = 0;
    for (int i = 0; i < 4; ++i) matchAvg += param.scoringMatrix[i][i];
    matchAvg /= 4.0;
    normScore /= ((seqNum * (seqNum-1)) / 2);
    normScore /= len;
    normScore /= matchAvg;
    return normScore;
}   

void calculateProfileFreq(float* profile, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, char type, int32_t profileLen, int32_t profileSize, std::pair<int, int> startPos) {
    int32_t refLen = util->seqsLen[nodes.first->identifier];
    int32_t qryLen = util->seqsLen[nodes.second->identifier];    
    int32_t refNum = tree->allNodes[nodes.first->identifier]->msaIdx.size();
    int32_t qryNum = tree->allNodes[nodes.second->identifier]->msaIdx.size();
    
    if (util->nowProcess < 2) {
        float refWeight = 0.0, qryWeight = 0.0;
        int32_t numThreshold = SEQNUM_TH;
        for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) qryWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
        bool storeFreq = (((refNum >= numThreshold || qryNum > numThreshold) || 
                           (!tree->allNodes[nodes.first->identifier]->msaFreq.empty() || !tree->allNodes[nodes.second->identifier]->msaFreq.empty())) &&
                           (startPos.first == 0 && startPos.second == 0));
        storeFreq = (storeFreq & (util->nowProcess != 1));

        if (util->nowProcess != 1) {
            for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx)  refWeight += tree->allNodes[util->seqsName[sIdx]]->weight;
            if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
                for (auto sIdx: tree->allNodes[nodes.first->identifier]->msaIdx) { 
                    int storage = util->seqsStorage[sIdx];
                    std::string name = util->seqsName[sIdx];
                    float w = tree->allNodes[name]->weight / refWeight * refNum;
                    tbb::this_task_arena::isolate( [&]{
                    tbb::parallel_for(tbb::blocked_range<int>(startPos.first, refLen), [&](tbb::blocked_range<int> r) {
                    for (int s = r.begin(); s < r.end(); ++s) {
                        int t = s - startPos.first;
                        int letterIndex = letterIdx(type, toupper(util->alnStorage[storage][sIdx][s]));
                        profile[profileSize*t+letterIndex] += 1.0 * w;
                    }
                    });
                    });
                }
                if (storeFreq) {
                    std::vector<std::vector<float>> freq (refLen, std::vector<float>(profileSize,0.0));
                    tree->allNodes[nodes.first->identifier]->msaFreq = freq;
                    for (int s = 0; s < refLen; ++s) for (int t = 0; t < profileSize; ++t) tree->allNodes[nodes.first->identifier]->msaFreq[s][t] = profile[profileSize*s+t] / refNum * refWeight;
                }
            }
            else {
                for (int s = startPos.first; s < int(tree->allNodes[nodes.first->identifier]->msaFreq.size()); ++s) {
                    int t = s - startPos.first;
                    for (int v = 0; v < profileSize; ++v) profile[profileSize*t+v] = tree->allNodes[nodes.first->identifier]->msaFreq[s][v] / refWeight * refNum;
                }
            }
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) { 
            for (auto sIdx: tree->allNodes[nodes.second->identifier]->msaIdx) { 
                int storage = util->seqsStorage[sIdx];
                std::string name = util->seqsName[sIdx];
                float w = tree->allNodes[name]->weight / qryWeight * qryNum;
                tbb::this_task_arena::isolate( [&]{
                tbb::parallel_for(tbb::blocked_range<int>(startPos.second, qryLen), [&](tbb::blocked_range<int> r) {
                for (int s = r.begin(); s < r.end(); ++s) {
                    int t = s - startPos.second;
                    int letterIndex = letterIdx(type, toupper(util->alnStorage[storage][sIdx][s]));
                    profile[profileSize*(profileLen+t)+letterIndex]+=1.0*w;
                }
                });
                });
            }
            if (storeFreq) {
                std::vector<std::vector<float>> freq (qryLen, std::vector<float>(profileSize,0.0));
                tree->allNodes[nodes.second->identifier]->msaFreq = freq;
                for (int s = 0; s < qryLen; ++s) for (int t = 0; t < profileSize; ++t) tree->allNodes[nodes.second->identifier]->msaFreq[s][t] = profile[profileSize*(profileLen+s)+t] / qryNum * qryWeight;
            }
        }
        else {
            for (int s = startPos.second; s < int(tree->allNodes[nodes.second->identifier]->msaFreq.size()); ++s) {
                int t = s - startPos.first;
                for (int v = 0; v < profileSize; ++v) profile[profileSize*(profileLen+t)+v] = tree->allNodes[nodes.second->identifier]->msaFreq[s][v] / qryWeight * qryNum;
            }
        }
        
    }
    else { 
        float refNum_f = 0, qryNum_f = 0;
        int subtreeRef = util->seqsIdx[nodes.first->identifier];
        int subtreeQry = util->seqsIdx[nodes.second->identifier]; 
        if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.first->identifier]->msaFreq = util->profileFreq[subtreeRef];
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.second->identifier]->msaFreq = util->profileFreq[subtreeQry];
        }
        for (int t = 0; t < profileSize; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < profileSize; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
        
        for (int s = 0; s < tree->allNodes[nodes.first->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < profileSize; ++t) profile[profileSize*s+t] = tree->allNodes[nodes.first->identifier]->msaFreq[s][t];
        }
        for (int s = 0; s < tree->allNodes[nodes.second->identifier]->msaFreq.size(); ++s) {
            for (int t = 0; t < profileSize; ++t) profile[profileSize*(profileLen+s)+t] = tree->allNodes[nodes.second->identifier]->msaFreq[s][t];
        } 
    }
    return;
}

void calculatePSGOP(float* hostFreq, float* hostGapOp, float* hostGapEx, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, int32_t seqLen, std::pair<int32_t,int32_t> offset, std::pair<int32_t,int32_t> lens, Params& param) {
    int32_t refLen = lens.first;
    int32_t qryLen = lens.second;
    int32_t offsetf = offset.first;
    int32_t offsetg = offset.second; 
    int32_t refNum = tree->allNodes[nodes.first->identifier]->msaIdx.size();
    int32_t qryNum = tree->allNodes[nodes.second->identifier]->msaIdx.size();
    int32_t profileSize = (option->type == 'n') ? 6 : 22;
    if (util->nowProcess == 2) {
        float refNum_f = 0, qryNum_f = 0;
        int subtreeRef = util->seqsIdx[nodes.first->identifier];
        int subtreeQry = util->seqsIdx[nodes.second->identifier]; 
        if (tree->allNodes[nodes.first->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.first->identifier]->msaFreq = util->profileFreq[subtreeRef];
        }
        if (tree->allNodes[nodes.second->identifier]->msaFreq.empty()) {
            tree->allNodes[nodes.second->identifier]->msaFreq = util->profileFreq[subtreeQry];
        }
        for (int t = 0; t < 6; ++t) refNum_f += tree->allNodes[nodes.first->identifier]->msaFreq[0][t]; 
        for (int t = 0; t < 6; ++t) qryNum_f += tree->allNodes[nodes.second->identifier]->msaFreq[0][t]; 
        refNum = static_cast<int32_t>(round(refNum_f));
        qryNum = static_cast<int32_t>(round(qryNum_f));
    }
    // Clustalw's method
    float minGOP = param.gapOpen * 0.2;
    float scale = (option->type == 'n') ? 0.3 : 1.0;
    if (option->psgop) {
        for (int s = 0; s < seqLen; ++s) {
            if (s < refLen) {
                float gapRatio = hostFreq[offsetf+profileSize*s+profileSize-1];
                if (gapRatio > 0) {
                    hostGapOp[offsetg+s] = param.gapOpen * scale * ((refNum-gapRatio)*1.0 / refNum);
                    hostGapEx[offsetg+s] = param.gapExtend * 0.5;
                    // hostGapEx[offsetg+s] = param.gapOpen * ((refNum-gapRatio)*1.0 / refNum);
                }
                else {
                    hostGapOp[offsetg+s] = param.gapOpen;
                    hostGapEx[offsetg+s] = param.gapExtend;
                }
            }
            else {
                hostGapOp[offsetg+s] = 0.0;
                hostGapEx[offsetg+s] = 0.0;
            }
        }
        for (int s = 0; s < seqLen; ++s) {
            if (s < qryLen) {
                float gapRatio = hostFreq[offsetf+profileSize*(seqLen+s)+profileSize-1];
                if (gapRatio > 0) {
                    hostGapOp[offsetg+seqLen+s] = param.gapOpen * scale * ((qryNum-gapRatio) * 1.0 / qryNum);
                    hostGapEx[offsetg+seqLen+s] = param.gapExtend * 0.5;
                    // hostGapEx[offsetg+seqLen+s] = param.gapOpen * ((qryNum-gapRatio) * 1.0 / qryNum);
                }
                else {
                    hostGapOp[offsetg+seqLen+s] = param.gapOpen;
                    hostGapEx[offsetg+seqLen+s] = param.gapExtend;
                }
            }
            else {
                hostGapOp[offsetg+seqLen+s] = 0.0;
                hostGapEx[offsetg+seqLen+s] = 0.0;
            }
        }   
    }
    else {
        for (int s = 0; s < seqLen; ++s) {
            hostGapOp[offsetg+s] = (s < refLen) ? param.gapOpen   : 0;
            hostGapEx[offsetg+s] = (s < refLen) ? param.gapExtend : 0;
        }
        for (int s = 0; s < seqLen; ++s) {
            hostGapOp[offsetg+seqLen+s] = (s < qryLen) ? param.gapOpen   : 0;
            hostGapEx[offsetg+seqLen+s] = (s < qryLen) ? param.gapExtend : 0;
        }
    }   
    return;
}

void removeGappyColumns(float* hostFreq, Tree* tree, std::pair<Node*, Node*>& nodes, msa::utility* util, msa::option* option, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, int32_t profileLen, int32_t maxProfileLen, std::pair<int,int>& lens, std::pair<int,int>& rawLens) {
    float gappyVertical = option->gappyVertical;
    int gappyHorizon = option->gappyHorizon, gappyLength;
    if ((gappyVertical == 1 || gappyHorizon < 1)) {
        lens.first = std::min(lens.first, maxProfileLen);
        lens.second = std::min(lens.second, maxProfileLen);
        rawLens.first = lens.first;
        rawLens.second = lens.second;
        return;
    }
    int32_t profileSize = (option->type == 'n') ? 6 : 22;
    int32_t rawIdx = 0, newIdx = 0;
    int32_t gapRef = 0, gapQry = 0;
    int32_t refLen = lens.first, qryLen = lens.second;
    float refNum_f = 0, qryNum_f = 0;
    for (int t = 0; t < profileSize; ++t) refNum_f += hostFreq[t]; 
    for (int t = 0; t < profileSize; ++t) qryNum_f += hostFreq[profileSize*profileLen+t];       
    int32_t refNum = (util->nowProcess < 2) ? tree->allNodes[nodes.first->identifier]->msaIdx.size() : static_cast<int32_t>(round(refNum_f));
    int32_t qryNum = (util->nowProcess < 2) ? tree->allNodes[nodes.second->identifier]->msaIdx.size(): static_cast<int32_t>(round(qryNum_f));
    int32_t numTH = static_cast<int32_t>(1/(1-gappyVertical));
    while (true) {
        if (rawIdx >= refLen) {
            for (int i = newIdx; i < refLen; ++i) for (int j = 0; j < profileSize; ++j) hostFreq[profileSize*i+j] = 0;
            break;
        }
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < refLen; ++gappyLength) {
            if ((hostFreq[profileSize*gappyLength+profileSize-1]+hostFreq[profileSize*gappyLength+profileSize-2])/refNum <= gappyVertical) break;
            if ((hostFreq[profileSize*gappyLength+profileSize-1]+hostFreq[profileSize*gappyLength+profileSize-2])/refNum > gappyVertical && refNum < numTH) break;
        }
        if (gappyLength - tempStart >= gappyHorizon) {
            gappyColumns.first.push(std::make_pair(tempStart, gappyLength-tempStart)); 
            gapRef += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, refLen); ++rawIdx) {
            if (newIdx == maxProfileLen) break;
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*newIdx+t] = hostFreq[profileSize*rawIdx+t];
            ++newIdx;
            
        } 
        if (newIdx == maxProfileLen) break;
    }
    lens.first = newIdx;
    rawLens.first = rawIdx;
    rawIdx = 0; newIdx = 0;
    while (true) {
        if (rawIdx >= qryLen) {
            for (int i = newIdx; i < qryLen; ++i) for (int j = 0; j < profileSize; ++j) hostFreq[profileSize*(profileLen+i)+j] = 0;
            break;
        }
        int tempStart = rawIdx;
        bool onlyN = false;
        for (gappyLength = rawIdx; gappyLength < qryLen; ++gappyLength) {
            if ((hostFreq[profileSize*(profileLen+gappyLength)+profileSize-1]+hostFreq[profileSize*(profileLen+gappyLength)+profileSize-2])/qryNum <= gappyVertical) break;
            if ((hostFreq[profileSize*(profileLen+gappyLength)+profileSize-1]+hostFreq[profileSize*(profileLen+gappyLength)+profileSize-2])/qryNum > gappyVertical && qryNum < numTH) break;
        }
        if (gappyLength - tempStart >= gappyHorizon) {
            gappyColumns.second.push(std::make_pair(tempStart, gappyLength-tempStart)); 
            gapQry += gappyLength-rawIdx;
            rawIdx += gappyLength-rawIdx;
        }
        for (rawIdx = rawIdx; rawIdx < std::min(gappyLength+1, qryLen); ++rawIdx) {
            if (newIdx == maxProfileLen) break;
            for (int t = 0; t < profileSize; ++t) hostFreq[profileSize*(profileLen+newIdx)+t] = hostFreq[profileSize*(profileLen+rawIdx)+t];
            ++newIdx;
            
        } 
        if (newIdx == maxProfileLen) break;
    }
    lens.second = newIdx;
    rawLens.second = rawIdx;
    return;
}

void addGappyColumnsBack(std::vector<int8_t>& aln_old, std::vector<int8_t>& aln, std::pair<std::queue<std::pair<int, int>>, std::queue<std::pair<int, int>>>& gappyColumns, std::pair<int, int>& debugIdx) {
    int rIdx = 0, qIdx = 0, j = 0; 
    int arIdx = 0, aqIdx = 0;
    bool endAln = (debugIdx.first <= 0);
    bool mergeCols = (debugIdx.second > 0);
    int maxrIdx = std::abs(debugIdx.first), maxqIdx = std::abs(debugIdx.second);
    bool preN = false;
    int offset = (aln.empty()) ? 0 : 1;
    int lastDir = aln_old.back();

    while ((j < aln_old.size() || (!gappyColumns.first.empty() || !gappyColumns.second.empty())) && (arIdx < maxrIdx || aqIdx < maxqIdx)) {
        bool gapR = (gappyColumns.first.empty())  ? false : ((rIdx == gappyColumns.first.front().first - offset) &&  (arIdx < maxrIdx));
        bool gapQ = (gappyColumns.second.empty()) ? false : ((qIdx == gappyColumns.second.front().first - offset) && (aqIdx < maxqIdx));
        int gapRLen = (!gapR) ? 0 : gappyColumns.first.front().second;
        int gapQLen = (!gapQ) ? 0 : gappyColumns.second.front().second;
        bool allNR = gapRLen < 0;
        bool allNQ = gapQLen < 0;
        gapRLen = std::abs(gapRLen);
        gapQLen = std::abs(gapQLen);
        if (gapR || gapQ) {
            if (gapR && !gapQ) {
                if (!allNR) {
                    for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                    gappyColumns.first.pop();
                    preN = allNR;
                }
                else {
                    int delLen = 0, qTerminal = (gappyColumns.second.empty()) ? aln_old.size() : gappyColumns.second.front().first; 
                    for (int k = j; k < aln_old.size(); ++k) {
                        if ((aln_old[k] & 0xFFFF) != 1 || (qIdx+delLen) >= qTerminal) break;
                        ++delLen; 
                    }
                    if (delLen > gapRLen) {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        gappyColumns.first.pop();
                        preN = false;
                    }
                    else {
                        for (int g = 0; g < delLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        for (int g = delLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                        gappyColumns.first.pop();
                        preN = true;
                    }
                }   
            }
            else if (!gapR && gapQ) {
                if (!allNQ) {
                    for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                    gappyColumns.second.pop();
                    preN = allNQ;
                }
                else {
                    int delLen = 0, rTerminal = (gappyColumns.first.empty()) ? aln_old.size() : gappyColumns.first.front().first; ; 
                    for (int k = j; k < aln_old.size(); ++k) {
                        if ((aln_old[k] & 0xFFFF) != 2 || (rIdx+delLen) >= rTerminal) break;
                        ++delLen; 
                    }
                    if (delLen > gapQLen) {
                        for (int g = 0; g < gapQLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        gappyColumns.second.pop();
                        preN = false;
                    }
                    else {
                        for (int g = 0; g < delLen; ++g) {++rIdx; ++qIdx; ++j; aln.push_back(0);}
                        for (int g = delLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                        gappyColumns.second.pop();
                        preN = true;
                    }
                }   
            }
            else if (gapR && gapQ) {
                if (!allNR && !allNQ) {
                    if (mergeCols) {
                        if (gapRLen >= gapQLen) {
                            for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                            for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                        }
                        else {
                            for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                            for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                        }
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                    }
                    gappyColumns.first.pop();
                    gappyColumns.second.pop();
                    preN = false;
                }
                else if (allNR && !allNQ) {
                    if (!preN) {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                        gappyColumns.second.pop();
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        gappyColumns.first.pop();
                    }   
                }
                else if (!allNR && allNQ) {
                    if (!preN) {
                        for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                        gappyColumns.first.pop();
                    }
                    else {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
                        gappyColumns.second.pop();
                    }
                }
                else {
                    if (gapRLen > gapQLen) {
                        for (int g = 0; g < gapQLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                        gappyColumns.first.front().first = rIdx;
                        gappyColumns.first.front().second = gapRLen - gapQLen;
                        gappyColumns.second.pop();
                    }
                    else {
                        for (int g = 0; g < gapRLen; ++g) {++qIdx; ++rIdx; aln.push_back(0);}
                        gappyColumns.second.front().first = qIdx;
                        gappyColumns.second.front().second = gapQLen - gapRLen;
                        gappyColumns.first.pop();
                    }
                    preN = true;
                }
            }
        }
        else {
            switch ((aln_old[j] & 0xFFFF)) {
                case 0: ++rIdx; ++qIdx; ++arIdx; ++aqIdx; aln.push_back(0); break;
                case 1: ++qIdx;         ++aqIdx;          aln.push_back(1); break;
                case 2: ++rIdx;         ++arIdx;          aln.push_back(2); break;
            }
            ++j;
        }
        if (gappyColumns.first.empty() && gappyColumns.second.empty()) {
            for (j = j; j < aln_old.size(); ++j) {
                switch ((aln_old[j] & 0xFFFF)) {
                    case 0: ++rIdx; ++qIdx; ++arIdx; ++aqIdx; aln.push_back(0); break;
                    case 1: ++qIdx;         ++aqIdx;          aln.push_back(1); break;
                    case 2: ++rIdx;         ++arIdx;          aln.push_back(2); break;
                }
            }
        }
    }
    if (endAln && (!gappyColumns.first.empty() || !gappyColumns.second.empty())) {
        bool gapR = !gappyColumns.first.empty();
        bool gapQ = !gappyColumns.second.empty();
        int gapRLen = (!gapR) ? 0 : gappyColumns.first.front().second;
        int gapQLen = (!gapQ) ? 0 : gappyColumns.second.front().second;
        if ((gapR && !gapQ)) {
            for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
            gappyColumns.first.pop();
        }
        else if ((gapQ && !gapR)) {
            for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            gappyColumns.second.pop();
        }
        else {
            if (mergeCols) {
                if (gapRLen >= gapQLen) {
                    for (int g = 0; g < gapQLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                    for (int g = gapQLen; g < gapRLen; ++g) {++rIdx;         aln.push_back(2);}
                }
                else {
                    for (int g = 0; g < gapRLen; ++g)       {++rIdx; ++qIdx; aln.push_back(0);}
                    for (int g = gapRLen; g < gapQLen; ++g) {++qIdx;         aln.push_back(1);}
                }
            }
            else {
                for (int g = 0; g < gapRLen; ++g) {++rIdx; aln.push_back(2);}
                for (int g = 0; g < gapQLen; ++g) {++qIdx; aln.push_back(1);}
            }
            gappyColumns.first.pop();
            gappyColumns.second.pop();
        }
    }
    debugIdx.first = rIdx; debugIdx.second = qIdx;
    return;
}

void mergeInsertions (Tree* tree, Node* nodeRef, std::vector<Node*>& nodes, msa::utility* util, std::vector<std::vector<int8_t>>& alnBad) {
    std::vector<int> alnIdx (alnBad.size(), 0);
    std::vector<std::vector<int8_t>> alnNew (alnBad.size());
    std::vector<int8_t> alnRef;
    while (true) {
        bool hasInsert = false;
        bool reachEnd = true;
        for (int i = 0; i < alnBad.size(); ++i) {
            if (!alnBad[i].empty()) {
                if (alnIdx[i] < alnBad[i].size()) {
                    reachEnd = false;
                    if (alnBad[i][alnIdx[i]] == 1) {
                        hasInsert = true;
                        break;
                    }
                }
            }
        }
        if (reachEnd) break;
        if (hasInsert) {
            // std::cout << "Insert at " << alnRef.size() << "\n";
            for (int i = 0; i < alnBad.size(); ++i) {
                if (!alnBad[i].empty()) {
                    if (alnBad[i][alnIdx[i]] == 1) {
                        alnNew[i].push_back(0);
                        alnIdx[i] += 1;
                        
                    } 
                    else {
                        alnNew[i].push_back(2);
                    }
                }
            }
            alnRef.push_back(2);
        }
        else {
            for (int i = 0; i < alnBad.size(); ++i) {
                if (!alnBad[i].empty()) {
                    alnNew[i].push_back(alnBad[i][alnIdx[i]]);
                    alnIdx[i] += 1;
                    
                }
            }
            alnRef.push_back(0);
        }
    } 
    std::cout << "Length after insertions: " << alnRef.size() << ", before insertions: " << util->seqsLen[nodeRef->identifier] << '\n';
    if (alnRef.size() != util->seqsLen[nodeRef->identifier]) { // Update large profile
        tbb::this_task_arena::isolate( [&]{
        tbb::parallel_for(tbb::blocked_range<int>(0, tree->allNodes[nodeRef->identifier]->msaIdx.size()), [&](tbb::blocked_range<int> range) {
        for (int idx = range.begin(); idx < range.end(); ++idx) {
            int sIdx = tree->allNodes[nodeRef->identifier]->msaIdx[idx];  
            util->memCheck(alnRef.size(), sIdx);
            int orgIdx = 0;
            int storeFrom = util->seqsStorage[sIdx];
            int storeTo = 1 - storeFrom;
            for (int k = 0; k < alnRef.size(); ++k) {
                if (alnRef[k] == 0) {
                    util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                    orgIdx++;
                }
                else if (alnRef[k] == 2) {
                    util->alnStorage[storeTo][sIdx][k] = '-';
                }
            }
            util->changeStorage(sIdx);
        }
        });
        });
        util->seqsLen[nodeRef->identifier] = alnRef.size();
    
    }
    tbb::this_task_arena::isolate( [&]{
    tbb::parallel_for(tbb::blocked_range<int>(0, nodes.size()), [&](tbb::blocked_range<int> range) {
    for (int i = range.begin(); i < range.end(); ++i) {
        if (!alnNew[i].empty()) {
            for (int idx = 0; idx < tree->allNodes[nodes[i]->identifier]->msaIdx.size(); ++idx) {
                int sIdx = tree->allNodes[nodes[i]->identifier]->msaIdx[idx];  
                util->memCheck(alnNew[i].size(), sIdx);
                int orgIdx = 0;
                int storeFrom = util->seqsStorage[sIdx];
                int storeTo = 1 - storeFrom;
                for (int k = 0; k < alnNew[i].size(); ++k) {
                    if (alnNew[i][k] == 0) {
                        util->alnStorage[storeTo][sIdx][k] = util->alnStorage[storeFrom][sIdx][orgIdx];
                        orgIdx++;
                    }
                    else if (alnNew[i][k] == 2) {
                        util->alnStorage[storeTo][sIdx][k] = '-';
                    }
                }
                util->changeStorage(sIdx);
            }
            util->seqsLen[nodes[i]->identifier] = alnNew[i].size();
        }
    }
    });
    });

    for (int i = 0; i < nodes.size(); ++i) {
        if (!alnBad[i].empty()) {
            for (auto q: tree->allNodes[nodes[i]->identifier]->msaIdx) {
                tree->allNodes[nodeRef->identifier]->msaIdx.push_back(q);
            }
            tree->allNodes[nodes[i]->identifier]->msaIdx.clear();
        }
    }
    tree->allNodes[nodeRef->identifier]->msaFreq.clear();
    
    return;
}

/*
void getReAlignRegion(std::vector<int8_t>& aln, std::vector<bool>& conserved, std::vector<reAlign*>& realign,  std::vector<std::pair<int,int>>& regions) {
    const int conservedSeed = 10; // minimum consecutive conserved columns to be consider as a seed
    const int minRealignLen = 5;
    int conservedLen = 0;
    bool whichRegion = 0; // 0: conserved, 1: non-conserved
    int noncon_start = 0;
    float noncon_num = 0;
    float noncon_ratio = 0.5;
    for (int i = 0; i < conserved.size(); ++i) {
        // std::cout << whichRegion << '\t' << conserved[i] << '\t' << conservedLen << '\n';
        if (!whichRegion) {
            if (!conserved[i]) {
                conservedLen = 0;
                whichRegion = 1;
                noncon_start = i;
                noncon_num = 1;
            }
        }
        else {
            if (conserved[i]) {
                conservedLen++;
                if (conservedLen >= conservedSeed) {
                    if (i-conservedLen-noncon_start+1 >= minRealignLen) {
                        // if (noncon_num / static_cast<float>(i-conservedLen-noncon_start+1) > noncon_ratio)
                        regions.push_back(std::make_pair(noncon_start, i-conservedLen));
                    }
                    whichRegion = 0;
                }
            }
            else {
                conservedLen = 0;
                noncon_num++;
            }
        }
    }
    // if (whichRegion && (conserved.size()-noncon_start >= minRealignLen)) {
    //     regions.push_back(std::make_pair(noncon_start, conserved.size()-1));
    // }
    int rIdx = 0, qIdx = 0, regionIdx = 0;
    int rs = 0, re = 0, qs = 0, qe = 0;
    for (int i = 0; i < aln.size(); ++i) {
        if (regionIdx == regions.size()) break;
        if (i == regions[regionIdx].first) {rs = rIdx; qs = qIdx;}
        if (i == regions[regionIdx].second) {
            re = rIdx; qe = qIdx;
            if (aln[i] == 2) qe--;
            if (aln[i] == 1) re--;
            reAlign* temp = new reAlign(rs, re, qs ,qe);
            realign.push_back(temp);
            ++regionIdx;
        }
        if (aln[i] == 0)     {rIdx++; qIdx++;}
        else if (aln[i] == 1) qIdx++;
        else if (aln[i] == 2) rIdx++;
    }
    // for (int i = 0; i < regions.size(); ++i) std::cout << '(' << regions[i].first << ',' << regions[i].second << ')' << '\n';
    for (auto r: realign) printf("(%d,%d)-(%d,%d)\n", r->rs, r->re, r->qs, r->qe);
    return;
}

void conservedColumns(std::vector<std::vector<float>>& freqs, std::vector<bool>& conserved, float weight) {
    float entropy = 0.0;
    float max_entropy = 1.0; // larger than this means not very conserved
    int count = 0;
    if (freqs[0].size() == 6) {
        float group_ratio [4] = {0.0};
        for (int k = 0; k < freqs.size(); ++k) {
            for (int i = 0; i < 4; ++i) group_ratio[i] += (freqs[k][i] / weight);
            for (int i = 0; i < 4; ++i) if (group_ratio[i] > 0) entropy -= (group_ratio[i] * std::log2(group_ratio[i]));
            for (int i = 0; i < 4; ++i) group_ratio[i] = 0;
            conserved[k] = (entropy < max_entropy);
        }
    }
    else {
        // Amino Acids             |  Group
        // A, F, I, L, M, V, Y, W  |   0
        // N, Q, S, T              |   1 
        // C, G, P                 |   2
        // H, K, R                 |   3
        // D, E                    |   4
        //                                  A      C      D      E      F      G      H      I      K     L       M
        std::map<int, int> amino_groups = {{0,0}, {1,2}, {2,4}, {3,4}, {4,0}, {5,2}, {6,3}, {7,0}, {8,3}, {9,0}, {10,0}, 
        //                                  N       P       Q       R       S       T       V       W       Y
                                           {11,1}, {12,2}, {13,1}, {14,3}, {15,1}, {16,1}, {17,0}, {18,0}, {19,0}};
        float group_ratio [5] = {0.0};
        for (int k = 0; k < freqs.size(); ++k) {
            entropy = 0;
            for (int i = 0; i < 20; ++i) group_ratio[amino_groups[i]] += (freqs[k][i] / weight);
            // for (int i = 0; i < 5; ++i) std::cout << group_ratio[i] << '\t';
            // std::cout << (group_ratio[0]+group_ratio[1]+group_ratio[2]+group_ratio[3]+group_ratio[4]) << '\n';
            for (int i = 0; i < 5; ++i) if (group_ratio[i] > 0) entropy -= (group_ratio[i] * std::log2(group_ratio[i]));
            for (int i = 0; i < 5; ++i) group_ratio[i] = 0;
            conserved[k] = (entropy < max_entropy);
            if (entropy >= max_entropy) count += 1;
            // std::cout << k << ": " << entropy << '\n';
        }
    }
    // if (count > conserved.size() / 10) {
    //     for (auto i: conserved) std::cout << i;
    //     std::cout << '\n';
    // }
    return;
}
*/
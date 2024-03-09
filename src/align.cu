#ifndef ALIGN_HPP
#include "align.cuh"
#define INF (1<<20)
#endif

int8_t updateState(STATE& currentHState, STATE& currentIState, STATE& currentDState)
{
    int8_t currentState = 0;
    switch (currentHState)
    {
    case STATE::HI:
        currentState = (currentState | 0x01); break;
    case STATE::HD:
        currentState = (currentState | 0x02); break;
    default: 
        currentState = currentState; break;
    }

    switch (currentIState)
    {
    case STATE::II:
        currentState = (currentState | 0x04);  break;
    default:
        currentState = currentState; break;
    }
    
    switch (currentDState)
    {
    case STATE::DD:
        currentState = (currentState | 0x08); break;
    default:
        currentState = currentState; break;
    }
    return currentState;
}

__device__ int8_t updateState_cuda(STATE& currentHState, STATE& currentIState, STATE& currentDState)
{
    int8_t currentState = 0;
    switch (currentHState)
    {
    case STATE::HI:
        currentState = (currentState | 0x01); break;
    case STATE::HD:
        currentState = (currentState | 0x02); break;
    default: 
        currentState = currentState; break;
    }

    switch (currentIState)
    {
    case STATE::II:
        currentState = (currentState | 0x04);  break;
    default:
        currentState = currentState; break;
    }
    
    switch (currentDState)
    {
    case STATE::DD:
        currentState = (currentState | 0x08); break;
    default:
        currentState = currentState; break;
    }
    return currentState;
}


void tracebackSeqtoSeq
(
    int8_t state,
    std::vector<int8_t>& TB,
    std::vector<int32_t>& wfLL,
    std::vector<int32_t>& wfLen,
    std::pair<std::string, std::string>& alignment,
    std::string& ref,
    std::string& query
){
    int32_t refLen = ref.size();
    int32_t queryLen = query.size();
    int32_t refIndex = refLen-1;
    int32_t queryIndex = queryLen-1;
    
    int32_t k=wfLL.size()-1;
    int32_t tbIndex = TB.size()-1;
    int8_t currentTB;
    int8_t dir;
    std::string refAlign, queryAlign;
    while (refIndex>=0 && queryIndex>=0)
    {

        currentTB = TB[tbIndex];
        if (state == 0) 
        { // Current State M
            state = currentTB & 0x03;
            if (state == 0) dir = 0;
            else if (state == 1)
            {
                dir = 1;
                if (currentTB & 0x04) state = 1;
                else state = 0;   
            }
            else
            {
                dir = 2;
                if (currentTB & 0x08) state = 2;
                else state = 0;
            }
        }
        else if (state == 1) 
        { // Current State I
            dir = 1;
            if (currentTB & 0x04) state = 1;
            else state = 0;
        } 
        else 
        { // Current State D
            dir = 2;
            if (currentTB & 0x08) state = 2;
            else state = 0;
        }

        tbIndex -= (refIndex - wfLL[k] + 1 + wfLen[k-1]);
        if (dir == 0) {tbIndex -= wfLen[k-2]; tbIndex += refIndex - wfLL[k-2]; refAlign+=ref[refIndex]; queryAlign+=query[queryIndex]; k--;k--;queryIndex--;refIndex--;}
        else if (dir == 1) {tbIndex += (refIndex - wfLL[k-1] + 1);refAlign+="-"; queryAlign+=query[queryIndex]; k--;queryIndex--;}
        else {tbIndex += (refIndex - wfLL[k-1]); refAlign+=ref[refIndex]; queryAlign+="-"; k--;refIndex--;}

    }

    while (refIndex>=0)
    {

        for (size_t i=0; i<ref.size(); i++) refAlign+=ref[refIndex];  
        for (size_t i=0; i<query.size(); i++) queryAlign+="-";  
        k--;refIndex--;
    }

    while (queryIndex>=0)
    {

        for (size_t i=0; i<ref.size(); i++) refAlign+="-"; 
        for (size_t i=0; i<query.size(); i++) queryAlign+=query[queryIndex]; 
        k--;queryIndex--;
    }

    refAlign = std::string(refAlign.rbegin(), refAlign.rend());
    queryAlign = std::string(queryAlign.rbegin(), queryAlign.rend());

    alignment = std::make_pair(refAlign, queryAlign);

}


void alignSeqToSeq
(
    std::string& ref,
    std::string& query,
    Params& param,
    std::pair<std::string, std::string>& alignment
){
    int32_t refLen = ref.size();
    int32_t queryLen = query.size();
    if (refLen<=0 || queryLen<=0) {fprintf(stderr, "Error: Ref/Query length <=0"); exit(1);}
    
    int maxWFLen = 0; //wavefront length
    if (param.marker != 0) maxWFLen = param.marker+2;
    else maxWFLen = ref.size() + query.size() + 2;

    int score = 0;
    int32_t *H[3], *I[2], *D[2];
    int32_t L[3], U[3];
    std::vector<int32_t> wfLL, wfLen;

    //Output
    int8_t state=0;
    std::vector<int8_t> TB;

    for(size_t i=0; i<3; i++)
    {
        H[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        if (i<2) 
        {
            I[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
            D[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        }
        L[i]=0;
        U[i]=0;
    }

    for (size_t i=0; i<3; i++)
    {
        for (size_t j=0; j<(size_t)maxWFLen; j++)
        {
            H[i][j] = -INF;
            if (i<2) {I[i][j] = 0; D[i][j] = -INF;} 
        }
    }

    for (int32_t k=0; k<refLen+queryLen+1; k++)
    {
        L[k%3] = (k<=queryLen)?0:k-queryLen;
        U[k%3] = (k<=refLen)?k:refLen;
        wfLL.push_back(L[k%3]);
        wfLen.push_back(U[k%3]-L[k%3]+1);
        // std::cout << L[k%3] << "-" << U[k%3] << "\n";
        for(int32_t i=L[k%3]; i<U[k%3]+1; i++) // i->Ref Index
        {   
            int32_t j=(k-i); //j->Query Index
            int32_t match = -INF, insOp = -INF, delOp = -INF, insExt = -INF, delExt = -INF;
            int32_t offset = i-L[k%3];
            int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1;
            int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
            int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
            STATE currentHState = STATE::HH;
            STATE currentIState = STATE::IH;
            STATE currentDState = STATE::DH;

            if (k==0) match = 0;

            if (offsetDiag>=0 && j>0)
            {
                if (ref[i-1]==query[j-1]) match = H[(k+1)%3][offsetDiag] + param.match;
                else match = H[(k+1)%3][offsetDiag] + param.mismatch;
            }
            
            if (offsetUp >= 0)
            {
                insOp = H[(k+2)%3][offsetUp] + param.gapOpen;
                insExt = I[(k+1)%2][offsetUp] + param.gapExtend;
            }

            if (offsetLeft >=0)
            {
                delOp = H[(k+2)%3][offsetLeft] + param.gapOpen;
                delExt = D[(k+1)%2][offsetLeft] + param.gapExtend;
            }

            I[k%2][offset] =  insOp;
            D[k%2][offset] =  delOp;
            
            if (insExt >= insOp) 
            {
                I[k%2][offset] = insExt;
                currentIState = STATE::II;
            }
            if (delExt >= delOp) 
            {
                D[k%2][offset] = delExt;
                currentDState = STATE::DD;
            }
            if (match > I[k%2][offset]) 
            {
                if (match > D[k%2][offset]) H[k%3][offset] = match;
                else 
                {
                    H[k%3][offset] = D[k%2][offset];
                    currentHState = STATE::HD;
                }
            }
            else if (I[k%2][offset] > D[k%2][offset]) 
            {
                H[k%3][offset] = I[k%2][offset];
                currentHState = STATE::HI;
            }
            else 
            {
                H[k%3][offset] = D[k%2][offset];
                currentHState = STATE::HD;
            }

            TB.push_back(updateState(currentHState, currentIState, currentDState));

            score = H[k%3][offset];
            state = currentHState;
        }
    }

    // Deallocate memory for scores
    for (size_t sIndx=0; sIndx<3; sIndx++) {
        std::free(H[sIndx]);
        if (sIndx < 2) {
            std::free(I[sIndx]);
            std::free(D[sIndx]);
        }
    }

    tracebackSeqtoSeq (state, TB, wfLL, wfLen, alignment, ref, query);
    printf("Score: %d\n", score);
    return;

}


void getCharFrequency
(
    std::vector<std::string>& ref,
    std::vector<std::vector<float>>& charFreq
){
    charFreq.resize(ref[0].size());
    for (size_t i=0; i<ref[0].size(); i++) 
    {
        charFreq[i].resize(5); //Only A,C,G,T,- supported
        for (int j=0; j<5; j++) charFreq[i][j]=0;

    }

    for (size_t refIndex=0; refIndex<ref[0].size(); refIndex++)
    {
        for (size_t i=0; i<ref.size(); i++)
        {
            if (ref[i][refIndex]=='A' || ref[i][refIndex]=='a') charFreq[refIndex][0]+=1;
            else if (ref[i][refIndex]=='C' || ref[i][refIndex]=='c') charFreq[refIndex][1]+=1;
            else if (ref[i][refIndex]=='G' || ref[i][refIndex]=='g') charFreq[refIndex][2]+=1;
            else if (ref[i][refIndex]=='T' || ref[i][refIndex]=='t') charFreq[refIndex][3]+=1;
            else charFreq[refIndex][4]+=1;
        }
    }
}

bool getSimilarityScore(int32_t refIndex, int32_t queryIdx, std::vector<std::vector<float>>& charFreqRef, std::vector<std::vector<float>>& charFreqQuery)
{
    float score=0;
    for (int i=0; i<5; i++) score+= std::sqrt(charFreqRef[refIndex][i]*charFreqQuery[queryIdx][i]);

    return (score>0.95); //ToDo: How to clerverly decide cut-off point
} 

void tracebackGrpToGrp
(
    int8_t state,
    std::vector<int8_t>& TB,
    std::vector<int32_t>& wfLL,
    std::vector<int32_t>& wfLen,
    std::pair<std::vector<std::string>, std::vector<std::string>>& alignment,
    std::vector<std::string>& ref,
    std::vector<std::string>& query
){
    int32_t refLen = ref[0].size();
    int32_t queryLen = query[0].size();
    int32_t refIndex = refLen-1;
    int32_t queryIndex = queryLen-1;
    
    int32_t k=wfLL.size()-1;
    int32_t tbIndex = TB.size()-1;
    int8_t currentTB;
    int8_t dir;
    std::vector<std::string> refAlign, queryAlign;
    refAlign.resize(ref.size()); queryAlign.resize(query.size());


    while (refIndex>=0 && queryIndex>=0)
    {
        currentTB = TB[tbIndex];
        if (state == 0) 
        { // Current State M
            state = currentTB & 0x03;
            if (state == 0) dir = 0;
            else if (state == 1)
            {
                dir = 1;
                if (currentTB & 0x04) state = 1;
                else state = 0;   
            }
            else
            {
                dir = 2;
                if (currentTB & 0x08) state = 2;
                else state = 0;
            }
        }
        else if (state == 1) 
        { // Current State I
            dir = 1;
            if (currentTB & 0x04) state = 1;
            else state = 0;
        } 
        else 
        { // Current State D
            dir = 2;
            if (currentTB & 0x08) state = 2;
            else state = 0;
        }

        tbIndex -= (refIndex - wfLL[k] + 1 + wfLen[k-1]);
        if (dir == 0) 
        {
            tbIndex -= wfLen[k-2]; tbIndex += refIndex - wfLL[k-2];
            for (size_t i=0; i<ref.size(); i++) refAlign[i]+=ref[i][refIndex]; 
            for (size_t i=0; i<query.size(); i++) queryAlign[i]+=query[i][queryIndex]; 
            k--;k--;queryIndex--;refIndex--;
        }
        else if (dir == 1) 
        {
            tbIndex += (refIndex - wfLL[k-1] + 1);
            for (size_t i=0; i<ref.size(); i++) refAlign[i]+="-"; 
            for (size_t i=0; i<query.size(); i++) queryAlign[i]+=query[i][queryIndex]; 
            k--;queryIndex--;
        }
        else 
        {
            tbIndex += (refIndex - wfLL[k-1]); 
            for (size_t i=0; i<ref.size(); i++) refAlign[i]+=ref[i][refIndex];  
            for (size_t i=0; i<query.size(); i++) queryAlign[i]+="-";  
            k--;refIndex--;
        }

    }


    while (refIndex>=0)
    {

        for (size_t i=0; i<ref.size(); i++) refAlign[i]+=ref[i][refIndex];  
        for (size_t i=0; i<query.size(); i++) queryAlign[i]+="-";  
        k--;refIndex--;
    }

    while (queryIndex>=0)
    {

        for (size_t i=0; i<ref.size(); i++) refAlign[i]+="-"; 
        for (size_t i=0; i<query.size(); i++) queryAlign[i]+=query[i][queryIndex]; 
        k--;queryIndex--;
    }

    // if (refIndex == 0)
    // {
    //     for (size_t i=0; i<ref.size(); i++) refAlign[i]+=ref[i][refIndex];  
    //     for (size_t i=0; i<query.size(); i++) queryAlign[i]+="-"; 
    //     k--;refIndex--;
    // }

    // if (queryIndex == 0)
    // {
    //     for (size_t i=0; i<ref.size(); i++) refAlign[i]+="-"; 
    //     for (size_t i=0; i<query.size(); i++) queryAlign[i]+=query[i][queryIndex]; 
    //     k--;queryIndex--;
    // }

    for (auto &s: refAlign) s = std::string(s.rbegin(), s.rend());
    for (auto &s: queryAlign) s = std::string(s.rbegin(), s.rend());

    alignment = std::make_pair(refAlign, queryAlign);

}


/*
void alignGrpToGrp
(
    std::vector<std::string>& ref,
    std::vector<std::string>& query,
    Params& param,
    std::pair<std::vector<std::string>, std::vector<std::string>>& alignment
){
    if (ref.size() < 1 || query.size() < 1) {fprintf(stderr, "Error: Number of Ref/Query <= 0\n"); exit(1);}
    int32_t refLen = ref[0].size();
    int32_t queryLen = query[0].size();
    if (refLen<=0 || queryLen<=0) {fprintf(stderr, "Error: Ref/Query length <= 0\n"); exit(1);}
    
    int maxWFLen = 0; //wavefront length
    if (param.marker != 0) maxWFLen = param.marker+2;
    else maxWFLen = ref[0].size() + query[0].size() + 2;


    // Get similarity score of input groups
    std::vector<std::vector<float>> charFreqRef, charFreqQuery; 
    getCharFrequency(ref, charFreqRef);
    getCharFrequency(query, charFreqQuery);    

    // std::cout << charFreqRef[0].size() << "-" << charFreqQuery[0].size() << "\n";
    int score = 0;
    int32_t *H[3], *I[2], *D[2];
    int32_t L[3], U[3];
    std::vector<int32_t> wfLL, wfLen;

    //Output
    int8_t state=0;
    std::vector<int8_t> TB;

    for(size_t i=0; i<3; i++)
    {
        H[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        if (i<2) 
        {
            I[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
            D[i] = (int32_t*) std::malloc(sizeof(int32_t)*maxWFLen);
        }
        L[i]=0;
        U[i]=0;
    }

    for (size_t i=0; i<3; i++)
    {
        for (size_t j=0; j<(size_t)maxWFLen; j++)
        {
            H[i][j] = -INF;
            if (i<2) {I[i][j] = 0; D[i][j] = -INF;} 
        }
    }

    for (int32_t k=0; k<refLen+queryLen+1; k++)
    {
        L[k%3] = (k<=queryLen)?0:k-queryLen;
        U[k%3] = (k<=refLen)?k:refLen;
        wfLL.push_back(L[k%3]);
        wfLen.push_back(U[k%3]-L[k%3]+1);
        for(int32_t i=L[k%3]; i<U[k%3]+1; i++) // i->Ref Index
        {   
            int32_t j=(k-i); //j->Query Index
            int32_t match = -INF, insOp = -INF, delOp = -INF, insExt = -INF, delExt = -INF;
            int32_t offset = i-L[k%3];
            int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1;
            int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
            int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
            STATE currentHState = STATE::HH;
            STATE currentIState = STATE::IH;
            STATE currentDState = STATE::DH;


            if (k==0) match=0;

            if (offsetDiag>=0 && j>0)
            {
                // std::cout << i << "-" << j << "-" << offsetDiag << "\n";
                bool groupMatch = getSimilarityScore(i-1, j-1, charFreqRef, charFreqQuery);
                if (groupMatch) match = H[(k+1)%3][offsetDiag] + param.match;
                else match = H[(k+1)%3][offsetDiag] + param.mismatch;
            }
            
            if (offsetUp >= 0)
            {
                insOp = H[(k+2)%3][offsetUp] + param.gapOpen;
                insExt = I[(k+1)%2][offsetUp] + param.gapExtend;
            }

            if (offsetLeft >=0)
            {
                delOp = H[(k+2)%3][offsetLeft] + param.gapOpen;
                delExt = D[(k+1)%2][offsetLeft] + param.gapExtend;
            }

            I[k%2][offset] =  insOp;
            D[k%2][offset] =  delOp;
            
            if (insExt >= insOp) 
            {
                I[k%2][offset] = insExt;
                currentIState = STATE::II;
            }
            if (delExt >= delOp) 
            {
                D[k%2][offset] = delExt;
                currentDState = STATE::DD;
            }
            if (match > I[k%2][offset]) 
            {
                if (match > D[k%2][offset]) H[k%3][offset] = match;
                else 
                {
                    H[k%3][offset] = D[k%2][offset];
                    currentHState = STATE::HD;
                }
            }
            else if (I[k%2][offset] > D[k%2][offset]) 
            {
                H[k%3][offset] = I[k%2][offset];
                currentHState = STATE::HI;
            }
            else 
            {
                H[k%3][offset] = D[k%2][offset];
                currentHState = STATE::HD;
            }

            // if (j==0) std::cout << (int)currentHState << "-" << (int)currentIState << "-" << (int)currentDState << "-" << (int)updateState(currentHState, currentIState, currentDState) << std::endl;

            TB.push_back(updateState(currentHState, currentIState, currentDState));

            score = H[k%3][offset];
            state = currentHState;
        }
    }

    // Deallocate memory for scores
    for (size_t sIndx=0; sIndx<3; sIndx++) {
        std::free(H[sIndx]);
        if (sIndx < 2) {
            std::free(I[sIndx]);
            std::free(D[sIndx]);
        }
    }

    tracebackGrpToGrp (state, TB, wfLL, wfLen, alignment, ref, query);
    // printf("Score: %d\n", score);
    return;

}
*/


__device__ char charFreqRef [10000*5];
__device__ char charFreqQry [10000*5];
__device__ int32_t globalH [3*10000*2];
__device__ int32_t globalI [2*10000*2];
__device__ int32_t globalD [2*10000*2];
__device__ int32_t globalWfLL  [10000*2];
__device__ int32_t globalWfLen [10000*2];
__device__ int8_t globalTB [10000*10000];

__global__ void alignGrpToGrp_cuda(char *ref, char *qry, int16_t* param, char *alignment, int32_t* seqInfo)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    // int gs = gridDim.x;
    int idx = bx*bs+tx;

    const size_t threadNum = 512;
    const int sharedMemSizeL = 4096; // Less Size
    const int sharedMemSize = 8192;
    
    __shared__ int32_t sharedWfLL  [sharedMemSizeL];
    __shared__ int32_t sharedWfLen [sharedMemSizeL];
    __shared__ int8_t  sharedTB    [sharedMemSize];
    __shared__ size_t tbSize;
    __shared__ int8_t state;
    __shared__ size_t wfLLSize;
    __shared__ int32_t leftGrid;
    // __shared__ int32_t sharedH     [sharedMemSize];
    // __shared__ int32_t sharedD     [sharedMemSize];
    // __shared__ int32_t sharedI     [sharedMemSize];
// 
    // __syncthreads();
    if (bx == 0) {
        int32_t seqLen = seqInfo[0];
        int32_t refLen = seqInfo[1];
        int32_t qryLen = seqInfo[2];
        int32_t refNum = seqInfo[3];
        int32_t qryNum = seqInfo[4];

        int16_t p_match = param[0];
        int16_t p_mismatch = param[1];
        int16_t p_gapOpen = param[2];
        int16_t p_gapExtend = param[3];
        int16_t p_xdrop = param[4];
        int16_t p_marker = param[5];

        int maxWFLen = 0; //wavefront length
        if (p_marker != 0) maxWFLen = p_marker + 2;
        else maxWFLen = refLen + qryLen + 2;

        
        // size_t tbSize = 0;

        // Traceback
        // char* refAlign, *qryAlign;

        
        int tbIdx = 0;
        
        int maxSeqLen = 20000;
        int maxLen = refLen > qryLen ? refLen : qryLen;
        __syncthreads();
        if (idx < threadNum) {
            for (int i = idx; i < maxLen; i += threadNum) {
                for (int j = i*5; j < i*5+5; j++) {
                    charFreqRef[j] = 0;
                    charFreqQry[j] = 0;
                }
                for (int j = 0; j < refNum; j++) {
                    int k = seqLen*j;
                    if (i < refLen) {
                        if      (ref[k+i]=='A' || ref[k+i]=='a') charFreqRef[i*5]+=1;
                        else if (ref[k+i]=='C' || ref[k+i]=='c') charFreqRef[i*5+1]+=1;
                        else if (ref[k+i]=='G' || ref[k+i]=='g') charFreqRef[i*5+2]+=1;
                        else if (ref[k+i]=='T' || ref[k+i]=='t') charFreqRef[i*5+3]+=1;
                        else charFreqRef[i*5+4]+=1;
                    }
                    if (i < qryLen) {
                        if      (qry[k+i]=='A' || qry[k+i]=='a') charFreqQry[i*5]+=1;
                        else if (qry[k+i]=='C' || qry[k+i]=='c') charFreqQry[i*5+1]+=1;
                        else if (qry[k+i]=='G' || qry[k+i]=='g') charFreqQry[i*5+2]+=1;
                        else if (qry[k+i]=='T' || qry[k+i]=='t') charFreqQry[i*5+3]+=1;
                        else charFreqQry[i*5+4]+=1;
                    }
                }
            }
        }
        __syncthreads();

        if (idx < threadNum) {
            for (int i = idx; i < 3*maxSeqLen; i = i + threadNum) {
                globalH[i] = -INF;
                if (i < 2*maxSeqLen) {
                    globalD[i] = -INF;
                    globalI[i] = 0;
                }
            }
        }
        __syncthreads();
        int32_t L[3], U[3];
        for(size_t i=0; i<3; i++) {
            L[i]=0;
            U[i]=0;
        }
        __syncthreads();
        // if (tx == 0) printf("tbSize:%d, ref:%d, qry:%d\n",  maxWFLen, refLen, qryLen);
        if (tx == 0) wfLLSize = 0;
        for (int k = 0; k < refLen+qryLen+1; ++k) {
            L[k%3] = (k<=qryLen)?0:k-qryLen;
            U[k%3] = (k<=refLen)?k:refLen;
            // TODO: Copy shared to global
            // if (k > 0 && k%sharedMemSizeL == 0) {
            //     if (idx < threadNum) {
            //         for (int i = idx; i < sharedMemSizeL; i += threadNum) {
            //             int q = k/sharedMemSizeL;
            //             // printf("%d, %d, %d\n", i, q, sharedMemSize);
            //             globalWfLen[sharedMemSizeL*(q-1)+i] = sharedWfLen[i];
            //             globalWfLL[sharedMemSizeL*(q-1)+i] = sharedWfLL[i];
            //         }
            //     }
            // }
            __syncthreads();
            if (tx == 0) {
                // sharedWfLL[k%sharedMemSizeL] = (L[k%3]);
                // sharedWfLen[k%sharedMemSizeL] = U[k%3]-L[k%3]+1;
                globalWfLen[wfLLSize] = L[k%3];
                globalWfLL[wfLLSize] = U[k%3]-L[k%3]+1;
                wfLLSize += 1;
            }
            __syncthreads();
            
            int32_t numRound = (U[k%3]+1-L[k%3]) / threadNum + 1;
            if (tx == 0) leftGrid = U[k%3]+1-L[k%3];
            __syncthreads();
            for (int round = 0; round < numRound; ++round) {
                if (idx < min(leftGrid, (int32_t)threadNum)) {
                    // if (idx == 0) printf("round: %d, tbIdx: %d, idx:%d\n",round, tbIdx, idx);
                    int32_t i=L[k%3]+round*threadNum+idx;
                    int32_t j=(k-i); //j->Query Index
                    int32_t match = -INF, insOp = -INF, delOp = -INF, insExt = -INF, delExt = -INF;
                    int32_t offset = i-L[k%3];
                    int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1;
                    int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                    int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                    STATE currentHState = STATE::HH;
                    STATE currentIState = STATE::IH;
                    STATE currentDState = STATE::DH;

                    int32_t tempI, tempD, tempH;

                    if (k==0) match=0;

                    if (offsetDiag>=0 && j>0) {
                        float similarScore=0;
                        for (int l=0; l<5; l++) similarScore+= sqrtf(charFreqRef[5*(i-1)+l]*charFreqQry[5*(j-1)+l]);
                        bool groupMatch = (similarScore>0.95);
                        if (groupMatch) match = globalH[(k+1)%3*maxSeqLen+offsetDiag] + p_match;
                        else match = globalH[(k+1)%3*maxSeqLen+offsetDiag] + p_mismatch;
                    }

                    if (offsetUp >= 0) {
                        insOp = globalH[(k+2)%3*maxSeqLen+offsetUp] + p_gapOpen;
                        insExt = globalI[(k+1)%2*maxSeqLen+offsetUp] + p_gapExtend;
                    }

                    if (offsetLeft >=0) {
                        delOp = globalH[(k+2)%3*maxSeqLen+offsetLeft] + p_gapOpen;
                        delExt = globalD[(k+1)%2*maxSeqLen+offsetLeft] + p_gapExtend;
                    }

                    tempI = insOp;
                    tempD = delOp;

                    
                    if (insExt >= insOp) {
                        tempI = insExt;
                        currentIState = STATE::II;
                    }
                    if (delExt >= delOp) {
                        tempD = delExt;
                        currentDState = STATE::DD;
                    }
                    if (match > tempI) {
                        if (match > tempD) tempH = match;
                        else 
                        {
                            tempH = tempD;
                            currentHState = STATE::HD;
                        }
                    }
                    else if (tempI > tempD) {
                        tempH = tempI;
                        currentHState = STATE::HI;
                    }
                    else {
                        tempH = tempD;
                        currentHState = STATE::HD;
                    }
                    // tempI : I[(k%2)*maxSeqLen+offset]
                    // tempD : D[(k%2)*maxSeqLen+offset]
                    // tempH : globalH[(k%3)*maxSeqLen+offset]
                    // sharedH[tbIdx+idx] = tempH;
                    // sharedI[tbIdx+idx] = tempI;
                    // sharedD[tbIdx+idx] = tempD;
                    globalI[(k%2)*maxSeqLen+offset] = tempI;
                    globalD[(k%2)*maxSeqLen+offset] = tempD;
                    globalH[(k%3)*maxSeqLen+offset] = tempH;
                    // sharedTB[tbIdx+idx] = updateState_cuda(currentHState, currentIState, currentDState);
                    globalTB[tbSize+idx] =  updateState_cuda(currentHState, currentIState, currentDState);
                    // score = tempH;
                    // state = currentHState;
                    // globalH[0] = tempH;
                    // if (tx == min(leftGrid, (int32_t)threadNum)) state = currentHState;
                    if (idx == min((int)leftGrid, (int)threadNum)) {
                        state = currentHState;
                    }
                }
                // __syncthreads();
                // if (tx < min(leftGrid, (int32_t)threadNum)) {
                //     printf("Writing to %d\n", tbIdx + tx);
                //     globalTB[tbIdx + tx] = sharedTB[tx];
                // }
                __syncthreads();
                if (tx == 0) {
                    tbIdx += min((size_t)leftGrid, (size_t)threadNum);
                    tbSize += min((size_t)leftGrid, (size_t)threadNum);
                    leftGrid -= min((int)leftGrid, (int)threadNum);
                }
                // __syncthreads();
                // if (tbIdx > (sharedMemSize/2)) {
                //     if (idx < threadNum) {
                //         for (int i = idx; i < sharedMemSize/2; i += threadNum) {
                //             int q = k/(sharedMemSize/2);
                //             // printf("%d, %d, %d\n", i, q, sharedMemSize);
                //             // globalTB[sharedMemSize/2*(q-1)+i] = sharedTB[i];
                //             // globalH [(k%3)*maxSeqLen+sharedMemSize/2*(q-1)+i] = sharedWfLL[i];
                //             // globalD [(k%2)*maxSeqLen+sharedMemSize/2*(q-1)+i] = sharedWfLL[i];
                //             // globalI [(k%2)*maxSeqLen+sharedMemSize/2*(q-1)+i] = sharedWfLL[i];
                //             // sharedTB[i] = sharedTB[sharedMemSize/2+i];
                //             // sharedH[i] = sharedH[sharedMemSize/2+i];
                //             // sharedD[i] = sharedD[sharedMemSize/2+i];
                //             // sharedI[i] = sharedI[sharedMemSize/2+i];
                //         }
                //     }
                //     if (idx == 0) tbIdx -= (sharedMemSize/2);
                // }
                __syncthreads();
                
            }
            __syncthreads();

        }
        __syncthreads();
    
        // Traceback
 
        
        // std::vector<std::string> refAlign, queryAlign;
        // refAlign.resize(ref.size()); queryAlign.resize(query.size());
        // if (tx == 0) printf("wSize: %d, tbSize:%d, ref:%d, qry:%d\n", wfLLSize, tbSize, refLen, qryLen);
        if (tx == 0) {
            int32_t tb = wfLLSize - 1;
            int32_t tbIndex = tbSize - 1; 
             
            int32_t refIndex = refLen-1;
            int32_t qryIndex = qryLen-1;

            int8_t dir;
            size_t alnIdx = 0;
            int8_t currentTB;

            // printf("wSize: %d, tbSize:%d, ref:%d, qry:%d\n", wfLLSize, tbSize, refLen, qryLen);

            while (refIndex>=0 && qryIndex>=0) {
                // printf("%d, %d, %d, %d\n", alnIdx, tbIndex, refIndex, qryIndex);
                currentTB = globalTB[tbIndex];
                if (state == 0) { 
                    // Current State M
                    state = currentTB & 0x03;
                    if (state == 0) dir = 0;
                    else if (state == 1)
                    {
                        dir = 1;
                        if (currentTB & 0x04) state = 1;
                        else state = 0;   
                    }
                    else
                    {
                        dir = 2;
                        if (currentTB & 0x08) state = 2;
                        else state = 0;
                    }
                }
                else if (state == 1) { 
                    // Current State I
                    dir = 1;
                    if (currentTB & 0x04) state = 1;
                    else state = 0;
                } 
                else { 
                    // Current State D
                    dir = 2;
                    if (currentTB & 0x08) state = 2;
                    else state = 0;
                }

                tbIndex -= (refIndex - globalWfLL[tb] + 1 + globalWfLen[tb-1]);
                if (dir == 0) {
                    tbIndex -= globalWfLen[tb-2]; tbIndex += refIndex - globalWfLL[tb-2];
                    for (size_t i=0; i<refNum; i++) alignment[i*seqLen+alnIdx] = ref[i*seqLen+refIndex]; 
                    for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
                    tb--;tb--;qryIndex--;refIndex--;
                }
                else if (dir == 1) 
                {
                    tbIndex += (refIndex - globalWfLL[tb-1] + 1);
                    for (size_t i=0; i<refNum; i++) alignment[i*seqLen+alnIdx] = '-'; 
                    for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
                    tb--;qryIndex--;
                }
                else 
                {
                    tbIndex += (refIndex - globalWfLL[tb-1]); 
                    for (size_t i=0; i<refNum; i++) alignment[i*seqLen+alnIdx] = ref[i*seqLen+refIndex];  
                    for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+alnIdx] = '-';  
                    tb--;refIndex--;
                }
                alnIdx++;
            }

            while (refIndex>=0) {
                for (size_t i=0; i<refNum; i++) alignment[i*seqLen+alnIdx] = ref[i*seqLen+refIndex];  
                for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+alnIdx] = '-';  
                tb--;refIndex--;alnIdx++;
            }

            while (qryIndex>=0) {

                for (size_t i=0; i<refNum; i++) alignment[i*seqLen+alnIdx] = '-'; 
                for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
                tb--;qryIndex--;alnIdx++;
            }
        
            __syncthreads();
        }
        
   

    }
    
    return;
}

/*
__global__ void alignGrpToGrp_cuda1(char *ref, char *qry, int16_t* param, char *alignment, int32_t* seqInfo)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    // int gs = gridDim.x;
    int idx = bx*bs+tx;

    int32_t seqLen = seqInfo[0];
    int32_t refLen = seqInfo[1];
    int32_t qryLen = seqInfo[2];
    int32_t refNum = seqInfo[3];
    int32_t qryNum = seqInfo[4];

    int16_t p_match = param[0];
    int16_t p_mismatch = param[1];
    int16_t p_gapOpen = param[2];
    int16_t p_gapExtend = param[3];
    int16_t p_xdrop = param[4];
    int16_t p_marker = param[5];

    int maxWFLen = 0; //wavefront length
    if (p_marker != 0) maxWFLen = p_marker + 2;
    else maxWFLen = refLen + qryLen + 2;

    // float* charFreqRef; 
    // float* charFreqQry;

    int score = 0;
    // int32_t *H_0, *H_1, *H_2, *I_0, *I_1, *D_0, *D_1;
    // __shared__ int32_t L[3], U[3];
    // int32_t *wfLL, *wfLen;
    int wfIdx = 0;
    int8_t state=0;
    // int8_t* TB;
    int tbIdx = 0;
    
    int globalMemSize = 10000;

    size_t batchSize = 512;

    int32_t refIndex = refLen-1;
    int32_t qryIndex = qryLen-1;
    int32_t k = wfIdx-1;
    int32_t tbIndex = tbIdx-1;
    int8_t currentTB;
    int8_t dir;
    size_t alnIdx = 0;
    char *refAlign, *qryAlign;

    const int memLength = 10000;
    
        
    
    if (bx == 0) {
        // if (ref.size() < 1 || query.size() < 1) {fprintf(stderr, "Error: Number of Ref/Query <= 0\n"); exit(1);}
        
        // printf("refLen %d, qryLen %d, refNum %d, qryNum %d", refLen, qryLen, refNum, qryNum);
        // if (refLen<=0 || queryLen<=0) {fprintf(stderr, "Error: Ref/Query length <= 0\n"); exit(1);}
        // if (tx == 0) {
        //    charFreqRef = new float [refLen*5];
        //    charFreqQry = new float [qryLen*5];
        // }
        
        if (tx < batchSize) {
            for (int i = idx; i < max(refLen, qryLen); i += batchSize) {
                for (int j = i*5; j < i*5+5; j++) {
                    charFreqRef[j] = 0;
                    charFreqQry[j] = 0;
                }
                for (int j = 0; j < refNum; j++) {
                    int k = seqLen*j;
                    if (i < refLen) {
                        if      (ref[k+i]=='A' || ref[k+i]=='a') charFreqRef[i*5]+=1;
                        else if (ref[k+i]=='C' || ref[k+i]=='c') charFreqRef[i*5+1]+=1;
                        else if (ref[k+i]=='G' || ref[k+i]=='g') charFreqRef[i*5+2]+=1;
                        else if (ref[k+i]=='T' || ref[k+i]=='t') charFreqRef[i*5+3]+=1;
                        else charFreqRef[i*5+4]+=1;
                    }
                    if (i < qryLen) {
                        if      (qry[k+i]=='A' || qry[k+i]=='a') charFreqQry[i*5]+=1;
                        else if (qry[k+i]=='C' || qry[k+i]=='c') charFreqQry[i*5+1]+=1;
                        else if (qry[k+i]=='G' || qry[k+i]=='g') charFreqQry[i*5+2]+=1;
                        else if (qry[k+i]=='T' || qry[k+i]=='t') charFreqQry[i*5+3]+=1;
                        else charFreqQry[i*5+4]+=1;
                    }
                }
            }
        }

        __syncthreads();

    
        
        const size_t localMemSize = 5000;
        const size_t exchangeMemSize = localMemSize * 6 / 10;
        // if (tx == 0) printf("Stage 1\n");
        if (tx == 0) {
            int32_t H     [3*localMemSize];
            int32_t D     [2*localMemSize];
            int32_t I     [2*localMemSize];
            for (int i = 0; i <3*localMemSize; i += 1) {
                H[i] = -INF;
                if ( i < 2*localMemSize) {
                    D[i] = -INF;
                    I[i] = 0;
                }
            }
            int32_t L[3], U[3];
            for(size_t i=0; i<3; i++) {
                L[i]=0;
                U[i]=0;
            }

            int memOffset = 0;

            for (int32_t k=0; k<refLen+qryLen+1; k++) {

                // Memory exchange between Local and Global
                if (k / exchangeMemSize > 0 & k % exchangeMemSize == 1000) {
                    // Testing function of Mem exchange
                    // for (int z=0; z<30; z++) printf("%d,",H[3*exchangeMemSize+z]);
                    // printf("\n%d\n", k);

                    for (int w = 0; w < exchangeMemSize; w++) {
                        for (int y = 0; y < 3; ++y) {
                            globalH[3*(w+memOffset)+y] = H[3*w+y];
                            if (y < 2) {
                                globalI[2*(w+memOffset)+y] = I[2*w+y];
                                globalD[2*(w+memOffset)+y] = D[2*w+y];
                            }
                        }
                        if (w < localMemSize - exchangeMemSize) {
                            for (int y = 0; y < 3; ++y) {
                                H[3*w+y] = H[3*(w+exchangeMemSize)+y];
                                if (y < 2) {
                                    I[2*w+y] = I[2*(w+exchangeMemSize)+y];
                                    D[2*w+y] = D[2*(w+exchangeMemSize)+y];
                                }
                            }
                        }
                    }
                    // for (int z=0; z<30; z++) printf("%d,",H[z]);
                    // printf("\n");
                    memOffset += exchangeMemSize;
                }
                // For testing 
                // for (int z=0; z <3; z++) H[3*(k-memOffset)+z] = k;
                
                // printf("%d\n", k);
                L[k%3] = (k<=qryLen)?0:k-qryLen;
                U[k%3] = (k<=refLen)?k:refLen;
                // wfLL[k] = (L[k%3]);
                // wfLen[k] = (U[k%3]-L[k%3]+1);
                ++wfIdx; 
                
                for(int32_t i=L[k%3]; i<U[k%3]+1; i++) { // i->Ref Index 
                    int32_t j=(k-i); //j->Query Index
                    int32_t match = -INF, insOp = -INF, delOp = -INF, insExt = -INF, delExt = -INF;
                    int32_t offset = i-L[k%3];
                    int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1;
                    int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                    int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                    STATE currentHState = STATE::HH;
                    STATE currentIState = STATE::IH;
                    STATE currentDState = STATE::DH;

                    if (k==0) match=0;
                    
                    if (offsetDiag>=0 && j>0) {
                        // std::cout << i << "-" << j << "-" << offsetDiag << "\n";
                        float score=0;
                        for (int l=0; l<5; l++) score+= sqrtf(charFreqRef[5*(i-1)+l]*charFreqQry[5*(j-1)+l]);
                        bool groupMatch = (score>0.95); 
                        //bool groupMatch = getSimilarityScore(i-1, j-1, charFreqRef, charFreqQuery);
                        if (groupMatch) {
                            match = H[3*(offsetDiag-memOffset) + ((k+1)%3)] + p_match;
                        }
                        else {
                            match = H[3*(offsetDiag-memOffset) + ((k+1)%3)] + p_mismatch;
                        }
                    }
                    
                    if (offsetUp >= 0) {
                        insOp =  H[3*(offsetUp-memOffset)+((k+2)%3)] + p_gapOpen;
                        insExt = I[2*(offsetUp-memOffset)+((k+1)%2)] + p_gapExtend;
                    }
                    
                    if (offsetLeft >=0) {
                        delOp =  H[3*(offsetLeft-memOffset)+((k+2)%3)] + p_gapOpen;
                        delExt = I[2*(offsetLeft-memOffset)+((k+1)%2)] + p_gapExtend;
                    }
                    
                    // printf("Stage 1\n");
                    I[2*(offset-memOffset)+(k%2)] = insOp;
                    D[2*(offset-memOffset)+(k%2)] = delOp;
                    
                    /*
                    if (insExt >= insOp) 
                    {
                        I[2*(offset-memOffset)+(k%2)] = insExt;
                        currentIState = STATE::II;
                    }
                    if (delExt >= delOp) 
                    {
                        D[2*(offset-memOffset)+(k%2)] = delExt;
                        currentDState = STATE::DD;
                    }
                    if (match > I[2*(offset-memOffset)+(k%2)]) 
                    {
                        if (match > D[2*(offset-memOffset)+(k%2)]) H[3*(offset-memOffset)+(k%3)] = match;
                        else 
                        {
                            H[3*(offset-memOffset)+(k%3)] = D[2*(offset-memOffset)+(k%2)];
                            currentHState = STATE::HD;
                        }
                    }
                    else if (I[2*(offset-memOffset)+(k%2)] > D[2*(offset-memOffset)+(k%2)]) 
                    {
                        H[3*(offset-memOffset)+(k%3)] = I[2*(offset-memOffset)+(k%2)];
                        currentHState = STATE::HI;
                    }
                    else 
                    {
                        H[3*(offset-memOffset)+(k%3)] = D[2*(offset-memOffset)+(k%2)];
                        currentHState = STATE::HD;
                    }
                    */


                    // if (j==0) std::cout << (int)currentHState << "-" << (int)currentIState << "-" << (int)currentDState << "-" << (int)updateState(currentHState, currentIState, currentDState) << std::endl;

                    // TB[tbIdx] = updateState_cuda(currentHState, currentIState, currentDState);
                    // ++tbIdx;
                    // score = H[3*offset+(k%3)];
                    // // score = H[k%3][offset];
                    // state = currentHState;
                    /*
                }
                
            
            }
        }
        
        // __syncthreads();
        
        // Deallocate memory for scores
        // cudaFree(H_0);
        // cudaFree(H_1);
        // cudaFree(H_2);
        // cudaFree(D_0);
        // cudaFree(D_1);
        // cudaFree(I_0);
        // cudaFree(I_1);
        

        /*
        // Traceback
        
        if (tx == 0) {
            refAlign = new char [refLen * refNum];
            qryAlign = new char [qryLen * qryNum];
        }
        
        // std::vector<std::string> refAlign, queryAlign;
        // refAlign.resize(ref.size()); queryAlign.resize(query.size());
        if (tx == 0) {
        while (refIndex>=0 && qryIndex>=0) {
            currentTB = TB[tbIndex];
            if (state == 0) { 
                // Current State M
                state = currentTB & 0x03;
                if (state == 0) dir = 0;
                else if (state == 1)
                {
                    dir = 1;
                    if (currentTB & 0x04) state = 1;
                    else state = 0;   
                }
                else
                {
                    dir = 2;
                    if (currentTB & 0x08) state = 2;
                    else state = 0;
                }
            }
            else if (state == 1) { 
                // Current State I
                dir = 1;
                if (currentTB & 0x04) state = 1;
                else state = 0;
            } 
            else { 
                // Current State D
                dir = 2;
                if (currentTB & 0x08) state = 2;
                else state = 0;
            }

            tbIndex -= (refIndex - wfLL[k] + 1 + wfLen[k-1]);
            if (dir == 0) {
                tbIndex -= wfLen[k-2]; tbIndex += refIndex - wfLL[k-2];
                for (size_t i=0; i<refNum; i++) refAlign[i*seqLen+alnIdx] = ref[i*seqLen+refIndex]; 
                for (size_t i=0; i<qryNum; i++) qryAlign[i*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
                k--;k--;qryIndex--;refIndex--;alnIdx++;
            }
            else if (dir == 1) 
            {
                tbIndex += (refIndex - wfLL[k-1] + 1);
                for (size_t i=0; i<refNum; i++) refAlign[i*seqLen+alnIdx] = '-'; 
                for (size_t i=0; i<qryNum; i++) qryAlign[i*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
                k--;qryIndex--;alnIdx++;
            }
            else 
            {
                tbIndex += (refIndex - wfLL[k-1]); 
                for (size_t i=0; i<refNum; i++) refAlign[i*seqLen+alnIdx] = ref[i*seqLen+refIndex];  
                for (size_t i=0; i<qryNum; i++) qryAlign[i*seqLen+alnIdx] = '-';  
                k--;refIndex--;alnIdx++;
            }
            __syncthreads();

        }

        while (refIndex>=0) {
            for (size_t i=0; i<refNum; i++) refAlign[i*seqLen+alnIdx] = ref[i*seqLen+refIndex];  
            for (size_t i=0; i<qryNum; i++) qryAlign[i*seqLen+alnIdx] = '-';  
            k--;refIndex--;alnIdx++;
            __syncthreads();
        }
        
        while (qryIndex>=0) {

            for (size_t i=0; i<refNum; i++) refAlign[i*seqLen+alnIdx] = '-'; 
            for (size_t i=0; i<qryNum; i++) qryAlign[i*seqLen+alnIdx] = qry[i*seqLen+qryIndex]; 
            k--;qryIndex--;alnIdx++;
            __syncthreads();
        }

        // tracebackGrpToGrp_cuda(state, TB, wfLL, wfLen, alignment, ref, qry);
        printf("Score: %d\n", score);
        // TODO: reverse 
        // for (auto &s: refAlign) s = std::string(s.rbegin(), s.rend());
        // for (auto &s: queryAlign) s = std::string(s.rbegin(), s.rend());

        for (size_t i = 0; i < refNum; ++i) {
            for (size_t j = 0; j < refLen; ++j) {
                alignment[i*seqLen+j] = refAlign[i*seqLen+j]; 
            }
        }
        for (size_t i = 0; i < qryNum; ++i) {
            for (size_t j = 0; j < qryLen; ++j) {
                alignment[(i+refNum)*seqLen+j] = qryAlign[i*seqLen+j]; 
            }
        }
        
        __syncthreads();
        delete refAlign;
        delete qryAlign;
        delete wfLL;
        delete wfLen;
        // cudaFree(refAlign);
        // cudaFree(qryAlign);
        // cudaFree(wfLL);
        // cudaFree(wfLen);
        // alignment = std::make_pair(refAlign, queryAlign);
        }
        
    }
    


    
    // if (bx == 0 && tx == 0) {
    //     delete charFreqRef;
    //     delete charFreqQry;
    // }
    
    return;
}
*/

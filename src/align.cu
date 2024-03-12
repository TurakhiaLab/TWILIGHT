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


__device__ float charFreqRef [40000*5];
__device__ float charFreqQry [40000*5];
__device__ int32_t globalH [3*40000*2];
__device__ int32_t globalI [2*40000*2];
__device__ int32_t globalD [2*40000*2];
__device__ int32_t globalWfLL  [40000*2];
__device__ int32_t globalWfLen [40000*2];
__device__ int8_t globalTB [40000*40000];

__global__ void alignGrpToGrp_cuda(
    char *ref, 
    char *qry,
    int16_t* param, 
    char *alignment, 
    int32_t* seqInfo
    // int8_t*  globalTB,
    // float*   charFreqRef,
    // float*   charFreqQry,
    // int32_t* globalH,
    // int32_t* globalD,
    // int32_t* globalI,
    // int32_t* globalWfLL,
    // int32_t* globalWfLen
    )
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
                }
                for (int j = 0; j < qryNum; j++) {
                    int k = seqLen*j;
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
                globalWfLL[wfLLSize] = L[k%3];
                globalWfLen[wfLLSize] = U[k%3]-L[k%3]+1;
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
        // if (tx == 0) {
        //     for (int i = 1100; i < 1200; ++i) {
        //         printf("%d,", globalWfLL[i]);
        //         if (i%10 == 9) printf("\n");
        //     }
        // }
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

            
            while (refIndex>=0 && qryIndex>=0) {
                // printf("%d, %d, %d, %d, %d\n", tbIndex, tb, refIndex, qryIndex, alnIdx);
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

__device__ int32_t sL[3];
__device__ int32_t sU[3];
__device__ int16_t S  [3*128];
__device__ int32_t max_score_list [128]; 

__global__ void alignGrpToGrp_talco(char *ref, char *qry, int16_t* param, char *alignment, int32_t* seqInfo)
{
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int bs = blockDim.x;
    // int gs = gridDim.x;
    int tidx = bx*bs+tx;

    const size_t threadNum = 128;
    const int fLen = 128; // frontier length (assuming anti-diagonal length cannot exceed 1024)

    __shared__ bool last_tile;
    __shared__ int8_t tile_aln[2*fLen];
    // __shared__ int16_t S  [3*fLen];
    __shared__ int32_t CS [3*fLen];
    __shared__ int16_t I  [2*fLen];
    __shared__ int32_t CI [2*fLen];
    __shared__ int16_t D  [2*fLen];
    __shared__ int32_t CD [2*fLen];
    __shared__ int8_t tb  [fLen*fLen];
    __shared__ float charFreqRef [fLen*5];
    __shared__ float charFreqQry [fLen*5];
    // __shared__ int32_t ftr_addr;
    __shared__ int16_t idx [7]; 
    // __shared__ int32_t sL[3];
    // __shared__ int32_t sU[3];

    // __shared__ int32_t max_score_list [fLen]; 
    // __shared__ int32_t max_score_prime; 
    __shared__ int32_t max_score_ref [fLen]; 
    __shared__ int32_t max_score_query [fLen];
    __shared__ int32_t max_score_start_addr; 
    __shared__ int32_t max_score_start_ftr;
    __shared__ int32_t max_score_ref_idx;    
    __shared__ int32_t max_score_query_idx;
    __shared__ int32_t max_score;

    __shared__ bool converged; 
    __shared__ bool conv_logic;

    
    // [0] ftr_length_idx
    // [1] ftr_lower_limit_idx
    // [2] tb_idx
    // [3] reference_idx
    // [4] query_idx
    // [5] globalAlnIdx
    // [6] tile
    __syncthreads();
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

        int16_t tile = 0;
        int32_t L[3], U[3];

        int32_t ftr_length [2*fLen];
        int32_t ftr_lower_limit [2*fLen];
        // initialize global values
        if (tidx == 0) {
            last_tile = false;
            for (int i = 0; i < 7; ++i) {
                idx[i] = 0;
            }
        }
        
        __syncthreads();


        while (!last_tile) {
            int32_t inf = p_xdrop + 1;
            int16_t reference_idx = idx[3];
            int16_t query_idx = idx[4];
            int32_t reference_length = refLen - reference_idx; 
            int32_t query_length = qryLen - query_idx;
            int32_t score = 0; 
            int32_t ftr_addr = 0;
            int16_t ftr_length_idx = 0;
            int16_t ftr_lower_limit_idx = 0;
            int32_t max_score_prime = -(p_xdrop+1);
            int16_t tb_idx = 0;
            for (size_t sIndx=0; sIndx<3; sIndx++) {
                L[sIndx] = sIndx;
                U[sIndx] = -sIndx;
            }
            for (size_t i = 0; i < 2*fLen; i++) {
                ftr_length[i] = 0;
                ftr_lower_limit[i] = 0;
            }

            if (tx == 0) {
                max_score_start_addr = 0; 
                max_score_start_ftr = 0;
                max_score_ref_idx = 0;    
                max_score_query_idx = 0;
                max_score = 0;
                converged = false;
                conv_logic = false;
                // ftr_length_idx = 0;
                // ftr_lower_limit_idx = 0;
                // printf("marker: %d, xdrop: %d\n", p_marker, p_xdrop);
            }
            __syncthreads();
            int32_t conv_score = 0; 
            int32_t conv_value = 0; 
            
            int8_t tb_state = 0;
            int8_t ptr = 0;  // Main pointer
            int8_t state = 0;
            int16_t aln_idx = 0;

            
            int32_t prev_conv_s = -1;

            bool Iptr = false;
            bool Dptr = false; // I and D states; xx00 -> M, x001 -> I (M) open, x101 -> I (I) extend, 0x10 -> D (M) open, 1x10 -> D (D) extend
            // if (tidx == 0) printf("Tile: %d, refIdx: %d, qryIdx: %d\n", tile, reference_idx, query_idx);
            __syncthreads();
            if (tidx == 0) {
                for (size_t sIndx=0; sIndx<3; sIndx++) {
                    sL[sIndx] = sIndx;
                    sU[sIndx] = -sIndx;
                }
            }
            __syncthreads();
            
            
            // Initialize shared memory
            if (tidx < threadNum) {
            // if (tidx == 0) {
                for (int i = 0; i < 5*fLen; i = i+threadNum) {
                    charFreqRef[i] = 0;
                    charFreqQry[i] = 0;
                    if (i < 3*fLen) {                        
                        S[i] = -1;
                        CS[i] = -1;
                    }
                    if (i < 2*fLen) {
                        I[i] = -1;
                        CI[i] = -1;
                        D[i] = -1;
                        CD[i] = -1;
                        tile_aln[i] = -1;
                        // ftr_length[i] = 0;
                        // ftr_lower_limit[i] = 0;
                    }
                    if (i < fLen) {
                        max_score_list [i] = -(p_xdrop+1); 
                        max_score_ref [i] = 0; 
                        max_score_query [i] = 0;
                    }
                }
            }
            // __syncthreads();
            // if (tidx == 0) printf("DEBUG: break 1\n");
            __syncthreads();
            for (int32_t k = 0; k < reference_length + query_length - 1; k++){
                if (tidx < fLen) {
                    max_score_list [tidx] = -(p_xdrop+1); 
                    max_score_ref [tidx] = 0; 
                    max_score_query [tidx] = 0;
                }
                __syncthreads();
                // if (tidx == 0) printf("k: %d\n", k);
                // if (tidx == 0) printf("ref+qry: %d\n", reference_length + query_length - 1);
                if (L[k%3] >= U[k%3]+1) { // No more cells to compute based on x-drop critieria
                    // std::cout << "No more cells to compute based on x-drop critieria" << std::endl;
                    // if (tidx == 0) printf("No more cells to compute based on x-drop critieria\n");
                    if (tidx == 0) printf("No more cells to compute based on x-drop critieria (L, U) = (%d, %d)\n", sL[k%3], sU[k%3]+1);
                    __syncthreads();
                    break;
                }
                // if (U[k%3]-L[k%3]+1 > fLen) { // Limit the size of the anti-diagonal
                //     fprintf(stderr, "ERROR: anti-diagonal larger than the max limit!\n");
                //     exit(1);
                // }
                // __syncthreads();
                // if (tidx == 0) printf("DEBUG: break 1.2\n");
                __syncthreads();
                // if (tidx == 0) {
                if (k <= p_marker) {
                    // printf("ftr length: %d ftr_addr: %d  ftr lower limit: %d idx0: %d, idx1: %d\n", U[k%3] - L[k%3] + 1, ftr_addr, L[k%3], ftr_length_idx, ftr_lower_limit_idx);
                    // printf("Before idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    ftr_length[ftr_length_idx] = (U[k%3] - L[k%3] + 1);
                    ftr_lower_limit[ftr_lower_limit_idx] = L[k%3];
                    // printf("before ftr_addr: %d\n", ftr_addr);
                    // printf("UL %d\n", sU[k%3] - sL[k%3] + 1);
                    ftr_addr += U[k%3] - L[k%3] + 1;
                    // printf("after ftr_addr: %d\n", ftr_addr);
                    ftr_length_idx += 1;
                    ftr_lower_limit_idx += 1;
                    // printf("After idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    // std::cout << "ftr length: " << U[k%3] - L[k%3] + 1 << " ftr_addr: " << ftr_addr << " ftr lower limit: " << L[k%3] << " tb len: " << tb.size() << std::endl;
                    
                }
                // }   
                __syncthreads();
                // if (tidx == 0) printf("DEBUG: break 1.4\n");
                // Calculate character frequency
                if (tidx < threadNum) {
                    for (int i = tidx; i < fLen; i += threadNum) {
                        for (int j = i*5; j < i*5+5; j++) {
                            charFreqRef[j] = 0;
                            charFreqQry[j] = 0;
                        }
                        for (int j = 0; j < refNum; j++) {
                            int k = seqLen*j;
                            if      (ref[k+reference_idx+i]=='A' || ref[k+reference_idx+i]=='a') charFreqRef[i*5]+=1;
                            else if (ref[k+reference_idx+i]=='C' || ref[k+reference_idx+i]=='c') charFreqRef[i*5+1]+=1;
                            else if (ref[k+reference_idx+i]=='G' || ref[k+reference_idx+i]=='g') charFreqRef[i*5+2]+=1;
                            else if (ref[k+reference_idx+i]=='T' || ref[k+reference_idx+i]=='t') charFreqRef[i*5+3]+=1;
                            else charFreqRef[i*5+4]+=1;
                        }
                        for (int j = 0; j < qryNum; j++) {
                            int k = seqLen*j;
                            if      (qry[k+query_idx+i]=='A' || qry[k+query_idx+i]=='a') charFreqQry[i*5]+=1;
                            else if (qry[k+query_idx+i]=='C' || qry[k+query_idx+i]=='c') charFreqQry[i*5+1]+=1;
                            else if (qry[k+query_idx+i]=='G' || qry[k+query_idx+i]=='g') charFreqQry[i*5+2]+=1;
                            else if (qry[k+query_idx+i]=='T' || qry[k+query_idx+i]=='t') charFreqQry[i*5+3]+=1;
                            else charFreqQry[i*5+4]+=1;
                        }
                    }
                }
                __syncthreads();
                // if (tidx == 0) printf("DEBUG1: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                int32_t numRound = (U[k%3]+1-L[k%3]) / threadNum + 1;
                int32_t leftGrid = U[k%3]+1-L[k%3];
                __syncthreads();
                // if (tidx == 0) printf("k: %d, Max Score: %d\n", k, max_score);
                // printf("%d, max: %d\n", tidx, max_score_list[tidx]);
                for (int round = 0; round < numRound; ++round) {
                    if (tidx < min(leftGrid, (int32_t)threadNum)) {
                        int16_t i=L[k%3]+round*threadNum+tidx;
                        int16_t Lprime = max(0, static_cast<int16_t>(k)-static_cast<int16_t>(reference_length) + 1);
                        int16_t j= min(static_cast<int16_t>(k), static_cast<int16_t>(reference_length - 1)) - (i-Lprime);; //j->Query Index
                        if (j < 0) if (tx == 0) { printf("ERROR: j less than 0.\n");}
                        int32_t match = -inf, insOp = -inf, delOp = -inf, insExt = -inf, delExt = -inf;
                        int32_t offset = i-L[k%3];
                        int32_t offsetDiag = L[k%3]-L[(k+1)%3]+offset-1; // L[0] - L[1] + 0 - 1
                        int32_t offsetUp = L[k%3]-L[(k+2)%3]+offset;
                        int32_t offsetLeft = L[k%3]-L[(k+2)%3]+offset-1;
                        int score_from_prev_tile = 0;
                        // if (tidx == 0) printf("O: %d, OD: %d OU: %d, OL: %d\n", offset, offsetDiag, offsetUp, offsetLeft);
                        if ((k==0) || ((offsetDiag >= 0) && (offsetDiag <= U[(k+1)%3]-L[(k+1)%3]))) {
                            if (k==0 && tile>0) {
                                score_from_prev_tile = tile*10;
                            }
                            float similarScore=0;
                            for (int l=0; l<5; l++) similarScore+= sqrtf(charFreqRef[5*(i)+l]*charFreqQry[5*(j)+l]);
                            bool groupMatch = (similarScore>0.95);
                            if (groupMatch) {
                                if (offsetDiag < 0) match = p_match + score_from_prev_tile; 
                                else match = S[(k+1)%3*fLen+offsetDiag] + p_match + score_from_prev_tile; 
                            }
                            else {
                                if (offsetDiag < 0) match = p_mismatch + score_from_prev_tile;
                                else match = S[(k+1)%3*fLen+offsetDiag] + p_mismatch + score_from_prev_tile;
                            }
                            // if (tidx == 0) printf("match: %d\n", match);
                        }
                        if ((offsetUp >= 0) && (offsetUp <= U[(k+2)%3]-L[(k+2)%3])) {
                            // insOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            // insExt = I[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                            delOp =  S[(k+2)%3*fLen+offsetUp] + p_gapOpen;
                            delExt = D[(k+1)%2*fLen+offsetUp] + p_gapExtend;
                            // if (tidx == 0) printf("del: %d, %d\n", delOp, delExt);
                        }
                        if ((offsetLeft >= 0) && (offsetLeft <= U[(k+2)%3]-L[(k+2)%3])) {
                            // delOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            // delExt = D[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                            insOp =  S[(k+2)%3*fLen+offsetLeft] + p_gapOpen;
                            insExt = I[(k+1)%2*fLen+offsetLeft] + p_gapExtend;
                        }
        
                        int16_t tempI, tempD, tempH;
                        tempI =  insOp;
                        tempD =  delOp;
                        Iptr = false;
                        Dptr = false;
                        // if (tile == 4 && k == 68) printf("id: %d, tempI: %d, tempD: %d, match: %d\n", tidx, tempI, tempD, match);
                        if (insExt >= insOp) {
                            tempI = insExt;
                            Iptr = true;
                        }
                        if (delExt >= delOp) {
                            tempD = delExt;
                            Dptr = true;
                        }
                        if (match > tempI) {
                            if (match > tempD) {
                                tempH = match;
                                ptr = 0;
                            }
                            else {
                                tempH = tempD;
                                ptr = 2;
                            }
                        }
                        else if (tempI > tempD) {
                            tempH = tempI;
                            ptr = 1;
                        }
                        else {
                            tempH = tempD;
                            ptr = 2;
                        }
                        // printf("%d: m: %d, i: %d, d:%d, h: %d, max: %d\n",tidx, match, tempI, tempD, tempH, max_score);
                        if (tempH < max_score-p_xdrop) {
                            tempH = -inf;
                        }
                        // printf("%d: m: %d, i: %d, d:%d, h: %d\n",tidx, match, tempI, tempD, tempH);
                        D[(k%2)*fLen+offset] = tempD; 
                        I[(k%2)*fLen+offset] = tempI;
                        S[(k%3)*fLen+offset] = tempH;

                        score = tempH;

                        max_score_query[tidx] = i;
                        max_score_ref[tidx] = j;
                        max_score_list[tidx] = score;

                        // printf("No.%d (r,q)=(%d,%d) pre_max:%d, score: %d\n", tidx, i, j, max_score_prime, score);
                        
                        if (k == p_marker - 1) { // Convergence algorithm
                            CS[(k%3)*fLen+offset] = (3 << 16) | (i & 0xFFFF); 
                            // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  tempCS << "\n";
                        } 
                        else if (k == p_marker) {
                            CS[(k%3)*fLen+offset] = (0 << 16) | (i & 0xFFFF);  // to extract value use (CS[k%3][offset] & 0xFFFF)
                            CI[(k%2)*fLen+offset] = (1 << 16) | (i & 0xFFFF);
                            CD[(k%2)*fLen+offset] = (2 << 16) | (i & 0xFFFF);
                            // if(DEBUG) std::cout << "Convergence Unique Id's: " <<  tempCS <<  " " << tempCI <<  " " << tempCD << "\n";
                        } 
                        else if (k >= p_marker + 1){
                            if (Iptr) {
                                CI[(k%2)*fLen+offset] = CI[((k+1)%2)*fLen+offsetLeft]; 
                            } else {
                                CI[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetLeft]; 
                            }

                            if (Dptr) {
                                CD[(k%2)*fLen+offset] = CD[((k+1)%2)*fLen+offsetUp];
                            } else {
                                CD[(k%2)*fLen+offset] = CS[((k+2)%3)*fLen+offsetUp];
                            }

                            if (ptr == 0) {
                                CS[(k%3)*fLen+offset] = CS[((k+1)%3)*fLen+offsetDiag];
                            } else if (ptr == 1) {
                                CS[(k%3)*fLen+offset] = CI[(k%2)*fLen+offset];
                            } else {
                                CS[(k%3)*fLen+offset] = CD[(k%2)*fLen+offset];
                            } 
                        }
                        // if (tile == 0 && tidx == 0) printf("k: %d, CS: %d, CI: %d, CD: %d, Iptr: %d, Dptr: %d\n",k, CS[k%3*fLen+offset], CI[k%2*fLen+offset], CD[k%2*fLen+offset], Iptr, Dptr);
                        if (Iptr) {
                            // std::cout << (ptr & 0xFF) << " ";
                            ptr |= 0x04; 
                            // std::cout << (ptr & 0xFF) << "\n";
                        }
                        if (Dptr) {
                            // std::cout << (ptr & 0xFF) << " ";
                            ptr |= 0x08;
                            // std::cout << (ptr & 0xFF) << "\n";
                        }
                        // if (tidx == 0) printf("Before KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);                      
                        if (k <= p_marker){
                            tb[tb_idx+tidx] = ptr;
                        }
                        // if (tidx == 0) printf("After KK idx0: %d, idx1: %d\n", ftr_length_idx, ftr_lower_limit_idx);
                    }
                    __syncthreads();
                    // Calculate Max
                    for (uint32_t r = min(leftGrid, (int32_t)threadNum)/2; r > 0; r >>= 1) {
                        if (tidx < r) {
                            // printf("Reduction max No.%d & %d, max_score: %d & %d, pre_max: %d\n", tidx, tidx+r, max_score_list[tidx], max_score_list[tidx+r], max_score_prime);
                            max_score_list[tidx]  = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_list[tidx+r] : max_score_list[tidx];
                            max_score_ref[tidx]   = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_ref[tidx+r] : max_score_ref[tidx];
                            max_score_query[tidx] = (max_score_list[tidx+r] > max_score_list[tidx]) ? max_score_query[tidx+r] : max_score_query[tidx];

                        }
                        __syncthreads();
                    }
                    // Update Max
                    if (max_score_prime < max_score_list[0]) {
                        max_score_prime = max_score_list[0];
                        if (tidx == 0) {
                            max_score_prime = max_score_list[0];
                            if (k <= p_marker) {
                                max_score_ref_idx = max_score_ref[0];
                                max_score_query_idx = max_score_query[0];
                                max_score_start_addr = ftr_addr - (U[k%3] - L[k%3] + 1)  + (max_score_query[0] - L[k%3]);
                                max_score_start_ftr = k;
                            }
                            // printf("max (%d, %d) = %d, %d, %d\n",max_score_ref_idx, max_score_query_idx, max_score_prime, max_score_start_addr, max_score_list[0]);
                        }
                    }
                    __syncthreads();
                    
                    __syncthreads();
                }
                __syncthreads();
                if (k <= p_marker){
                    tb_idx += min((int)leftGrid, (int)threadNum);
                }
                __syncthreads();
                // if (tidx == 0) printf("DEBUG2: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                int32_t newL = L[k%3];
                int32_t newU = U[k%3];
                // if (tidx == 0) printf("k: %d,newL: %d, newU: %d\n",k, L[k%3], U[k%3]);
                while (newL <= U[k%3]) {         
                    int32_t offset = newL - L[k%3];
                    // if (tidx == 0) printf("DDDDD: %d, %d, %d, %d\n", newU, L[k%3], offset, S[(k%3)*fLen+offset]);
                    // printf("k: %d,newL: %d, offset: %d\n",S[(k%3)*fLen+offset], L[k%3], offset);
                    if (S[(k%3)*fLen+offset] <= -inf) {
                        newL++;
                    }
                    else {
                        break;
                    }
                }
                while (newU >= L[k%3]) {
                    int32_t offset = newU - L[k%3];
                    // if (tidx == 0) printf("DDDDD: %d, %d, %d, %d\n", newU, L[k%3], offset, S[(k%3)*fLen+offset]);
                    if (S[(k%3)*fLen+offset] <= -inf) {
                        newU--;
                    }
                    else {
                        break;
                    }
                }
                int32_t v1 = static_cast<int32_t>(query_length)-1;
                int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                int32_t v3 = newU+1;
                int32_t Lprime = max(static_cast<int32_t>(0), v2);
                __syncthreads();
                // if (tidx == 0) printf("newU: %d, + %d = %d ", newU, 1, v3);
                // if (tidx == 0) printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                // if (tidx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
        
                L[(k+1)%3] = max(newL, Lprime);
                U[(k+1)%3] = min(v1, v3); 
        
                // if (tidx == 0) printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
        
                
                // sL[(k+1)%3] = L[(k+1)%3];
                // sU[(k+1)%3] = U[(k+1)%3];

                if (tidx == 0) {
                    // printf("k: %d,newL: %d, newU: %d\n",k, newL, newU);
                    if ((!converged) && (k < reference_length + query_length - 2)) {
                        int32_t start = newL - L[k%3];
                        int32_t length = newU - newL;
                        int32_t conv_I = CI[(k%2)*fLen+start];
                        int32_t conv_D = CD[(k%2)*fLen+start];
                        int32_t conv_S = CS[(k%3)*fLen+start];
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_I != CI[(k%2)*fLen+i]){
                                conv_I = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_D != CD[(k%2)*fLen+i]){
                                conv_D = -1;
                                break;
                            } 
                        }
                        for (int32_t i = start + 1; i <= start + length; i++ ){
                            if (conv_S != CS[(k%3)*fLen+i]){
                                conv_S = -1;
                                break;
                            } 
                        }
                        // printf("k: %d, %d, %d, %d\n", k, conv_I, conv_D, conv_S);
                        if ((conv_I == conv_D) && (conv_I == conv_S) && (prev_conv_s == conv_S) && (conv_I != -1)){
                            converged = true; 
                            conv_value = prev_conv_s;
                            conv_score = max_score_prime;
                            // if (DEBUG)  std::cout << "Converged at: " << conv_value << "\n";
                            // printf("Converged at: %d\n", conv_value);
                        }
                        prev_conv_s = conv_S;
                    }
                    // int32_t v1 = static_cast<int32_t>(query_length)-1;
                    // int32_t v2 = static_cast<int32_t>(k)+2-static_cast<int32_t>(reference_length);
                    // int32_t v3 = newU+1;

                    // int32_t Lprime = max(static_cast<int32_t>(0), v2);

                    // printf("v1: %d, v2: %d, v3: %d\n", v1, v2, v3);
                    // // if (tidx == 0) printf("DEBUG3: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
            
                    // L[(k+1)%3] = max(newL, Lprime);
                    // U[(k+1)%3] = min(v1, v3); 
            
                    // printf("DEBUG4: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                    // printf("DEBUG5 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
            
                    
                    // sL[(k+1)%3] = L[(k+1)%3];
                    // sU[(k+1)%3] = U[(k+1)%3];
    

                    // printf("DEBUG44: L: %d, %d, %d, U:%d, %d, %d\n", L[k%3], L[(k+1)%3], L[(k+2)%3], U[k%3], U[(k+1)%3], U[(k+2)%3]);
                    // printf("DEBUG6 L: %d, %d, %d, U:%d, %d, %d\n", sL[k%3], sL[(k+1)%3], sL[(k+2)%3], sU[k%3], sU[(k+1)%3], sU[(k+2)%3]);
                    // Update max_score
                    max_score = max_score_prime;
                    // printf("Max score: %d\n", max_score);

                    if ((converged) && (max_score > conv_score)){
                        conv_logic = true;
                        // printf("Tile %d Convergence logic found: %d\n", tile, max_score);
                        // if (DEBUG) std::cout << "Convergence logic found: ";
                        // break;
                    }
                    
                }
                __syncthreads();
                if (conv_logic) break;
            }
            __syncthreads();
            // if (tidx == 0) printf("DEBUG: break 2\n");
            if (tidx == 0) {
                // printf("Frontier addr: %d \ntb_start_ftr: %d\nmarker: %d\n", ftr_addr, ftr_length_idx, p_marker);
                int32_t conv_ref_idx = 0; 
                int32_t conv_query_idx = 0;
                int32_t tb_start_addr = 0; 
                int32_t tb_start_ftr = 0;  
                // printf("ftr_addr: %d\n", ftr_addr);
                // printf("ftr_l_idx: %d\n", ftr_length_idx);
                // printf("ftr[]: %d\n", ftr_length[ftr_length_idx - 1]);
                // printf("limit_idx: %d\n", ftr_lower_limit_idx);
                // printf("limit[]: %d\n", ftr_lower_limit[ftr_lower_limit_idx - 1]);
                if (conv_logic) {       
                    conv_query_idx = conv_value & 0xFFFF;
                    tb_state = (conv_value >> 16) & 0xFFFF;
                    conv_ref_idx = p_marker - conv_query_idx; 
                    conv_ref_idx -= (tb_state == 3) ? 1: 0;
                    tb_start_addr = ftr_addr - ftr_length[ftr_length_idx - 1];
                    tb_start_addr = (tb_state == 3) ? tb_start_addr - ftr_length[ftr_length_idx - 2] + (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 2]) : 
                                                      tb_start_addr +  (conv_query_idx - ftr_lower_limit[ftr_lower_limit_idx - 1]);
                    tb_start_ftr = (tb_state == 3) ? ftr_length_idx - 2: ftr_length_idx - 1;
                    // if (DEBUG) std::cout <<  " conv query idx: " << conv_query_idx << " " << (tb_state&0xFFFF) << " " << conv_ref_idx << " " << conv_value << std::endl;
                    // printf(" conv query idx: %d, %d, %d, %d\n", conv_query_idx ,(tb_state&0xFFFF) ,conv_ref_idx , conv_value );
                
                } else {
                    conv_query_idx = max_score_query_idx;
                    conv_ref_idx = max_score_ref_idx;
                    tb_start_addr = max_score_start_addr;
                    tb_start_ftr = max_score_start_ftr;
                    tb_state = 0;
                    last_tile = true;
                }
                // printf("Before: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                reference_idx += conv_ref_idx;
                query_idx += conv_query_idx;
                // printf("After: convref: %d, convqry:%d\n", conv_ref_idx, conv_query_idx);
                // printf("After: refidx: %d, qryidx:%d\n", reference_idx, query_idx);
                // Start Traceback
                // printf("DEBUG: break 3\n");
                // Traceback(ftr_length, ftr_lower_limit, tb_start_addr, tb_start_ftr, (tb_state%3), conv_query_idx, conv_ref_idx, tb, aln);
                int32_t addr = tb_start_addr; 
                int16_t ftr = tb_start_ftr;
                int16_t traceback_idx = conv_query_idx;
                int16_t qry_idx = conv_query_idx;
                int16_t ref_idx = conv_ref_idx;
                int8_t  tb_value = 0;
                state = tb_state%3;

                int8_t  dir = 0;
                // bool checkpoint = false;
                // int count = 0;
                // if(tile == 155) printf("ftr: %d\n", ftr);
                // if(tile == 155) printf("state: %d\n", tb_state);
                // if(tile == 155) printf("addr: %d\n", addr);
                // if(tile == 155) printf("aln_idx: %d\n", aln_idx);
                // if(tile == 155) printf("dir: %d\n", dir);
                // if(tile == 155) printf("traceback_idx: %d\n", traceback_idx);
                while (ftr >= 0) {
                    // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
                    if (addr < 0) {
                        printf("ERROR: tb addr < 0 (%d)!\n", addr);
                        // exit(1);
                        break;
                    }
                    // count++;
                    tb_value = tb[addr];
                    // if (DEBUG)
                    // {
                    // printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, idx , ftr_lower_limit[ftr]);
                    // printf(" fL[ftr - 1]: %d, ll[ftr-2]: %d\n", ftr_length[ftr - 1], ftr_lower_limit[ftr-2]);
                    // printf(" fL[ftr - 2]: %d", ftr_length[ftr - 2]);
                    // printf(" Tb: %d", ( tb_value&0xFFFF));
                    // }
                    // if (tile < 10) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d, tb_value: %d\n", addr, (state & 0xFFFF), ftr, traceback_idx , ftr_lower_limit[ftr], tb_value);
                    if (state == 0) { // Current State M
                        state = tb_value & 0x03;
                        if (state == 0) {
                            dir = 0;
                        }
                        else if (state == 1) {
                            dir = 1;
                            if (tb_value & 0x04) {
                                state = 1;
                            } else {
                                state = 0;
                            }   
                        }
                        else {
                            dir = 2;
                            if (tb_value & 0x08) {
                                state = 2;
                            } else {
                                state = 0;
                            }
                        }
                    } 
                    else if (state == 1) { // Current State I
                        dir = 1;
                        if (tb_value & 0x04) {
                            state = 1;
                        } else {
                            state = 0;
                        }
                    } 
                    else { // Current State D
                        dir = 2;
                        if (tb_value & 0x08) {
                            state = 2;
                        } else {
                            state = 0;
                        }
                    }
                    // if (tile == 155 || tile == 140) printf("Start Addr:%d state:%d, ftr:%d, idx:%d, ll[ftr-1]:%d\n", addr, (state & 0xFFFF), ftr, traceback_idx , ftr_lower_limit[ftr]);
                    addr = addr - (traceback_idx  - ftr_lower_limit[ftr] + 1) - (ftr_length[ftr - 1]);

                    if (dir == 0) {
                        addr = addr - (ftr_length[ftr - 2]) + (traceback_idx - ftr_lower_limit[ftr  - 2]);
                        ftr -= 2;
                        traceback_idx -= 1;
                        qry_idx--;
                        ref_idx--;
                    }
                    else if (dir == 1) {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1]);
                        ftr -= 1;
                        traceback_idx -=1;
                        qry_idx--;
                    }
                    else {
                        addr = addr + (traceback_idx - ftr_lower_limit[ftr  - 1] + 1);
                        ftr -= 1;
                        ref_idx--;
                    }

                    
                    tile_aln[aln_idx] = dir;
                    ++aln_idx;    
                    // state = next_state;
                    // if (DEBUG)  std::cout << " Final State: " << (state&0xFF) << " End Addr: " << addr << " index:" << ref_idx << ", " << query_idx << std::endl;
                }
                // std::cout << count << std::endl;   

                // if (DEBUG) std::cout << ref_idx << " " << qry_idx << std::endl; 
                
                state = tb_state % 3;
                // if (DEBUG) {
                //     std::cout << "tb_state: " <<  (tb_state & 0xFFFF) << std::endl;
                //     int count = 0;
                //     for (auto &a: aln){
                //         std::cout << count << ": " << (a & 0xFFFF) << "\n";
                //         count += 1;
                //     }
                // }

                // Write global memory
                int32_t refIndex = reference_idx;
                int32_t qryIndex = query_idx;
                int16_t globalAlnIdx = idx[5];
                // printf("TB: refidx: %d, qryidx:%d, globalAln: %d\n", refIndex, qryIndex, globalAlnIdx);
                for (int j = aln_idx -1; j >= 0; --j) {
                    // printf("i: %d, dir:%d, globalAln: %d\n", j, tile_aln[j], globalAlnIdx);
                    if (tile_aln[j] == 0) {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = ref[i*seqLen+refIndex]; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = qry[i*seqLen+qryIndex]; 
                        qryIndex++;refIndex++;
                    }
                    else if (tile_aln[j] == 1) {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = '-'; 
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = qry[i*seqLen+qryIndex]; 
                        qryIndex++;
                    }
                    else {
                        for (size_t i=0; i<refNum; i++) alignment[i*seqLen+globalAlnIdx] = ref[i*seqLen+refIndex];  
                        for (size_t i=0; i<qryNum; i++) alignment[(i+refNum)*seqLen+globalAlnIdx] = '-'; 
                        refIndex++;
                    }
                    globalAlnIdx++;
                }
                // printf("Finish TB of tile %d!!!!\n", tile);
                idx[3] = reference_idx;
                idx[4] = query_idx;
                idx[5] = globalAlnIdx;
                // tile++;
            }
           
            __syncthreads();
            if (tidx < threadNum) tile++;
        }
        __syncthreads();
        
        // __syncthreads();
        // TODO: Add score to debug
        // score = Score(params, aln, reference[n], query[n], reference_idx, query_idx);
    }


}

#ifndef ALIGN_HPP
#include "align.hpp"
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



void alignGrpToGrp
(
    std::vector<std::string>& ref,
    std::vector<std::string>& query,
    Params& param,
    std::pair<std::vector<std::string>, std::vector<std::string>>& alignment
){
    if (ref.size() < 1 || query.size() < 1) {fprintf(stderr, "Error: Number of Ref/Query <= 0"); exit(1);}
    int32_t refLen = ref[0].size();
    int32_t queryLen = query[0].size();
    if (refLen<=0 || queryLen<=0) {fprintf(stderr, "Error: Ref/Query length <= 0"); exit(1);}
    
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


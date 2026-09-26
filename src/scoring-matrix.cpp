#ifndef MSA_HPP
#include "msa.hpp"
#endif

#ifndef BLOSUM_HPP
#include "blosum.hpp"
#endif

#include <iostream>

char checkOnly(char inChar) {
    switch (inChar) {
        case 'E': return 'p';
        case 'F': return 'p';
        case 'I': return 'p';
        case 'J': return 'p';
        case 'L': return 'p';
        case 'P': return 'p';
        case 'Q': return 'p';
        case 'Z': return 'p';
        case 'U': return 'n';
        default:  return 'x';
    }
}

constexpr std::array<int, 256> make_aa_lookup()
{
    std::array<int, 256> table{};
    table['A'] = 1;
    table['C'] = 2;
    table['D'] = 3;
    table['E'] = 4;
    table['F'] = 5;
    table['G'] = 6;
    table['H'] = 7;
    table['I'] = 8;
    table['K'] = 9;
    table['L'] = 10;
    table['M'] = 11;
    table['N'] = 12;
    table['P'] = 13;
    table['Q'] = 14;
    table['R'] = 15;
    table['S'] = 16;
    table['T'] = 17;
    table['V'] = 18;
    table['W'] = 19;
    table['Y'] = 20;
    table['-'] = 22;
    table['.'] = 22;
    return table;
}

constexpr std::array<int, 256> make_nb_lookup()
{
    std::array<int, 256> table{};
    table['A'] = 1;
    table['C'] = 2;
    table['G'] = 3;
    table['T'] = 4;
    table['U'] = 4;
    table['-'] = 6;
    table['.'] = 6;
    return table;
}

constexpr auto aa_lookup = make_aa_lookup();
constexpr auto nb_lookup = make_nb_lookup();

int letterIdx(char type, char c) {
    if (type == 'p') {
        int val = aa_lookup[static_cast<unsigned char>(c)];
        return (val == 0) ? 20 : val - 1;
    }
    else {
        int val = nb_lookup[static_cast<unsigned char>(c)];
        return (val == 0) ? 4 : val - 1;
    }   
}

msa::Params::Params(po::variables_map& vm, char type) {
    bool userDefine = vm.count("matrix");
    float userGOP = vm["gap-open"].as<float>();
    float userGEP = vm["gap-extend"].as<float>();
    float userGBoundary = vm.count("gap-ends") ? vm["gap-ends"].as<float>() : userGEP;
    float userOffset = vm.count("offset") ? vm["offset"].as<float>() : 0.0f;
    float xdrop = round(vm["xdrop"].as<float>());
    if (userGOP > 0 || userGEP > 0 || userGBoundary > 0)  {
        std::cerr << "ERROR: Gap penalties must be less than or equal to 0.\n";
        exit(1);
    }
    if (xdrop <= 0)  {
        std::cerr << "ERROR: XDrop value should be larger than 0.\n";
        exit(1);
    }
    float targetScore = vm.count("scale") ? vm["scale"].as<float>() : 60.0f;
    float userLocalGOP = vm.count("local-gap-open") ? vm["local-gap-open"].as<float>() : -2.00f;
    float userLocalGEP = vm.count("local-gap-extend") ? vm["local-gap-extend"].as<float>() : -0.10f;
    if (userLocalGOP > 0 || userLocalGEP > 0) {
        std::cerr << "ERROR: Local gap penalties must be less than or equal to 0.\n";
        exit(1);
    }
    this->offset = userOffset * targetScore;
    this->gapOpen = userGOP * targetScore;
    this->gapExtend = userGEP * targetScore;
    this->gapBoundary = userGBoundary * targetScore;
    this->localGapOpen = userLocalGOP * targetScore;
    this->localGapExtend = userLocalGEP * targetScore;
    this->gapTerminal = 0.0f;
    this->xdrop =  (this->gapExtend == 0) ? xdrop : -1*xdrop*this->gapExtend;

    this->matrixSize = (type == 'n') ? 5 : 21;
    this->scoringMatrix = new float* [this->matrixSize];
    for (int i = 0; i < this->matrixSize; ++i) this->scoringMatrix[i] = new float[this->matrixSize];
        
    int blosumType = vm["blosum"].as<int>();
    if (!userDefine) {
        if (type == 'n') {
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    if (i == 4 || j == 4)        this->scoringMatrix[i][j] = vm.count("wildcard") ? vm["match"].as<float>() : 0.0;
                    else if (i == j)             this->scoringMatrix[i][j] = vm["match"].as<float>();
                    else if (std::abs(i-j) == 2) this->scoringMatrix[i][j] = vm["transition"].as<float>();
                    else                         this->scoringMatrix[i][j] = vm["mismatch"].as<float>();
                }
            }
        }
        else if (type == 'p') {
            if (blosumType != 43 && blosumType != 45 && blosumType != 62 && blosumType != 80 && blosumType != 99) {
                std::cerr << "WARNING: Invalid substitution matrix \"" << blosumType
                          << "\". Please choose from 43 (PFASUM43), 45, 62, 80, or 99 (MIQS).\n";
                std::cerr << "Using default: BLOSUM62.\n";
                blosumType = 62;
            }
            float Nscore = 0.0;
            if (blosumType == 43) {
                for (int i = 0; i < 20; ++i) Nscore += PFASUM43[i][i];
            } else if (blosumType == 99) {
                for (int i = 0; i < 20; ++i) Nscore += MIQS[i][i];
            } else if (blosumType == 45) {
                for (int i = 0; i < 20; ++i) Nscore += BLOSUM45[i][i];
            } else if (blosumType == 80) {
                for (int i = 0; i < 20; ++i) Nscore += BLOSUM80[i][i];
            } else {
                for (int i = 0; i < 20; ++i) Nscore += BLOSUM62[i][i];
            }
            Nscore /= 20;
            for (int i = 0; i < 21; ++i) {
                this->scoringMatrix[i][20] = vm.count("wildcard") ? 5 * Nscore : 0.0;
                this->scoringMatrix[20][i] = vm.count("wildcard") ? 5 * Nscore : 0.0;
            }
            if (blosumType == 43) {
                for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) this->scoringMatrix[i][j] = 5 * PFASUM43[i][j];
            }
            else if (blosumType == 62) {
                for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) this->scoringMatrix[i][j] = 5 * BLOSUM62[i][j];
            }
            else if (blosumType == 45) {
                for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) this->scoringMatrix[i][j] = 5 * BLOSUM45[i][j];
            }
            else if (blosumType == 80) {
                for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) this->scoringMatrix[i][j] = 5 * BLOSUM80[i][j];
            }
            else if (blosumType == 99) {
                for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) this->scoringMatrix[i][j] = 5 * MIQS[i][j];
            }
            
        }
    }
    else {
        std::string matrixFileName = vm["matrix"].as<std::string>();
        std::ifstream matrixFile(matrixFileName);
        if (!matrixFile) {
            fprintf(stderr, "ERROR: can't open %s\n", matrixFileName.c_str());
            exit(1);
        }
        std::string word;
        std::vector<int> charVec;
        int readCount = 0, charNum = this->matrixSize-1;
        while (matrixFile >> word) {
            if (readCount == charNum) {
                bool isNumber = true;
                try {
                    size_t pos;
                    std::stod(word, &pos);
                    if (pos != word.size()) isNumber = false; // extra non-numeric chars
                } catch (...) {
                    isNumber = false;
                }
                if (!isNumber) {
                    charNum = this->matrixSize;
                }
            }
            if (readCount < charNum) {
                char letter = toupper(word[0]);
                int ambig = (type == 'n') ? 4 : 20;
                if (letterIdx(type, letter) == ambig && charNum == this->matrixSize-1) {
                    std::string seqType = (type == 'n') ? " for nucleotide sequences.\n" : " for protein sequences.\n";
                    std::cerr << "Unrecognized letter \"" << letter << "\"" << seqType;
                    exit(1);
                }
                charVec.push_back(letterIdx(type, letter));
                readCount++;
            }
            else {
                int x = (readCount-charNum) / charNum;
                int y = (readCount-charNum) % charNum;
                int i = charVec[x];
                int j = charVec[y];
                this->scoringMatrix[i][j] = std::stof(word.c_str());
                readCount++;
            }
        }
        matrixFile.close();
        if (charNum == this->matrixSize-1) {
            float Nscore = 0;
            for (int i = 0; i < charNum; ++i) Nscore += this->scoringMatrix[i][i];
            Nscore = vm.count("wildcard") ? (Nscore / charNum) : 0.0;
            for (int i = 0; i < this->matrixSize; ++i) {
                this->scoringMatrix[i][this->matrixSize-1] = Nscore;
                this->scoringMatrix[this->matrixSize-1][i] = Nscore;
            }
        }
        
    }
    if (!userDefine) this->normalizeScoringMatrix(targetScore, this->offset);

    std::map<char, int> letterMap;
    if (type == 'n') {
        letterMap = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'U', 3}, {'-', 5}};
    }
    else {
        letterMap = {{'A',0}, {'C',1}, {'D',2}, {'E',3}, {'F',4}, {'G',5}, {'H',6}, {'I',7}, {'K',8}, {'L',9}, {'M',10}, {'N',11}, {'P',12}, {'Q',13}, {'R',14}, {'S',15}, {'T',16}, {'V',17}, {'W',18}, {'Y',19}, {'-',21}}; // 20 for all other characters (ambiguous)
    }

    this->gapTerminal = 0.0f;
    
    if (vm.count("verbose")) {
        std::cerr << "======== Parameters ========\n";
        if (type == 'p') {
            if      (blosumType == 43) std::cerr << "PFASUM43\n";
            else if (blosumType == 45) std::cerr << "BLOSUM45\n";
            else if (blosumType == 62) std::cerr << "BLOSUM62\n";
            else if (blosumType == 80) std::cerr << "BLOSUM80\n";
            else if (blosumType == 99) std::cerr << "MIQS\n";
        }
        std::cerr << std::setw(5) << " ";
        for (size_t i = 0; i < this->matrixSize-1; ++i) {
            auto letter = letterMap.begin();
            std::advance(letter, i);
            letter++;
            std::cerr << std::setw(5) << letter->first;
        }
        std::cerr << std::setw(5) << ((type == 'n') ? 'N' : 'X');
        std::cerr << "\n";
        for (size_t i = 0; i < this->matrixSize; ++i) {
            auto letter = letterMap.begin();
            std::advance(letter, i);
            letter++;
            if (i < this->matrixSize-1) std::cerr << std::setw(5) << letter->first;
            else std::cerr << std::setw(5) << ((type == 'n') ? 'N' : 'X');
            for (size_t j = 0; j < this->matrixSize; ++j) {
                std::cerr << std::setw(5) << static_cast<int>(this->scoringMatrix[i][j]);
            }
            std::cerr << "\n";
        }
        std::cerr << "Gap-Open:       " << this->gapOpen << " (normalized: " << userGOP << ")\n"
                  << "Gap-Extend:     " << this->gapExtend << " (normalized: " << userGEP << ")\n"
                  << "Gap-Terminal:   " << this->gapTerminal << "\n"
                  << "Offset:         " << this->offset << " (normalized: " << userOffset << ")\n"
                  << "Xdrop:          " << this->xdrop << '\n';
    }
}

void msa::Params::normalizeScoringMatrix(float targetScore, float offset) {
    
    int standardSize = matrixSize - 1; 
    std::vector<float> freq(standardSize, 0.0f);

    if (standardSize == 4) {
        for (int i = 0; i < standardSize; ++i) {
            freq[i] = 0.25f;
        }
    } else if (standardSize == 20) {
        freq = {
            0.077f, 0.020f, 0.052f, 0.062f, 0.040f, 0.074f, 0.023f, 0.052f, 0.059f, 0.091f, 
            0.024f, 0.043f, 0.051f, 0.041f, 0.051f, 0.069f, 0.059f, 0.066f, 0.014f, 0.032f
        };
    } else {
        for (int i = 0; i < standardSize; ++i) freq[i] = 1.0f / standardSize;
    }

    // Step 1: Zero-Expectation Centering on standard alphabet
    float bg_average = 0.0f;
    for (int i = 0; i < standardSize; ++i) {
        for (int j = 0; j < standardSize; ++j) {
            bg_average += scoringMatrix[i][j] * freq[i] * freq[j];
        }
    }

    for (int i = 0; i < standardSize; ++i) {
        for (int j = 0; j < standardSize; ++j) {
            scoringMatrix[i][j] -= bg_average;
        }
    }

    // Step 2: Rescaling to targetScore (60-scale) on diagonal average
    float id_average = 0.0f;
    for (int i = 0; i < standardSize; ++i) {
        id_average += scoringMatrix[i][i] * freq[i];
    }

    float scale = (id_average != 0.0f) ? (targetScore / id_average) : 1.0f;
    for (int i = 0; i < standardSize; ++i) {
        for (int j = 0; j < standardSize; ++j) {
            scoringMatrix[i][j] *= scale;
        }
    }

    // Step 3: Offset subtraction on standard alphabet
    if (offset != 0.0f) {
        for (int i = 0; i < standardSize; ++i) {
            for (int j = 0; j < standardSize; ++j) {
                scoringMatrix[i][j] -= offset;
            }
        }
    }

    // Step 4: Wildcard/ambiguous character reset to strictly 0.0
    int wildcardIdx = matrixSize - 1;
    for (int i = 0; i < matrixSize; ++i) {
        scoringMatrix[i][wildcardIdx] = 0.0f;
        scoringMatrix[wildcardIdx][i] = 0.0f;
    }
}


msa::Params::~Params() {
    for (int i = 0; i < this->matrixSize; ++i) delete[] this->scoringMatrix[i];
    delete [] scoringMatrix;
}
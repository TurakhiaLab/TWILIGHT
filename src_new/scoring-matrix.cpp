#ifndef MSA_HPP
#include "msa.hpp"
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
    this->gapOpen = vm["gap-open"].as<float>();
    this->gapExtend = vm["gap-extend"].as<float>();
    this->gapBoundary = vm.count("gap-ends") ? vm["gap-ends"].as<float>() : this->gapExtend;
    float xdrop = round(vm["xdrop"].as<float>());
    if (this->gapOpen > 0 || this->gapExtend > 0 || this->gapBoundary > 0)  {
        std::cerr << "ERROR: Gap penalties must be less than or equal to 0.\n";
        exit(1);
    }
    if (xdrop <= 0)  {
        std::cerr << "ERROR: XDrop value should be larger than 0.\n";
        exit(1);
    }
    this->xdrop =  (this->gapExtend == 0) ? xdrop : -1*xdrop*this->gapExtend;

    this->matrixSize = (type == 'n') ? 5 : 21;
    this->scoringMatrix = new float* [this->matrixSize];
    for (int i = 0; i < this->matrixSize; ++i) this->scoringMatrix[i] = new float[this->matrixSize];
        
            
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
            float Nscore = 0.0;
            for (int i = 0; i < 20; ++i) Nscore += BLOSUM62[i][i];
            Nscore /= 20;
            for (int i = 0; i < 21; ++i) {
                for (int j = 0; j < 21; ++j) {
                    if (i == 20 || j == 20) this->scoringMatrix[i][j] = vm.count("wildcard") ? 5 * Nscore : 0.0;
                    else                    this->scoringMatrix[i][j] = 5 * BLOSUM62[i][j];
                }
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

    std::map<char, int> letterMap;
    if (type == 'n') {
        letterMap = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'U', 3}, {'-', 5}};
    }
    else {
        letterMap = {{'A',0}, {'C',1}, {'D',2}, {'E',3}, {'F',4}, {'G',5}, {'H',6}, {'I',7}, {'K',8}, {'L',9}, {'M',10}, {'N',11}, {'P',12}, {'Q',13}, {'R',14}, {'S',15}, {'T',16}, {'V',17}, {'W',18}, {'Y',19}, {'-',21}}; // 20 for all other characters (ambiguous)
    }
    
    if (vm.count("verbose")) {
        std::cout << "======== Parameters ========\n";
        std::cout << std::setw(5) << " ";
        for (size_t i = 0; i < this->matrixSize-1; ++i) {
            auto letter = letterMap.begin();
            std::advance(letter, i);
            letter++;
            std::cout << std::setw(5) << letter->first;
        }
        std::cout << std::setw(5) << ((type == 'n') ? 'N' : 'X');
        std::cout << "\n";
        for (size_t i = 0; i < this->matrixSize; ++i) {
            auto letter = letterMap.begin();
            std::advance(letter, i);
            letter++;
            if (i < this->matrixSize-1) std::cout << std::setw(5) << letter->first;
            else std::cout << std::setw(5) << ((type == 'n') ? 'N' : 'X');
            for (size_t j = 0; j < this->matrixSize; ++j) {
                std::cout << std::setw(5) << this->scoringMatrix[i][j];
            }
            std::cout << "\n";
        }
        std::cout << "Gap-Open:   " << this->gapOpen << "\n"
                  << "Gap-Extend: " << this->gapExtend << "\n"
                  << "Gap-Ends:   " << this->gapBoundary << "\n"
                  << "Xdrop:      " << this->xdrop << '\n';
        std::cout << "============================\n";
    }
}

msa::Params::~Params() {
    for (int i = 0; i < this->matrixSize; ++i) delete[] this->scoringMatrix[i];
    delete [] scoringMatrix;
}
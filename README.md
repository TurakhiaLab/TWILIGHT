# TWILIGHT

## Introduction
TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It utilizes the TALCO method to efficiently perform global alignments with high accuracy. TWILIGHT is compatible with Linux systems and is able to take advantage of CUDA-capable GPUs for further acceleration.
## Build Instructions
#### Please make sure if you have the below libraries installed
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
#### Install and compile
```
git clone https://github.com/TurakhiaLab/TWILIGHT.git
mkdir TWILIGHT/build && cd TWILIGHT/build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make twilight
```

## Run Instructions
#### See Help for more details
```
./twilight -h
```
#### Build MSA from raw sequences
```
./twilight -t <tree file> -i <sequence file> -o <output file>
```
#### Merge multiple MSA files
```
./twilight -f <directory containing all MSA files> -o <output file>
```
#### For large dataset, divide into multiple subalignments and align sequentially to reduce memory usage
```
./twilight -t <tree file> -i <sequence file> -o <output file> -d <temporary directory> -a <max subalignment size> --merge-subtree <merger method>
```


## Sample commands for the provided test data
#### Build MSA on difficult alignment (more gappy) 
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln --psgop y
```
#### Build MSA on short-branched sequences
```
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln --gappy 1 --psgop n
```
#### Divide into multiple subalignments and align sequentially
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -d RNASim_temp -a 200 --psgop y
```
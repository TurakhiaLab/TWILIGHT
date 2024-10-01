# MSA-Accel

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make msa-gpu
```

## Run Instructions
### Help
```
./msa-gpu -h
```
### You can create an output directory for the following example commands
```
mkdir ../output
```
### Run with default settings
```
./msa-gpu -t ../dataset/sars_2000.nwk -i ../dataset/sars_2000.fa -o ../output/sars_2000.aln
```
### Run with transitivity merger
```
./msa-gpu -t ../dataset/RNASim_10000.tre -i ../dataset/RNASim_10000.fa -o ../output/RNASim_10000_merge.aln -l 500
```
### Run with the removing-gappy-colum feature
```
./msa-gpu -t ../dataset/RNASim_10000.tre -i ../dataset/RNASim_10000.fa -o ../output/RNASim_10000.aln --gappy-vertical 0.9 --gappy-horizon 1
```

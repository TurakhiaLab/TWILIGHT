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
### Run
```
mkdir output
./msa-gpu -t ../dataset/sars_2000.nwk -i ../dataset/sars_2000.fa -o ../output/sars_2000.aln
```
### Run with transitivity merger
```
mkdir output
./msa-gpu -t ../dataset/16384-taxon-trial-1.nwk -i ../dataset/16384-taxon-trial-1.fa -o ../output/16384-taxon-trial-1_merge_mod.aln  -s t -l 50
```




# MSA-Accel

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make test-msa
make msa-gpu
```

## Run Instructions
### Help
'''
./msa-gpu -h
'''
### Run baseline
* To enable output, uncomment test_msa.cpp line 275-282
'''
./test-msa -t ../dataset/gene1.best.treefile -s ../dataset/gene1.aln
'''
### Test traditional approach
* To enable output, uncomment test_msa.cu line 557-564
* test_msa.cu line 291 (blockSize) has to match align.cu line 633 (threadNum)
```
./msa-gpu -t ../dataset/gene1.best.treefile -s ../dataset/gene1.aln
```
### Test TALCO approach
* modify test_msa.cu line 291 (blockSize) to 128
* comment test_msa.cu line 303-309
* uncomment test_msa.cu line 295-301
* There are some bugs
'''
./msa-gpu -t ../dataset/sars_20.nwk -s ../dataset/sars_20.fa 
'''



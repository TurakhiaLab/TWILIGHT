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
```
./msa-gpu -h
```
### Run
```
./msa-gpu -t ../dataset/sars_20.nwk -s ../dataset/sars_20.fa 
```



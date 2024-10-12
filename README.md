# MSA-Accel

## Build Instructions
```
mkdir build
cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make msa-accel
```

## Run Instructions
### Please see 'Help' for more details
```
./msa-accel -h
```
### Run with default settings
```
./msa-accel -t <tree file> -i <sequence file> -o <output file>
```
### Run with transitivity merger
```
./msa-accel -t <tree file> -i <sequence file> -o <output file> -l <maximum subtree size>
```
### Run with the removing-gappy-column feature (Recommanded)
```
./msa-accel -t <tree file> -i <sequence file> -o <output file> --gappy-vertical 0.9 --gappy-horizon 1
```
### For large dataset, align subtrees sequentially to reduce memory usage
```
./msa-accel -t <tree file> -i <sequence file> -o <output file> --gappy-vertical 0.9 --gappy-horizon 1 --temp-dir <directory for stoting temporary files> -m <maximum subtree size> --merge-subtree <merger method>
```

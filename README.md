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
### Run with the default setting
```
./msa-accel -t <tree file> -i <sequence file> -o <output file>
```
### Run with transitivity merger
```
./msa-accel -t <tree file> -i <sequence file> -o <output file> -l <maximum subtree size>
```

### For large dataset, align subtrees sequentially to reduce memory usage
```
./msa-accel -t <tree file> -i <sequence file> -o <output file> --temp-dir <directory for storing temporary files> -m <maximum subtree size> --merge-subtree <merger method>
```

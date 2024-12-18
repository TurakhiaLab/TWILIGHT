# TWILIGHT

## Build Instructions
#### Please make sure if you have the below libraries installed
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
#### Install and compile
```
git clone https://github.com/TurakhiaLab/msa-accel.git
mkdir msa-accel/build && cd msa-accel/build
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
./twilight -t ../dataset/RNASim_10000.nwk -i ../dataset/RNASim_10000.fa -o RNASim_10000.aln --psgop y
```
#### Build MSA on short-branched sequences
```
./twilight -t ../dataset/sars_2000.nwk -i ../dataset/sars_2000.fa -o sars_2000.aln --gappy 1 --psgop n
```
#### Divide into multiple subalignments and align sequentially
```
./twilight -t ../dataset/RNASim_10000.nwk -i ../dataset/RNASim_10000.fa -o RNASim_10000.aln -d RNASim_10000_temp -a 3000 --psgop y
```
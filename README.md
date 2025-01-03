# TWILIGHT

## Introduction
TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration.
## Install
### From source code (only tested on Linux)
#### Please make sure if you have the below libraries installed
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
#### Clone the repository and compile
```
git clone https://github.com/TurakhiaLab/TWILIGHT.git
mkdir TWILIGHT/build && cd TWILIGHT/build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```
### Using Docker locally (recommended for MAC users)
This docker image only supports the CPU version of TWILIGHT. If you would like to use the GPU-accelerated version, please build it from source code.
#### Clone the repository
```
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
```
#### Build the image and run the container
```
docker build -t twilight_image .
docker run -it twilight_image
```
## Run Instructions
#### See Help for more details
```
./twilight -h
```
#### Fast mode
```
./twilight -t <tree file> -i <sequence file> -o <output file>
```
#### Memory-friendly mode
```
./twilight -t <tree file> -i <sequence file> -o <output file> -d <temporary directory> -a <maximum subalignment size>
```
#### Merge mode
```
./twilight -f <directory containing all MSA files> -o <output file>
```

## Sample commands for the provided test data
#### Fast mode on gappy alignments
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y
```
#### Fast mode on short-branched sequences
```
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n
```
#### Memory-friendly mode
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -d RNASim_temp -a 200 -p y
```

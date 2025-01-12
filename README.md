# TWILIGHT

## Table of Contents
- [Overview](#overview)
- [Clone the repository](#clone)
- [Install dependencies](#install) 
- [Build TWILIGHT](#custom)
- [Docker](#docker)
- [Run TWILIGHT](#run)
- [Example Command](#example)

<br>

## <a name="overview"></a> Overview

TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration. 

TWILIGHT accepts unaligned sequences in FASTA format and an optional input guide tree in Newick format to generate output alignments in FASTA format. When a guide tree is unavailable, users can utilize the provided snakemake workflow to estimate trees using external tools.

## <a name="clone"></a> Clone the TWILIGHT repository.

```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
```

## <a name="install"></a> Install the required dependencies.

For TWILIGHT, please make sure the following libraries are installed in your system.
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
For iterative mode, please install snakemake `pip install snakemake` and tools you would like to use and configure the executable path in `workflow/config.yaml`,
- Initial guide tree: MashTree, MAFFT(PartTree)
- Subsequent iterations: FastTree, IQ-Tree, RAxML

## <a name="custom"></a> Build TWILIGHT.

```bash
cd TWILIGHT
mkdir build && cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```
## <a name="docker"></a> An alternative way: Use Docker.
Clone the repository and go into the TWILIGHT directory. Docker currently only supports vanilla mode.
```
docker build -t twilight_image .
docker run -it twilight_image
```

## <a name="run"></a> Run TWILIGHT
### TWILIGHT vanilla mode
See Help for more details
```
./twilight -h
```
Default mode
```
./twilight -t <tree file> -i <sequence file> -o <output file>
```
To reduce memory usage
```
./twilight -t <tree file> -i <sequence file> -o <output file> -d <temporary directory> -m <maximum subtree size>
```
Merge multiple MSAs
```
./twilight -f <directory containing all MSA files> -o <output file>
```
### TWILIGHT iterative mode
Go into the TWILIGHT directory and update the `workflow/config.yaml` file to include the path to your input unaligned sequence file, tree inference tools, the number of cpu threads, and other parameters.

Dry run to see the workflow
```
snakemake -n
```
Run iterative mode
```
snakemake
```
## <a name="example"></a> Example commands
Default mode
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y
```
Align short-branched sequences
```
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n
```
To reduce memory usage
```
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -d RNASim_temp -a 200 -p y
```

# TWILIGHT: Tall and WIde aLIGnments at High Throughput

## Table of Contents
- [Overview](#overview)
- [Quick start](#start)
  - [Installation](#install)
  - [Execution](#run)
  - [Examples](#example)
  - [Iterative](#iterative)
- [Citation](#cite)


<br>

## <a name="overview"></a> Overview

TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration. 

TWILIGHT default mode requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate output alignments in FASTA format. When a guide tree is unavailable, users can utilize the iterative mode which provides a snakemake workflow to estimate guide trees using external tools.

## <a name="start"></a> Quick start
### <a name="install"></a> Installation
#### Using Unix commands (requires sudo access to install dependant libraries)
Dependent libraries
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
Build TWILIGHT
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
mkdir build && cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```
#### Using Docker
Docker currently only supports default mode.
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
docker build -t twilight_image .
docker run -it twilight_image
```
### <a name="run"></a> Execution
For more information about TWILIGHT's options and instructions, see Help.
```bash
./twilight -h
```
To run TWILIGHT with defualt configuration:
```bash
./twilight -t <guide_tree> -i <input_fasta> -o <output_fasts>
```
Specify a maximum subtree size to reduce memory usage when memory usage is the major concern.
```bash
./twilight -t <guide_tree> -i <input_fasta> -o <output_fasts> -m <maximum_subtree_size>
```
To merge multiple MSAs, please move all MSA files into a folder.
```bash
./twilight -f <folder> -o <output_fasts>
```
### <a name="iterative"></a> Iterative mode
Iterative mode provides a snakemake workflow to estimate guide trees using external tools. Please first install snakemake `pip install snakemake` and tools you would like to use and configure the executable path in `workflow/config.yaml`.
- Initial guide tree: MashTree, MAFFT(PartTree)
- Subsequent iterations: FastTree, IQ-Tree, RAxML

Update the `workflow/config.yaml` file to include the path to your input sequence file, tree inference tools, the number of cpu threads, and other parameters.

Dry run to see the workflow
```bash
snakemake -n
```
Run iterative mode
```bash
snakemake
```
### <a name="Examples"></a> Example commands
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln # default configuration
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y # for gappy and divergent alignments
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y -d RNASim_temp -a 200 # reduce memory usage
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n # for short-branched sequences
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n -x ../dataset/subsitution.txt --gap-open -20 --gap-extend -4 # using user-defined scoring system
```
## <a name="cite"></a> Citation
TBA.
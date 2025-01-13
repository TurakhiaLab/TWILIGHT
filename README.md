# TWILIGHT: Tall and WIde aLIGnments at High Throughput

## Table of Contents
- [Overview](#overview)
- [Quick start](#start)
  - [Install TWILIGHT](#install)
    - [Using Unix commands](#unix)
    - [Using Docker](#docker)
    - [Example command](#example) 
  - [TWILIGHT modes](#mode)
    - [Default mode](#default)
    - [Iterative mode](#iterative)
- [Other example commands](#o_example)
- [Citation](#cite)


<br>

## <a name="overview"></a> Overview

TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration. 

By default, TWILIGHT requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate the output alignment in FASTA format. When a guide tree is unavailable, users can utilize the iterative mode, which provides a snakemake workflow to estimate guide trees using external tools.

## <a name="start"></a> Quick start
### <a name="install"></a> Install TWILIGHT
#### <a name="unix"></a> Using Unix commands (requires sudo access to install dependant libraries)
Dependent libraries.
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
Build TWILIGHT.
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
mkdir build && cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```
#### <a name="docker"></a> Using Docker
Docker currently only supports default mode.
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
docker build -t twilight_image .
docker run -it twilight_image
```
### <a name="example"></a> Example command
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln
```
## <a name="mode"></a> TWILIGHT modes
### <a name="default"></a> Default mode
For more information about TWILIGHT's options and instructions, see Help.
```bash
./twilight -h
```
Run TWILIGHT with default configuration.
```bash
./twilight -t guide_tree -i input_fasta -o output_fasta
```
Specify a maximum subtree size to reduce memory usage when memory usage is the major concern.
```bash
./twilight -t guide_tree -i input_fasta -o output_fasta -m maximum_subtree_size
```
To merge multiple MSAs, please move all MSA files into a folder.
```bash
./twilight -f folder -o output_fasta
```
### <a name="iterative"></a> Iterative mode

The TWILIGHT iterative mode provides a snakemake workflow to estimate guide trees using external tools. Please first install snakemake `pip install snakemake` and tools you would like to use and configure the executable path in `workflow/config.yaml`.
- Initial guide tree: MashTree, MAFFT(PartTree)
- Subsequent iterations: FastTree, IQ-Tree, RAxML

Update the `workflow/config.yaml` file to include the path to your input sequence file, tree inference tools, the number of cpu threads, and other parameters.

Dry run to check the workflow.
```bash
snakemake -n
```
Run the TWILIGHT iterative mode.
```bash
snakemake
```
## <a name="o_example"></a> Other example commands
Align divergent sequences.
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y
```
Reduce memory usage.
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -p y -d RNASim_temp -m 200
```
Align short-branched sequences.
```bash
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n
```
Specify a user-defined scoring matrix.
```bash
./twilight -t ../dataset/sars_20.nwk -i ../dataset/sars_20.fa -o sars_20.aln -r 1 -p n -x ../dataset/substitution.txt --gap-open -20 --gap-extend -4
```
## <a name="cite"></a> Citation
TBA.
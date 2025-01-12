# TWILIGHT: Tall and WIde aLIGnments at High Throughput

## Table of Contents
- [Overview](#overview)
- [Quick start](#start)
  - [Unix commands](#unix)
  - [Docker](#docker)
- [TWILIGHT modes](#mode)
  - [Default](#default)
  - [Iterative](#iterative)
- [Citation](#cite)


<br>

## <a name="overview"></a> Overview

TWILIGHT is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration. 

TWILIGHT accepts unaligned sequences in FASTA format and an optional input guide tree in Newick format to generate output alignments in FASTA format. When a guide tree is unavailable, users can utilize the provided snakemake workflow to estimate trees using external tools.

## <a name="start"></a> Quick start

### <a name="unix"></a> Unix commands (requires sudo access to install dependant libraries)
#### Dependent libraries
```
sudo apt-get install libboost-all-dev # BOOST
sudo apt-get install libtbb2 # TBB 2.0
```
#### Build TWILIGHT
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
mkdir build && cd build
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz
tar -xvzf 2019_U9.tar.gz
cmake  -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake  ..
make
```
### <a name="docker"></a> Use Docker
Docker currently only supports default mode.
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
docker build -t twilight_image .
docker run -it twilight_image
```
## <a name="mode"></a> TWILIGHT modes
### <a name="default"></a> Default mode
See Help for more details
```
./twilight -h
```
Default
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
### <a name="iterative"></a> Iterative mode
Please install snakemake `pip install snakemake` and tools you would like to use and configure the executable path in `workflow/config.yaml`.
- Initial guide tree: MashTree, MAFFT(PartTree)
- Subsequent iterations: FastTree, IQ-Tree, RAxML

Update the `workflow/config.yaml` file to include the path to your input sequence file, tree inference tools, the number of cpu threads, and other parameters.

Dry run to see the workflow
```
snakemake -n
```
Run iterative mode
```
snakemake
```
## <a name="cite"></a> Citation
TBA.
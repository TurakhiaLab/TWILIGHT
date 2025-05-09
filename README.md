



<div align="center">
    
# TWILIGHT: Tall and Wide Alignments at High Throughput

[license-badge]: https://img.shields.io/badge/License-MIT-yellow.svg 
[license-link]: https://github.com/TurakhiaLab/TWILIGHT/blob/main/LICENSE

[![License][license-badge]][license-link]
[<img src="https://img.shields.io/badge/Build with-CMake-green.svg?logo=snakemake">](https://cmake.org)
[<img src="https://img.shields.io/badge/Made with-Snakemake-aquamarine.svg?logo=snakemake">](https://snakemake.readthedocs.io/en/v7.19.1/index.html)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/twilight/README.html)

<div align="center">
  <img src="docs/images/logo.png" width="800"/>
</div>

</div>

## Table of Contents
- [Introduction](#intro) ([Wiki](https://turakhia.ucsd.edu/TWILIGHT/))
- [Installation](#install)
  - [Using Installation Script](#script)
  - [Using Conda](#conda)
  - [Using Dockerfile](#docker)
- [Run TWILIGHT](#run)
  - [Default mode](#default)
  - [Iterative mode](#iterative)
- [Contributions](#contribution)
- [Citing TWILIGHT](#cite)

<br>

## <a name="intro"></a> Introduction

TWILIGHT (**T**all and **Wi**de A**lig**nments at **H**igh **T**hroughput) is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration.

By default, TWILIGHT requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate the output alignment in FASTA format (Fig. 1a). When a guide tree is unavailable, users can utilize the iterative mode, which provides a Snakemake workflow to estimate guide trees using external tools (Fig 1b.).

TWILIGHT adopts the progressive alignment algorithm (Fig. 1c) and employs tiling strategies to band alignments (Fig. 1e). Combined with a divide-and-conquer technique (Fig. 1a), a novel heuristic dealing with gappy columns (Fig. 1d) and support for GPU acceleration (Fig. 1f), TWILIGHT demonstrates exceptional speed and memory efficiency.

<div align="center">
    <div><b>Figure 1: Overview of TWILIGHT alogorithm</b></div>
    <img src="docs/images/overview.png" width="800"/>
</div>

## <a name="install"></a> Installation
### <a name="script"></a> Using installation script (requires sudo access and only for Ubuntu)  

This has been tested only on Ubuntu. Users without sudo access are advised to install TWILIGHT via [Conda](#conda). For users on other platforms or operating systems, please refer to the [Docker](#docker) section below for installation instructions.

**Step 0:** Dependencies  
- Git: `sudo apt install -y git` 
- Conda: Optional, for iterative mode  

**Step 1:** Clone the repository
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
```
**Step 2:** Install dependencies (requires sudo access)

Skip this step if the below libraries are already installed.
```bash
- wget
- build-essential 
- cmake 
- libboost-all-dev 
- libtbb2 
- protobuf-compiler
```
Otherwise,
```bash
bash ./install/installDependencies.sh
```
**Step 3:** Install TWILIGHT

If CUDA-capable GPUs are detected, the GPU version will be built; otherwise, the CPU version will be used.
```bash
bash ./install/buildTWILIGHT.sh
```
**Step 4:** Enter `build` directory and run TWILIGHT
```bash
cd build
./twilight --help
```
**Step 5 (optional):** Install TWILIGHT iterative mode

**Step 5.1** Create and activate a Conda environment (ensure Conda is installed first)
```bash
cd ../ # Return to TWILIGHT home directory
conda create -n twilight -y
conda activate twilight
```
**Step 5.2** Install Snakemake and tree inference tools
```bash
bash ./install/installIterative.sh
```

### <a name="conda"></a> Using Conda
TWILIGHT is currently available for installation via Conda on the linux/amd64 platform only. Users on other platforms please refer to the [Docker](#docker) section.

**Step 1:** Create and activate a Conda environment (ensure Conda is installed first)
```bash
conda create -n twilight -y
conda activate twilight
```
**Step 2:** Install TWILIGHT
```bash
conda install bioconda::twilight
```
**Step 3 (optional):** Install TWILIGHT iterative mode
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
bash ./install/installIterative.sh
```


### <a name="docker"></a> Using Dockerfile
The Dockerfile installed all the dependencies and tools for TWILIGHT default/iterative mode. 

**Step 0:** Dependencies
- Git: `sudo apt install -y git` 
- Docker

**Step 1:** Clone the repository
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
```
**Step 2:** Build a docker image

CPU version
```bash
cd docker/cpu
docker build -t twilight .
```
GPU version (using nvidia/cuda as base image)
```bash
cd docker/gpu
docker build -t twilight .
```
**Step 3:** Build and run docker container

CPU version
```bash
docker run --platform=linux/amd64 -it twilight
```
GPU version  
Since the prebuilt version included in the Docker image is for CPU only, please rebuild TWILIGHT within the container to enable GPU support.
```bash
docker run --platform=linux/amd64 --gpus all -it twilight
bash ./install/installTWILIGHT.sh # Within the container
```
**Step 4:** Enter `build` directory and run TWILIGHT
```bash
cd build
./twilight -h
```


## <a name="run"></a> Run TWILIGHT
### <a name="default"></a> Default Mode
For more information about TWILIGHT's options and instructions, see [wiki](https://turakhia.ucsd.edu/TWILIGHT/) or *Help* for more details. 
```bash
cd build
./twilight -h
```
#### Default Configuration
Usage syntax
```bash
./twilight -t <path to tree file> -i <path to sequence file> -o <path to output file>
```
Example
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln
```
#### Divide-and-Conquer Method
TWILIGHT divides tree into subtrees with at most *m* leaves, which is specified by the user, and align subtrees sequentially to reduce the CPU’s main memory usage.
Usage syntax
```bash
./twilight -t <path to tree file> -i <path to sequence file> -o <path to output file> -m <maximum subtree size>
```
Example
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -m 200
```
#### Merge Multiple MSA Files
To merge multiple MSAs, please move the MSA files into a folder.  
Usage syntax
```bash
./twilight -f <path to the folder> -o <path to output file>
```
Example
```bash
./twilight -f ../dataset/RNASim_subalignments/ -o RNASim.aln
```
### <a name="iterative"></a> Iterative Mode
TWILIGHT iterative mode provides a Snakemake workflow to estimate guide trees using external tools. 

Options for tree inference tools:
- Initial guide tree: `parttree`, `maffttree`, `mashtree`
- Subsequent iterations: `fasttree`, `iqtree`, `raxml`

**Step 1:** Enter `workflow` directory
```bash
cd workflow
```
**Step 2:** See [wiki](https://turakhia.ucsd.edu/TWILIGHT/) for more details for the configurations. For people who install TWILIGHT via Conda, please replace the executable path `"../build/twilight"` with `"twilight"` in `config.yaml`.

**Step 3:** Run TWILIGHT iterative mode.  
Usage syntax
```bash
snakemake --cores [num threads] --config SEQ=[sequence] OUT=[output] DIR=[directory] ITER=[iterations] INITTREE=[tree method] ITERTREE=[tree method] OUTTREE=[tree] GETTREE=[yes/no]
```
Example  
- Using default configurations
```bash
snakemake --cores 8 --config SEQ=../dataset/RNASim.fa OUT=RNASim.aln DIR=tempDir
```
- Specifying all command line options
```bash
snakemake --cores 8 --config SEQ=../dataset/RNASim.fa OUT=RNASim.aln DIR=tempDir ITER=2 INITTREE=maffttree ITERTREE=raxml OUTTREE=RNASim.tree GETTREE=yes
```
##  <a name="contribution"></a> Contributions
We welcome contributions from the community to enhance the capabilities of **TWILIGHT**. If you encounter any issues or have suggestions for improvement, please open an issue on [TWILIGHT GitHub page](https://github.com/TurakhiaLab/TWILIGHT). For general inquiries and support, reach out to our team.

##  <a name="cite"></a> Citing TWILIGHT
TBA.

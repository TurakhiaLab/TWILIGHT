



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
  - [Summary](#summary) 
  - [Using Conda](#conda)
  - [Using Install Script](#script)
  - [Using Dockerfile](#docker)
- [Run TWILIGHT](#run)
  - [Default mode](#default)
  - [Iterative mode](#iterative)
- [Contributions](#contribution)
- [Citing TWILIGHT](#cite)

<br>

## <a name="intro"></a> Introduction

TWILIGHT (**T**all and **Wi**de A**lig**nments at **H**igh **T**hroughput) is a tool designed for ultrafast and ultralarge multiple sequence alignment. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration.

By default, TWILIGHT requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate the output alignment in FASTA format (Fig. 1a, <a name="default"></a>**default mode**). When a guide tree is unavailable, TWILIGHT provides a Snakemake workflow to estimate guide trees using external tools (Fig 1b, <a name="iter"></a>**iterative mode**).

TWILIGHT adopts the progressive alignment algorithm (Fig. 1c) and employs tiling strategies to band alignments (Fig. 1e). Combined with a divide-and-conquer technique (Fig. 1a), a novel heuristic dealing with gappy columns (Fig. 1d) and support for GPU acceleration (Fig. 1f), TWILIGHT demonstrates exceptional speed and memory efficiency.

<div align="center">
    <div><b>Figure 1: Overview of TWILIGHT alogorithm</b></div>
    <img src="docs/images/overview.png" width="800"/>
</div>

## <a name="install"></a> Installation
### <a name="summary"></a> Installation summary (choose your installation method)

TWILIGHT offers multiple installation methods for different platforms and hardware setups:
- Conda is recommended for most users needing the [default mode](#default) and *partial* [iterative mode](#iter) support, as some tree tools may be unavailable on certain platforms.
- Install script is required for AMD GPU support.
- Docker (built from the provided Dockerfile) is recommended for full support for [iterative mode](#iter).

| Platform / Setup       | [Conda](#conda) | [Script](#script) | [Docker](#docker) |
|------------------------|-----------------|-------------------|-------------------|
| Linux (x86_64)         | ✅               | ✅                | ✅                |
| Linux (aarch64)        | ✅               | ✅                | 🟡                |
| macOS (Intel Chip)     | ✅               | ✅                | ✅                |
| macOS (Apple Silicon)  | ✅               | ✅                | 🟡                |
| NVIDIA GPU             | ✅               | ✅                | ✅                |
| AMD GPU                | ❌               | ✅                | ❌                |

🟡 The Docker image is currently built for the `linux/amd64` platform. While it can run on `arm64` systems (e.g., Apple Silicon or Linux aarch64) via emulation, this may lead to reduced performance.

⚠️ To enable **GPU support**, the appropriate GPU drivers must be installed on the host system. This applies to all installation methods (Installation script, Conda, and Docker). The CUDA toolkits and libraries are included in Conda and Docker setups, but must be installed manually when using the installation script.  

### <a name="conda"></a> Using Conda
TWILIGHT is available on multiple platforms via Conda. See [TWILIGHT Bioconda Page](https://anaconda.org/bioconda/twilight) for details.  

**Step 1:** Create and activate a Conda environment (ensure Conda is installed first)
```bash
conda create -n twilight -y
conda activate twilight
# Set up channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
# Install TWILIGHT
conda install bioconda::twilight
```
**Step 2 (optional):** Install TWILIGHT iterative mode
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
bash ./install/installIterative.sh
```

### <a name="script"></a> Using installation script (requires sudo access if certain common libraries are not already installed)  

Users without sudo access are advised to install TWILIGHT via [Conda](#conda) or [Docker](#docker).

**Step 1:** Clone the repository
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
```
**Step 2:** Install dependencies (requires sudo access)

TWILIGHT depends on the following common system libraries, which are typically pre-installed on most development environments:
```bash
- wget
- build-essential 
- cmake 
- libboost-all-dev 
```
It also requires `libtbb-dev`, which is not always pre-installed on all systems. For users who do not have sudo access and are missing **only** `libtbb-dev`, our script builds and installs TBB from source in the local user environment, with **no sudo access required**.

For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
```bash
sudo apt install -y wget build-essential libboost-all-dev cmake libtbb-dev
```
For Mac users, install dependencies using **Homebrew**:
```bash
xcode-select --install # if not already installed
brew install wget boost cmake tbb
```

**Step 3:** Build TWILIGHT

Our build script automatically detects the best available compute backend **(CPU, NVIDIA GPU, or AMD GPU)** and builds TWILIGHT accordingly. Alternatively, users can manually specify the desired target platform.

Automatic build:
```bash
bash ./install/buildTWILIGHT.sh
```
Build for a specific platform:
```bash
bash ./install/buildTWILIGHT.sh cuda # For NVIDIA GPUs
bash ./install/buildTWILIGHT.sh hip  # For AMD GPUs
```
**Step 4:** The TWILIGHT executable is located in the `bin` directory and can be run as follows:
```bash
cd bin
./twilight --help
```
**Step 5 (optional)** Install TWILIGHT iterative mode (ensure Conda is installed first)
```bash
# Create and activate a Conda environment 
conda create -n twilight -y
conda activate twilight
# Install Snakemake and tree inference tools
bash ./install/installIterative.sh
```

### <a name="docker"></a> Using Dockerfile
The Dockerfile installed all the dependencies and tools for TWILIGHT default/iterative mode. 

**Step 1:** Clone the repository
```bash
git clone https://github.com/TurakhiaLab/TWILIGHT.git
cd TWILIGHT
```
**Step 2:** Build a docker image (ensure Docker is installed first)

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
**Step 3:** Start and run docker container

CPU version
```bash
docker run --platform=linux/amd64 -it twilight
```
GPU version  
```bash
docker run --platform=linux/amd64 --gpus all -it twilight
```
**Step 4:** Run TWILIGHT
```bash
cd bin
./twilight -h
```


## <a name="run"></a> Run TWILIGHT
### <a name="default"></a> Default Mode
For more information about TWILIGHT's options and instructions, see [wiki](https://turakhia.ucsd.edu/TWILIGHT/) or *Help* for more details. 
```bash
cd bin
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
**Step 2:** See [wiki](https://turakhia.ucsd.edu/TWILIGHT/) for more details for the configurations. For users who install TWILIGHT via Conda, please replace the executable path `"../bin/twilight"` with `"twilight"` in `config.yaml`. Feel free to switch to a more powerful tree tool if available, such as replacing `"raxmlHPC"` with `"raxmlHPC-PTHREADS-AVX2"` for better performance. 

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

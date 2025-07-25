# <b>Welcome to TWILIGHT Wiki</b>
<div align="center">
    <img src="images/logo.png"/>
</div>

## <b>TWILIGHT Video Tutorial</b>

<iframe width="1000" height="600" src="https://www.youtube.com/embed/bn8N_Y-_g_A?si=unLXo8xNFAiNDa-K" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

## <b>Introduction</b> 
### <b>Overview</b><a name="overview"></a>
TWILIGHT (**T**all and **Wi**de A**lig**nments at **H**igh **T**hroughput) is an innovative MSA tool designed to leverage massive parallelism for efficient handling of tall and wide alignment tasks. It is able to scale to millions of long nucleotide sequences (>10000 bases). TWILIGHT can run on CPU-only platforms (Linux/Mac) or take advantage of CUDA-capable GPUs for further acceleration.

By default, TWILIGHT requires an unaligned sequence file in FASTA format and an input guide tree in Newick format to generate the output alignment in FASTA format (Fig. 1a, see [**TWILIGHT Default Mode**](#default) for more details). When a guide tree is unavailable, users can utilize the iterative mode, which provides a Snakemake workflow to estimate guide trees using external tools (Fig. 1b, see [**TWILIGHT Iterative Mode**](#iterative) for more details).

TWILIGHT adopts the **progressive alignment algorithm (Fig. 1c)** and employs **tiling strategies to band alignments (Fig. 1e)**. Combined with a **divide-and-conquer technique (Fig. 1a)**, **a novel heuristic dealing with gappy columns (Fig. 1d)** and **support for GPU acceleration (Fig. 1f)**, TWILIGHT demonstrates exceptional speed and memory efficiency.  

<div align="center">
    <div><b>Figure 1: Overview of TWILIGHT algorithm</b></div>
    <img src="images/overview.png" width="1000"/>
</div>

### <b>Key Features</b>

#### <b>Divide-and-Conquer Technique</b>  
TWILIGHT's ***divide-and-conquer technique (Fig. 1a)*** starts by dividing the initial guide tree into smaller subtrees with at most *m* leaves, which controls memory usage by limiting the number of sequences processed at once. TWILIGHT sequentially iterates over each subtree, loading its corresponding sequences into memory, and computing the *subalignment* for the subtree. After all *subalignments* are computed, TWILIGHT merges *subalignments* with transitivity merger algorithm of PASTA (Mirarab et al., 2015) or progressive alignment.  

#### <b>Progressive Alignment and Scheduling</b>
TWILIGHT performs progressive alignment based on the guide tree structure, first
aligning the most similar sequences or sequence groups and then progressively adding less similar sequences or groups. The ***scheduling algorithm (Fig. 1c)*** initializes the leaf nodes and then performs a post-order traversal to set the alignment order for internal nodes. This approach allows certain child nodes (e.g., a, b, d, f) to be aligned in parallel, while others (e.g., c, e, g) are processed sequentially, ensuring efficient parallel processing while respecting tree dependencies.  

#### <b>Removing gappy columns</b>  
As more sequences are added, the alignment tends to grow longer, often with many gap-dominated columns, leading to increased computational time. To improve efficiency, TWILIGHT ***excludes gappy columns (Fig. 1d)***—where gaps exceed a specified threshold—before the alignment step, resulting in shorter profiles and reducing computational overhead. After the alignment step, TWILIGHT restores the excluded gappy columns. If gappy columns from both profiles align at the same position, they are merged, reducing the overall alignment length. Otherwise, a gappy column from one profile is aligned with a new gap introduced in the other
profile.  

#### <b>Banded Alignment and Tiling Strategy</b>
TWILIGHT employs a modified X-Drop algorithm to control memory usage in pairwise profile alignment, ensuring linear scaling with sequence length. However, storing traceback pointers in limited GPU shared memory still poses challenges. To address this, TWILIGHT integrates the ***TALCO tiling algorithm  (Walia et al., 2024) (Fig. 1e)***, which limits traceback storage to a constant memory usage by using convergence pointers.  

#### <b>Parallelization Techniques</b>
TWILIGHT employs parallelization techniques on both CPU and GPU to maximize efficiency. On CPUs, it is implemented in C++ with Threading Building Blocks (TBB) to exploit thread-level parallelism, processing independent alignments in parallel while also parallelizing profile calculations and sequence updates. On GPUs, TWILIGHT utilizes ***three levels of parallelism (Fig. 1f)***: multi-GPU parallelism distributes alignment tasks into batches and sends to multiple GPUs, inter-alignment parallelism assigns each GPU thread block to a pairwise profile alignment, and intra-alignment parallelism processes multiple cells in the dynamic programming matrix wavefront simultaneously.  






<a name="install"></a>
## <b>Installation Methods</b>

### **Installation summary (choose your installation method)**

TWILIGHT offers multiple installation methods for different platforms and hardware setups:

- Conda is recommended for most users needing the [default mode](#overview) and *partial* [iterative mode](#overview) support, as some tree tools may be unavailable on certain platforms.
- Install script is required for AMD GPU support.
- Docker (built from the provided Dockerfile) is recommended for full support for [iterative mode](#overview).

| Platform / Setup       | [Conda](#conda) | [Script](#script) | [Docker](#docker) |
|------------------------|-----------------|-------------------|-------------------|
| Linux (x86_64)         | ✅               | ✅                | ✅                |
| Linux (aarch64)        | ✅               | ✅                | 🟡                |
| macOS (Intel Chip)     | ✅               | ✅                | ✅                |
| macOS (Apple Silicon)  | ✅               | ✅                | 🟡                |
| NVIDIA GPU             | ✅               | ✅                | ✅                |
| AMD GPU                | ❌               | ✅                | ❌                |

!!!Note
    🟡 The Docker image is currently built for the `linux/amd64` platform. While it can run on `arm64` systems (e.g., Apple Silicon or Linux aarch64) via emulation, this may lead to reduced performance.

!!!Note
    ⚠️ To enable **GPU support**, the appropriate GPU drivers must be installed on the host system. This applies to all installation methods (Installation script, Conda, and Docker). The CUDA toolkits and libraries are included in Conda and Docker setups, but must be installed manually when using the installation script.  


### **Using Conda** <a name=conda></a>

TWILIGHT is available on multiple platforms via Conda. See [TWILIGHT Bioconda Page](https://anaconda.org/bioconda/twilight) for details.  

0. <a name=installconda></a> Install Conda (if not installed)
    ```bash
    # Replace <PLATFORM> with the actual platform of your system
    # See https://repo.anaconda.com/miniconda/ for details
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-<PLATFORM>.sh
    chmod +x Miniconda3-latest-<PLATFORM>.sh
    ./Miniconda3-latest-<PLATFORM>.sh
    export PATH="$HOME/miniconda3/bin:$PATH"
    # Reload shell configuration
    source ~/.bashrc  # For Bash (Most Linux)
    source ~/.zshrc   # For Zsh (Default in Mac)
    ```
1. Create and activate a Conda environment
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
2. Install TWILIGHT iterative mode
    ```bash
    git clone https://github.com/TurakhiaLab/TWILIGHT.git
    cd TWILIGHT
    bash ./install/installIterative.sh
    ```

### **Using installation script (requires sudo access)** <a name=script></a>
Users without sudo access are advised to install TWILIGHT via [Conda](#conda) or [Docker](#docker).

1. Clone the repository
    ```bash
    git clone https://github.com/TurakhiaLab/TWILIGHT.git
    cd TWILIGHT
    ```
2. Install dependencies (requires sudo access)  
    TWILIGHT depends on the following common system libraries, which are typically pre-installed on most development    environments:
    ```bash
    - wget
    - build-essential 
    - cmake 
    - libboost-all-dev 
    ```
    It also requires `libtbb-dev`, which is not always pre-installed on all systems. For users who do not have sudo     access and are missing **only** `libtbb-dev`, our script builds and installs TBB from source in the local user  environment, with **no sudo access required**.

    For Ubuntu users with sudo access, if any of the required libraries are missing, you can install them with:
    ```bash
    sudo apt install -y wget build-essential libboost-all-dev cmake libtbb-dev
    ```
    For Mac users, install dependencies using **Homebrew**:
    ```bash
    xcode-select --install # if not already installed
    brew install wget boost cmake tbb
    ```
3. Build TWILIGHT  
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
4. The TWILIGHT executable is located in the `bin` directory and can be run as follows:
    ```bash
    cd bin
    ./twilight --help
    ```
5. (optional) Install TWILIGHT iterative mode (ensure [Conda is installed](#installconda) first)
    ```bash
    # Create and activate a Conda environment 
    conda create -n twilight -y
    conda activate twilight
    # Install Snakemake and tree inference tools
    bash ./install/installIterative.sh
    ```

### **Using Dockerfile** <a name=docker></a> 
The Dockerfile installed all the dependencies and tools for TWILIGHT default/iterative mode.  

1. Clone the repository  
    ```bash
    git clone https://github.com/TurakhiaLab/TWILIGHT.git
    cd TWILIGHT
    ```
2. Build the docker image  
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
3. Start and run docker container   
    CPU version
    ```bash
    docker run --platform=linux/amd64 -it twilight
    ```
    GPU version  
    ```bash
    docker run --platform=linux/amd64 --gpus all -it twilight
    ```
4. Run TWILIGHT
    ```bash
    ./twilight --help
    ```

## <a name="Run TWILIGHT"></a> **Run TWILIGHT**

### <a name="default"></a> **Default Mode**  

#### Functionalities
<div name="table1" align="center"> <b>Table 1:</b> List of functionalities supported by TWILIGHT </div>

| **Option**                       | **Description**                                                                                                      |
|----------------------------------|----------------------------------------------------------------------------------------------------------------------| 
| `-t`, `--tree`                   | Input guide tree file path (Newick format).                                                                          |
| `-i`, `--sequences`              | Input unaligned sequences path (FASTA format).                                                                       |
| `-f`, `--files`                  | Path to the directory containing all MSA files (FASTA format) to be merged.                                          |
| `-o`, `--output`                 | Output MSA file path.                                                                                                |
| `-C`, `--cpu`                    | Number of CPU cores. Default: all available cores.                                                                   |
| `-m`, `--max-subtree`            | Maximum number of leaves in a subtree. Used when divide-and-conquer method is applied.                               |
| `-g`, `--max-group`              | Maximum number of leaves in a group within a subtree. Groups are merged using transitivity merger.                   |
| `-d`, `--temp-dir`               | Directory for storing temporary files. Used when `-f` or `-m` is specified.                                          |
| `-r`, `--remove-gappy`           | Threshold for removing gappy columns. Set to 1 to disable this feature. Default: 0.95. <br> It is recommended to set a higher threshold (i.e. 0.999) or to 1 when the alignment is known to be less gappy.                                           |
| `--type`               | Data type. `n`: nucleotide, `p`: protein. Will be automatically inferred if not provided.                                       |
| `-w`, `--wildcard`               | Treat unknown or ambiguous bases as wildcards and align them to usual letters.                                       |
| `-v`, `--verbose`                | Print out every detail process.                                                                                      | 
| `-p`, `--pgsop`                  | Position-specific gap open penalty. `y` for enable and `n` for disable. Default: `y`| 
| `-x`, `--matrix`                 | User-specified substitution matrix path.                                                                             |
| `-k`, `--keep-temp`              | Keep the temporary directory.                                                                                        |
| `-s`, `--sum-of-pairs-score`     | Calculate the sum-of-pairs-score after alignment.                                                                    |
| `--type`                         | Data type. `n` for nucleotide sequences and `p` for protein sequences. Will be automatically inferred if not provided.                          |
| `--no-align-gappy`               | Do not align gappy columns. This will create a longer MSA (larger file).                                             |
| `--length-deviation`               | Sequences whose lengths deviate from the average by more than the specified fraction will be deferred or excluded. Default: disabled. |
| `--max-ambig`               | Sequences with an ambiguous character proportion exceeding the specified threshold will be deferred or excluded. Default: 0.1.  |
| `--filter`               | Exclude sequences with high ambiguity or length deviation. Default: disabled. |
| `--merge`                        | Method to merge subtrees. `t` for Transitivity Merger and `p` for Pregressive Alignment. Default: `p`.               |
| `--match`                        | Match score. Default: 18.                                                                                            | 
| `--mismatch`                     | Mismatch penalty (transversions). Default: -8.                                                                       | 
| `--transition`                   | Transition score. Default: -4.                                                                                        |
| `--gap-open`                     | Gap-open penalty. Default: -50.                                                                                      |
| `--gap-extend`                   | Gap-extend penalty. Default: -5.                                                                                     |
| `--xdrop`                        | X-drop value (scale). The actual X-drop will be multiplied by the gap-extend penalty. Default: 600.                  |
| `--output-type`                  | Way to present MSA. `FASTA` for FASTA format and `CIGAR` for CIGAR-like compressed format. Default: `FASTA`.         |
| `--check`                  | Check the final alignment. Sequences with no legal alignment will be displayed.       |
| `-h`, `--help`                   | Print help messages.                                                                                                 |
| `--version`                   | Show TWILIGHT version.                                                                                                 |
| **Options only for GPU version**                                                                                                                        |
| `-G`, `--gpu`                    | Number of GPU. Default: all available GPUs.                                                                          |
| `--gpu-index`                    | Specify the GPU to be used, separated by commas. For example, `0,2,3` uses 3 GPUs with the index 0, 2 and 3.         |
| `--cpu-only`                     | Run the program only using CPU.                                                                                      |


!!!Note
    All files used for the examples below can be found in the `TWILIGHT/dataset`.

Enter into the build directory (assuming `$TWILIGHT_HOME` directs to the TWILIGHT repository directory)  
```bash
cd $TWILIGHT_HOME/bin
./twilight -h
```
#### **Default Configuration**  
Generate MSA by giving an unaligned sequence file and a guide tree file.  

* Usage syntax
```bash
./twilight -t <path to tree file> -i <path to sequence file> -o <path to output file>
```
* Example
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln
```
#### **Divide-and-Conquer Method**  
Devide tree into subtrees with at most *m* leaves and align them sequentially to reduce the CPU’s main memory usage.

* Usage syntax
```bash
./twilight -t <path to tree file> -i <path to sequence file> -o <path to output file> -m <maximum subtree size>
```
* Example
```bash
./twilight -t ../dataset/RNASim.nwk -i ../dataset/RNASim.fa -o RNASim.aln -m 200
```
#### **Merge Multiple MSA Files**
Assume a star-tree structure to merge multiple MSAs. Please move all MSA files to be merged into one folder.  

* Usage syntax
```bash
./twilight -f <path to the folder> -o <path to output file>
```
* Example
```bash
./twilight -f ../dataset/RNASim_subalignments/ -o RNASim.aln
```
### <a name="iterative"></a> **Iterative Mode**
TWILIGHT iterative mode estimate guide trees using external tools.  

#### **Command Line Configurations**
<div name="table2" align="center"> <b>Table 2:</b> List of command line configurations in TWILIGHT iterative mode</div>  

| **Configuration <br>(used with `--config`)**    | **Description and Options**                                                                                   |
|-----------------------------------------------------|---------------------------------------------------------------------------------------------------| 
| `TYPE`                            | Input sequence type, options: `n` for nucleotide or `p` for protein sequences. (required)    |                                                               
| `SEQ`                            | Input unaligned sequences path (FASTA format).                                                                       |
| `OUT`                            | Output MSA file path.                                                                                                |
| `DIR`                            | Directory for storing temporary files.                                                                               |
| `ITER`                           | Number of iterations.                                                                                                |
| `INITTREE`                       | Tree method for initial guide tree, options: `parttree`, `maffttree` or `mashtree`.                                  |
| `ITERTREE`                       | Tree method for subsequent iterations, options: `fasttree`, `iqtree` or `raxml`.                                      | 
| `GETTREE`                        | Estimate tree after final alignment, options: `yes` or `no`.                                                          |
| `OUTTREE`                        | Output Tree file path, used when `GETTREE=yes`.                                                                      |  
| **Snakemake Options**                                                                                                                                   |
| `-n`                             | Dry run, check workflow.                                                                                             |
| `--cores`                        | Number of cores. If you would like to use more than 32 cores, please also modify `num_threads` in `workflow/config.yaml`. |

!!!Note
    Default options in `workflow/config.yaml` will be used if the configuration is not specified. You can also modify `workflow/config.yaml` to set your options.
    
!!!Note
    For users who install TWILIGHT via Conda, please replace the executable path `"../bin/twilight"` with `"twilight"` in `config.yaml`. Feel free to switch to a more powerful tree tool if available, such as replacing `"raxmlHPC"` with `"raxmlHPC-PTHREADS-AVX2"` for better performance. 

#### **Run TWILIGHT iterative mode**
1. Enter into `workflow` directory
```bash
cd $TWILIGHT_HOME/workflow
```
2. (optional) Update the `workflow/config.yaml` file  
In addition to command line options, `workflow/config.yaml` provides additional options for functionalities in TWILIGHT and substitution models for tree inference tools. See `workflow/config.yaml` for more details.  
3. Run
    * Usage syntax
    ```bash
    snakemake --cores [num threads] --config TYPE=[n/p] SEQ=[sequence] OUT=[output] DIR=[directory] ITER=[iterations] INITTREE=[tree method] ITERTREE=[tree method] OUTTREE=[tree] GETTREE=[yes/no]
    ```
    * Example  
    Using default configurations  
    ```bash
    snakemake --cores 8 --config TYPE=n SEQ=../dataset/RNASim.fa OUT=RNASim.aln DIR=tempDir
    ```
    Specifying all command line options
    ```bash
    snakemake --cores 8 --config TYPE=n SEQ=../dataset/RNASim.fa OUT=RNASim.aln DIR=tempDir ITER=2 INITTREE=maffttree ITERTREE=raxml OUTTREE=RNASim.tree GETTREE=yes
    ```

## <b>Contributions</b>
We welcome contributions from the community to enhance the capabilities of TWILIGHT. If you encounter any issues or have suggestions for improvement, please open an issue on [TWILIGHT GitHub page](https://github.com/TurakhiaLab/TWILIGHT). For general inquiries and support, reach out to our team.

## <b>Citing TWILIGHT</b>
If you use the <b>TWILIGHT</b> in your research or publications, we kindly request that you cite the following paper:  

Yu-Hsiang Tseng, Sumit Walia, Yatish Turakhia, "<i>Ultrafast and ultralarge multiple sequence alignments using TWILIGHT</i>", Bioinformatics, Volume 41, Issue Supplement_1, July 2025, Pages i332–i341, doi: [10.1093/bioinformatics/btaf212](https://doi.org/10.1093/bioinformatics/btaf212)
# TWILIGHT Scripts & Pipelines Quick Guide

This directory contains utility scripts and the iterative refinement pipeline for **TWILIGHT**.

---

## Table of Contents
1. [Overview](#1-overview)
2. [Prerequisites & Dependencies](#2-prerequisites--dependencies)
3. [Quick Start](#3-quick-start)
4. [twilight_dipper_iterative.py Guide](#4-twilight_dipper_iterativepy-guide)
   - [Command Line Parameters](#command-line-parameters)
   - [Advanced Features & Flags](#advanced-features--flags)
5. [References](#5-references)

---

## 1. Overview

**`twilight_dipper_iterative.py`** is an automated iterative refinement engine for high-accuracy multiple sequence alignment. It alternates TWILIGHT accurate consistency-based progressive alignment (`-A 2` / `-A 3`) with ultra-fast DIPPER profile guide tree reconstruction, guided by objective score tracking and early stopping.

> [!NOTE]
> **Nucleotide Specificity**: The iterative pipeline is currently designed specifically for **NUCLEOTIDE** sequences (DNA and RNA) because upstream DIPPER currently supports nucleotide sequence guide tree reconstruction. Protein support hooks are present in the script and will be enabled upon the official upstream release of protein-enabled DIPPER.

---

## 2. Prerequisites & Dependencies

- **Python**: Version 3.8 or higher.
- **Python Libraries**: `numpy` (optional, recommended for `--trim-gappy-cols`).
- **TWILIGHT**: Compiled binary at `../bin/twilight` (automatically detected).
- **DIPPER**: Compiled binary (`dipper` for GPU or `dipper_cpu` for CPU) in your `PATH` or standard installation directory.
  - Source Repository: [https://github.com/TurakhiaLab/DIPPER](https://github.com/TurakhiaLab/DIPPER)

---

## 3. Quick Start

### Minimal One-Liner (Automatic Initial Tree via DIPPER)
```bash
python3 scripts/twilight_dipper_iterative.py \
    -i input_sequences.fasta \
    -o aligned_output.fasta
```

### Standard High-Accuracy Production Run (32 Threads)
```bash
python3 scripts/twilight_dipper_iterative.py \
    -i input_sequences.fasta \
    -o aligned_output.fasta \
    -A 2 \
    --max-iters 10 \
    --min-iters 3 \
    --patience 2 \
    --trim-gappy-cols 0.99 \
    --threads 32
```

---

## 4. `twilight_dipper_iterative.py` Guide

### Command Line Parameters

> [!TIP]
> Run `python3 scripts/twilight_dipper_iterative.py --help` for the complete list of options and descriptions.

| Parameter | Default | Description |
| :--- | :---: | :--- |
| **`-i, --input`** | *Required* | Path to input unaligned sequences (FASTA). |
| **`-o, --output`** | *Required* | Path for output aligned MSA (FASTA). |
| **`-t, --tree`** | `None` | Path to initial guide tree (Newick). If omitted, DIPPER generates one automatically from raw sequences. |
| **`-A, --acc-mode`** | `2` | TWILIGHT accuracy mode: `2` (Pairwise consistency) or `3` (Transitive consistency / TCS). |
| **`-c, --threads`** | `16` | Number of CPU threads for both TWILIGHT and DIPPER. |
| **`--max-iters`** | `10` | Maximum number of refinement iterations. |
| **`--min-iters`** | `3` | Minimum iterations to execute before early stopping can trigger. |
| **`--patience`** | `2` | Stop after $N$ consecutive iterations without score improvement. |
| **`--min-delta`** | `0.001` | Minimum relative score improvement (0.1%) to count as progress. |
| **`--trim-gappy-cols`** | `0.0` | Trim columns with gap proportion $\ge$ threshold from MSA before passing to DIPPER (e.g. `0.99` for $\ge 99\%$ gaps; `0.0` disables trimming). |
| **`--keep-tree`** | `True` | Save the best reconstructed guide tree alongside MSA as `<output>.tree`. |
| **`--keep-temp`** | `False` | Keep intermediate iteration files (`iter*.aln`, `iter*.tree`, `iter*.score`) for debugging and inspection. |
| **`--temp-dir`** | `None` | Custom directory for temporary iteration files (default: system temp directory). |
| **`--score-type`** | `1` | Scoring function for convergence: `1` (Consistency Objective Score) or `2` (WSP). |
| **`--sample-max-group`**| `200` | Maximum group size for consistency partitioning in TWILIGHT. |
| **`--refine-method`** | `upgma`| Algorithm for Layer 1 subtree refinement: `upgma` (default) or `nj`. (Subtree refinement is automatically enabled in accurate mode). |
| **`--twilight-bin`** | `None` | Custom path to the `twilight` executable (default: auto-detected as `../bin/twilight`). |
| **`--dipper-bin`** | `None` | Custom path to the `dipper` or `dipper_cpu` executable (default: auto-detected). |

---

### Advanced Features & Flags

#### A. Gappy Column Trimming (`--trim-gappy-cols`)
In long or highly divergent alignments, sparse columns with high gap ratios (e.g. $\ge 99\%$ gaps) slow down tree reconstruction without providing phylogenetic signal.
```bash
python3 scripts/twilight_dipper_iterative.py \
    -i sequences.fa -o out.aln \
    --trim-gappy-cols 0.99
```
*Note: Trimming is applied **only** to the temporary MSA fed into DIPPER. The exported final MSA preserves all original columns without truncation.*

#### B. Consistency Level Selection (`-A 2` vs `-A 3`)
- **`-A 2` (Pairwise Consistency)**:
  - Fast, memory-efficient consistency transformation.
  - Recommended for large-scale datasets ($>10,000$ sequences).
- **`-A 3` (Transitive Consistency Score - TCS)**:
  - Higher topological accuracy and SP-score.
  - Recommended for challenging benchmark datasets or divergent sequences.

#### C. Objective Score Function (`--score-type`)
- `--score-type 1` (Default): **Consistency Objective Score**. Measures the agreement of alignment columns against the local pairwise consistency library.
- `--score-type 2`: **WSP (Weighted Sum-of-Pairs)**.

#### D. Subtree Refinement Control
In TWILIGHT accurate mode (`-A 2` / `-A 3`), Layer 1 guide subtree refinement is automatically enabled by default to boost topological and alignment accuracy without significant runtime overhead. You can select the refinement algorithm:
- `--refine-method upgma` (Default, recommended): Fast UPGMA guide subtree refinement based on intra-group pairwise alignment distances.
- `--refine-method nj`: Neighbor-Joining guide subtree refinement.

#### E. Retaining Tree and Intermediate Files
- `--keep-tree`: Enabled by default. Automatically exports the guide tree corresponding to the highest-scoring alignment as `<output>.tree`.
- `--keep-temp`: Preserves all per-iteration alignments and trees in `--temp-dir` for fine-grained convergence analysis.

---

## 5. References

If you use DIPPER or this iterative pipeline, please cite:
- **DIPPER Paper**:  
  Sumit Walia, Zexing Chen, Yu-Hsiang Tseng, Yatish Turakhia, *"Ultrafast and Ultralarge Distance-Based Phylogenetics Using DIPPER"*, Nature Computational Science (2026).  
  DOI: [10.1038/s43588-026-01015-8](https://doi.org/10.1038/s43588-026-01015-8)
- **DIPPER Source Repository**:  
  [https://github.com/TurakhiaLab/DIPPER](https://github.com/TurakhiaLab/DIPPER)

# TWILIGHT C API Documentation

TWILIGHT provides a C-compatible library interface ([`include/twilight.h`](../include/twilight.h)) for embedding in-memory, tree-guided progressive multiple sequence alignment directly into C or C++.

> **Supported Scope**
> - Supports standard guide-tree progressive alignment (equivalent to CLI `-t <tree> -i <seqs>`), including `--accuracy-level` `1` (standard), `2` (pairwise consistency), and `3` (transitive consistency / TCS).
> - **Not supported:** Divide-and-conquer subtree partitioning (`-m / --max-subtree`) and sequence placement modes (`-a / --alignment`).

---

## Table of Contents
- [Build the Static Library](#build-the-static-library)
- [Include & Link in Your Project](#include--link-in-your-project)
- [Quick Start Example](#quick-start-example)
- [API Reference](#api-reference)
  - [Options (`tw_opt_t`)](#options-tw_opt_t)
  - [Input Sequences (`tw_seq_t`)](#input-sequences-tw_seq_t)
  - [Output MSA (`tw_msa_t`)](#output-msa-tw_msa_t)
  - [Functions & Status Codes](#functions--status-codes)

---

## Build the Static Library

Building TWILIGHT with CMake produces:
- Static library: `lib/libtwilight.a` (CMake target: `twilight_lib`)
- Public C header: `include/twilight.h`
- CLI executable: `bin/twilight`
- C example executable: `bin/c_api_example`

```bash
cmake -S . -B build
cmake --build build -j
```

---

## Include & Link in Your Project

### Option 1: CMake (Recommended)
```cmake
add_subdirectory(path/to/TWILIGHT)

add_executable(my_app main.c)
target_link_libraries(my_app PRIVATE twilight_lib)
```

### Option 2: Direct Compilation (GCC / Clang)
Compile your C source with `-Iinclude`, then link against `lib/libtwilight.a`, TBB, zlib, and the C++ standard library (add `-L/usr/local/cuda/lib64 -lcudart` if `libtwilight.a` was built with CUDA):

```bash
# 1. Compile C source
gcc -O2 -Iinclude -c main.c -o main.o

# 2a. Link (CPU build)
g++ main.o lib/libtwilight.a -ltbb -lz -lpthread -o my_app

# 2b. Link (CUDA build)
g++ main.o lib/libtwilight.a -ltbb -lz -lpthread -L/usr/local/cuda/lib64 -lcudart -o my_app
```

---

## Quick Start Example

```c
#include "twilight.h"
#include <stdio.h>

int main(void) {
    /* 1. Initialize options with TWILIGHT defaults */
    tw_opt_t opt;
    tw_opt_init(&opt);
    opt.accuracy_level = TW_ACC_PAIRWISE; /* 1: standard, 2: pairwise consistency, 3: TCS */
    opt.num_threads = 4;                  /* 0 = auto (all available CPU cores) */
    opt.quiet = 1;                        /* 1 = silent (default for C API), 0 = log to stderr */

    /* 2. Prepare guide tree (Newick string) and input sequences (non-owning pointers) */
    const char *newick = "((seq1:0.1,seq2:0.2):0.1,seq3:0.3);";
    tw_seq_t seqs[3] = {
        { "seq1", 0, "ACGTACGTACGT", 12 },
        { "seq2", 0, "ACGTTCGTACGT", 12 },
        { "seq3", 0, "ACGTACGACGT",  11 }
    };

    /* 3. Run in-memory alignment */
    tw_msa_t *msa = NULL;
    tw_status_t st = tw_align(&opt, newick, seqs, 3, &msa);
    if (st != TW_OK) {
        fprintf(stderr, "TWILIGHT error: %s\n", tw_strerror(st));
        return 1;
    }

    /* 4. Access aligned sequences via 2D pointer msa->seqs[i] */
    for (size_t i = 0; i < msa->n_seqs; ++i) {
        printf(">%s\n%s\n", msa->names[i], msa->seqs[i]);
    }

    /* 5. Release memory owned by TWILIGHT */
    tw_msa_free(msa);
    return 0;
}
```

A complete runnable program supporting both inline demos and FASTA/Newick file inputs is provided in [`c_api_example.c`](c_api_example.c):

```bash
# Run built-in demo & error handling checks
./bin/c_api_example

# Run on FASTA + Newick files (accuracy level 2, 8 threads)
./bin/c_api_example dataset/RNASim.nwk dataset/RNASim.fa out.aln 2 8
```

---

## API Reference

### Options (`tw_opt_t`)

Always call `tw_opt_init(&opt)` first to populate defaults (or pass `NULL` to `tw_align()` to use all defaults). You can also apply a preset via `tw_opt_set_preset(&opt, "fast" | "accurate" | "tcs")`.

| Field | Default | Description |
| :--- | :--- | :--- |
| `num_threads` | `0` | Number of CPU threads (`0` = all available cores) |
| `cpu_only` | `0` | `1` forces CPU execution; `0` uses GPU if built with CUDA and `accuracy_level == 1` |
| `quiet` | `1` | `1` suppresses all stderr output (including errors); `0` logs to stderr |
| `verbose` | `0` | `1` enables `[DEBUG]` logs when `quiet == 0` |
| `seq_type` | `TW_SEQ_AUTO` | `TW_SEQ_AUTO` (`0`), `TW_SEQ_NUCL` (`1`), or `TW_SEQ_PROTEIN` (`2`) |
| `reroot` | `1` | `1` auto-reroots tree for balanced parallelism; `0` keeps original root |
| `prune_tree` | `1` | `1` prunes tree leaves missing from `seqs` |
| `remove_gappy` | `0.95f` | Threshold in `(0, 1]` for removing gappy columns (`1.0f` disables) |
| `max_ambig` | `0.1f` | Max ambiguous character fraction in `(0, 1]` |
| `len_dev` | `0.0f` | Allowed length deviation from median (`0.0f` disables) |
| `min_len` / `max_len` | `0` / `0` | Sequence length bounds (`0` = disabled/unlimited) |
| `filter_low_qual` | `0` | `0` defers low-quality sequences; `1` excludes them from output |
| `wildcard` | `0` | `1` treats ambiguous characters as wildcards |
| `match` / `mismatch` / `transition` | `18` / `-8` / `-4` | Nucleotide substitution scores |
| `gap_open` / `gap_extend` / `gap_ends` | auto | Normalized gap penalties `<= 0` (positive values select type-specific defaults) |
| `scale` / `xdrop` | auto / `600.0f` | Normalization scale (`<= 0` = auto) and X-drop multiplier |
| `blosum` | `62` | Protein matrix: `43` (PFASUM43), `45`, `62`, `80`, or `99` (MIQS) |
| `accuracy_level` | `1` | `TW_ACC_STANDARD` (`1`), `TW_ACC_PAIRWISE` (`2`), `TW_ACC_TRANSITIVE` (`3`) |
| `sample_mode` | `0` | `0`: auto, `1`: all-to-all, `2`: high rate, `3`: low rate |
| `sample_max_group` | `0` | Max consistency group size (`0` = auto) |
| `sample_rate` | `0.0f` | Hierarchical representative sample rate in `(0, 1)` (`0.0f` = auto) |
| `sample_top_centrality` | `0.0f` | Top-centrality quota in `(0, 1]` (`0.0f` = auto) |
| `consistency_weight` | `5.0f` | Multiplier for consistency bonus in `-A 2` / `-A 3` |
| `refine_subtree` | `-1` | `-1`: auto (ON for level 2/3), `0`: off, `1`: on |
| `refine_method` | `0` | `TW_REFINE_UPGMA` (`0`) or `TW_REFINE_NJ` (`1`) |
| `iterative_refine` | `-1` | `-1`: auto (ON for level 2/3), `0`: off, `1`: on |
| `refine_iters` / `refine_batch_size` | `1` / `4` | Iterative refinement iterations and batch size |
| `compute_score` | `0` | `0`: off, `1`: consistency objective score, `2`: WSP score |

### Input Sequences (`tw_seq_t`)

`tw_seq_t` is a non-owning view of caller memory:

```c
typedef struct {
    const char *name;  /* Sequence identifier (must match guide tree leaf names) */
    size_t name_len;   /* Length of name in bytes; 0 = compute via strlen(name) */
    const char *seq;   /* Unaligned sequence characters */
    size_t seq_len;    /* Length of seq in bytes; 0 = compute via strlen(seq) */
} tw_seq_t;
```

- Setting `name_len` or `seq_len` to `0` instructs TWILIGHT to compute the length via `strlen()`.
- Passing explicit lengths avoids `strlen()` overhead and supports non-null-terminated buffers (e.g., slices in a sequence database).

### Output MSA (`tw_msa_t`)

`tw_msa_t` is allocated by `tw_align()` and must be freed with `tw_msa_free()`:

```c
typedef struct {
    size_t n_seqs;       /* Number of aligned sequences in output */
    size_t aln_len;      /* Number of columns in the alignment */
    char **names;        /* names[i]: null-terminated sequence name */
    char **seqs;         /* seqs[i]: null-terminated aligned sequence of length aln_len */
    int *input_indices;  /* input_indices[i]: original index in input tw_seq_t[] */
    double score;        /* Alignment score if compute_score > 0, else 0.0 */
    ...
} tw_msa_t;
```

### Functions & Status Codes

- `void tw_opt_init(tw_opt_t *opt);`
- `tw_status_t tw_opt_set_preset(tw_opt_t *opt, const char *preset);`
- `tw_status_t tw_align(const tw_opt_t *opt, const char *newick_tree, const tw_seq_t *seqs, size_t n_seqs, tw_msa_t **out_msa);`
- `void tw_msa_free(tw_msa_t *msa);`
- `const char *tw_strerror(tw_status_t status);`
- `const char *tw_version(void);`

Return codes (`tw_status_t`): `TW_OK` (`0`), `TW_ERR_INVALID_ARG` (`-1`), `TW_ERR_EMPTY_INPUT` (`-2`), `TW_ERR_TREE_PARSE` (`-3`), `TW_ERR_IO` (`-4`), `TW_ERR_OOM` (`-5`), `TW_ERR_UNSUPPORTED` (`-6`), `TW_ERR_INTERNAL` (`-7`).

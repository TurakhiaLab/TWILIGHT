#ifndef TWILIGHT_H
#define TWILIGHT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * TWILIGHT C API
 *
 * Supported scope:
 *   - Guide-tree progressive multiple sequence alignment (equivalent to
 *     CLI `-t <tree> -i <seqs>`, supporting accuracy levels 1, 2, and 3).
 *   - Divide-and-conquer subtree partitioning (`-m / --max-subtree`) and
 *     placement modes (`-a / --alignment`) are NOT supported by this API.
 * ========================================================================= */

/* Status / error codes returned by TWILIGHT C API functions */
typedef enum {
    TW_OK              =  0,
    TW_ERR_INVALID_ARG = -1,  /* Invalid option value or NULL pointer */
    TW_ERR_EMPTY_INPUT = -2,  /* No valid sequences provided or matched */
    TW_ERR_TREE_PARSE  = -3,  /* Failed to parse Newick guide tree */
    TW_ERR_IO          = -4,  /* I/O failure */
    TW_ERR_OOM         = -5,  /* Memory allocation failure */
    TW_ERR_UNSUPPORTED = -6,  /* Unsupported mode or option combination */
    TW_ERR_INTERNAL    = -7   /* Unexpected internal error */
} tw_status_t;

/* Sequence alphabet type */
typedef enum {
    TW_SEQ_AUTO    = 0,  /* Auto-detect from input sequences */
    TW_SEQ_NUCL    = 1,  /* Nucleotide (DNA / RNA) */
    TW_SEQ_PROTEIN = 2   /* Protein (amino acids) */
} tw_seq_type_t;

/* Accuracy level (matches CLI -A / --accuracy-level) */
typedef enum {
    TW_ACC_STANDARD   = 1,  /* Traditional progressive sum-of-pairs */
    TW_ACC_PAIRWISE   = 2,  /* + Pairwise consistency */
    TW_ACC_TRANSITIVE = 3   /* + Transitive consistency (TCS) */
} tw_acc_level_t;

/* Layer 1 guide subtree refinement algorithm */
typedef enum {
    TW_REFINE_UPGMA = 0,
    TW_REFINE_NJ    = 1
} tw_refine_method_t;

/* =========================================================================
 * Options Struct (equivalent to minimap2's mm_mapopt_t)
 * Always initialize with tw_opt_init(&opt) before modifying fields.
 * ========================================================================= */
typedef struct {
    /* Hardware & Logging */
    int num_threads;             /* 0 = auto (all available CPU cores) */
    int cpu_only;                /* 0 = auto (use GPU if built with CUDA and -A 1), 1 = force CPU */
    int quiet;                   /* 1 = completely silent (default for C API), 0 = log to stderr */
    int verbose;                 /* 0 = normal, 1 = emit [DEBUG] messages (requires quiet == 0) */

    /* Data Type & Tree */
    int seq_type;                /* tw_seq_type_t: TW_SEQ_AUTO (default), NUCL, or PROTEIN */
    int reroot;                  /* 1 = auto reroot for parallelism (default), 0 = keep root */
    int prune_tree;              /* 1 = prune tree leaves not present in input seqs (default), 0 = off */

    /* Column & Sequence Filtering */
    float remove_gappy;          /* Threshold for removing gappy columns in (0, 1]; default 0.95f (1.0f = disable) */
    float max_ambig;             /* Max ambiguous character fraction in (0, 1]; default 0.1f */
    float len_dev;               /* Allowed length deviation from median (>= 0); default 0.0f (disabled) */
    int min_len;                 /* Minimum sequence length; default 0 (disabled) */
    int max_len;                 /* Maximum sequence length; default 0 (unlimited) */
    int filter_low_qual;         /* 0 = defer low quality seqs (default), 1 = exclude them from output */
    int wildcard;                /* 0 = off (default), 1 = treat ambiguous characters as wildcards */

    /* Scoring Parameters */
    float match;                 /* Nucleotide match score (default: 18.0f) */
    float mismatch;              /* Nucleotide transversion mismatch penalty (default: -8.0f) */
    float transition;            /* Nucleotide transition score (default: -4.0f) */
    float gap_open;              /* Normalized gap-open penalty <= 0; > 0 means auto (-2.5 nucl / -1.8 prot) */
    float gap_extend;            /* Normalized gap-extend penalty <= 0; > 0 means auto (-0.25 nucl / -0.18 prot) */
    float gap_ends;              /* End gap penalty <= 0; > 0 means same as gap_extend */
    float scale;                 /* Normalization scale > 0; <= 0 means auto (20.0 nucl / 600.0 prot) */
    float xdrop;                 /* X-drop scale multiplier (default: 600.0f) */
    int blosum;                  /* Protein matrix: 43 (PFASUM43), 45, 62 (default), 80, or 99 (MIQS) */

    /* Accurate / Consistency Mode (-A 2 / -A 3) */
    int accuracy_level;          /* tw_acc_level_t: 1 (default), 2, or 3 */
    int sample_mode;             /* 0: auto (default), 1: all-to-all, 2: high, 3: low */
    int sample_max_group;        /* 0 = auto; > 1 sets max group size */
    float sample_rate;           /* 0.0f = auto; (0, 1) sets hierarchical sample rate */
    float sample_top_centrality; /* 0.0f = auto; (0, 1] sets top-centrality quota */
    float consistency_weight;    /* Consistency bonus weight multiplier (default: 5.0f) */
    int refine_subtree;          /* -1 = auto (ON for level 2/3), 0 = off, 1 = on */
    int refine_method;           /* tw_refine_method_t: TW_REFINE_UPGMA (default) or TW_REFINE_NJ */
    int iterative_refine;        /* -1 = auto (ON for level 2/3), 0 = off, 1 = on */
    int refine_iters;            /* Max iterative refinement iterations (default: 1) */
    int refine_batch_size;       /* Iterative refinement batch size (default: 4) */

    /* Optional Alignment Score Computation */
    int compute_score;           /* 0 = off (default), 1 = consistency objective, 2 = WSP */
} tw_opt_t;

/* Initialize options struct with TWILIGHT defaults */
void tw_opt_init(tw_opt_t *opt);

/*
 * Apply a named preset to an initialized options struct.
 * Supported presets:
 *   "fast" / "standard"   -> accuracy_level = 1
 *   "accurate" / "mafft"  -> accuracy_level = 2
 *   "tcs"                 -> accuracy_level = 3
 * Returns TW_OK on success, TW_ERR_INVALID_ARG if preset is unknown.
 */
tw_status_t tw_opt_set_preset(tw_opt_t *opt, const char *preset);

/* =========================================================================
 * Input & Output Structures
 * ========================================================================= */

/* Non-owning view of a single unaligned input sequence in caller memory */
typedef struct {
    const char *name;  /* Sequence identifier (must match guide tree leaf names) */
    size_t name_len;   /* Length of name in bytes; 0 means compute via strlen(name) */
    const char *seq;   /* Unaligned sequence characters */
    size_t seq_len;    /* Length of seq in bytes; 0 means compute via strlen(seq) */
} tw_seq_t;

/*
 * Output Multiple Sequence Alignment owned by TWILIGHT.
 * Layout: contiguous backing buffers with 2D row pointers (seqs[i][j]).
 * Always release with tw_msa_free().
 */
typedef struct {
    size_t n_seqs;       /* Number of aligned sequences in output */
    size_t aln_len;      /* Length of each aligned sequence (number of MSA columns) */
    char **names;        /* Array of n_seqs null-terminated sequence names */
    char **seqs;         /* 2D pointer: seqs[i] is null-terminated aligned seq of length aln_len */
    int *input_indices;  /* Maps output row i -> original index in input tw_seq_t[] */
    double score;        /* Alignment score if opt->compute_score > 0, else 0.0 */

    /* Internal contiguous backing buffers; do not free or modify directly */
    char *_seq_buf;
    char *_name_buf;
} tw_msa_t;

/* =========================================================================
 * Core Alignment & Utility Functions
 * ========================================================================= */

/*
 * Perform tree-guided progressive multiple sequence alignment entirely in memory.
 *
 * Parameters:
 *   opt         - Options struct (NULL uses tw_opt_init defaults).
 *   newick_tree - Null-terminated Newick guide tree string (required).
 *   seqs        - Array of n_seqs unaligned input sequences (required).
 *   n_seqs      - Number of input sequences (must be >= 2).
 *   out_msa     - Receives newly allocated tw_msa_t* on TW_OK (caller frees via tw_msa_free).
 *
 * Returns TW_OK on success, or a negative tw_status_t error code on failure.
 */
tw_status_t tw_align(
    const tw_opt_t *opt,
    const char *newick_tree,
    const tw_seq_t *seqs,
    size_t n_seqs,
    tw_msa_t **out_msa
);

/* Free all memory owned by a tw_msa_t returned by tw_align() (safe on NULL) */
void tw_msa_free(tw_msa_t *msa);

/* Return a human-readable string describing a tw_status_t error code */
const char *tw_strerror(tw_status_t status);

/* Return TWILIGHT library version string */
const char *tw_version(void);

#ifdef __cplusplus
}
#endif

#endif /* TWILIGHT_H */

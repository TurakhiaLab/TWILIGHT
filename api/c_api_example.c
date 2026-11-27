#include "twilight.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Simple FASTA reader for testing the C API with benchmark files */
typedef struct {
    tw_seq_t *items;
    char **owned_names;
    char **owned_seqs;
    size_t count;
    size_t cap;
} fasta_list_t;

static void fasta_list_free(fasta_list_t *list) {
    size_t i;
    for (i = 0; i < list->count; ++i) {
        free(list->owned_names[i]);
        free(list->owned_seqs[i]);
    }
    free(list->owned_names);
    free(list->owned_seqs);
    free(list->items);
}

static char *read_entire_file(const char *path) {
    FILE *fp = fopen(path, "rb");
    long sz;
    char *buf;
    if (!fp) return NULL;
    fseek(fp, 0, SEEK_END);
    sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(fp); return NULL; }
    if (fread(buf, 1, (size_t)sz, fp) != (size_t)sz) {
        free(buf);
        fclose(fp);
        return NULL;
    }
    buf[sz] = '\0';
    fclose(fp);
    return buf;
}

static int read_fasta(const char *path, fasta_list_t *out) {
    FILE *fp = fopen(path, "r");
    char line[8192];
    char *cur_seq = NULL;
    size_t cur_len = 0, cur_cap = 0;
    char *cur_name = NULL;

    memset(out, 0, sizeof(*out));
    if (!fp) return -1;

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;
        if (line[0] == '>') {
            if (cur_name) {
                if (out->count == out->cap) {
                    out->cap = (out->cap == 0) ? 64 : out->cap * 2;
                    out->items = (tw_seq_t *)realloc(out->items, out->cap * sizeof(tw_seq_t));
                    out->owned_names = (char **)realloc(out->owned_names, out->cap * sizeof(char *));
                    out->owned_seqs = (char **)realloc(out->owned_seqs, out->cap * sizeof(char *));
                }
                out->owned_names[out->count] = cur_name;
                out->owned_seqs[out->count] = cur_seq ? cur_seq : (char *)calloc(1, 1);
                out->items[out->count].name = cur_name;
                out->items[out->count].name_len = strlen(cur_name);
                out->items[out->count].seq = out->owned_seqs[out->count];
                out->items[out->count].seq_len = cur_len;
                out->count++;
                cur_seq = NULL;
                cur_len = cur_cap = 0;
            }
            cur_name = (char *)malloc(len);
            memcpy(cur_name, line + 1, len);
        } else {
            if (cur_len + len + 1 > cur_cap) {
                cur_cap = (cur_cap == 0) ? (len + 256) : (cur_cap * 2 + len);
                cur_seq = (char *)realloc(cur_seq, cur_cap);
            }
            memcpy(cur_seq + cur_len, line, len);
            cur_len += len;
            cur_seq[cur_len] = '\0';
        }
    }
    if (cur_name) {
        if (out->count == out->cap) {
            out->cap = (out->cap == 0) ? 64 : out->cap * 2;
            out->items = (tw_seq_t *)realloc(out->items, out->cap * sizeof(tw_seq_t));
            out->owned_names = (char **)realloc(out->owned_names, out->cap * sizeof(char *));
            out->owned_seqs = (char **)realloc(out->owned_seqs, out->cap * sizeof(char *));
        }
        out->owned_names[out->count] = cur_name;
        out->owned_seqs[out->count] = cur_seq ? cur_seq : (char *)calloc(1, 1);
        out->items[out->count].name = cur_name;
        out->items[out->count].name_len = strlen(cur_name);
        out->items[out->count].seq = out->owned_seqs[out->count];
        out->items[out->count].seq_len = cur_len;
        out->count++;
    }
    fclose(fp);
    return 0;
}

static int run_inline_demo(void) {
    tw_opt_t opt;
    tw_msa_t *msa = NULL;
    tw_status_t st;
    size_t i;

    const char *tree = "((seq1:0.1,seq2:0.2):0.1,seq3:0.3);";
    tw_seq_t inputs[3] = {
        { "seq1", 0, "ACGTACGTACGT", 12 },
        { "seq2", 0, "ACGTTCGTACGT", 12 },
        { "seq3", 0, "ACGTACGACGT",  11 }
    };

    printf("TWILIGHT C API version: %s\n", tw_version());

    tw_opt_init(&opt);
    opt.accuracy_level = TW_ACC_PAIRWISE;
    opt.num_threads = 2;
    opt.quiet = 1;

    st = tw_align(&opt, tree, inputs, 3, &msa);
    if (st != TW_OK) {
        fprintf(stderr, "Inline demo failed: %s\n", tw_strerror(st));
        return 1;
    }

    printf("Inline demo aligned %zu sequences (aln_len = %zu):\n", msa->n_seqs, msa->aln_len);
    for (i = 0; i < msa->n_seqs; ++i) {
        printf(">%s (input_idx=%d)\n%s\n", msa->names[i], msa->input_indices[i], msa->seqs[i]);
    }
    tw_msa_free(msa);

    /* Test error handling without process exit */
    st = tw_align(&opt, "((broken_newick", inputs, 3, &msa);
    if (st != TW_ERR_TREE_PARSE || msa != NULL) {
        fprintf(stderr, "Expected TW_ERR_TREE_PARSE, got %d (%s)\n", (int)st, tw_strerror(st));
        return 1;
    }
    opt.remove_gappy = 2.0f;
    st = tw_align(&opt, tree, inputs, 3, &msa);
    if (st != TW_ERR_INVALID_ARG || msa != NULL) {
        fprintf(stderr, "Expected TW_ERR_INVALID_ARG, got %d (%s)\n", (int)st, tw_strerror(st));
        return 1;
    }
    printf("Error handling checks passed!\n");
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        return run_inline_demo();
    }

    if (argc < 5) {
        fprintf(stderr, "Usage: %s <tree.nwk> <seqs.fa> <out.aln> <accuracy_level:1|2|3> [threads]\n", argv[0]);
        return 1;
    }

    const char *tree_path = argv[1];
    const char *seq_path = argv[2];
    const char *out_path = argv[3];
    int acc_level = atoi(argv[4]);
    int threads = (argc >= 6) ? atoi(argv[5]) : 4;

    char *tree_str = read_entire_file(tree_path);
    fasta_list_t seqs;
    tw_opt_t opt;
    tw_msa_t *msa = NULL;
    tw_status_t st;
    FILE *out_fp;
    size_t i;

    if (!tree_str) {
        fprintf(stderr, "Failed to read tree file: %s\n", tree_path);
        return 1;
    }
    if (read_fasta(seq_path, &seqs) != 0) {
        fprintf(stderr, "Failed to read FASTA file: %s\n", seq_path);
        free(tree_str);
        return 1;
    }

    tw_opt_init(&opt);
    opt.accuracy_level = acc_level;
    opt.num_threads = threads;
    opt.quiet = 1;

    st = tw_align(&opt, tree_str, seqs.items, seqs.count, &msa);
    if (st != TW_OK) {
        fprintf(stderr, "tw_align failed: %s\n", tw_strerror(st));
        fasta_list_free(&seqs);
        free(tree_str);
        return 1;
    }

    out_fp = fopen(out_path, "w");
    if (!out_fp) {
        fprintf(stderr, "Failed to open output file: %s\n", out_path);
        tw_msa_free(msa);
        fasta_list_free(&seqs);
        free(tree_str);
        return 1;
    }
    for (i = 0; i < msa->n_seqs; ++i) {
        fprintf(out_fp, ">%s\n%s\n", msa->names[i], msa->seqs[i]);
    }
    fclose(out_fp);

    tw_msa_free(msa);
    fasta_list_free(&seqs);
    free(tree_str);
    return 0;
}

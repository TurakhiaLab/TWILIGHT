#!/usr/bin/env python3
"""
================================================================================
TWILIGHT + DIPPER Iterative Alignment Refinement Pipeline
================================================================================

DISCLAIMER:
  This iterative refinement pipeline is specifically designed for NUCLEOTIDE
  sequences (DNA / RNA). DIPPER currently does NOT support protein sequences.

Description:
  Alternates progressive multiple sequence alignment (TWILIGHT accurate mode)
  with ultra-fast guide tree reconstruction from alignment profiles (DIPPER).

  In each iteration:
    1. TWILIGHT aligns sequences guided by the current tree (-A 2 or -A 3)
       and computes an objective score (Consistency score or WSP).
    2. The alignment score is evaluated against previous iterations.
    3. If convergence criteria or early stopping conditions are met (patience),
       the loop terminates and rolls back to the highest-scoring checkpoint.
    4. Otherwise, DIPPER reconstructs a refined guide tree from the current
       alignment, and the next iteration begins.

Usage Examples:
  # Basic run with auto guide tree generation:
  python3 scripts/twilight_dipper_iterative.py -i input.fasta -o output.aln

  # Custom iterations, patience, and CPU threads:
  python3 scripts/twilight_dipper_iterative.py -i input.fasta -o output.aln \\
      --max-iters 10 --min-iters 3 --patience 2 --threads 16

  # Using custom guide tree and transitive consistency (-A 3):
  python3 scripts/twilight_dipper_iterative.py -i input.fasta -o output.aln \\
      -t initial.tree -A 3 --score-type 1

Parameters:
  -i, --input            Input unaligned FASTA file (Required).
  -o, --output           Output aligned MSA file in FASTA format (Required).
  -t, --tree             Initial guide tree (Newick format, optional). If not
                         provided, DIPPER generates one from input sequences.
  --max-iters            Maximum number of iterations (Default: 10).
  --min-iters            Minimum iterations before early stopping (Default: 3).
  --patience             Patience: stop after N iterations without score
                         improvement (Default: 2).
  --min-delta            Minimum relative score improvement threshold to count
                         as progress (Default: 0.001, i.e. 0.1%).
  -A, --acc-mode         TWILIGHT accuracy mode: 2 (+ pairwise consistency, default)
                         or 3 (+ transitive consistency).
  --score-type           Scoring method for convergence:
                           1: Consistency Objective Score (Default, recommended)
                           2: WSP (Weighted Profile Sum-of-Pairs)
  --sample-max-group     Max sequences per consistency group in TWILIGHT (Default: 100).
  --refine-subtree       Enable Layer 1 guide subtree refinement (Default: True).
  --refine-method        Subtree refinement algorithm: 'upgma' (Default) or 'nj'.
  --twilight-bin         Path to TWILIGHT executable (Default: <repo>/bin/twilight).
  --dipper-bin           Path to DIPPER executable (Default: dipper_cpu).
  -c, --threads          Number of CPU threads to use (Default: 16).
  --keep-temp            Keep intermediate files for debugging (Default: False).
  --keep-tree            Export the best guide tree alongside MSA as <output>.tree.
================================================================================
"""

import os
import sys
import argparse
import subprocess
import shutil
import tempfile
import re
import time

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="TWILIGHT + DIPPER Iterative Alignment Refinement Pipeline (for Nucleotides)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Disclaimer:
  This iterative script is designed for NUCLEOTIDE sequences (DNA/RNA)
  because DIPPER currently does not support protein sequence guide trees.
"""
    )

    # I/O Arguments
    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument("-i", "--input", required=True, type=str,
                          help="Input unaligned sequences in FASTA format (Required).")
    io_group.add_argument("-o", "--output", required=True, type=str,
                          help="Output aligned MSA file in FASTA format (Required).")
    io_group.add_argument("-t", "--tree", default=None, type=str,
                          help="Initial guide tree in Newick format (Optional). "
                               "If omitted, DIPPER will automatically generate a tree from unaligned sequences.")

    # Iteration & Stopping Criteria
    iter_group = parser.add_argument_group("Iteration & Stopping Criteria")
    iter_group.add_argument("--max-iters", default=10, type=int,
                            help="Maximum number of refinement iterations (Default: 10).")
    iter_group.add_argument("--min-iters", default=3, type=int,
                            help="Minimum iterations to run before early stopping can trigger (Default: 3).")
    iter_group.add_argument("--patience", default=2, type=int,
                            help="Number of iterations without score improvement before early stopping (Default: 2).")
    iter_group.add_argument("--min-delta", default=0.001, type=float,
                            help="Minimum relative score improvement required to qualify as progress (Default: 0.001 = 0.1%%).")

    # TWILIGHT Alignment Options
    tw_group = parser.add_argument_group("TWILIGHT Alignment Options")
    tw_group.add_argument("-A", "--acc-mode", default=2, type=int, choices=[2, 3],
                          help="TWILIGHT accuracy level: 2 (+ pairwise consistency, default) or 3 (+ transitive consistency).")
    tw_group.add_argument("--score-type", default=1, type=int, choices=[1, 2],
                          help="Convergence scoring objective: 1: Consistency Score (Default, recommended), 2: WSP.")
    tw_group.add_argument("--sample-max-group", default=100, type=int,
                          help="Maximum group size for consistency partitioning in TWILIGHT (Default: 100).")
    tw_group.add_argument("--refine-subtree", dest="refine_subtree", action="store_true", default=True,
                          help="Enable Layer 1 guide subtree refinement using pairwise local alignment distances (Default: True).")
    tw_group.add_argument("--no-refine-subtree", dest="refine_subtree", action="store_false",
                          help="Disable Layer 1 guide subtree refinement.")
    tw_group.add_argument("--refine-method", default="upgma", choices=["upgma", "nj"],
                          help="Algorithm for Layer 1 subtree refinement: 'upgma' (Default, recommended) or 'nj'.")

    # Executables & Hardware
    exec_group = parser.add_argument_group("Executables & Hardware Options")
    exec_group.add_argument("--twilight-bin", default=None, type=str,
                            help="Path to the 'twilight' executable (Default: auto-detected as ../bin/twilight).")
    exec_group.add_argument("--dipper-bin", default=None, type=str,
                            help="Path to the 'dipper' or 'dipper_cpu' executable (Default: auto-detected).")
    exec_group.add_argument("-c", "--threads", default=16, type=int,
                            help="Number of CPU threads to use for both TWILIGHT and DIPPER (Default: 16).")
    exec_group.add_argument("--temp-dir", default=None, type=str,
                            help="Custom directory for temporary iteration files (Default: auto-generated in system temp).")
    exec_group.add_argument("--keep-temp", action="store_true", default=False,
                            help="Keep intermediate iteration files after completion (Default: False).")
    exec_group.add_argument("--keep-tree", action="store_true", default=True,
                            help="Save the best reconstructed guide tree as <output>.tree (Default: True).")
    exec_group.add_argument("-v", "--verbose", action="store_true", default=False,
                            help="Enable verbose output during execution.")

    return parser.parse_args()


def resolve_executables(twilight_arg, dipper_arg):
    """Auto-detect and validate TWILIGHT and DIPPER executable binaries."""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Resolve TWILIGHT
    tw_candidates = []
    if twilight_arg:
        tw_candidates.append(os.path.abspath(twilight_arg))
    tw_candidates.extend([
        os.path.abspath(os.path.join(script_dir, "..", "bin", "twilight")),
        os.path.abspath(os.path.join(script_dir, "bin", "twilight")),
        shutil.which("twilight")
    ])

    twilight_bin = None
    for cand in tw_candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            twilight_bin = cand
            break

    if not twilight_bin:
        sys.stderr.write("ERROR: Could not locate executable 'twilight' binary.\n"
                         "Please specify it via --twilight-bin /path/to/twilight\n")
        sys.exit(1)

    # 2. Resolve DIPPER
    dip_candidates = []
    if dipper_arg:
        dip_candidates.append(os.path.abspath(dipper_arg))
    dip_candidates.extend([
        "/home/y3tseng@AD.UCSD.EDU/DIPPER_main/DIPPER/bin/dipper_cpu",
        "/home/y3tseng@AD.UCSD.EDU/DIPPER_main/DIPPER/bin/dipper",
        shutil.which("dipper_cpu"),
        shutil.which("dipper")
    ])

    dipper_bin = None
    for cand in dip_candidates:
        if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
            dipper_bin = cand
            break

    if not dipper_bin:
        sys.stderr.write("ERROR: Could not locate executable DIPPER binary ('dipper_cpu' or 'dipper').\n"
                         "Please specify it via --dipper-bin /path/to/dipper_cpu\n")
        sys.exit(1)

    return twilight_bin, dipper_bin


def sanitize_newick_tree(in_tree_path, out_tree_path):
    """
    Sanitize Newick tree by replacing invalid root / branch lengths like
    ':-nan' or ':nan' (which can occur in DIPPER outputs) with ':0.0001'.
    """
    with open(in_tree_path, 'r') as f:
        content = f.read()
    sanitized = re.sub(r':[-+]?(?:nan|inf)', ':0.0001', content, flags=re.IGNORECASE)
    with open(out_tree_path, 'w') as f:
        f.write(sanitized)


def parse_score_file(score_file_path):
    """
    Parse TWILIGHT .score file to extract:
      - primary score (line 1)
      - consistency score
      - WSP score
      - alignment length
      - gap ratio
    """
    stats = {
        'score': 0.0,
        'consistency': 0.0,
        'wsp': 0.0,
        'aln_len': 0,
        'gap_ratio': 0.0
    }
    if not os.path.isfile(score_file_path):
        return stats

    with open(score_file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        return stats

    try:
        stats['score'] = float(lines[0])
    except ValueError:
        pass

    for line in lines[1:]:
        if line.startswith("CONSISTENCY_SCORE="):
            try: stats['consistency'] = float(line.split("=")[1])
            except ValueError: pass
        elif line.startswith("WSP_SCORE="):
            try: stats['wsp'] = float(line.split("=")[1])
            except ValueError: pass
        elif line.startswith("ALIGNMENT_LENGTH="):
            try: stats['aln_len'] = int(line.split("=")[1])
            except ValueError: pass
        elif line.startswith("GAP_RATIO="):
            try: stats['gap_ratio'] = float(line.split("=")[1])
            except ValueError: pass

    return stats


def print_banner():
    """Print pipeline banner and disclaimer."""
    print("=" * 80)
    print("  TWILIGHT + DIPPER Iterative Alignment Refinement Pipeline")
    print("=" * 80)
    print("  [DISCLAIMER] Designed specifically for NUCLEOTIDE sequences (DNA/RNA).")
    print("               DIPPER guide tree reconstruction currently does not support protein.")
    print("=" * 80)
    sys.stdout.flush()


def main():
    args = parse_arguments()
    print_banner()

    # Validate input FASTA
    input_fasta = os.path.abspath(args.input)
    if not os.path.isfile(input_fasta) or os.path.getsize(input_fasta) == 0:
        sys.stderr.write(f"ERROR: Input sequence file '{input_fasta}' does not exist or is empty.\n")
        sys.exit(1)

    output_aln = os.path.abspath(args.output)
    output_dir = os.path.dirname(output_aln)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Resolve executable binaries
    twilight_bin, dipper_bin = resolve_executables(args.twilight_bin, args.dipper_bin)
    print(f"[*] TWILIGHT Binary : {twilight_bin}")
    print(f"[*] DIPPER Binary   : {dipper_bin}")
    print(f"[*] Input Sequences : {input_fasta}")
    print(f"[*] Output MSA      : {output_aln}")
    print(f"[*] Max Iterations  : {args.max_iters} (Min: {args.min_iters}, Patience: {args.patience})")
    print(f"[*] Score Mode      : {'1: Consistency Objective Score' if args.score_type == 1 else '2: WSP'}")
    print(f"[*] Subtree Refine  : {args.refine_subtree} (Method: {args.refine_method}, Group Size: {args.sample_max_group})")
    print(f"[*] CPU Threads     : {args.threads}")
    print("-" * 80)
    sys.stdout.flush()

    # Create temporary working directory
    if args.temp_dir:
        temp_dir = os.path.abspath(args.temp_dir)
        os.makedirs(temp_dir, exist_ok=True)
        cleanup_temp = False
    else:
        temp_dir = tempfile.mkdtemp(prefix="twilight_iter_")
        cleanup_temp = not args.keep_temp

    pipeline_start_time = time.time()

    try:
        # Step 1: Initial Guide Tree
        if args.tree:
            init_tree = os.path.abspath(args.tree)
            if not os.path.isfile(init_tree):
                sys.stderr.write(f"ERROR: Provided initial tree '{init_tree}' not found.\n")
                sys.exit(1)
            clean_curr_tree = os.path.join(temp_dir, "iter0_clean.tree")
            sanitize_newick_tree(init_tree, clean_curr_tree)
            print(f"[*] Using provided guide tree: {init_tree}")
        else:
            print("[*] No guide tree provided. Generating initial guide tree using DIPPER from raw sequences...")
            raw_init_tree = os.path.join(temp_dir, "init_dipper_raw.tree")
            clean_curr_tree = os.path.join(temp_dir, "iter0_clean.tree")
            cmd_init_dip = [
                dipper_bin,
                "-i", "r",
                "-o", "t",
                "-m", "0",
                "-I", input_fasta,
                "-O", raw_init_tree,
                "-T", str(args.threads)
            ]
            t_dip0 = time.time()
            res_dip0 = subprocess.run(cmd_init_dip, capture_output=True, text=True)
            if res_dip0.returncode != 0:
                sys.stderr.write(f"ERROR: DIPPER failed to generate initial guide tree:\n{res_dip0.stderr}\n")
                sys.exit(1)
            sanitize_newick_tree(raw_init_tree, clean_curr_tree)
            print(f"[*] Initial guide tree built in {time.time() - t_dip0:.2f} s")

        print("=" * 105)
        print(f"{'Iter':<5} | {'Score':<11} | {'Consistency':<11} | {'WSP Score':<11} | {'Length':<7} | {'Gap %':<7} | {'TW(s)':<6} | {'DIP(s)':<6} | Status")
        print("-" * 105)
        sys.stdout.flush()

        best_score = -float('inf')
        best_iter = 1
        best_aln_path = os.path.join(temp_dir, "best_checkpoint.aln")
        best_tree_path = os.path.join(temp_dir, "best_checkpoint.tree")
        patience_counter = 0

        total_iters_run = 0

        # Iteration Loop
        for it in range(1, args.max_iters + 1):
            total_iters_run = it
            curr_aln = os.path.join(temp_dir, f"iter{it}.aln")
            curr_score_file = os.path.join(temp_dir, f"iter{it}.score")

            # 1. Run TWILIGHT
            cmd_tw = [
                twilight_bin,
                "-t", clean_curr_tree,
                "-i", input_fasta,
                "-o", curr_aln,
                "-A", str(args.acc_mode),
                "--sample-max-group", str(args.sample_max_group),
                "--cpu", str(args.threads),
                "--overwrite",
                "--score", str(args.score_type),
                "--score-file", curr_score_file
            ]
            if args.refine_subtree:
                cmd_tw.extend(["--refine-subtree", "--refine-method", args.refine_method])

            t_tw = time.time()
            res_tw = subprocess.run(cmd_tw, capture_output=True, text=True)
            tw_time = time.time() - t_tw
            if res_tw.returncode != 0 or not os.path.isfile(curr_aln):
                sys.stderr.write(f"\nERROR: TWILIGHT failed at iteration {it}:\n{res_tw.stderr}\n")
                sys.exit(1)

            # 2. Parse Score
            stats = parse_score_file(curr_score_file)
            curr_score = stats['score']
            gap_pct = stats['gap_ratio'] * 100.0

            # 3. Checkpoint Tracking & Early Stopping Evaluation
            improved = False
            if it == 1:
                improved = True
            else:
                # Require relative improvement > min_delta
                threshold = best_score + (abs(best_score) * args.min_delta if abs(best_score) > 1e-6 else args.min_delta)
                if curr_score > threshold:
                    improved = True

            status_str = ""
            if improved:
                best_score = curr_score
                best_iter = it
                patience_counter = 0
                shutil.copyfile(curr_aln, best_aln_path)
                shutil.copyfile(clean_curr_tree, best_tree_path)
                status_str = f"★ New Best (Iter {it})"
            else:
                patience_counter += 1
                status_str = f"Patience {patience_counter}/{args.patience}"

            # 4. DIPPER Tree Reconstruction (if not the last iteration)
            dip_time = 0.0
            if it < args.max_iters and not (it >= args.min_iters and patience_counter >= args.patience):
                raw_next_tree = os.path.join(temp_dir, f"iter{it}_dipper_raw.tree")
                clean_next_tree = os.path.join(temp_dir, f"iter{it}_clean.tree")
                cmd_dip = [
                    dipper_bin,
                    "-i", "m",
                    "-o", "t",
                    "-m", "1",
                    "-d", "1",
                    "-I", curr_aln,
                    "-O", raw_next_tree,
                    "-T", str(args.threads)
                ]
                t_dip = time.time()
                res_dip = subprocess.run(cmd_dip, capture_output=True, text=True)
                dip_time = time.time() - t_dip
                if res_dip.returncode != 0 or not os.path.isfile(raw_next_tree):
                    sys.stderr.write(f"\nWARNING: DIPPER failed to reconstruct tree from iteration {it} alignment:\n{res_dip.stderr}\n")
                    break
                sanitize_newick_tree(raw_next_tree, clean_next_tree)
                clean_curr_tree = clean_next_tree

            print(f"{it:<5} | {curr_score:<11.6f} | {stats['consistency']:<11.6f} | {stats['wsp']:<11.6f} | {stats['aln_len']:<7} | {gap_pct:<6.2f}% | {tw_time:<6.2f} | {dip_time:<6.2f} | {status_str}")
            sys.stdout.flush()

            # Early stopping check
            if it >= args.min_iters and patience_counter >= args.patience:
                print("-" * 105)
                print(f"[*] Early stopping triggered: Score did not improve for {patience_counter} consecutive iterations.")
                break

        print("=" * 105)

        # Final Rollback & Output Export
        if not os.path.isfile(best_aln_path):
            sys.stderr.write("ERROR: No valid alignment checkpoint was produced.\n")
            sys.exit(1)

        shutil.copyfile(best_aln_path, output_aln)
        print(f"[✓] Final MSA exported to: {output_aln}")
        print(f"[✓] Best iteration selected: Iteration {best_iter} (Score: {best_score:.6f})")

        if args.keep_tree and os.path.isfile(best_tree_path):
            best_out_tree = output_aln + ".tree"
            shutil.copyfile(best_tree_path, best_out_tree)
            print(f"[✓] Best guide tree exported to: {best_out_tree}")

        total_elapsed = time.time() - pipeline_start_time
        print(f"[*] Total Execution Time: {total_elapsed:.2f} s across {total_iters_run} iterations.")

    finally:
        if cleanup_temp and os.path.isdir(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()


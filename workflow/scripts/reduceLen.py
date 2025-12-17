import argparse
import time
import sys
import gzip
import numpy as np
from numba import jit, prange, set_num_threads

# --- Helper: Efficient FASTA Reader (Supports .gz) ---
def read_fasta_to_numpy(filename):
    """
    Reads a FASTA (or .gz FASTA) file and converts it to a 2D NumPy character array.
    """
    headers = []
    sequences = []
    current_seq = []
    
    # Determine if we need to open with gzip or standard open
    if filename.endswith('.gz'):
        # 'rt' mode opens it as text, handling the decompression automatically
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')

    try:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
                headers.append(line)
            else:
                current_seq.append(line)
        # Add the last sequence
        if current_seq:
            sequences.append("".join(current_seq))
    finally:
        f.close()
            
    if not sequences:
        raise ValueError("No sequences found in input file.")

    # Convert to NumPy array of characters (S1 = 1-byte string)
    try:
        # Use 'S1' (bytes) for performance equivalent to C++ char
        seq_matrix = np.array([list(s) for s in sequences], dtype='S1')
    except ValueError:
        raise ValueError("Sequences must all be the same length for this operation.")
        
    return headers, seq_matrix

# --- Core Logic: Parallel Filter (Numba) ---
# nopython=True: Compile to machine code (no Python interpreter slow-down)
# parallel=True: Enable automatic parallelization (OpenMP/TBB backend)
@jit(nopython=True, parallel=True)
def get_column_mask(seq_matrix, threshold):
    n_seqs, n_cols = seq_matrix.shape
    max_allowed_gaps = int(np.floor(threshold * n_seqs))
    
    keep_mask = np.zeros(n_cols, dtype=np.bool_)
    
    for col_idx in prange(n_cols):
        gap_count = 0
        for seq_idx in range(n_seqs):
            # b'-' is the byte representation of a dash
            if seq_matrix[seq_idx, col_idx] == b'-':
                gap_count += 1
        
        if gap_count <= max_allowed_gaps:
            keep_mask[col_idx] = True
            
    return keep_mask

# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description='Reduce alignment length to speedup tree inference process')
    parser.add_argument('inaln', help='Input alignment (FASTA or .gz)')
    parser.add_argument('outaln', help='Output alignment (Uncompressed FASTA)')
    parser.add_argument('threshold', type=float, help='Minimum gap proportion for a column be removed')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    args = parser.parse_args()

    # 1. Configure Threads
    if args.threads > 1:
        set_num_threads(args.threads)
        print(f"Using {args.threads} threads.")

    try:
        print(f"Reading alignment from {args.inaln}...")
        headers, seq_matrix = read_fasta_to_numpy(args.inaln)
        
        num_seqs, num_cols = seq_matrix.shape
        print(f"Original dimensions: {num_seqs} sequences, {num_cols} columns")

        # Start Timing
        start_time = time.perf_counter()

        # 2. Parallel Analysis
        keep_mask = get_column_mask(seq_matrix, args.threshold)

        # 3. Filtering (Slicing)
        filtered_matrix = seq_matrix[:, keep_mask]

        end_time = time.perf_counter()
        elapsed_ms = (end_time - start_time) * 1000

        new_cols = filtered_matrix.shape[1]
        print(f"Original length: {num_cols}, length after removing gappy columns: {new_cols}")
        print(f"Remove gappy columns in {elapsed_ms:.2f} ms")

        # 4. Write Output (Uncompressed)
        print(f"Writing output to {args.outaln}...")
        with open(args.outaln, 'w') as f:
            for i, header in enumerate(headers):
                seq_str = filtered_matrix[i].tobytes().decode('utf-8')
                f.write(f"{header}\n{seq_str}\n")
        
        print("Done.")

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
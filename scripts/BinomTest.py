#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Optional ISA-L acceleration for gzip I/O (SIMD decompress + multithreaded
# compress). Falls back to stdlib gzip if not installed, so the script stays
# portable and output is standard gzip either way.
try:
    from isal import igzip as _igzip
    from isal import igzip_threaded as _igzip_threaded
    _HAVE_ISAL = True
except ImportError:
    _HAVE_ISAL = False


def _open_read(path):
    """Open a possibly-gzipped input for reading bytes, via ISA-L when available."""
    path = str(path)
    if path.endswith(".gz"):
        if _HAVE_ISAL:
            return _igzip.open(path, "rb")
        import gzip
        return gzip.open(path, "rb")
    return open(path, "rb")


def _write_tsv_gz(df, path, threads, compresslevel=3):
    """Write a standard-gzip TSV. Uses multithreaded ISA-L compression when
    available, else pandas' single-threaded gzip. Decompressed content is
    identical either way (compression level changes only the .gz byte size)."""
    if _HAVE_ISAL:
        level = max(0, min(3, compresslevel))  # ISA-L supports levels 0-3 only
        with _igzip_threaded.open(path, "wt", compresslevel=level,
                                  threads=max(1, threads), encoding="utf-8",
                                  newline="") as fh:
            df.to_csv(fh, sep="\t", index=False, lineterminator="\n")
    else:
        df.to_csv(path, sep="\t", index=False, compression="gzip")


def binom_pval(kn, null_prob):
    """Exact two-sided binomial test p-value for one (k, n) pair.

    Identical call as the original per-row version; only the number of calls
    changes (once per *unique* (meth, coverage) pair instead of per site).
    """
    k, n = kn
    return binomtest(int(k), n=int(n), p=null_prob, alternative='two-sided').pvalue


def main():
    parser = argparse.ArgumentParser(description="Binomial test for methylation sites with FDR correction")
    parser.add_argument("-i", "--input", required=True, help="Input summarized file (tsv.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output filtered file (tsv.gz)")
    parser.add_argument("--nonconv_chr", type=str, default=None, help="Non-conversion control chromosome (recommended when available, e.g., 'chloroplast')")
    parser.add_argument("--null_prob", type=float, default=None, help="Null hypothesis probability (use precomputed value when nonconv_chr is unavailable; default: estimated from nonconv_chr or 0.5)")
    parser.add_argument("--fdr_threshold", type=float, default=0.05, help="FDR threshold (default: 0.05)")
    parser.add_argument("--min_coverage", type=int, default=0, help="Minimum coverage to retain site (default: 0)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    args = parser.parse_args()

    # Validate input
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    # Create output directory if needed
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading input: {args.input}")
    with _open_read(args.input) as _in_fh:
        df = pd.read_csv(_in_fh, sep='\t', compression=None, low_memory=False, dtype={'chr': str})

    if df.empty:
        raise ValueError("Input file is empty")

    print(f"Loaded {len(df):,} sites")

    # Ensure required columns
    if 'start' not in df.columns and 'pos' in df.columns:
        df['start'] = df['pos']
    df['end'] = df['start'] + 1
    df['coverage'] = df['meth'] + df['unmeth']
    # Filter by coverage (single pass)
    min_cov = max(1, args.min_coverage)  # Ensure at least 1 to avoid division by zero
    df = df[df['coverage'] >= min_cov]
    print(f"After coverage filter (>={min_cov}): {len(df):,} sites")

    if df.empty:
        raise ValueError(f"No sites remaining after coverage filter (min={min_cov})")

    # Estimate or use null probability
    null_prob = args.null_prob
    if null_prob is None and args.nonconv_chr:
        subset = df[df['chr'] == args.nonconv_chr]
        if not subset.empty:
            null_prob = subset['meth'].sum() / (subset['meth'].sum() + subset['unmeth'].sum())
            print(f"Estimated null_prob from {args.nonconv_chr}: {null_prob:.4f}")
        else:
            print(f"Warning: No data found for nonconv_chr '{args.nonconv_chr}', using default 0.5")
            null_prob = 0.5
    elif null_prob is None:
        null_prob = 0.5
        print(f"Using default null_prob: {null_prob}")
    else:
        print(f"Using specified null_prob: {null_prob}")

    # Binomial test.
    #
    # The p-value of a two-sided binomial test depends only on (k, n) = (meth,
    # coverage) for a fixed null_prob, not on the row. Across ~1.2e8 sites there
    # are only tens of thousands of distinct (k, n) pairs, so we compute the
    # exact scipy.stats.binomtest once per unique pair and map the result back
    # to every site. This is bit-for-bit identical to the original per-row
    # apply(), but calls binomtest 3-4 orders of magnitude fewer times.
    print(f"Running binomial tests (null_prob={null_prob:.4f}, threads={args.threads})...")
    uniq = df[['meth', 'coverage']].drop_duplicates().reset_index(drop=True)
    print(f"Unique (k, n) pairs: {len(uniq):,} (vs {len(df):,} sites)")

    pairs = list(zip(uniq['meth'].to_numpy(), uniq['coverage'].to_numpy()))
    if args.threads > 1 and len(pairs) > 1000:
        with mp.Pool(args.threads) as pool:
            pvals = list(tqdm(pool.imap(partial(binom_pval, null_prob=null_prob), pairs, chunksize=256),
                              total=len(pairs), desc="Binom tests (unique pairs)"))
    else:
        pvals = [binom_pval(p, null_prob) for p in tqdm(pairs, desc="Binom tests (unique pairs)")]
    uniq['pval'] = pvals

    # Map p-values back to all sites (vectorized left join on integer keys,
    # preserves the row order of the filtered frame).
    print("Mapping p-values back to all sites...")
    df = df.merge(uniq, on=['meth', 'coverage'], how='left')

    # FDR correction
    print("Applying FDR correction...")
    df['FDR'] = multipletests(df['pval'], method='fdr_bh')[1]

    n_sig = (df['FDR'] <= args.fdr_threshold).sum()
    print(f"Significant sites (FDR <= {args.fdr_threshold}): {n_sig:,} / {len(df):,} ({100*n_sig/len(df):.2f}%)")

    # Set meth=0 for non-significant sites (Weighted methylation level approach)
    df.loc[df['FDR'] > args.fdr_threshold, 'meth'] = 0
    print(f"Non-significant sites set to meth=0 for Weighted methylation level calculation")

    # Output
    df = df[['chr', 'pos', 'strand', 'meth', 'unmeth', 'context']]
    _write_tsv_gz(df, args.output, threads=args.threads)
    print(f"\nOutput written to: {args.output}")

if __name__ == "__main__":
    main()

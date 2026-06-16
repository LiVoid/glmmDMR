#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


# modkit pileup BED (0-based BED-like) expected columns:
# 1: chr, 2: start, 3: end, 4: name (e.g., m,CHH,0), 6: strand,
# 12: n_mod, 13: n_canonical
COL_INDEX = {
    "chr": 0,
    "start": 1,
    "end": 2,
    "name": 3,
    "strand": 5,
    "n_mod": 11,
    "n_canonical": 12,
}


def normalize_context(ctx: str) -> str:
    if ctx == "CG":
        return "CpG"
    return ctx


def parse_context(name_series: pd.Series) -> pd.Series:
    # name field is typically "m,CHH,0". Use second token as context.
    tokens = name_series.fillna("").astype(str).str.split(",")
    context = tokens.str[1].fillna("NA")
    return context.map(normalize_context)


def load_modkit_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        compression="infer",
        usecols=[
            COL_INDEX["chr"],
            COL_INDEX["start"],
            COL_INDEX["name"],
            COL_INDEX["strand"],
            COL_INDEX["n_mod"],
            COL_INDEX["n_canonical"],
        ],
        names=["chr", "start", "name", "strand", "n_mod", "n_canonical"],
        dtype={"chr": str, "start": "int64", "name": str, "strand": str},
    )

    if df.empty:
        return pd.DataFrame(columns=["chr", "pos", "strand", "meth", "unmeth", "context"])

    df["context"] = parse_context(df["name"])
    df["pos"] = df["start"] + 1  # convert BED start (0-based) to 1-based cytosine position
    df["meth"] = pd.to_numeric(df["n_mod"], errors="coerce").fillna(0).astype("int64")
    df["unmeth"] = pd.to_numeric(df["n_canonical"], errors="coerce").fillna(0).astype("int64")

    out = df[["chr", "pos", "strand", "meth", "unmeth", "context"]].copy()
    out = out.groupby(["chr", "pos", "strand", "context"], as_index=False)[["meth", "unmeth"]].sum()
    return out[["chr", "pos", "strand", "meth", "unmeth", "context"]]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert modkit pileup BED to BinomTest.py input format"
    )
    parser.add_argument("-i", "--input", required=True, help="Input modkit pileup BED(.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output tsv[.gz]")
    parser.add_argument(
        "--min_coverage",
        type=int,
        default=0,
        help="Keep rows with meth+unmeth >= this value (default: 0)",
    )
    parser.add_argument(
        "--contexts",
        nargs="*",
        default=None,
        help="Optional context filter (e.g., CpG CHG CHH)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    out = load_modkit_bed(input_path)

    if args.min_coverage > 0:
        cov = out["meth"] + out["unmeth"]
        out = out[cov >= args.min_coverage]

    if args.contexts:
        keep_contexts = {normalize_context(x) for x in args.contexts}
        out = out[out["context"].isin(keep_contexts)]

    compression = "gzip" if str(output_path).endswith(".gz") else None
    out.to_csv(output_path, sep="\t", index=False, compression=compression)

    print(f"Input rows: {len(out):,}")
    print(f"Output written: {output_path}")
    if not out.empty:
        print("Context counts:")
        print(out["context"].value_counts().sort_index())


if __name__ == "__main__":
    main()

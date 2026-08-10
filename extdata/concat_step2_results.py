#!/usr/bin/env python3
"""Combine CASTIE Step 2 files into the wide table consumed by Step 3."""

import argparse
import os
import re
import sys
from pathlib import Path

os.environ.setdefault("POLARS_MAX_THREADS", "1")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate Step 2 results, split dynamic p-values, and filter by MAF."
    )
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--contexts",
        required=True,
        help="Comma-separated dynamic-covariate names in the order used by Step 2.",
    )
    parser.add_argument("--file-pattern", default="*", help="Input glob [default: *]")
    parser.add_argument(
        "--gene-regex",
        required=True,
        help="Regex applied to each filename; it must contain a named group '(?P<gene>...)'.",
    )
    parser.add_argument("--maf-min", type=float, default=0.05)
    parser.add_argument("--maf-max", type=float, default=0.95)
    parser.add_argument(
        "--threads", type=int, default=1, help="Polars worker threads [default: 1]"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl

    contexts = [name.strip() for name in args.contexts.split(",") if name.strip()]
    if not contexts:
        raise SystemExit("--contexts must contain at least one name")
    if not 0 <= args.maf_min < args.maf_max <= 1:
        raise SystemExit("Require 0 <= --maf-min < --maf-max <= 1")

    gene_regex = re.compile(args.gene_regex)
    if "gene" not in gene_regex.groupindex:
        raise SystemExit("--gene-regex must contain a named '(?P<gene>...)' group")

    files = sorted(path for path in args.input_dir.glob(args.file_pattern) if path.is_file())
    if not files:
        raise SystemExit(f"No files matched {args.input_dir / args.file_pattern}")

    frames = []
    for path in files:
        match = gene_regex.search(path.name)
        if match is None:
            print(f"[SKIP] filename does not match --gene-regex: {path.name}", file=sys.stderr)
            continue
        gene = match.group("gene")
        frame = (
            pl.scan_csv(
                path,
                separator="\t",
                null_values=["NA", "NaN", ""],
                schema_overrides={
                    "MarkerID": pl.String,
                    "AF_Allele2": pl.Float64,
                    "p.value": pl.Float64,
                    "pval_ge": pl.String,
                    "pval_ge_CCT": pl.Float64,
                },
            )
            .select(["MarkerID", "AF_Allele2", "p.value", "pval_ge", "pval_ge_CCT"])
            .filter(
                (pl.col("AF_Allele2") > args.maf_min)
                & (pl.col("AF_Allele2") < args.maf_max)
            )
            .with_columns(
                pl.lit(gene).alias("Gene"),
                pl.col("p.value").alias("pval_main"),
                pl.col("pval_ge").str.split(",").alias("dynamic_pvalues"),
            )
            .with_columns(
                [
                    pl.col("dynamic_pvalues")
                    .list.get(index, null_on_oob=True)
                    .cast(pl.Float64, strict=False)
                    .alias(context)
                    for index, context in enumerate(contexts)
                ]
            )
            .select(
                ["Gene", "MarkerID", "AF_Allele2", "pval_main"]
                + contexts
                + ["pval_ge_CCT"]
            )
        )
        frames.append(frame)

    if not frames:
        raise SystemExit("No filenames matched --gene-regex")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined = pl.concat(frames, how="vertical").collect(engine="streaming")
    combined.write_csv(args.output, separator="\t")
    print(f"Loaded files: {len(frames)}")
    print(f"Output rows: {combined.height}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()

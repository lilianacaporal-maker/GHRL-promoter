
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Epigenetic profiling of the putative GHRL promoter using ENCODE BigWig tracks.
- Loads DNase-seq, H3K4me3, H3K27ac and optional CTCF signals.
- Extracts values in a region and converts genomic coordinates to positions relative to a given TSS.
- Optionally saves merged CSV and a plot (disabled with --dry-run).

Usage (dry-run):
  python src/python/03_epigenetics_encode_tracks.py \
    --dnase   /path/to/ENCODE_dnase.bigWig \
    --h3k4me3 /path/to/ENCODE_H3K4me3.bigWig \
    --h3k27ac /path/to/ENCODE_H3K27ac.bigWig \
    --ctcf    /path/to/ENCODE_CTCF.bigWig \
    --chrom chr3 --start 10288734 --end 10292734 \
    --tss 10290734 \
    --dry-run
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pybigwig

def read_bw(bw_path, chrom, start, end, label):
    """Read BigWig values into a DataFrame with columns: pos, <label>."""
    if bw_path is None or not os.path.exists(bw_path):
        return None
    bw = pybigwig.open(bw_path)
    vals = bw.values(chrom, start, end)
    bw.close()
    vals = np.array([0.0 if v is None else float(v) for v in vals], dtype=float)
    return pd.DataFrame({"pos": np.arange(start, end), label: vals})

def add_rel(df, tss):
    """Add column 'rel' = position relative to TSS."""
    if df is None: return None
    df["rel"] = df["pos"] - tss
    return df

def plot_epigenetic(df_list, labels, out_png):
    """Plot multiple epigenetic signals over relative position."""
    plt.figure(figsize=(12, 6))
    colors = ["#2ca02c", "#1f77b4", "#ff7f0e", "#9467bd"]
    for df, lab, color in zip(df_list, labels, colors):
        if df is None: continue
        plt.plot(df["rel"], df[lab], label=lab, color=color)
    plt.axvline(0, color="black", linestyle="--", alpha=0.6, label="TSS (0)")
    plt.title("Epigenetic signals around the GHRL TSS")
    plt.xlabel("Position relative to TSS (bp)")
    plt.ylabel("Signal intensity (BigWig)")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dnase",   required=True, help="ENCODE BigWig DNase-seq.")
    ap.add_argument("--h3k4me3", required=True, help="ENCODE BigWig H3K4me3.")
    ap.add_argument("--h3k27ac", required=True, help="ENCODE BigWig H3K27ac.")
    ap.add_argument("--ctcf",    default=None,  help="ENCODE BigWig CTCF (optional).")
    ap.add_argument("--chrom",   default="chr3")
    ap.add_argument("--start",   type=int, default=10288734)
    ap.add_argument("--end",     type=int, default=10292734)
    ap.add_argument("--tss",     type=int, default=10290734)
    ap.add_argument("--out",     default="results/epigenetics_signals.csv")
    ap.add_argument("--fig",     default="results/epigenetics_region.png")
    ap.add_argument("--dry-run", action="store_true", help="Run without writing CSV/PNG.")
    args = ap.parse_args()

    dnase = read_bw(args.dnase,   args.chrom, args.start, args.end, "DNase")
    k4me3 = read_bw(args.h3k4me3, args.chrom, args.start, args.end, "H3K4me3")
    k27ac = read_bw(args.h3k27ac, args.chrom, args.start, args.end, "H3K27ac")
    ctcf  = read_bw(args.ctcf,    args.chrom, args.start, args.end, "CTCF") if args.ctcf else None

    dnase = add_rel(dnase, args.tss)
    k4me3 = add_rel(k4me3, args.tss)
    k27ac = add_rel(k27ac, args.tss)
    ctcf  = add_rel(ctcf,  args.tss) if ctcf is not None else None

    # Merge on 'pos'
    dfs = [d for d in [dnase, k4me3, k27ac, ctcf] if d is not None]
    merged = dfs[0][["pos", "rel"]].copy()
    for d in dfs:
        merged = merged.merge(d.drop(columns=["rel"]), on="pos", how="left")

    if args.dry_run:
        print("Dry-run: no files were written.")
    else:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        merged.to_csv(args.out, index=False)
        plot_epigenetic([dnase, k4me3, k27ac, ctcf], ["DNase", "H3K4me3", "H3K27ac", "CTCF"], args.fig)
        print(f"CSV -> {args.out}")
        print(f"Figure -> {args.fig}")

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Evolutionary conservation in the putative GHRL promoter region.
- Reads UCSC BigWig tracks: phastCons100way and phyloP100way.
- Extracts values in a given window (default: chr3:10288734-10292734).
- Detects conserved peaks using threshold >= 0.8.
- Computes Pearson and Spearman correlations.
- Optionally plots the profiles and saves CSV (disabled with --dry-run).

Usage (dry-run, no files written):
  python src/python/02_conservation_phastcons_phylop.py \
    --phastcons /path/to/hg38.phastCons100way.bw \
    --phylop    /path/to/hg38.phyloP100way.bw \
    --chrom chr3 --start 10288734 --end 10292734 \
    --threshold 0.8 \
    --dry-run
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pybigwig

def read_bigwig_values(bw_path, chrom, start, end):
    """Return a numpy array of BigWig values for [start, end) on a chromosome."""
    if not os.path.exists(bw_path):
        raise FileNotFoundError(f"BigWig not found: {bw_path}")
    bw = pybigwig.open(bw_path)
    vals = bw.values(chrom, start, end)
    bw.close()
    # Replace None with 0.0 for robust processing
    return np.array([0.0 if v is None else float(v) for v in vals], dtype=float)

def find_peak_intervals(values, threshold):
    """
    Find contiguous intervals where values >= threshold.
    Returns list of tuples: (start_idx, end_idx, mean_value, max_value).
    """
    peaks = []
    in_peak = False
    s = 0
    for i, v in enumerate(values):
        if v >= threshold and not in_peak:
            in_peak = True; s = i
        elif v < threshold and in_peak:
            e = i
            seg = values[s:e]
            peaks.append((s, e, float(np.mean(seg)), float(np.max(seg))))
            in_peak = False
    if in_peak:
        e = len(values)
        seg = values[s:e]
        peaks.append((s, e, float(np.mean(seg)), float(np.max(seg))))
    return peaks

def intersect_peaks(phast, phyl, threshold):
    """
    Compute intervals where BOTH phastCons and phyloP are >= threshold simultaneously.
    Returns list of (start_idx, end_idx).
    """
    mask = (phast >= threshold) & (phyl >= threshold)
    intervals = []
    in_run = False
    s = 0
    for i, flag in enumerate(mask):
        if flag and not in_run:
            in_run = True; s = i
        elif not flag and in_run:
            intervals.append((s, i))
            in_run = False
    if in_run:
        intervals.append((s, len(mask)))
    return intervals

def correlations(a, b):
    """Pearson and Spearman correlations, ignoring NaNs/Infs."""
    df = pd.DataFrame({"a": a, "b": b}).replace([np.inf, -np.inf], np.nan).dropna()
    return float(df["a"].corr(df["b"], method="pearson")), float(df["a"].corr(df["b"], method="spearman"))

def plot_profiles(phastcons, phylop, chrom, start, end, out_png, threshold):
    """Plot phastCons and phyloP signals over genomic coordinates."""
    x = np.arange(start, end)
    plt.figure(figsize=(12, 4))
    plt.plot(x, phastcons, label="phastCons100way", color="#1f77b4")
    plt.plot(x, phylop,    label="phyloP100way",   color="#ff7f0e", alpha=0.85)
    plt.axhline(threshold, color="gray", linestyle="--", label=f"Threshold {threshold}")
    plt.title(f"Evolutionary conservation {chrom}:{start}-{end}")
    plt.xlabel("Genomic coordinate (bp)")
    plt.ylabel("Conservation score")
    plt.legend()
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=200)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phastcons", required=True, help="UCSC BigWig phastCons100way (hg38/hg19).")
    ap.add_argument("--phylop",    required=True, help="UCSC BigWig phyloP100way (hg38/hg19).")
    ap.add_argument("--chrom",     default="chr3")
    ap.add_argument("--start",     type=int, default=10288734)
    ap.add_argument("--end",       type=int, default=10292734)
    ap.add_argument("--threshold", type=float, default=0.8)
    ap.add_argument("--out",       default="results/conservation_summary.csv")
    ap.add_argument("--fig",       default="results/conservation_region.png")
    ap.add_argument("--dry-run",   action="store_true", help="Run without writing CSV/PNG.")
    args = ap.parse_args()

    phast = read_bigwig_values(args.phastcons, args.chrom, args.start, args.end)
    phyl  = read_bigwig_values(args.phylop,    args.chrom, args.start, args.end)

    pearson, spearman = correlations(phast, phyl)
    peaks_phast = find_peak_intervals(phast, args.threshold)
    peaks_phyl  = find_peak_intervals(phyl,  args.threshold)
    both_intervals = intersect_peaks(phast, phyl, args.threshold)

    print(f"Pearson: {pearson:.3f} | Spearman: {spearman:.3f}")
    print(f"Peaks >= {args.threshold} (phastCons): {len(peaks_phast)} | (phyloP): {len(peaks_phyl)}")
    print(f"Intervals with BOTH metrics >= {args.threshold}: {len(both_intervals)}")

    # Assemble summary rows
    rows = []
    for source, peaks in [("phastCons", peaks_phast), ("phyloP", peaks_phyl)]:
        for (s_i, e_i, mean_v, max_v) in peaks:
            rows.append({
                "metric": source,
                "peak_start": args.start + s_i,
                "peak_end":   args.start + e_i,
                "length_bp":  e_i - s_i,
                "mean_score": mean_v,
                "max_score":  max_v
            })
    for (s_i, e_i) in both_intervals:
        rows.append({
            "metric": "both>=threshold",
            "peak_start": args.start + s_i,
            "peak_end":   args.start + e_i,
            "length_bp":  e_i - s_i,
            "mean_score": float(np.mean([phast[s_i:e_i], phyl[s_i:e_i]])),
            "max_score":  float(np.max([np.max(phast[s_i:e_i]), np.max(phyl[s_i:e_i])]))
        })

    df = pd.DataFrame(rows)

    if args.dry_run:
        print("Dry-run: no files were written.")
    else:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        df.to_csv(args.out, index=False)
        plot_profiles(phast, phyl, args.chrom, args.start, args.end, args.fig, args.threshold)
        print(f"CSV -> {args.out}")
        print(f"Figure -> {args.fig}")

if __name__ == "__main__":
    main()

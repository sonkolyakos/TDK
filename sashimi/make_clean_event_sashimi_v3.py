#!/usr/bin/env python3
"""
make_clean_event_sashimi_v3.py

Compact publication-style event-centric sashimi plots.

Changes:
  - Arcs (junction counts) are extracted directly from rMATS IJC/SJC to guarantee 
    visibility of the skipping junction (avoids the BAM "top 4" noise limit).
  - Maintained the original v2 title/subtitle formatting (no explicit PSI % breakdown).
  - Kept the depth filter (`--min-reads`) to drop phantom events.
  - Increased font sizes for numbers and labels for better readability (Title size kept original).
"""

from pathlib import Path
import argparse
import math

import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Arc

EVENT_TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]


def clean_gene(x):
    return str(x).replace('"', '').strip().upper()


def get_mean_sum(row, col1, col2):
    """Kiszámolja a replikátumok átlagos IJC és SJC olvasatszámának összegét."""
    def parse_counts(val_str):
        vals = [float(i) for i in str(val_str).split(',') if i.strip() and i.strip() != 'NA']
        return sum(vals) / len(vals) if vals else 0.0

    val1 = parse_counts(row.get(col1, '0'))
    val2 = parse_counts(row.get(col2, '0'))
    return val1 + val2


def read_best_event(rmats_dir, gene, preferred_type=None, min_reads_threshold=20.0):
    rows = []

    types = EVENT_TYPES
    if preferred_type:
        types = [preferred_type] + [x for x in EVENT_TYPES if x != preferred_type]

    for typ in types:
        f = Path(rmats_dir) / f"{typ}.MATS.JC.txt"
        if not f.exists():
            continue

        d = pd.read_csv(f, sep="\t", low_memory=False)

        gene_cols = [
            c for c in ["geneSymbol", "GeneID", "gene", "Gene_Symbol"]
            if c in d.columns
        ]

        if not gene_cols or "FDR" not in d.columns or "IncLevelDifference" not in d.columns:
            continue

        mask = np.zeros(len(d), dtype=bool)
        for c in gene_cols:
            mask |= d[c].map(clean_gene).eq(gene.upper()).to_numpy()

        x = d.loc[mask].copy()
        if x.empty:
            continue
            
        # --- Mélység kritérium (Depth filter) ---
        valid_mask = []
        for idx, row in x.iterrows():
            ctr_reads = get_mean_sum(row, 'IJC_SAMPLE_1', 'SJC_SAMPLE_1')
            tis_reads = get_mean_sum(row, 'IJC_SAMPLE_2', 'SJC_SAMPLE_2')
            valid_mask.append((ctr_reads >= min_reads_threshold) and (tis_reads >= min_reads_threshold))
            
        x = x[valid_mask].copy()
        
        if x.empty:
            continue
        # ----------------------------------------

        x["FDR_num"] = pd.to_numeric(x["FDR"], errors="coerce")
        x["dPSI_num"] = pd.to_numeric(x["IncLevelDifference"], errors="coerce")

        x = x[
            x["FDR_num"].le(0.05)
            & x["dPSI_num"].abs().ge(0.10)
        ].copy()

        if x.empty:
            continue

        x["event_type"] = typ
        rows.append(x)

        if preferred_type and typ == preferred_type:
            break

    if not rows:
        return None

    x = pd.concat(rows, ignore_index=True)
    x["abs_dpsi"] = x["dPSI_num"].abs()

    x = x.sort_values(
        ["FDR_num", "abs_dpsi"],
        ascending=[True, False]
    )

    return x.iloc[0]


def as_int(row, name):
    if name not in row.index:
        return None
    v = pd.to_numeric(pd.Series([row[name]]), errors="coerce").iloc[0]
    return None if pd.isna(v) else int(v)


def event_blocks(row):
    typ = row["event_type"]
    blocks = []

    def add(a, b, label):
        s = as_int(row, a)
        e = as_int(row, b)
        if s is not None and e is not None and e > s:
            blocks.append((s, e, label))

    if typ == "SE":
        add("upstreamES", "upstreamEE", "upstream")
        add("exonStart_0base", "exonEnd", "alternative")
        add("downstreamES", "downstreamEE", "downstream")
    elif typ in ("A5SS", "A3SS"):
        add("longExonStart_0base", "longExonEnd", "long")
        add("shortES", "shortEE", "short")
        add("flankingES", "flankingEE", "flanking")
    elif typ == "MXE":
        add("upstreamES", "upstreamEE", "upstream")
        add("1stExonStart_0base", "1stExonEnd", "alternative1")
        add("2ndExonStart_0base", "2ndExonEnd", "alternative2")
        add("downstreamES", "downstreamEE", "downstream")
    elif typ == "RI":
        add("upstreamES", "upstreamEE", "upstream")
        add("riExonStart_0base", "riExonEnd", "retained_region")
        add("downstreamES", "downstreamEE", "downstream")

    seen = set()
    out = []
    for b in sorted(blocks):
        key = (b[0], b[1])
        if key not in seen:
            out.append(b)
            seen.add(key)
    return out


def event_region(blocks, flank=150):
    if not blocks:
        raise ValueError("Could not define event region.")
    start = max(0, min(x[0] for x in blocks) - flank)
    end = max(x[1] for x in blocks) + flank
    return start, end


def bam_chrom(bam, chrom):
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        refs = set(fh.references)
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append("chr" + chrom)
    for c in candidates:
        if c in refs:
            return c
    raise ValueError(f"{chrom} not found in {bam}")


def read_cov(bam, chrom, start, end):
    cov = np.zeros(end - start, dtype=float)
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for read in fh.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            pos = read.reference_start
            for op, length in read.cigartuples or []:
                if op in (0, 7, 8):
                    a = max(start, pos)
                    b = min(end, pos + length)
                    if b > a:
                        cov[a-start:b-start] += 1
                    pos += length
                elif op in (2, 3):
                    pos += length
    return cov


def aggregate_coverage(bams, chrom0, start, end):
    coverages = []
    for bam in bams:
        chrom = bam_chrom(bam, chrom0)
        cov = read_cov(bam, chrom, start, end)
        coverages.append(cov)
    if not coverages:
        return np.zeros(end - start, dtype=float)
    return np.mean(np.vstack(coverages), axis=0)


def get_rmats_junctions(row, sample_id):
    """
    Közvetlenül az rMATS táblázatból nyeri ki a nyers IJC/SJC számokat.
    Ezzel elkerülhető a BAM olvasás miatti exklúziós ív elvesztése.
    """
    if row["event_type"] != "SE":
        return []

    ijc_str = str(row.get(f"IJC_SAMPLE_{sample_id}", "0"))
    sjc_str = str(row.get(f"SJC_SAMPLE_{sample_id}", "0"))

    ijc_vals = [float(x) for x in ijc_str.split(",") if x.strip()]
    sjc_vals = [float(x) for x in sjc_str.split(",") if x.strip()]

    ijc_mean = sum(ijc_vals) / len(ijc_vals) if ijc_vals else 0
    sjc_mean = sum(sjc_vals) / len(sjc_vals) if sjc_vals else 0

    up_ee = as_int(row, "upstreamEE")
    ex_st = as_int(row, "exonStart_0base")
    ex_en = as_int(row, "exonEnd")
    dn_es = as_int(row, "downstreamES")

    junctions = []
    if None not in (up_ee, ex_st, ex_en, dn_es):
        junctions.append((up_ee, ex_st, ijc_mean / 2.0))
        junctions.append((ex_en, dn_es, ijc_mean / 2.0))
        junctions.append((up_ee, dn_es, sjc_mean))

    return junctions


def build_transform(start, end, blocks, intron_display=90.0):
    exon_intervals = sorted((s, e) for s, e, _ in blocks)
    merged = []
    for s, e in exon_intervals:
        if not merged or s > merged[-1][1]:
            merged.append([s, e])
        else:
            merged[-1][1] = max(merged[-1][1], e)

    segments = []
    cursor = start
    display = 0.0

    for s, e in merged:
        if s > cursor:
            gap_display = min(float(s - cursor), intron_display)
            segments.append((cursor, s, display, display + gap_display))
            display += gap_display
        exon_width = max(1, e - s)
        segments.append((s, e, display, display + exon_width))
        display += exon_width
        cursor = e

    if cursor < end:
        gap_display = min(float(end - cursor), intron_display)
        segments.append((cursor, end, display, display + gap_display))
        display += gap_display

    def transform(vals):
        vals = np.asarray(vals, dtype=float)
        out = np.empty_like(vals)
        for idx, val in np.ndenumerate(vals):
            out[idx] = np.nan
            for gs, ge, ds, de in segments:
                if gs <= val <= ge:
                    frac = 0 if ge == gs else (val - gs) / (ge - gs)
                    out[idx] = ds + frac * (de - ds)
                    break
        return out
    return transform, display


def assign_label_levels(items, total_width):
    if not items:
        return []
    items = sorted(items, key=lambda x: x["center"])
    min_sep = total_width * 0.095
    levels = []
    for item in items:
        assigned = False
        for lev in range(4):
            conflict = False
            for old in levels:
                if (old["level"] == lev and abs(old["center"] - item["center"]) < min_sep):
                    conflict = True
                    break
            if not conflict:
                item["level"] = lev
                assigned = True
                break
        if not assigned:
            item["level"] = 3
        levels.append(item)
    return levels


def draw_track(ax, x_display, coverage, junctions, transform, label, ymax, total_width):
    ax.fill_between(x_display, 0, coverage, linewidth=0, alpha=0.95)
    ax.set_ylim(0, ymax * 1.48)
    # Megnövelt y-tengely felirat méret (9 -> 12)
    ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=12, fontweight="bold")
    ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    # Megnövelt y-tengely érték méret (7 -> 10)
    ax.tick_params(axis="y", labelsize=10, length=3, width=0.8)

    prepared = []
    for donor, acceptor, count in junctions:
        if count <= 0.0:
            continue
        xd, xa = transform([donor, acceptor])
        if not np.isfinite(xd) or not np.isfinite(xa) or xa <= xd:
            continue
        width = xa - xd
        span_fraction = width / max(total_width, 1)
        height = ymax * (0.20 + 0.24 * min(1.0, span_fraction))
        prepared.append({
            "xd": xd, "xa": xa, "width": width, 
            "center": (xd + xa) / 2, "height": height, "count": count
        })

    prepared = assign_label_levels(prepared, total_width)

    for item in prepared:
        arc = Arc((item["center"], 0), width=item["width"], height=2 * item["height"],
                  theta1=0, theta2=180, linewidth=0.65, edgecolor="0.38", alpha=0.72, zorder=3)
        ax.add_patch(arc)
        label_y = item["height"] + ymax * (0.035 + 0.075 * item["level"])
        # Megnövelt csomóponti szám méret (6.1 -> 10)
        ax.text(item["center"], label_y, f"{item['count']:.0f}", ha="center", va="bottom",
                fontsize=10, color="0.12", zorder=5, 
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.88, pad=0.5))


def draw_schematic(ax, blocks, transform, event_type, start, end):
    ax.set_ylim(-0.45, 0.65)
    ax.axis("off")
    xs, xe = transform([start, end])
    ax.plot([xs, xe], [0, 0], linewidth=0.65)

    for s, e, label in blocks:
        x1, x2 = transform([s, e])
        if not np.isfinite(x1) or not np.isfinite(x2):
            continue
        rect = Rectangle((x1, -0.15), max(1, x2 - x1), 0.30, linewidth=0, fill=True)
        ax.add_patch(rect)
        if label.startswith("alternative") or label in ("short", "retained_region"):
            # Megnövelt exon felirat méret (6.3 -> 9)
            ax.text((x1 + x2) / 2, 0.27, label.replace("_", " "), 
                    ha="center", va="bottom", fontsize=9)

    # Megnövelt eseménytípus felirat méret (7.4 -> 11)
    ax.text(xs, -0.31, f"{event_type} event", ha="left", va="top", fontsize=11, fontweight="bold")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellline", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--rmats", required=True)
    parser.add_argument("--ctr", nargs="+", required=True)
    parser.add_argument("--tis", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--preferred-type", default=None)
    parser.add_argument("--flank", type=int, default=150)
    
    # Mélység kritérium
    parser.add_argument("--min-reads", type=float, default=20.0, 
                        help="Minimum rMATS IJC+SJC average reads required.")
    
    # Dummy argumentumok a bash parancsok visszafelé kompatibilitása érdekében
    parser.add_argument("--min-junction", type=float, default=5, help="Ignored")
    parser.add_argument("--max-junctions", type=int, default=4, help="Ignored")

    args = parser.parse_args()

    row = read_best_event(args.rmats, args.gene, args.preferred_type, min_reads_threshold=args.min_reads)

    if row is None:
        print(f"[{args.cellline}] No strict rMATS event found with sufficient reads (>{args.min_reads}) for {args.gene}")
        raise SystemExit(1)

    blocks = event_blocks(row)
    start, end = event_region(blocks, args.flank)
    chrom = str(row["chr"]).replace('"', '').strip()

    ctr_cov = aggregate_coverage([Path(x) for x in args.ctr], chrom, start, end)
    tis_cov = aggregate_coverage([Path(x) for x in args.tis], chrom, start, end)

    transform, total_width = build_transform(start, end, blocks)
    genomic_x = np.arange(start, end)
    display_x = transform(genomic_x)

    ymax = max(float(np.nanpercentile(ctr_cov, 99.5)), float(np.nanpercentile(tis_cov, 99.5)), 1.0)

    # RMATS alapú ív kinyerés a biztos megjelenésért
    ctr_show = get_rmats_junctions(row, "1")
    tis_show = get_rmats_junctions(row, "2")

    fig = plt.figure(figsize=(6.3, 3.15))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.46], hspace=0.07)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)

    draw_track(ax1, display_x, ctr_cov, ctr_show, transform, "CTR", ymax, total_width)
    draw_track(ax2, display_x, tis_cov, tis_show, transform, "TIS", ymax, total_width)
    draw_schematic(ax3, blocks, transform, row["event_type"], start, end)
    ax1.set_xlim(0, total_width)

    # EREDETI V2 CÍM ÉS LÁBLÉC FORMÁTUM VISSZAÁLLÍTÁSA
    dpsi = -float(row["dPSI_num"])
    fdr = float(row["FDR_num"])

    # Eredeti cím méret (10.8)
    fig.suptitle(
        (
            f"{args.cellline}  {args.gene}  |  "
            f"{row['event_type']}   "
            f"ΔPSI={dpsi:+.3f}   "
            f"FDR={fdr:.2g}"
        ),
        fontsize=10.8,
        fontweight="bold",
        y=0.98
    )

    # Megnövelt lábléc méret (6.2 -> 9)
    fig.text(
        0.995,
        0.01,
        (
            f"{chrom}:{start:,}-{end:,}  |  "
            f"introns compressed"
        ),
        ha="right",
        va="bottom",
        fontsize=9
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(
        f"{args.cellline}\t"
        f"{args.gene}\t"
        f"{row['event_type']}\t"
        f"{chrom}:{start}-{end}\t"
        f"{dpsi:+.4f}\t"
        f"{fdr:.4g}\t"
        f"{out}"
    )

if __name__ == "__main__":
    main()

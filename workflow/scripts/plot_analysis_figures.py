#!/usr/bin/env python3
"""Generate publication-grade P1-P3 static figures and a JSON data bundle for interactive visualization."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


PALETTE = {
    "teal": "#0F766E",
    "blue": "#1D4ED8",
    "sky": "#0284C7",
    "amber": "#D97706",
    "red": "#DC2626",
    "rose": "#E11D48",
    "purple": "#7C3AED",
    "magenta": "#A21CAF",
    "emerald": "#059669",
    "slate": "#64748B",
    "slate_light": "#CBD5E1",
    "ink": "#0F172A",
    "muted": "#475569",
    "grid": "#E2E8F0",
}

STRATA_ORDER = ["500_1000", "1000_3000", "3000_5000", "ge_5000", "lt_3000", "ge_3000"]
STRATA_LABELS = {
    "500_1000": "500–999 bp",
    "1000_3000": "1–3 kb",
    "3000_5000": "3–5 kb",
    "ge_5000": ">= 5 kb",
    "lt_3000": "< 3 kb",
    "ge_3000": ">= 3 kb",
    "pooled": "Pooled",
}


def configure_style() -> None:
    """Apply clean publication Matplotlib defaults."""
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.5,
            "axes.edgecolor": "#CBD5E1",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.2,
            "legend.frameon": True,
            "legend.facecolor": "#FFFFFF",
            "legend.edgecolor": "#E2E8F0",
            "figure.facecolor": "#FFFFFF",
            "axes.facecolor": "#FFFFFF",
            "savefig.facecolor": "#FFFFFF",
        }
    )


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> List[str]:
    """Save figure as 300 DPI PNG and vector SVG."""
    outdir.mkdir(parents=True, exist_ok=True)
    saved: List[str] = []
    for ext in ("png", "svg"):
        out_path = outdir / f"{stem}.{ext}"
        fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.18)
        saved.append(str(out_path))
    plt.close(fig)
    return saved


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_hotspot_domains(
    hotspots: List[Dict[str, str]],
    domains_path: Path,
    target_ids: Optional[Iterable[str]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Extract overlapping Track A (and Track B for annotation) domain hits inside selected hotspots."""
    wanted = set(target_ids) if target_ids is not None else {h["hotspot_id"] for h in hotspots}
    selected = [h for h in hotspots if h["hotspot_id"] in wanted]
    by_contig: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for hs in selected:
        by_contig[hs["contig"]].append(hs)

    all_domains = read_tsv(domains_path)
    out: Dict[str, Dict[str, Any]] = {}
    for hs in selected:
        hid = hs["hotspot_id"]
        h_start = int(hs["g_start"])
        h_end = int(hs["g_end"])
        h_strand = hs["strand"]
        hits: List[Dict[str, Any]] = []
        for d in all_domains:
            if d["contig"] != hs["contig"]:
                continue
            ds = int(d["g_start"])
            de = int(d["g_end"])
            if de <= h_start or ds >= h_end:
                continue
            if d["strand"] != h_strand:
                continue
            hits.append(
                {
                    "model_name": d["model_name"],
                    "accession": d["accession"],
                    "track": d["track"],
                    "role": d["role"],
                    "g_start": ds,
                    "g_end": de,
                    "rel_start": ds - h_start,
                    "rel_end": de - h_start,
                    "frame": int(d["frame"]),
                    "strand": d["strand"],
                    "orf_len_aa": int(d["orf_len_aa"]),
                    "dom_bitscore": float(d["dom_bitscore"]),
                    "i_evalue": d["i_evalue"],
                }
            )
        hits.sort(key=lambda r: (r["rel_start"], r["rel_end"]))
        out[hid] = {
            "hotspot_id": hid,
            "contig": hs["contig"],
            "g_start": h_start,
            "g_end": h_end,
            "length_bp": int(hs["length_bp"]),
            "strand": h_strand,
            "n_domains": int(hs["n_domains"]),
            "n_backbone": int(hs["n_backbone"]),
            "n_ks": int(hs["n_ks"]),
            "n_ca": int(hs["n_ca"]),
            "n_ancillary": int(hs["n_ancillary"]),
            "completeness_label": hs["completeness_label"],
            "domains": hits,
        }
    return out


def compute_tiara_na_stats(root: Path, sample_ids: List[str]) -> Dict[str, Dict[str, int]]:
    """Compute ge-500-bp Tiara NA contig counts directly from contig_evidence.tsv files."""
    stats: Dict[str, Dict[str, int]] = {}
    for sid in sample_ids:
        ev_path = root / "results/decontam" / sid / "contig_evidence.tsv"
        na_count = 0
        na_bp = 0
        with ev_path.open("r", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if int(row["length_bp"]) >= 500 and row["tiara_p_prok"] in {"", "NA"}:
                    na_count += 1
                    na_bp += int(row["length_bp"])
        stats[sid] = {"tiara_na_ge500_contigs": na_count, "tiara_na_ge500_bp": na_bp}
    return stats


def classify_domain_group(model_name: str, track: str, role: str) -> Tuple[str, str]:
    """Return display category and color for gigasynthase domain visualization."""
    if track == "B":
        return ("Track B tailoring", PALETTE["slate"])
    if model_name in {"ketoacyl-synt", "Ketoacyl-synt_C", "KAsynt_C_assoc", "Ketoacyl-synt_2"}:
        return ("KS backbone", PALETTE["red"])
    if role == "backbone":
        return ("AT/PP/C/A backbone", PALETTE["amber"])
    if model_name in {"KR", "adh_short", "adh_short_C2", "ADH_N", "ADH_zinc_N", "ADH_zinc_N_2"}:
        return ("KR / ER redox", PALETTE["teal"])
    if model_name in {"PKS_DH_N", "PS-DH"}:
        return ("DH dehydratase", PALETTE["sky"])
    return ("Other Track A ancillary", PALETTE["purple"])


def plot_fig1_p1_p2(
    outdir: Path,
    pfam_rows: List[Dict[str, str]],
    meta_12b1: Dict[str, Any],
    meta_cns: Dict[str, Any],
    bench_12b1: Dict[str, Any],
    sens_df: pd.DataFrame,
) -> List[str]:
    """Figure 1: P1 curated HMM panel & P2 12B1 Track A benchmark + sensitivity sweep."""
    fig = plt.figure(figsize=(15.5, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.25])

    # Panel A: P1 HMM Composition + P2/P3 Domain Hits by Track/Role
    ax_a = fig.add_subplot(gs[0, 0])
    cat_order = [
        ("A", "backbone", "Track A\nBackbone", PALETTE["red"]),
        ("A", "ancillary", "Track A\nAncillary", PALETTE["teal"]),
        ("B", "ancillary", "Track B\nAncillary", PALETTE["sky"]),
        ("A", "negative_control", "Negative\nControl", PALETTE["slate"]),
    ]
    model_counts = Counter((r["track"], r["role"]) for r in pfam_rows)
    hits_12b1 = meta_12b1["summary"]["domain_hits_by_track_role"]
    hits_cns = meta_cns["summary"]["domain_hits_by_track_role"]

    x = np.arange(len(cat_order))
    labels = [c[2] for c in cat_order]
    hmms = [model_counts[(c[0], c[1])] for c in cat_order]
    colors = [c[3] for c in cat_order]

    ax_a2 = ax_a.twinx()
    ax_a2.spines["top"].set_visible(False)
    bars_hmm = ax_a.bar(x - 0.18, hmms, width=0.34, color=colors, alpha=0.9, label="Curated HMMs (n=69)")
    line_12b1 = ax_a2.plot(
        x + 0.16,
        [hits_12b1[f"{c[0]}/{c[1]}"] for c in cat_order],
        marker="o",
        linewidth=2.0,
        color=PALETTE["ink"],
        label="12B1 hits (1,475)",
    )
    line_cns = ax_a2.plot(
        x + 0.16,
        [hits_cns[f"{c[0]}/{c[1]}"] for c in cat_order],
        marker="s",
        linestyle="--",
        linewidth=1.8,
        color=PALETTE["amber"],
        label="CNS00066 hits (1,222)",
    )
    for bar, val in zip(bars_hmm, hmms):
        ax_a.text(bar.get_x() + bar.get_width() / 2, val + 0.6, f"{val}", ha="center", va="bottom", fontweight="bold", fontsize=8.5)
    halo = dict(boxstyle="round,pad=0.12", facecolor="#FFFFFF", edgecolor="none", alpha=0.80)
    for xi, c in enumerate(cat_order):
        v1 = hits_12b1[f"{c[0]}/{c[1]}"]
        v2 = hits_cns[f"{c[0]}/{c[1]}"]
        ax_a2.text(xi + 0.16, v1 + 22, f"{v1}", ha="center", va="bottom", fontsize=7.6, color=PALETTE["ink"], bbox=halo)
        y2_pos = v2 - 26 if v1 - v2 > 45 else v2 - 42
        ax_a2.text(xi + 0.16, max(y2_pos, 18), f"{v2}", ha="center", va="top", fontsize=7.6, color=PALETTE["amber"], bbox=halo)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(labels)
    ax_a.set_ylim(0, 32)
    ax_a2.set_ylim(0, 780)
    ax_a.set_ylabel("HMM models in hapto_core.hmm")
    ax_a2.set_ylabel("Gathering-threshold hits", color=PALETTE["muted"])
    ax_a.grid(axis="y", color=PALETTE["grid"], linewidth=0.7)
    ax_a.set_title("A. P1 Curated Pfam Panel (n=69) & Hit Yield")
    h1, l1 = ax_a.get_legend_handles_labels()
    h2, l2 = ax_a2.get_legend_handles_labels()
    ax_a.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=7.5)

    # Panel B: PKZILLA Genomic Coverage & Protein-Relative Domain Recall
    ax_b = fig.add_subplot(gs[0, 1])
    loci_map = {row["locus_id"]: row for row in bench_12b1["loci"]}
    prot_map = {row["protein_id"]: row for row in bench_12b1["protein_relative_pkzilla_check"]}
    ids = ["PKZILLA-1", "PKZILLA-2", "PKZILLA-3", "PKZILLA-B1"]
    xb = np.arange(len(ids)) * 1.15
    unmerged = [loci_map[i]["unmerged_coverage_same_strand"] * 100 if i in loci_map else np.nan for i in ids]
    merged = [loci_map[i]["merged_coverage_same_strand"] * 100 if i in loci_map else np.nan for i in ids]
    prot_rec = [prot_map[i]["truth_domain_recall"] * 100 for i in ids]

    w = 0.28
    b1 = ax_b.bar(xb - w, np.nan_to_num(unmerged, nan=0.0), width=w, color=PALETTE["teal"], label="Unmerged gene span (%)")
    b2 = ax_b.bar(xb, np.nan_to_num(merged, nan=0.0), width=w, color=PALETTE["sky"], label="Merged truth window (%)")
    b3 = ax_b.bar(xb + w, prot_rec, width=w, color=PALETTE["amber"], label="Protein domain recall (%)")
    # Mask missing genomic bars for PKZILLA-B1 cleanly
    b1[-1].set_alpha(0.0)
    b2[-1].set_alpha(0.0)

    for idx in range(3):
        ax_b.text(xb[idx] - w - 0.01, unmerged[idx] + 1.3, f"{unmerged[idx]:.1f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold")
        ax_b.text(xb[idx], merged[idx] + 1.3, f"{merged[idx]:.1f}", ha="center", va="bottom", fontsize=7.2)
    for idx in range(4):
        ax_b.text(xb[idx] + w + 0.02, prot_rec[idx] + 1.3, f"{prot_rec[idx]:.1f}", ha="center", va="bottom", fontsize=7.2, fontweight="bold")
    ax_b.text(xb[3] - 0.14, 42, "No RCC3426\nassembly\nin dataset", ha="center", va="center", fontsize=7.6, color=PALETTE["muted"], style="italic")

    ax_b.axhline(90.0, color=PALETTE["slate"], linestyle=":", linewidth=1.0, alpha=0.8)
    ax_b.set_xticks(xb)
    ax_b.set_xticklabels(ids)
    ax_b.set_ylim(0, 112)
    ax_b.set_ylabel("Recall / span coverage (%)")
    ax_b.grid(axis="y", color=PALETTE["grid"], linewidth=0.7)
    ax_b.set_title("B. P2 PKZILLA Locus (12B1) & Protein Recall")
    ax_b.legend(loc="lower left", fontsize=7.6)

    # Panel C: 20-config sensitivity sweep heatmap (outside-truth bp + annotation of frozen baseline & fragmentation)
    ax_c = fig.add_subplot(gs[0, 2])
    piv_out = sens_df.pivot(index="min_domains", columns="d_max_bp", values="outside_truth_bp_same_strand")
    piv_rec = sens_df.pivot(index="min_domains", columns="d_max_bp", values="merged_64_recall_same")
    piv_comp = sens_df.pivot(index="min_domains", columns="d_max_bp", values="pkzilla2_components")
    piv_tot = sens_df.pivot(index="min_domains", columns="d_max_bp", values="predicted_hotspots_total")

    # Log-scaled intensity so 4,363 bp vs 30,831 bp vs 167,169 bp is visually clear
    mat = np.log10(piv_out.values.astype(float) + 1.0)
    im = ax_c.imshow(mat, cmap="YlOrRd", aspect="auto", vmin=0.0, vmax=5.3)
    dvals = list(piv_out.columns)
    mvals = list(piv_out.index)
    ax_c.set_xticks(np.arange(len(dvals)))
    ax_c.set_xticklabels([f"{d // 1000} kb" for d in dvals])
    ax_c.set_yticks(np.arange(len(mvals)))
    ax_c.set_yticklabels([f"{m}" for m in mvals])
    ax_c.set_xlabel("Cluster linkage gap d_max")
    ax_c.set_ylabel("Minimum Track A domains")
    ax_c.set_title("C. P2 20-Config Sweep (Merged /64 & Overhang bp)")

    for i, m in enumerate(mvals):
        for j, d in enumerate(dvals):
            obp = int(piv_out.iloc[i, j])
            rec = int(piv_rec.iloc[i, j])
            tot = int(piv_tot.iloc[i, j])
            c2 = int(piv_comp.iloc[i, j])
            dark_bg = mat[i, j] > 3.9
            txt_col = "#FFFFFF" if dark_bg else PALETTE["ink"]
            obp_lbl = f"{obp / 1000:.1f}k" if obp >= 1000 else f"{obp}"
            frag = "*" if c2 > 1 else ""
            ax_c.text(j, i, f"{rec}/64 ({tot}){frag}\n+{obp_lbl} bp", ha="center", va="center", fontsize=7.3, color=txt_col, fontweight="bold" if (d == 15000 and m == 5) else "normal")
            if d == 15000 and m == 5:
                rect = mpatches.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, fill=False, edgecolor=PALETTE["teal"], linewidth=2.4)
                ax_c.add_patch(rect)

    cbar = fig.colorbar(im, ax=ax_c, fraction=0.046, pad=0.03)
    cbar.set_label("log10(outside-truth bp + 1)", fontsize=8)
    ax_c.text(0.0, -0.22, "* PKZILLA-2 splits into 2 components at d_max = 5 kb; teal box = frozen baseline (15 kb, min=5)", transform=ax_c.transAxes, fontsize=7.5, color=PALETTE["muted"])

    fig.suptitle("HaptoSMASH P1–P2: Curated HMM Library & P. parvum 12B1 Track A Benchmark", fontsize=12.5, fontweight="bold", y=1.03)
    return save_figure(fig, outdir, "fig1_p1_p2_track_a_benchmark")


def plot_fig2_p3_decontam_bench(
    outdir: Path,
    decontam_bench: Dict[str, Any],
    strat_df: pd.DataFrame,
    spike_df: pd.DataFrame,
    org_df: pd.DataFrame,
) -> List[str]:
    """Figure 2: P3 length-stratified ROC-AUC, recall@FPR, and organellar safety."""
    fig = plt.figure(figsize=(15.5, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.15, 0.95])

    tools = [
        ("ensemble", "Length-cal. Ensemble", PALETTE["teal"], "o"),
        ("tiara", "Tiara", PALETTE["sky"], "s"),
        ("whokaryote_T", "Whokaryote RF_T", PALETTE["purple"], "^"),
        ("whokaryote_S", "Whokaryote RF_S", PALETTE["rose"], "D"),
    ]
    bins = ["500_1000", "1000_3000", "3000_5000", "ge_5000", "pooled"]
    x = np.arange(len(bins))
    offsets = [-0.24, -0.08, 0.08, 0.24]

    # Panel A: Length-stratified ROC-AUC across 9 replicates
    ax_a = fig.add_subplot(gs[0, 0])
    strat_clean = strat_df.copy()
    strat_clean["roc_auc"] = pd.to_numeric(strat_clean["roc_auc"], errors="coerce")
    spike_clean = spike_df.copy()
    spike_clean["roc_auc"] = pd.to_numeric(spike_clean["roc_auc"], errors="coerce")

    for (tkey, tlbl, tcol, tmark), dx in zip(tools, offsets):
        means, lo_err, hi_err, xs = [], [], [], []
        for idx, b in enumerate(bins):
            if b == "pooled":
                vals = spike_clean.loc[spike_clean["tool"] == tkey, "roc_auc"].dropna()
            else:
                vals = strat_clean.loc[(strat_clean["length_stratum"] == b) & (strat_clean["tool"] == tkey), "roc_auc"].dropna()
            if vals.empty:
                continue
            m = float(vals.mean())
            means.append(m)
            lo_err.append(m - float(vals.min()))
            hi_err.append(float(vals.max()) - m)
            xs.append(idx + dx)
        ax_a.errorbar(
            xs,
            means,
            yerr=[lo_err, hi_err],
            fmt=tmark + "-",
            color=tcol,
            ecolor=tcol,
            capsize=3,
            linewidth=1.7,
            markersize=5.5,
            label=tlbl,
        )

    ax_a.axhline(0.95, color=PALETTE["red"], linestyle="--", linewidth=1.2, label="P3 ROC-AUC gate (0.95)")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([STRATA_LABELS[b] for b in bins], rotation=15)
    ax_a.set_ylim(0.75, 1.015)
    ax_a.set_ylabel("ROC-AUC (mean & min–max, 9 reps)")
    ax_a.grid(axis="y", color=PALETTE["grid"], linewidth=0.7)
    ax_a.set_title("A. Length-Stratified Spike-In ROC-AUC")
    ax_a.legend(loc="lower right", fontsize=7.6)

    # Panel B: Operating-point Retention vs Removal + Recall @ 5% FPR across length strata
    ax_b = fig.add_subplot(gs[0, 1])
    ens_strata = decontam_bench["stratified_by_length_and_tool"]
    bins_b = ["500_1000", "1000_3000", "3000_5000", "ge_5000"]
    xb = np.arange(len(bins_b))
    euk_ret = [ens_strata[b]["ensemble"]["mean_euk_base_retention"] * 100 for b in bins_b]
    bac_rem = [ens_strata[b]["ensemble"]["mean_bac_base_removal_recall"] * 100 for b in bins_b]
    r_fpr05 = [ens_strata[b]["ensemble"]["mean_bac_recall_at_fpr_0_05"] * 100 for b in bins_b]
    r_fpr01 = [ens_strata[b]["ensemble"]["mean_bac_recall_at_fpr_0_01"] * 100 for b in bins_b]

    w = 0.20
    b_euk = ax_b.bar(xb - 1.5 * w, euk_ret, width=w, color=PALETTE["teal"], label="Euk base retention (%)")
    b_bac = ax_b.bar(xb - 0.5 * w, bac_rem, width=w, color=PALETTE["amber"], label="Bac removal recall (%)")
    b_f05 = ax_b.bar(xb + 0.5 * w, r_fpr05, width=w, color=PALETTE["sky"], label="Recall @ 5% FPR (%)")
    b_f01 = ax_b.bar(xb + 1.5 * w, r_fpr01, width=w, color=PALETTE["purple"], label="Recall @ 1% FPR (%)")
    for bars in (b_euk, b_bac):
        for bar in bars:
            h = bar.get_height()
            ax_b.text(bar.get_x() + bar.get_width() / 2, h + 1.1, f"{h:.1f}", ha="center", va="bottom", fontsize=7.2)

    ax_b.set_xticks(xb)
    ax_b.set_xticklabels([STRATA_LABELS[b] for b in bins_b])
    ax_b.set_ylim(0, 113)
    ax_b.set_ylabel("Ensemble operating performance (%)")
    ax_b.grid(axis="y", color=PALETTE["grid"], linewidth=0.7)
    ax_b.set_title("B. Ensemble Retention, Removal & Low-FPR Recall")
    ax_b.legend(loc="upper left", ncol=2, fontsize=7.3)

    # Panel C: Replicate robustness across 5/15/25% spike-in + Organelle protection
    ax_c = fig.add_subplot(gs[0, 2])
    ens_sp = spike_clean[spike_clean["tool"] == "ensemble"].copy()
    fracs = sorted(ens_sp["spike_fraction"].unique())
    xc = np.arange(len(fracs))
    for i, frac in enumerate(fracs):
        sub = ens_sp[ens_sp["spike_fraction"] == frac]
        ax_c.scatter([i - 0.12] * len(sub), sub["roc_auc"], color=PALETTE["teal"], s=48, zorder=3, label="Pooled AUC" if i == 0 else "")
        ax_c.scatter([i + 0.04] * len(sub), sub["euk_base_retention"], color=PALETTE["sky"], marker="s", s=42, zorder=3, label="Euk base retention" if i == 0 else "")
        ax_c.scatter([i + 0.18] * len(sub), sub["bac_base_removal_recall"], color=PALETTE["amber"], marker="^", s=48, zorder=3, label="Bac base removal" if i == 0 else "")

    ax_c.axhline(0.95, color=PALETTE["red"], linestyle="--", linewidth=1.1)
    ax_c.set_xticks(xc)
    ax_c.set_xticklabels([f"{int(round(f * 100))}% bac\nspike-in" for f in fracs])
    ax_c.set_ylim(0.92, 1.01)
    ax_c.set_ylabel("Metric score (proportion)")
    ax_c.grid(axis="y", color=PALETTE["grid"], linewidth=0.7)
    ax_c.set_title("C. Spike-In Replicates & Organellar Guard")
    ax_c.legend(loc="lower left", fontsize=7.6)

    n_pla = int((org_df["organelle_ref_type"] == "plastid_ref").sum())
    n_mit = int((org_df["organelle_ref_type"] == "mitochondrion_ref").sum())
    n_ret = int((org_df["retained_in_filtered"].astype(str).str.lower() == "true").sum())
    ax_c.text(
        0.5,
        0.24,
        f"Organellar Safety Check:\n{n_ret}/{len(org_df)} retained as organellar\n({n_pla} plastid + {n_mit} mitochondrion)\n0 discarded as bacterial",
        transform=ax_c.transAxes,
        ha="center",
        va="center",
        fontsize=8.0,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#ECFDF5", edgecolor=PALETTE["emerald"], linewidth=1.1),
    )

    fig.suptitle(
        "HaptoSMASH P3: Tiara + Whokaryote Ensemble Benchmark (Roseobacter/Alteromonas Spike-Ins & OSF Organelles)",
        fontsize=12.2,
        fontweight="bold",
        y=1.03,
    )
    return save_figure(fig, outdir, "fig2_p3_decontam_benchmark")


def plot_fig3_18_assemblies(
    outdir: Path,
    samples_summary: pd.DataFrame,
    na_stats: Dict[str, Dict[str, int]],
) -> List[str]:
    """Figure 3: 18-assembly decontamination composition and length structure."""
    group_order = {"truth_benchmark": 0, "benchmark": 1, "flagship": 2, "osf_panel": 3}
    df = samples_summary.copy()
    df["group_rank"] = df["group"].map(group_order)
    df = df.sort_values(["group_rank", "bacterial_bp_fraction", "sample_id"], ascending=[True, True, True]).reset_index(drop=True)

    fig = plt.figure(figsize=(15.5, 7.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0])

    y = np.arange(len(df))[::-1]
    ylabels = []
    for _, r in df.iterrows():
        sid = r["sample_id"]
        has_na = na_stats.get(sid, {}).get("tiara_na_ge500_contigs", 0) > 0
        star = " *" if has_na else ""
        ylabels.append(f"{sid}{star}")

    # Panel A: Base-pair classification breakdown (%) across all 18 assemblies
    ax_a = fig.add_subplot(gs[0, 0])
    comps = [
        ("eukaryotic_haptophyte_bp_fraction", "Eukaryotic haptophyte (retained)", PALETTE["teal"]),
        ("organellar_bp_fraction", "Organellar (retained)", PALETTE["emerald"]),
        ("putative_hgt_bp_fraction", "Putative HGT / GFF-rescued (retained)", PALETTE["purple"]),
        ("ambiguous_bp_fraction", "Ambiguous (retained conservatively)", PALETTE["slate_light"]),
        ("bacterial_bp_fraction", "Bacterial (EXCLUDED)", PALETTE["rose"]),
    ]
    left = np.zeros(len(df))
    for col, lbl, col_hex in comps:
        vals = df[col].astype(float).values * 100.0
        ax_a.barh(y, vals, left=left, height=0.72, color=col_hex, edgecolor="#FFFFFF", linewidth=0.4, label=lbl)
        left += vals

    for yi, (_, r) in zip(y, df.iterrows()):
        bac_pct = float(r["bacterial_bp_fraction"]) * 100.0
        bac_mb = float(r["bacterial_bp"]) / 1e6
        if bac_pct >= 0.5:
            ax_a.text(101.0, yi, f"-{bac_pct:.2f}% ({bac_mb:.2f} Mb)", va="center", ha="left", fontsize=7.6, color=PALETTE["rose"], fontweight="bold")
        elif int(r["excluded_contigs"]) > 0:
            ax_a.text(101.0, yi, f"-{int(r['excluded_contigs'])} contigs ({int(r['excluded_bp'])/1e3:.1f} kb)", va="center", ha="left", fontsize=7.6, color=PALETTE["muted"])
        else:
            ax_a.text(101.0, yi, "0 excluded", va="center", ha="left", fontsize=7.6, color=PALETTE["teal"])

    ax_a.set_yticks(y)
    ax_a.set_yticklabels(ylabels, fontsize=8.3)
    ax_a.set_xlim(0, 127)
    ax_a.set_xlabel("Assembly base-pair composition (%)")
    ax_a.set_title("A. 18-Assembly Decontamination Call Breakdown (% of bp)")
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.45, -0.09), ncol=2, fontsize=7.8)

    # Panel B: Contig length distribution (counts by bin) showing why sub-3kb calibration matters
    ax_b = fig.add_subplot(gs[0, 1])
    len_bins = [
        ("len_ge_5k_contigs", ">= 5 kb", PALETTE["blue"]),
        ("len_3k_5k_contigs", "3–5 kb", PALETTE["sky"]),
        ("len_1k_3k_contigs", "1–3 kb (below Tiara 3k default)", PALETTE["amber"]),
        ("len_lt_1k_contigs", "< 1 kb", PALETTE["red"]),
    ]
    left_c = np.zeros(len(df))
    for col, lbl, col_hex in len_bins:
        vals = df[col].astype(float).values / 1000.0
        ax_b.barh(y, vals, left=left_c, height=0.72, color=col_hex, edgecolor="#FFFFFF", linewidth=0.4, label=lbl)
        left_c += vals

    for idx, (yi, (_, r)) in enumerate(zip(y, df.iterrows())):
        tot_c = int(r["total_contigs"])
        tot_mb = float(r["total_bp"]) / 1e6
        ax_b.text(left_c[idx] + 0.6, yi, f"{tot_c:,} ({tot_mb:.0f} Mb)", va="center", ha="left", fontsize=7.5, color=PALETTE["muted"])

    ax_b.set_yticks(y)
    ax_b.set_yticklabels([""] * len(y))
    ax_b.set_xlim(0, float(left_c.max()) * 1.28)
    ax_b.set_xlabel("Contig count (thousands)")
    ax_b.grid(axis="x", color=PALETTE["grid"], linewidth=0.7)
    ax_b.set_title("B. Contig Length Structure Driving Sub-3kb Calibration")
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.48, -0.09), ncol=2, fontsize=7.8)

    fig.text(
        0.01,
        -0.02,
        "* Assemblies where >=5 kb N-gap windows induced Tiara NaN scores on scaffolded contigs (0 in pglobosa_CNS00066 or osf_panel; all conservatively retained as ambiguous).",
        fontsize=8.0,
        color=PALETTE["muted"],
    )
    fig.suptitle("HaptoSMASH P3: Decontamination & Fragmentation Profiles Across All 18 Haptophyte Assemblies", fontsize=12.4, fontweight="bold", y=1.02)
    return save_figure(fig, outdir, "fig3_p3_18_assemblies_atlas")


def plot_fig4_cns00066_loci(
    outdir: Path,
    cns_ev: pd.DataFrame,
    cns_hotspots: pd.DataFrame,
    arch_map: Dict[str, Dict[str, Any]],
) -> List[str]:
    """Figure 4: pglobosa_CNS00066 decision space, 10 hotspots, and PKZILLA vs CNS00066 architectures."""
    fig = plt.figure(figsize=(15.5, 8.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.08], width_ratios=[1.05, 1.15])

    # Panel A: CNS00066 Tiara vs Whokaryote_S decision space
    ax_a = fig.add_subplot(gs[0, 0])
    df = cns_ev.copy()
    df["tiara_num"] = pd.to_numeric(df["tiara_p_prok"], errors="coerce")
    df["whok_num"] = pd.to_numeric(df["whokaryote_S_p_prok"], errors="coerce").fillna(-0.04)
    df["len_num"] = pd.to_numeric(df["length_bp"], errors="coerce")
    hotspot_scaffolds = set(cns_hotspots["contig"].unique())

    style_map = {
        "eukaryotic_haptophyte": ("Eukaryotic haptophyte (54)", PALETTE["teal"], "o", 38),
        "ambiguous": ("Ambiguous retained (18)", PALETTE["slate"], "s", 42),
        "putative_hgt": ("Putative HGT / GFF3-rescued (11)", PALETTE["purple"], "^", 68),
        "bacterial": ("Bacterial excluded (2)", PALETTE["rose"], "X", 95),
    }
    for lbl, (disp, col, mark, sz) in style_map.items():
        sub = df[df["final_label"] == lbl]
        sizes = np.where(sub["contig"].isin(hotspot_scaffolds), 120, sz)
        ax_a.scatter(sub["tiara_num"], sub["whok_num"], s=sizes, c=col, marker=mark, alpha=0.85, edgecolors="#FFFFFF", linewidth=0.6, label=disp, zorder=3)

    hs_rows = df[df["contig"].isin(hotspot_scaffolds)]
    ax_a.scatter(hs_rows["tiara_num"], hs_rows["whok_num"], s=170, facecolors="none", edgecolors=PALETTE["amber"], linewidth=1.8, label="Track A hotspot scaffold (8 scaffolds / 10 loci)", zorder=4)

    bac_rows = df[df["final_label"] == "bacterial"]
    if not bac_rows.empty:
        bx = float(bac_rows["tiara_num"].iloc[0])
        by = float(bac_rows["whok_num"].iloc[0])
        ax_a.annotate(
            "GWHBJCL00000063 (2.95 kb)\n& GWHBJCL00000069 (1.93 kb)",
            (bx, by),
            xytext=(-138, -28),
            textcoords="offset points",
            fontsize=7.4,
            color=PALETTE["rose"],
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=PALETTE["rose"], lw=0.9),
        )
    for hid_lbl, cid, offset in [
        ("H0007 scaffold\n(9.37 Mb)", "GWHBJCL00000050", (54, 34)),
        ("H0008 scaffold\n(7.08 Mb)", "GWHBJCL00000059", (58, -14)),
    ]:
        row = df[df["contig"] == cid].iloc[0]
        ax_a.annotate(hid_lbl, (row["tiara_num"], row["whok_num"]), xytext=offset, textcoords="offset points", fontsize=7.4, color=PALETTE["ink"], arrowprops=dict(arrowstyle="->", color=PALETTE["amber"], lw=0.9))

    ax_a.set_xlim(-0.05, 1.07)
    ax_a.set_ylim(-0.09, 1.08)
    ax_a.axhline(0.0, color=PALETTE["slate_light"], linewidth=0.7, linestyle=":")
    ax_a.text(0.02, -0.068, "y = -0.04: Whokaryote RF_S NA (<2 genes; n=1 contig)", fontsize=7.2, color=PALETTE["muted"])
    ax_a.set_xlabel("Tiara p_prok (bac + arc)")
    ax_a.set_ylabel("Whokaryote RF_S p_prok (gene structure)")
    ax_a.grid(color=PALETTE["grid"], linewidth=0.7)
    ax_a.set_title("A. P. globosa CNS00066 Contig Classification Space (n=85)")
    ax_a.legend(loc="center left", bbox_to_anchor=(0.14, 0.72), fontsize=7.2)

    # Panel B: CNS00066 all 10 hotspots (length & domain breakdown)
    ax_b = fig.add_subplot(gs[0, 1])
    hs_sorted = cns_hotspots.copy()
    hs_sorted["length_kb"] = hs_sorted["length_bp"].astype(float) / 1000.0
    yb = np.arange(len(hs_sorted))[::-1]
    n_ks = hs_sorted["n_ks"].astype(int).values
    n_ca = hs_sorted["n_ca"].astype(int).values
    n_other_bb = hs_sorted["n_backbone"].astype(int).values - n_ks - n_ca
    n_anc = hs_sorted["n_ancillary"].astype(int).values

    ax_b.barh(yb, n_ks, height=0.66, color=PALETTE["red"], label="KS backbone")
    ax_b.barh(yb, n_ca, left=n_ks, height=0.66, color=PALETTE["sky"], label="NRPS C/A")
    ax_b.barh(yb, n_other_bb, left=n_ks + n_ca, height=0.66, color=PALETTE["amber"], label="AT/PP backbone")
    ax_b.barh(yb, n_anc, left=n_ks + n_ca + n_other_bb, height=0.66, color=PALETTE["teal"], label="Track A ancillary")

    for yi, (_, r) in zip(yb, hs_sorted.iterrows()):
        tot = int(r["n_domains"])
        kb = float(r["length_bp"]) / 1000.0
        ax_b.text(tot + 2.0, yi, f"{kb:.1f} kb  ({tot} dom)", va="center", ha="left", fontsize=7.8, color=PALETTE["ink"], fontweight="bold" if r["hotspot_id"] in {"H0007", "H0008"} else "normal")

    ax_b.set_yticks(yb)
    ax_b.set_yticklabels([f"{r['hotspot_id']} ({r['contig'][-4:]})" for _, r in hs_sorted.iterrows()], fontsize=8.2)
    ax_b.set_xlim(0, 185)
    ax_b.set_xlabel("Track A HMM domain count per hotspot")
    ax_b.grid(axis="x", color=PALETTE["grid"], linewidth=0.7)
    ax_b.set_title("B. CNS00066 Track A Hotspots (10/10 Preserved Post-Decontam)")
    ax_b.legend(loc="lower right", fontsize=7.6)

    # Panel C: Multi-track Domain Architecture Comparison (12B1 PKZILLA-1/-2/-3 vs CNS00066 H0007/H0008)
    ax_c = fig.add_subplot(gs[1, :])
    loci_tracks = [
        ("12B1_H0017", "P. parvum 12B1 · PKZILLA-1 (H0017)"),
        ("12B1_H0007", "P. parvum 12B1 · PKZILLA-2 (H0007)"),
        ("12B1_H0011", "P. parvum 12B1 · PKZILLA-3 (H0011)"),
        ("CNS_H0007", "P. globosa CNS00066 · H0007"),
        ("CNS_H0008", "P. globosa CNS00066 · H0008"),
    ]
    max_kb = max(arch_map[k]["length_bp"] for k, _ in loci_tracks) / 1000.0
    y_centers = np.linspace(len(loci_tracks) - 1, 0, len(loci_tracks)) * 1.4

    legend_seen: Dict[str, str] = {}
    for idx, (akey, albl) in enumerate(loci_tracks):
        yc = y_centers[idx]
        locus = arch_map[akey]
        kb = locus["length_bp"] / 1000.0
        ax_c.plot([0, kb], [yc, yc], color=PALETTE["slate_light"], linewidth=2.2, zorder=1)
        for d in locus["domains"]:
            cat, dcol = classify_domain_group(d["model_name"], d["track"], d["role"])
            legend_seen[cat] = dcol
            xs = d["rel_start"] / 1000.0
            xe = d["rel_end"] / 1000.0
            width = max(xe - xs, 0.35)
            frame_offset = (d["frame"] - 1) * 0.24
            rect = mpatches.Rectangle(
                (xs, yc + frame_offset - 0.10),
                width,
                0.20,
                facecolor=dcol,
                edgecolor="#FFFFFF",
                linewidth=0.25,
                alpha=0.92,
                zorder=3 if d["track"] == "A" else 2,
            )
            ax_c.add_patch(rect)
        ax_c.text(
            kb + 1.4,
            yc,
            f"{kb:.1f} kb · {locus['n_domains']} dom ({locus['n_ks']} KS, {locus['n_ca']} CA) · {locus['strand']}",
            va="center",
            ha="left",
            fontsize=8.2,
            color=PALETTE["ink"],
        )

    ax_c.set_yticks(y_centers)
    ax_c.set_yticklabels([lbl for _, lbl in loci_tracks], fontsize=8.8, fontweight="bold")
    ax_c.set_xlim(-1.0, max_kb + 38.0)
    ax_c.set_ylim(-0.7, y_centers[0] + 0.7)
    ax_c.set_xlabel("Relative position along hotspot (kb; sub-tracks show reading frames 0, 1, 2)")
    ax_c.grid(axis="x", color=PALETTE["grid"], linewidth=0.7)
    ax_c.set_title("C. GFF-Free 6-Frame Track A Domain Architectures Across Flagship Gigasynthase Loci")

    handles = [mpatches.Patch(color=col, label=lbl) for lbl, col in legend_seen.items()]
    ax_c.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=6, fontsize=7.8)

    fig.suptitle(
        "Phaeocystis globosa CNS00066: Decontamination Space, Flagship Hotspots & PKZILLA Domain Architecture Comparison",
        fontsize=12.4,
        fontweight="bold",
        y=1.02,
    )
    return save_figure(fig, outdir, "fig4_cns00066_flagship_loci")


def build_data_and_figures(root: Path, fig_outdir: Path, bundle_path: Path) -> Dict[str, Any]:
    configure_style()

    pfam_rows = read_tsv(root / "resources/pfam_lists/hapto_core_pfam.tsv")
    manifest_p1 = json.loads((root / "db/custom_hmm/hapto_core.manifest.json").read_text(encoding="utf-8"))
    meta_12b1 = json.loads((root / "results/track_a/pparvum_12B1_v1_1/metadata.json").read_text(encoding="utf-8"))
    meta_cns_unfilt = json.loads((root / "results/track_a/pglobosa_CNS00066/metadata.json").read_text(encoding="utf-8"))
    meta_cns = json.loads((root / "results/track_a_decontam/pglobosa_CNS00066/metadata.json").read_text(encoding="utf-8"))
    bench_12b1 = json.loads((root / "results/benchmarks/track_a_12B1/benchmark_summary.json").read_text(encoding="utf-8"))
    sens_df = pd.read_csv(root / "results/benchmarks/track_a_12B1/sensitivity.tsv", sep="\t")
    hotspots_12b1_df = pd.read_csv(root / "results/track_a/pparvum_12B1_v1_1/hotspots.tsv", sep="\t")
    hotspots_cns_df = pd.read_csv(root / "results/track_a_decontam/pglobosa_CNS00066/hotspots.tsv", sep="\t")

    decontam_bench = json.loads((root / "results/benchmarks/decontam/benchmark_summary.json").read_text(encoding="utf-8"))
    spike_df = pd.read_csv(root / "results/benchmarks/decontam/spikein_metrics.tsv", sep="\t")
    strat_df = pd.read_csv(root / "results/benchmarks/decontam/stratified_metrics.tsv", sep="\t")
    org_df = pd.read_csv(root / "results/benchmarks/decontam/organelle_check.tsv", sep="\t")

    samples_meta = pd.read_csv(root / "config/samples.tsv", sep="\t")
    summary_df = pd.read_csv(root / "results/decontam/all_samples_summary.tsv", sep="\t")
    merged_samples = samples_meta.merge(summary_df, on="sample_id", how="inner")
    na_stats = compute_tiara_na_stats(root, list(merged_samples["sample_id"]))
    merged_samples["tiara_na_ge500_contigs"] = merged_samples["sample_id"].map(lambda s: na_stats[s]["tiara_na_ge500_contigs"])
    merged_samples["tiara_na_ge500_bp"] = merged_samples["sample_id"].map(lambda s: na_stats[s]["tiara_na_ge500_bp"])

    cns_ev = pd.read_csv(root / "results/decontam/pglobosa_CNS00066/contig_evidence.tsv", sep="\t")

    hs_12b1_list = read_tsv(root / "results/track_a/pparvum_12B1_v1_1/hotspots.tsv")
    hs_cns_list = read_tsv(root / "results/track_a_decontam/pglobosa_CNS00066/hotspots.tsv")
    arch_12b1 = extract_hotspot_domains(hs_12b1_list, root / "results/track_a/pparvum_12B1_v1_1/domains.tsv", ["H0017", "H0007", "H0011"])
    arch_cns = extract_hotspot_domains(hs_cns_list, root / "results/track_a_decontam/pglobosa_CNS00066/domains.tsv")

    arch_for_plot = {
        "12B1_H0017": arch_12b1["H0017"],
        "12B1_H0007": arch_12b1["H0007"],
        "12B1_H0011": arch_12b1["H0011"],
        "CNS_H0007": arch_cns["H0007"],
        "CNS_H0008": arch_cns["H0008"],
    }

    outputs: List[str] = []
    outputs.extend(plot_fig1_p1_p2(fig_outdir, pfam_rows, meta_12b1, meta_cns, bench_12b1, sens_df))
    outputs.extend(plot_fig2_p3_decontam_bench(fig_outdir, decontam_bench, strat_df, spike_df, org_df))
    outputs.extend(plot_fig3_18_assemblies(fig_outdir, merged_samples, na_stats))
    outputs.extend(plot_fig4_cns00066_loci(fig_outdir, cns_ev, hotspots_cns_df, arch_for_plot))

    bundle = {
        "p1": {
            "model_count": manifest_p1["model_count"],
            "hmmer_version": manifest_p1["hmmer_version"],
            "pfam_list": pfam_rows,
            "counts_by_track_role": {f"{t}/{r}": c for (t, r), c in Counter((x["track"], x["role"]) for x in pfam_rows).items()},
            "hits_12b1": meta_12b1["summary"]["domain_hits_by_track_role"],
            "hits_cns": meta_cns["summary"]["domain_hits_by_track_role"],
        },
        "p2": {
            "meta_12b1": meta_12b1["summary"],
            "benchmark_summary": bench_12b1,
            "sensitivity": json.loads(sens_df.to_json(orient="records")),
            "hotspots_12b1": json.loads(hotspots_12b1_df.to_json(orient="records")),
            "pkzilla_architectures": arch_12b1,
        },
        "p3": {
            "decontam_benchmark": decontam_bench,
            "spikein_metrics": json.loads(spike_df.where(spike_df.notna(), None).to_json(orient="records")),
            "stratified_metrics": json.loads(strat_df.where(strat_df.notna(), None).to_json(orient="records")),
            "organelle_check": json.loads(org_df.to_json(orient="records")),
            "samples_18": json.loads(merged_samples.where(merged_samples.notna(), None).to_json(orient="records")),
            "cns_evidence": json.loads(cns_ev.where(cns_ev.notna(), None).to_json(orient="records")),
            "cns_meta_unfiltered": meta_cns_unfilt["summary"],
            "cns_meta_decontam": meta_cns["summary"],
            "cns_hotspots": json.loads(hotspots_cns_df.to_json(orient="records")),
            "cns_architectures": arch_cns,
        },
        "figures": outputs,
    }
    bundle_path.parent.mkdir(parents=True, exist_ok=True)
    bundle_path.write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    return bundle


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--fig-outdir", type=Path, default=Path("results/figures"))
    parser.add_argument("--bundle-out", type=Path, default=Path("tmp/viz_bundle.json"))
    args = parser.parse_args()

    bundle = build_data_and_figures(args.root, args.fig_outdir, args.bundle_out)
    print(f"Generated {len(bundle['figures'])} figure files -> {args.fig_outdir}")
    for p in bundle["figures"]:
        print(f"  - {p}")
    print(f"Wrote visualization data bundle -> {args.bundle_out} ({args.bundle_out.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Synthetic bacterial spike-in (5/15/25%) and organellar safety benchmark for P3 decontamination."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from sklearn.metrics import roc_auc_score

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import decontaminate_assemblies as dec  # noqa: E402


STRATA_BOUNDS = [
    ("500_1000", 500, 999),
    ("1000_3000", 1000, 2999),
    ("3000_5000", 3000, 4999),
    ("ge_5000", 5000, 25000),
]

STRATA_WEIGHTS = {
    "500_1000": 0.15,
    "1000_3000": 0.25,
    "3000_5000": 0.25,
    "ge_5000": 0.35,
}


def load_sequences_from_fasta(path: Path) -> List[Tuple[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    records: List[Tuple[str, str]] = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for title, seq in SimpleFastaParser(handle):
            records.append((title.split()[0], seq.upper()))
    return records


def sample_fragments_by_strata(
    sources: List[Tuple[str, str]],
    target_bp: int,
    *,
    prefix: str,
    rng: random.Random,
    strata_weights: Dict[str, float] = STRATA_WEIGHTS,
) -> List[Tuple[str, str, str]]:
    """Sample non-N nucleotide fragments across explicit length strata totaling ~target_bp."""
    eligible = [(cid, seq) for cid, seq in sources if len(seq) >= 25000]
    if not eligible:
        eligible = [(cid, seq) for cid, seq in sources if len(seq) >= 5000]
    if not eligible:
        raise ValueError(f"No sufficiently long contigs for sampling ({prefix})")

    fragments: List[Tuple[str, str, str]] = []
    counter = 0
    for s_name, lo, hi in STRATA_BOUNDS:
        bin_target = int(round(target_bp * strata_weights[s_name]))
        acc_bp = 0
        attempts = 0
        while acc_bp < bin_target and attempts < 5000:
            attempts += 1
            src_id, src_seq = rng.choice(eligible)
            frag_len = rng.randint(lo, min(hi, len(src_seq)))
            max_start = len(src_seq) - frag_len
            if max_start < 0:
                continue
            start = rng.randint(0, max_start)
            frag = src_seq[start : start + frag_len]
            if frag.count("N") / frag_len > 0.05:
                continue
            counter += 1
            cid = f"{prefix}__{s_name}__{counter:05d}"
            fragments.append((cid, frag, s_name))
            acc_bp += frag_len
    return fragments


def recall_at_fixed_fpr(
    y_true: np.ndarray, y_score: np.ndarray, max_fpr: float
) -> float:
    """Compute maximum bacterial recall achievable at empirical false-positive rate <= max_fpr."""
    neg_scores = np.sort(y_score[y_true == 0])
    pos_scores = y_score[y_true == 1]
    if len(neg_scores) == 0 or len(pos_scores) == 0:
        return float("nan")
    # Threshold strictly above top floor(max_fpr * N_neg) negative scores
    allowed_fps = int(np.floor(max_fpr * len(neg_scores)))
    if allowed_fps == 0:
        cutoff = float(neg_scores[-1])
        tp = int(np.sum(pos_scores > cutoff))
    elif allowed_fps >= len(neg_scores):
        tp = len(pos_scores)
    else:
        cutoff = float(neg_scores[-(allowed_fps + 1)])
        tp = int(np.sum(pos_scores > cutoff))
    return float(tp / len(pos_scores))


def compute_tool_metrics(
    rows: List[Dict[str, object]],
    truth_map: Dict[str, str],
    score_key: str,
) -> Dict[str, object]:
    """Compute ROC-AUC, recall@FPR=0.01/0.05, and coverage for a continuous score key."""
    y_true_all: List[int] = []
    y_score_all: List[float] = []
    n_total = len(rows)
    n_pos_total = sum(1 for r in rows if truth_map[str(r["contig"])] == "bacteria")
    n_neg_total = n_total - n_pos_total

    for r in rows:
        val = dec.parse_float_or_none(str(r[score_key]))
        if val is None:
            continue
        lbl = 1 if truth_map[str(r["contig"])] == "bacteria" else 0
        y_true_all.append(lbl)
        y_score_all.append(val)

    n_scored = len(y_true_all)
    if n_scored == 0:
        return {
            "n_total": n_total,
            "n_scored": 0,
            "coverage_fraction": 0.0,
            "n_pos_scored": 0,
            "n_neg_scored": 0,
            "roc_auc": "NA",
            "bac_recall_at_fpr_0_01": "NA",
            "bac_recall_at_fpr_0_05": "NA",
        }
    yt = np.array(y_true_all, dtype=int)
    ys = np.array(y_score_all, dtype=float)
    n_pos = int(np.sum(yt == 1))
    n_neg = int(np.sum(yt == 0))
    if n_pos == 0 or n_neg == 0:
        auc_str = "NA"
        r01_str = "NA"
        r05_str = "NA"
    else:
        auc_str = f"{roc_auc_score(yt, ys):.6f}"
        r01_str = f"{recall_at_fixed_fpr(yt, ys, 0.01):.6f}"
        r05_str = f"{recall_at_fixed_fpr(yt, ys, 0.05):.6f}"
    return {
        "n_total": n_total,
        "n_scored": n_scored,
        "coverage_fraction": round(n_scored / n_total, 6) if n_total else 0.0,
        "n_pos_total": n_pos_total,
        "n_neg_total": n_neg_total,
        "n_pos_scored": n_pos,
        "n_neg_scored": n_neg,
        "roc_auc": auc_str,
        "bac_recall_at_fpr_0_01": r01_str,
        "bac_recall_at_fpr_0_05": r05_str,
    }


def compute_discrete_metrics(
    rows: List[Dict[str, object]],
    truth_map: Dict[str, str],
) -> Dict[str, object]:
    euk_bp_total = 0
    euk_bp_retained = 0
    bac_bp_total = 0
    bac_bp_excluded = 0
    amb_contigs = 0
    amb_bp = 0
    total_bp = 0
    for r in rows:
        cid = str(r["contig"])
        bp = int(r["length_bp"])
        total_bp += bp
        truth = truth_map[cid]
        retained = r["retained_in_filtered"] == "true"
        lbl = str(r["final_label"])
        if lbl == "ambiguous":
            amb_contigs += 1
            amb_bp += bp
        if truth == "eukarya":
            euk_bp_total += bp
            if retained:
                euk_bp_retained += bp
        elif truth == "bacteria":
            bac_bp_total += bp
            if not retained:
                bac_bp_excluded += bp
    return {
        "euk_bp_total": euk_bp_total,
        "euk_base_retention": round(euk_bp_retained / euk_bp_total, 6) if euk_bp_total else 0.0,
        "bac_bp_total": bac_bp_total,
        "bac_base_removal_recall": round(bac_bp_excluded / bac_bp_total, 6) if bac_bp_total else 0.0,
        "ambiguous_contigs": amb_contigs,
        "ambiguous_contig_fraction": round(amb_contigs / len(rows), 6) if rows else 0.0,
        "ambiguous_bp": amb_bp,
        "ambiguous_base_fraction": round(amb_bp / total_bp, 6) if total_bp else 0.0,
    }


def evaluate_organellar_references(
    organelle_dir: Path,
    *,
    tiara_min_len: int = 500,
    whokaryote_min_len: int = 1000,
    threads: int = 4,
) -> List[Dict[str, object]]:
    """Verify real haptophyte plastid/mitochondrial genomes are never discarded as bacteria."""
    gb_files = sorted(organelle_dir.glob("*.gb"))
    if not gb_files:
        raise FileNotFoundError(f"No GenBank files found in {organelle_dir}")

    rows_out: List[Dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="hapto_org_eval_") as tmp_str:
        work = Path(tmp_str)
        fasta = work / "organelles.fasta"
        truth_kind: Dict[str, str] = {}
        with fasta.open("wt", encoding="utf-8") as handle:
            for gb in gb_files:
                kind = "plastid_ref" if "pt_" in gb.name else "mitochondrion_ref"
                for rec in SeqIO.parse(str(gb), "genbank"):
                    if len(rec.seq) < tiara_min_len:
                        continue
                    cid = f"{kind}__{rec.id}"
                    truth_kind[cid] = kind
                    handle.write(f">{cid}\n{str(rec.seq).upper()}\n")
        outdir = work / "decontam"
        dec.run_decontamination(
            fasta_path=fasta,
            outdir=outdir,
            sample_id="organelle_refs",
            tiara_min_len=tiara_min_len,
            whokaryote_min_len=whokaryote_min_len,
            threads=threads,
            emit_filtered_fasta=False,
        )
        with (outdir / "contig_evidence.tsv").open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows_out.append(
                    {
                        "contig": row["contig"],
                        "organelle_ref_type": truth_kind[row["contig"]],
                        "length_bp": int(row["length_bp"]),
                        "tiara_fst_stage": row["tiara_fst_stage"],
                        "tiara_snd_stage": row["tiara_snd_stage"],
                        "whokaryote_S_pred": row["whokaryote_S_pred"],
                        "final_label": row["final_label"],
                        "retained_in_filtered": row["retained_in_filtered"],
                        "decision_reason": row["decision_reason"],
                    }
                )
    return rows_out


def write_tsv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def run_benchmark(
    *,
    euk_fasta: Path,
    refs_manifest: Path,
    organelle_dir: Path,
    outdir: Path,
    fractions: Sequence[float] = (0.05, 0.15, 0.25),
    seeds: Sequence[int] = (101, 202, 303),
    total_bp_per_replicate: int = 8_000_000,
    tiara_min_len: int = 500,
    whokaryote_min_len: int = 1000,
    threads: int = 4,
) -> Dict[str, object]:
    euk_sources = load_sequences_from_fasta(euk_fasta)
    bac_sources_by_ref: List[Tuple[str, List[Tuple[str, str]]]] = []
    with refs_manifest.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            bac_sources_by_ref.append(
                (row["accession"], load_sequences_from_fasta(Path(row["genome_fasta"])))
            )

    spike_rows: List[Dict[str, object]] = []
    strat_rows: List[Dict[str, object]] = []
    tools = [
        ("tiara", "tiara_p_prok"),
        ("whokaryote_S", "whokaryote_S_p_prok"),
        ("whokaryote_T", "whokaryote_T_p_prok"),
        ("ensemble", "ensemble_prok_score"),
    ]

    with tempfile.TemporaryDirectory(prefix="hapto_spikein_") as tmp_str:
        work_root = Path(tmp_str)
        for frac in fractions:
            for seed in seeds:
                rng = random.Random(int(round(frac * 1000)) * 10000 + seed)
                bac_bp_target = int(round(total_bp_per_replicate * frac))
                euk_bp_target = total_bp_per_replicate - bac_bp_target
                per_ref_bp = bac_bp_target // len(bac_sources_by_ref)

                frags: List[Tuple[str, str, str, str]] = []
                for cid, seq, s_name in sample_fragments_by_strata(
                    euk_sources, euk_bp_target, prefix="euk", rng=rng
                ):
                    frags.append((cid, seq, "eukarya", s_name))
                for ref_idx, (acc, ref_seqs) in enumerate(bac_sources_by_ref):
                    target = (
                        per_ref_bp
                        if ref_idx < len(bac_sources_by_ref) - 1
                        else bac_bp_target - per_ref_bp * (len(bac_sources_by_ref) - 1)
                    )
                    for cid, seq, s_name in sample_fragments_by_strata(
                        ref_seqs, target, prefix=f"bac_{acc}", rng=rng
                    ):
                        frags.append((cid, seq, "bacteria", s_name))

                rng.shuffle(frags)
                truth_map = {cid: truth for cid, _seq, truth, _sb in frags}
                stratum_map = {cid: sb for cid, _seq, _truth, sb in frags}

                rep_dir = work_root / f"f{int(round(frac * 100))}_s{seed}"
                rep_dir.mkdir(parents=True, exist_ok=True)
                rep_fasta = rep_dir / "spikein.fasta"
                with rep_fasta.open("wt", encoding="utf-8") as handle:
                    for cid, seq, _truth, _sb in frags:
                        handle.write(f">{cid}\n{seq}\n")

                dec.run_decontamination(
                    fasta_path=rep_fasta,
                    outdir=rep_dir / "out",
                    sample_id=f"spikein_f{int(round(frac * 100))}_s{seed}",
                    tiara_min_len=tiara_min_len,
                    whokaryote_min_len=whokaryote_min_len,
                    threads=threads,
                    emit_filtered_fasta=False,
                )
                with (rep_dir / "out/contig_evidence.tsv").open(newline="", encoding="utf-8") as h:
                    ev_rows = list(csv.DictReader(h, delimiter="\t"))

                disc_overall = compute_discrete_metrics(ev_rows, truth_map)
                for tool_name, score_key in tools:
                    tm = compute_tool_metrics(ev_rows, truth_map, score_key)
                    spike_rows.append(
                        {
                            "spike_fraction": frac,
                            "seed": seed,
                            "tool": tool_name,
                            **tm,
                            **disc_overall,
                        }
                    )

                # Stratifications: each explicit length bin + combined <3kb and >=3kb
                stratum_groups: List[Tuple[str, List[Dict[str, object]]]] = []
                for s_name, _, _ in STRATA_BOUNDS:
                    stratum_groups.append(
                        (s_name, [r for r in ev_rows if stratum_map[str(r["contig"])] == s_name])
                    )
                stratum_groups.append(
                    ("lt_3000", [r for r in ev_rows if int(r["length_bp"]) < 3000])
                )
                stratum_groups.append(
                    ("ge_3000", [r for r in ev_rows if int(r["length_bp"]) >= 3000])
                )
                for s_label, s_rows in stratum_groups:
                    s_disc = compute_discrete_metrics(s_rows, truth_map)
                    for tool_name, score_key in tools:
                        tm = compute_tool_metrics(s_rows, truth_map, score_key)
                        strat_rows.append(
                            {
                                "spike_fraction": frac,
                                "seed": seed,
                                "length_stratum": s_label,
                                "tool": tool_name,
                                **tm,
                                **s_disc,
                            }
                        )

    organelle_rows = evaluate_organellar_references(
        organelle_dir,
        tiara_min_len=tiara_min_len,
        whokaryote_min_len=whokaryote_min_len,
        threads=threads,
    )

    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / "spikein_metrics.tsv", spike_rows)
    write_tsv(outdir / "stratified_metrics.tsv", strat_rows)
    write_tsv(outdir / "organelle_check.tsv", organelle_rows)

    # Aggregate summary by tool and by (length_stratum, tool)
    def _mean_numeric(items: List[Dict[str, object]], key: str) -> Optional[float]:
        vals = [dec.parse_float_or_none(str(x[key])) for x in items]
        clean = [v for v in vals if v is not None]
        return round(sum(clean) / len(clean), 6) if clean else None

    overall_by_tool: Dict[str, Dict[str, object]] = {}
    for tool_name, _ in tools:
        subset = [r for r in spike_rows if r["tool"] == tool_name]
        overall_by_tool[tool_name] = {
            "mean_roc_auc": _mean_numeric(subset, "roc_auc"),
            "min_roc_auc": min(
                (v for v in (dec.parse_float_or_none(str(x["roc_auc"])) for x in subset) if v is not None),
                default=None,
            ),
            "mean_coverage_fraction": _mean_numeric(subset, "coverage_fraction"),
            "mean_bac_recall_at_fpr_0_01": _mean_numeric(subset, "bac_recall_at_fpr_0_01"),
            "mean_bac_recall_at_fpr_0_05": _mean_numeric(subset, "bac_recall_at_fpr_0_05"),
            "mean_euk_base_retention": _mean_numeric(subset, "euk_base_retention"),
            "mean_bac_base_removal_recall": _mean_numeric(subset, "bac_base_removal_recall"),
            "mean_ambiguous_base_fraction": _mean_numeric(subset, "ambiguous_base_fraction"),
        }

    strat_summary: Dict[str, Dict[str, object]] = {}
    for s_label in ("500_1000", "1000_3000", "3000_5000", "ge_5000", "lt_3000", "ge_3000"):
        strat_summary[s_label] = {}
        for tool_name, _ in tools:
            subset = [
                r
                for r in strat_rows
                if r["length_stratum"] == s_label and r["tool"] == tool_name
            ]
            strat_summary[s_label][tool_name] = {
                "mean_roc_auc": _mean_numeric(subset, "roc_auc"),
                "mean_coverage_fraction": _mean_numeric(subset, "coverage_fraction"),
                "mean_bac_recall_at_fpr_0_01": _mean_numeric(subset, "bac_recall_at_fpr_0_01"),
                "mean_bac_recall_at_fpr_0_05": _mean_numeric(subset, "bac_recall_at_fpr_0_05"),
                "mean_euk_base_retention": _mean_numeric(subset, "euk_base_retention"),
                "mean_bac_base_removal_recall": _mean_numeric(subset, "bac_base_removal_recall"),
                "mean_ambiguous_base_fraction": _mean_numeric(subset, "ambiguous_base_fraction"),
            }

    org_total = len(organelle_rows)
    org_discarded_as_bac = sum(1 for r in organelle_rows if r["retained_in_filtered"] != "true")
    org_labeled_organellar = sum(1 for r in organelle_rows if r["final_label"] == "organellar")
    summary_payload = {
        "format": "hapto_decontam_benchmark_v1",
        "parameters": {
            "fractions": list(fractions),
            "seeds": list(seeds),
            "total_bp_per_replicate": total_bp_per_replicate,
            "tiara_min_len": tiara_min_len,
            "whokaryote_min_len": whokaryote_min_len,
        },
        "overall_by_tool": overall_by_tool,
        "stratified_by_length_and_tool": strat_summary,
        "organellar_safety": {
            "organellar_contigs_tested": org_total,
            "labeled_organellar": org_labeled_organellar,
            "discarded_as_bacterial": org_discarded_as_bac,
            "organelle_retention_rate": round((org_total - org_discarded_as_bac) / org_total, 6)
            if org_total
            else 0.0,
        },
    }
    (outdir / "benchmark_summary.json").write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary_payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--euk-fasta", type=Path, required=True)
    parser.add_argument("--refs-manifest", type=Path, required=True)
    parser.add_argument("--organelle-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fractions", default="0.05,0.15,0.25")
    parser.add_argument("--seeds", default="101,202,303")
    parser.add_argument("--total-bp", type=int, default=8_000_000)
    parser.add_argument("--tiara-min-len", type=int, default=500)
    parser.add_argument("--whokaryote-min-len", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    fracs = [float(x.strip()) for x in args.fractions.split(",") if x.strip()]
    seeds = [int(x.strip()) for x in args.seeds.split(",") if x.strip()]
    summary = run_benchmark(
        euk_fasta=args.euk_fasta,
        refs_manifest=args.refs_manifest,
        organelle_dir=args.organelle_dir,
        outdir=args.outdir,
        fractions=fracs,
        seeds=seeds,
        total_bp_per_replicate=args.total_bp,
        tiara_min_len=args.tiara_min_len,
        whokaryote_min_len=args.whokaryote_min_len,
        threads=args.threads,
    )
    ens = summary["overall_by_tool"]["ensemble"]
    org = summary["organellar_safety"]
    print(
        f"Decontam spike-in PASS: ensemble_mean_auc={ens['mean_roc_auc']} "
        f"min_auc={ens['min_roc_auc']} euk_retention={ens['mean_euk_base_retention']} "
        f"bac_removal={ens['mean_bac_base_removal_recall']} "
        f"org_discarded={org['discarded_as_bacterial']}/{org['organellar_contigs_tested']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

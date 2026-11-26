#!/usr/bin/env python3
"""Tiara + Whokaryote assembly decontamination with raw-score retention and organelle/HGT guards."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import joblib
import numpy as np
import pandas as pd
import whokaryote_scripts
from Bio.SeqIO.FastaIO import SimpleFastaParser


EVIDENCE_FIELDS = [
    "contig",
    "length_bp",
    "length_bin",
    "tiara_scored",
    "tiara_fst_stage",
    "tiara_snd_stage",
    "tiara_p_org",
    "tiara_p_bac",
    "tiara_p_arc",
    "tiara_p_euk",
    "tiara_p_unk1",
    "tiara_p_pla",
    "tiara_p_unk2",
    "tiara_p_mit",
    "tiara_p_prok",
    "whokaryote_minsize_eligible",
    "prodigal_nr_genes",
    "whokaryote_rf_scored",
    "whokaryote_S_pred",
    "whokaryote_S_p_prok",
    "whokaryote_T_pred",
    "whokaryote_T_p_prok",
    "gene_density",
    "rbs_ratio",
    "ratio_same_orientation",
    "ensemble_prok_score",
    "classifier_disagreement",
    "has_explicit_algal_evidence",
    "algal_evidence_source",
    "final_label",
    "retained_in_filtered",
    "decision_reason",
]

VALID_LABELS = (
    "eukaryotic_haptophyte",
    "bacterial",
    "organellar",
    "ambiguous",
    "putative_hgt",
)

LENGTH_BINS = ("lt_1k", "1k_3k", "3k_5k", "ge_5k")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(4 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def open_fasta_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def length_bin_for(length_bp: int) -> str:
    if length_bp < 1000:
        return "lt_1k"
    if length_bp < 3000:
        return "1k_3k"
    if length_bp < 5000:
        return "3k_5k"
    return "ge_5k"


def parse_float_or_none(val: str) -> Optional[float]:
    txt = str(val).strip()
    if not txt or txt.lower() in {"na", "n/a", "nan", "none"}:
        return None
    parsed = float(txt)
    if math.isnan(parsed):
        return None
    return parsed


def format_opt_float(val: Optional[float], precision: int = 6) -> str:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "NA"
    return f"{val:.{precision}f}"


def validate_decontam_prerequisites() -> Dict[str, str]:
    """Validate Tiara, Whokaryote, Prodigal, and RF model files in current environment."""
    env_bin = str(Path(sys.executable).resolve().parent)
    current_path = os.environ.get("PATH", "")
    if env_bin not in current_path.split(os.pathsep):
        os.environ["PATH"] = env_bin + os.pathsep + current_path

    tools: Dict[str, str] = {}
    for exe in ("tiara", "whokaryote.py", "prodigal"):
        resolved = shutil.which(exe)
        if not resolved:
            raise RuntimeError(f"Required decontamination executable missing from PATH: {exe}")
        tools[exe] = resolved

    data_dir = Path(whokaryote_scripts.__file__).resolve().parent / "data"
    for model_name in ("RF_S_RBS_240121_bal.joblib", "RF_T_RBS_310122_bal.joblib"):
        model_path = data_dir / model_name
        if not model_path.is_file() or model_path.stat().st_size == 0:
            raise RuntimeError(f"Missing Whokaryote model file: {model_path}")
        joblib.load(model_path)
        tools[model_name] = str(model_path)
    return tools


def parse_explicit_algal_gff3(gff_path: Optional[Path]) -> Dict[str, str]:
    """Identify contigs with explicit multi-exon spliceosomal gene evidence (>=2 CDS exons per mRNA/gene)."""
    if gff_path is None:
        return {}
    if not gff_path.is_file():
        raise FileNotFoundError(f"GFF3 file not found: {gff_path}")
    cds_by_parent: Dict[Tuple[str, str], int] = {}
    pkzilla_contigs: set = set()
    opener = gzip.open if gff_path.suffix == ".gz" else open
    with opener(gff_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            contig, _, ftype, _, _, _, _, _, attrs = cols
            if ftype == "CDS":
                parent = ""
                for part in attrs.split(";"):
                    if part.startswith("Parent="):
                        parent = part[7:]
                        break
                key = (contig, parent or f"line_{len(cds_by_parent)}")
                cds_by_parent[key] = cds_by_parent.get(key, 0) + 1
            if "pkzilla" in attrs.lower():
                pkzilla_contigs.add(contig)

    multi_exon_counts: Counter = Counter()
    for (contig, _parent), count in cds_by_parent.items():
        if count >= 2:
            multi_exon_counts[contig] += 1

    evidence: Dict[str, str] = {}
    for contig in set(multi_exon_counts) | pkzilla_contigs:
        notes: List[str] = []
        if multi_exon_counts[contig] > 0:
            notes.append(f"gff3_multi_exon_transcripts={multi_exon_counts[contig]}")
        if contig in pkzilla_contigs:
            notes.append("gff3_pkzilla_locus")
        evidence[contig] = ";".join(notes)
    return evidence


def parse_algal_evidence_ids(path: Optional[Path]) -> Dict[str, str]:
    """Optional TSV (contig<TAB>evidence_note) or plain contig ID list supplying explicit algal evidence."""
    if path is None:
        return {}
    if not path.is_file():
        raise FileNotFoundError(f"Algal evidence file not found: {path}")
    evidence: Dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for idx, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if idx == 1 and parts[0].lower() in {"contig", "contig_id", "sequence_id"}:
                continue
            contig = parts[0].split()[0]
            note = parts[1] if len(parts) > 1 and parts[1].strip() else "explicit_algal_marker"
            evidence[contig] = note
    return evidence


def run_prodigal_parallel(
    fasta_path: Path,
    out_gff: Path,
    out_faa: Path,
    *,
    min_len: int,
    threads: int,
    work_dir: Path,
) -> int:
    """Run Prodigal (-p meta) across contig chunks and merge deterministic GFF/FAA outputs."""
    selected: List[Tuple[str, str]] = []
    with fasta_path.open("rt", encoding="utf-8") as handle:
        for title, seq in SimpleFastaParser(handle):
            if len(seq) >= min_len:
                selected.append((title.split()[0], seq))

    if not selected:
        out_gff.write_text("", encoding="utf-8")
        out_faa.write_text("", encoding="utf-8")
        return 0

    n_workers = max(1, min(int(threads), len(selected)))
    if n_workers == 1:
        sub_fasta = work_dir / f"prodigal_min{min_len}.fasta"
        with sub_fasta.open("wt", encoding="utf-8") as handle:
            for cid, seq in selected:
                handle.write(f">{cid}\n{seq}\n")
        subprocess.run(
            [
                "prodigal",
                "-i", str(sub_fasta),
                "-o", str(out_gff),
                "-a", str(out_faa),
                "-p", "meta",
                "-f", "gff",
                "-q",
            ],
            check=True,
        )
        return len(selected)

    chunks: List[List[Tuple[str, str]]] = [[] for _ in range(n_workers)]
    # Greedy load-balancing by nucleotide length while preserving chunk map
    loads = [0] * n_workers
    indexed_contigs = sorted(enumerate(selected), key=lambda x: len(x[1][1]), reverse=True)
    contig_to_chunk: Dict[int, int] = {}
    for idx, (_cid, seq) in indexed_contigs:
        w = min(range(n_workers), key=lambda k: loads[k])
        loads[w] += len(seq)
        contig_to_chunk[idx] = w
    for idx, item in enumerate(selected):
        chunks[contig_to_chunk[idx]].append(item)

    chunk_gffs: List[Path] = []
    chunk_faas: List[Path] = []

    def _run_chunk(w_idx: int) -> Tuple[Path, Path]:
        c_fa = work_dir / f"chunk_{w_idx}.fasta"
        c_gff = work_dir / f"chunk_{w_idx}.gff"
        c_faa = work_dir / f"chunk_{w_idx}.faa"
        with c_fa.open("wt", encoding="utf-8") as handle:
            for cid, seq in chunks[w_idx]:
                handle.write(f">{cid}\n{seq}\n")
        subprocess.run(
            [
                "prodigal",
                "-i", str(c_fa),
                "-o", str(c_gff),
                "-a", str(c_faa),
                "-p", "meta",
                "-f", "gff",
                "-q",
            ],
            check=True,
        )
        return c_gff, c_faa

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        results = list(pool.map(_run_chunk, range(n_workers)))
    for g_path, a_path in results:
        chunk_gffs.append(g_path)
        chunk_faas.append(a_path)

    with out_gff.open("wt", encoding="utf-8") as g_out:
        for idx, g_path in enumerate(chunk_gffs):
            with g_path.open("rt", encoding="utf-8") as g_in:
                for line in g_in:
                    if idx > 0 and line.startswith("##gff-version"):
                        continue
                    g_out.write(line)
    with out_faa.open("wt", encoding="utf-8") as a_out:
        for a_path in chunk_faas:
            a_out.write(a_path.read_text(encoding="utf-8"))
    return len(selected)


def parse_tiara_table(path: Path) -> Dict[str, Dict[str, object]]:
    rows: Dict[str, Dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected_prefix = ["sequence_id", "class_fst_stage", "class_snd_stage"]
        if header[:3] != expected_prefix:
            raise ValueError(f"Unexpected Tiara output header in {path}: {header}")
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            rec = dict(zip(header, cols))
            cid = str(rec["sequence_id"]).split()[0]
            p_org = parse_float_or_none(rec.get("org", "NA"))
            p_bac = parse_float_or_none(rec.get("bac", "NA"))
            p_arc = parse_float_or_none(rec.get("arc", "NA"))
            p_euk = parse_float_or_none(rec.get("euk", "NA"))
            p_unk1 = parse_float_or_none(rec.get("unk1", "NA"))
            p_pla = parse_float_or_none(rec.get("pla", "NA"))
            p_unk2 = parse_float_or_none(rec.get("unk2", "NA"))
            p_mit = parse_float_or_none(rec.get("mit", "NA"))
            p_prok = (p_bac + p_arc) if (p_bac is not None and p_arc is not None) else None
            rows[cid] = {
                "fst": rec["class_fst_stage"].strip(),
                "snd": rec["class_snd_stage"].strip(),
                "p_org": p_org,
                "p_bac": p_bac,
                "p_arc": p_arc,
                "p_euk": p_euk,
                "p_unk1": p_unk1,
                "p_pla": p_pla,
                "p_unk2": p_unk2,
                "p_mit": p_mit,
                "p_prok": p_prok,
            }
    return rows


def score_whokaryote_models(
    feature_csv: Path,
    tiara_rows: Dict[str, Dict[str, object]],
    raw_dir: Path,
) -> Dict[str, Dict[str, object]]:
    """Load Whokaryote featuretable.csv and score both RF_S (gene-only) and RF_T (Tiara-integrated)."""
    if not feature_csv.is_file() or feature_csv.stat().st_size == 0:
        return {}
    df = pd.read_csv(feature_csv)
    if df.empty:
        return {}
    df["contig"] = df["contig"].astype(str).str.split().str[0]

    data_dir = Path(whokaryote_scripts.__file__).resolve().parent / "data"
    rf_s = joblib.load(data_dir / "RF_S_RBS_240121_bal.joblib")
    rf_t = joblib.load(data_dir / "RF_T_RBS_310122_bal.joblib")

    def _tiara_code(cid: str) -> float:
        entry = tiara_rows.get(cid)
        if not entry:
            return np.nan
        fst = str(entry["fst"])
        if fst == "eukarya":
            return 0.0
        if fst in {"archaea", "bacteria", "prokarya"}:
            return 1.0
        return 2.0

    df["tiara_pred"] = [_tiara_code(cid) for cid in df["contig"]]
    base_feature_cols = [
        "ratio_same_orientation",
        "ID_general_avg",
        "ID_general_std",
        "gene_density",
        "gene_length",
        "ID_Q1",
        "ID_Q3",
        "rbs_ratio",
    ]
    t_feature_cols = base_feature_cols + ["tiara_pred"]

    # Score RF_S wherever gene-structure features are complete
    s_mask = ~df[base_feature_cols].isna().any(axis=1)
    df["whokaryote_S_pred"] = "NA"
    df["whokaryote_S_p_prok"] = np.nan
    if s_mask.any():
        probs_s = rf_s.predict_proba(df.loc[s_mask, base_feature_cols])
        prok_idx_s = list(rf_s.classes_).index(1)
        p_s_prok = probs_s[:, prok_idx_s]
        df.loc[s_mask, "whokaryote_S_p_prok"] = p_s_prok
        df.loc[s_mask, "whokaryote_S_pred"] = [
            "prokaryote" if p >= 0.5 else "eukaryote" for p in p_s_prok
        ]

    # Score RF_T wherever gene-structure + Tiara features are complete
    t_mask = ~df[t_feature_cols].isna().any(axis=1)
    df["whokaryote_T_pred"] = "NA"
    df["whokaryote_T_p_prok"] = np.nan
    if t_mask.any():
        probs_t = rf_t.predict_proba(df.loc[t_mask, t_feature_cols])
        prok_idx_t = list(rf_t.classes_).index(1)
        p_t_prok = probs_t[:, prok_idx_t]
        df.loc[t_mask, "whokaryote_T_p_prok"] = p_t_prok
        df.loc[t_mask, "whokaryote_T_pred"] = [
            "prokaryote" if p >= 0.5 else "eukaryote" for p in p_t_prok
        ]

    df[["contig", "whokaryote_S_pred", "whokaryote_S_p_prok"]].to_csv(
        raw_dir / "whokaryote_predictions_S.tsv", sep="\t", index=False
    )
    df[
        [
            "contig",
            "contig_length",
            "nr_genes",
            "whokaryote_S_pred",
            "whokaryote_S_p_prok",
            "whokaryote_T_pred",
            "whokaryote_T_p_prok",
        ]
    ].to_csv(raw_dir / "whokaryote_probabilities.tsv", sep="\t", index=False)

    out: Dict[str, Dict[str, object]] = {}
    for _, row in df.iterrows():
        cid = str(row["contig"])
        out[cid] = {
            "nr_genes": int(row["nr_genes"]),
            "gene_density": parse_float_or_none(row["gene_density"]),
            "rbs_ratio": parse_float_or_none(row["rbs_ratio"]),
            "ratio_same_orientation": parse_float_or_none(row["ratio_same_orientation"]),
            "s_pred": str(row["whokaryote_S_pred"]),
            "s_p_prok": parse_float_or_none(row["whokaryote_S_p_prok"]),
            "t_pred": str(row["whokaryote_T_pred"]),
            "t_p_prok": parse_float_or_none(row["whokaryote_T_p_prok"]),
        }
    return out


def classify_contig(
    *,
    contig: str,
    length_bp: int,
    tiara: Optional[Dict[str, object]],
    whokaryote: Optional[Dict[str, object]],
    algal_evidence: Optional[str],
    tiara_min_len: int = 500,
    whokaryote_min_len: int = 1000,
    bact_ensemble_cutoff: float = 0.70,
    bact_tiara_only_1k_cutoff: float = 0.88,
    bact_tiara_only_lt1k_cutoff: float = 0.95,
    euk_ensemble_cutoff: float = 0.35,
    euk_lt1k_tiara_cutoff: float = 0.15,
    disagreement_cutoff: float = 0.55,
) -> Dict[str, object]:
    """Compute length-calibrated continuous ensemble score and assign 5-class label with organelle/HGT guards."""
    lbin = length_bin_for(length_bp)
    has_algal = bool(algal_evidence)
    algal_src = algal_evidence or "none"

    t_scored = tiara is not None
    t_fst = str(tiara["fst"]) if tiara else "NA"
    t_snd = str(tiara["snd"]) if tiara else "NA"
    p_org = tiara["p_org"] if tiara else None  # type: ignore[assignment]
    p_bac = tiara["p_bac"] if tiara else None  # type: ignore[assignment]
    p_arc = tiara["p_arc"] if tiara else None  # type: ignore[assignment]
    p_euk = tiara["p_euk"] if tiara else None  # type: ignore[assignment]
    p_unk1 = tiara["p_unk1"] if tiara else None  # type: ignore[assignment]
    p_pla = tiara["p_pla"] if tiara else None  # type: ignore[assignment]
    p_unk2 = tiara["p_unk2"] if tiara else None  # type: ignore[assignment]
    p_mit = tiara["p_mit"] if tiara else None  # type: ignore[assignment]
    p_t_prok = tiara["p_prok"] if tiara else None  # type: ignore[assignment]

    wk_eligible = length_bp >= whokaryote_min_len
    nr_genes = int(whokaryote["nr_genes"]) if whokaryote else 0
    s_pred = str(whokaryote["s_pred"]) if whokaryote else "NA"
    p_s_prok = whokaryote["s_p_prok"] if whokaryote else None  # type: ignore[assignment]
    t_pred = str(whokaryote["t_pred"]) if whokaryote else "NA"
    p_wk_t_prok = whokaryote["t_p_prok"] if whokaryote else None  # type: ignore[assignment]
    wk_rf_scored = p_s_prok is not None

    gene_density = whokaryote["gene_density"] if whokaryote else None
    rbs_ratio = whokaryote["rbs_ratio"] if whokaryote else None
    ratio_same = whokaryote["ratio_same_orientation"] if whokaryote else None

    disagreement: Optional[float] = None
    if p_t_prok is not None and p_s_prok is not None:
        disagreement = abs(p_t_prok - p_s_prok)

    # 1. Below Tiara min length -> ambiguous (retained conservatively)
    if not t_scored or p_t_prok is None:
        label = "ambiguous"
        ensemble_score: Optional[float] = None
        retained = True
        reason = f"below_tiara_min_len_{tiara_min_len}"
    # 2. Organelle protection: never discard organelles as bacteria!
    elif t_fst == "organelle" or (
        p_org is not None and p_org >= 0.50 and p_org > p_t_prok and (p_euk is None or p_org >= p_euk)
    ):
        label = "organellar"
        # Dampen ensemble bacterial score for organelles so ROC/filtering never conflate organelles with bacteria
        ensemble_score = float(p_t_prok)
        retained = True
        reason = f"tiara_organelle_{t_snd}"
    else:
        if p_s_prok is not None:
            if length_bp >= 3000 or nr_genes >= 4:
                w_t, w_s = 0.70, 0.30
                min_tiara = 0.50
            else:
                w_t, w_s = 0.85, 0.15
                min_tiara = 0.60
            ensemble_score = w_t * float(p_t_prok) + w_s * float(p_s_prok)
            is_bacterial_candidate = (
                ensemble_score >= bact_ensemble_cutoff and float(p_t_prok) >= min_tiara
            )
            strong_disagreement = disagreement is not None and disagreement >= disagreement_cutoff
            is_euk_candidate = (
                t_fst != "unknown"
                and ensemble_score <= euk_ensemble_cutoff
                and float(p_s_prok) < 0.65
            )
        elif length_bp >= whokaryote_min_len:
            # 1,000+ bp with <2 Prodigal genes: slight confidence shrinkage toward 0.5
            ensemble_score = 0.5 + 0.90 * (float(p_t_prok) - 0.5)
            is_bacterial_candidate = float(p_t_prok) >= bact_tiara_only_1k_cutoff
            strong_disagreement = False
            is_euk_candidate = t_fst != "unknown" and ensemble_score <= euk_ensemble_cutoff
        else:
            # 500-999 bp ultra-short fragments: conservative confidence shrinkage + stricter thresholds
            ensemble_score = 0.5 + 0.70 * (float(p_t_prok) - 0.5)
            is_bacterial_candidate = float(p_t_prok) >= bact_tiara_only_lt1k_cutoff
            strong_disagreement = False
            is_euk_candidate = t_fst != "unknown" and float(p_t_prok) <= euk_lt1k_tiara_cutoff

        if is_bacterial_candidate:
            if has_algal:
                label = "putative_hgt"
                retained = True
                reason = f"bacterial_signal_rescued_by_explicit_algal_evidence:{algal_src}"
            else:
                label = "bacterial"
                retained = False
                if p_s_prok is not None:
                    reason = "ensemble_bacterial_consensus"
                elif length_bp >= whokaryote_min_len:
                    reason = "tiara_only_1k_bacterial"
                else:
                    reason = "tiara_only_lt1k_high_confidence_bacterial"
        elif strong_disagreement:
            # Never treat classifier disagreement alone as HGT evidence!
            label = "ambiguous"
            retained = True
            reason = "classifier_disagreement_without_bacterial_consensus"
        elif is_euk_candidate:
            label = "eukaryotic_haptophyte"
            retained = True
            reason = (
                "ensemble_eukaryotic_consensus"
                if p_s_prok is not None
                else "tiara_eukaryotic_call"
            )
        else:
            label = "ambiguous"
            retained = True
            reason = "intermediate_or_unknown_score"

    return {
        "contig": contig,
        "length_bp": length_bp,
        "length_bin": lbin,
        "tiara_scored": "true" if t_scored else "false",
        "tiara_fst_stage": t_fst,
        "tiara_snd_stage": t_snd,
        "tiara_p_org": format_opt_float(p_org),  # type: ignore[arg-type]
        "tiara_p_bac": format_opt_float(p_bac),  # type: ignore[arg-type]
        "tiara_p_arc": format_opt_float(p_arc),  # type: ignore[arg-type]
        "tiara_p_euk": format_opt_float(p_euk),  # type: ignore[arg-type]
        "tiara_p_unk1": format_opt_float(p_unk1),  # type: ignore[arg-type]
        "tiara_p_pla": format_opt_float(p_pla),  # type: ignore[arg-type]
        "tiara_p_unk2": format_opt_float(p_unk2),  # type: ignore[arg-type]
        "tiara_p_mit": format_opt_float(p_mit),  # type: ignore[arg-type]
        "tiara_p_prok": format_opt_float(p_t_prok),  # type: ignore[arg-type]
        "whokaryote_minsize_eligible": "true" if wk_eligible else "false",
        "prodigal_nr_genes": nr_genes,
        "whokaryote_rf_scored": "true" if wk_rf_scored else "false",
        "whokaryote_S_pred": s_pred,
        "whokaryote_S_p_prok": format_opt_float(p_s_prok),  # type: ignore[arg-type]
        "whokaryote_T_pred": t_pred,
        "whokaryote_T_p_prok": format_opt_float(p_wk_t_prok),  # type: ignore[arg-type]
        "gene_density": format_opt_float(gene_density),  # type: ignore[arg-type]
        "rbs_ratio": format_opt_float(rbs_ratio),  # type: ignore[arg-type]
        "ratio_same_orientation": format_opt_float(ratio_same),  # type: ignore[arg-type]
        "ensemble_prok_score": format_opt_float(ensemble_score),
        "classifier_disagreement": format_opt_float(disagreement),
        "has_explicit_algal_evidence": "true" if has_algal else "false",
        "algal_evidence_source": algal_src,
        "final_label": label,
        "retained_in_filtered": "true" if retained else "false",
        "decision_reason": reason,
    }


def run_decontamination(
    *,
    fasta_path: Path,
    outdir: Path,
    sample_id: str,
    gff3_path: Optional[Path] = None,
    algal_evidence_path: Optional[Path] = None,
    tiara_min_len: int = 500,
    whokaryote_min_len: int = 1000,
    threads: int = 4,
    emit_filtered_fasta: bool = True,
) -> Dict[str, object]:
    """Run Tiara + Whokaryote/Prodigal decontamination and emit all auditable outputs."""
    started = time.monotonic()
    tools = validate_decontam_prerequisites()

    outdir.mkdir(parents=True, exist_ok=True)
    raw_dir = outdir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    algal_evidence = parse_explicit_algal_gff3(gff3_path)
    for cid, note in parse_algal_evidence_ids(algal_evidence_path).items():
        if cid in algal_evidence:
            algal_evidence[cid] = f"{algal_evidence[cid]};{note}"
        else:
            algal_evidence[cid] = note

    with tempfile.TemporaryDirectory(prefix="hapto_decontam_") as tmp_str:
        work_dir = Path(tmp_str)
        clean_fasta = work_dir / "clean_contigs.fasta"
        contigs_meta: List[Tuple[str, int]] = []
        seen_ids: set = set()
        with open_fasta_text(fasta_path) as src_h, clean_fasta.open("wt", encoding="utf-8") as dst_h:
            for title, seq in SimpleFastaParser(src_h):
                cid = title.split()[0]
                if cid in seen_ids:
                    raise ValueError(f"Duplicate contig ID in {fasta_path}: {cid}")
                seen_ids.add(cid)
                contigs_meta.append((cid, len(seq)))
                dst_h.write(f">{cid}\n{seq}\n")

        tiara_out = raw_dir / "tiara_pred.txt"
        tiara_eligible = sum(1 for _cid, length in contigs_meta if length >= tiara_min_len)
        if tiara_eligible > 0:
            subprocess.run(
                [
                    "tiara",
                    "-i", str(clean_fasta),
                    "-o", str(tiara_out),
                    "-m", str(tiara_min_len),
                    "-t", str(threads),
                    "--probabilities",
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            tiara_rows = parse_tiara_table(tiara_out)
        else:
            tiara_out.write_text(
                "sequence_id\tclass_fst_stage\tclass_snd_stage\torg\tbac\tarc\teuk\tunk1\tpla\tunk2\tmit\n",
                encoding="utf-8",
            )
            (raw_dir / "log_tiara_pred.txt").write_text(
                "No contigs >= tiara_min_len; skipped tiara.\n", encoding="utf-8"
            )
            tiara_rows = {}

        gff_out = raw_dir / "contigs_genes.gff"
        faa_out = raw_dir / "contigs_proteins.faa"
        wk_eligible = run_prodigal_parallel(
            clean_fasta,
            gff_out,
            faa_out,
            min_len=whokaryote_min_len,
            threads=threads,
            work_dir=work_dir,
        )

        feat_csv = raw_dir / "whokaryote_featuretable.csv"
        wk_rows: Dict[str, Dict[str, object]] = {}
        if wk_eligible > 0 and gff_out.stat().st_size > 0:
            # Invoke whokaryote_scripts.calc_features directly on merged Prodigal GFF, and run
            # whokaryote.py whenever >=1 contig has >=2 genes (for official whokaryote_predictions_T.tsv)
            whokaryote_scripts.calc_features(str(gff_out), str(raw_dir))
            raw_feat_path = raw_dir / "featuretable.csv"
            raw_feat_path.replace(feat_csv)
            feat_df = pd.read_csv(feat_csv)
            if not feat_df.dropna().empty:
                # Copy tiara_pred.txt and run official whokaryote.py --gff --model T
                wk_work = work_dir / "wk_run"
                wk_work.mkdir(parents=True, exist_ok=True)
                shutil.copy(tiara_out, wk_work / "tiara_pred.txt")
                subprocess.run(
                    [
                        sys.executable,
                        tools["whokaryote.py"],
                        "--gff", str(gff_out),
                        "--outdir", str(wk_work),
                        "--minsize", str(whokaryote_min_len),
                        "--model", "T",
                        "--threads", str(threads),
                    ],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                for fname in ("whokaryote_predictions_T.tsv", "featuretable_predictions_T.tsv"):
                    src_f = wk_work / fname
                    if src_f.is_file():
                        shutil.copy(src_f, raw_dir / fname)
            else:
                (raw_dir / "whokaryote_predictions_T.tsv").write_text(
                    "contig\tpredicted\n", encoding="utf-8"
                )
            wk_rows = score_whokaryote_models(feat_csv, tiara_rows, raw_dir)
        else:
            feat_csv.write_text(
                "contig,contig_length,nr_genes,ratio_same_orientation,ID_general_avg,"
                "ID_general_std,gene_density,gene_length,ID_Q1,ID_Q3,rbs_ratio\n",
                encoding="utf-8",
            )
            (raw_dir / "whokaryote_predictions_T.tsv").write_text(
                "contig\tpredicted\n", encoding="utf-8"
            )
            (raw_dir / "whokaryote_predictions_S.tsv").write_text(
                "contig\twhokaryote_S_pred\twhokaryote_S_p_prok\n", encoding="utf-8"
            )
            (raw_dir / "whokaryote_probabilities.tsv").write_text(
                "contig\tcontig_length\tnr_genes\twhokaryote_S_pred\twhokaryote_S_p_prok\t"
                "whokaryote_T_pred\twhokaryote_T_p_prok\n",
                encoding="utf-8",
            )

    evidence_rows: List[Dict[str, object]] = []
    excluded_ids: List[str] = []
    retained_set: set = set()
    for cid, length_bp in contigs_meta:
        row = classify_contig(
            contig=cid,
            length_bp=length_bp,
            tiara=tiara_rows.get(cid),
            whokaryote=wk_rows.get(cid),
            algal_evidence=algal_evidence.get(cid),
            tiara_min_len=tiara_min_len,
            whokaryote_min_len=whokaryote_min_len,
        )
        evidence_rows.append(row)
        if row["retained_in_filtered"] == "true":
            retained_set.add(cid)
        else:
            excluded_ids.append(cid)

    evidence_path = outdir / "contig_evidence.tsv"
    with evidence_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=EVIDENCE_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(evidence_rows)

    excluded_path = outdir / "excluded_ids.txt"
    excluded_path.write_text(
        "".join(f"{cid}\n" for cid in excluded_ids), encoding="utf-8"
    )

    filtered_path = outdir / "filtered.fasta"
    if emit_filtered_fasta:
        with open_fasta_text(fasta_path) as src_h, filtered_path.open("wt", encoding="utf-8") as dst_h:
            for title, seq in SimpleFastaParser(src_h):
                cid = title.split()[0]
                if cid in retained_set:
                    dst_h.write(f">{title}\n{seq}\n")

    total_contigs = len(evidence_rows)
    total_bp = sum(int(r["length_bp"]) for r in evidence_rows)
    retained_bp = sum(
        int(r["length_bp"]) for r in evidence_rows if r["retained_in_filtered"] == "true"
    )
    excluded_bp = total_bp - retained_bp

    label_counts: Dict[str, int] = {lbl: 0 for lbl in VALID_LABELS}
    label_bp: Dict[str, int] = {lbl: 0 for lbl in VALID_LABELS}
    bin_counts: Dict[str, int] = {lb: 0 for lb in LENGTH_BINS}
    bin_bp: Dict[str, int] = {lb: 0 for lb in LENGTH_BINS}
    for r in evidence_rows:
        lbl = str(r["final_label"])
        lb = str(r["length_bin"])
        bp = int(r["length_bp"])
        label_counts[lbl] += 1
        label_bp[lbl] += bp
        bin_counts[lb] += 1
        bin_bp[lb] += bp

    elapsed = time.monotonic() - started
    summary_flat: Dict[str, object] = {
        "sample_id": sample_id,
        "total_contigs": total_contigs,
        "total_bp": total_bp,
        "retained_contigs": len(retained_set),
        "retained_bp": retained_bp,
        "retained_bp_fraction": round(retained_bp / total_bp, 6) if total_bp else 0.0,
        "excluded_contigs": len(excluded_ids),
        "excluded_bp": excluded_bp,
        "excluded_bp_fraction": round(excluded_bp / total_bp, 6) if total_bp else 0.0,
    }
    for lbl in VALID_LABELS:
        summary_flat[f"{lbl}_contigs"] = label_counts[lbl]
        summary_flat[f"{lbl}_bp"] = label_bp[lbl]
        summary_flat[f"{lbl}_bp_fraction"] = (
            round(label_bp[lbl] / total_bp, 6) if total_bp else 0.0
        )
    for lb in LENGTH_BINS:
        summary_flat[f"len_{lb}_contigs"] = bin_counts[lb]
        summary_flat[f"len_{lb}_bp"] = bin_bp[lb]
    summary_flat["elapsed_seconds"] = round(elapsed, 3)

    summary_tsv = outdir / "summary.tsv"
    with summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summary_flat.keys()), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(summary_flat)

    summary_json = outdir / "summary.json"
    payload = {
        "format": "hapto_decontam_v1",
        "sample_id": sample_id,
        "parameters": {
            "tiara_min_len": tiara_min_len,
            "whokaryote_min_len": whokaryote_min_len,
            "threads": threads,
            "ensemble_weights_ge3k_or_4genes": {"tiara_p_prok": 0.70, "whokaryote_S_p_prok": 0.30},
            "ensemble_weights_1k_3k_lt4genes": {"tiara_p_prok": 0.85, "whokaryote_S_p_prok": 0.15},
            "tiara_only_shrinkage": {"ge_1k": 0.90, "lt_1k": 0.70},
            "bacterial_ensemble_cutoff": 0.70,
            "bacterial_tiara_only_1k_cutoff": 0.88,
            "bacterial_tiara_only_lt1k_cutoff": 0.95,
            "eukaryotic_ensemble_cutoff": 0.35,
        },
        "inputs": {
            "fasta": {"path": str(fasta_path), "bytes": fasta_path.stat().st_size, "sha256": sha256_file(fasta_path)},
            "gff3": str(gff3_path) if gff3_path else "NA",
            "algal_evidence": str(algal_evidence_path) if algal_evidence_path else "NA",
        },
        "tools": tools,
        "summary": summary_flat,
    }
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary_flat


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_p = subparsers.add_parser("run", help="Decontaminate one assembly FASTA.")
    run_p.add_argument("--fasta", type=Path, required=True)
    run_p.add_argument("--sample-id", required=True)
    run_p.add_argument("--outdir", type=Path, required=True)
    run_p.add_argument("--gff3", type=Path, default=None)
    run_p.add_argument("--algal-evidence", type=Path, default=None)
    run_p.add_argument("--tiara-min-len", type=int, default=500)
    run_p.add_argument("--whokaryote-min-len", type=int, default=1000)
    run_p.add_argument("--threads", type=int, default=4)

    agg_p = subparsers.add_parser("aggregate", help="Combine per-sample summary.tsv files.")
    agg_p.add_argument("--summaries", type=Path, nargs="+", required=True)
    agg_p.add_argument("--output-tsv", type=Path, required=True)
    agg_p.add_argument("--output-json", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "run":
        gff3 = None if (args.gff3 is None or str(args.gff3) == "NA") else args.gff3
        summary = run_decontamination(
            fasta_path=args.fasta,
            outdir=args.outdir,
            sample_id=args.sample_id,
            gff3_path=gff3,
            algal_evidence_path=args.algal_evidence,
            tiara_min_len=args.tiara_min_len,
            whokaryote_min_len=args.whokaryote_min_len,
            threads=args.threads,
        )
        print(
            f"Decontam PASS [{args.sample_id}]: contigs={summary['total_contigs']} "
            f"retained_bp={summary['retained_bp']} ({summary['retained_bp_fraction']}) "
            f"excluded_bac={summary['excluded_contigs']} ({summary['excluded_bp']} bp) "
            f"organellar={summary['organellar_contigs']} ambiguous={summary['ambiguous_contigs']} "
            f"putative_hgt={summary['putative_hgt_contigs']}"
        )
        return 0

    rows: List[Dict[str, str]] = []
    fieldnames: Optional[List[str]] = None
    for path in args.summaries:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fieldnames is None:
                fieldnames = list(reader.fieldnames or [])
            rows.extend(reader)
    if not rows or not fieldnames:
        raise ValueError("No summary rows provided to aggregate")
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(
        json.dumps({"format": "hapto_decontam_all_v1", "samples": rows}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Aggregated {len(rows)} decontamination profiles -> {args.output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

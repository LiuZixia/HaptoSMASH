#!/usr/bin/env python3
"""Evaluate Track A gigasynthase predictions against 12B1 PKZILLA truth and controls."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pyhmmer
from Bio import SeqIO

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import detect_gigasynthase_hotspots as det  # noqa: E402


ORIG_COORD_RE = re.compile(
    r"^(?:(?P<label>.+?)__)?(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+):[^_]*_(?P<strand>[+-])$"
)
MERGED_COORD_RE = re.compile(
    r"^(?P<contig>[^:]+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)$"
)
MANIFEST_LOCUS_RE = re.compile(
    r"^[^;]+;(?P<contig>[^:]+):(?P<start1>\d+)-(?P<end1>\d+)\((?P<strand>[+-])\)$"
)

CODON_FOR_AA = {
    "A": "GCT", "C": "TGC", "D": "GAC", "E": "GAA", "F": "TTC",
    "G": "GGT", "H": "CAC", "I": "ATC", "K": "AAG", "L": "CTG",
    "M": "ATG", "N": "AAC", "P": "CCT", "Q": "CAG", "R": "CGT",
    "S": "TCC", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAC",
    "X": "NNN",
}


@dataclass(frozen=True)
class GenomicInterval:
    label: str
    contig: str
    start0: int
    end0: int
    strand: str

    @property
    def length_bp(self) -> int:
        return self.end0 - self.start0


@dataclass(frozen=True)
class TruthTsvRow:
    hotspot_coord_raw: str
    size: int
    original_label: Optional[str]
    original_contig: str
    original_start_raw: int
    original_end_raw: int
    original_strand: str
    merged_contig: str
    merged_start0: int
    merged_end0: int
    merged_strand: str


def parse_truth_bed(path: Path) -> list[GenomicInterval]:
    """Parse authoritative 0-based half-open genomic BED6 truth intervals."""
    intervals: list[GenomicInterval] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue
            cols = stripped.split("\t")
            if len(cols) < 6:
                raise ValueError(f"{path}:{line_number}: expected BED6 columns, got {len(cols)}")
            contig = cols[0]
            start0, end0 = int(cols[1]), int(cols[2])
            name, strand = cols[3], cols[5]
            if not (0 <= start0 < end0):
                raise ValueError(f"{path}:{line_number}: invalid genomic BED interval {start0}-{end0}")
            if strand not in {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid strand {strand!r}")
            intervals.append(
                GenomicInterval(label=name, contig=contig, start0=start0, end0=end0, strand=strand)
            )
    if not intervals:
        raise ValueError(f"No truth intervals loaded from {path}")
    return intervals


def parse_truth_tsv(path: Path) -> list[TruthTsvRow]:
    """Parse truth TSV keeping original v0.3 coordinates and lifted v1.1 merged coordinates strictly separate."""
    rows: list[TruthTsvRow] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["hotspot_coord", "size"]:
            raise ValueError(f"Unexpected truth TSV header in {path}: {reader.fieldnames}")
        for line_number, row in enumerate(reader, 2):
            raw = row["hotspot_coord"]
            size = int(row["size"])
            parts = raw.split("::")
            if len(parts) != 2:
                raise ValueError(f"{path}:{line_number}: expected original::merged format, got {raw!r}")
            orig_match = ORIG_COORD_RE.fullmatch(parts[0])
            merged_match = MERGED_COORD_RE.fullmatch(parts[1])
            if not orig_match or not merged_match:
                raise ValueError(f"{path}:{line_number}: could not parse coordinate parts in {raw!r}")
            m_start0 = int(merged_match.group("start"))
            m_end0 = int(merged_match.group("end"))
            if m_end0 - m_start0 != size:
                raise ValueError(
                    f"{path}:{line_number}: merged span {m_end0 - m_start0} != size {size}"
                )
            rows.append(
                TruthTsvRow(
                    hotspot_coord_raw=raw,
                    size=size,
                    original_label=orig_match.group("label"),
                    original_contig=orig_match.group("contig"),
                    original_start_raw=int(orig_match.group("start")),
                    original_end_raw=int(orig_match.group("end")),
                    original_strand=orig_match.group("strand"),
                    merged_contig=merged_match.group("contig"),
                    merged_start0=m_start0,
                    merged_end0=m_end0,
                    merged_strand=merged_match.group("strand"),
                )
            )
    return rows


def verify_bed_and_tsv_consistency(
    bed_intervals: list[GenomicInterval], tsv_rows: list[TruthTsvRow]
) -> None:
    if len(bed_intervals) != len(tsv_rows):
        raise ValueError(
            f"Truth BED count ({len(bed_intervals)}) != truth TSV count ({len(tsv_rows)})"
        )
    bed_set = {(b.contig, b.start0, b.end0, b.strand) for b in bed_intervals}
    tsv_set = {
        (t.merged_contig, t.merged_start0, t.merged_end0, t.merged_strand) for t in tsv_rows
    }
    if bed_set != tsv_set:
        raise ValueError("Merged coordinates in truth TSV do not match authoritative BED6 intervals")


def parse_unmerged_pkzilla_loci(manifest_path: Path) -> list[GenomicInterval]:
    """Load 12B1 PKZILLA-1/-2/-3 unmerged genomic spans (converted from 1-based inclusive to 0-based half-open)."""
    loci: list[GenomicInterval] = []
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sid = row["sequence_id"]
            if sid not in {"PKZILLA-1", "PKZILLA-2", "PKZILLA-3"}:
                continue
            match = MANIFEST_LOCUS_RE.fullmatch(row["source_identifier_or_coordinates"])
            if not match:
                raise ValueError(f"Could not parse genomic locus for {sid} in {manifest_path}")
            start1 = int(match.group("start1"))
            end1 = int(match.group("end1"))
            if not (1 <= start1 <= end1):
                raise ValueError(f"Invalid 1-based coordinates for {sid}: {start1}-{end1}")
            loci.append(
                GenomicInterval(
                    label=sid,
                    contig=match.group("contig"),
                    start0=start1 - 1,
                    end0=end1,
                    strand=match.group("strand"),
                )
            )
    if len(loci) != 3:
        raise ValueError(f"Expected 3 12B1 PKZILLA loci in {manifest_path}, found {len(loci)}")
    return loci


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted((s, e) for s, e in intervals if e > s)
    if not ordered:
        return []
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_s, prev_e = merged[-1]
        if start <= prev_e:
            merged[-1] = (prev_s, max(prev_e, end))
        else:
            merged.append((start, end))
    return merged


def union_bp(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(e - s for s, e in merge_intervals(intervals))


def intersection_bp(
    a_intervals: Iterable[tuple[int, int]],
    b_intervals: Iterable[tuple[int, int]],
) -> int:
    a = merge_intervals(a_intervals)
    b = merge_intervals(b_intervals)
    i = j = total = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if hi > lo:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def read_hotspots_tsv(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                {
                    **row,
                    "g_start": int(row["g_start"]),
                    "g_end": int(row["g_end"]),
                    "length_bp": int(row["length_bp"]),
                    "n_domains": int(row["n_domains"]),
                    "n_ks_ca_backbone": int(row["n_ks_ca_backbone"]),
                    "n_neg_control_overlaps": int(row["n_neg_control_overlaps"]),
                }
            )
    return rows


def evaluate_predictions(
    hotspots: list[dict[str, object]],
    unmerged_loci: list[GenomicInterval],
    merged_bed: list[GenomicInterval],
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    pred_by_contig: dict[str, list[dict[str, object]]] = {}
    pred_by_cs: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in hotspots:
        contig = str(row["contig"])
        strand = str(row["strand"])
        pred_by_contig.setdefault(contig, []).append(row)
        pred_by_cs.setdefault((contig, strand), []).append(row)

    merged_by_prefix: dict[str, GenomicInterval] = {}
    for item in merged_bed:
        if "__" in item.label:
            merged_by_prefix[item.label.split("__", 1)[0]] = item

    locus_rows: list[dict[str, object]] = []
    for locus in unmerged_loci:
        same_comps = [
            r
            for r in pred_by_cs.get((locus.contig, locus.strand), [])
            if int(r["g_end"]) > locus.start0 and int(r["g_start"]) < locus.end0
        ]
        any_comps = [
            r
            for r in pred_by_contig.get(locus.contig, [])
            if int(r["g_end"]) > locus.start0 and int(r["g_start"]) < locus.end0
        ]
        same_intervals = [(int(r["g_start"]), int(r["g_end"])) for r in same_comps]
        any_intervals = [(int(r["g_start"]), int(r["g_end"])) for r in any_comps]
        same_cov_bp = intersection_bp(same_intervals, [(locus.start0, locus.end0)])
        any_cov_bp = intersection_bp(any_intervals, [(locus.start0, locus.end0)])

        merged_win = merged_by_prefix.get(locus.label)
        merged_same_bp = (
            intersection_bp(same_intervals, [(merged_win.start0, merged_win.end0)])
            if merged_win
            else 0
        )
        merged_any_bp = (
            intersection_bp(any_intervals, [(merged_win.start0, merged_win.end0)])
            if merged_win
            else 0
        )
        locus_rows.append(
            {
                "locus_id": locus.label,
                "contig": locus.contig,
                "unmerged_start0": locus.start0,
                "unmerged_end0": locus.end0,
                "unmerged_length_bp": locus.length_bp,
                "strand": locus.strand,
                "detected_same_strand": len(same_comps) > 0,
                "detected_any_strand": len(any_comps) > 0,
                "components_same_strand": len(same_comps),
                "components_any_strand": len(any_comps),
                "unmerged_covered_bp_same_strand": same_cov_bp,
                "unmerged_coverage_same_strand": round(same_cov_bp / locus.length_bp, 6),
                "unmerged_covered_bp_any_strand": any_cov_bp,
                "unmerged_coverage_any_strand": round(any_cov_bp / locus.length_bp, 6),
                "merged_window_bp": merged_win.length_bp if merged_win else 0,
                "merged_overlap_bp_same_strand": merged_same_bp,
                "merged_coverage_same_strand": (
                    round(merged_same_bp / merged_win.length_bp, 6) if merged_win else 0.0
                ),
                "merged_overlap_bp_any_strand": merged_any_bp,
                "merged_coverage_any_strand": (
                    round(merged_any_bp / merged_win.length_bp, 6) if merged_win else 0.0
                ),
                "component_ids_same_strand": ";".join(str(r["hotspot_id"]) for r in same_comps) or "none",
                "component_domains_same_strand": sum(int(r["n_domains"]) for r in same_comps),
                "edge_truncated_components": sum(
                    1 for r in same_comps if r["edge_truncated"] == "true"
                ),
            }
        )

    merged_rows: list[dict[str, object]] = []
    truth_intervals_by_contig: dict[str, list[tuple[int, int]]] = {}
    truth_intervals_by_cs: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for truth in merged_bed:
        truth_intervals_by_contig.setdefault(truth.contig, []).append((truth.start0, truth.end0))
        truth_intervals_by_cs.setdefault((truth.contig, truth.strand), []).append(
            (truth.start0, truth.end0)
        )
        same_comps = [
            r
            for r in pred_by_cs.get((truth.contig, truth.strand), [])
            if int(r["g_end"]) > truth.start0 and int(r["g_start"]) < truth.end0
        ]
        any_comps = [
            r
            for r in pred_by_contig.get(truth.contig, [])
            if int(r["g_end"]) > truth.start0 and int(r["g_start"]) < truth.end0
        ]
        same_bp = intersection_bp(
            [(int(r["g_start"]), int(r["g_end"])) for r in same_comps],
            [(truth.start0, truth.end0)],
        )
        any_bp = intersection_bp(
            [(int(r["g_start"]), int(r["g_end"])) for r in any_comps],
            [(truth.start0, truth.end0)],
        )
        merged_rows.append(
            {
                "truth_name": truth.label,
                "contig": truth.contig,
                "start0": truth.start0,
                "end0": truth.end0,
                "strand": truth.strand,
                "length_bp": truth.length_bp,
                "detected_same_strand": len(same_comps) > 0,
                "detected_any_strand": len(any_comps) > 0,
                "components_same_strand": len(same_comps),
                "overlap_bp_same_strand": same_bp,
                "coverage_same_strand": round(same_bp / truth.length_bp, 6),
                "overlap_bp_any_strand": any_bp,
                "coverage_any_strand": round(any_bp / truth.length_bp, 6),
            }
        )

    total_pred_bp_any = sum(
        union_bp((int(r["g_start"]), int(r["g_end"])) for r in rows)
        for rows in pred_by_contig.values()
    )
    inside_truth_bp_any = sum(
        intersection_bp(
            [(int(r["g_start"]), int(r["g_end"])) for r in rows],
            truth_intervals_by_contig.get(contig, []),
        )
        for contig, rows in pred_by_contig.items()
    )
    outside_truth_bp_any = total_pred_bp_any - inside_truth_bp_any

    total_pred_bp_same = sum(
        union_bp((int(r["g_start"]), int(r["g_end"])) for r in rows)
        for rows in pred_by_cs.values()
    )
    inside_truth_bp_same = sum(
        intersection_bp(
            [(int(r["g_start"]), int(r["g_end"])) for r in rows],
            truth_intervals_by_cs.get(cs, []),
        )
        for cs, rows in pred_by_cs.items()
    )
    outside_truth_bp_same = total_pred_bp_same - inside_truth_bp_same

    novel_hotspots_any = 0
    novel_hotspots_same = 0
    for row in hotspots:
        contig = str(row["contig"])
        strand = str(row["strand"])
        iv = [(int(row["g_start"]), int(row["g_end"]))]
        if intersection_bp(iv, truth_intervals_by_contig.get(contig, [])) == 0:
            novel_hotspots_any += 1
        if intersection_bp(iv, truth_intervals_by_cs.get((contig, strand), [])) == 0:
            novel_hotspots_same += 1

    truth_total_bp = sum(
        union_bp(ivs) for ivs in truth_intervals_by_contig.values()
    )
    summary = {
        "pkzilla_unmerged_loci_total": len(locus_rows),
        "pkzilla_recall_same_strand": sum(1 for r in locus_rows if r["detected_same_strand"]),
        "pkzilla_recall_any_strand": sum(1 for r in locus_rows if r["detected_any_strand"]),
        "pkzilla_min_unmerged_coverage_same_strand": min(
            float(r["unmerged_coverage_same_strand"]) for r in locus_rows
        ),
        "pkzilla_mean_unmerged_coverage_same_strand": round(
            sum(float(r["unmerged_coverage_same_strand"]) for r in locus_rows) / len(locus_rows),
            6,
        ),
        "pkzilla_components_per_locus": {
            str(r["locus_id"]): int(r["components_same_strand"]) for r in locus_rows
        },
        "merged_truth_regions_total": len(merged_rows),
        "merged_truth_recall_same_strand": sum(
            1 for r in merged_rows if r["detected_same_strand"]
        ),
        "merged_truth_recall_any_strand": sum(
            1 for r in merged_rows if r["detected_any_strand"]
        ),
        "merged_truth_total_bp": truth_total_bp,
        "predicted_hotspots_total": len(hotspots),
        "predicted_union_bp_any_strand": total_pred_bp_any,
        "predicted_inside_truth_bp_any_strand": inside_truth_bp_any,
        "predicted_outside_truth_bp_any_strand": outside_truth_bp_any,
        "predicted_outside_truth_fraction_any_strand": (
            round(outside_truth_bp_any / total_pred_bp_any, 6) if total_pred_bp_any else 0.0
        ),
        "predicted_union_bp_same_strand": total_pred_bp_same,
        "predicted_inside_truth_bp_same_strand": inside_truth_bp_same,
        "predicted_outside_truth_bp_same_strand": outside_truth_bp_same,
        "predicted_outside_truth_fraction_same_strand": (
            round(outside_truth_bp_same / total_pred_bp_same, 6) if total_pred_bp_same else 0.0
        ),
        "predicted_hotspots_with_zero_truth_overlap_any_strand": novel_hotspots_any,
        "predicted_hotspots_with_zero_truth_overlap_same_strand": novel_hotspots_same,
        "predicted_12b1_hotspots_with_neg_control_domains": sum(
            1 for r in hotspots if int(r["n_neg_control_overlaps"]) > 0
        ),
    }
    return summary, locus_rows, merged_rows


def evaluate_pufa_fas_controls(
    pufa_faa: Path,
    hmm_path: Path,
    pfam_table: dict[str, det.PfamModelInfo],
    *,
    ehuxleyi_fasta: Optional[Path] = None,
    d_max_bp: int = 15000,
    min_domains: int = 5,
    min_orf_aa: int = 150,
    cpus: int = 2,
) -> list[dict[str, object]]:
    """Evaluate negative-control PUFA/FAS seeds (and real E. huxleyi FAS/KAS III scaffolds if provided)."""
    hmms, ga_cutoffs = det.load_hmms(hmm_path, pfam_table)
    pad = "TAGTGA" * 20
    results: list[dict[str, object]] = []

    for rec in SeqIO.parse(str(pufa_faa), "fasta"):
        aa = str(rec.seq).upper().rstrip("*")
        dna = pad + "".join(CODON_FOR_AA[c] for c in aa) + pad
        orfs = list(det.extract_six_frame_orfs(rec.id, dna, min_orf_aa=min_orf_aa))
        hits = det.scan_orf_batch(orfs, hmms, pfam_table, ga_cutoffs, cpus=cpus) if orfs else []
        hotspots = det.cluster_track_a_hits(hits, d_max_bp=d_max_bp, min_domains=min_domains)
        results.append(
            {
                "control_id": rec.id,
                "source_type": "seed_backtranslated_dna",
                "length_bp": len(dna),
                "domain_hits_total": len(hits),
                "track_a_cluster_hits": sum(
                    1 for h in hits if h.track == "A" and h.role in {"backbone", "ancillary"}
                ),
                "neg_control_hmm_hits": sum(
                    1 for h in hits if h.role == "negative_control"
                ),
                "false_hotspot_calls": len(hotspots),
                "models_hit": ";".join(sorted({h.model_name for h in hits})) or "none",
            }
        )

    if ehuxleyi_fasta is not None and ehuxleyi_fasta.is_file():
        target_scaffolds = {"NW_005201381.1", "NW_005202303.1"}
        with det.open_text_maybe_gzip(ehuxleyi_fasta) as handle:
            for contig_idx, rec in enumerate(SeqIO.parse(handle, "fasta")):
                if rec.id not in target_scaffolds:
                    continue
                orfs = list(
                    det.extract_six_frame_orfs(
                        rec.id, str(rec.seq), contig_index=contig_idx, min_orf_aa=min_orf_aa
                    )
                )
                hits = (
                    det.scan_orf_batch(orfs, hmms, pfam_table, ga_cutoffs, cpus=cpus)
                    if orfs
                    else []
                )
                hotspots = det.cluster_track_a_hits(
                    hits, d_max_bp=d_max_bp, min_domains=min_domains
                )
                results.append(
                    {
                        "control_id": rec.id,
                        "source_type": "ehuxleyi_CCMP1516_genomic_scaffold",
                        "length_bp": len(rec.seq),
                        "domain_hits_total": len(hits),
                        "track_a_cluster_hits": sum(
                            1 for h in hits if h.track == "A" and h.role in {"backbone", "ancillary"}
                        ),
                        "neg_control_hmm_hits": sum(
                            1 for h in hits if h.role == "negative_control"
                        ),
                        "false_hotspot_calls": len(hotspots),
                        "models_hit": ";".join(sorted({h.model_name for h in hits})) or "none",
                    }
                )
    return results


def evaluate_pkzilla_protein_domains(
    pkzilla_faa: Path,
    domains_dir: Path,
    hmm_path: Path,
    pfam_table: dict[str, det.PfamModelInfo],
    *,
    cpus: int = 2,
) -> list[dict[str, object]]:
    """Protein-relative domain check for PKZILLA-1/-2/-3/-B1 (never mixed with genomic intervals)."""
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(pkzilla_faa, digital=True, alphabet=alphabet) as sf:
        seqs = list(sf)
    hmms, _ = det.load_hmms(hmm_path, pfam_table)

    hits_by_prot: dict[str, list[tuple[int, int, str, str]]] = {}
    for top_hits in pyhmmer.hmmscan(seqs, hmms, bit_cutoffs="gathering", cpus=cpus):
        prot_id = det.as_str(top_hits.query.name)
        prot_hits: list[tuple[int, int, str, str]] = []
        for hit in top_hits:
            if not hit.included:
                continue
            acc = det.as_str(hit.accession)
            info = pfam_table[acc]
            if info.track != "A" or info.role not in {"backbone", "ancillary"}:
                continue
            for dom in hit.domains.included:
                aln = dom.alignment
                prot_hits.append((int(aln.target_from) - 1, int(aln.target_to), acc, info.model_name))
        hits_by_prot[prot_id] = prot_hits

    rows: list[dict[str, object]] = []
    for prot_id in ("PKZILLA-1", "PKZILLA-2", "PKZILLA-3", "PKZILLA-B1"):
        bed_path = domains_dir / f"{prot_id}_PKS_domains.bed"
        truth_domains: list[tuple[int, int, str]] = []
        skipped_placeholder = 0
        for line in bed_path.read_text(encoding="utf-8").splitlines():
            cols = line.split("\t")
            s, e, label = int(cols[1]), int(cols[2]), cols[3]
            if prot_id == "PKZILLA-3" and s == e == 24346 and label == "C1":
                skipped_placeholder += 1
                continue
            truth_domains.append((s, e, label))

        prot_hits = hits_by_prot.get(prot_id, [])
        pred_ivs = [(s, e) for s, e, _, _ in prot_hits]
        recalled_truth_domains = 0
        for ts, te, _ in truth_domains:
            if intersection_bp(pred_ivs, [(ts, te)]) > 0:
                recalled_truth_domains += 1
        truth_union = union_bp([(ts, te) for ts, te, _ in truth_domains])
        pred_in_truth = intersection_bp(pred_ivs, [(ts, te) for ts, te, _ in truth_domains])
        rows.append(
            {
                "protein_id": prot_id,
                "coordinate_space": "protein_relative_aa",
                "truth_domains": len(truth_domains),
                "skipped_zero_length_placeholder": skipped_placeholder,
                "truth_domains_recalled": recalled_truth_domains,
                "truth_domain_recall": round(recalled_truth_domains / len(truth_domains), 6),
                "truth_domain_aa_union": truth_union,
                "pred_domain_aa_overlap": pred_in_truth,
                "pred_track_a_domain_calls": len(prot_hits),
                "pred_ks_ca_calls": sum(
                    1 for _, _, acc, _ in prot_hits if acc in (det.KS_ACCESSIONS | det.CA_ACCESSIONS)
                ),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    ev = subparsers.add_parser("evaluate", help="Evaluate frozen Track A predictions against 12B1 truth.")
    ev.add_argument("--hotspots-tsv", type=Path, required=True)
    ev.add_argument("--truth-bed", type=Path, required=True)
    ev.add_argument("--truth-tsv", type=Path, required=True)
    ev.add_argument("--seeds-manifest", type=Path, required=True)
    ev.add_argument("--pufa-faa", type=Path, required=True)
    ev.add_argument("--pkzilla-faa", type=Path, required=True)
    ev.add_argument("--pkzilla-domains-dir", type=Path, required=True)
    ev.add_argument("--hmm", type=Path, required=True)
    ev.add_argument("--pfam-list", type=Path, required=True)
    ev.add_argument("--ehuxleyi-fasta", type=Path, default=None)
    ev.add_argument("--outdir", type=Path, required=True)
    ev.add_argument("--cpus", type=int, default=2)

    sens = subparsers.add_parser("sensitivity", help="Sweep d_max and min_domains on existing domains.tsv.")
    sens.add_argument("--domains-tsv", type=Path, required=True)
    sens.add_argument("--truth-bed", type=Path, required=True)
    sens.add_argument("--truth-tsv", type=Path, required=True)
    sens.add_argument("--seeds-manifest", type=Path, required=True)
    sens.add_argument("--d-max-values", default="5000,10000,15000,25000,50000")
    sens.add_argument("--min-domains-values", default="3,5,8,10")
    sens.add_argument("--output", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    merged_bed = parse_truth_bed(args.truth_bed)
    truth_tsv_rows = parse_truth_tsv(args.truth_tsv)
    verify_bed_and_tsv_consistency(merged_bed, truth_tsv_rows)
    unmerged_loci = parse_unmerged_pkzilla_loci(args.seeds_manifest)

    if args.command == "evaluate":
        hotspots = read_hotspots_tsv(args.hotspots_tsv)
        summary, locus_rows, merged_rows = evaluate_predictions(
            hotspots, unmerged_loci, merged_bed
        )
        pfam_table = det.load_pfam_table(args.pfam_list)
        control_rows = evaluate_pufa_fas_controls(
            args.pufa_faa,
            args.hmm,
            pfam_table,
            ehuxleyi_fasta=args.ehuxleyi_fasta,
            cpus=args.cpus,
        )
        protein_rows = evaluate_pkzilla_protein_domains(
            args.pkzilla_faa,
            args.pkzilla_domains_dir,
            args.hmm,
            pfam_table,
            cpus=args.cpus,
        )

        args.outdir.mkdir(parents=True, exist_ok=True)
        write_tsv(args.outdir / "locus_metrics.tsv", locus_rows)
        write_tsv(args.outdir / "merged_bed_metrics.tsv", merged_rows)
        write_tsv(args.outdir / "control_eval.tsv", control_rows)
        write_tsv(args.outdir / "protein_domain_recall.tsv", protein_rows)

        payload = {
            "format": "hapto_track_a_benchmark_v1",
            "evaluation_scope": "frozen_baseline_pparvum_12B1_v1_1",
            "pkzilla_b1_note": (
                "PKZILLA-B1 originates from Prymnesium parvum RCC3426, which has no genome "
                "assembly in config/samples.tsv or data/. Therefore genomic benchmark covers "
                "12B1 PKZILLA-1/-2/-3 (plus all 64 merged 12B1 truth BED regions), while "
                "PKZILLA-B1 is evaluated separately at the protein-relative domain level."
            ),
            "summary": summary,
            "loci": locus_rows,
            "pufa_fas_controls_separate": control_rows,
            "protein_relative_pkzilla_check": protein_rows,
        }
        summary_path = args.outdir / "benchmark_summary.json"
        summary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(
            f"12B1 benchmark PASS: recall={summary['pkzilla_recall_same_strand']}/3 "
            f"min_cov={summary['pkzilla_min_unmerged_coverage_same_strand']:.4f} "
            f"merged_64_recall_same={summary['merged_truth_recall_same_strand']}/64 "
            f"merged_64_recall_any={summary['merged_truth_recall_any_strand']}/64 "
            f"outside_truth_bp={summary['predicted_outside_truth_bp_same_strand']} "
            f"summary={summary_path}"
        )
        return 0

    hits = det.read_domains_tsv(args.domains_tsv)
    d_values = [int(x.strip()) for x in args.d_max_values.split(",") if x.strip()]
    n_values = [int(x.strip()) for x in args.min_domains_values.split(",") if x.strip()]
    sens_rows: list[dict[str, object]] = []
    for d_max in d_values:
        for min_dom in n_values:
            hotspots = det.cluster_track_a_hits(hits, d_max_bp=d_max, min_domains=min_dom)
            summary, locus_rows, _ = evaluate_predictions(hotspots, unmerged_loci, merged_bed)
            by_id = {str(r["locus_id"]): r for r in locus_rows}
            sens_rows.append(
                {
                    "is_frozen_baseline": d_max == 15000 and min_dom == 5,
                    "d_max_bp": d_max,
                    "min_domains": min_dom,
                    "pkzilla_recall_same_strand": summary["pkzilla_recall_same_strand"],
                    "pkzilla1_cov": by_id["PKZILLA-1"]["unmerged_coverage_same_strand"],
                    "pkzilla1_components": by_id["PKZILLA-1"]["components_same_strand"],
                    "pkzilla2_cov": by_id["PKZILLA-2"]["unmerged_coverage_same_strand"],
                    "pkzilla2_components": by_id["PKZILLA-2"]["components_same_strand"],
                    "pkzilla3_cov": by_id["PKZILLA-3"]["unmerged_coverage_same_strand"],
                    "pkzilla3_components": by_id["PKZILLA-3"]["components_same_strand"],
                    "merged_64_recall_same": summary["merged_truth_recall_same_strand"],
                    "predicted_hotspots_total": summary["predicted_hotspots_total"],
                    "outside_truth_bp_same_strand": summary["predicted_outside_truth_bp_same_strand"],
                }
            )
    write_tsv(args.output, sens_rows)
    print(f"Sensitivity sweep PASS: configs={len(sens_rows)} output={args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

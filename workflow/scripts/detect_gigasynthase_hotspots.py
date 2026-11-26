#!/usr/bin/env python3
"""GFF-free Track A gigasynthase hotspot detector using 6-frame ORFs and PyHMMER."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import re
import time
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional

import Bio
import pyhmmer
from Bio import SeqIO
from Bio.Seq import Seq


KS_ACCESSIONS = frozenset({"PF00109.33", "PF02801.29", "PF16197.13", "PF13723.14"})
CA_ACCESSIONS = frozenset({"PF00668.26", "PF00501.35"})
CLUSTER_TRACK = "A"
CLUSTER_ROLES = frozenset({"backbone", "ancillary"})
VALID_DNA_RE = re.compile(r"[^ACGTURYKMSWBDHVN]")

DOMAINS_FIELDS = [
    "contig",
    "g_start",
    "g_end",
    "strand",
    "frame",
    "orf_id",
    "orf_g_start",
    "orf_g_end",
    "orf_len_aa",
    "orf_at_contig_left",
    "orf_at_contig_right",
    "aa_start",
    "aa_end",
    "env_start",
    "env_end",
    "accession",
    "model_name",
    "track",
    "role",
    "i_evalue",
    "c_evalue",
    "dom_bitscore",
    "hit_evalue",
    "hit_bitscore",
    "ga_seq_threshold",
    "ga_dom_threshold",
    "threshold_mode",
]

HOTSPOTS_FIELDS = [
    "hotspot_id",
    "contig",
    "g_start",
    "g_end",
    "length_bp",
    "strand",
    "n_domains",
    "n_backbone",
    "n_ks_ca_backbone",
    "n_ks",
    "n_ca",
    "n_ancillary",
    "n_unique_models",
    "n_orfs",
    "frames",
    "edge_truncated",
    "completeness_label",
    "edge_evidence",
    "max_dom_bitscore",
    "sum_dom_bitscore",
    "min_i_evalue",
    "n_neg_control_overlaps",
    "neg_control_models",
    "n_track_b_overlaps",
    "models_summary",
]


@dataclass(frozen=True)
class PfamModelInfo:
    accession: str
    model_name: str
    track: str
    role: str
    rationale: str


@dataclass(frozen=True)
class OrfRecord:
    orf_id: str
    contig: str
    contig_index: int
    contig_length: int
    strand: str
    frame: int
    g_start: int
    g_end: int
    len_aa: int
    at_contig_left: bool
    at_contig_right: bool
    aa_seq: str


@dataclass(frozen=True)
class DomainHit:
    contig: str
    contig_index: int
    g_start: int
    g_end: int
    strand: str
    frame: int
    orf_id: str
    orf_g_start: int
    orf_g_end: int
    orf_len_aa: int
    orf_at_contig_left: bool
    orf_at_contig_right: bool
    aa_start: int
    aa_end: int
    env_start: int
    env_end: int
    accession: str
    model_name: str
    track: str
    role: str
    i_evalue: float
    c_evalue: float
    dom_bitscore: float
    hit_evalue: float
    hit_bitscore: float
    ga_seq_threshold: float
    ga_dom_threshold: float
    threshold_mode: str = "gathering"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_metadata(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def as_str(value: object) -> str:
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8")
    return str(value)


def normalize_dna(sequence: str) -> str:
    upper = sequence.upper().replace("U", "T")
    return VALID_DNA_RE.sub("N", upper)


def load_pfam_table(path: Path) -> dict[str, PfamModelInfo]:
    models: dict[str, PfamModelInfo] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"accession", "model_name", "track", "role", "rationale"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        for row in reader:
            acc = row["accession"].strip()
            if acc in models:
                raise ValueError(f"Duplicate accession in {path}: {acc}")
            models[acc] = PfamModelInfo(
                accession=acc,
                model_name=row["model_name"].strip(),
                track=row["track"].strip(),
                role=row["role"].strip(),
                rationale=row["rationale"].strip(),
            )
    missing_anchors = (KS_ACCESSIONS | CA_ACCESSIONS) - set(models)
    if missing_anchors:
        raise ValueError(f"Required KS/C/A anchor accessions absent from {path}: {sorted(missing_anchors)}")
    return models


def aa_interval_to_genomic(
    orf_g_start: int,
    orf_g_end: int,
    strand: str,
    aa_start: int,
    aa_end: int,
) -> tuple[int, int]:
    """Map 0-based half-open ORF AA coordinates [aa_start, aa_end) to 0-based half-open genomic coordinates."""
    if aa_start < 0 or aa_end <= aa_start:
        raise ValueError(f"Invalid AA interval [{aa_start}, {aa_end})")
    orf_bp = orf_g_end - orf_g_start
    if 3 * aa_end > orf_bp:
        raise ValueError(
            f"AA interval [{aa_start}, {aa_end}) exceeds ORF span {orf_bp} bp ({orf_bp // 3} aa)"
        )
    if strand == "+":
        return orf_g_start + 3 * aa_start, orf_g_start + 3 * aa_end
    if strand == "-":
        return orf_g_end - 3 * aa_end, orf_g_end - 3 * aa_start
    raise ValueError(f"Invalid strand: {strand!r}")


def extract_six_frame_orfs(
    contig: str,
    sequence: str,
    *,
    contig_index: int = 0,
    min_orf_aa: int = 150,
    genetic_code: int = 1,
) -> Iterator[OrfRecord]:
    """Yield lossless stop-to-stop 6-frame ORFs (>= min_orf_aa, no ATG required)."""
    if min_orf_aa < 1:
        raise ValueError(f"min_orf_aa must be >= 1, got {min_orf_aa}")
    clean = normalize_dna(sequence)
    contig_len = len(clean)
    if contig_len < 3 * min_orf_aa:
        return

    fwd_seq = Seq(clean)
    rc_seq = fwd_seq.reverse_complement()

    for strand, strand_seq in (("+", fwd_seq), ("-", rc_seq)):
        for frame in (0, 1, 2):
            n_codons = (contig_len - frame) // 3
            if n_codons < min_orf_aa:
                continue
            coding_end = frame + 3 * n_codons
            translated = str(
                strand_seq[frame:coding_end].translate(table=genetic_code, stop_symbol="*")
            )
            pos = 0
            n_aa = len(translated)
            while pos < n_aa:
                stop_idx = translated.find("*", pos)
                end_idx = n_aa if stop_idx == -1 else stop_idx
                seg_len = end_idx - pos
                if seg_len >= min_orf_aa:
                    seg_aa = translated[pos:end_idx]
                    s_start = frame + 3 * pos
                    s_end = frame + 3 * end_idx
                    if strand == "+":
                        g_start = s_start
                        g_end = s_end
                        at_left = pos == 0
                        at_right = end_idx == n_codons
                    else:
                        g_start = contig_len - s_end
                        g_end = contig_len - s_start
                        at_left = end_idx == n_codons
                        at_right = pos == 0
                    orf_id = f"{contig}|{strand}f{frame}|{g_start}-{g_end}"
                    yield OrfRecord(
                        orf_id=orf_id,
                        contig=contig,
                        contig_index=contig_index,
                        contig_length=contig_len,
                        strand=strand,
                        frame=frame,
                        g_start=g_start,
                        g_end=g_end,
                        len_aa=seg_len,
                        at_contig_left=at_left,
                        at_contig_right=at_right,
                        aa_seq=seg_aa,
                    )
                pos = end_idx + 1


def iter_fasta_orfs(
    fasta_path: Path,
    *,
    min_orf_aa: int = 150,
    genetic_code: int = 1,
) -> tuple[Iterator[OrfRecord], dict[str, int]]:
    """Return a streaming ORF iterator and mutable contig-stats accumulator."""
    stats = {"contigs": 0, "total_bp": 0, "orfs_extracted": 0}

    def _generator() -> Iterator[OrfRecord]:
        with open_text_maybe_gzip(fasta_path) as handle:
            for contig_idx, record in enumerate(SeqIO.parse(handle, "fasta")):
                seq_str = str(record.seq)
                stats["contigs"] += 1
                stats["total_bp"] += len(seq_str)
                for orf in extract_six_frame_orfs(
                    record.id,
                    seq_str,
                    contig_index=contig_idx,
                    min_orf_aa=min_orf_aa,
                    genetic_code=genetic_code,
                ):
                    stats["orfs_extracted"] += 1
                    yield orf

    return _generator(), stats


def load_hmms(hmm_path: Path, pfam_table: dict[str, PfamModelInfo]):
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        hmms = list(hmm_file)
    if not hmms:
        raise ValueError(f"No HMMs loaded from {hmm_path}")
    ga_cutoffs: dict[str, tuple[float, float]] = {}
    for hmm in hmms:
        acc = as_str(hmm.accession)
        name = as_str(hmm.name)
        if acc not in pfam_table:
            raise ValueError(f"HMM accession {acc} ({name}) missing from Pfam selection TSV")
        if pfam_table[acc].model_name != name:
            raise ValueError(
                f"HMM name mismatch for {acc}: hmm={name} table={pfam_table[acc].model_name}"
            )
        ga = hmm.cutoffs.gathering
        if ga is None:
            raise ValueError(f"HMM {acc} ({name}) lacks gathering (GA) cutoffs")
        ga_cutoffs[acc] = (float(ga[0]), float(ga[1]))
    return hmms, ga_cutoffs


def scan_orf_batch(
    orfs: list[OrfRecord],
    hmms,
    pfam_table: dict[str, PfamModelInfo],
    ga_cutoffs: dict[str, tuple[float, float]],
    *,
    cpus: int,
) -> list[DomainHit]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    by_id = {orf.orf_id: orf for orf in orfs}
    digital_seqs = [
        pyhmmer.easel.TextSequence(name=orf.orf_id, sequence=orf.aa_seq).digitize(alphabet)
        for orf in orfs
    ]
    hits_out: list[DomainHit] = []
    for top_hits in pyhmmer.hmmscan(digital_seqs, hmms, bit_cutoffs="gathering", cpus=cpus):
        orf = by_id[as_str(top_hits.query.name)]
        for hit in top_hits:
            if not hit.included:
                continue
            acc = as_str(hit.accession)
            model_name = as_str(hit.name)
            info = pfam_table[acc]
            ga_seq, ga_dom = ga_cutoffs[acc]
            for dom in hit.domains.included:
                aln = dom.alignment
                aa_start = int(aln.target_from) - 1
                aa_end = int(aln.target_to)
                env_start = int(dom.env_from) - 1
                env_end = int(dom.env_to)
                g_start, g_end = aa_interval_to_genomic(
                    orf.g_start, orf.g_end, orf.strand, aa_start, aa_end
                )
                hits_out.append(
                    DomainHit(
                        contig=orf.contig,
                        contig_index=orf.contig_index,
                        g_start=g_start,
                        g_end=g_end,
                        strand=orf.strand,
                        frame=orf.frame,
                        orf_id=orf.orf_id,
                        orf_g_start=orf.g_start,
                        orf_g_end=orf.g_end,
                        orf_len_aa=orf.len_aa,
                        orf_at_contig_left=orf.at_contig_left,
                        orf_at_contig_right=orf.at_contig_right,
                        aa_start=aa_start,
                        aa_end=aa_end,
                        env_start=env_start,
                        env_end=env_end,
                        accession=acc,
                        model_name=model_name,
                        track=info.track,
                        role=info.role,
                        i_evalue=float(dom.i_evalue),
                        c_evalue=float(dom.c_evalue),
                        dom_bitscore=float(dom.score),
                        hit_evalue=float(hit.evalue),
                        hit_bitscore=float(hit.score),
                        ga_seq_threshold=ga_seq,
                        ga_dom_threshold=ga_dom,
                    )
                )
    return hits_out


def scan_fasta_to_hits(
    fasta_path: Path,
    hmm_path: Path,
    pfam_table: dict[str, PfamModelInfo],
    *,
    min_orf_aa: int = 150,
    genetic_code: int = 1,
    batch_size: int = 2000,
    cpus: int = 1,
) -> tuple[list[DomainHit], dict[str, int]]:
    hmms, ga_cutoffs = load_hmms(hmm_path, pfam_table)
    orf_iter, stats = iter_fasta_orfs(
        fasta_path, min_orf_aa=min_orf_aa, genetic_code=genetic_code
    )
    all_hits: list[DomainHit] = []
    batch: list[OrfRecord] = []
    for orf in orf_iter:
        batch.append(orf)
        if len(batch) >= batch_size:
            all_hits.extend(
                scan_orf_batch(batch, hmms, pfam_table, ga_cutoffs, cpus=cpus)
            )
            batch.clear()
    if batch:
        all_hits.extend(
            scan_orf_batch(batch, hmms, pfam_table, ga_cutoffs, cpus=cpus)
        )
    all_hits.sort(
        key=lambda h: (
            h.contig_index,
            h.contig,
            h.g_start,
            h.g_end,
            h.strand,
            h.frame,
            h.accession,
            h.aa_start,
            h.aa_end,
        )
    )
    return all_hits, stats


def cluster_track_a_hits(
    hits: Iterable[DomainHit],
    *,
    d_max_bp: int = 15000,
    min_domains: int = 5,
    require_ks_or_ca: bool = True,
    hotspot_prefix: str = "H",
) -> list[dict[str, object]]:
    """Cluster Track A (backbone + ancillary) domain hits on same contig & strand."""
    hit_list = list(hits)
    track_a_hits = [
        h for h in hit_list if h.track == CLUSTER_TRACK and h.role in CLUSTER_ROLES
    ]
    neg_controls_by_cs: dict[tuple[str, str], list[DomainHit]] = {}
    track_b_by_cs: dict[tuple[str, str], list[DomainHit]] = {}
    for h in hit_list:
        if h.track == "A" and h.role == "negative_control":
            neg_controls_by_cs.setdefault((h.contig, h.strand), []).append(h)
        elif h.track == "B":
            track_b_by_cs.setdefault((h.contig, h.strand), []).append(h)

    grouped: dict[tuple[str, str], list[DomainHit]] = {}
    for h in track_a_hits:
        grouped.setdefault((h.contig, h.strand), []).append(h)

    raw_clusters: list[tuple[int, str, int, int, str, list[DomainHit]]] = []
    for (contig, strand), cs_hits in grouped.items():
        cs_hits.sort(key=lambda h: (h.g_start, h.g_end, h.frame, h.accession, h.aa_start))
        cluster: list[DomainHit] = [cs_hits[0]]
        c_contig_idx = cs_hits[0].contig_index
        c_start = cs_hits[0].g_start
        c_max_end = cs_hits[0].g_end
        for hit in cs_hits[1:]:
            if hit.g_start - c_max_end <= d_max_bp:
                cluster.append(hit)
                if hit.g_end > c_max_end:
                    c_max_end = hit.g_end
            else:
                raw_clusters.append((c_contig_idx, contig, c_start, c_max_end, strand, cluster))
                cluster = [hit]
                c_contig_idx = hit.contig_index
                c_start = hit.g_start
                c_max_end = hit.g_end
        raw_clusters.append((c_contig_idx, contig, c_start, c_max_end, strand, cluster))

    raw_clusters.sort(key=lambda item: (item[0], item[1], item[2], item[3], item[4]))

    hotspots: list[dict[str, object]] = []
    for _, contig, g_start, g_end, strand, cluster in raw_clusters:
        n_domains = len(cluster)
        if n_domains < min_domains:
            continue
        n_ks = sum(1 for h in cluster if h.accession in KS_ACCESSIONS)
        n_ca = sum(1 for h in cluster if h.accession in CA_ACCESSIONS)
        n_ks_ca = n_ks + n_ca
        if require_ks_or_ca and n_ks_ca < 1:
            continue

        n_backbone = sum(1 for h in cluster if h.role == "backbone")
        n_ancillary = sum(1 for h in cluster if h.role == "ancillary")
        edge_items: list[str] = []
        seen_edge: set[str] = set()
        for h in cluster:
            if h.orf_at_contig_left:
                tag = f"contig_left:{h.orf_id}"
                if tag not in seen_edge:
                    seen_edge.add(tag)
                    edge_items.append(tag)
            if h.orf_at_contig_right:
                tag = f"contig_right:{h.orf_id}"
                if tag not in seen_edge:
                    seen_edge.add(tag)
                    edge_items.append(tag)
        edge_truncated = bool(edge_items)

        neg_overlaps = [
            nh
            for nh in neg_controls_by_cs.get((contig, strand), [])
            if not (nh.g_end <= g_start or nh.g_start >= g_end)
        ]
        b_overlaps = [
            bh
            for bh in track_b_by_cs.get((contig, strand), [])
            if not (bh.g_end <= g_start or bh.g_start >= g_end)
        ]

        model_counts = Counter(h.model_name for h in cluster)
        models_summary = ";".join(
            f"{name}:{count}" for name, count in sorted(model_counts.items())
        )
        neg_models = sorted({nh.model_name for nh in neg_overlaps})
        hotspot_id = f"{hotspot_prefix}{len(hotspots) + 1:04d}"
        hotspots.append(
            {
                "hotspot_id": hotspot_id,
                "contig": contig,
                "g_start": g_start,
                "g_end": g_end,
                "length_bp": g_end - g_start,
                "strand": strand,
                "n_domains": n_domains,
                "n_backbone": n_backbone,
                "n_ks_ca_backbone": n_ks_ca,
                "n_ks": n_ks,
                "n_ca": n_ca,
                "n_ancillary": n_ancillary,
                "n_unique_models": len(model_counts),
                "n_orfs": len({h.orf_id for h in cluster}),
                "frames": ",".join(str(f) for f in sorted({h.frame for h in cluster})),
                "edge_truncated": "true" if edge_truncated else "false",
                "completeness_label": "Edge-truncated" if edge_truncated else "Internal",
                "edge_evidence": ";".join(edge_items) if edge_items else "none",
                "max_dom_bitscore": f"{max(h.dom_bitscore for h in cluster):.2f}",
                "sum_dom_bitscore": f"{sum(h.dom_bitscore for h in cluster):.2f}",
                "min_i_evalue": f"{min(h.i_evalue for h in cluster):.3e}",
                "n_neg_control_overlaps": len(neg_overlaps),
                "neg_control_models": ",".join(neg_models) if neg_models else "none",
                "n_track_b_overlaps": len(b_overlaps),
                "models_summary": models_summary,
            }
        )
    return hotspots


def write_domains_tsv(path: Path, hits: Iterable[DomainHit]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DOMAINS_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for hit in hits:
            row = asdict(hit)
            del row["contig_index"]
            row["orf_at_contig_left"] = "true" if hit.orf_at_contig_left else "false"
            row["orf_at_contig_right"] = "true" if hit.orf_at_contig_right else "false"
            row["i_evalue"] = f"{hit.i_evalue:.3e}"
            row["c_evalue"] = f"{hit.c_evalue:.3e}"
            row["dom_bitscore"] = f"{hit.dom_bitscore:.2f}"
            row["hit_evalue"] = f"{hit.hit_evalue:.3e}"
            row["hit_bitscore"] = f"{hit.hit_bitscore:.2f}"
            row["ga_seq_threshold"] = f"{hit.ga_seq_threshold:.2f}"
            row["ga_dom_threshold"] = f"{hit.ga_dom_threshold:.2f}"
            writer.writerow(row)


def read_domains_tsv(path: Path) -> list[DomainHit]:
    hits: list[DomainHit] = []
    contig_order: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != DOMAINS_FIELDS:
            raise ValueError(f"Unexpected domains.tsv columns in {path}: {reader.fieldnames}")
        for row in reader:
            contig = row["contig"]
            if contig not in contig_order:
                contig_order[contig] = len(contig_order)
            hits.append(
                DomainHit(
                    contig=contig,
                    contig_index=contig_order[contig],
                    g_start=int(row["g_start"]),
                    g_end=int(row["g_end"]),
                    strand=row["strand"],
                    frame=int(row["frame"]),
                    orf_id=row["orf_id"],
                    orf_g_start=int(row["orf_g_start"]),
                    orf_g_end=int(row["orf_g_end"]),
                    orf_len_aa=int(row["orf_len_aa"]),
                    orf_at_contig_left=row["orf_at_contig_left"] == "true",
                    orf_at_contig_right=row["orf_at_contig_right"] == "true",
                    aa_start=int(row["aa_start"]),
                    aa_end=int(row["aa_end"]),
                    env_start=int(row["env_start"]),
                    env_end=int(row["env_end"]),
                    accession=row["accession"],
                    model_name=row["model_name"],
                    track=row["track"],
                    role=row["role"],
                    i_evalue=float(row["i_evalue"]),
                    c_evalue=float(row["c_evalue"]),
                    dom_bitscore=float(row["dom_bitscore"]),
                    hit_evalue=float(row["hit_evalue"]),
                    hit_bitscore=float(row["hit_bitscore"]),
                    ga_seq_threshold=float(row["ga_seq_threshold"]),
                    ga_dom_threshold=float(row["ga_dom_threshold"]),
                    threshold_mode=row["threshold_mode"],
                )
            )
    return hits


def write_hotspots_outputs(
    outdir: Path, hotspots: list[dict[str, object]]
) -> tuple[Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / "hotspots.tsv"
    bed_path = outdir / "hotspots.bed"
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=HOTSPOTS_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(hotspots)

    with bed_path.open("w", encoding="utf-8") as handle:
        for row in hotspots:
            score = min(1000, max(0, int(round(float(str(row["sum_dom_bitscore"]))))))
            name = (
                f"{row['hotspot_id']}|N={row['n_domains']}"
                f"|KS_CA={row['n_ks_ca_backbone']}|{row['completeness_label']}"
            )
            handle.write(
                f"{row['contig']}\t{row['g_start']}\t{row['g_end']}\t{name}\t{score}\t{row['strand']}\n"
            )
    return tsv_path, bed_path


def write_metadata_json(
    path: Path,
    *,
    mode: str,
    parameters: dict[str, object],
    inputs: dict[str, object],
    hits: list[DomainHit],
    hotspots: list[dict[str, object]],
    scan_stats: Optional[dict[str, int]],
    elapsed_seconds: float,
    unfiltered_warning: Optional[str] = None,
) -> None:
    role_counts = Counter(f"{h.track}/{h.role}" for h in hits)
    payload: dict[str, object] = {
        "format": "hapto_gigasynthase_track_a_v1",
        "mode": mode,
        "script": "workflow/scripts/detect_gigasynthase_hotspots.py",
        "versions": {
            "python": platform.python_version(),
            "biopython": Bio.__version__,
            "pyhmmer": pyhmmer.__version__,
        },
        "parameters": parameters,
        "inputs": inputs,
        "summary": {
            "domain_hits_total": len(hits),
            "domain_hits_by_track_role": dict(sorted(role_counts.items())),
            "hotspots_total": len(hotspots),
            "hotspots_edge_truncated": sum(
                1 for row in hotspots if row["edge_truncated"] == "true"
            ),
            "hotspots_with_neg_control_overlaps": sum(
                1 for row in hotspots if int(str(row["n_neg_control_overlaps"])) > 0
            ),
            "elapsed_seconds": round(elapsed_seconds, 3),
        },
    }
    if scan_stats is not None:
        payload["summary"]["contigs"] = scan_stats["contigs"]  # type: ignore[index]
        payload["summary"]["total_bp"] = scan_stats["total_bp"]  # type: ignore[index]
        payload["summary"]["orfs_ge_min_aa"] = scan_stats["orfs_extracted"]  # type: ignore[index]
    if unfiltered_warning:
        payload["status_note"] = unfiltered_warning
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    tmp.replace(path)


def add_cluster_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--d-max", type=int, default=15000)
    parser.add_argument("--min-domains", type=int, default=5)
    parser.add_argument(
        "--no-require-ks-or-ca",
        action="store_true",
        help="Disable the baseline requirement of >=1 KS or C/A backbone domain hit.",
    )
    parser.add_argument("--hotspot-prefix", default="H")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser("scan", help="Extract 6-frame ORFs, scan HMMs, and cluster Track A hotspots.")
    scan.add_argument("--fasta", type=Path, required=True)
    scan.add_argument("--hmm", type=Path, required=True)
    scan.add_argument("--pfam-list", type=Path, required=True)
    scan.add_argument("--outdir", type=Path, required=True)
    scan.add_argument("--min-orf-aa", type=int, default=150)
    scan.add_argument("--genetic-code", type=int, default=1)
    scan.add_argument("--batch-size", type=int, default=2000)
    scan.add_argument("--cpus", type=int, default=4)
    scan.add_argument("--unfiltered-warning", default="")
    add_cluster_args(scan)

    recluster = subparsers.add_parser("recluster", help="Re-cluster an existing domains.tsv without re-scanning.")
    recluster.add_argument("--domains", type=Path, required=True)
    recluster.add_argument("--outdir", type=Path, required=True)
    add_cluster_args(recluster)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    started = time.monotonic()
    require_ks_or_ca = not args.no_require_ks_or_ca

    if args.command == "scan":
        pfam_table = load_pfam_table(args.pfam_list)
        hits, scan_stats = scan_fasta_to_hits(
            args.fasta,
            args.hmm,
            pfam_table,
            min_orf_aa=args.min_orf_aa,
            genetic_code=args.genetic_code,
            batch_size=args.batch_size,
            cpus=args.cpus,
        )
        domains_path = args.outdir / "domains.tsv"
        write_domains_tsv(domains_path, hits)
        hotspots = cluster_track_a_hits(
            hits,
            d_max_bp=args.d_max,
            min_domains=args.min_domains,
            require_ks_or_ca=require_ks_or_ca,
            hotspot_prefix=args.hotspot_prefix,
        )
        hotspots_tsv, hotspots_bed = write_hotspots_outputs(args.outdir, hotspots)
        elapsed = time.monotonic() - started
        metadata_path = args.outdir / "metadata.json"
        write_metadata_json(
            metadata_path,
            mode="scan",
            parameters={
                "min_orf_aa": args.min_orf_aa,
                "genetic_code": args.genetic_code,
                "batch_size": args.batch_size,
                "cpus": args.cpus,
                "d_max_bp": args.d_max,
                "min_domains": args.min_domains,
                "require_ks_or_ca": require_ks_or_ca,
                "threshold_mode": "gathering",
                "clustered_track": CLUSTER_TRACK,
                "clustered_roles": sorted(CLUSTER_ROLES),
                "ks_accessions": sorted(KS_ACCESSIONS),
                "ca_accessions": sorted(CA_ACCESSIONS),
            },
            inputs={
                "fasta": file_metadata(args.fasta),
                "hmm": file_metadata(args.hmm),
                "pfam_list": file_metadata(args.pfam_list),
            },
            hits=hits,
            hotspots=hotspots,
            scan_stats=scan_stats,
            elapsed_seconds=elapsed,
            unfiltered_warning=args.unfiltered_warning or None,
        )
        print(
            f"Track A scan PASS: contigs={scan_stats['contigs']} bp={scan_stats['total_bp']} "
            f"orfs={scan_stats['orfs_extracted']} domain_hits={len(hits)} "
            f"hotspots={len(hotspots)} elapsed={elapsed:.2f}s "
            f"outputs=({domains_path}, {hotspots_tsv}, {hotspots_bed}, {metadata_path})"
        )
        return 0

    hits = read_domains_tsv(args.domains)
    hotspots = cluster_track_a_hits(
        hits,
        d_max_bp=args.d_max,
        min_domains=args.min_domains,
        require_ks_or_ca=require_ks_or_ca,
        hotspot_prefix=args.hotspot_prefix,
    )
    hotspots_tsv, hotspots_bed = write_hotspots_outputs(args.outdir, hotspots)
    elapsed = time.monotonic() - started
    metadata_path = args.outdir / "metadata.json"
    write_metadata_json(
        metadata_path,
        mode="recluster",
        parameters={
            "d_max_bp": args.d_max,
            "min_domains": args.min_domains,
            "require_ks_or_ca": require_ks_or_ca,
            "threshold_mode": "gathering",
            "clustered_track": CLUSTER_TRACK,
            "clustered_roles": sorted(CLUSTER_ROLES),
        },
        inputs={"domains_tsv": file_metadata(args.domains)},
        hits=hits,
        hotspots=hotspots,
        scan_stats=None,
        elapsed_seconds=elapsed,
    )
    print(
        f"Track A recluster PASS: domain_hits={len(hits)} hotspots={len(hotspots)} "
        f"outputs=({hotspots_tsv}, {hotspots_bed}, {metadata_path})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

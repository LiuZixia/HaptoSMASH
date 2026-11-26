#!/usr/bin/env python3
"""Unit and synthetic integration tests for Track A 6-frame ORF and hotspot detection."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import detect_gigasynthase_hotspots as det  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
PFAM_LIST = REPO_ROOT / "resources/pfam_lists/hapto_core_pfam.tsv"

CODON_FOR_AA = {
    "A": "GCT",
    "C": "TGC",
    "D": "GAC",
    "E": "GAA",
    "F": "TTC",
    "G": "GGT",
    "H": "CAC",
    "I": "ATC",
    "K": "AAG",
    "L": "CTG",
    "M": "ATG",
    "N": "AAC",
    "P": "CCT",
    "Q": "CAG",
    "R": "CGT",
    "S": "TCC",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAC",
    "X": "NNN",
}


def backtranslate(aa: str) -> str:
    return "".join(CODON_FOR_AA[ch] for ch in aa)


def resolve_hmm_path() -> Path | None:
    for candidate in (
        REPO_ROOT / "db/custom_hmm/hapto_core.hmm",
        REPO_ROOT / "tmp/hmm_repro/hapto_core.hmm",
    ):
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def make_hit(
    *,
    contig: str = "ctg1",
    g_start: int,
    g_end: int,
    strand: str = "+",
    frame: int = 0,
    accession: str = "PF00109.33",
    model_name: str = "ketoacyl-synt",
    track: str = "A",
    role: str = "backbone",
    orf_id: str | None = None,
    orf_at_left: bool = False,
    orf_at_right: bool = False,
) -> det.DomainHit:
    oid = orf_id or f"{contig}|{strand}f{frame}|{g_start}-{g_end}"
    return det.DomainHit(
        contig=contig,
        contig_index=0,
        g_start=g_start,
        g_end=g_end,
        strand=strand,
        frame=frame,
        orf_id=oid,
        orf_g_start=g_start,
        orf_g_end=g_end,
        orf_len_aa=max(150, (g_end - g_start) // 3),
        orf_at_contig_left=orf_at_left,
        orf_at_contig_right=orf_at_right,
        aa_start=0,
        aa_end=max(1, (g_end - g_start) // 3),
        env_start=0,
        env_end=max(1, (g_end - g_start) // 3),
        accession=accession,
        model_name=model_name,
        track=track,
        role=role,
        i_evalue=1e-20,
        c_evalue=1e-21,
        dom_bitscore=100.0,
        hit_evalue=1e-20,
        hit_bitscore=100.0,
        ga_seq_threshold=25.0,
        ga_dom_threshold=25.0,
    )


def test_six_frame_forward_reverse_lossless_and_no_atg_required():
    # Construct forward frame 1 ORF without any ATG start codon (160 aa of Ala/Val/Leu)
    fwd_aa = ("AVL" * 54)[:160]
    assert "M" not in fwd_aa
    # Prefix 1 nt + stop TAA (4 nt), ORF (480 nt), stop TAG (3 nt), spacer,
    # then a reverse-strand ORF in frame 2 on RC.
    rev_aa = ("GDE" * 55)[:165]
    assert "M" not in rev_aa
    rev_dna = str(Seq("TAA" + backtranslate(rev_aa) + "TGA").reverse_complement())
    # Assemble contig: [0:1]="A", [1:4]="TAA" (stop in frame 1), [4:484]=fwd_dna, [484:487]="TAG",
    # spacer "ACGT"*5 (20 nt), rev_dna (501 nt), trailing "CC" (2 nt)
    fwd_dna = backtranslate(fwd_aa)
    contig = "A" + "TAA" + fwd_dna + "TAG" + ("ACGT" * 5) + rev_dna + "CC"

    orfs = list(det.extract_six_frame_orfs("ctg_test", contig, min_orf_aa=150))
    fwd_matches = [o for o in orfs if o.strand == "+" and o.aa_seq == fwd_aa]
    rev_matches = [o for o in orfs if o.strand == "-" and o.aa_seq == rev_aa]
    assert len(fwd_matches) == 1
    assert len(rev_matches) == 1

    f_orf = fwd_matches[0]
    assert f_orf.frame == 1
    assert (f_orf.g_start, f_orf.g_end) == (4, 484)
    assert not f_orf.at_contig_left
    assert not f_orf.at_contig_right
    assert str(Seq(contig[f_orf.g_start : f_orf.g_end]).translate()) == fwd_aa

    r_orf = rev_matches[0]
    assert not r_orf.at_contig_left
    assert not r_orf.at_contig_right
    assert (
        str(Seq(contig[r_orf.g_start : r_orf.g_end]).reverse_complement().translate())
        == rev_aa
    )

    # Lossless AA-subinterval mapping for both strands
    f_gs, f_ge = det.aa_interval_to_genomic(f_orf.g_start, f_orf.g_end, "+", 10, 45)
    assert str(Seq(contig[f_gs:f_ge]).translate()) == fwd_aa[10:45]

    r_gs, r_ge = det.aa_interval_to_genomic(r_orf.g_start, r_orf.g_end, "-", 12, 57)
    assert (
        str(Seq(contig[r_gs:r_ge]).reverse_complement().translate()) == rev_aa[12:57]
    )


def test_stop_splitting_min_length_and_contig_edges():
    aa1 = ("S" * 150)
    aa_short = ("T" * 149)
    aa2 = ("V" * 155)
    # Contig starts directly in frame 0 without stop (left edge truncated!),
    # has internal stop, short 149-aa ORF (filtered), stop, then 155-aa ORF running off 3' end (right edge truncated!).
    contig = backtranslate(aa1) + "TAA" + backtranslate(aa_short) + "TAG" + backtranslate(aa2) + "A"

    orfs = [o for o in det.extract_six_frame_orfs("edge_ctg", contig, min_orf_aa=150) if o.strand == "+" and o.frame == 0]
    assert [(o.len_aa, o.g_start, o.g_end, o.at_contig_left, o.at_contig_right) for o in orfs] == [
        (150, 0, 450, True, False),
        (155, 450 + 3 + 447 + 3, len(contig) - 1, False, True),
    ]

    # Reverse strand edge truncation: place the same open-ended ORF on reverse strand at left/right boundaries
    rc_contig = str(Seq(contig).reverse_complement())
    rc_orfs = [o for o in det.extract_six_frame_orfs("edge_ctg_rc", rc_contig, min_orf_aa=150) if o.strand == "-" and o.frame == 0]
    assert len(rc_orfs) == 2
    assert rc_orfs[0].aa_seq == aa1
    assert rc_orfs[0].at_contig_right is True and rc_orfs[0].at_contig_left is False
    assert rc_orfs[1].aa_seq == aa2
    assert rc_orfs[1].at_contig_left is True and rc_orfs[1].at_contig_right is False


def test_clustering_rules_gap_anchor_track_b_exclusion_and_edge_label():
    # 1. Exact d_max=15000 boundary, running max_end with nested interval, and preserved overlapping hits
    hits = [
        make_hit(g_start=1000, g_end=5000, accession="PF00109.33", model_name="ketoacyl-synt", role="backbone"),
        make_hit(g_start=1200, g_end=2500, accession="PF02801.29", model_name="Ketoacyl-synt_C", role="backbone"),
        make_hit(g_start=20000, g_end=21000, accession="PF08659.18", model_name="KR", role="ancillary"),  # gap=15000 -> merges
        make_hit(g_start=20100, g_end=20900, accession="PF00106.32", model_name="adh_short", role="ancillary"),  # overlap preserved
        make_hit(g_start=22000, g_end=23000, accession="PF00550.33", model_name="PP-binding", role="backbone"),
        # Opposite strand at same coordinates must NOT merge
        make_hit(g_start=1000, g_end=23000, strand="-", accession="PF00109.33", model_name="ketoacyl-synt", role="backbone"),
        # Beyond 15000 bp gap (23000 + 15001 = 38001) -> separate, <5 hits so dropped
        make_hit(g_start=38001, g_end=39000, accession="PF00109.33", model_name="ketoacyl-synt", role="backbone"),
        # Track B hits adjacent to cluster must NOT be included
        make_hit(g_start=24000, g_end=25000, accession="PF00305.26", model_name="Lipoxygenase", track="B", role="ancillary"),
        # Negative control hit overlapping cluster: excluded from N_dom, recorded in overlap count
        make_hit(g_start=1500, g_end=2000, accession="PF17828.7", model_name="FAS_N", track="A", role="negative_control"),
    ]
    hotspots = det.cluster_track_a_hits(hits, d_max_bp=15000, min_domains=5)
    assert len(hotspots) == 1
    hs = hotspots[0]
    assert (hs["g_start"], hs["g_end"], hs["strand"]) == (1000, 23000, "+")
    assert hs["n_domains"] == 5
    assert hs["n_ks_ca_backbone"] == 2
    assert hs["completeness_label"] == "Internal"
    assert hs["n_neg_control_overlaps"] == 1
    assert hs["neg_control_models"] == "FAS_N"
    assert hs["n_track_b_overlaps"] == 0

    # 2. Ancillary-only cluster (5 hits without KS or C/A) must NOT be called
    anc_only = [
        make_hit(g_start=1000 + i * 1000, g_end=1500 + i * 1000, accession="PF08659.18", model_name="KR", role="ancillary")
        for i in range(5)
    ]
    assert det.cluster_track_a_hits(anc_only) == []

    # 3. Track B cluster (6 hits) must NOT be called
    track_b_only = [
        make_hit(g_start=1000 + i * 1000, g_end=1500 + i * 1000, accession="PF00305.26", model_name="Lipoxygenase", track="B", role="ancillary")
        for i in range(6)
    ]
    assert det.cluster_track_a_hits(track_b_only) == []

    # 4. Edge-truncated label requires actual terminal ORF boundary evidence
    edge_hits = [
        make_hit(g_start=0, g_end=600, accession="PF00109.33", model_name="ketoacyl-synt", role="backbone", orf_at_left=True),
        *[
            make_hit(g_start=1000 + i * 800, g_end=1500 + i * 800, accession="PF08659.18", model_name="KR", role="ancillary")
            for i in range(4)
        ],
    ]
    edge_hs = det.cluster_track_a_hits(edge_hits)
    assert len(edge_hs) == 1
    assert edge_hs[0]["edge_truncated"] == "true"
    assert edge_hs[0]["completeness_label"] == "Edge-truncated"
    assert "contig_left:" in str(edge_hs[0]["edge_evidence"])


def test_end_to_end_synthetic_scan_and_recluster(tmp_path: Path):
    hmm_path = resolve_hmm_path()
    if hmm_path is None:
        pytest.skip("hapto_core.hmm not accessible in current test environment")

    pk_records = SeqIO.to_dict(SeqIO.parse(REPO_ROOT / "resources/seeds/pkzilla_seeds.faa", "fasta"))
    ox_records = SeqIO.to_dict(SeqIO.parse(REPO_ROOT / "resources/seeds/halogenase_oxylipin_seeds.faa", "fasta"))
    neg_records = SeqIO.to_dict(SeqIO.parse(REPO_ROOT / "resources/seeds/pufa_neg_controls.faa", "fasta"))

    # Contig 1 (+): multi-exon PKZILLA-1 slice (first 6,000 aa -> multiple KS/AT/KR/DH/ACP domains) split by an 80-bp frame-shifting intron
    pk_aa = str(pk_records["PKZILLA-1"].seq)[:6000]
    exon1 = backtranslate(pk_aa[:3000])
    exon2 = backtranslate(pk_aa[3000:])
    intron = "TAACTGACG" * 9  # 81 bp containing stops
    pad = "TAGTGA" * 20
    ctg1 = pad + exon1 + intron + exon2 + pad

    # Contig 2 (-, edge-truncated at left end of contig): 3,500 aa slice of PKZILLA-2 placed on reverse strand touching contig left boundary
    pk2_aa = str(pk_records["PKZILLA-2"].seq)[:3500]
    ctg2 = str(Seq(backtranslate(pk2_aa)).reverse_complement()) + pad

    # Contig 3: Track B only (concatenated spacer-separated LOX_EMIHU copies -> 6 Track B hits, 0 Track A hotspots)
    lox_aa = str(ox_records["LOX_EMIHU"].seq)
    ctg3 = pad + ("TAA" + backtranslate(lox_aa)) * 6 + "TAA" + pad

    # Contig 4: E. huxleyi FAS negative control (XP_005790240.1)
    fas_aa = str(neg_records["XP_005790240.1"].seq)
    ctg4 = pad + backtranslate(fas_aa) + pad

    # Contig 5: E. huxleyi KAS III negative control (XP_005793321.1 -> A/negative_control only)
    kas3_aa = str(neg_records["XP_005793321.1"].seq)
    ctg5 = pad + backtranslate(kas3_aa) + pad

    fasta_path = tmp_path / "synthetic.fna"
    with fasta_path.open("w", encoding="utf-8") as handle:
        handle.write(
            f">ctg_pk1\n{ctg1}\n>ctg_pk2_edge\n{ctg2}\n"
            f">ctg_track_b\n{ctg3}\n>ctg_fas\n{ctg4}\n>ctg_kas3\n{ctg5}\n"
        )

    pfam_table = det.load_pfam_table(PFAM_LIST)
    hits, stats = det.scan_fasta_to_hits(
        fasta_path, hmm_path, pfam_table, min_orf_aa=150, batch_size=50, cpus=2
    )
    assert stats["contigs"] == 5
    hotspots = det.cluster_track_a_hits(hits, d_max_bp=15000, min_domains=5)

    by_contig = {}
    for row in hotspots:
        by_contig.setdefault(row["contig"], []).append(row)

    assert "ctg_pk1" in by_contig
    assert len(by_contig["ctg_pk1"]) == 1  # intron bridged into single hotspot!
    assert by_contig["ctg_pk1"][0]["strand"] == "+"
    assert by_contig["ctg_pk1"][0]["completeness_label"] == "Internal"

    assert "ctg_pk2_edge" in by_contig
    assert by_contig["ctg_pk2_edge"][0]["strand"] == "-"
    assert by_contig["ctg_pk2_edge"][0]["completeness_label"] == "Edge-truncated"

    assert "ctg_track_b" not in by_contig  # Track B excluded from Track A clustering!
    assert "ctg_kas3" not in by_contig  # A/negative_control hits excluded from Track A clustering!

    kas3_hits = [h for h in hits if h.contig == "ctg_kas3"]
    assert len(kas3_hits) >= 2 and all(h.role == "negative_control" for h in kas3_hits)

    fas_hits = [h for h in hits if h.contig == "ctg_fas"]
    assert len(fas_hits) >= 5  # Baseline FAS domains are recorded for separate false-call reporting

    # Write outputs and verify recluster reproduces identical hotspots
    outdir = tmp_path / "out_scan"
    domains_tsv = outdir / "domains.tsv"
    det.write_domains_tsv(domains_tsv, hits)
    det.write_hotspots_outputs(outdir, hotspots)
    det.write_metadata_json(
        outdir / "metadata.json",
        mode="scan",
        parameters={"d_max_bp": 15000, "min_domains": 5},
        inputs={"fasta": det.file_metadata(fasta_path)},
        hits=hits,
        hotspots=hotspots,
        scan_stats=stats,
        elapsed_seconds=0.1,
    )
    reread_hits = det.read_domains_tsv(domains_tsv)
    reclustered = det.cluster_track_a_hits(reread_hits, d_max_bp=15000, min_domains=5)
    assert reclustered == hotspots

    meta = json.loads((outdir / "metadata.json").read_text(encoding="utf-8"))
    assert meta["summary"]["hotspots_total"] == len(hotspots)
    with (outdir / "hotspots.bed").open(encoding="utf-8") as handle:
        bed_rows = [line.strip().split("\t") for line in handle if line.strip()]
    assert len(bed_rows) == len(hotspots)
    with (outdir / "hotspots.tsv").open(encoding="utf-8") as handle:
        tsv_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(tsv_rows) == len(hotspots)

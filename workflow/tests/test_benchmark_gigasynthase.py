#!/usr/bin/env python3
"""Unit tests for PKZILLA truth coordinate parsing and Track A benchmark metrics."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import benchmark_gigasynthase as bm  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_truth_bed_and_tsv_separate_coordinate_parsing(tmp_path: Path):
    bed_path = tmp_path / "truth.bed"
    tsv_path = tmp_path / "truth.tsv"
    bed_path.write_text(
        "12B1-Scaf17\t2284380\t2427553\tPKZILLA-1__12B1-Scaf17:2284381-2427444:._+\t0\t+\n"
        "12B1-Scaf14\t378171\t415553\t12B1-Scaf14:378172-415553:._-\t0\t-\n",
        encoding="utf-8",
    )
    tsv_path.write_text(
        "hotspot_coord\tsize\n"
        "PKZILLA-1__12B1-Scaf17:2284381-2427444:._+::12B1-Scaf17:2284380-2427553(+)\t143173\n"
        "12B1-Scaf14:378172-415553:._-::12B1-Scaf14:378171-415553(-)\t37382\n",
        encoding="utf-8",
    )

    bed = bm.parse_truth_bed(bed_path)
    tsv = bm.parse_truth_tsv(tsv_path)
    bm.verify_bed_and_tsv_consistency(bed, tsv)

    assert tsv[0].original_label == "PKZILLA-1"
    assert (tsv[0].original_start_raw, tsv[0].original_end_raw) == (2284381, 2427444)
    assert (tsv[0].merged_start0, tsv[0].merged_end0) == (2284380, 2427553)
    assert (tsv[1].original_start_raw, tsv[1].original_end_raw) == (378172, 415553)
    assert (tsv[1].merged_start0, tsv[1].merged_end0) == (378171, 415553)

    # Mismatch between BED and TSV merged coordinates must fail loudly
    bad_tsv = tmp_path / "bad.tsv"
    bad_tsv.write_text(
        "hotspot_coord\tsize\n"
        "PKZILLA-1__12B1-Scaf17:2284381-2427444:._+::12B1-Scaf17:2284381-2427444(+)\t143063\n"
        "12B1-Scaf14:378172-415553:._-::12B1-Scaf14:378171-415553(-)\t37382\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="do not match"):
        bm.verify_bed_and_tsv_consistency(bed, bm.parse_truth_tsv(bad_tsv))


def test_parse_unmerged_pkzilla_loci_from_manifest(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sequence_id\tsource_identifier_or_coordinates\n"
        "PKZILLA-1\tPKZILLA-1;12B1-Scaf17:2286087-2423367(+)\n"
        "PKZILLA-2\tPKZILLA-2;12B1-Scaf7:1023621-1116873(+)\n"
        "PKZILLA-3\tPKZILLA-3;12B1-Scaf10:233715-308247(+)\n"
        "PKZILLA-B1\tPKZILLA-B1;rc_scaffold_1158-RCC3426:9886-127934(+)\n",
        encoding="utf-8",
    )
    loci = bm.parse_unmerged_pkzilla_loci(manifest)
    by_id = {locus.label: locus for locus in loci}
    assert set(by_id) == {"PKZILLA-1", "PKZILLA-2", "PKZILLA-3"}
    assert (by_id["PKZILLA-1"].contig, by_id["PKZILLA-1"].start0, by_id["PKZILLA-1"].end0, by_id["PKZILLA-1"].strand) == (
        "12B1-Scaf17", 2286086, 2423367, "+"
    )

    real_manifest = REPO_ROOT / "resources/seeds/manifest.tsv"
    if real_manifest.is_file():
        assert bm.parse_unmerged_pkzilla_loci(real_manifest) == loci


def test_benchmark_metrics_unmerged_fragmentation_and_outside_truth():
    unmerged = [
        bm.GenomicInterval("PKZILLA-1", "Scaf1", 1000, 11000, "+"),
        bm.GenomicInterval("PKZILLA-2", "Scaf2", 2000, 12000, "+"),
        bm.GenomicInterval("PKZILLA-3", "Scaf3", 5000, 15000, "+"),
    ]
    merged = [
        bm.GenomicInterval("PKZILLA-1__...", "Scaf1", 500, 11500, "+"),
        bm.GenomicInterval("PKZILLA-2__...", "Scaf2", 1500, 12500, "+"),
        bm.GenomicInterval("PKZILLA-3__...", "Scaf3", 4500, 15500, "+"),
    ]
    hotspots = [
        # PKZILLA-1 split into 2 components covering 9000 / 10000 bp (0.9)
        {"hotspot_id": "H0001", "contig": "Scaf1", "g_start": 1000, "g_end": 6000, "length_bp": 5000, "strand": "+", "n_domains": 10, "n_ks_ca_backbone": 2, "edge_truncated": "false", "n_neg_control_overlaps": 0},
        {"hotspot_id": "H0002", "contig": "Scaf1", "g_start": 7000, "g_end": 11000, "length_bp": 4000, "strand": "+", "n_domains": 8, "n_ks_ca_backbone": 1, "edge_truncated": "false", "n_neg_control_overlaps": 0},
        # PKZILLA-2 opposite strand only -> 0 same-strand coverage, 1 any-strand coverage
        {"hotspot_id": "H0003", "contig": "Scaf2", "g_start": 2000, "g_end": 12000, "length_bp": 10000, "strand": "-", "n_domains": 12, "n_ks_ca_backbone": 2, "edge_truncated": "false", "n_neg_control_overlaps": 0},
        # Novel call outside all truth
        {"hotspot_id": "H0004", "contig": "Scaf9", "g_start": 0, "g_end": 3000, "length_bp": 3000, "strand": "+", "n_domains": 6, "n_ks_ca_backbone": 1, "edge_truncated": "true", "n_neg_control_overlaps": 1},
    ]
    summary, locus_rows, _ = bm.evaluate_predictions(hotspots, unmerged, merged)
    assert summary["pkzilla_recall_same_strand"] == 1
    assert summary["pkzilla_recall_any_strand"] == 2
    by_loc = {r["locus_id"]: r for r in locus_rows}
    assert by_loc["PKZILLA-1"]["components_same_strand"] == 2
    assert by_loc["PKZILLA-1"]["unmerged_coverage_same_strand"] == 0.9
    assert by_loc["PKZILLA-2"]["detected_same_strand"] is False
    assert by_loc["PKZILLA-2"]["detected_any_strand"] is True
    assert summary["predicted_outside_truth_bp_any_strand"] == 3000
    assert summary["predicted_outside_truth_bp_same_strand"] == 13000
    assert summary["predicted_12b1_hotspots_with_neg_control_domains"] == 1

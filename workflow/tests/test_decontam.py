#!/usr/bin/env python3
"""Unit and integration tests for P3 Tiara + Whokaryote decontamination."""

from __future__ import annotations

import csv
import gzip
import json
import random
import subprocess
import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
DECONTAM_PY = Path(__file__).resolve().parents[2] / "tmp/env_smokes/decontam/bin/python"
REPO_ROOT = Path(__file__).resolve().parents[2]


def run_decontam_py(code: str) -> str:
    proc = subprocess.run(
        [str(DECONTAM_PY), "-c", code],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout


def test_python38_syntax_and_decontam_prerequisites():
    # Compile in Python 3.8 decontam env and verify tools/models validate
    subprocess.run(
        [
            str(DECONTAM_PY),
            "-m",
            "py_compile",
            "workflow/scripts/decontaminate_assemblies.py",
            "workflow/scripts/benchmark_decontam.py",
        ],
        cwd=str(REPO_ROOT),
        check=True,
    )
    out = run_decontam_py(
        "import sys; sys.path.insert(0, 'workflow/scripts'); "
        "import decontaminate_assemblies as d; "
        "tools = d.validate_decontam_prerequisites(); "
        "assert 'tiara' in tools and 'whokaryote.py' in tools and 'prodigal' in tools; "
        "assert 'RF_S_RBS_240121_bal.joblib' in tools and 'RF_T_RBS_310122_bal.joblib' in tools; "
        "print('PREREQS_OK')"
    )
    assert "PREREQS_OK" in out


def test_five_class_labeling_organelle_and_hgt_guards():
    code = """
import json, sys
sys.path.insert(0, 'workflow/scripts')
import decontaminate_assemblies as d

# 1. Organelle protection: even if Whokaryote_S says prokaryote (0.95), Tiara organelle stays organellar & retained
org_row = d.classify_contig(
    contig='pt1', length_bp=40000,
    tiara={'fst': 'organelle', 'snd': 'plastid', 'p_org': 0.92, 'p_bac': 0.05, 'p_arc': 0.01, 'p_euk': 0.02, 'p_unk1': 0.0, 'p_pla': 0.95, 'p_unk2': 0.0, 'p_mit': 0.05, 'p_prok': 0.06},
    whokaryote={'nr_genes': 30, 'gene_density': 0.85, 'rbs_ratio': 0.7, 'ratio_same_orientation': 0.8, 's_pred': 'prokaryote', 's_p_prok': 0.95, 't_pred': 'eukaryote', 't_p_prok': 0.1},
    algal_evidence=None,
)
assert org_row['final_label'] == 'organellar' and org_row['retained_in_filtered'] == 'true'

# 2. Disagreement alone must NEVER produce putative_hgt -> becomes ambiguous & retained
dis_row = d.classify_contig(
    contig='dis1', length_bp=4000,
    tiara={'fst': 'bacteria', 'snd': 'n/a', 'p_org': 0.01, 'p_bac': 0.85, 'p_arc': 0.02, 'p_euk': 0.12, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.87},
    whokaryote={'nr_genes': 3, 'gene_density': 0.5, 'rbs_ratio': 0.2, 'ratio_same_orientation': 0.5, 's_pred': 'eukaryote', 's_p_prok': 0.10, 't_pred': 'prokaryote', 't_p_prok': 0.6},
    algal_evidence=None,
)
assert dis_row['final_label'] == 'ambiguous' and dis_row['retained_in_filtered'] == 'true'

# 3. Clear bacterial consensus without algal evidence -> bacterial & excluded
bac_row = d.classify_contig(
    contig='bac1', length_bp=12000,
    tiara={'fst': 'bacteria', 'snd': 'n/a', 'p_org': 0.0, 'p_bac': 0.94, 'p_arc': 0.02, 'p_euk': 0.04, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.96},
    whokaryote={'nr_genes': 10, 'gene_density': 0.88, 'rbs_ratio': 0.75, 'ratio_same_orientation': 0.8, 's_pred': 'prokaryote', 's_p_prok': 0.92, 't_pred': 'prokaryote', 't_p_prok': 0.95},
    algal_evidence=None,
)
assert bac_row['final_label'] == 'bacterial' and bac_row['retained_in_filtered'] == 'false'

# 4. Same bacterial signal WITH explicit multi-exon algal evidence -> putative_hgt & rescued
hgt_row = d.classify_contig(
    contig='hgt1', length_bp=12000,
    tiara={'fst': 'bacteria', 'snd': 'n/a', 'p_org': 0.0, 'p_bac': 0.94, 'p_arc': 0.02, 'p_euk': 0.04, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.96},
    whokaryote={'nr_genes': 10, 'gene_density': 0.88, 'rbs_ratio': 0.75, 'ratio_same_orientation': 0.8, 's_pred': 'prokaryote', 's_p_prok': 0.92, 't_pred': 'prokaryote', 't_p_prok': 0.95},
    algal_evidence='gff3_multi_exon_transcripts=2',
)
assert hgt_row['final_label'] == 'putative_hgt' and hgt_row['retained_in_filtered'] == 'true'

# 5. Clear eukaryotic consensus -> eukaryotic_haptophyte & retained
euk_row = d.classify_contig(
    contig='euk1', length_bp=20000,
    tiara={'fst': 'eukarya', 'snd': 'n/a', 'p_org': 0.0, 'p_bac': 0.01, 'p_arc': 0.0, 'p_euk': 0.99, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.01},
    whokaryote={'nr_genes': 4, 'gene_density': 0.25, 'rbs_ratio': 0.1, 'ratio_same_orientation': 0.5, 's_pred': 'eukaryote', 's_p_prok': 0.05, 't_pred': 'eukaryote', 't_p_prok': 0.02},
    algal_evidence=None,
)
assert euk_row['final_label'] == 'eukaryotic_haptophyte' and euk_row['retained_in_filtered'] == 'true'

# 6. Below Tiara min length -> ambiguous & retained
tiny_row = d.classify_contig(
    contig='tiny1', length_bp=300,
    tiara=None, whokaryote=None, algal_evidence=None,
)
assert tiny_row['final_label'] == 'ambiguous' and tiny_row['retained_in_filtered'] == 'true'

# 7. Length-calibrated <1kb protection: p_tiara_prok=0.85 stays ambiguous (<0.95), p_tiara_prok=0.97 -> bacterial
lt1k_amb = d.classify_contig(
    contig='short_amb', length_bp=700,
    tiara={'fst': 'bacteria', 'snd': 'n/a', 'p_org': 0.0, 'p_bac': 0.85, 'p_arc': 0.0, 'p_euk': 0.15, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.85},
    whokaryote=None, algal_evidence=None,
)
assert lt1k_amb['final_label'] == 'ambiguous' and float(lt1k_amb['ensemble_prok_score']) < 0.80
lt1k_bac = d.classify_contig(
    contig='short_bac', length_bp=700,
    tiara={'fst': 'bacteria', 'snd': 'n/a', 'p_org': 0.0, 'p_bac': 0.97, 'p_arc': 0.0, 'p_euk': 0.03, 'p_unk1': 0.0, 'p_pla': None, 'p_unk2': None, 'p_mit': None, 'p_prok': 0.97},
    whokaryote=None, algal_evidence=None,
)
assert lt1k_bac['final_label'] == 'bacterial' and lt1k_bac['retained_in_filtered'] == 'false'
print('LABELS_OK')
"""
    assert "LABELS_OK" in run_decontam_py(code)


def test_end_to_end_gzipped_decontam_and_gff_evidence(tmp_path: Path):
    rng = random.Random(123)
    in_gz = tmp_path / "sample.fasta.gz"
    gff3 = tmp_path / "annot.gff3"
    outdir = tmp_path / "out"

    with gzip.open(in_gz, "wt", encoding="utf-8") as h:
        for name, length in [("ctg_short_300", 300), ("ctg_1500", 1500), ("ctg_6000", 6000)]:
            seq = "".join(rng.choices("ACGT", k=length))
            h.write(f">{name} desc=test\n{seq}\n")
    gff3.write_text(
        "ctg_6000\tTest\tmRNA\t1\t2000\t.\t+\t.\tID=rna1\n"
        "ctg_6000\tTest\tCDS\t10\t500\t.\t+\t0\tID=cds1;Parent=rna1\n"
        "ctg_6000\tTest\tCDS\t700\t1200\t.\t+\t0\tID=cds2;Parent=rna1\n",
        encoding="utf-8",
    )
    orig_bytes = in_gz.read_bytes()

    subprocess.run(
        [
            str(DECONTAM_PY),
            "workflow/scripts/decontaminate_assemblies.py",
            "run",
            "--fasta",
            str(in_gz),
            "--sample-id",
            "syn_test",
            "--outdir",
            str(outdir),
            "--gff3",
            str(gff3),
            "--tiara-min-len",
            "500",
            "--whokaryote-min-len",
            "1000",
            "--threads",
            "2",
        ],
        cwd=str(REPO_ROOT),
        check=True,
    )

    assert in_gz.read_bytes() == orig_bytes  # Source assembly untouched
    assert (outdir / "raw/tiara_pred.txt").is_file()
    assert (outdir / "raw/whokaryote_featuretable.csv").is_file()
    assert (outdir / "filtered.fasta").is_file()
    assert (outdir / "excluded_ids.txt").is_file()

    with (outdir / "contig_evidence.tsv").open(newline="", encoding="utf-8") as h:
        ev = {r["contig"]: r for r in csv.DictReader(h, delimiter="\t")}
    assert len(ev) == 3
    assert ev["ctg_short_300"]["final_label"] == "ambiguous"
    assert ev["ctg_6000"]["has_explicit_algal_evidence"] == "true"

    summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
    assert summary["summary"]["total_contigs"] == 3


def test_benchmark_metric_functions():
    code = """
import sys, numpy as np
sys.path.insert(0, 'workflow/scripts')
import benchmark_decontam as bm

yt = np.array([0, 0, 0, 0, 1, 1, 1, 1])
ys = np.array([0.1, 0.2, 0.3, 0.8, 0.5, 0.6, 0.85, 0.95])
assert bm.recall_at_fixed_fpr(yt, ys, 0.01) == 0.5  # > 0.8 -> 2/4
assert bm.recall_at_fixed_fpr(yt, ys, 0.25) == 1.0  # > 0.3 -> 4/4
print('METRICS_OK')
"""
    assert "METRICS_OK" in run_decontam_py(code)

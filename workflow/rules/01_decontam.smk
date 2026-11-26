import csv
from pathlib import Path


def _load_samples_meta(manifest_path):
    mapping = {}
    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["sample_id"]] = row
    return mapping


SAMPLES_META = _load_samples_meta(config["samples_tsv"])
ALL_SAMPLE_IDS = list(SAMPLES_META.keys())
DECONTAM_CFG = config.get("decontam", {})
TIARA_MIN_LEN = int(DECONTAM_CFG.get("tiara_min_len", 500))
WHOKARYOTE_MIN_LEN = int(DECONTAM_CFG.get("whokaryote_min_len", 1000))
DECONTAM_PY = config.get("decontam_python", "python")


rule decontam_sample:
    input:
        fasta=lambda wc: SAMPLES_META[wc.sample_id]["genome_fasta"],
        gff3=lambda wc: (
            SAMPLES_META[wc.sample_id]["gff3"]
            if SAMPLES_META[wc.sample_id]["gff3"] != "NA"
            else []
        ),
        script="workflow/scripts/decontaminate_assemblies.py",
    output:
        filtered="results/decontam/{sample_id}/filtered.fasta",
        excluded="results/decontam/{sample_id}/excluded_ids.txt",
        evidence="results/decontam/{sample_id}/contig_evidence.tsv",
        summary_json="results/decontam/{sample_id}/summary.json",
        summary_tsv="results/decontam/{sample_id}/summary.tsv",
        tiara_raw="results/decontam/{sample_id}/raw/tiara_pred.txt",
        whokaryote_feat="results/decontam/{sample_id}/raw/whokaryote_featuretable.csv",
    params:
        py=DECONTAM_PY,
        outdir="results/decontam/{sample_id}",
        gff3_arg=lambda wc: SAMPLES_META[wc.sample_id]["gff3"],
        tiara_min_len=TIARA_MIN_LEN,
        whokaryote_min_len=WHOKARYOTE_MIN_LEN,
    threads: 4
    conda:
        "../../envs/decontam.yaml"
    shell:
        r"""
        {params.py} {input.script} run \
            --fasta {input.fasta} \
            --sample-id {wildcards.sample_id} \
            --outdir {params.outdir} \
            --gff3 {params.gff3_arg} \
            --tiara-min-len {params.tiara_min_len} \
            --whokaryote-min-len {params.whokaryote_min_len} \
            --threads {threads}
        """


rule decontam_all_summary:
    input:
        summaries=expand("results/decontam/{sample_id}/summary.tsv", sample_id=ALL_SAMPLE_IDS),
        script="workflow/scripts/decontaminate_assemblies.py",
    output:
        tsv="results/decontam/all_samples_summary.tsv",
        json="results/decontam/all_samples_summary.json",
    params:
        py=DECONTAM_PY,
    conda:
        "../../envs/decontam.yaml"
    shell:
        r"""
        {params.py} {input.script} aggregate \
            --summaries {input.summaries} \
            --output-tsv {output.tsv} \
            --output-json {output.json}
        """


rule benchmark_decontam:
    input:
        euk_fasta=SAMPLES_META["pparvum_12B1_v1_1"]["genome_fasta"],
        refs_manifest=DECONTAM_CFG.get("references_manifest", "resources/decontam_refs/manifest.tsv"),
        organelle_dir=DECONTAM_CFG.get(
            "organelle_dir", "data/genomes/phaeocystis_osf_vka93/organellar genomes"
        ),
        script="workflow/scripts/benchmark_decontam.py",
        core_script="workflow/scripts/decontaminate_assemblies.py",
    output:
        summary="results/benchmarks/decontam/benchmark_summary.json",
        spikein="results/benchmarks/decontam/spikein_metrics.tsv",
        stratified="results/benchmarks/decontam/stratified_metrics.tsv",
        organelle="results/benchmarks/decontam/organelle_check.tsv",
    params:
        py=DECONTAM_PY,
        outdir="results/benchmarks/decontam",
        tiara_min_len=TIARA_MIN_LEN,
        whokaryote_min_len=WHOKARYOTE_MIN_LEN,
    threads: 4
    conda:
        "../../envs/decontam.yaml"
    shell:
        r"""
        {params.py} {input.script} \
            --euk-fasta {input.euk_fasta} \
            --refs-manifest {input.refs_manifest} \
            --organelle-dir "{input.organelle_dir}" \
            --outdir {params.outdir} \
            --tiara-min-len {params.tiara_min_len} \
            --whokaryote-min-len {params.whokaryote_min_len} \
            --threads {threads}
        """


rule track_a_detect_decontam:
    input:
        fasta="results/decontam/{sample_id}/filtered.fasta",
        hmm=config["databases"]["custom_hmm"],
        pressed=multiext(config["databases"]["custom_hmm"], ".h3f", ".h3i", ".h3m", ".h3p"),
        manifest=config["databases"]["custom_hmm_manifest"],
        pfam_list="resources/pfam_lists/hapto_core_pfam.tsv",
        script="workflow/scripts/detect_gigasynthase_hotspots.py",
    output:
        domains="results/track_a_decontam/{sample_id}/domains.tsv",
        hotspots="results/track_a_decontam/{sample_id}/hotspots.tsv",
        bed="results/track_a_decontam/{sample_id}/hotspots.bed",
        metadata="results/track_a_decontam/{sample_id}/metadata.json",
    params:
        outdir="results/track_a_decontam/{sample_id}",
        min_orf_aa=150,
        d_max=15000,
        min_domains=5,
    threads: config.get("threads", 16)
    conda:
        "../../envs/core.yaml"
    shell:
        r"""
        python {input.script} scan \
            --fasta {input.fasta} \
            --hmm {input.hmm} \
            --pfam-list {input.pfam_list} \
            --outdir {params.outdir} \
            --min-orf-aa {params.min_orf_aa} \
            --d-max {params.d_max} \
            --min-domains {params.min_domains} \
            --cpus {threads} \
            --unfiltered-warning ""
        """

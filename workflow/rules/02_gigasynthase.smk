import csv
from pathlib import Path


def _load_sample_fastas(manifest_path):
    mapping = {}
    with Path(manifest_path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["sample_id"]] = row["genome_fasta"]
    return mapping


SAMPLE_FASTAS = _load_sample_fastas(config["samples_tsv"])
CUSTOM_HMM = config["databases"]["custom_hmm"]
CUSTOM_HMM_MANIFEST = config["databases"]["custom_hmm_manifest"]
PFAM_LIST = "resources/pfam_lists/hapto_core_pfam.tsv"


rule track_a_detect:
    input:
        fasta=lambda wc: SAMPLE_FASTAS[wc.sample_id],
        hmm=CUSTOM_HMM,
        pressed=multiext(CUSTOM_HMM, ".h3f", ".h3i", ".h3m", ".h3p"),
        manifest=CUSTOM_HMM_MANIFEST,
        pfam_list=PFAM_LIST,
        script="workflow/scripts/detect_gigasynthase_hotspots.py",
    output:
        domains="results/track_a/{sample_id}/domains.tsv",
        hotspots="results/track_a/{sample_id}/hotspots.tsv",
        bed="results/track_a/{sample_id}/hotspots.bed",
        metadata="results/track_a/{sample_id}/metadata.json",
    params:
        outdir="results/track_a/{sample_id}",
        min_orf_aa=150,
        d_max=15000,
        min_domains=5,
        unfiltered_warning=lambda wc: (
            "UNFILTERED_EXPLORATORY_SCAN: Assembly has not been decontaminated; do not interpret biologically."
            if wc.sample_id == "pglobosa_CNS00066"
            else ""
        ),
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
            --unfiltered-warning "{params.unfiltered_warning}"
        """


rule benchmark_track_a_12b1:
    input:
        hotspots="results/track_a/pparvum_12B1_v1_1/hotspots.tsv",
        truth_bed=config["benchmarks"]["pkzilla_hotspots_bed"],
        truth_tsv="data/benchmarks/pkzilla_truth/zenodo_10309063_hotspot_liftover/2023-12-08_Hotspot_liftover/v1.1_liftover_PKS_hotspots_5000.merge.tsv",
        seeds_manifest="resources/seeds/manifest.tsv",
        pufa_faa="resources/seeds/pufa_neg_controls.faa",
        pkzilla_faa="resources/seeds/pkzilla_seeds.faa",
        pkzilla_domains_dir=config["benchmarks"]["pkzilla_domains_dir"],
        hmm=CUSTOM_HMM,
        pfam_list=PFAM_LIST,
        ehuxleyi_fasta=SAMPLE_FASTAS["ehuxleyi_CCMP1516"],
        script="workflow/scripts/benchmark_gigasynthase.py",
    output:
        summary="results/benchmarks/track_a_12B1/benchmark_summary.json",
        loci="results/benchmarks/track_a_12B1/locus_metrics.tsv",
        merged="results/benchmarks/track_a_12B1/merged_bed_metrics.tsv",
        controls="results/benchmarks/track_a_12B1/control_eval.tsv",
        protein_check="results/benchmarks/track_a_12B1/protein_domain_recall.tsv",
    params:
        outdir="results/benchmarks/track_a_12B1",
    threads: 2
    conda:
        "../../envs/core.yaml"
    shell:
        r"""
        python {input.script} evaluate \
            --hotspots-tsv {input.hotspots} \
            --truth-bed {input.truth_bed} \
            --truth-tsv {input.truth_tsv} \
            --seeds-manifest {input.seeds_manifest} \
            --pufa-faa {input.pufa_faa} \
            --pkzilla-faa {input.pkzilla_faa} \
            --pkzilla-domains-dir {input.pkzilla_domains_dir} \
            --hmm {input.hmm} \
            --pfam-list {input.pfam_list} \
            --ehuxleyi-fasta {input.ehuxleyi_fasta} \
            --outdir {params.outdir} \
            --cpus {threads}
        """


rule sensitivity_track_a_12b1:
    input:
        domains="results/track_a/pparvum_12B1_v1_1/domains.tsv",
        truth_bed=config["benchmarks"]["pkzilla_hotspots_bed"],
        truth_tsv="data/benchmarks/pkzilla_truth/zenodo_10309063_hotspot_liftover/2023-12-08_Hotspot_liftover/v1.1_liftover_PKS_hotspots_5000.merge.tsv",
        seeds_manifest="resources/seeds/manifest.tsv",
        script="workflow/scripts/benchmark_gigasynthase.py",
    output:
        sensitivity="results/benchmarks/track_a_12B1/sensitivity.tsv",
    conda:
        "../../envs/core.yaml"
    shell:
        r"""
        python {input.script} sensitivity \
            --domains-tsv {input.domains} \
            --truth-bed {input.truth_bed} \
            --truth-tsv {input.truth_tsv} \
            --seeds-manifest {input.seeds_manifest} \
            --output {output.sensitivity}
        """

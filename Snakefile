configfile: "config/config.yaml"


include: "workflow/rules/00_prep_db.smk"
include: "workflow/rules/01_decontam.smk"
include: "workflow/rules/02_gigasynthase.smk"


rule all:
    input:
        config["databases"]["custom_hmm_manifest"],
        "tmp/hmm_smoke/scan.calls.tsv",
        "results/track_a/pparvum_12B1_v1_1/hotspots.tsv",
        "results/benchmarks/track_a_12B1/benchmark_summary.json",
        "results/benchmarks/decontam/benchmark_summary.json",
        "results/decontam/all_samples_summary.tsv",
        "results/track_a_decontam/pglobosa_CNS00066/hotspots.tsv",

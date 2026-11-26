PFAM_SELECTION = "resources/pfam_lists/hapto_core_pfam.tsv"
CUSTOM_HMM = config["databases"]["custom_hmm"]
CUSTOM_HMM_MANIFEST = config["databases"]["custom_hmm_manifest"]


rule prepare_hapto_core:
    input:
        selection=PFAM_SELECTION,
        pfam_hmm=config["databases"]["pfam_hmm"],
        pfam_dat=config["databases"]["pfam_dat"],
        builder="workflow/scripts/build_hmm_database.py",
    output:
        hmm=CUSTOM_HMM,
        pressed=multiext(CUSTOM_HMM, ".h3f", ".h3i", ".h3m", ".h3p"),
        manifest=CUSTOM_HMM_MANIFEST,
    params:
        release=config["pfam_release"],
        expected_count=config["pfam_expected_count"],
    conda:
        "../../envs/core.yaml"
    shell:
        r"""
        python {input.builder} extract \
            --selection {input.selection} \
            --pfam-dat {input.pfam_dat} \
            --pfam-hmm {input.pfam_hmm} \
            --pfam-release {params.release} \
            --expected-pfam-count {params.expected_count} \
            --output {output.hmm}
        hmmpress -f {output.hmm}
        python {input.builder} manifest \
            --selection {input.selection} \
            --pfam-dat {input.pfam_dat} \
            --pfam-hmm {input.pfam_hmm} \
            --pfam-release {params.release} \
            --expected-pfam-count {params.expected_count} \
            --output-hmm {output.hmm} \
            --pressed {output.pressed} \
            --output-manifest {output.manifest}
        """


rule hmm_scan_smoke:
    input:
        fixture="resources/test_data/hmm_smoke.faa",
        expected="resources/test_data/hmm_smoke_expected.json",
        checker="workflow/scripts/check_hmm_smoke.py",
        hmm=CUSTOM_HMM,
        pressed=multiext(CUSTOM_HMM, ".h3f", ".h3i", ".h3m", ".h3p"),
        manifest=CUSTOM_HMM_MANIFEST,
    output:
        domtbl="tmp/hmm_smoke/scan.domtbl",
        calls="tmp/hmm_smoke/scan.calls.tsv",
    log:
        "tmp/hmm_smoke/hmmscan.log",
    conda:
        "../../envs/core.yaml"
    shell:
        r"""
        hmmscan --cut_ga --domtblout {output.domtbl} {input.hmm} {input.fixture} > {log}
        python {input.checker} validate \
            --fixture {input.fixture} \
            --domtbl {output.domtbl} \
            --expected {input.expected} \
            --calls {output.calls}
        """

#!/usr/bin/env python3
"""Build the small, provenance-tracked seed collections used by HaptoSMASH."""

import argparse
import csv
import hashlib
import io
import time
import urllib.request
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


FIELDS = [
    "sequence_id",
    "gene_family",
    "species_strain",
    "role",
    "source_file_or_accession",
    "source_identifier_or_coordinates",
    "length_aa",
    "sequence_sha256",
    "source_sha256",
    "evidence_notes",
]


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_index(path):
    with path.open("rt") as handle:
        records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def local_spec(output, sequence_id, family, species, role, path, source_id, evidence):
    return {
        "output": output,
        "sequence_id": sequence_id,
        "gene_family": family,
        "species_strain": species,
        "role": role,
        "path": path,
        "source_id": source_id,
        "evidence_notes": evidence,
    }


def remote_spec(output, sequence_id, family, species, role, accession, evidence):
    return local_spec(
        output, sequence_id, family, species, role, None, accession, evidence
    ) | {"accession": accession}


def fetch_uniprot(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = urllib.request.Request(url, headers={"User-Agent": "HaptoSMASH-seed-curation/1.0"})
    last_error = None
    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
            records = list(SeqIO.parse(io.StringIO(payload.decode("utf-8")), "fasta"))
            if len(records) != 1 or accession not in records[0].description:
                raise ValueError(f"Unexpected UniProt FASTA response for {accession}")
            return records[0], payload, url
        except Exception as error:  # Network errors are reported after bounded retries.
            last_error = error
            if attempt < 3:
                time.sleep(2 * attempt)
    raise RuntimeError(f"Unable to retrieve {url}: {last_error}")


def write_fasta(path, records):
    with path.open("wt") as handle:
        for record, row in records:
            handle.write(
                f">{row['sequence_id']} family={row['gene_family']} "
                f"role={row['role']} length_aa={len(record.seq)}\n"
            )
            sequence = str(record.seq)
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = root / "resources/seeds"
    data = root / "data"
    pk_dir = data / "benchmarks/pkzilla_truth/zenodo_11068523_pks_domain_to_module"
    cns = data / "genomes/phaeocystis_globosa_CNS00066/proteins.faa"
    emi = data / "genomes/phaeocystis_osf_vka93"
    emi_proteins = data / "genomes/emiliania_huxleyi_CCMP1516/proteins.faa"
    p12 = pk_dir / "12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa"
    pb1 = pk_dir / "derived_files/renamed_PKZILLA-1-on-B.loci.nt.rename.consensus.fa.pep.fa"
    osf_dsyb = emi / "code/gene_mapping/parse_gene_mapping/dsyB_coords.bed"
    osf_genes = emi / "code/gene_mapping/Phaglo1_manual.bed"

    osf_dsyb_text = osf_dsyb.read_text()
    osf_gene_text = osf_genes.read_text()
    if "Phaglo1_scaffold_3\t69032\t71993\tPhaglo1_dsyB" not in osf_dsyb_text:
        raise ValueError("Expected OSF Phaglo1 DSYB candidate coordinates were not found")
    if "Phaglo1_scaffold_3\t69032\t72230\tPhaglo1_16308" not in osf_gene_text:
        raise ValueError("Expected OSF Phaglo1_16308 gene model was not found")

    pk_evidence = (
        "Full-length source record; protein-relative domain files were not concatenated. "
        "Fallon et al., Science 2024, DOI:10.1126/science.ado3290."
    )
    alma_evidence = (
        "CNS00066 source identifier verified in proteins.faa and annotation.gff3; "
        "OSF Alma1_filtered.bed supplies independent Phaeocystis homolog coordinates; "
        "Alcolombri et al., Science 2015, DOI:10.1126/science.aab1586."
    )
    dsyb_evidence = (
        "Functionally characterized bacterial DSYB reference (Curson et al., Nature Microbiology "
        "2018, DOI:10.1038/s41564-018-0119-5). OSF dsyB_coords.bed identifies "
        "Phaglo1_scaffold_3:69032-71993 overlapping Phaglo1_16308, but that candidate is not "
        "extracted because the available BED is a gene span rather than verified protein evidence."
    )
    specs = [
        local_spec("pkzilla_seeds.faa", "PKZILLA-1", "PKZILLA-1", "Prymnesium parvum 12B1", "positive", p12, "PKZILLA-1;12B1-Scaf17:2286087-2423367(+)", pk_evidence),
        local_spec("pkzilla_seeds.faa", "PKZILLA-2", "PKZILLA-2", "Prymnesium parvum 12B1", "positive", p12, "PKZILLA-2;12B1-Scaf7:1023621-1116873(+)", pk_evidence),
        local_spec("pkzilla_seeds.faa", "PKZILLA-3", "PKZILLA-3", "Prymnesium parvum 12B1", "positive", p12, "PKZILLA-3;12B1-Scaf10:233715-308247(+)", pk_evidence),
        local_spec("pkzilla_seeds.faa", "PKZILLA-B1", "PKZILLA-B1", "Prymnesium parvum RCC3426", "positive", pb1, "PKZILLA-B1;rc_scaffold_1158-RCC3426_pilon_152285-290527cns_:9886-127934(+)", pk_evidence),
        local_spec("dmsp_seeds.faa", "GWHPBJCL017043", "PgMT2", "Phaeocystis globosa CNS00066", "positive", cns, "GWHPBJCL017043;Pgl17043;GWHBJCL00000058:2559386-2560934(+)", "CNS00066 PgMT2 mapping verified from local gene/protein identifiers; DMSP biosynthesis gene-family evidence from PMID:39193189."),
        local_spec("dmsp_seeds.faa", "GWHPBJCL009637", "Alma1", "Phaeocystis globosa CNS00066", "positive", cns, "GWHPBJCL009637;Pgl09637;GWHBJCL00000046:3086654-3088037(+)", alma_evidence),
        local_spec("dmsp_seeds.faa", "GWHPBJCL008660", "Alma2", "Phaeocystis globosa CNS00066", "positive", cns, "GWHPBJCL008660;Pgl08660;GWHBJCL00000047:6675615-6676559(+)", alma_evidence),
        local_spec("dmsp_seeds.faa", "GWHPBJCL017348", "Alma3", "Phaeocystis globosa CNS00066", "positive", cns, "GWHPBJCL017348;Pgl17348;GWHBJCL00000058:3709013-3710450(+)", alma_evidence),
        remote_spec("dmsp_seeds.faa", "DSYB_ROSEA", "DsyB", "Roseibium aggregatum", "positive", "A0A1P8C677", dsyb_evidence),
        remote_spec("dmsp_seeds.faa", "ALMA1_EMIHU", "Alma1", "Emiliania huxleyi", "positive", "P0DN21", "Reviewed UniProt P0DN21; experimentally characterized DMSP lyase (Alcolombri et al., Science 2015, DOI:10.1126/science.aab1586)."),
        remote_spec("pufa_neg_controls.faa", "PFAA_SHEON", "PUFA synthase PfaA", "Shewanella oneidensis MR-1", "negative_control", "Q8EGK0", "UniProt Q8EGK0; non-prymnesin PUFA synthase control for Track A backbone false positives."),
        local_spec("pufa_neg_controls.faa", "XP_005790240.1", "Fatty acid synthase", "Emiliania huxleyi CCMP1516", "negative_control", emi_proteins, "XP_005790240.1;NW_005201381.1:249218-261374(+)", "RefSeq protein and product annotation verified in local proteins.faa and annotation.gff3; FAS negative control."),
        local_spec("pufa_neg_controls.faa", "XP_005793321.1", "Beta-ketoacyl-ACP synthase", "Emiliania huxleyi CCMP1516", "negative_control", emi_proteins, "XP_005793321.1;NW_005202303.1:305615-306990(-)", "RefSeq protein and product annotation verified in local proteins.faa and annotation.gff3; KAS/FAS negative control."),
        remote_spec("halogenase_oxylipin_seeds.faa", "PRNA_PSEFL", "Tryptophan 7-halogenase PrnA", "Pseudomonas fluorescens", "positive", "P95480", "Reviewed UniProt P95480; experimentally characterized flavin-dependent tryptophan halogenase."),
        remote_spec("halogenase_oxylipin_seeds.faa", "LOX_EMIHU", "Arachidonate 15-lipoxygenase-like", "Emiliania huxleyi", "positive", "Q0MYU5", "UniProt Q0MYU5; haptophyte lipoxygenase reference for oxylipin-related Track B profiles."),
    ]

    indexes = {}
    grouped = {}
    manifest = []
    seen_ids = set()
    seen_sequences = {}
    for spec in specs:
        sequence_id = spec["sequence_id"]
        if sequence_id in seen_ids:
            raise ValueError(f"Duplicate output sequence ID: {sequence_id}")
        seen_ids.add(sequence_id)
        if spec["path"] is not None:
            path = spec["path"]
            if path not in indexes:
                indexes[path] = read_index(path)
            source_key = spec["source_id"].split(";", 1)[0]
            if source_key not in indexes[path]:
                raise ValueError(f"Source identifier {source_key} not found in {path}")
            record = indexes[path][source_key]
            source_name = str(path.relative_to(root))
            source_hash = sha256_file(path)
        else:
            record, payload, source_name = fetch_uniprot(spec["accession"])
            source_hash = sha256_bytes(payload)
        sequence = str(record.seq).upper()
        if sequence.endswith("*"):
            sequence = sequence[:-1]
            record.seq = Seq(sequence)
        if not sequence or "*" in sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise ValueError(f"Invalid or incomplete protein sequence for {sequence_id}")
        sequence_hash = sha256_bytes(sequence.encode("ascii"))
        if sequence_hash in seen_sequences:
            raise ValueError(f"Duplicate sequence: {sequence_id} and {seen_sequences[sequence_hash]}")
        seen_sequences[sequence_hash] = sequence_id
        row = {key: spec[key] for key in ["sequence_id", "gene_family", "species_strain", "role", "evidence_notes"]}
        row.update({
            "source_file_or_accession": source_name,
            "source_identifier_or_coordinates": spec["source_id"],
            "length_aa": len(sequence),
            "sequence_sha256": sequence_hash,
            "source_sha256": source_hash,
        })
        record.id = sequence_id
        record.name = sequence_id
        grouped.setdefault(spec["output"], []).append((record, row))
        manifest.append(row)
        print(f"selected\t{sequence_id}\t{len(sequence)} aa\t{source_name}")

    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, records in grouped.items():
        write_fasta(output_dir / filename, records)
    with (output_dir / "manifest.tsv").open("wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)

    for filename in ["pkzilla_seeds.faa", "dmsp_seeds.faa", "pufa_neg_controls.faa", "halogenase_oxylipin_seeds.faa"]:
        if filename not in grouped:
            raise ValueError(f"Required output was not generated: {filename}")
    print(f"Validated {len(manifest)} unique full-length seed records in 4 collections")
    print("Phaglo1_16308 remains an OSF DSYB candidate, not an extracted full-length seed")


if __name__ == "__main__":
    main()

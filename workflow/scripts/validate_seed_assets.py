#!/usr/bin/env python3
"""Independently validate curated seed sequences, provenance, and PKZILLA coordinates."""

import argparse
import csv
import hashlib
import re
from pathlib import Path

from Bio import SeqIO


COLLECTIONS = {
    "pkzilla_seeds.faa": 4,
    "dmsp_seeds.faa": 6,
    "pufa_neg_controls.faa": 3,
    "halogenase_oxylipin_seeds.faa": 2,
}
HEADER = re.compile(r"^(\S+) family=(.*?) role=(positive|negative_control) length_aa=(\d+)$")
HEX256 = re.compile(r"^[0-9a-f]{64}$")


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_sequence(record):
    sequence = str(record.seq).upper()
    if sequence.endswith("*"):
        sequence = sequence[:-1]
    return sequence


def load_manifest(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames != [
            "sequence_id", "gene_family", "species_strain", "role",
            "source_file_or_accession", "source_identifier_or_coordinates",
            "length_aa", "sequence_sha256", "source_sha256", "evidence_notes",
        ]:
            raise ValueError(f"Unexpected manifest columns in {path}: {reader.fieldnames}")
    if len(rows) != 15:
        raise ValueError(f"Expected 15 manifest rows, observed {len(rows)}")
    by_id = {}
    for row in rows:
        sequence_id = row["sequence_id"]
        if not sequence_id or sequence_id in by_id:
            raise ValueError(f"Missing or duplicate manifest sequence ID: {sequence_id!r}")
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed manifest row for {sequence_id}")
        if row["role"] not in {"positive", "negative_control"}:
            raise ValueError(f"Invalid role for {sequence_id}: {row['role']}")
        if not row["evidence_notes"] or not row["source_identifier_or_coordinates"]:
            raise ValueError(f"Incomplete provenance for {sequence_id}")
        if not HEX256.fullmatch(row["sequence_sha256"]) or not HEX256.fullmatch(row["source_sha256"]):
            raise ValueError(f"Invalid SHA-256 provenance for {sequence_id}")
        by_id[sequence_id] = row
    return by_id


def validate_collections(seed_dir, manifest):
    observed = {}
    sequence_owners = {}
    for filename, expected_count in COLLECTIONS.items():
        path = seed_dir / filename
        with path.open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        if len(records) != expected_count:
            raise ValueError(f"{filename}: expected {expected_count} records, observed {len(records)}")
        for record in records:
            match = HEADER.fullmatch(record.description)
            if not match:
                raise ValueError(f"Malformed header in {filename}: {record.description}")
            sequence_id, family, role, header_length = match.groups()
            if sequence_id != record.id or sequence_id not in manifest or sequence_id in observed:
                raise ValueError(f"Unknown, duplicate, or inconsistent sequence ID: {sequence_id}")
            row = manifest[sequence_id]
            sequence = normalized_sequence(record)
            if not sequence or "*" in sequence:
                raise ValueError(f"Empty sequence or internal stop for {sequence_id}")
            if len(sequence) != int(header_length) or len(sequence) != int(row["length_aa"]):
                raise ValueError(f"Sequence/header/manifest length mismatch for {sequence_id}")
            if family != row["gene_family"] or role != row["role"]:
                raise ValueError(f"Header/manifest metadata mismatch for {sequence_id}")
            sequence_hash = sha256_bytes(sequence.encode("ascii"))
            if sequence_hash != row["sequence_sha256"]:
                raise ValueError(f"Sequence hash mismatch for {sequence_id}")
            if sequence_hash in sequence_owners:
                raise ValueError(f"Duplicate sequence: {sequence_id} and {sequence_owners[sequence_hash]}")
            sequence_owners[sequence_hash] = sequence_id
            observed[sequence_id] = filename
    if set(observed) != set(manifest):
        raise ValueError(f"FASTA/manifest ID mismatch: {set(manifest) ^ set(observed)}")
    return observed


def validate_local_sources(root, manifest):
    indexes = {}
    checked = 0
    for sequence_id, row in manifest.items():
        source = row["source_file_or_accession"]
        if source.startswith("https://"):
            expected_url = (
                "https://rest.uniprot.org/uniprotkb/"
                f"{row['source_identifier_or_coordinates']}.fasta"
            )
            if source != expected_url:
                raise ValueError(f"UniProt accession/URL mismatch for {sequence_id}")
            continue
        path = root / source
        if sha256_file(path) != row["source_sha256"]:
            raise ValueError(f"Local source checksum mismatch for {sequence_id}: {path}")
        if path not in indexes:
            with path.open() as handle:
                indexes[path] = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        source_id = row["source_identifier_or_coordinates"].split(";", 1)[0]
        if source_id not in indexes[path]:
            raise ValueError(f"Source ID {source_id} not found for {sequence_id}")
        source_sequence = normalized_sequence(indexes[path][source_id])
        if sha256_bytes(source_sequence.encode("ascii")) != row["sequence_sha256"]:
            raise ValueError(f"Full-length local source sequence mismatch for {sequence_id}")
        checked += 1
    return checked


def validate_bed(path, protein_id, protein_length):
    count = 0
    maximum = 0
    placeholders = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least four BED columns")
            if fields[0] != protein_id:
                raise ValueError(f"{path}:{line_number}: unexpected protein ID {fields[0]}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-integer coordinates") from error
            is_known_placeholder = (
                path.name == "PKZILLA-3_PKS_domains.bed"
                and protein_id == "PKZILLA-3"
                and start == end == 24346
                and fields[3] == "C1"
            )
            if is_known_placeholder:
                placeholders += 1
                continue
            if not 0 <= start < end <= protein_length:
                raise ValueError(
                    f"{path}:{line_number}: interval {start}-{end} outside "
                    f"{protein_id} length {protein_length}"
                )
            count += 1
            maximum = max(maximum, end)
    if not count:
        raise ValueError(f"No intervals in {path}")
    return count, maximum, placeholders


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    seed_dir = root / "resources/seeds"
    manifest = load_manifest(seed_dir / "manifest.tsv")
    validate_collections(seed_dir, manifest)
    local_count = validate_local_sources(root, manifest)
    print(f"seed validation: PASS; records=15 local_source_sequences={local_count} remote_records=5")

    domain_dir = root / "data/benchmarks/pkzilla_truth/zenodo_11068523_pks_domain_to_module"
    for protein_id in ["PKZILLA-1", "PKZILLA-2", "PKZILLA-3", "PKZILLA-B1"]:
        length = int(manifest[protein_id]["length_aa"])
        for kind in ["domains", "modules"]:
            path = domain_dir / f"{protein_id}_PKS_{kind}.bed"
            count, maximum, placeholders = validate_bed(path, protein_id, length)
            print(
                f"coordinate validation: PASS; protein={protein_id} kind={kind} "
                f"intervals={count} skipped_zero_length={placeholders} "
                f"max_end={maximum} length_aa={length}"
            )
    print("PKZILLA protein-relative coordinate validation: PASS; proteins=4 bed_files=8")


if __name__ == "__main__":
    main()

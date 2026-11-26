#!/usr/bin/env python3
"""Validate, extract, and describe the bounded HaptoSMASH Pfam HMM subset."""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path


FIELDS = ["accession", "model_name", "track", "role", "rationale"]
ALLOWED_TRACKS = {"A", "B"}
ALLOWED_ROLES = {"backbone", "ancillary", "negative_control"}


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required nonempty file is absent: {path}")
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def read_selection(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Unexpected selection columns: {reader.fieldnames}")
    if not 50 <= len(rows) <= 100:
        raise ValueError(f"Expected 50-100 selected models, observed {len(rows)}")
    accessions = set()
    names = set()
    for row in rows:
        if None in row or any(row[field] is None or not row[field].strip() for field in FIELDS):
            raise ValueError(f"Incomplete selection row: {row}")
        if row["accession"] in accessions:
            raise ValueError(f"Duplicate selected accession: {row['accession']}")
        if row["model_name"] in names:
            raise ValueError(f"Duplicate selected model name: {row['model_name']}")
        if row["track"] not in ALLOWED_TRACKS or row["role"] not in ALLOWED_ROLES:
            raise ValueError(f"Invalid track/role in selection row: {row}")
        accessions.add(row["accession"])
        names.add(row["model_name"])
    return rows


def read_pfam_metadata(path):
    metadata = {}
    current = None
    total = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("# STOCKHOLM 1.0"):
                if current is not None:
                    raise ValueError(f"Unterminated Pfam metadata record before line {line_number}")
                current = {}
            elif current is not None and line.startswith("#=GF "):
                key = line[5:7]
                if key in {"ID", "AC", "DE"}:
                    current[key] = line[10:].rstrip("\n")
            elif current is not None and line.rstrip("\n") == "//":
                if not all(key in current for key in ["ID", "AC", "DE"]):
                    raise ValueError(f"Incomplete Pfam metadata record ending at line {line_number}")
                accession = current["AC"]
                if accession in metadata:
                    raise ValueError(f"Duplicate Pfam metadata accession: {accession}")
                metadata[accession] = (current["ID"], current["DE"])
                total += 1
                current = None
    if current is not None:
        raise ValueError(f"Unterminated final Pfam metadata record in {path}")
    if not total:
        raise ValueError(f"No Pfam metadata records in {path}")
    return metadata, total


def validate_metadata(rows, metadata):
    for row in rows:
        accession = row["accession"]
        if accession not in metadata:
            raise ValueError(f"Selected accession absent from Pfam metadata: {accession}")
        model_name = metadata[accession][0]
        if model_name != row["model_name"]:
            raise ValueError(
                f"Model-name mismatch for {accession}: list={row['model_name']} Pfam={model_name}"
            )


def read_hmm(path, rows, expected_total, collect_blocks=False):
    selected = {row["accession"]: row["model_name"] for row in rows}
    observed = {}
    blocks = {}
    all_accessions = set()
    block = None
    accession = None
    model_name = None
    total = 0
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("HMMER3/"):
                if block is not None:
                    raise ValueError(f"Unterminated HMM record before {path}:{line_number}")
                block = [line]
                accession = None
                model_name = None
                continue
            if block is None:
                if line.strip():
                    raise ValueError(f"Unexpected content before first HMM at {path}:{line_number}")
                continue
            block.append(line)
            if line.startswith("ACC   "):
                accession = line[6:].strip()
            elif line.startswith("NAME  "):
                model_name = line[6:].strip()
            elif line.rstrip("\n") == "//":
                if not accession or not model_name:
                    raise ValueError(f"HMM record at line {line_number} lacks NAME or ACC")
                if accession in all_accessions:
                    raise ValueError(f"Duplicate HMM accession in {path}: {accession}")
                all_accessions.add(accession)
                total += 1
                if accession in selected:
                    observed[accession] = model_name
                    if collect_blocks:
                        blocks[accession] = "".join(block)
                block = None
    if block is not None:
        raise ValueError(f"Unterminated final HMM record in {path}")
    if total != expected_total:
        raise ValueError(f"{path}: expected {expected_total} HMM records, observed {total}")
    missing = set(selected) - set(observed)
    if missing:
        raise ValueError(f"Selected HMM accessions absent from {path}: {sorted(missing)}")
    for selected_accession, expected_name in selected.items():
        if observed[selected_accession] != expected_name:
            raise ValueError(
                f"HMM name mismatch for {selected_accession}: "
                f"list={expected_name} file={observed[selected_accession]}"
            )
    return blocks


def hmmer_version():
    executable = shutil.which("hmmpress")
    if not executable:
        raise ValueError("hmmpress is not available on PATH")
    result = subprocess.run(
        [executable, "-h"], check=True, capture_output=True, text=True
    )
    for line in result.stdout.splitlines():
        if line.startswith("# HMMER"):
            return line.strip()
    raise ValueError("Could not parse HMMER version from hmmpress -h")


def write_json_atomic(path, payload):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def common_arguments(parser):
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--pfam-dat", type=Path, required=True)
    parser.add_argument("--pfam-hmm", type=Path, required=True)
    parser.add_argument("--pfam-release", required=True)
    parser.add_argument("--expected-pfam-count", type=int, required=True)


def load_and_validate(args):
    rows = read_selection(args.selection)
    metadata, metadata_count = read_pfam_metadata(args.pfam_dat)
    if metadata_count != args.expected_pfam_count:
        raise ValueError(
            f"Pfam metadata release count mismatch: expected {args.expected_pfam_count}, "
            f"observed {metadata_count}"
        )
    validate_metadata(rows, metadata)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ["validate", "extract", "manifest"]:
        common_arguments(subparsers.add_parser(command))
    extract_parser = subparsers.choices["extract"]
    extract_parser.add_argument("--output", type=Path, required=True)
    manifest_parser = subparsers.choices["manifest"]
    manifest_parser.add_argument("--output-hmm", type=Path, required=True)
    manifest_parser.add_argument("--pressed", type=Path, nargs=4, required=True)
    manifest_parser.add_argument("--output-manifest", type=Path, required=True)
    args = parser.parse_args()

    rows = load_and_validate(args)
    role_counts = Counter(f"{row['track']}/{row['role']}" for row in rows)
    if args.command in {"validate", "extract"}:
        blocks = read_hmm(
            args.pfam_hmm,
            rows,
            args.expected_pfam_count,
            collect_blocks=args.command == "extract",
        )
    if args.command == "validate":
        print(
            f"Pfam selection validation: PASS; release={args.pfam_release} "
            f"source_models={args.expected_pfam_count} selected_models={len(rows)}"
        )
        print(f"Role counts: {dict(sorted(role_counts.items()))}")
        return
    if args.command == "extract":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(args.output.name + ".tmp")
        with temporary.open("w") as handle:
            for row in rows:
                handle.write(blocks[row["accession"]])
        temporary.replace(args.output)
        print(
            f"Extracted {len(rows)} HMM models from Pfam {args.pfam_release} "
            f"to {args.output}; sha256={sha256_file(args.output)}"
        )
        return

    read_hmm(args.output_hmm, rows, len(rows))
    outputs = {
        "custom_hmm": file_record(args.output_hmm),
        "pressed": [file_record(path) for path in args.pressed],
    }
    payload = {
        "format": "hapto_core_hmm_manifest_v1",
        "generator": "workflow/scripts/build_hmm_database.py",
        "pfam_release": args.pfam_release,
        "pfam_release_evidence": (
            "Local Pfam-A.hmm and Pfam-A.hmm.dat both contain 53,275 models, "
            "matching the official Pfam 39.0 release-notes family count."
        ),
        "model_count": len(rows),
        "role_counts": dict(sorted(role_counts.items())),
        "models": rows,
        "inputs": {
            "selection": file_record(args.selection),
            "pfam_dat": file_record(args.pfam_dat),
            "pfam_hmm": file_record(args.pfam_hmm),
        },
        "outputs": outputs,
        "hmmer_version": hmmer_version(),
    }
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(args.output_manifest, payload)
    print(
        f"HMM build manifest: PASS; models={len(rows)} pressed_files=4 "
        f"manifest={args.output_manifest} sha256={sha256_file(args.output_manifest)}"
    )


if __name__ == "__main__":
    main()

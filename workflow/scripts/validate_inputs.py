#!/usr/bin/env python3
"""Validate HaptoSMASH sample and reference manifests without modifying inputs."""

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Iterable, Optional

import yaml


REQUIRED_COLUMNS = (
    "sample_id",
    "species",
    "strain",
    "assembly_accession",
    "source",
    "group",
    "genome_fasta",
    "gff3",
    "proteins_faa",
)
ALLOWED_GROUPS = {"truth_benchmark", "benchmark", "flagship", "osf_panel"}
OPTIONAL_COLUMNS = ("gff3", "proteins_faa")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="config/samples.tsv")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--root", default=".", help="Repository root for relative paths")
    parser.add_argument("--report", help="Optional per-sample TSV inventory report")
    parser.add_argument("--expected-count", type=int, default=18)
    return parser.parse_args()


def input_path(root: Path, value: str, context: str, errors: list[str]) -> Optional[Path]:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        errors.append(f"{context}: path must be relative and remain lexically below the repository: {value}")
        return None
    return root / path


def fasta_starts_correctly(path: Path) -> bool:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    return line.startswith(">")
    except (OSError, UnicodeError, gzip.BadGzipFile):
        return False
    return False


def check_file(
    root: Path,
    value: str,
    context: str,
    errors: list[str],
    *,
    fasta: bool = False,
) -> tuple[Optional[Path], int]:
    path = input_path(root, value, context, errors)
    if path is None:
        return None, 0
    if not path.is_file():
        errors.append(f"{context}: file does not exist: {value}")
        return path, 0
    size = path.stat().st_size
    if size == 0:
        errors.append(f"{context}: file is empty: {value}")
    elif fasta and not fasta_starts_correctly(path):
        errors.append(f"{context}: not a readable FASTA file: {value}")
    return path, size


def iter_config_paths(config: dict) -> Iterable[tuple[str, str]]:
    for section in ("databases", "benchmarks", "seeds"):
        for key, value in config.get(section, {}).items():
            if isinstance(value, str):
                yield f"config.{section}.{key}", value


def validate(args: argparse.Namespace) -> tuple[list[str], list[dict[str, str]]]:
    root = Path(args.root).resolve()
    manifest_path = input_path(root, args.manifest, "manifest", errors := [])
    config_path = input_path(root, args.config, "config", errors)
    rows: list[dict[str, str]] = []
    report_rows: list[dict[str, str]] = []

    if manifest_path is None or not manifest_path.is_file():
        errors.append(f"Manifest does not exist: {args.manifest}")
        return errors, report_rows
    if config_path is None or not config_path.is_file():
        errors.append(f"Config does not exist: {args.config}")
        return errors, report_rows

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            errors.append(f"Manifest is missing required columns: {', '.join(missing)}")
            return errors, report_rows
        rows = list(reader)

    if args.expected_count >= 0 and len(rows) != args.expected_count:
        errors.append(f"Expected {args.expected_count} samples, found {len(rows)}")

    seen_ids: set[str] = set()
    seen_accessions: set[str] = set()
    for line_number, row in enumerate(rows, start=2):
        sample_id = (row.get("sample_id") or "").strip()
        context = f"line {line_number} ({sample_id or 'missing sample_id'})"
        if not sample_id:
            errors.append(f"line {line_number}: empty sample_id")
        elif sample_id in seen_ids:
            errors.append(f"{context}: duplicate sample_id")
        seen_ids.add(sample_id)

        for column in REQUIRED_COLUMNS:
            if not (row.get(column) or "").strip():
                errors.append(f"{context}: empty required field {column}")
        group = (row.get("group") or "").strip()
        if group and group not in ALLOWED_GROUPS:
            errors.append(f"{context}: unsupported group {group!r}")
        species = (row.get("species") or "").strip()
        if species.endswith("_sp"):
            errors.append(f"{context}: unresolved species identifier {species!r}")
        strain = (row.get("strain") or "").strip()
        if strain.upper() in {"NA", "UNKNOWN", "UNRESOLVED"}:
            errors.append(f"{context}: unresolved strain identifier")
        accession = (row.get("assembly_accession") or "").strip()
        if accession in seen_accessions:
            errors.append(f"{context}: duplicate assembly_accession {accession!r}")
        seen_accessions.add(accession)

        report = {"sample_id": sample_id, "species": species, "strain": strain, "group": group}
        for column in ("genome_fasta", "gff3", "proteins_faa"):
            value = (row.get(column) or "").strip()
            if column in OPTIONAL_COLUMNS and value == "NA":
                report[f"{column}_bytes"] = "0"
                continue
            if not value:
                report[f"{column}_bytes"] = "0"
                continue
            _, size = check_file(
                root, value, f"{context}.{column}", errors, fasta=column != "gff3"
            )
            report[f"{column}_bytes"] = str(size)
        report_rows.append(report)

    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    if config.get("samples_tsv") != args.manifest:
        errors.append(
            f"config.samples_tsv is {config.get('samples_tsv')!r}, expected {args.manifest!r}"
        )
    for context, value in iter_config_paths(config):
        path = input_path(root, value, context, errors)
        if path is None:
            continue
        if not path.exists():
            errors.append(f"{context}: path does not exist: {value}")
        elif path.is_file() and path.stat().st_size == 0:
            errors.append(f"{context}: file is empty: {value}")
        elif path.is_dir() and not any(path.iterdir()):
            errors.append(f"{context}: directory is empty: {value}")

    return errors, report_rows


def write_report(path_value: str, rows: list[dict[str, str]]) -> None:
    path = Path(path_value)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "species", "strain", "group", "genome_fasta_bytes", "gff3_bytes", "proteins_faa_bytes"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    errors, rows = validate(args)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    if args.report:
        write_report(args.report, rows)
    total_bytes = sum(int(row["genome_fasta_bytes"]) for row in rows)
    print(f"Validated {len(rows)} samples ({total_bytes:,} compressed/uncompressed genome bytes).")
    print("All configured database, benchmark, and seed paths are present and non-empty.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate and validate the tiny HaptoSMASH custom-HMM scan fixture."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    return sha256_bytes(path.read_bytes())


def load_records(path):
    with path.open() as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def make_fixture(args):
    pkzilla = load_records(args.seeds_dir / "pkzilla_seeds.faa")
    oxylipin = load_records(args.seeds_dir / "halogenase_oxylipin_seeds.faa")
    ks_sequence = str(pkzilla["PKZILLA-1"].seq)[1350:1786]
    if len(ks_sequence) != 436:
        raise ValueError("Unexpected PKZILLA-1 KS1 fixture length")
    records = [
        SeqRecord(
            Seq(ks_sequence),
            id="positive_KS",
            description="source=PKZILLA-1 protein_start=1350 protein_end=1786 coordinates=0-based-half-open",
        ),
        SeqRecord(
            Seq(str(oxylipin["LOX_EMIHU"].seq)),
            id="positive_LOX",
            description="source=LOX_EMIHU full_length=true",
        ),
        SeqRecord(
            Seq("ACDEFGHIKLMNPQRSTVWY" * 6),
            id="negative",
            description="synthetic=true expected_calls=0",
        ),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        SeqIO.write(records, handle, "fasta")
    expected = {
        "format": "hapto_core_hmm_smoke_v1",
        "fixture_sha256": sha256_file(args.output),
        "expectations": {
            "positive_KS": {"required_accessions": ["PF00109.33"]},
            "positive_LOX": {"required_accessions": ["PF00305.26"]},
            "negative": {"required_accessions": [], "maximum_calls": 0},
        },
        "records": {
            record.id: {
                "length_aa": len(record.seq),
                "sequence_sha256": sha256_bytes(str(record.seq).encode("ascii")),
                "description": record.description,
            }
            for record in records
        },
    }
    write_json(args.expected, expected)
    print(
        f"Smoke fixture: records=3 fasta={args.output} "
        f"sha256={expected['fixture_sha256']} expected={args.expected}"
    )


def read_domtbl(path):
    calls = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split(None, 22)
            if len(fields) < 22:
                raise ValueError(f"Malformed domtblout row at {path}:{line_number}")
            calls.append({
                "query_id": fields[3],
                "model_name": fields[0],
                "accession": fields[1],
                "independent_evalue": fields[12],
                "bitscore": fields[13],
                "ali_from": fields[17],
                "ali_to": fields[18],
            })
    return sorted(
        calls,
        key=lambda row: (
            row["query_id"], row["accession"], int(row["ali_from"]), int(row["ali_to"])
        ),
    )


def validate(args):
    expected = json.loads(args.expected.read_text())
    if sha256_file(args.fixture) != expected["fixture_sha256"]:
        raise ValueError("Smoke FASTA checksum does not match expected JSON")
    fixture_records = load_records(args.fixture)
    if set(fixture_records) != set(expected["expectations"]):
        raise ValueError("Smoke FASTA IDs do not match expected JSON")
    calls = read_domtbl(args.domtbl)
    by_query = {}
    for call in calls:
        by_query.setdefault(call["query_id"], set()).add(call["accession"])
    for query_id, expectation in expected["expectations"].items():
        observed = by_query.get(query_id, set())
        missing = set(expectation["required_accessions"]) - observed
        if missing:
            raise ValueError(f"Missing required calls for {query_id}: {sorted(missing)}")
        maximum = expectation.get("maximum_calls")
        if maximum is not None and len(by_query.get(query_id, set())) > maximum:
            raise ValueError(f"Unexpected calls for {query_id}: {sorted(observed)}")
        if maximum == 0 and any(call["query_id"] == query_id for call in calls):
            raise ValueError(f"Negative fixture {query_id} has domain calls")
    args.calls.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "query_id", "model_name", "accession", "independent_evalue",
        "bitscore", "ali_from", "ali_to",
    ]
    with args.calls.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(calls)
    print(
        f"HMM scan smoke: PASS; records=3 domain_calls={len(calls)} "
        f"positive_KS=PF00109.33 positive_LOX=PF00305.26 negative_calls=0 "
        f"calls={args.calls}"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    fixture = subparsers.add_parser("make-fixture")
    fixture.add_argument("--seeds-dir", type=Path, required=True)
    fixture.add_argument("--output", type=Path, required=True)
    fixture.add_argument("--expected", type=Path, required=True)
    scan = subparsers.add_parser("validate")
    scan.add_argument("--fixture", type=Path, required=True)
    scan.add_argument("--domtbl", type=Path, required=True)
    scan.add_argument("--expected", type=Path, required=True)
    scan.add_argument("--calls", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "make-fixture":
        make_fixture(args)
    else:
        validate(args)


if __name__ == "__main__":
    main()

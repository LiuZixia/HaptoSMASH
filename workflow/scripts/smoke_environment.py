#!/usr/bin/env python3
"""Check imports and command-line entry points for one HaptoSMASH environment."""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path


def check_imports(modules: tuple[str, ...]) -> None:
    for module in modules:
        imported = importlib.import_module(module)
        version = getattr(imported, "__version__", None)
        print(f"import {module}: OK" + (f" ({version})" if version else ""))


def check_executables(names: tuple[str, ...]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError(f"missing executables: {', '.join(missing)}")
    for name in names:
        print(f"executable {name}: {shutil.which(name)}")


def run_command(command: list[str]) -> None:
    print("+ " + " ".join(command))
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode:
        raise RuntimeError(f"command exited {result.returncode}: {' '.join(command)}")


def check_whokaryote_models() -> None:
    import joblib
    import whokaryote_scripts

    data_dir = Path(whokaryote_scripts.__file__).resolve().parent / "data"
    for name in ("RF_S_RBS_240121_bal.joblib", "RF_T_RBS_310122_bal.joblib"):
        path = data_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing Whokaryote model: {path}")
        joblib.load(path)
        print(f"Whokaryote model {name}: load OK ({path.stat().st_size} bytes)")


def checks_for(mode: str):
    if mode == "driver":
        yield "driver imports", lambda: check_imports(("snakemake", "conda_lock", "yaml"))
        yield "driver executables", lambda: check_executables(
            ("snakemake", "conda-lock", "datasets", "dataformat", "aria2c", "gh")
        )
        yield "Snakemake version", lambda: run_command(["snakemake", "--version"])
        yield "conda-lock version", lambda: run_command(["conda-lock", "--version"])
    elif mode == "core":
        yield "core imports", lambda: check_imports(
            ("pyhmmer", "Bio", "polars", "pandas", "numpy", "scipy", "pyranges", "yaml", "pytest")
        )
        yield "core executables", lambda: check_executables(
            (
                "hmmpress", "hmmscan", "seqkit", "bedtools", "diamond", "mmseqs",
                "miniprot", "mafft", "trimal", "iqtree", "hisat2", "samtools", "stringtie",
            )
        )
        yield "HMMER version", lambda: run_command(["hmmpress", "-h"])
        yield "seqkit version", lambda: run_command(["seqkit", "version"])
        yield "bedtools version", lambda: run_command(["bedtools", "--version"])
    elif mode == "decontam":
        yield "decontamination imports", lambda: check_imports(
            ("joblib", "pandas", "numpy", "sklearn", "whokaryote_scripts")
        )
        yield "decontamination executables", lambda: check_executables(
            ("tiara", "whokaryote.py", "prodigal")
        )
        yield "Tiara help", lambda: run_command(["tiara", "--help"])
        yield "Whokaryote help", lambda: run_command(["whokaryote.py", "--help"])
        yield "Prodigal version", lambda: run_command(["prodigal", "-v"])
        yield "Whokaryote models", check_whokaryote_models
    elif mode == "engine":
        yield "engine imports", lambda: check_imports(
            (
                "antismash", "Bio", "brawn", "helperlibs", "jsonschema", "BCBio",
                "nrpys", "MOODS", "orjson", "sklearn",
            )
        )
        yield "engine executables", lambda: check_executables(
            ("hmmscan", "hmmpress", "hmmsearch2", "diamond", "fasttree", "prodigal", "blastp")
        )
        yield "antiSMASH version", lambda: run_command([sys.executable, "-m", "antismash", "--version"])
    else:
        raise ValueError(mode)


def main() -> int:
    env_bin = str(Path(sys.executable).resolve().parent)
    os.environ["PATH"] = env_bin + os.pathsep + os.environ.get("PATH", "")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("driver", "core", "decontam", "engine"))
    args = parser.parse_args()
    failures = []
    for label, check in checks_for(args.mode):
        print(f"\n== {label} ==")
        try:
            check()
        except Exception as error:
            failures.append((label, str(error)))
            print(f"FAIL: {error}", file=sys.stderr)
    if failures:
        print(f"\n{args.mode} smoke test failed: {len(failures)} check group(s)", file=sys.stderr)
        return 1
    print(f"\n{args.mode} smoke test: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

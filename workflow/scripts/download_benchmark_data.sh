#!/usr/bin/env bash
set -euo pipefail

# Retry each failed download 10 times by default.
# Use --task-timeout SECONDS to retry until the shared download deadline.
MAX_RETRIES=10
TASK_TIMEOUT=0
RETRY_INTERVAL=5
CHECK_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-retries)
            [[ $# -ge 2 && "$2" =~ ^[0-9]+$ ]] || { echo "Error: --max-retries requires a non-negative integer." >&2; exit 2; }
            MAX_RETRIES="$2"
            shift 2
            ;;
        --task-timeout)
            [[ $# -ge 2 && "$2" =~ ^[1-9][0-9]*$ ]] || { echo "Error: --task-timeout requires a positive integer in seconds." >&2; exit 2; }
            TASK_TIMEOUT="$2"
            shift 2
            ;;
        --check-only)
            CHECK_ONLY=1
            shift
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            exit 2
            ;;
    esac
done

TASK_DEADLINE=$(( $(date +%s) + TASK_TIMEOUT ))

run_with_retry() {
    local attempt=1
    local status now remaining sleep_time

    while true; do
        if (( TASK_TIMEOUT > 0 )); then
            now=$(date +%s)
            remaining=$((TASK_DEADLINE - now))
            if (( remaining <= 0 )); then
                echo "Download deadline exceeded: $*" >&2
                return 124
            fi
            echo "Attempt ${attempt}; ${remaining} seconds remain: $*"
            if timeout "${remaining}s" "$@"; then
                return 0
            fi
            status=$?
        else
            echo "Attempt ${attempt}: $*"
            if "$@"; then
                return 0
            fi
            status=$?
            if (( attempt > MAX_RETRIES )); then
                echo "Download failed after ${MAX_RETRIES} retries: $*" >&2
                return "${status}"
            fi
        fi

        echo "Retrying in ${RETRY_INTERVAL} seconds: $*" >&2
        if (( TASK_TIMEOUT > 0 )); then
            now=$(date +%s)
            remaining=$((TASK_DEADLINE - now))
            if (( remaining <= 0 )); then
                echo "Download deadline exceeded." >&2
                return 124
            fi
            sleep_time=$(( remaining < RETRY_INTERVAL ? remaining : RETRY_INTERVAL ))
            sleep "${sleep_time}"
        else
            sleep "${RETRY_INTERVAL}"
        fi
        attempt=$((attempt + 1))
    done
}

# ==============================================================================
# HaptoSMASH automated download script for benchmark data, standalone databases
# and the supplementary benchmark sets (NGDC GWH, Zenodo PKZILLA, OSF vka93).
# Runtime environment: conda activate haptosmash-driver
# ==============================================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DB_DIR="${REPO_ROOT}/db"
DATA_DIR="${REPO_ROOT}/data"
SEED_DIR="${REPO_ROOT}/resources/seeds"
TMP_DIR="${DATA_DIR}/_tmp_downloads"

# Supplementary benchmark data (formerly download_extra_benchmarks.sh)
GWH_DIR="${DATA_DIR}/genomes/phaeocystis_globosa_CNS00066"
PKZ_DIR="${DATA_DIR}/benchmarks/pkzilla_truth"
OSF_DIR="${DATA_DIR}/genomes/phaeocystis_osf_vka93"

DECOMPRESS_CMD="gzip -dc"
command -v pigz >/dev/null 2>&1 && DECOMPRESS_CMD="pigz -dc"

# OSF rate-limit tuning (see step [8/8]). The files.osf.io backend throttles by
# request RATE, not just per-connection bandwidth, and 403s aria2's default
# User-Agent. The bulk is fetched fast with aria2c (browser UA); the long tail of
# tiny files is then cleared by a rate-paced curl resume that backs off on 429.
OSF_WORKERS="${OSF_WORKERS:-2}"        # parallel curl workers for the resume pass
OSF_SLEEP="${OSF_SLEEP:-2}"            # seconds each worker pauses between files
OSF_COOLDOWN="${OSF_COOLDOWN:-45}"     # seconds to rest between resume passes
OSF_MAX_PASSES="${OSF_MAX_PASSES:-6}"  # capped passes so a run always terminates

if (( CHECK_ONLY )); then
    python3 - "${DB_DIR}" "${DATA_DIR}" "${SEED_DIR}" <<'PYEOF'
import sys
from pathlib import Path


db_dir, data_dir, seed_dir = map(Path, sys.argv[1:4])
required_files = {
    "NCBI haptophyte catalog": data_dir / "genomes/catalog/ncbi_haptophyta_all.tsv",
    "NCBI P. parvum catalog": data_dir / "genomes/catalog/ncbi_prymnesium_parvum.tsv",
    "Pfam-A HMM": db_dir / "pfam/Pfam-A.hmm",
    "Pfam-A metadata": db_dir / "pfam/Pfam-A.hmm.dat",
    "MIBiG 4.0 proteins": db_dir / "mibig/mibig_prot_seqs_4.0.fasta",
    "PKZILLA-1 seed": seed_dir / "pkzilla1_uniprot_A0AB34IYJ6.faa",
}
for sample in (
    "emiliania_huxleyi_CCMP1516",
    "diacronema_lutheri_CCMP1325",
    "chrysochromulina_tobinii_CCMP291",
):
    required_files[f"{sample} genome"] = data_dir / f"genomes/{sample}/genome.fna"
for name in ("genome.fna", "annotation.gff3", "proteins.faa", "cds.fna", "rna.fna"):
    required_files[f"CNS00066 {name}"] = data_dir / f"genomes/phaeocystis_globosa_CNS00066/{name}"

missing = 0
for label, path in required_files.items():
    if path.is_file() and path.stat().st_size > 0:
        print(f"    [skip] {label}: {path}")
    else:
        print(f"    [pending] {label}: {path}")
        missing += 1
pp_dir = data_dir / "genomes/prymnesium_parvum_ncbi"
pp_complete = (pp_dir / ".done").is_file() or any(
    path.stat().st_size > 0 for path in pp_dir.glob("GC*/genome.fna")
)
print(f"    [{'skip' if pp_complete else 'pending'}] P. parvum collection: {pp_dir}")
missing += 0 if pp_complete else 1
for label, path in {
    "Zenodo domain/module content": data_dir / "benchmarks/pkzilla_truth/zenodo_11068523_pks_domain_to_module/12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa",
    "Zenodo hotspot content": data_dir / "benchmarks/pkzilla_truth/zenodo_10309063_hotspot_liftover/2023-12-08_Hotspot_liftover/v1.1_liftover_PKS_hotspots_5000.merge.bed",
}.items():
    complete = path.is_file() and path.stat().st_size > 0
    print(f"    [{'skip' if complete else 'pending'}] {label}: {path}")
    missing += 0 if complete else 1
for label, path in {
    "MIBiG JSON": db_dir / "mibig/json/mibig_json_4.0",
    "MIBiG GenBank": db_dir / "mibig/gbk/mibig_gbk_4.0",
}.items():
    if path.is_dir() and any(path.iterdir()):
        print(f"    [skip] {label}: {path}")
    else:
        print(f"    [pending] {label}: {path}")
        missing += 1

manifest = data_dir / "genomes/phaeocystis_osf_vka93/osf_vka93_manifest.tsv"
complete = pending = 0
if not manifest.is_file() or manifest.stat().st_size == 0:
    print(f"    [pending] OSF manifest: {manifest}")
    missing += 1
else:
    with manifest.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                print(f"    [pending] OSF manifest line {line_number}: expected three columns")
                pending += 1
                continue
            relative, expected_text, _url = fields
            target = data_dir / "genomes/phaeocystis_osf_vka93" / relative
            try:
                expected = int(expected_text)
            except ValueError:
                expected = -1
            if expected >= 0 and target.is_file() and target.stat().st_size == expected:
                complete += 1
            else:
                pending += 1
    print(f"    OSF resume accounting: complete={complete} pending={pending} total={complete + pending}")

print(f"Offline inventory: required_pending={missing} osf_pending={pending}")
raise SystemExit(1 if missing or pending else 0)
PYEOF
    exit $?
fi

mkdir -p "${DB_DIR}/pfam" "${DB_DIR}/mibig" \
         "${DATA_DIR}/genomes/catalog" \
         "${DATA_DIR}/benchmarks/pkzilla_truth" \
         "${GWH_DIR}" "${OSF_DIR}" \
         "${SEED_DIR}" "${TMP_DIR}"

echo ">>> [1/8] Exporting NCBI assembly summaries for the full Haptophyta phylum (TaxID: 2830) and Prymnesium parvum (TaxID: 97485)..."
datasets summary genome taxon 2830 --as-json-lines | \
  dataformat tsv genome \
    --fields accession,organism-name,organism-infraspecific-strain,assminfo-level,assmstats-total-sequence-len,assminfo-release-date \
  > "${DATA_DIR}/genomes/catalog/ncbi_haptophyta_all.tsv"

datasets summary genome taxon 97485 --as-json-lines | \
  dataformat tsv genome \
    --fields accession,organism-name,organism-infraspecific-strain,assminfo-level,assmstats-total-sequence-len,assminfo-release-date \
  > "${DATA_DIR}/genomes/catalog/ncbi_prymnesium_parvum.tsv"

echo "    Generated catalog table: data/genomes/catalog/ncbi_haptophyta_all.tsv"

# Generic unpack/normalise function: extract a datasets zip into standard genome.fna / annotation.gff3 / proteins.faa
unpack_ncbi_dataset() {
    local zip_file="$1"
    local target_dir="$2"
    local acc="$3"

    mkdir -p "${target_dir}"
    local work_dir="${TMP_DIR}/unpack_${acc}"
    rm -rf "${work_dir}"
    unzip -q "${zip_file}" -d "${work_dir}"

    local src_dir="${work_dir}/ncbi_dataset/data/${acc}"
    if [[ -d "${src_dir}" ]]; then
        find "${src_dir}" -maxdepth 1 -name "*_genomic.fna" -exec cp -f {} "${target_dir}/genome.fna" \;
        [[ -f "${src_dir}/genomic.gff" ]] && cp -f "${src_dir}/genomic.gff" "${target_dir}/annotation.gff3"
        [[ -f "${src_dir}/protein.faa" ]] && cp -f "${src_dir}/protein.faa" "${target_dir}/proteins.faa"
        cp -f "${work_dir}/ncbi_dataset/data/assembly_data_report.jsonl" "${target_dir}/" 2>/dev/null || true
    fi
    rm -rf "${work_dir}"
}

echo ">>> [2/8] Downloading NCBI core haptophyte reference genomes (E. huxleyi / D. lutheri / C. tobinii)..."
declare -A CORE_GENOMES=(
    ["emiliania_huxleyi_CCMP1516"]="GCF_000372725.1"
    ["diacronema_lutheri_CCMP1325"]="GCA_019693415.1"
    ["chrysochromulina_tobinii_CCMP291"]="GCA_001275005.1"
)

for sample in "${!CORE_GENOMES[@]}"; do
    acc="${CORE_GENOMES[$sample]}"
    out_dir="${DATA_DIR}/genomes/${sample}"
    if [[ -s "${out_dir}/genome.fna" ]]; then
        echo "    [skip] ${sample} (${acc}) already exists"
        continue
    fi
    echo "    [download] ${sample} (${acc})..."
    zip_path="${TMP_DIR}/${acc}.zip"
    run_with_retry datasets download genome accession "${acc}" \
        --include genome,gff3,protein \
        --filename "${zip_path}"
    unpack_ncbi_dataset "${zip_path}" "${out_dir}" "${acc}"
    rm -f "${zip_path}"
done

echo ">>> [3/8] Downloading publicly available NCBI Prymnesium parvum genomes (including the 12B1 v1.1 / Type B strain)..."
PP_DIR="${DATA_DIR}/genomes/prymnesium_parvum_ncbi"
PP_HAS_GENOMES=0
for genome_path in "${PP_DIR}"/GC*/genome.fna; do
    if [[ -s "${genome_path}" ]]; then
        PP_HAS_GENOMES=1
        break
    fi
done
if [[ ! -f "${PP_DIR}/.done" && "${PP_HAS_GENOMES}" -eq 0 ]]; then
    zip_path="${TMP_DIR}/prymnesium_parvum_all.zip"
    run_with_retry datasets download genome taxon 97485 \
        --include genome,gff3,protein \
        --dehydrated \
        --filename "${zip_path}"
    rehydrate_dir="${TMP_DIR}/prymnesium_parvum_rehydrate"
    if [[ ! -d "${rehydrate_dir}/ncbi_dataset" ]]; then
        unzip -q "${zip_path}" -d "${rehydrate_dir}"
    fi
    run_with_retry datasets rehydrate --directory "${rehydrate_dir}"
    mkdir -p "${PP_DIR}"
    for acc_dir in "${rehydrate_dir}/ncbi_dataset/data/"GC*; do
        [[ -d "${acc_dir}" ]] || continue
        acc="$(basename "${acc_dir}")"
        sub_dir="${PP_DIR}/${acc}"
        mkdir -p "${sub_dir}"
        find "${acc_dir}" -maxdepth 1 -name "*_genomic.fna" -exec cp -f {} "${sub_dir}/genome.fna" \;
        [[ -f "${acc_dir}/genomic.gff" ]] && cp -f "${acc_dir}/genomic.gff" "${sub_dir}/annotation.gff3"
        [[ -f "${acc_dir}/protein.faa" ]] && cp -f "${acc_dir}/protein.faa" "${sub_dir}/proteins.faa"
        echo "    Normalised P. parvum assembly: ${sub_dir}"
    done
    cp -f "${rehydrate_dir}/ncbi_dataset/data/assembly_data_report.jsonl" "${PP_DIR}/" 2>/dev/null || true
    touch "${PP_DIR}/.done"
    rm -rf "${rehydrate_dir}" "${zip_path}"
else
    echo "    [skip] Prymnesium parvum NCBI collection already exists"
fi

echo ">>> [4/8] Downloading standalone latest Pfam-A and MIBiG 4.0 databases..."
# 4.1 Standalone latest Pfam-A (used by the external PyHMMER extra-large hotspot engine and custom-rule database extraction; antiSMASH bundles Pfam 35.0)
if [[ ! -s "${DB_DIR}/pfam/Pfam-A.hmm" ]]; then
    run_with_retry aria2c -c --max-tries=1 -x 8 -s 8 -d "${DB_DIR}/pfam" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    run_with_retry aria2c -c --max-tries=1 -x 8 -s 8 -d "${DB_DIR}/pfam" \
        https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    pigz -d "${DB_DIR}/pfam/Pfam-A.hmm.gz"
    pigz -d "${DB_DIR}/pfam/Pfam-A.hmm.dat.gz"
else
    echo "    [skip] Pfam-A.hmm already exists"
fi

# 4.2 Standalone MIBiG 4.0 (JSON + GBK + protein FASTA)
MIBIG_BASE="https://dl.secondarymetabolites.org/mibig"
if [[ ! -s "${DB_DIR}/mibig/mibig_prot_seqs_4.0.fasta" ]]; then
    run_with_retry aria2c -c --max-tries=1 -x 8 -s 8 -d "${DB_DIR}/mibig" "${MIBIG_BASE}/mibig_json_4.0.tar.gz"
    run_with_retry aria2c -c --max-tries=1 -x 8 -s 8 -d "${DB_DIR}/mibig" "${MIBIG_BASE}/mibig_gbk_4.0.tar.gz"
    run_with_retry aria2c -c --max-tries=1 -x 8 -s 8 -d "${DB_DIR}/mibig" "${MIBIG_BASE}/mibig_prot_seqs_4.0.fasta"
    mkdir -p "${DB_DIR}/mibig/json" "${DB_DIR}/mibig/gbk"
    tar -xzf "${DB_DIR}/mibig/mibig_json_4.0.tar.gz" -C "${DB_DIR}/mibig/json"
    tar -xzf "${DB_DIR}/mibig/mibig_gbk_4.0.tar.gz" -C "${DB_DIR}/mibig/gbk"
    rm -f "${DB_DIR}/mibig/mibig_json_4.0.tar.gz" "${DB_DIR}/mibig/mibig_gbk_4.0.tar.gz"
else
    echo "    [skip] MIBiG 4.0 already exists"
fi

echo ">>> [5/8] Automatically fetching the UniProt PKZILLA-1 positive-control seed protein (A0AB34IYJ6, 45,212 aa)..."
if [[ ! -s "${SEED_DIR}/pkzilla1_uniprot_A0AB34IYJ6.faa" ]]; then
    run_with_retry curl -fSL -C - "https://rest.uniprot.org/uniprotkb/A0AB34IYJ6.fasta" \
        -o "${SEED_DIR}/pkzilla1_uniprot_A0AB34IYJ6.faa"
fi

# ==============================================================================
# Supplementary benchmark sets (merged from the former download_extra_benchmarks.sh)
# ==============================================================================

echo ">>> [6/8] Downloading NGDC GWH Phaeocystis globosa CNS00066 (GWHBJCL00000000)..."
GWH_BASE="https://download.cncb.ac.cn/gwh/Others/Phaeocystis_globosa_Pgl_GWHBJCL00000000"

declare -A GWH_FILES=(
    ["GWHBJCL00000000.genome.fasta.gz"]="genome.fna"
    ["GWHBJCL00000000.gff.gz"]="annotation.gff3"
    ["GWHBJCL00000000.Protein.faa.gz"]="proteins.faa"
    ["GWHBJCL00000000.CDS.fasta.gz"]="cds.fna"
    ["GWHBJCL00000000.RNA.fasta.gz"]="rna.fna"
)

for remote_name in "${!GWH_FILES[@]}"; do
    local_name="${GWH_FILES[$remote_name]}"
    out_file="${GWH_DIR}/${local_name}"
    if [[ -s "${out_file}" ]]; then
        echo "    [SKIP] ${local_name} already exists"
        continue
    fi
    echo "    [FETCH] ${remote_name} -> ${local_name}"
    run_with_retry aria2c -c -x 4 -s 4 --retry-wait=5 --max-tries=5 \
        -d "${TMP_DIR}" -o "${remote_name}" \
        "${GWH_BASE}/${remote_name}"
    ${DECOMPRESS_CMD} "${TMP_DIR}/${remote_name}" > "${out_file}.tmp"
    mv "${out_file}.tmp" "${out_file}"
    rm -f "${TMP_DIR}/${remote_name}"
done

echo ">>> [7/8] Downloading Zenodo PKZILLA ground-truth datasets (11068523 & 10309063)..."

download_and_unzip() {
    local url="$1"
    local zip_name="$2"
    local dest_dir="$3"
    local sentinel="$4"

    if [[ -f "${dest_dir}/.done" || -s "${dest_dir}/${sentinel}" ]]; then
        echo "    [SKIP] $(basename "${dest_dir}") already extracted"
        return
    fi
    echo "    [FETCH] ${zip_name}"
    run_with_retry aria2c -c -x 4 -s 4 --retry-wait=5 --max-tries=5 \
        -d "${TMP_DIR}" -o "${zip_name}" "${url}"
    mkdir -p "${dest_dir}"
    unzip -qo "${TMP_DIR}/${zip_name}" -d "${dest_dir}"
    touch "${dest_dir}/.done"
    rm -f "${TMP_DIR}/${zip_name}"
}

download_and_unzip \
    "https://zenodo.org/records/11068523/files/PKS%20domain%20to%20module%20conversion%20script.zip?download=1" \
    "zenodo_11068523_pks_domain_to_module.zip" \
    "${PKZ_DIR}/zenodo_11068523_pks_domain_to_module" \
    "12B1_v1.1.gff3__12B1_scaffolds_v1.1.fasta.pep.fa"

download_and_unzip \
    "https://zenodo.org/records/10309063/files/Hotspot_liftover.zip?download=1" \
    "zenodo_10309063_hotspot_liftover.zip" \
    "${PKZ_DIR}/zenodo_10309063_hotspot_liftover" \
    "2023-12-08_Hotspot_liftover/v1.1_liftover_PKS_hotspots_5000.merge.bed"

echo ">>> [8/8] Indexing and downloading OSF vka93 (Phaeocystis genomes & annotations)..."
MANIFEST="${OSF_DIR}/osf_vka93_manifest.tsv"

# (a) Index every file in the vka93 project tree via the paginated OSF API v2.
python3 - "${MANIFEST}" << 'PYEOF'
import json
import os
import re
import sys
import urllib.request

manifest_path = sys.argv[1]
HEADERS = {"User-Agent": "HaptoSMASH-Downloader/1.0"}


def fetch_json(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def iter_paginated(url):
    while url:
        payload = fetch_json(url)
        for item in payload.get("data", []):
            yield item
        url = payload.get("links", {}).get("next")


def sanitize(name):
    return re.sub(r"[^\w.\-]+", "_", name.strip()).strip("_") or "unnamed"


records = []


def walk_files(url, prefix):
    for item in iter_paginated(url):
        attrs = item.get("attributes", {})
        kind = attrs.get("kind")
        name = attrs.get("name", "unnamed")
        rel_path = os.path.join(prefix, name) if prefix else name
        if kind == "folder":
            sub_url = item["relationships"]["files"]["links"]["related"]["href"]
            walk_files(sub_url, rel_path)
        elif kind == "file":
            dl_url = item.get("links", {}).get("download")
            size = attrs.get("size") or 0
            if dl_url:
                records.append((rel_path, str(size), dl_url))


def walk_node(node_id, prefix=""):
    # Traverse all storage providers attached to this node (typically osfstorage)
    providers_url = f"https://api.osf.io/v2/nodes/{node_id}/files/"
    for prov in iter_paginated(providers_url):
        files_url = prov["relationships"]["files"]["links"]["related"]["href"]
        walk_files(files_url, prefix)

    # Recursively traverse child components
    children_url = f"https://api.osf.io/v2/nodes/{node_id}/children/"
    for child in iter_paginated(children_url):
        child_id = child["id"]
        title = sanitize(child.get("attributes", {}).get("title", child_id))
        sub_prefix = os.path.join(prefix, title) if prefix else title
        walk_node(child_id, sub_prefix)


walk_node("vka93")

with open(manifest_path, "w", encoding="utf-8") as fh:
    for rel_path, size, dl_url in records:
        fh.write(f"{rel_path}\t{size}\t{dl_url}\n")

total_mb = sum(int(r[1]) for r in records) / (1024 * 1024)
print(f"    OSF vka93 manifest indexed: {len(records)} files, estimated total {total_mb:.1f} MB")
PYEOF

# Helper: how many manifest entries are not yet present at their expected size.
osf_count_pending() {
    local pend=0 rel sz _u target
    while IFS=$'\t' read -r rel sz _u; do
        [[ -n "${rel}" ]] || continue
        target="${OSF_DIR}/${rel}"
        if [[ -f "${target}" && "$(stat -c%s "${target}" 2>/dev/null)" == "${sz}" ]]; then
            :
        else
            pend=$((pend + 1))
        fi
    done < "${MANIFEST}"
    echo "${pend}"
}

# Helper: rebuild a pending list (url<TAB>size<TAB>rel) of files still missing/short.
osf_build_pending() {
    local rel sz url target
    : > "${OSF_PENDING}"
    while IFS=$'\t' read -r rel sz url; do
        [[ -n "${rel}" ]] || continue
        target="${OSF_DIR}/${rel}"
        if [[ -f "${target}" && "${sz}" -gt 0 && "$(stat -c%s "${target}" 2>/dev/null)" == "${sz}" ]]; then
            continue
        fi
        printf '%s\t%s\t%s\n' "${url}" "${sz}" "${rel}" >> "${OSF_PENDING}"
    done < "${MANIFEST}"
    wc -l < "${OSF_PENDING}"
}

# Helper: fetch ONE OSF file with browser UA and rate-limit-aware curl back-off.
# Skips files already at their expected size so it is safe to re-run any number of times.
osf_fetch_one() {
    local line="$1" url size rel target
    IFS=$'\t' read -r url size rel <<< "${line}"
    [[ -n "${url:-}" && -n "${rel:-}" ]] || return 0
    target="${OSF_DIR}/${rel}"
    if [[ -n "${size:-}" && "${size}" != "0" && -f "${target}" && "$(stat -c%s "${target}" 2>/dev/null)" == "${size}" ]]; then
        return 0
    fi
    mkdir -p "$(dirname "${target}")"
    curl -fsSL -A "Mozilla/5.0" --connect-timeout 20 --max-time 180 \
        --retry 10 --retry-delay 20 --retry-max-time 600 --retry-all-errors \
        -C - -o "${target}" "${url}" \
        || echo "      [FAIL] ${rel}"
    sleep "${OSF_SLEEP}"
}

# (b) Fast bulk pass: aria2c -i fetches every pending file concurrently (-j 12 files,
#     -x 16 -s 16 connections each). This streams the large assemblies quickly.
#     Failures are tolerated (|| true); the paced resume below clears whatever the
#     OSF rate limiter rejects, so the whole step still converges.
A2_INPUT="${TMP_DIR}/osf_vka93.aria.txt"
: > "${A2_INPUT}"
while IFS=$'\t' read -r rel_path expected_size dl_url; do
    [[ -n "${rel_path}" ]] || continue
    target="${OSF_DIR}/${rel_path}"
    target_dir="$(dirname "${target}")"
    target_name="$(basename "${target}")"
    mkdir -p "${target_dir}"

    if [[ -f "${target}" && "${expected_size}" -gt 0 ]]; then
        actual_size=$(stat -c\%s "${target}")
        if [[ "${actual_size}" -eq "${expected_size}" ]]; then
            echo "    [SKIP] ${rel_path}"
            continue
        fi
    fi

    printf '%s\n  dir=%s\n  out=%s\n' "${dl_url}" "${target_dir}" "${target_name}" >> "${A2_INPUT}"
done < "${MANIFEST}"

aria2c -i "${A2_INPUT}" -c -j 12 -x 16 -s 16 \
    --user-agent="Mozilla/5.0" \
    --retry-wait=5 --max-tries=5 --auto-file-renaming=false \
    --file-allocation=none --console-log-level=warn --summary-interval=15 || true

# (c) Rate-paced resume: OSF 429s the high-rate bulk on the long tail of tiny files,
#     so finish the remainder with a few curl workers that sleep between requests and
#     back off on 429. Re-runs after a cooldown until nothing is left (capped passes).
OSF_PENDING="${TMP_DIR}/osf_pending.tsv"
for pass in $(seq 1 "${OSF_MAX_PASSES}"); do
    remaining=$(osf_count_pending)
    if (( remaining == 0 )); then
        echo "    OSF vka93 complete."
        break
    fi
    queued=$(osf_build_pending)
    echo "    [resume pass ${pass}] ${queued} file(s) @ ${OSF_WORKERS} curl workers (429 back-off)"
    while IFS= read -r line; do
        [[ -n "${line}" ]] || continue
        osf_fetch_one "${line}" &
        while [[ "$(jobs -rp | wc -l)" -ge "${OSF_WORKERS}" ]]; do wait -n || true; done
    done < "${OSF_PENDING}"
    wait || true
    remaining=$(osf_count_pending)
    echo "    [resume pass ${pass}] done; remaining=${remaining}"
    if (( remaining == 0 )); then
        echo "    OSF vka93 complete."
        break
    fi
    echo "    cooldown ${OSF_COOLDOWN}s before next pass..."
    sleep "${OSF_COOLDOWN}"
done

rm -rf "${TMP_DIR}"

echo "=============================================================================="
echo "All automated + supplementary downloads complete!"
echo "--- P. parvum strain comparison ---"
cat "${DATA_DIR}/genomes/catalog/ncbi_prymnesium_parvum.tsv"
echo "--- CNS00066 (NGDC GWH) ---"
ls -lh "${GWH_DIR}"
echo "--- PKZILLA Truth ---"
find "${PKZ_DIR}" -maxdepth 2 -not -name ".done"
echo "--- OSF vka93 ---"
echo "files: $(find "${OSF_DIR}" -type f ! -name '*manifest*' ! -name '*.aria2' | wc -l)"; du -sh "${OSF_DIR}"
echo "=============================================================================="

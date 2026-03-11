#!/bin/bash
# =============================================================================
# feedback3_blast_to_braker.sh
#
# FEEDBACK LOOP 3: BLAST hits → BRAKER re-annotation with enriched hints
#
# WORKFLOW:
#   1. Parse the BLAST results (blast_results.txt) produced by Step 7
#   2. Extract the subject sequences of all high-confidence hits
#      (e-value <= evalue_cutoff, pident >= pident_min) using efetch
#      or a local BLAST database if available
#   3. De-duplicate and combine these sequences with the original OrthoDB
#      protein hint file used in Step 3
#   4. Re-run BRAKER using the enriched protein hint set on the same
#      masked genome, placing output in a new braker_enriched/ directory
#   5. Re-run Step 4 protein extraction on the new annotation
#   6. Compare gene counts and mean protein length between the original
#      and enriched annotations; write a diff report
#
# WHY THIS HELPS:
#   The initial OrthoDB hints cover the organism's broader taxonomic group.
#   BLAST hits from Step 7 reveal which specific protein families are actually
#   present in the designed sequences — and those families may be poorly
#   represented in the generic OrthoDB partition.  Adding them as direct
#   hints can recover missed exons, fix gene merges/splits, and improve the
#   annotation of exactly the protein classes we care about designing.
#
# USAGE:
#   ./feedback3_blast_to_braker.sh [config.yaml|config.json|config.txt]
#   ./feedback3_blast_to_braker.sh \
#       --run_dir=./output/latest \
#       --evalue_cutoff=1e-10 \
#       --pident_min=40 \
#       --max_hits_per_query=5 \
#       --genemark_path=/path/to/genemark/bin
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir           Pipeline run directory (default: ./output/latest)
#   evalue_cutoff     Maximum e-value for a BLAST hit to be included (default: 1e-10)
#   pident_min        Minimum percent identity (default: 40)
#   max_hits_per_query Max subject sequences fetched per query (default: 5)
#   genemark_path     Path to GeneMark bin dir (default: ~/genemark-etp-full/...)
#   is_eukaryote      true/false, controls BRAKER vs Prokka (default: true)
#   bam_file          Path to RNA-Seq BAM to carry forward into re-annotation
#   output_dir        Where to write loop outputs (default: <run_dir>/feedback3_loop)
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# Shared helpers (mirror genome_to_design.sh style)
# ─────────────────────────────────────────────
log() { echo "[$(date '+%F %T')] [$1] ${*:2}" | tee -a "$master_log"; }

run_step() {
    local name=$1 critical=$2 desc=$3 fn=$4
    local done_f="$checkpoint_dir/.${name}.done"
    local fail_f="$checkpoint_dir/.${name}.failed"
    local step_log="$log_dir/${name}.log"
    if [ -f "$done_f" ]; then
        log SKIP "$name — already done. Delete $done_f to re-run."
        return 0
    fi
    log START "$name — $desc"
    local t0; t0=$(date +%s)
    "$fn" 2>&1 | tee -a "$step_log" | while IFS= read -r line; do
        echo "[$(date '+%F %T')] [$name] $line" >> "$master_log"
    done
    local rc=${PIPESTATUS[0]} elapsed=$(( $(date +%s) - t0 ))
    if [ "$rc" -eq 0 ]; then
        echo "Completed $(date '+%F %T') | elapsed ${elapsed}s" > "$done_f"
        rm -f "$fail_f"
        log DONE "$name — ${elapsed}s"
    else
        echo "Failed $(date '+%F %T') | exit=$rc | elapsed ${elapsed}s" > "$fail_f"
        log FAIL "$name — exit $rc. Log: $step_log"
        [ "$critical" = true ] && { log HALT "$name is critical — stopping."; exit "$rc"; }
        log WARN "$name is non-critical — continuing."
        return "$rc"
    fi
}

# ─────────────────────────────────────────────
# Config parsing — supports yaml, json, txt, CLI flags
# ─────────────────────────────────────────────
parse_config() {
    local file=$1 ext="${1##*.}"
    case "$ext" in
        yaml|yml) python3 -c "
import sys, yaml
cfg = yaml.safe_load(open('$file')) or {}
for k, v in cfg.items(): print(f'CFG_{k.upper()}={v}')
" ;;
        json) python3 -c "
import sys, json
cfg = json.load(open('$file'))
for k, v in cfg.items(): print(f'CFG_{k.upper()}={v}')
" ;;
        *) grep -v '^\s*#' "$file" | grep ':' | while IFS=: read -r k v; do
               echo "CFG_$(echo "$k" | tr '[:lower:]' '[:upper:]' | tr -d ' ')=$(echo "$v" | sed 's/^ //')"
           done ;;
    esac
}

# ── Defaults ──
RUN_DIR=""
EVALUE_CUTOFF="1e-10"
PIDENT_MIN=40
MAX_HITS=5
GENEMARK_PATH="${HOME}/genemark-etp-full/gmetp_linux_64/bin"
IS_EUKARYOTE="true"
BAM_FILE=""
OUTPUT_DIR=""

# ── Load config file if first arg looks like a file ──
if [ $# -ge 1 ] && [ -f "$1" ]; then
    eval "$(parse_config "$1")"
    [ -n "${CFG_RUN_DIR:-}"         ] && RUN_DIR="$CFG_RUN_DIR"
    [ -n "${CFG_EVALUE_CUTOFF:-}"   ] && EVALUE_CUTOFF="$CFG_EVALUE_CUTOFF"
    [ -n "${CFG_PIDENT_MIN:-}"      ] && PIDENT_MIN="$CFG_PIDENT_MIN"
    [ -n "${CFG_MAX_HITS_PER_QUERY:-}" ] && MAX_HITS="$CFG_MAX_HITS_PER_QUERY"
    [ -n "${CFG_GENEMARK_PATH:-}"   ] && GENEMARK_PATH="$CFG_GENEMARK_PATH"
    [ -n "${CFG_IS_EUKARYOTE:-}"    ] && IS_EUKARYOTE="$CFG_IS_EUKARYOTE"
    [ -n "${CFG_BAM_FILE:-}"        ] && BAM_FILE="$CFG_BAM_FILE"
    [ -n "${CFG_OUTPUT_DIR:-}"      ] && OUTPUT_DIR="$CFG_OUTPUT_DIR"
    shift
fi

# ── CLI flags override config ──
for arg in "$@"; do
    case "$arg" in
        --run_dir=*)           RUN_DIR="${arg#*=}" ;;
        --evalue_cutoff=*)     EVALUE_CUTOFF="${arg#*=}" ;;
        --pident_min=*)        PIDENT_MIN="${arg#*=}" ;;
        --max_hits_per_query=*) MAX_HITS="${arg#*=}" ;;
        --genemark_path=*)     GENEMARK_PATH="${arg#*=}" ;;
        --is_eukaryote=*)      IS_EUKARYOTE="${arg#*=}" ;;
        --bam=*)               BAM_FILE="${arg#*=}" ;;
        --output_dir=*)        OUTPUT_DIR="${arg#*=}" ;;
    esac
done

# ── Resolve paths ──
if [ -z "$RUN_DIR" ]; then
    RUN_DIR=$(realpath ./output/latest 2>/dev/null || echo "")
fi
[ -d "$RUN_DIR" ] || { echo "ERROR: run_dir '$RUN_DIR' not found."; exit 1; }
RUN_DIR=$(realpath "$RUN_DIR")
OUTPUT_DIR="${OUTPUT_DIR:-${RUN_DIR}/feedback3_loop}"
mkdir -p "$OUTPUT_DIR"
checkpoint_dir="$OUTPUT_DIR"
log_dir="$OUTPUT_DIR/logs"
mkdir -p "$log_dir"
master_log="$log_dir/feedback3.log"
num_threads=$(nproc)

log INFO "=========================================="
log INFO "Feedback Loop 3: BLAST hits → BRAKER re-annotation"
log INFO "run_dir        : $RUN_DIR"
log INFO "evalue_cutoff  : $EVALUE_CUTOFF"
log INFO "pident_min     : $PIDENT_MIN"
log INFO "max_hits       : $MAX_HITS"
log INFO "is_eukaryote   : $IS_EUKARYOTE"
log INFO "output_dir     : $OUTPUT_DIR"
log INFO "=========================================="

. ~/miniconda3/etc/profile.d/conda.sh

BLAST_RESULTS="$RUN_DIR/blast_results.txt"
[ -f "$BLAST_RESULTS" ] || { log HALT "blast_results.txt not found in $RUN_DIR"; exit 1; }

# ─────────────────────────────────────────────
# STEP FB3-1 — Extract subject accessions from BLAST results [CRITICAL]
# blast_results.txt is tabular fmt 6:
#   qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle
# We filter by evalue and pident, collect unique subject accessions,
# and write them to a file for efetch in the next step.
# ─────────────────────────────────────────────
ACCS_FILE="$OUTPUT_DIR/blast_subject_accs.txt"

_fb3_extract_accs() {
    conda activate braker_env || return 1
    python3 - <<PYEOF
import csv, math

blast_file   = "$BLAST_RESULTS"
out_file     = "$ACCS_FILE"
evalue_cut   = float("$EVALUE_CUTOFF")
pident_min   = float("$PIDENT_MIN")
max_per_q    = int("$MAX_HITS")

accs_by_query = {}
with open(blast_file) as f:
    for row in csv.reader(f, delimiter="\t"):
        if len(row) < 12:
            continue
        qseqid, sseqid = row[0], row[1]
        try:
            pident = float(row[2])
            evalue = float(row[10])
        except ValueError:
            continue
        if evalue <= evalue_cut and pident >= pident_min:
            accs_by_query.setdefault(qseqid, set()).add(sseqid)

all_accs = set()
for q, accs in accs_by_query.items():
    selected = list(accs)[:max_per_q]
    all_accs.update(selected)

with open(out_file, "w") as f:
    for acc in sorted(all_accs):
        # Strip version suffix for efetch compatibility (NP_123456.1 -> NP_123456)
        f.write(acc.split(".")[0] + "\n")

print(f"Queries with hits : {len(accs_by_query)}")
print(f"Unique accessions : {len(all_accs)}")
print(f"Written to        : {out_file}")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb3_extract_accs" true \
    "Extract subject accessions from BLAST results (evalue<=$EVALUE_CUTOFF, pident>=$PIDENT_MIN)" \
    _fb3_extract_accs

[ -s "$ACCS_FILE" ] || { log HALT "No accessions extracted — check BLAST results and thresholds"; exit 1; }

# ─────────────────────────────────────────────
# STEP FB3-2 — Fetch protein sequences for all accessions [CRITICAL]
# Uses NCBI Entrez efetch in batches of 200 to avoid rate limiting.
# Outputs to blast_hint_proteins.faa in the output directory.
# ─────────────────────────────────────────────
BLAST_HINTS_FAA="$OUTPUT_DIR/blast_hint_proteins.faa"

_fb3_fetch_sequences() {
    conda activate braker_env || return 1
    python3 - <<PYEOF
from Bio import Entrez, SeqIO
import time, math

Entrez.email = "PricklyPearEnterprises@gmail.com"

accs_file = "$ACCS_FILE"
out_file  = "$BLAST_HINTS_FAA"

with open(accs_file) as f:
    accs = [line.strip() for line in f if line.strip()]

print(f"Fetching {len(accs)} sequences from NCBI...")
batch_size = 200
records    = []
n_batches  = math.ceil(len(accs) / batch_size)

for i in range(0, len(accs), batch_size):
    batch = accs[i:i + batch_size]
    batch_num = i // batch_size + 1
    print(f"Batch {batch_num}/{n_batches}: {len(batch)} accessions")
    try:
        handle = Entrez.efetch(
            db="protein", id=",".join(batch),
            rettype="fasta", retmode="text"
        )
        batch_records = list(SeqIO.parse(handle, "fasta"))
        records.extend(batch_records)
        print(f"  Fetched {len(batch_records)} sequences")
    except Exception as e:
        print(f"  WARNING: efetch failed for batch {batch_num}: {e}")
    # Respect NCBI rate limit (3 requests/sec without API key)
    time.sleep(0.4)

# Sanitize headers (spaces/pipes → underscores, matching pipeline convention)
for rec in records:
    rec.id          = rec.id.replace(" ", "_").replace("|", "_")
    rec.description = ""

SeqIO.write(records, out_file, "fasta")
print(f"Total sequences fetched : {len(records)}")
print(f"Written to              : {out_file}")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb3_fetch_sequences" true \
    "Fetch protein sequences for BLAST hit accessions via Entrez" \
    _fb3_fetch_sequences

# ─────────────────────────────────────────────
# STEP FB3-3 — Merge BLAST hints with original OrthoDB hints [non-critical]
# Concatenate blast_hint_proteins.faa with the OrthoDB file used in Step 3.
# De-duplicate by sequence ID to avoid redundant hints inflating BRAKER's
# hint weight.  The merged file becomes the new --prot_seq for BRAKER.
# ─────────────────────────────────────────────
ENRICHED_HINTS="$OUTPUT_DIR/enriched_hints.faa"

_fb3_merge_hints() {
    conda activate braker_env || return 1
    # Find the OrthoDB file that was used — it lives under ~/orthodb_*.fasta
    local orthodb_file
    orthodb_file=$(ls "${HOME}"/orthodb_*.fasta 2>/dev/null | head -1)
    if [ -z "$orthodb_file" ]; then
        log WARN "No orthodb_*.fasta found in HOME — using BLAST hints only"
        cp "$BLAST_HINTS_FAA" "$ENRICHED_HINTS"
        conda deactivate
        return 0
    fi
    log INFO "Merging $BLAST_HINTS_FAA + $orthodb_file"
    python3 - <<PYEOF
from Bio import SeqIO

files    = ["$orthodb_file", "$BLAST_HINTS_FAA"]
seen_ids = set()
merged   = []
for fn in files:
    for rec in SeqIO.parse(fn, "fasta"):
        if rec.id not in seen_ids:
            seen_ids.add(rec.id)
            merged.append(rec)

SeqIO.write(merged, "$ENRICHED_HINTS", "fasta")
print(f"Merged {len(merged)} unique protein sequences -> $ENRICHED_HINTS")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb3_merge_hints" false \
    "Merge BLAST hints with original OrthoDB hints, de-duplicate" \
    _fb3_merge_hints

[ -f "$ENRICHED_HINTS" ] || { log HALT "Enriched hints file missing"; exit 1; }

# ─────────────────────────────────────────────
# STEP FB3-4 — BRAKER re-annotation with enriched hints [CRITICAL]
# Uses the same masked genome and GENEMARK_PATH as Step 3, but swaps
# in the enriched --prot_seq.  Outputs to braker_enriched/ under run_dir.
# ─────────────────────────────────────────────
BRAKER_ENRICHED_DIR="$RUN_DIR/braker_enriched"
# Resolve genome — prefer masked.fna, fall back to clean.fna
GENOME=""
for g in "$RUN_DIR/masked.fna" "$RUN_DIR/clean.fna"; do
    [ -f "$g" ] && { GENOME="$g"; break; }
done
[ -n "$GENOME" ] || { log HALT "No masked.fna or clean.fna found in $RUN_DIR"; exit 1; }

ANNO_PREFIX=$(basename "$RUN_DIR" | sed 's/run_[0-9_]*/enriched/')

_fb3_braker_enriched() {
    conda activate braker_env || return 1
    export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config"
    export AUGUSTUS_BIN_PATH="$CONDA_PREFIX/bin"
    export AUGUSTUS_SCRIPTS_PATH="$CONDA_PREFIX/bin"
    export GENEMARK_PATH="$GENEMARK_PATH"
    export PROTHINT_PATH="${GENEMARK_PATH}/ProtHint/bin"

    # Ensure Augustus config dir is writable
    if [ ! -w "$AUGUSTUS_CONFIG_PATH/species" ]; then
        cp -r "$AUGUSTUS_CONFIG_PATH" "${HOME}/augustus_config_writable"
        export AUGUSTUS_CONFIG_PATH="${HOME}/augustus_config_writable"
    fi
    # Remove any stale species dir
    rm -rf "$AUGUSTUS_CONFIG_PATH/species/$ANNO_PREFIX"
    mkdir -p "$BRAKER_ENRICHED_DIR"

    local braker_cmd="$CONDA_PREFIX/bin/braker.pl \
        --AUGUSTUS_CONFIG_PATH=\"$AUGUSTUS_CONFIG_PATH\" \
        --AUGUSTUS_BIN_PATH=\"$AUGUSTUS_BIN_PATH\" \
        --AUGUSTUS_SCRIPTS_PATH=\"$AUGUSTUS_SCRIPTS_PATH\" \
        --GENEMARK_PATH=\"$GENEMARK_PATH\" \
        --PROTHINT_PATH=\"$GENEMARK_PATH/ProtHint/bin\" \
        --genome=\"$GENOME\" \
        --species=\"$ANNO_PREFIX\" \
        --prot_seq=\"$ENRICHED_HINTS\" \
        --workingdir=\"$BRAKER_ENRICHED_DIR\" \
        --threads=$num_threads \
        --softmasking \
        --nocleanup \
        --augustus_args=\"--min_intron_len=10 --max_intron_len=25000\""

    # Carry forward RNA-Seq BAM if one was used originally or passed via config
    if [ -n "$BAM_FILE" ] && [ -f "$BAM_FILE" ]; then
        braker_cmd="$braker_cmd --bam=\"$BAM_FILE\""
    elif [ -f "$RUN_DIR/rnaseq_merged.bam" ]; then
        braker_cmd="$braker_cmd --bam=\"$RUN_DIR/rnaseq_merged.bam\""
    fi

    eval "$braker_cmd"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "BRAKER enriched failed — last 50 lines of braker.log:"
        tail -50 "$BRAKER_ENRICHED_DIR/braker.log" 2>/dev/null
    fi
    conda deactivate
    return $rc
}

run_step "fb3_braker_enriched" true \
    "BRAKER re-annotation with enriched protein hints" \
    _fb3_braker_enriched

# ─────────────────────────────────────────────
# STEP FB3-5 — Protein extraction on enriched annotation [CRITICAL]
# Mirrors Step 4 of genome_to_design.sh exactly.
# ─────────────────────────────────────────────
ENRICHED_FAA="$OUTPUT_DIR/proteins_enriched.faa"

_fb3_protein_extract() {
    conda activate braker_env || return 1
    local src_faa="$BRAKER_ENRICHED_DIR/${ANNO_PREFIX}.faa"
    [ -f "$src_faa" ] || { echo "ERROR: $src_faa not found"; return 1; }

    python3 - <<PYEOF
from Bio import SeqIO
import re

faa_file = "$src_faa"
out_file = "$ENRICHED_FAA"
min_len  = 200  # eukaryote default; lower to 100 for prokaryotes

all_records = list(SeqIO.parse(faa_file, "fasta"))
by_gene = {}
for r in all_records:
    m = re.match(r'^(.+)\.t\d+$', r.id)
    gene_id = m.group(1) if m else r.id
    if gene_id not in by_gene or len(r.seq) > len(by_gene[gene_id].seq):
        by_gene[gene_id] = r

proteins = sorted(
    [r for r in by_gene.values() if len(r.seq) > min_len],
    key=lambda r: r.id
)
SeqIO.write(proteins, out_file, "fasta")
print(f"Input isoforms   : {len(all_records)}")
print(f"Unique genes     : {len(by_gene)}")
print(f"Written (>{min_len}aa): {len(proteins)} -> {out_file}")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb3_protein_extract" true \
    "Extract and filter proteins from enriched annotation" \
    _fb3_protein_extract

# ─────────────────────────────────────────────
# STEP FB3-6 — Diff report: original vs enriched annotation [non-critical]
# ─────────────────────────────────────────────
DIFF_REPORT="$OUTPUT_DIR/annotation_diff.txt"
ORIG_FAA="$RUN_DIR/proteins.faa"

_fb3_diff_report() {
    conda activate braker_env || return 1
    python3 - <<PYEOF
from Bio import SeqIO

def stats(faa):
    recs = list(SeqIO.parse(faa, "fasta"))
    lengths = [len(r.seq) for r in recs]
    return {
        "count":    len(recs),
        "mean_len": sum(lengths) / len(lengths) if lengths else 0,
        "min_len":  min(lengths) if lengths else 0,
        "max_len":  max(lengths) if lengths else 0,
    }

orig = stats("$ORIG_FAA")
enr  = stats("$ENRICHED_FAA")

report = [
    "Annotation comparison: original vs BLAST-enriched hints",
    "=" * 56,
    f"{'Metric':<25} {'Original':>12} {'Enriched':>12} {'Delta':>10}",
    "-" * 56,
    f"{'Protein count':<25} {orig['count']:>12} {enr['count']:>12} {enr['count'] - orig['count']:>+10}",
    f"{'Mean length (aa)':<25} {orig['mean_len']:>12.1f} {enr['mean_len']:>12.1f} {enr['mean_len'] - orig['mean_len']:>+10.1f}",
    f"{'Min length (aa)':<25} {orig['min_len']:>12} {enr['min_len']:>12} {enr['min_len'] - orig['min_len']:>+10}",
    f"{'Max length (aa)':<25} {orig['max_len']:>12} {enr['max_len']:>12} {enr['max_len'] - orig['max_len']:>+10}",
    "=" * 56,
]

text = "\n".join(report)
print(text)
with open("$DIFF_REPORT", "w") as f:
    f.write(text + "\n")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb3_diff_report" false \
    "Compare original vs enriched annotation protein sets" \
    _fb3_diff_report

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────
log INFO "=========================================="
log INFO "Feedback Loop 3 complete."
log INFO "Enriched proteins  : $ENRICHED_FAA"
log INFO "Diff report        : $DIFF_REPORT"
log INFO "Enriched BRAKER    : $BRAKER_ENRICHED_DIR"
log INFO ""
log INFO "To continue the design pipeline from the enriched proteins:"
log INFO "  Copy $ENRICHED_FAA to $RUN_DIR/proteins.faa"
log INFO "  Then re-run genome_to_design.sh from step 5 (rm .step5_rfdiffusion.done)"
log INFO "=========================================="
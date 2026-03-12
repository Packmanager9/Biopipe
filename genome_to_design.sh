#!/bin/bash
# Usage: ./genome_to_design.sh "organism_name" output_dir [eukaryote] [--GENEMARK_PATH=/path] [--bam=/path/to/bam or comma-separated] [--auto_rnaseq] [--force]
# Example: ./genome_to_design.sh "Taraxacum officinale" ./output true --GENEMARK_PATH=${HOME}/genemark-etp-full/gmetp_linux_64/bin/gmes --auto_rnaseq
# Requires: Conda environments 'braker_env' (annotation, now with repeatmodeler repeatmasker), 'SE3nv' (RFdiffusion), 'colabfold' (structure prediction)
# For --auto_rnaseq: Needs sra-tools, star, samtools installed in braker_env.
#
# RESUMABILITY:
# Re-running the script will skip any step that already completed successfully.
# To force a single step to re-run: rm $run_dir/.step4_protein_extract.done
# To restart the entire pipeline: pass --force flag
# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────

# Logging helper (your original, assuming it's before this)
log() {
  local level="$1"; shift
  echo "[$(date '+%F %T')] [$level] $*" | tee -a "$master_log"
}

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate braker_env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

organism="$1"
out_dir="$2"
is_eukaryote="${3:-false}"
bam_file=""
auto_rnaseq=false
force_rerun=false
GENEMARK_PATH=""  # Initialize empty to allow default later

for arg in "$@"; do
  case "$arg" in
    --GENEMARK_PATH=*) GENEMARK_PATH="${arg#*=}" ;;
    --bam=*) bam_file="${arg#*=}" ;;
    --auto_rnaseq) auto_rnaseq=true ;;
    --force) force_rerun=true ;;
  esac
done

# Set default GENEMARK_PATH if not provided via flag
GENEMARK_PATH="${GENEMARK_PATH:-${HOME}/genemark-etp-full/gmetp_linux_64/bin}"
export GENEMARK_PATH  # Export for use in subshells

out_dir=$(realpath "$out_dir") # Resolve to absolute path to avoid relative dir issues

# ─────────────────────────────────────────────
# Detect number of available threads for parallelism
# ─────────────────────────────────────────────
num_threads=$(nproc) # Dynamically set to number of CPU threads (24 on this machine)
echo "INFO: Detected $num_threads available CPU threads for parallel processing."
# ─────────────────────────────────────────────
# Directory and log setup
# ─────────────────────────────────────────────
mkdir -p "$out_dir" || { echo "Error: Cannot create $out_dir"; exit 1; }
latest_symlink="$out_dir/latest"
if [ "$force_rerun" = true ]; then
timestamp=$(date +%Y%m%d_%H%M)
run_dir="$out_dir/run_$timestamp"
mkdir -p "$run_dir"
cd "$run_dir" || exit 1
checkpoint_dir="$(pwd)"
else
if [ -d "$latest_symlink" ]; then
run_dir=$(realpath "$latest_symlink")
cd "$run_dir" || exit 1
checkpoint_dir="$(pwd)"
if [ -f "$checkpoint_dir/.pipeline_complete.done" ]; then
timestamp=$(date +%Y%m%d_%H%M)
run_dir="$out_dir/run_$timestamp"
mkdir -p "$run_dir"
cd "$run_dir" || exit 1
checkpoint_dir="$(pwd)"
else
log "INFO" "Resuming incomplete run in $run_dir"
fi
else
timestamp=$(date +%Y%m%d_%H%M)
run_dir="$out_dir/run_$timestamp"
mkdir -p "$run_dir"
cd "$run_dir" || exit 1
checkpoint_dir="$(pwd)"
fi
fi
# Check if organism changed
organism_file="$checkpoint_dir/.pipeline_organism"
if [ -f "$organism_file" ] && [ "$(cat "$organism_file")" != "$organism" ]; then
log "INFO" "Organism changed from $(cat "$organism_file") to $organism: starting new run"
timestamp=$(date +%Y%m%d_%H%M)
run_dir="$out_dir/run_$timestamp"
mkdir -p "$run_dir"
cd "$run_dir" || exit 1
checkpoint_dir="$(pwd)"
organism_file="$checkpoint_dir/.pipeline_organism"
fi
# Clear all checkpoints if --force
if [ "$force_rerun" = true ]; then
echo "INFO: --force passed: clearing all checkpoints"
rm -f "$checkpoint_dir"/.step*.done "$checkpoint_dir"/.step*.failed
fi
# Save current organism
echo "$organism" > "$organism_file"
log_dir="$(pwd)/logs"
mkdir -p "$log_dir"
master_log="$log_dir/pipeline.log"
# ─────────────────────────────────────────────
# Progress bar helper
# ─────────────────────────────────────────────
progress() {
local current=$1 total=$2 last_elapsed=$3 avg=$4 remaining=$5
local percent=$((current * 100 / total))
local bar_length=50
local filled=$((percent * bar_length / 100))
local bar="\r["
for ((i=0; i<filled; i++)); do bar+="="; done
for ((i=filled; i<bar_length; i++)); do bar+=" "; done
bar+="] %d%% (%d/%d) | Last: %ds | Avg: %ds | ETA: %ds"
printf "$bar" $percent $current $total $last_elapsed $avg $remaining > /dev/tty
}
# ─────────────────────────────────────────────
# Checkpoint runner
#
# run_step STEP_NAME IS_CRITICAL "Description" shell_function
#
# IS_CRITICAL=true → failure halts the pipeline immediately
# IS_CRITICAL=false → failure is logged as a warning; pipeline continues
#
# Outputs:
# logs/<step_name>.log per-step stdout+stderr
# .step_name.done written on success (contains timestamp + elapsed)
# .step_name.failed written on failure (contains exit code + timestamp)
# ─────────────────────────────────────────────
run_step() {
local step_name="$1"
local is_critical="$2"
local description="$3"
local fn="$4"
local done_file="$checkpoint_dir/.${step_name}.done"
local fail_file="$checkpoint_dir/.${step_name}.failed"
local step_log="$log_dir/${step_name}.log"
# ── Skip if already done ──
if [ -f "$done_file" ]; then
log "SKIP" "$step_name — already completed ($(cat "$done_file"))"
log "SKIP" " Delete $done_file to re-run this step."
return 0
fi
# ── Run ──
log "START" "$step_name — $description"
local start_ts
start_ts=$(date +%s)
# Capture output to step log; also mirror each line into master log with prefix
"$fn" 2>&1 \
    | tee -a "$step_log" \
    | tr '\r' '\n' \
    | grep -v '^[[:space:]]*$' \
    | while IFS= read -r line; do
        echo "[$(date '+%F %T')] [$step_name] $line" >> "$master_log"
      done
local exit_code="${PIPESTATUS[0]}"
local elapsed=$(( $(date +%s) - start_ts ))
if [ "$exit_code" -eq 0 ]; then
echo "Completed $(date '+%F %T') | elapsed ${elapsed}s" > "$done_file"
rm -f "$fail_file"
log "DONE" "$step_name — finished in ${elapsed}s"
return 0
else
echo "Failed $(date '+%F %T') | exit_code=$exit_code | elapsed ${elapsed}s" > "$fail_file"
log "FAIL" "$step_name — exit code $exit_code after ${elapsed}s. Full log: $step_log"
if [ "$is_critical" = true ]; then
log "HALT" "$step_name is critical — pipeline cannot continue."
log "HALT" "Fix the error, then re-run (the step will resume from here)."
exit "$exit_code"
else
log "WARN" "$step_name is non-critical — continuing to next step."
return "$exit_code"
fi
fi
}
# ─────────────────────────────────────────────
# Conda setup
# ─────────────────────────────────────────────
. ~/miniconda3/etc/profile.d/conda.sh
email="Example@mail.com"
anno_prefix=$(echo "$organism" | sed 's/ /_/g')
log "INFO" "Pipeline started | organism: $organism | eukaryote: $is_eukaryote | force: $force_rerun"
log "INFO" "Working dir: $(pwd)"
# ─────────────────────────────────────────────
# STEP 1 — Genome fetch [CRITICAL]
# ─────────────────────────────────────────────
# _step1_genome_fetch() { 
# # Wrap genome_fetch with 'script' to simulate a TTY, preventing hangs in non-interactive environments
# # Assumes 'script' (from util-linux) is available; if not, install via apt/yum or conda install -c conda-forge util-linux
# log "DEBUG" "Running genome_fetch with simulated TTY to avoid potential hangs"
# script -q -c "${HOME}/Desktop/favbooks/genome_fetch.sh \"$organism\"" /dev/null # Use absolute path here
# }
_step1_genome_fetch() {
    # Resolve genome_fetch.sh: prefer /usr/local/bin (installed), then same dir as this script
    local fetch_script
    if command -v genome_fetch.sh >/dev/null 2>&1; then
        fetch_script="genome_fetch.sh"
    else
        local script_dir
        script_dir="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
        fetch_script="$script_dir/genome_fetch.sh"
    fi
    if [ ! -f "$fetch_script" ] && ! command -v "$fetch_script" >/dev/null 2>&1; then
        log "HALT" "genome_fetch.sh not found in PATH or alongside genome_to_design.sh"
        return 1
    fi
    log "DEBUG" "Running genome_fetch with simulated TTY to avoid potential hangs"
    if command -v script >/dev/null 2>&1; then
        script -q -c "$fetch_script \"$organism\"" /dev/null
    else
        log "WARN" "'script' not found — running genome_fetch.sh directly (no TTY simulation)"
        bash "$fetch_script" "$organism"
    fi
} 
run_step "step1_genome_fetch" true "Fetch reference genome for '$organism'" _step1_genome_fetch
# Decompress — runs every time (idempotent, fast, not worth checkpointing)
# Note: gunzip is standard bash, but keeping conda activate if bioenv has specific version/dependencies
conda activate braker_env
find . -type f -name "*.gz" -exec gunzip -f {} \; 2>>"$log_dir/step1_genome_fetch.log" \
|| log "WARN" "gunzip failed for some files — may already be decompressed"
conda deactivate
# Broader search for organism_dir (up to maxdepth 5, to handle deeper nesting like refseq/plant/GCF_xxx or human_readable variants)
organism_pattern="*$(echo "$organism" | sed 's/ /_/g')*"
organism_dir=$(find . -maxdepth 5 -type d -iname "$organism_pattern" | head -1)
# If found, flatten structure to root-level organism_dir for consistency (use exact organism name without _dir suffix)
if [ -n "$organism_dir" ]; then
    new_dir="./$(echo "$organism" | sed 's/ /_/g')"
    mkdir -p "$new_dir"
    find "$organism_dir" -type f \( -name "*.fna" -o -name "*.fasta" \) -exec mv {} "$new_dir/" \;
    organism_dir="$new_dir"
    log "INFO" "Flattened genome files to $organism_dir"
fi
# Fallback: If no matching dir, search for any dir with .fna/.fasta and assume it's the genome dir
if [ -z "$organism_dir" ]; then
    organism_dir=$(find . -maxdepth 5 -type f \( -name "*.fna" -o -name "*.fasta" \) | head -1 | xargs dirname)
    if [ -n "$organism_dir" ]; then
        log "WARN" "No exact organism dir match – using fallback dir: $organism_dir"
        # Flatten as above
        new_dir="./$(echo "$organism" | sed 's/ /_/g')"
        mkdir -p "$new_dir"
        find "$organism_dir" -type f \( -name "*.fna" -o -name "*.fasta" \) -exec mv {} "$new_dir/" \;
        organism_dir="$new_dir"
        log "INFO" "Flattened fallback genome files to $organism_dir"
    fi
fi
# Existing check/re-run (keep, but now with broader find)
if [ -z "$organism_dir" ] && [ -f "$checkpoint_dir/.step1_genome_fetch.done" ]; then
    log "WARN" "Genome directory not found after skipping step1_genome_fetch. Forcing fresh re-run of step1."
    rm -f "$checkpoint_dir/.step1_genome_fetch.done" "$checkpoint_dir/.step1_genome_fetch.failed"
    run_step "step1_genome_fetch" true "Fetch reference genome for '$organism' (forced re-run)" _step1_genome_fetch
    # Re-decompress if needed
    conda activate braker_env
    find . -type f -name "*.gz" -exec gunzip -f {} \; 2>>"$log_dir/step1_genome_fetch.log" \
    || log "WARN" "gunzip failed for some files — may already be decompressed"
    conda deactivate
    # Re-find the directory
    organism_dir=$(find . -maxdepth 5 -type d -iname "$organism_pattern" | head -1)
    # Re-flatten if found
    if [ -n "$organism_dir" ]; then
        new_dir="./$(echo "$organism" | sed 's/ /_/g')"
        mkdir -p "$new_dir"
        find "$organism_dir" -type f \( -name "*.fna" -o -name "*.fasta" \) -exec mv {} "$new_dir/" \;
        organism_dir="$new_dir"
        log "INFO" "Flattened genome files to $organism_dir"
    fi
    # Re-fallback if still no match
    if [ -z "$organism_dir" ]; then
        organism_dir=$(find . -maxdepth 5 -type f \( -name "*.fna" -o -name "*.fasta" \) | head -1 | xargs dirname)
        if [ -n "$organism_dir" ]; then
            log "WARN" "No exact organism dir match after re-fetch – using fallback dir: $organism_dir"
            new_dir="./$(echo "$organism" | sed 's/ /_/g')"
            mkdir -p "$new_dir"
            find "$organism_dir" -type f \( -name "*.fna" -o -name "*.fasta" \) -exec mv {} "$new_dir/" \;
            organism_dir="$new_dir"
            log "INFO" "Flattened fallback genome files to $organism_dir"
        fi
    fi
fi
# Now proceed with genome_file check
if [ -z "$organism_dir" ]; then
log "HALT" "No genome directory found even after attempting re-fetch — cannot continue"
exit 1
fi
# Broader genome_file search (recursive under organism_dir)
genome_file=$(find "$organism_dir" -type f \( -name "*.fna" -o -name "*.fasta" \) | head -1)
[ -f "$genome_file" ] || { log "HALT" "No genome file found in $organism_dir — cannot continue"; exit 1; }
log "INFO" "Genome file resolved: $genome_file"
# Compute genome length for STAR --genomeSAindexNbases
genome_len=$(awk '!/^>/ {len += length($0)} END {print len}' "$genome_file")
sa_index=$(python3 -c "import math; l = $genome_len if $genome_len > 0 else 1; val = math.log2(l) / 2 - 1; print(min(14, int(val)))")
log "INFO" "Computed genome length: $genome_len bp, STAR genomeSAindexNbases: $sa_index"
# ─────────────────────────────────────────────
# STEP 1b — Sanitize FASTA headers [non-critical: falls back to original]
# ─────────────────────────────────────────────
_step1b_sanitize_headers() {
conda activate braker_env || return 1
local sanitized_genome="${genome_file%.*}_sanitized.fna"
python3 -c "
from Bio import SeqIO
import sys
in_file = '$genome_file'
out_file = '$sanitized_genome'
records = []
for rec in SeqIO.parse(in_file, 'fasta'):
    rec.id = rec.id.replace(' ', '_').replace('|', '_')
    rec.description = ''
    records.append(rec)
SeqIO.write(records, out_file, 'fasta')
print(f'Sanitized headers in {in_file} -> {out_file}')
" || { echo "Header sanitization failed"; return 1; }
local rc=$?
conda deactivate
return $rc
}
run_step "step1b_sanitize_headers" false "Sanitize FASTA headers in genome (replace spaces/| with _)" _step1b_sanitize_headers
if [ -f "${genome_file%.*}_sanitized.fna" ]; then
genome_file="${genome_file%.*}_sanitized.fna"
log "INFO" "Using sanitized genome: $genome_file"
fi
# ─────────────────────────────────────────────
# STEP 2 — QC with BBTools [non-critical: falls back to raw genome]
# ─────────────────────────────────────────────
_step2_qc() {
conda activate braker_env || return 1
bbduk.sh in="$genome_file" out=clean.fna ref=adapters stats=qc_stats.txt Xmx=8g
local rc=$?
conda deactivate
return $rc
}
run_step "step2_qc" false "QC genome with bbduk (adapter trimming)" _step2_qc
if [ -f "clean.fna" ]; then
clean_genome="$(pwd)/clean.fna"
log "INFO" "Using QC-cleaned genome: $clean_genome"
else
clean_genome="$(realpath "$genome_file")"
log "WARN" "clean.fna not found — falling back to raw genome: $clean_genome"
fi
# ─────────────────────────────────────────────
# STEP 2b — Repeat masking [non-critical: falls back to clean genome]
# ─────────────────────────────────────────────
_step2b_repeatmask() {
conda activate braker_env || return 1
mkdir -p repeat_out

# ── BuildDatabase ──────────────────────────────────────────────────────────
echo "[REPEATMASK] Building repeat database from $(basename "$clean_genome")"
BuildDatabase -name "$anno_prefix" "$clean_genome"
local db_rc=$?
if [ $db_rc -ne 0 ]; then
    echo "[REPEATMASK] BuildDatabase failed ($db_rc) — skipping masking"
    cp "$clean_genome" masked.fna
    conda deactivate
    return 0
fi

# ── RepeatModeler ──────────────────────────────────────────────────────────
# Tee full (unfiltered) output to a dedicated detail log so nothing is lost.
# Strip only the per-second countdown lines from stdout so pipeline.log stays
# readable while still showing round headers, errors and summary lines.
local rm_detail_log="$log_dir/repeatmodeler_detail.log"
echo "[REPEATMASK] Starting RepeatModeler (full log → $rm_detail_log)"
echo "[REPEATMASK] This step takes several hours on large genomes — round progress will appear below"

# RepeatModeler -database "$anno_prefix" -pa $num_threads 2>&1 | \
#     tee -a "$rm_detail_log" | \
#     grep -vE '^\s+[0-9]+%\s+completed,\s+[0-9]' | \
#     grep -vE '^\s*$'

RepeatModeler -database "$anno_prefix" -pa $num_threads 2>&1 | \
    tee -a "$rm_detail_log" | \
    perl -pe 's/\b(\d+):(\d+):(\d+)\b/sprintf("%02d:%02d:%02d",$1,$2,$3)/ge' | \
    grep -vE '^\s+[0-9]+%\s+completed,\s+[0-9]' | \
    grep -vE '^\s*$'


local rm_rc=${PIPESTATUS[0]}   # RepeatModeler exit — not grep's

if [ $rm_rc -ne 0 ]; then
    echo "[REPEATMASK] RepeatModeler exited $rm_rc"
    if [ "$is_eukaryote" = false ]; then
        echo "[REPEATMASK] Non-eukaryote: treating as no-repeats, continuing without masking"
        cp "$clean_genome" masked.fna
        conda deactivate
        return 0
    else
        echo "[REPEATMASK] Eukaryote: masking required — last 20 lines of detail log:"
        tail -20 "$rm_detail_log"
        conda deactivate
        return $rm_rc
    fi
fi

# ── Resolve repeat library ─────────────────────────────────────────────────
# [ -f "RM_*/glob" ] does NOT expand in bash — must use find.
local repeat_lib
repeat_lib=$(find . -maxdepth 2 -name "consensi.fa.classified" | sort | tail -1)
if [ -z "$repeat_lib" ] || [ ! -s "$repeat_lib" ]; then
    # Fall back to unclassified consensi if classified was not produced
    repeat_lib=$(find . -maxdepth 2 -name "consensi.fa" | sort | tail -1)
fi
if [ -z "$repeat_lib" ] || [ ! -s "$repeat_lib" ]; then
    echo "[REPEATMASK] RepeatModeler produced no repeat library (consensi.fa is empty or missing)"
    echo "[REPEATMASK] Continuing without masking"
    cp "$clean_genome" masked.fna
    conda deactivate
    return 0
fi
echo "[REPEATMASK] Using repeat library: $repeat_lib ($(wc -l < "$repeat_lib") lines)"

# ── RepeatMasker ───────────────────────────────────────────────────────────
echo "[REPEATMASK] Starting RepeatMasker"
RepeatMasker -pa $num_threads -lib "$repeat_lib" -xsmall -dir repeat_out "$clean_genome"
local rc=$?
if [ $rc -eq 0 ] && [ -f "repeat_out/clean.fna.masked" ]; then
    mv "repeat_out/clean.fna.masked" masked.fna
    local masked_pct
    masked_pct=$(grep -i 'bases masked' repeat_out/clean.fna.tbl 2>/dev/null | grep -oE '[0-9]+\.[0-9]+%' | head -1)
    echo "[REPEATMASK] RepeatMasker complete — masked ${masked_pct:-unknown} of genome"
elif [ $rc -eq 0 ]; then
    echo "[REPEATMASK] RepeatMasker succeeded but masked file not found — falling back to clean genome"
    cp "$clean_genome" masked.fna
fi
conda deactivate
return $rc
}
repeat_critical=false
if [ "$is_eukaryote" = true ]; then repeat_critical=true; fi
# FIX: Skip repeat masking if genome is too small
if [ $genome_len -lt 1000000 ]; then # Arbitrary threshold for "small" genome (1 MB)
log "WARN" "Genome too small ($genome_len bp) for reliable repeat masking - skipping step2b"
cp "$clean_genome" masked.fna
touch "$checkpoint_dir/.step2b_repeatmask.done"
anno_genome="masked.fna"
else
run_step "step2b_repeatmask" "$repeat_critical" "Mask repeats with RepeatModeler/RepeatMasker" _step2b_repeatmask
if [ -f "masked.fna" ]; then
anno_genome="masked.fna"
log "INFO" "Using masked genome for annotation: $anno_genome"
else
anno_genome="$clean_genome"
log "WARN" "masked.fna not found — falling back to $anno_genome"
fi
fi
# ─────────────────────────────────────────────
# STEP 3 — Genome annotation [CRITICAL]
# ─────────────────────────────────────────────
_step3_braker() {
  conda activate braker_env || return 1
  export AUGUSTUS_CONFIG_PATH="$CONDA_PREFIX/config" # Point to env's configs
  export AUGUSTUS_BIN_PATH="$CONDA_PREFIX/bin"
  export AUGUSTUS_SCRIPTS_PATH="$CONDA_PREFIX/bin"
  export BAMTOOLS_PATH="$CONDA_PREFIX/bin/bamtools"
  export PATH="$AUGUSTUS_BIN_PATH:$PATH"
  export GENEMARK_PATH="${GENEMARK_PATH:-${HOME}/genemark-etp-full/gmetp_linux_64/bin}"
  export PROTHINT_PATH="${GENEMARK_PATH}/ProtHint/bin" # Adjust if ProtHint is not in subdir
  echo "GENEMARK_PATH=$GENEMARK_PATH" # For debugging

  # Ensure Augustus config dir is writable
  if [ ! -w "$AUGUSTUS_CONFIG_PATH/species" ]; then
    echo "Augustus config not writable -- copying to ${HOME}/augustus_config_writable"
    cp -r "$AUGUSTUS_CONFIG_PATH" "${HOME}/augustus_config_writable"
    export AUGUSTUS_CONFIG_PATH="${HOME}/augustus_config_writable"
  fi

# Remove stale species dir so BRAKER trains fresh
local species_dir="$AUGUSTUS_CONFIG_PATH/species/$anno_prefix"
if [ -d "$species_dir" ]; then
echo "Removing stale Augustus species dir: $species_dir"
rm -rf "$species_dir"
fi
local braker_out="$(pwd)/braker_out"
mkdir -p "$braker_out"
# Auto-query and download RNA-Seq if --auto_rnaseq and no --bam
if [ "$auto_rnaseq" = true ] && [ -z "$bam_file" ]; then
echo "Auto-querying RNA-Seq for $organism..."
local srr_list=$(python3 - <<EOF
from Bio import Entrez
Entrez.email = "$email"
try:
    # Broaden strategy terms and remove date filter
    query = f'$organism[Organism] AND ("RNA-Seq"[Strategy] OR "transcriptome sequencing"[Strategy] OR "rna seq"[Strategy] OR "mRNA-Seq"[Strategy] OR "ncRNA-Seq"[Strategy])'
    handle = Entrez.esearch(db='sra', term=query, retmax=20, sort='relevance')
    record = Entrez.read(handle)
    ids = record['IdList']
    srrs = []
    if ids:
        for id in ids:
            handle = Entrez.efetch(db='sra', id=id, retmode='xml')
            xml = handle.read().decode('utf-8')
            import re
            matches = re.findall(r'<RUN acc="(SRR|ERR|SRA)\d+"/>', xml)
            srrs.extend([m.replace('<RUN acc="', '').replace('"/>', '') for m in matches])
    if not srrs:
        # Fallback to genus
        genus = '$organism'.split()[0]
        query_genus = f'{genus}[Organism] AND ("RNA-Seq"[Strategy] OR "transcriptome sequencing"[Strategy] OR "rna seq"[Strategy])'
        handle = Entrez.esearch(db='sra', term=query_genus, retmax=20, sort='relevance')
        record = Entrez.read(handle)
        ids = record['IdList']
        if ids:
            for id in ids:
                handle = Entrez.efetch(db='sra', id=id, retmode='xml')
                xml = handle.read().decode('utf-8')
                import re
                matches = re.findall(r'<RUN acc="(SRR|ERR|SRA)\d+"/>', xml)
                srrs.extend([m.replace('<RUN acc="', '').replace('"/>', '') for m in matches])
    print(' '.join(srrs[:10])) # Limit to first 10 SRR total for efficiency
except Exception as e:
    print('') # Fallback on error
EOF
)
if [ -n "$srr_list" ]; then
echo "Found SRA runs: $srr_list – downloading with prefetch (HTTPS transport)..."
for srr in $srr_list; do
prefetch --transport http -O . $srr || echo "Warning: prefetch failed for $srr – skipping"
done
local successful_srr=$(ls *.sra 2>/dev/null | sed 's/\.sra$//' | tr '\n' ' ')
if [ -n "$successful_srr" ]; then
srr_list="$successful_srr"
num_srr=$(echo "$srr_list" | wc -w)
# Fasterq-dump loop - parallelize with GNU parallel if available
if command -v parallel >/dev/null; then
echo "Using GNU parallel for fasterq-dump"
parallel -j $((num_threads / 2)) fasterq-dump --split-files {} ::: $srr_list
else
local i=0
local total_elapsed=0
for srr in $srr_list; do
i=$((i + 1))
echo "Downloading FASTQ $i/$num_srr: $srr"
local start_time=$(date +%s)
fasterq-dump --split-files $srr
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_srr - i)))
progress $i $num_srr $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed FASTQ download $i/$num_srr in ${elapsed}s"
if [ $i -lt $num_srr ]; then
echo "Estimated time remaining for downloads: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
# Build STAR index if not exists
if [ ! -d "star_index" ]; then
STAR --genomeDir star_index --runMode genomeGenerate --genomeFastaFiles "$anno_genome" --runThreadN $num_threads --genomeSAindexNbases $sa_index
fi
# STAR alignment loop - parallelize with GNU parallel if available
if command -v parallel >/dev/null; then
echo "Using GNU parallel for STAR alignments"
parallel -j $((num_threads / 4)) 'STAR --genomeDir star_index --readFilesIn {}_1.fastq {}_2.fastq --outFileNamePrefix {}_ --outSAMtype BAM SortedByCoordinate --runThreadN $((num_threads / 4))' ::: $srr_list
else
i=0
total_elapsed=0
for srr in $srr_list; do
i=$((i + 1))
echo "Aligning $i/$num_srr: $srr"
local start_time=$(date +%s)
if [ -f "${srr}_2.fastq" ]; then
STAR --genomeDir star_index --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outFileNamePrefix ${srr}_ --outSAMtype BAM SortedByCoordinate --runThreadN $num_threads
else
STAR --genomeDir star_index --readFilesIn ${srr}.fastq --outFileNamePrefix ${srr}_ --outSAMtype BAM SortedByCoordinate --runThreadN $num_threads
fi
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_srr - i)))
progress $i $num_srr $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed alignment $i/$num_srr in ${elapsed}s"
if [ $i -lt $num_srr ]; then
echo "Estimated time remaining for alignments: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
samtools merge -@ $num_threads rnaseq_merged.bam *_Aligned.sortedByCoord.out.bam
samtools index rnaseq_merged.bam
bam_file="$(pwd)/rnaseq_merged.bam"
else
echo "Warning: No SRA files downloaded – proceeding with proteins only"
fi
else
echo "Warning: No suitable RNA-Seq found – proceeding with proteins only"
fi
fi
# Determine OrthoDB group based on organism taxonomy using Biopython Entrez
local orthodb_group=$(python3 - <<EOF
from Bio import Entrez
Entrez.email = "$email" # Required for Entrez queries
try:
    handle = Entrez.esearch(db="taxonomy", term="$organism")
    record = Entrez.read(handle)
    taxid = record['IdList'][0] if record['IdList'] else None
    if taxid:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
        lineage = [tax['ScientificName'].lower() for tax in record[0]['LineageEx']]
        # Map to OrthoDB partition
        if 'fungi' in lineage:
            print('fungi')
        elif 'metazoa' in lineage:
            if 'vertebrata' in lineage:
                print('vertebrata')
            elif 'arthropoda' in lineage:
                print('arthropoda')
            else:
                print('metazoa')
        elif 'viridiplantae' in lineage:
            print('viridiplantae')
        elif 'alveolata' in lineage:
            print('alveolata')
        elif 'stramenopiles' in lineage:
            print('stramenopiles')
        elif 'amoebozoa' in lineage:
            print('amoebozoa')
        elif 'euglenozoa' in lineage:
            print('euglenozoa')
        else:
            print('eukaryota')
    else:
        print('eukaryota')
except Exception as e:
    print('eukaryota') # Fallback on error
EOF
) || { echo "Failed to determine OrthoDB group, falling back to eukaryota"; orthodb_group="eukaryota"; }
echo "Determined OrthoDB group: $orthodb_group"
# Set OrthoDB file path
local orthodb_file="$HOME/orthodb_${orthodb_group}.fasta"
# Download if missing
if [ ! -f "$orthodb_file" ]; then
local taxa_cap=$(echo "$orthodb_group" | sed 's/./\u&/')
echo "Downloading OrthoDB v12 ${orthodb_group} proteins for BRAKER hints from alternative source..."
echo "Determined OrthoDB group: $orthodb_group"
echo "Expected OrthoDB file: $HOME/orthodb_${orthodb_group}.fasta"
if [ -f "$HOME/orthodb_${orthodb_group}.fasta" ]; then
  echo "File found — skipping download"
else
  echo "File not found — will download"
fi
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/${taxa_cap}.fa.gz -P "${HOME}"
gunzip ~/${taxa_cap}.fa.gz || { echo "gunzip failed"; return 1; }
mv ~/${taxa_cap}.fa "${orthodb_file}.tmp"
# Sanitize OrthoDB headers (replace spaces/| with _)
python3 -c "
from Bio import SeqIO
in_file = '${orthodb_file}.tmp'
out_file = '$orthodb_file'
records = []
for rec in SeqIO.parse(in_file, 'fasta'):
    rec.id = rec.id.replace(' ', '_').replace('|', '_')
    rec.description = ''
    records.append(rec)
SeqIO.write(records, out_file, 'fasta')
print(f'Sanitized OrthoDB headers -> {out_file}')
" || { echo "OrthoDB sanitization failed"; rm "${orthodb_file}.tmp"; return 1; }
rm "${orthodb_file}.tmp"
fi
[ -f "$orthodb_file" ] || { echo "OrthoDB download failed"; return 1; }
# Check if the file is empty or too small
if [ ! -s "$orthodb_file" ]; then
echo "OrthoDB file is empty - redownloading..."
rm -f "$orthodb_file"
local taxa_cap=$(echo "$orthodb_group" | sed 's/./\u&/')
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/${taxa_cap}.fa.gz -P "${HOME}" #pattern matched?
gunzip ~/${taxa_cap}.fa.gz || { echo "gunzip failed"; return 1; }
mv ~/${taxa_cap}.fa "${orthodb_file}.tmp"
# Sanitize as above
python3 -c "
from Bio import SeqIO
in_file = '${orthodb_file}.tmp'
out_file = '$orthodb_file'
records = []
for rec in SeqIO.parse(in_file, 'fasta'):
    rec.id = rec.id.replace(' ', '_').replace('|', '_')
    rec.description = ''
    records.append(rec)
SeqIO.write(records, out_file, 'fasta')
print(f'Sanitized OrthoDB headers -> {out_file}')
" || { echo "OrthoDB sanitization failed"; rm "${orthodb_file}.tmp"; return 1; }
rm "${orthodb_file}.tmp"
  [ -s "$orthodb_file" ] || { echo "Redownload failed - OrthoDB file still empty"; return 1; }
fi
# --nocleanup prevents BRAKER deleting genome.fa (its header-sanitized genome
# copy). We need it for gffread so sequence IDs match the GTF headers.
braker_cmd="$CONDA_PREFIX/bin/braker.pl \
--AUGUSTUS_CONFIG_PATH=\"$AUGUSTUS_CONFIG_PATH\" \
--AUGUSTUS_BIN_PATH=\"$AUGUSTUS_BIN_PATH\" \
--AUGUSTUS_SCRIPTS_PATH=\"$AUGUSTUS_SCRIPTS_PATH\" \
--GENEMARK_PATH=\"$GENEMARK_PATH\" \
--PROTHINT_PATH=\"$GENEMARK_PATH/ProtHint/bin\" \
--genome=\"$anno_genome\" \
--species=\"$anno_prefix\" \
--prot_seq=\"$orthodb_file\" \
--workingdir=\"$braker_out\" \
--threads=$num_threads \
--softmasking \
--nocleanup \
--augustus_args=\"--min_intron_len=10 --max_intron_len=25000\""
if [ -n "$bam_file" ]; then
braker_cmd="$braker_cmd --bam=\"$bam_file\""
fi
eval "$braker_cmd"
local rc=$?
# FIX: Check for low hints and fallback to --esmode
if [ $rc -ne 0 ] && grep -q "less than 1000 introns" "$braker_out/braker.log"; then
log "WARN" "Low hints detected - retrying BRAKER in ES mode"
braker_cmd="${braker_cmd/--prot_seq=\"$orthodb_file\"/} --esmode"
eval "$braker_cmd"
rc=$?
fi
conda deactivate
if [ $rc -ne 0 ]; then
echo "BRAKER exited $rc -- last 50 lines of braker.log:"
tail -50 "$braker_out/braker.log" 2>/dev/null
return $rc
fi
local gtf_file="$braker_out/braker.gtf"
[ -f "$gtf_file" ] || { echo "braker.gtf not found in $braker_out"; return 1; }
# Resolve the genome file for gffread.
# BRAKER sanitizes FASTA headers (strips whitespace/pipes) and saves the
# result as genome.fa in --workingdir. The GTF uses those sanitized IDs.
# --nocleanup keeps genome.fa alive. Fall back to Python remap if missing.
local braker_genome="$braker_out/genome.fa"
if [ ! -f "$braker_genome" ]; then
echo "genome.fa not found -- attempting remap via genome_header.map"
local header_map="$braker_out/genome_header.map"
if [ -f "$header_map" ]; then
python3 -c "
map_file = '$header_map'
in_fasta = '$anno_genome'
out_fasta = '$braker_out/genome_remapped.fa'
mapping = {}
with open(map_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            mapping[parts[0]] = parts[1]
with open(in_fasta) as fin, open(out_fasta, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            orig = line[1:].split()[0]
            fout.write('>' + mapping.get(orig, orig) + '\n')
        else:
            fout.write(line)
print(f'Remapped {len(mapping)} headers -> {out_fasta}')
"
braker_genome="$braker_out/genome_remapped.fa"
else
echo "Warning: genome_header.map also missing -- falling back to anno_genome"
braker_genome="$anno_genome"
fi
else
echo "Using BRAKER genome.fa for gffread: $braker_genome"
fi
# Repair GTF: GeneMark sometimes emits records without transcript_id.
# gffread requires it -- inject from gene_id (e.g. g1 -> g1.t1).
local fixed_gtf="$braker_out/braker_fixed.gtf"
python3 -c "
import re
in_path = '$gtf_file'
out_path = '$fixed_gtf'
ok = repaired = 0
with open(in_path) as fin, open(out_path, 'w') as fout:
    for line in fin:
        if line.startswith('#') or not line.strip():
            fout.write(line); continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 9:
            fout.write(line); continue
        attrs = parts[8]
        if 'transcript_id' in attrs:
            ok += 1; fout.write(line); continue
        m = re.search(r'gene_id\s+\"{0,1}([^\";\s]+)\"{0,1}', attrs) or re.match(r'(\w+)', attrs)
        if m:
            tid = m.group(1).strip() + '.t1'
            parts[8] = attrs.rstrip('; ') + '; transcript_id \"' + tid + '\";'
            repaired += 1
        fout.write('\t'.join(parts) + '\n')
print(f'GTF repair: {ok} OK, {repaired} transcript_id fields injected -> {out_path}')
"
[ -f "$fixed_gtf" ] || { echo "GTF repair script failed"; return 1; }
# Use full path to gffread -- conda activate inside a bash function does not
# always update PATH reliably in all shell environments.
local gffread_bin="${HOME}/miniconda3/envs/braker_env/bin/gffread"  # Already good, minor tweak
[ -x "$gffread_bin" ] || {
echo "Error: gffread not found at $gffread_bin"
echo "Fix: mamba install -n braker_env -c bioconda gffread"
return 1
}
echo "Using gffread: $gffread_bin"
"$gffread_bin" "$fixed_gtf" -g "$braker_genome" -y "$anno_prefix.faa" || return 1
return 0
}
_step3_prokka() {
conda activate braker_env || return 1
prokka --outdir anno --prefix "$anno_prefix" "$anno_genome" --cpus $num_threads
local rc=$?
conda deactivate
return $rc
}
if [ "$is_eukaryote" = "true" ]; then
run_step "step3_annotate" true "Annotate eukaryote genome with BRAKER" _step3_braker
anno_dir="$(pwd)/braker_out"
else
run_step "step3_annotate" true "Annotate prokaryote genome with Prokka" _step3_prokka
anno_dir="anno"
fi
# ─────────────────────────────────────────────
# STEP 4 — Protein extraction and filtering [CRITICAL]
# ─────────────────────────────────────────────
_step4_protein_extract() {
conda activate braker_env || return 1
local min_len=100
if [ "$is_eukaryote" = "true" ]; then min_len=200; fi
local faa_file="$anno_dir/${anno_prefix}.faa"
FAA_FILE="${faa_file}" python3 - <<EOF
import sys, re, os
from Bio import SeqIO
faa_file = os.environ['FAA_FILE']
all_records = list(SeqIO.parse(faa_file, 'fasta'))
by_gene = {}
for r in all_records:
    import re as _re
    m = _re.match(r'^(.+)\.t\d+$', r.id)
    gene_id = m.group(1) if m else r.id
    if gene_id not in by_gene or len(r.seq) > len(by_gene[gene_id].seq):
        by_gene[gene_id] = r
proteins = sorted([r for r in by_gene.values() if len(r.seq) > $min_len], key=lambda r: r.id)
SeqIO.write(proteins, 'proteins.faa', 'fasta')
dropped = sum(1 for r in by_gene.values() if len(r.seq) <= $min_len)
print(f'Input isoforms : {len(all_records)}')
print(f'Unique genes : {len(by_gene)}')
print(f'Dropped <=${min_len}aa: {dropped}')
print(f'Written : {len(proteins)} proteins -> proteins.faa')
EOF
local rc=$?
conda deactivate
return $rc
}
run_step "step4_protein_extract" true "Extract and filter proteins >100aa" _step4_protein_extract
[ -f "proteins.faa" ] || { log "HALT" "proteins.faa missing after step 4 — cannot continue"; exit 1; }
# ─────────────────────────────────────────────
# STEP 5 — RFdiffusion backbone design [non-critical]
# ─────────────────────────────────────────────
_step5_rfdiffusion() {
# RFdiffusion requires the SE3nv environment (installed from ~/RFdiffusion/env/SE3nv.yml)
conda activate SE3nv || { echo "Error: SE3nv conda env not found."; echo "Run: conda env create -f ~/RFdiffusion/env/SE3nv.yml && conda activate SE3nv && pip install ~/RFdiffusion/env/SE3Transformer && pip install -e ~/RFdiffusion --no-deps"; return 1; }
mkdir -p designs
local num_designs=$(grep -c '^>' proteins.faa)
# Use a ColabFold-predicted structure as motif if one exists (from a pre-run).
# Fall back to de novo unconditional design otherwise.
local motif_pdb=""
if [ -d "colabfold_pre" ]; then
motif_pdb=$(find colabfold_pre -name "*rank_001*.pdb" 2>/dev/null | head -1)
fi
if [ -n "$motif_pdb" ] && [ -f "$motif_pdb" ]; then
echo "Using ColabFold structure as motif: $motif_pdb"
local contig="'contigmap.contigs=[10-40/A1-50/10-40]'"
local input_pdb="inference.input_pdb=\"$motif_pdb\""
else
echo "No motif PDB found in colabfold_pre/ -- running de novo unconditional design"
echo "Tip: pre-run ColabFold on proteins.faa and place results in colabfold_pre/ to enable motif-scaffolding"
local contig="'contigmap.contigs=[100-200]'"
local input_pdb=""
fi
# Parallelize RFdiffusion designs using GNU parallel if available
if command -v parallel >/dev/null; then
echo "Using GNU parallel for RFdiffusion designs"
seq 0 $((num_designs - 1)) | parallel -j $((num_threads / 4)) 'python "${HOME}/RFdiffusion/scripts/run_inference.py" \
  inference.output_prefix=designs/'"${anno_prefix}"'_design_{} \
  inference.model_directory_path="${HOME}/models" \
  '"$input_pdb"' \
  inference.num_designs=1 \
  '"$contig"' \
  inference.ckpt_override_path="${HOME}/models/Base_ckpt.pt"'
else
local i=0
local total_elapsed=0
for ((design_num=0; design_num<num_designs; design_num++)); do
i=$((i + 1))
echo "Generating design $i/$num_designs"
local start_time=$(date +%s)
python "${HOME}/RFdiffusion/scripts/run_inference.py" \
inference.output_prefix=designs/"${anno_prefix}"_design_${design_num} \
inference.model_directory_path="${HOME}/models" \
$input_pdb \
inference.num_designs=1 \
$contig \
inference.ckpt_override_path="${HOME}/models/Base_ckpt.pt"
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_designs - i)))
progress $i $num_designs $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed design $i/$num_designs in ${elapsed}s"
if [ $i -lt $num_designs ]; then
echo "Estimated time remaining: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
local rc=$?
conda deactivate
return $rc
}
run_step "step5_rfdiffusion" false "Design backbones with RFdiffusion" _step5_rfdiffusion
# ─────────────────────────────────────────────
# STEP 5b — ProteinMPNN sequence design [non-critical]
# ─────────────────────────────────────────────
_step5b_proteinmpnn() {
conda activate SE3nv || return 1
mkdir -p proteinmpnn_out/seqs proteinmpnn_out/split_seqs
local pdb_list=$(ls designs/*.pdb 2>/dev/null)
local num_pdbs=$(echo "$pdb_list" | wc -l)
# Parallelize ProteinMPNN if GNU parallel available
if command -v parallel >/dev/null; then
echo "Using GNU parallel for ProteinMPNN"
parallel -j $((num_threads / 4)) 'python "${HOME}/ProteinMPNN/protein_mpnn_run.py" \
  --pdb_path {} \
  --out_folder "proteinmpnn_out/" \
  --num_seq_per_target 8 \
  --sampling_temp "0.1" \
  --batch_size 8' ::: $pdb_list
else
local i=0
local total_elapsed=0
for pdb in $pdb_list; do
i=$((i + 1))
echo "Processing PDB $i/$num_pdbs: $pdb"
local start_time=$(date +%s)
python "${HOME}/ProteinMPNN/protein_mpnn_run.py" \
--pdb_path "$pdb" \
--out_folder "proteinmpnn_out/" \
--num_seq_per_target 8 \
--sampling_temp "0.1" \
--batch_size 8
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_pdbs - i)))
progress $i $num_pdbs $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed PDB $i/$num_pdbs in ${elapsed}s"
if [ $i -lt $num_pdbs ]; then
echo "Estimated time remaining: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
# Split multi-sequence .fa into single-sequence .fasta files
conda activate braker_env || return 1
for fa in proteinmpnn_out/seqs/*.fa; do
  [ -f "$fa" ] || continue
design_basename=$(basename "$fa" .fa)
mkdir -p "proteinmpnn_out/split_seqs/$design_basename"
python3 -c "
from Bio import SeqIO
records = list(SeqIO.parse('$fa', 'fasta'))
for idx, rec in enumerate(records):
    output_file = 'proteinmpnn_out/split_seqs/$design_basename/${design_basename}_seq{idx}.fasta'.format(idx=idx)
    SeqIO.write([rec], output_file, 'fasta')
print(f'Split {len(records)} sequences from $fa')
" || { echo "Splitting failed for $fa"; }
done
conda deactivate
local rc=$?
conda deactivate
return $rc
}
run_step "step5b_proteinmpnn" false "Design sequences for RFdiffusion backbones with ProteinMPNN" _step5b_proteinmpnn
# ─────────────────────────────────────────────
# STEP 6 — ColabFold structure prediction [non-critical]
# ─────────────────────────────────────────────
_step6_colabfold() {
local seqs_dir="proteinmpnn_out/split_seqs"
if [ ! -d "$seqs_dir" ] || [ $(find "$seqs_dir" -name "*.fasta" 2>/dev/null | wc -l) -eq 0 ]; then
echo "No split seqs found — nothing to fold"
return 0
fi
conda activate colabfold || return 1
local fasta_list=$(find "$seqs_dir" -name "*.fasta" 2>/dev/null)
local num_fasta=$(echo "$fasta_list" | wc -l)
echo "Starting batch processing of $num_fasta FASTA files with ColabFold..."
# Parallelize ColabFold if GNU parallel available (ColabFold supports GPU, but for CPU, limit jobs)
if command -v parallel >/dev/null; then
echo "Using GNU parallel for ColabFold"
parallel -j $((num_threads / 6)) 'colabfold_batch {} colabfold_out --num-recycle 3 --use-gpu-relax' ::: $fasta_list
else
local i=0
local total_elapsed=0
for fasta in $fasta_list; do
i=$((i + 1))
echo "Processing FASTA $i/$num_fasta: $fasta"
local start_time=$(date +%s)
colabfold_batch "$fasta" colabfold_out --num-recycle 3 --use-gpu-relax # Comment out --use-gpu-relax if no GPU or issues
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_fasta - i)))
progress $i $num_fasta $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed FASTA $i/$num_fasta in ${elapsed}s"
if [ $i -lt $num_fasta ]; then
echo "Estimated time remaining: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
local rc=$?
conda deactivate
return $rc
}
run_step "step6_colabfold" false "Predict/refine structures with ColabFold" _step6_colabfold
# ─────────────────────────────────────────────
# STEP 7 — BLAST validation [non-critical]
# ─────────────────────────────────────────────
_step7_blast() {
local seqs_dir="proteinmpnn_out/split_seqs"
if [ $(find "$seqs_dir" -name "*.fasta" 2>/dev/null | wc -l) -eq 0 ]; then
echo "No design FASTA found — skipping BLAST"
return 0
fi
conda activate braker_env || return 1
mkdir -p blast_results
local fasta_list=$(find "$seqs_dir" -name "*.fasta" 2>/dev/null)
local num_fasta=$(echo "$fasta_list" | wc -l)
echo "Starting batch BLAST of $num_fasta FASTA files..."
# Parallelize BLAST loop with GNU parallel if available (blastp is single-threaded, so parallel jobs)
if command -v parallel >/dev/null; then
echo "Using GNU parallel for BLAST"
parallel -j $((num_threads / 2)) '
  basename=$(basename {} .fasta)
  blastp -query {} -db nr -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -max_target_seqs 1 -remote > "blast_results/${basename}.blast"
  ' ::: $fasta_list
else
local i=0
local total_elapsed=0
for seq_fa in $fasta_list; do
i=$((i + 1))
echo "Processing FASTA $i/$num_fasta: $seq_fa"
local start_time=$(date +%s)
basename=$(basename "$seq_fa" .fasta)
blastp -query "$seq_fa" -db nr -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -max_target_seqs 1 -remote > "blast_results/${basename}.blast" # Added -remote for remote BLAST if local db not available
local elapsed=$(( $(date +%s) - start_time ))
total_elapsed=$((total_elapsed + elapsed))
local avg=$((total_elapsed / i))
local remaining=$((avg * (num_fasta - i)))
progress $i $num_fasta $elapsed $avg $remaining
printf "\n" > /dev/tty
echo "Completed FASTA $i/$num_fasta in ${elapsed}s"
if [ $i -lt $num_fasta ]; then
echo "Estimated time remaining: ${remaining}s"
fi
done
printf "\n" > /dev/tty
fi
# Combine all .blast into blast_results.txt for reference
cat blast_results/*.blast > blast_results.txt
local rc=$?
conda deactivate
return $rc
}
run_step "step7_blast" false "Validate designs with remote BLAST" _step7_blast
# ─────────────────────────────────────────────
# STEP 8 — Rename based on putative function [non-critical]
# ─────────────────────────────────────────────
_step8_rename() {
conda activate braker_env || return 1
python3 -c "
import os
import re
anno_prefix = '$anno_prefix'
blast_dir = 'blast_results'
colab_dir = 'colabfold_out'
for blast_file in os.listdir(blast_dir):
    if not blast_file.endswith('.blast'):
        continue
    basename = blast_file.replace('.blast', '')
    function = 'unknown'
    path = os.path.join(blast_dir, blast_file)
    with open(path, 'r') as f:
        line = f.readline().strip()
        if line:
            parts = line.split('\t')
            if len(parts) >= 13:
                stitle = parts[12]
                function = re.sub(r'[^a-zA-Z0-9_]', '_', stitle)[:50]
    # Find all files in colabfold_out starting with basename_
    files = [f for f in os.listdir(colab_dir) if f.startswith(basename + '_')]
    for old_name in files:
        suffix = old_name[len(basename) + 1:]
        new_name = f'{anno_prefix}_{function}_{suffix}'
        os.rename(os.path.join(colab_dir, old_name), os.path.join(colab_dir, new_name))
print('Renaming completed based on BLAST results')
" || { echo "Renaming failed"; return 1; }
local rc=$?
conda deactivate
return $rc
}
run_step "step8_rename" false "Rename final structures based on putative functions from BLAST" _step8_rename
# ─────────────────────────────────────────────
# Pipeline summary
# ─────────────────────────────────────────────
echo "" | tee -a "$master_log"
log "INFO" "════════════════════════════════════════"
log "INFO" "Pipeline complete. Step summary:"
for sentinel in \
"$checkpoint_dir/.step1_genome_fetch.done" "$checkpoint_dir/.step1_genome_fetch.failed" \
"$checkpoint_dir/.step1b_sanitize_headers.done" "$checkpoint_dir/.step1b_sanitize_headers.failed" \
"$checkpoint_dir/.step2_qc.done" "$checkpoint_dir/.step2_qc.failed" \
"$checkpoint_dir/.step2b_repeatmask.done" "$checkpoint_dir/.step2b_repeatmask.failed" \
"$checkpoint_dir/.step3_annotate.done" "$checkpoint_dir/.step3_annotate.failed" \
"$checkpoint_dir/.step4_protein_extract.done" "$checkpoint_dir/.step4_protein_extract.failed" \
"$checkpoint_dir/.step5_rfdiffusion.done" "$checkpoint_dir/.step5_rfdiffusion.failed" \
"$checkpoint_dir/.step5b_proteinmpnn.done" "$checkpoint_dir/.step5b_proteinmpnn.failed" \
"$checkpoint_dir/.step6_colabfold.done" "$checkpoint_dir/.step6_colabfold.failed" \
"$checkpoint_dir/.step7_blast.done" "$checkpoint_dir/.step7_blast.failed" \
"$checkpoint_dir/.step8_rename.done" "$checkpoint_dir/.step8_rename.failed"; do
    [ -e "$sentinel" ] || continue
label=$(basename "$sentinel")
if [[ "$label" == *.done ]]; then
status="✓ DONE "
else
status="✗ FAILED"
fi
log "INFO" " $status $label — $(cat "$sentinel")"
done
log "INFO" "════════════════════════════════════════"
log "INFO" "Annotations : $anno_dir/"
log "INFO" "Designs : designs/"
log "INFO" "Structures : colabfold_out/"
log "INFO" "Logs : $log_dir/"
log "INFO" "Master log : $master_log"
touch "$checkpoint_dir/.pipeline_complete.done"
# Update latest symlink
ln -sfn "$(basename "$(pwd)")" "$out_dir/latest"
# ─────────────────────────────────────────────
# Generate auto-documentation report
# ─────────────────────────────────────────────
report_file="$(pwd)/README.md"
echo "# Genome-to-Design Pipeline Run Report" > "$report_file"
echo "" >> "$report_file"
echo "## Run Parameters" >> "$report_file"
echo "- **Organism**: $organism" >> "$report_file"
echo "- **Output Directory**: $out_dir" >> "$report_file"
echo "- **Eukaryote Mode**: $is_eukaryote" >> "$report_file"
echo "- **BAM File**: ${bam_file:-None}" >> "$report_file"
echo "- **Auto RNA-Seq**: $auto_rnaseq" >> "$report_file"
echo "- **Force Rerun**: $force_rerun" >> "$report_file"
echo "- **Run Start Time**: $(date '+%F %T')" >> "$report_file"
echo "" >> "$report_file"
echo "## Step Summary" >> "$report_file"
for sentinel in \
"$checkpoint_dir/.step1_genome_fetch.done" "$checkpoint_dir/.step1_genome_fetch.failed" \
"$checkpoint_dir/.step1b_sanitize_headers.done" "$checkpoint_dir/.step1b_sanitize_headers.failed" \
"$checkpoint_dir/.step2_qc.done" "$checkpoint_dir/.step2_qc.failed" \
"$checkpoint_dir/.step2b_repeatmask.done" "$checkpoint_dir/.step2b_repeatmask.failed" \
"$checkpoint_dir/.step3_annotate.done" "$checkpoint_dir/.step3_annotate.failed" \
"$checkpoint_dir/.step4_protein_extract.done" "$checkpoint_dir/.step4_protein_extract.failed" \
"$checkpoint_dir/.step5_rfdiffusion.done" "$checkpoint_dir/.step5_rfdiffusion.failed" \
"$checkpoint_dir/.step5b_proteinmpnn.done" "$checkpoint_dir/.step5b_proteinmpnn.failed" \
"$checkpoint_dir/.step6_colabfold.done" "$checkpoint_dir/.step6_colabfold.failed" \
"$checkpoint_dir/.step7_blast.done" "$checkpoint_dir/.step7_blast.failed" \
"$checkpoint_dir/.step8_rename.done" "$checkpoint_dir/.step8_rename.failed"; do
    [ -e "$sentinel" ] || continue
label=$(basename "$sentinel" .done | basename - .failed)
if [[ "$sentinel" == *.done ]]; then
status="✓ Done"
else
status="✗ Failed"
fi
echo "- **$label**: $status - $(cat "$sentinel")" >> "$report_file"
done
echo "" >> "$report_file"
echo "## Key Outputs" >> "$report_file"
echo "- Annotations: $anno_dir/" >> "$report_file"
echo "- Designs: designs/" >> "$report_file"
echo "- Structures: colabfold_out/" >> "$report_file"
echo "- Logs: $log_dir/" >> "$report_file"
echo "- Master Log: $master_log" >> "$report_file"
if [ -f "blast_results.txt" ]; then
echo "- BLAST Results: blast_results.txt" >> "$report_file"
fi
echo "" >> "$report_file"
echo "## Citations" >> "$report_file"
echo "This pipeline uses the following tools (please cite accordingly):" >> "$report_file"
echo "- BRAKER (for eukaryotic annotation)" >> "$report_file"
echo "- Prokka (for prokaryotic annotation)" >> "$report_file"
echo "- RFdiffusion (for backbone design)" >> "$report_file"
echo "- ProteinMPNN (for sequence design)" >> "$report_file"
echo "- ColabFold (for structure prediction)" >> "$report_file"
echo "- BLAST (for validation)" >> "$report_file"
echo "- Biopython, Entrez, etc. (for various utilities)" >> "$report_file"
log "INFO" "Generated report: $report_file"
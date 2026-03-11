#!/bin/bash
# =============================================================================
# feedback5_designed_proteins_to_annotation.sh
#
# FEEDBACK LOOP 5: Validated designed proteins → enriched annotation hints
#
# WORKFLOW:
#   1. Filter ColabFold-predicted designed proteins by pLDDT threshold to
#      keep only high-confidence designs (they are effectively validated
#      structural predictions at this point)
#   2. Back-translate the amino acid sequences to codon-optimized DNA
#      (optional, only needed if running a new genome_to_design.sh later)
#   3. Combine the passing designed protein sequences with the original
#      OrthoDB hint file
#   4. Write the combined hint set as designed_hints.faa for use with either:
#       a) The current organism: re-run BRAKER to potentially recover genes
#          in exactly the protein families the designed sequences represent
#       b) A RELATED organism: pass as --prot_seq when running a fresh
#          genome_to_design.sh on a new species to bootstrap annotation
#          quality using our own curated, high-confidence sequences
#   5. Write a manifest listing which designed sequences were included,
#      their mean pLDDT, source backbone, and the BLAST function label
#
# USAGE:
#   ./feedback5_designed_proteins_to_annotation.sh [config.yaml|config.json|config.txt]
#   ./feedback5_designed_proteins_to_annotation.sh \
#       --run_dir=./output/latest \
#       --plddt_min=75 \
#       --target_organism="Nicotiana benthamiana" \
#       --target_output=./nicotiana_output
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir            Source pipeline run directory (default: ./output/latest)
#   plddt_min          Minimum ColabFold mean pLDDT to include a design (default: 75)
#   target_organism    If set, automatically launch genome_to_design.sh on a
#                      new organism using the enriched hints (optional)
#   target_output      Output directory for the target organism run (optional)
#   target_eukaryote   true/false for the target organism (default: true)
#   include_blast_func true  = annotate sequences with BLAST function labels
#                      false = use generic design IDs (default: true)
#   output_dir         Where to write loop outputs (default: <run_dir>/feedback5_loop)
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# Shared helpers
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
# Config parsing
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
PLDDT_MIN=75
TARGET_ORGANISM=""
TARGET_OUTPUT=""
TARGET_EUKARYOTE="true"
INCLUDE_BLAST_FUNC="true"
OUTPUT_DIR=""

if [ $# -ge 1 ] && [ -f "$1" ]; then
    eval "$(parse_config "$1")"
    [ -n "${CFG_RUN_DIR:-}"           ] && RUN_DIR="$CFG_RUN_DIR"
    [ -n "${CFG_PLDDT_MIN:-}"         ] && PLDDT_MIN="$CFG_PLDDT_MIN"
    [ -n "${CFG_TARGET_ORGANISM:-}"   ] && TARGET_ORGANISM="$CFG_TARGET_ORGANISM"
    [ -n "${CFG_TARGET_OUTPUT:-}"     ] && TARGET_OUTPUT="$CFG_TARGET_OUTPUT"
    [ -n "${CFG_TARGET_EUKARYOTE:-}"  ] && TARGET_EUKARYOTE="$CFG_TARGET_EUKARYOTE"
    [ -n "${CFG_INCLUDE_BLAST_FUNC:-}"] && INCLUDE_BLAST_FUNC="$CFG_INCLUDE_BLAST_FUNC"
    [ -n "${CFG_OUTPUT_DIR:-}"        ] && OUTPUT_DIR="$CFG_OUTPUT_DIR"
    shift
fi

for arg in "$@"; do
    case "$arg" in
        --run_dir=*)           RUN_DIR="${arg#*=}" ;;
        --plddt_min=*)         PLDDT_MIN="${arg#*=}" ;;
        --target_organism=*)   TARGET_ORGANISM="${arg#*=}" ;;
        --target_output=*)     TARGET_OUTPUT="${arg#*=}" ;;
        --target_eukaryote=*)  TARGET_EUKARYOTE="${arg#*=}" ;;
        --include_blast_func=*) INCLUDE_BLAST_FUNC="${arg#*=}" ;;
        --output_dir=*)        OUTPUT_DIR="${arg#*=}" ;;
    esac
done

if [ -z "$RUN_DIR" ]; then
    RUN_DIR=$(realpath ./output/latest 2>/dev/null || echo "")
fi
[ -d "$RUN_DIR" ] || { echo "ERROR: run_dir '$RUN_DIR' not found."; exit 1; }
RUN_DIR=$(realpath "$RUN_DIR")
OUTPUT_DIR="${OUTPUT_DIR:-${RUN_DIR}/feedback5_loop}"
mkdir -p "$OUTPUT_DIR"
checkpoint_dir="$OUTPUT_DIR"
log_dir="$OUTPUT_DIR/logs"
mkdir -p "$log_dir"
master_log="$log_dir/feedback5.log"

log INFO "=========================================="
log INFO "Feedback Loop 5: Designed proteins → annotation hints"
log INFO "run_dir          : $RUN_DIR"
log INFO "plddt_min        : $PLDDT_MIN"
log INFO "target_organism  : ${TARGET_ORGANISM:-'(none — generate hints only)'}"
log INFO "include_blast    : $INCLUDE_BLAST_FUNC"
log INFO "output_dir       : $OUTPUT_DIR"
log INFO "=========================================="

. ~/miniconda3/etc/profile.d/conda.sh

# ─────────────────────────────────────────────
# STEP FB5-1 — Filter validated designs by pLDDT [CRITICAL]
# Reads ColabFold JSON scores, cross-references with ProteinMPNN split_seqs
# to get the amino acid sequences, and selects those above the threshold.
# If include_blast_func=true, also looks up the BLAST function label
# (from blast_results/) and uses it as the sequence description.
# ─────────────────────────────────────────────
VALIDATED_FAA="$OUTPUT_DIR/validated_designs.faa"
MANIFEST="$OUTPUT_DIR/feedback5_manifest.tsv"

_fb5_filter_validated() {
    conda activate braker_env || return 1
    python3 - <<PYEOF
import json, glob, os, csv
from Bio import SeqIO

colab_dir     = "$RUN_DIR/colabfold_out"
split_dir     = "$RUN_DIR/proteinmpnn_out/split_seqs"
blast_dir     = "$RUN_DIR/blast_results"
out_faa       = "$VALIDATED_FAA"
manifest_path = "$MANIFEST"
plddt_min     = float("$PLDDT_MIN")
use_blast     = "$INCLUDE_BLAST_FUNC".lower() == "true"

# ── 1. Build pLDDT map: design_stem → mean pLDDT ──
plddt_map = {}
for sj in glob.glob(f"{colab_dir}/**/*scores*rank_001*.json", recursive=True):
    try:
        data   = json.load(open(sj))
        plddt  = data.get("plddt", [])
        if not plddt:
            continue
        mean = sum(plddt) / len(plddt)
        stem = os.path.basename(sj).split("_scores")[0]
        plddt_map[stem] = max(plddt_map.get(stem, 0), mean)
    except Exception as e:
        print(f"WARN: {sj}: {e}")

print(f"Scored designs: {len(plddt_map)}")

# ── 2. Build BLAST function map (optional) ──
blast_map = {}
if use_blast and os.path.isdir(blast_dir):
    for bf in glob.glob(f"{blast_dir}/*.blast"):
        stem = os.path.basename(bf).replace(".blast", "")
        with open(bf) as f:
            line = f.readline().strip()
            if line:
                parts = line.split("\t")
                if len(parts) >= 13:
                    import re
                    fn = re.sub(r'[^a-zA-Z0-9_\s]', '', parts[12])[:60].strip()
                    blast_map[stem] = fn

# ── 3. Load FASTA sequences for passing designs ──
records  = []
manifest = []
for fasta_file in glob.glob(f"{split_dir}/**/*.fasta", recursive=True):
    stem = os.path.basename(fasta_file).replace(".fasta", "")
    # Match stem to pLDDT map (may have _seqN suffix)
    score = None
    for key in plddt_map:
        if stem.startswith(key) or key.startswith(stem):
            score = plddt_map[key]
            break
    if score is None:
        # Try partial match
        for key in plddt_map:
            if key in stem or stem in key:
                score = plddt_map[key]
                break
    if score is None or score < plddt_min:
        continue

    try:
        recs = list(SeqIO.parse(fasta_file, "fasta"))
        if not recs:
            continue
        rec = recs[0]
        # Translate DNA → AA if the source is a nucleotide FASTA
        seq = str(rec.seq)
        if set(seq.upper()) - set("ACGTN-"):
            # Likely already amino acid
            aa_seq = seq
        else:
            from Bio.Seq import Seq
            aa_seq = str(Seq(seq).translate(to_stop=True))

        func_label = blast_map.get(stem, "designed_protein")
        rec.id          = f"{stem}|pLDDT={score:.1f}|{func_label}"
        rec.description = ""
        rec.seq         = type(rec.seq)(aa_seq)
        records.append(rec)
        manifest.append({
            "design":      stem,
            "mean_plddt":  f"{score:.1f}",
            "blast_func":  func_label,
            "fasta_src":   fasta_file,
        })
    except Exception as e:
        print(f"WARN: could not process {fasta_file}: {e}")

SeqIO.write(records, out_faa, "fasta")
print(f"Sequences passing pLDDT >= {plddt_min}: {len(records)}")
print(f"Written to: {out_faa}")

with open(manifest_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["design", "mean_plddt", "blast_func", "fasta_src"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(manifest)
print(f"Manifest: {manifest_path}")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb5_filter_validated" true \
    "Filter designed proteins by pLDDT >= $PLDDT_MIN" \
    _fb5_filter_validated

[ -s "$VALIDATED_FAA" ] || {
    log HALT "No sequences passed the pLDDT filter — lower plddt_min or run more designs"
    exit 1
}

# ─────────────────────────────────────────────
# STEP FB5-2 — Merge with OrthoDB hints [non-critical]
# Produces an enriched hint file for use as --prot_seq in BRAKER.
# De-duplicates by sequence ID.
# ─────────────────────────────────────────────
ENRICHED_HINTS="$OUTPUT_DIR/enriched_with_designs.faa"

_fb5_merge_hints() {
    conda activate braker_env || return 1
    local orthodb_file
    orthodb_file=$(ls "${HOME}"/orthodb_*.fasta 2>/dev/null | head -1)

    if [ -z "$orthodb_file" ]; then
        log WARN "No orthodb_*.fasta found — enriched hints will be designs only"
        cp "$VALIDATED_FAA" "$ENRICHED_HINTS"
        conda deactivate
        return 0
    fi

    log INFO "Merging with $orthodb_file"
    python3 - <<PYEOF
from Bio import SeqIO

seen = set()
merged = []
for fn in ["$orthodb_file", "$VALIDATED_FAA"]:
    for rec in SeqIO.parse(fn, "fasta"):
        if rec.id not in seen:
            seen.add(rec.id)
            merged.append(rec)

SeqIO.write(merged, "$ENRICHED_HINTS", "fasta")
print(f"Merged {len(merged)} unique sequences -> $ENRICHED_HINTS")
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb5_merge_hints" false \
    "Merge validated designs with OrthoDB hints" \
    _fb5_merge_hints

# ─────────────────────────────────────────────
# STEP FB5-3 — Optionally launch genome_to_design.sh on a new organism
# If target_organism is set, this step passes the enriched hints as an
# additional protein hint file by copying it to the expected OrthoDB
# location for the target organism's taxonomic group.
# The user is responsible for having genome_to_design.sh in PATH or $HOME.
# ─────────────────────────────────────────────
_fb5_launch_target_run() {
    if [ -z "$TARGET_ORGANISM" ]; then
        echo "No target_organism set — skipping target run"
        return 0
    fi
    [ -n "$TARGET_OUTPUT" ] || TARGET_OUTPUT="./output_${TARGET_ORGANISM// /_}"
    mkdir -p "$TARGET_OUTPUT"

    # Place the enriched hints where BRAKER will find them.
    # genome_to_design.sh downloads to ~/orthodb_<group>.fasta;
    # we cannot know the group in advance, so we write a "custom" file
    # and pass it via the BRAKER command-line through a wrapper env var.
    local hints_dest="${HOME}/orthodb_designed_hints.fasta"
    cp "$ENRICHED_HINTS" "$hints_dest"
    log INFO "Enriched hints copied to $hints_dest"
    log INFO "Launching genome_to_design.sh for '$TARGET_ORGANISM'"
    log INFO "  Target output: $TARGET_OUTPUT"
    log INFO "  Eukaryote    : $TARGET_EUKARYOTE"

    # Find genome_to_design.sh — check PATH then common locations
    local gd_script
    gd_script=$(command -v genome_to_design.sh 2>/dev/null || echo "${HOME}/genome_to_design.sh")
    [ -f "$gd_script" ] || {
        echo "ERROR: genome_to_design.sh not found. Place it in PATH or \$HOME."
        return 1
    }

    bash "$gd_script" \
        "$TARGET_ORGANISM" \
        "$TARGET_OUTPUT" \
        "$TARGET_EUKARYOTE" \
        --auto_rnaseq \
        2>&1 | tee -a "$log_dir/fb5_target_run.log"
}

run_step "fb5_launch_target" false \
    "Optionally launch genome_to_design.sh for target organism with enriched hints" \
    _fb5_launch_target_run

# ─────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────
n_validated=$(grep -c '^>' "$VALIDATED_FAA" 2>/dev/null || echo 0)

log INFO "=========================================="
log INFO "Feedback Loop 5 complete."
log INFO "Validated designs (pLDDT >= $PLDDT_MIN) : $n_validated"
log INFO "Validated FAA   : $VALIDATED_FAA"
log INFO "Enriched hints  : $ENRICHED_HINTS"
log INFO "Manifest        : $MANIFEST"
log INFO ""
log INFO "To use enriched hints for current organism re-annotation:"
log INFO "  See feedback3_blast_to_braker.sh with the --enriched_hints flag"
log INFO ""
log INFO "To use as --prot_seq for a new genome_to_design.sh run:"
log INFO "  export ORTHODB_OVERRIDE=$ENRICHED_HINTS"
log INFO "  ./genome_to_design.sh 'New Organism' ./new_output true"
log INFO "=========================================="
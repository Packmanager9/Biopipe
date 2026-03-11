#!/bin/bash
# =============================================================================
# feedback1_colabfold_to_rfdiffusion.sh
#
# FEEDBACK LOOP 1: ColabFold outputs → RFdiffusion motif-scaffolded redesign
#
# WORKFLOW:
#   1. Parse ColabFold output JSONs to rank all predicted structures by mean
#      pLDDT score across all residues
#   2. Select the top-N highest-confidence structures (default N=3)
#   3. Feed each selected PDB back to RFdiffusion as a motif, using
#      contigmap.contigs to scaffold variable-length flanking regions around
#      the fixed motif residues
#   4. Run ProteinMPNN on the new backbones to generate sequences
#   5. Run ColabFold again on those sequences (one more cycle)
#   6. Write a summary comparing original vs redesigned pLDDT scores
#
# This is the "exploit the best hits" loop: instead of treating every
# RFdiffusion design equally, we pick the structures that ColabFold already
# validated well and ask RFdiffusion to build new scaffolds around them.
# Each iteration should converge toward designs with higher structural
# confidence.  Run this script after genome_to_design.sh completes Step 6.
#
# USAGE:
#   ./feedback1_colabfold_to_rfdiffusion.sh [config.yaml|config.json|config.txt]
#   ./feedback1_colabfold_to_rfdiffusion.sh \
#       --run_dir=./output/latest \
#       --top_n=3 \
#       --plddt_min=70 \
#       --motif_residues="A1-50" \
#       --flank_min=10 --flank_max=40 \
#       --iterations=2 \
#       --num_designs=4
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir          Path to the pipeline run directory (default: ./output/latest)
#   top_n            Number of top structures to use as motifs (default: 3)
#   plddt_min        Minimum mean pLDDT to qualify a structure (default: 70)
#   motif_residues   Residue range for the fixed motif e.g. "A1-50" (default: A1-50)
#   flank_min        Minimum length of scaffolded flanking regions (default: 10)
#   flank_max        Maximum length of scaffolded flanking regions (default: 40)
#   iterations       Number of ColabFold→RFdiffusion cycles to run (default: 2)
#   num_designs      RFdiffusion designs per motif per iteration (default: 4)
#   num_seqs         ProteinMPNN sequences per backbone (default: 8)
#   output_dir       Where to write loop outputs (default: <run_dir>/feedback1_loop)
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# Shared helpers (mirror genome_to_design.sh)
# ─────────────────────────────────────────────
log() { echo "[$(date '+%F %T')] [$1] ${*:2}" | tee -a "$master_log"; }

progress() {
    local cur=$1 tot=$2 elapsed=$3 avg=$4 rem=$5
    printf "\r[%-50s] %d%% (%d/%d) | Last: %ds | Avg: %ds | ETA: %ds" \
        "$(printf '=%.0s' $(seq 1 $((cur * 50 / tot))))" \
        $((cur * 100 / tot)) "$cur" "$tot" "$elapsed" "$avg" "$rem" > /dev/tty
}

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
    local file=$1
    local ext="${file##*.}"
    case "$ext" in
        yaml|yml)
            # Requires python3 + pyyaml (already a pipeline dependency)
            python3 - "$file" <<'EOF'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
for k, v in cfg.items():
    print(f"CFG_{k.upper()}={v}")
EOF
            ;;
        json)
            python3 - "$file" <<'EOF'
import sys, json
cfg = json.load(open(sys.argv[1]))
for k, v in cfg.items():
    print(f"CFG_{k.upper()}={v}")
EOF
            ;;
        *)
            # Plain key: value text file
            grep -v '^\s*#' "$file" | grep ':' | while IFS=: read -r k v; do
                k=$(echo "$k" | tr '[:lower:]' '[:upper:]' | tr -d ' ')
                v=$(echo "$v" | sed 's/^ //')
                echo "CFG_${k}=${v}"
            done
            ;;
    esac
}

# ── Defaults ──
RUN_DIR=""
TOP_N=3
PLDDT_MIN=70
MOTIF_RESIDUES="A1-50"
FLANK_MIN=10
FLANK_MAX=40
ITERATIONS=2
NUM_DESIGNS=4
NUM_SEQS=8
OUTPUT_DIR=""

# ── Load config file if first arg looks like a file ──
if [ $# -ge 1 ] && [ -f "$1" ]; then
    eval "$(parse_config "$1")"
    [ -n "${CFG_RUN_DIR:-}"        ] && RUN_DIR="$CFG_RUN_DIR"
    [ -n "${CFG_TOP_N:-}"          ] && TOP_N="$CFG_TOP_N"
    [ -n "${CFG_PLDDT_MIN:-}"      ] && PLDDT_MIN="$CFG_PLDDT_MIN"
    [ -n "${CFG_MOTIF_RESIDUES:-}" ] && MOTIF_RESIDUES="$CFG_MOTIF_RESIDUES"
    [ -n "${CFG_FLANK_MIN:-}"      ] && FLANK_MIN="$CFG_FLANK_MIN"
    [ -n "${CFG_FLANK_MAX:-}"      ] && FLANK_MAX="$CFG_FLANK_MAX"
    [ -n "${CFG_ITERATIONS:-}"     ] && ITERATIONS="$CFG_ITERATIONS"
    [ -n "${CFG_NUM_DESIGNS:-}"    ] && NUM_DESIGNS="$CFG_NUM_DESIGNS"
    [ -n "${CFG_NUM_SEQS:-}"       ] && NUM_SEQS="$CFG_NUM_SEQS"
    [ -n "${CFG_OUTPUT_DIR:-}"     ] && OUTPUT_DIR="$CFG_OUTPUT_DIR"
    shift
fi

# ── CLI flags override config file ──
for arg in "$@"; do
    case "$arg" in
        --run_dir=*)        RUN_DIR="${arg#*=}" ;;
        --top_n=*)          TOP_N="${arg#*=}" ;;
        --plddt_min=*)      PLDDT_MIN="${arg#*=}" ;;
        --motif_residues=*) MOTIF_RESIDUES="${arg#*=}" ;;
        --flank_min=*)      FLANK_MIN="${arg#*=}" ;;
        --flank_max=*)      FLANK_MAX="${arg#*=}" ;;
        --iterations=*)     ITERATIONS="${arg#*=}" ;;
        --num_designs=*)    NUM_DESIGNS="${arg#*=}" ;;
        --num_seqs=*)       NUM_SEQS="${arg#*=}" ;;
        --output_dir=*)     OUTPUT_DIR="${arg#*=}" ;;
    esac
done

# ── Resolve run_dir ──
if [ -z "$RUN_DIR" ]; then
    RUN_DIR=$(realpath ./output/latest 2>/dev/null || echo "")
fi
[ -d "$RUN_DIR" ] || { echo "ERROR: run_dir '$RUN_DIR' not found."; exit 1; }
RUN_DIR=$(realpath "$RUN_DIR")

OUTPUT_DIR="${OUTPUT_DIR:-${RUN_DIR}/feedback1_loop}"
mkdir -p "$OUTPUT_DIR"
checkpoint_dir="$OUTPUT_DIR"
log_dir="$OUTPUT_DIR/logs"
mkdir -p "$log_dir"
master_log="$log_dir/feedback1.log"
num_threads=$(nproc)

log INFO "=========================================="
log INFO "Feedback Loop 1: ColabFold → RFdiffusion"
log INFO "run_dir:        $RUN_DIR"
log INFO "top_n:          $TOP_N"
log INFO "plddt_min:      $PLDDT_MIN"
log INFO "motif_residues: $MOTIF_RESIDUES"
log INFO "flanks:         $FLANK_MIN–$FLANK_MAX"
log INFO "iterations:     $ITERATIONS"
log INFO "num_designs:    $NUM_DESIGNS"
log INFO "output_dir:     $OUTPUT_DIR"
log INFO "=========================================="

. ~/miniconda3/etc/profile.d/conda.sh

# ─────────────────────────────────────────────
# STEP FB1-1 — Rank ColabFold structures by mean pLDDT [CRITICAL]
# Parse all *scores*.json files produced by ColabFold; each contains
# per-residue pLDDT arrays for up to 5 ranked models.  We take the
# mean pLDDT of model rank_001 as the primary quality signal.
# ─────────────────────────────────────────────
COLAB_DIR="$RUN_DIR/colabfold_out"
TOP_PDBS_FILE="$OUTPUT_DIR/top_pdbs.txt"

_fb1_rank_structures() {
    conda activate braker_env || return 1
    python3 - <<PYEOF
import json, glob, os

colab_dir  = "$COLAB_DIR"
out_file   = "$TOP_PDBS_FILE"
top_n      = $TOP_N
plddt_min  = $PLDDT_MIN

scores_jsons = glob.glob(f"{colab_dir}/**/*scores*rank_001*.json", recursive=True)
if not scores_jsons:
    # Fallback: any scores json
    scores_jsons = glob.glob(f"{colab_dir}/**/*scores*.json", recursive=True)

results = []
for sj in scores_jsons:
    try:
        data  = json.load(open(sj))
        plddt = data.get("plddt", [])
        if not plddt:
            continue
        mean_plddt = sum(plddt) / len(plddt)
        # Derive PDB path — ColabFold names them identically except .pdb
        pdb = sj.replace("scores", "unrelaxed").replace(".json", ".pdb")
        if not os.path.isfile(pdb):
            pdb = sj.replace("_scores_", "_relaxed_").replace(".json", ".pdb")
        if not os.path.isfile(pdb):
            # Try relaxed variant
            pdb = sj.replace(".json", ".pdb").replace("scores", "relaxed")
        if not os.path.isfile(pdb):
            print(f"WARN: no PDB found for {sj}")
            continue
        results.append((mean_plddt, pdb))
    except Exception as e:
        print(f"WARN: could not parse {sj}: {e}")

# Sort descending by pLDDT, filter by threshold
results = sorted(results, reverse=True)
qualified = [(s, p) for s, p in results if s >= plddt_min]
selected  = qualified[:top_n]

print(f"Total structures found : {len(results)}")
print(f"Above pLDDT {plddt_min}        : {len(qualified)}")
print(f"Selected (top {top_n})        : {len(selected)}")

with open(out_file, "w") as f:
    for score, pdb in selected:
        f.write(f"{score:.2f}\t{pdb}\n")
        print(f"  {score:.2f}  {pdb}")

if not selected:
    print("ERROR: no structures passed the pLDDT threshold — lower plddt_min or run more designs")
    exit(1)
PYEOF
    local rc=$?
    conda deactivate
    return $rc
}

run_step "fb1_rank_structures" true \
    "Rank ColabFold structures by mean pLDDT, select top $TOP_N above $PLDDT_MIN" \
    _fb1_rank_structures

[ -f "$TOP_PDBS_FILE" ] || { log HALT "top_pdbs.txt missing — cannot continue"; exit 1; }

# ─────────────────────────────────────────────
# MAIN ITERATION LOOP
# Each iteration:
#   a) Take the current set of top PDBs as motifs
#   b) Run RFdiffusion with motif scaffolding
#   c) Run ProteinMPNN on new backbones
#   d) Run ColabFold on new sequences
#   e) Re-rank; keep only structures that beat the previous best pLDDT
#   f) Write per-iteration summary to $OUTPUT_DIR/iter_N_summary.tsv
# ─────────────────────────────────────────────
current_pdbs_file="$TOP_PDBS_FILE"

for iter in $(seq 1 "$ITERATIONS"); do
    log INFO "──────────────────────────────────────────"
    log INFO "Starting iteration $iter / $ITERATIONS"

    iter_dir="$OUTPUT_DIR/iter_${iter}"
    mkdir -p "$iter_dir/designs" "$iter_dir/mpnn_out/seqs" \
             "$iter_dir/mpnn_out/split_seqs" "$iter_dir/colabfold_out"

    iter_checkpoint="$iter_dir"
    iter_log_dir="$iter_dir/logs"
    mkdir -p "$iter_log_dir"

    # ── (a) RFdiffusion motif-scaffolded redesign ───────────────────────────
    # Read PDB list written by the ranking step (or previous iteration)
    mapfile -t pdb_entries < <(awk '{print $2}' "$current_pdbs_file")
    num_pdbs=${#pdb_entries[@]}
    log INFO "Iter $iter: running RFdiffusion on $num_pdbs motif PDBs"

    _fb1_rfdiffusion_iter() {
        conda activate SE3nv || return 1
        local i=0 total_elapsed=0
        for pdb_path in "${pdb_entries[@]}"; do
            i=$((i + 1))
            local t0; t0=$(date +%s)
            local stem; stem=$(basename "$pdb_path" .pdb)
            echo "RFdiffusion motif design $i/$num_pdbs: $stem"

            # contigmap: fixed motif flanked by variable scaffolding regions
            # e.g. "10-40/A1-50/10-40" means 10-40 aa flank, residues A1-50
            # fixed from the motif PDB, then 10-40 aa flank on the C-terminus
            local contig="${FLANK_MIN}-${FLANK_MAX}/${MOTIF_RESIDUES}/${FLANK_MIN}-${FLANK_MAX}"

            for d in $(seq 1 "$NUM_DESIGNS"); do
                python "${HOME}/RFdiffusion/scripts/run_inference.py" \
                    inference.output_prefix="$iter_dir/designs/${stem}_iter${iter}_d${d}" \
                    inference.model_directory_path="${HOME}/models" \
                    inference.input_pdb="$pdb_path" \
                    "contigmap.contigs=[$contig]" \
                    inference.num_designs=1 \
                    inference.ckpt_override_path="${HOME}/models/Base_ckpt.pt"
            done

            local elapsed=$(( $(date +%s) - t0 ))
            total_elapsed=$((total_elapsed + elapsed))
            local avg=$((total_elapsed / i))
            local rem=$((avg * (num_pdbs - i)))
            progress "$i" "$num_pdbs" "$elapsed" "$avg" "$rem"
            printf "\n" > /dev/tty
        done
        conda deactivate
    }

    run_step "fb1_iter${iter}_rfdiffusion" true \
        "Iter $iter: motif-scaffolded RFdiffusion" _fb1_rfdiffusion_iter

    # ── (b) ProteinMPNN on new backbones ────────────────────────────────────
    _fb1_mpnn_iter() {
        conda activate SE3nv || return 1
        local pdb_list; pdb_list=$(find "$iter_dir/designs" -name "*.pdb" 2>/dev/null)
        [ -n "$pdb_list" ] || { echo "No PDBs found in $iter_dir/designs"; return 1; }
        local num_pdbs_mpnn; num_pdbs_mpnn=$(echo "$pdb_list" | wc -l)
        local i=0 total_elapsed=0

        for pdb in $pdb_list; do
            i=$((i + 1))
            local t0; t0=$(date +%s)
            python "${HOME}/ProteinMPNN/protein_mpnn_run.py" \
                --pdb_path "$pdb" \
                --out_folder "$iter_dir/mpnn_out/" \
                --num_seq_per_target "$NUM_SEQS" \
                --sampling_temp "0.1" \
                --batch_size 8
            local elapsed=$(( $(date +%s) - t0 ))
            total_elapsed=$((total_elapsed + elapsed))
            local avg=$((total_elapsed / i))
            local rem=$((avg * (num_pdbs_mpnn - i)))
            progress "$i" "$num_pdbs_mpnn" "$elapsed" "$avg" "$rem"
            printf "\n" > /dev/tty
        done

        # Split multi-sequence .fa into per-sequence .fasta for ColabFold
        conda activate braker_env
        for fa in "$iter_dir/mpnn_out/seqs/"*.fa; do
            [ -f "$fa" ] || continue
            local base; base=$(basename "$fa" .fa)
            mkdir -p "$iter_dir/mpnn_out/split_seqs/$base"
            python3 -c "
from Bio import SeqIO
records = list(SeqIO.parse('$fa', 'fasta'))
for idx, rec in enumerate(records):
    out = '$iter_dir/mpnn_out/split_seqs/$base/${base}_seq{}.fasta'.format(idx)
    SeqIO.write([rec], out, 'fasta')
print(f'Split {len(records)} seqs from $fa')
"
        done
        conda deactivate
    }

    run_step "fb1_iter${iter}_mpnn" true \
        "Iter $iter: ProteinMPNN sequence design" _fb1_mpnn_iter

    # ── (c) ColabFold on new sequences ──────────────────────────────────────
    _fb1_colabfold_iter() {
        conda activate colabfold || return 1
        local fasta_list; fasta_list=$(find "$iter_dir/mpnn_out/split_seqs" -name "*.fasta" 2>/dev/null)
        [ -n "$fasta_list" ] || { echo "No FASTA files found"; return 1; }
        local num_f; num_f=$(echo "$fasta_list" | wc -l)
        local i=0 total_elapsed=0

        for fasta in $fasta_list; do
            i=$((i + 1))
            local t0; t0=$(date +%s)
            colabfold_batch "$fasta" "$iter_dir/colabfold_out" \
                --num-recycle 3 --use-gpu-relax
            local elapsed=$(( $(date +%s) - t0 ))
            total_elapsed=$((total_elapsed + elapsed))
            local avg=$((total_elapsed / i))
            local rem=$((avg * (num_f - i)))
            progress "$i" "$num_f" "$elapsed" "$avg" "$rem"
            printf "\n" > /dev/tty
        done
        conda deactivate
    }

    run_step "fb1_iter${iter}_colabfold" false \
        "Iter $iter: ColabFold structure prediction" _fb1_colabfold_iter

    # ── (d) Re-rank; write summary; update current_pdbs_file ────────────────
    iter_pdbs_file="$iter_dir/top_pdbs.txt"

    _fb1_rerank_iter() {
        conda activate braker_env || return 1
        python3 - <<PYEOF
import json, glob, os

colab_dir  = "$iter_dir/colabfold_out"
orig_file  = "$current_pdbs_file"
out_file   = "$iter_pdbs_file"
summary_f  = "$iter_dir/iter_${iter}_summary.tsv"
top_n      = $TOP_N
plddt_min  = $PLDDT_MIN

# Load previous best scores for comparison
prev_scores = {}
with open(orig_file) as f:
    for line in f:
        score, pdb = line.strip().split("\t")
        prev_scores[os.path.basename(pdb)] = float(score)

scores_jsons = glob.glob(f"{colab_dir}/**/*scores*rank_001*.json", recursive=True)
results = []
for sj in scores_jsons:
    try:
        data       = json.load(open(sj))
        plddt_vals = data.get("plddt", [])
        if not plddt_vals:
            continue
        mean_plddt = sum(plddt_vals) / len(plddt_vals)
        pdb = sj.replace("scores", "unrelaxed").replace(".json", ".pdb")
        if not os.path.isfile(pdb):
            pdb = sj.replace(".json", ".pdb").replace("scores", "relaxed")
        if not os.path.isfile(pdb):
            continue
        results.append((mean_plddt, pdb))
    except Exception as e:
        print(f"WARN: {e}")

results = sorted(results, reverse=True)
selected = [r for r in results if r[0] >= plddt_min][:top_n]

with open(out_file, "w") as f:
    for score, pdb in selected:
        f.write(f"{score:.2f}\t{pdb}\n")

# Write iteration summary TSV
with open(summary_f, "w") as f:
    f.write("pdb\tplddt_iter${iter}\tplddt_prev\tdelta\n")
    for score, pdb in results[:top_n * 2]:
        base = os.path.basename(pdb)
        prev = prev_scores.get(base, float("nan"))
        delta = score - prev if prev == prev else float("nan")
        f.write(f"{base}\t{score:.2f}\t{prev}\t{delta:+.2f}\n")

print(f"Iter $iter summary -> {summary_f}")
print(f"New top {top_n} PDBs -> {out_file}")
for s, p in selected:
    print(f"  {s:.2f}  {p}")
PYEOF
        local rc=$?
        conda deactivate
        return $rc
    }

    run_step "fb1_iter${iter}_rerank" false \
        "Iter $iter: re-rank new structures, update top PDB list" _fb1_rerank_iter

    # If re-ranking produced a new file, use it next iteration
    [ -f "$iter_pdbs_file" ] && current_pdbs_file="$iter_pdbs_file"

    log INFO "Iteration $iter complete. Best PDBs: $current_pdbs_file"
done

# ─────────────────────────────────────────────
# Final report
# ─────────────────────────────────────────────
log INFO "=========================================="
log INFO "Feedback Loop 1 complete."
log INFO "Final top structures: $current_pdbs_file"
log INFO "Per-iteration summaries:"
for s in "$OUTPUT_DIR"/iter_*/iter_*_summary.tsv; do
    [ -f "$s" ] && log INFO "  $s"
done
log INFO "=========================================="
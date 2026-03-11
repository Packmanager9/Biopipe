#!/usr/bin/env python3
# =============================================================================
# feedback2_plddt_mpnn_resample.py
#
# FEEDBACK LOOP 2: ColabFold pLDDT scores → ProteinMPNN resampling
#
# WORKFLOW:
#   1. Parse all ColabFold output JSONs and extract mean pLDDT per sequence
#   2. Split sequences into two bins:
#        PASS  — mean pLDDT >= plddt_pass  (exit loop, keep these)
#        RETRY — mean pLDDT < plddt_pass   (send backbone back to ProteinMPNN)
#   3. For RETRY backbones, look up the corresponding RFdiffusion PDB and
#      re-run ProteinMPNN with a higher temperature or more sequences to
#      diversify the search
#   4. Re-run ColabFold on the new sequences
#   5. Repeat up to max_iterations; sequences that never pass are flagged
#      as "unconverged" in the final report
#   6. Write a TSV summary: design ID, iteration reached, final pLDDT,
#      pass/fail status
#
# The key insight is that a poor ColabFold prediction means the sequence
# doesn't fold to the intended backbone — not that the backbone is bad.
# Re-running ProteinMPNN on the SAME backbone with more diversity often
# finds a sequence that folds cleanly.
#
# USAGE:
#   python feedback2_plddt_mpnn_resample.py config.yaml
#   python feedback2_plddt_mpnn_resample.py config.json
#   python feedback2_plddt_mpnn_resample.py config.txt
#   python feedback2_plddt_mpnn_resample.py \
#       --run_dir=./output/latest \
#       --plddt_pass=75 \
#       --plddt_warn=60 \
#       --max_iterations=3 \
#       --resample_temp=0.2 \
#       --resample_n=16
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir         Pipeline run directory (default: ./output/latest symlink)
#   plddt_pass      Mean pLDDT threshold to accept a design (default: 75)
#   plddt_warn      Below this, print a warning even if passing (default: 60)
#   max_iterations  Maximum resample rounds before marking unconverged (default: 3)
#   resample_temp   ProteinMPNN sampling temperature on retry (default: 0.2,
#                   slightly higher than production 0.1 to explore sequence space)
#   resample_n      Number of new sequences per retry backbone (default: 16)
#   output_dir      Where to write loop outputs (default: <run_dir>/feedback2_loop)
# =============================================================================

import os
import sys
import json
import glob
import shutil
import argparse
import subprocess
import textwrap
from pathlib import Path
from datetime import datetime

# ── optional yaml/json config loading ────────────────────────────────────────
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# ─────────────────────────────────────────────
# Logging — mirrors shell pipeline style
# ─────────────────────────────────────────────
_log_file = None

def log(level: str, msg: str):
    ts  = datetime.now().strftime("%F %T")
    line = f"[{ts}] [{level}] {msg}"
    print(line)
    if _log_file:
        with open(_log_file, "a") as f:
            f.write(line + "\n")

# ─────────────────────────────────────────────
# Checkpoint helpers
# ─────────────────────────────────────────────
def is_done(checkpoint_dir: Path, name: str) -> bool:
    return (checkpoint_dir / f".{name}.done").exists()

def mark_done(checkpoint_dir: Path, name: str, elapsed: float):
    ts = datetime.now().strftime("%F %T")
    (checkpoint_dir / f".{name}.done").write_text(
        f"Completed {ts} | elapsed {elapsed:.0f}s\n"
    )

def mark_failed(checkpoint_dir: Path, name: str, elapsed: float, rc: int):
    ts = datetime.now().strftime("%F %T")
    (checkpoint_dir / f".{name}.failed").write_text(
        f"Failed {ts} | exit={rc} | elapsed {elapsed:.0f}s\n"
    )

# ─────────────────────────────────────────────
# Config loading — yaml / json / txt / argparse
# ─────────────────────────────────────────────
def load_config_file(path: str) -> dict:
    ext = Path(path).suffix.lower()
    with open(path) as f:
        if ext in (".yaml", ".yml"):
            if not HAS_YAML:
                raise ImportError("pyyaml required: pip install pyyaml")
            return yaml.safe_load(f) or {}
        elif ext == ".json":
            return json.load(f)
        else:
            # Plain key: value text
            cfg = {}
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    k, v = line.split(":", 1)
                    cfg[k.strip()] = v.strip()
            return cfg

def parse_args():
    p = argparse.ArgumentParser(
        description="Feedback Loop 2: pLDDT-gated ProteinMPNN resampling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Config file (yaml/json/txt) may be passed as the first positional arg.
            CLI flags always override config file values.
        """)
    )
    p.add_argument("config", nargs="?", help="Path to yaml/json/txt config file")
    p.add_argument("--run_dir",        default=None)
    p.add_argument("--plddt_pass",     type=float, default=None)
    p.add_argument("--plddt_warn",     type=float, default=None)
    p.add_argument("--max_iterations", type=int,   default=None)
    p.add_argument("--resample_temp",  type=float, default=None)
    p.add_argument("--resample_n",     type=int,   default=None)
    p.add_argument("--output_dir",     default=None)
    return p.parse_args()

def resolve_config(args) -> dict:
    # Start with hard defaults
    cfg = {
        "run_dir":        str(Path("./output/latest").resolve()),
        "plddt_pass":     75.0,
        "plddt_warn":     60.0,
        "max_iterations": 3,
        "resample_temp":  0.2,
        "resample_n":     16,
        "output_dir":     None,
    }
    # Layer config file on top
    if args.config and Path(args.config).is_file():
        file_cfg = load_config_file(args.config)
        for k in cfg:
            if k in file_cfg:
                cfg[k] = file_cfg[k]
    # Layer CLI flags on top (only if explicitly provided)
    for k in cfg:
        val = getattr(args, k, None)
        if val is not None:
            cfg[k] = val
    return cfg

# ─────────────────────────────────────────────
# pLDDT parsing
# ─────────────────────────────────────────────
def parse_colabfold_scores(colabfold_dir: Path) -> dict:
    """
    Returns {design_id: {"mean_plddt": float, "pdb": Path, "fasta": Path}}
    Searches recursively for *scores*rank_001*.json files.
    Falls back to any scores json if rank_001 not found.
    """
    pattern = str(colabfold_dir / "**" / "*scores*rank_001*.json")
    jsons   = glob.glob(pattern, recursive=True)
    if not jsons:
        jsons = glob.glob(str(colabfold_dir / "**" / "*scores*.json"), recursive=True)

    results = {}
    for sj in jsons:
        sj = Path(sj)
        try:
            data   = json.loads(sj.read_text())
            plddt  = data.get("plddt", [])
            if not plddt:
                continue
            mean   = sum(plddt) / len(plddt)
            # Infer PDB path (ColabFold may write unrelaxed or relaxed)
            stem   = sj.stem.replace("_scores_rank_001", "").replace("_scores", "")
            pdb    = sj.parent / (sj.name.replace("scores", "unrelaxed").replace(".json", ".pdb"))
            if not pdb.exists():
                pdb = sj.parent / (sj.name.replace("scores", "relaxed").replace(".json", ".pdb"))
            if not pdb.exists():
                pdb = None
            results[stem] = {"mean_plddt": mean, "pdb": pdb, "json": sj}
        except Exception as e:
            log("WARN", f"Could not parse {sj}: {e}")
    return results

# ─────────────────────────────────────────────
# Backbone lookup
# Attempt to find the RFdiffusion PDB that generated a given ProteinMPNN
# sequence by matching design stem names (design_X_dY_seqZ → design_X_dY).
# ─────────────────────────────────────────────
def find_backbone_pdb(design_stem: str, designs_dir: Path) -> Path | None:
    # Strip trailing _seqN suffix if present
    parts = design_stem.rsplit("_seq", 1)
    backbone_stem = parts[0]
    matches = list(designs_dir.glob(f"{backbone_stem}*.pdb"))
    if matches:
        return matches[0]
    # Broader search
    matches = list(designs_dir.rglob(f"*{backbone_stem}*.pdb"))
    return matches[0] if matches else None

# ─────────────────────────────────────────────
# ProteinMPNN runner (subprocess into SE3nv env)
# ─────────────────────────────────────────────
def run_proteinmpnn(backbone_pdb: Path, out_dir: Path,
                    num_seqs: int, temp: float, log_path: Path) -> bool:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", "SE3nv",
        "python", str(Path.home() / "ProteinMPNN" / "protein_mpnn_run.py"),
        "--pdb_path",           str(backbone_pdb),
        "--out_folder",         str(out_dir),
        "--num_seq_per_target", str(num_seqs),
        "--sampling_temp",      str(temp),
        "--batch_size",         "8",
    ]
    log("INFO", f"ProteinMPNN: {backbone_pdb.name} → {out_dir}")
    with open(log_path, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    return rc == 0

# ─────────────────────────────────────────────
# ColabFold runner (subprocess into colabfold env)
# ─────────────────────────────────────────────
def run_colabfold(fasta_dir: Path, out_dir: Path, log_path: Path) -> bool:
    out_dir.mkdir(parents=True, exist_ok=True)
    fastas = list(fasta_dir.rglob("*.fasta"))
    if not fastas:
        log("WARN", f"No FASTA files in {fasta_dir}")
        return False
    all_ok = True
    for fa in fastas:
        cmd = [
            "conda", "run", "-n", "colabfold",
            "colabfold_batch", str(fa), str(out_dir),
            "--num-recycle", "3", "--use-gpu-relax",
        ]
        log("INFO", f"ColabFold: {fa.name}")
        with open(log_path, "a") as lf:
            rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            log("WARN", f"ColabFold exited {rc} for {fa.name}")
            all_ok = False
    return all_ok

# ─────────────────────────────────────────────
# Sequence splitting helper (multi-seq .fa → one .fasta per seq)
# ─────────────────────────────────────────────
def split_fasta(fa_path: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        from Bio import SeqIO
        records = list(SeqIO.parse(str(fa_path), "fasta"))
        for i, rec in enumerate(records):
            out = out_dir / f"{fa_path.stem}_seq{i}.fasta"
            SeqIO.write([rec], str(out), "fasta")
        log("INFO", f"Split {len(records)} seqs from {fa_path.name}")
    except ImportError:
        log("WARN", "Biopython not available for splitting — copying whole file")
        shutil.copy(fa_path, out_dir / fa_path.name)

# ─────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────
def main():
    global _log_file

    args = parse_args()
    cfg  = resolve_config(args)

    run_dir       = Path(cfg["run_dir"]).expanduser().resolve()
    plddt_pass    = float(cfg["plddt_pass"])
    plddt_warn    = float(cfg["plddt_warn"])
    max_iter      = int(cfg["max_iterations"])
    resample_temp = float(cfg["resample_temp"])
    resample_n    = int(cfg["resample_n"])
    output_dir    = Path(cfg["output_dir"]) if cfg["output_dir"] else run_dir / "feedback2_loop"
    output_dir.mkdir(parents=True, exist_ok=True)

    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    _log_file = str(log_dir / "feedback2.log")

    log("INFO", "==========================================")
    log("INFO", "Feedback Loop 2: pLDDT-gated ProteinMPNN resampling")
    log("INFO", f"run_dir        : {run_dir}")
    log("INFO", f"plddt_pass     : {plddt_pass}")
    log("INFO", f"plddt_warn     : {plddt_warn}")
    log("INFO", f"max_iterations : {max_iter}")
    log("INFO", f"resample_temp  : {resample_temp}")
    log("INFO", f"resample_n     : {resample_n}")
    log("INFO", f"output_dir     : {output_dir}")
    log("INFO", "==========================================")

    colab_dir   = run_dir / "colabfold_out"
    designs_dir = run_dir / "designs"

    # ── Initial scoring of Step 6 ColabFold outputs ──────────────────────────
    log("INFO", "Scoring initial ColabFold structures...")
    all_scores = parse_colabfold_scores(colab_dir)
    if not all_scores:
        log("HALT", f"No ColabFold score JSONs found in {colab_dir}")
        sys.exit(1)
    log("INFO", f"Found {len(all_scores)} scored structures")

    # Track per-design state across iterations
    # state: {stem: {"plddt_history": [...], "status": "retry"|"pass"|"unconverged"}}
    state = {}
    for stem, info in all_scores.items():
        m = info["mean_plddt"]
        if m >= plddt_pass:
            status = "pass"
        else:
            status = "retry"
        if m < plddt_warn and status == "pass":
            log("WARN", f"{stem}: pLDDT {m:.1f} passed threshold but below warn level {plddt_warn}")
        state[stem] = {
            "plddt_history": [m],
            "status":        status,
            "backbone_pdb":  find_backbone_pdb(stem, designs_dir),
            "latest_pdb":    info["pdb"],
        }

    n_pass  = sum(1 for v in state.values() if v["status"] == "pass")
    n_retry = sum(1 for v in state.values() if v["status"] == "retry")
    log("INFO", f"Initial pass: {n_pass}  retry: {n_retry}")

    # ─────────────────────────────────────────────
    # Resample iterations
    # ─────────────────────────────────────────────
    for iteration in range(1, max_iter + 1):
        retry_list = [s for s, v in state.items() if v["status"] == "retry"]
        if not retry_list:
            log("INFO", f"All designs passed — stopping at iteration {iteration}")
            break

        log("INFO", f"──────────────────────────────────────────")
        log("INFO", f"Iteration {iteration}: resampling {len(retry_list)} designs")

        iter_dir     = output_dir / f"iter_{iteration}"
        mpnn_dir     = iter_dir / "mpnn_out"
        split_dir    = iter_dir / "split_seqs"
        colab_iter   = iter_dir / "colabfold_out"

        for stem in retry_list:
            backbone = state[stem]["backbone_pdb"]
            if backbone is None or not backbone.exists():
                log("WARN", f"{stem}: backbone PDB not found — marking unconverged")
                state[stem]["status"] = "unconverged"
                continue

            mpnn_stem_dir = mpnn_dir / stem
            mpnn_log      = log_dir / f"iter{iteration}_{stem}_mpnn.log"

            # Run ProteinMPNN at a slightly elevated temperature to diversify
            ok = run_proteinmpnn(backbone, mpnn_stem_dir, resample_n,
                                 resample_temp, mpnn_log)
            if not ok:
                log("WARN", f"{stem}: ProteinMPNN failed at iter {iteration}")
                continue

            # Split sequences for ColabFold
            for fa in mpnn_stem_dir.rglob("*.fa"):
                split_fasta(fa, split_dir / stem)

        # Run ColabFold on all retry sequences from this iteration
        colab_log = log_dir / f"iter{iteration}_colabfold.log"
        run_colabfold(split_dir, colab_iter, colab_log)

        # Re-score
        new_scores = parse_colabfold_scores(colab_iter)
        for stem in retry_list:
            # Find new scores for this stem (may match multiple seqN variants)
            matching = {k: v for k, v in new_scores.items() if stem in k}
            if not matching:
                log("WARN", f"{stem}: no new ColabFold output at iter {iteration}")
                if iteration == max_iter:
                    state[stem]["status"] = "unconverged"
                continue
            # Take the best new pLDDT across all resampled sequences
            best = max(matching.values(), key=lambda x: x["mean_plddt"])
            m    = best["mean_plddt"]
            state[stem]["plddt_history"].append(m)
            if m >= plddt_pass:
                state[stem]["status"]     = "pass"
                state[stem]["latest_pdb"] = best["pdb"]
                log("INFO", f"{stem}: PASSED at iter {iteration} — pLDDT {m:.1f}")
            else:
                log("INFO", f"{stem}: still failing — pLDDT {m:.1f}")
                if iteration == max_iter:
                    state[stem]["status"] = "unconverged"

    # ─────────────────────────────────────────────
    # Final summary report
    # ─────────────────────────────────────────────
    summary_path = output_dir / "feedback2_summary.tsv"
    with open(summary_path, "w") as f:
        f.write("design\tstatus\tfinal_plddt\tplddt_history\tbackbone_pdb\n")
        for stem, v in sorted(state.items()):
            final  = v["plddt_history"][-1]
            hist   = ",".join(f"{x:.1f}" for x in v["plddt_history"])
            bpdb   = str(v["backbone_pdb"]) if v["backbone_pdb"] else "not_found"
            f.write(f"{stem}\t{v['status']}\t{final:.1f}\t{hist}\t{bpdb}\n")

    n_pass  = sum(1 for v in state.values() if v["status"] == "pass")
    n_unc   = sum(1 for v in state.values() if v["status"] == "unconverged")
    n_retry = sum(1 for v in state.values() if v["status"] == "retry")

    log("INFO", "==========================================")
    log("INFO", "Feedback Loop 2 complete.")
    log("INFO", f"  Passed        : {n_pass}")
    log("INFO", f"  Unconverged   : {n_unc}")
    log("INFO", f"  Still retrying: {n_retry}")
    log("INFO", f"  Summary       : {summary_path}")
    log("INFO", "==========================================")

    # Exit non-zero if any designs never converged (useful for CI / snakemake)
    if n_unc > 0:
        log("WARN", f"{n_unc} design(s) never reached pLDDT >= {plddt_pass} "
                    f"after {max_iter} iterations — see summary for details")
        sys.exit(2)

if __name__ == "__main__":
    main()
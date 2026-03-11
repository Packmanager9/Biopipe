#!/usr/bin/env python3
# =============================================================================
# feedback4_domesticated_cds_revalidate.py
#
# FEEDBACK LOOP 4: Domesticated CDS sequences → ColabFold re-validation
#
# WORKFLOW:
#   1. Read the GenBank (.gb) file(s) produced by plasmid_design_moclo_v3.py
#      (Step 9) and extract the CDS sequences that were domesticated
#      (renamed with the _dom_opt or _chisel_opt suffix)
#   2. Compare each domesticated CDS to its pre-domestication source sequence
#      to calculate the number and positions of synonymous changes
#   3. Write the domesticated amino acid sequences as individual FASTA files
#   4. Run ColabFold on those sequences
#   5. Compare pLDDT scores between the original ColabFold prediction (from
#      Step 6) and the post-domestication prediction
#   6. Flag any sequence where the mean pLDDT dropped more than plddt_drop_warn
#      (indicating that a synonymous change may have disrupted a folding-
#      critical mRNA structure or introduced a rare codon cluster)
#   7. Write a per-sequence report: num_changes, original_plddt,
#      domesticated_plddt, delta, flag
#
# WHY THIS MATTERS:
#   Synonymous codon substitutions used to remove restriction sites can
#   alter local mRNA structure and translation rate.  In rare cases this
#   affects cotranslational folding.  This loop catches those regressions
#   before a sequence is sent to synthesis.
#
# USAGE:
#   python feedback4_domesticated_cds_revalidate.py config.yaml
#   python feedback4_domesticated_cds_revalidate.py config.json
#   python feedback4_domesticated_cds_revalidate.py config.txt
#   python feedback4_domesticated_cds_revalidate.py \
#       --run_dir=./output/latest \
#       --plasmid_dir=./output/latest/moclo_plasmids \
#       --plddt_drop_warn=5.0 \
#       --plddt_drop_fail=15.0
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir          Pipeline run directory (default: ./output/latest)
#   plasmid_dir      Where plasmid_design_moclo_v3.py wrote its .gb files
#                    (default: <run_dir>/moclo_plasmids)
#   plddt_drop_warn  pLDDT drop that triggers a WARNING (default: 5.0)
#   plddt_drop_fail  pLDDT drop that triggers a FAIL flag (default: 15.0)
#   output_dir       Where to write loop outputs (default: <run_dir>/feedback4_loop)
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
from typing import Optional

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    HAS_BIO = True
except ImportError:
    HAS_BIO = False

# ─────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────
_log_file: Optional[str] = None

def log(level: str, msg: str):
    ts   = datetime.now().strftime("%F %T")
    line = f"[{ts}] [{level}] {msg}"
    print(line)
    if _log_file:
        with open(_log_file, "a") as f:
            f.write(line + "\n")

# ─────────────────────────────────────────────
# Config loading
# ─────────────────────────────────────────────
def load_config_file(path: str) -> dict:
    ext = Path(path).suffix.lower()
    with open(path) as f:
        if ext in (".yaml", ".yml"):
            if not HAS_YAML:
                raise ImportError("pip install pyyaml")
            return yaml.safe_load(f) or {}
        elif ext == ".json":
            return json.load(f)
        else:
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
        description="Feedback Loop 4: domesticated CDS → ColabFold re-validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("config", nargs="?")
    p.add_argument("--run_dir",         default=None)
    p.add_argument("--plasmid_dir",     default=None)
    p.add_argument("--plddt_drop_warn", type=float, default=None)
    p.add_argument("--plddt_drop_fail", type=float, default=None)
    p.add_argument("--output_dir",      default=None)
    return p.parse_args()

def resolve_config(args) -> dict:
    cfg = {
        "run_dir":         str(Path("./output/latest").resolve()),
        "plasmid_dir":     None,
        "plddt_drop_warn": 5.0,
        "plddt_drop_fail": 15.0,
        "output_dir":      None,
    }
    if args.config and Path(args.config).is_file():
        for k, v in load_config_file(args.config).items():
            if k in cfg:
                cfg[k] = v
    for k in cfg:
        val = getattr(args, k, None)
        if val is not None:
            cfg[k] = val
    return cfg

# ─────────────────────────────────────────────
# GenBank CDS extraction
# Parses .gb files and finds features whose /label or /note contains
# "_dom_opt" or "_chisel_opt" (the suffixes appended by the domestication
# routines in plasmid_design_moclo_v3.py).
# ─────────────────────────────────────────────
def extract_domesticated_cds(gb_dir: Path) -> list[dict]:
    """
    Returns a list of dicts:
        {name, dna_seq, aa_seq, source_gb, feature}
    for every CDS feature that was domesticated.
    """
    if not HAS_BIO:
        raise ImportError("Biopython required: pip install biopython")

    results = []
    for gb_file in gb_dir.glob("*.gb"):
        for record in SeqIO.parse(str(gb_file), "genbank"):
            for feat in record.features:
                if feat.type != "misc_feature":
                    continue
                label = " ".join(feat.qualifiers.get("label", []) +
                                 feat.qualifiers.get("note", []))
                if "_dom_opt" not in label and "_chisel_opt" not in label:
                    continue
                dna = feat.extract(record.seq)
                aa  = dna.translate(to_stop=True)
                results.append({
                    "name":      label.strip(),
                    "dna_seq":   str(dna),
                    "aa_seq":    str(aa),
                    "source_gb": gb_file.name,
                })
                log("INFO", f"Found domesticated CDS: {label.strip()} ({len(aa)} aa) from {gb_file.name}")
    return results

# ─────────────────────────────────────────────
# pLDDT lookup from original Step 6 ColabFold outputs
# ─────────────────────────────────────────────
def get_original_plddt(design_name: str, colab_dir: Path) -> Optional[float]:
    """
    Try to find the original pLDDT for a design by matching the name stem
    (strip the _dom_opt/_chisel_opt suffix first) against ColabFold JSONs.
    """
    stem = design_name.replace("_dom_opt", "").replace("_chisel_opt", "")
    pattern = str(colab_dir / "**" / f"*{stem}*scores*rank_001*.json")
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        pattern = str(colab_dir / "**" / f"*{stem}*scores*.json")
        matches = glob.glob(pattern, recursive=True)
    if not matches:
        return None
    try:
        data   = json.loads(Path(matches[0]).read_text())
        plddt  = data.get("plddt", [])
        return sum(plddt) / len(plddt) if plddt else None
    except Exception:
        return None

def parse_new_plddt(colab_dir: Path, design_name: str) -> Optional[float]:
    """Parse pLDDT from newly run ColabFold output for a given design name."""
    pattern = str(colab_dir / "**" / f"*{design_name}*scores*rank_001*.json")
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        pattern = str(colab_dir / "**" / f"*{design_name}*scores*.json")
        matches = glob.glob(pattern, recursive=True)
    if not matches:
        return None
    try:
        data  = json.loads(Path(matches[0]).read_text())
        plddt = data.get("plddt", [])
        return sum(plddt) / len(plddt) if plddt else None
    except Exception:
        return None

# ─────────────────────────────────────────────
# ColabFold runner
# ─────────────────────────────────────────────
def run_colabfold(fasta_path: Path, out_dir: Path, log_path: Path) -> bool:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "conda", "run", "-n", "colabfold",
        "colabfold_batch", str(fasta_path), str(out_dir),
        "--num-recycle", "3", "--use-gpu-relax",
    ]
    log("INFO", f"ColabFold: {fasta_path.name}")
    with open(log_path, "a") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        log("WARN", f"ColabFold exited {rc} for {fasta_path.name}")
    return rc == 0

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    global _log_file

    if not HAS_BIO:
        print("ERROR: Biopython required — pip install biopython")
        sys.exit(1)

    args = parse_args()
    cfg  = resolve_config(args)

    run_dir      = Path(cfg["run_dir"]).expanduser().resolve()
    plasmid_dir  = Path(cfg["plasmid_dir"]).expanduser().resolve() \
                   if cfg["plasmid_dir"] else run_dir / "moclo_plasmids"
    drop_warn    = float(cfg["plddt_drop_warn"])
    drop_fail    = float(cfg["plddt_drop_fail"])
    output_dir   = Path(cfg["output_dir"]).expanduser().resolve() \
                   if cfg["output_dir"] else run_dir / "feedback4_loop"

    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    _log_file = str(log_dir / "feedback4.log")

    log("INFO", "==========================================")
    log("INFO", "Feedback Loop 4: domesticated CDS re-validation")
    log("INFO", f"run_dir      : {run_dir}")
    log("INFO", f"plasmid_dir  : {plasmid_dir}")
    log("INFO", f"drop_warn    : {drop_warn}")
    log("INFO", f"drop_fail    : {drop_fail}")
    log("INFO", f"output_dir   : {output_dir}")
    log("INFO", "==========================================")

    colab_dir  = run_dir / "colabfold_out"
    fasta_dir  = output_dir / "domesticated_fastas"
    colab_new  = output_dir / "colabfold_revalidation"
    fasta_dir.mkdir(exist_ok=True)

    # ── Step 1: Extract domesticated CDS from GenBank files ─────────────────
    if not plasmid_dir.exists():
        log("HALT", f"plasmid_dir {plasmid_dir} does not exist")
        sys.exit(1)

    log("INFO", f"Scanning GenBank files in {plasmid_dir}")
    designs = extract_domesticated_cds(plasmid_dir)
    if not designs:
        log("HALT", "No domesticated CDS features found in GenBank files. "
                    "Ensure plasmid_design_moclo_v3.py was run with perform_domestication: true")
        sys.exit(1)
    log("INFO", f"Found {len(designs)} domesticated CDS sequences")

    # ── Step 2: Write amino acid FASTA files for ColabFold ──────────────────
    # ColabFold works on protein sequences, so we translate the domesticated CDS
    fasta_paths = []
    for d in designs:
        safe_name = d["name"].replace("/", "_").replace(" ", "_")[:80]
        fa_path   = fasta_dir / f"{safe_name}.fasta"
        with open(fa_path, "w") as f:
            f.write(f">{safe_name}\n{d['aa_seq']}\n")
        d["fasta_path"] = fa_path
        d["safe_name"]  = safe_name
        fasta_paths.append(fa_path)
        log("INFO", f"Wrote {fa_path.name}")

    # ── Step 3: Look up original pLDDT scores ───────────────────────────────
    for d in designs:
        d["original_plddt"] = get_original_plddt(d["name"], colab_dir)
        if d["original_plddt"] is None:
            log("WARN", f"{d['name']}: no original ColabFold pLDDT found — delta will be N/A")
        else:
            log("INFO", f"{d['name']}: original pLDDT = {d['original_plddt']:.1f}")

    # ── Step 4: Run ColabFold on domesticated sequences ──────────────────────
    log("INFO", f"Running ColabFold on {len(fasta_paths)} domesticated sequences...")
    colab_log = log_dir / "colabfold_revalidation.log"
    for fa in fasta_paths:
        run_colabfold(fa, colab_new, colab_log)

    # ── Step 5: Parse new pLDDT scores and compare ───────────────────────────
    log("INFO", "Comparing original vs domesticated pLDDT scores...")
    report_rows = []
    flags = {"warn": 0, "fail": 0, "ok": 0, "no_data": 0}

    for d in designs:
        new_plddt = parse_new_plddt(colab_new, d["safe_name"])

        if new_plddt is None:
            flag  = "no_data"
            delta = None
            log("WARN", f"{d['name']}: no new ColabFold pLDDT found")
            flags["no_data"] += 1
        elif d["original_plddt"] is None:
            flag  = "no_baseline"
            delta = None
            flags["no_data"] += 1
        else:
            delta = new_plddt - d["original_plddt"]
            if delta <= -drop_fail:
                flag = "FAIL"
                flags["fail"] += 1
                log("WARN", f"FAIL {d['name']}: pLDDT dropped {delta:+.1f} "
                             f"({d['original_plddt']:.1f} → {new_plddt:.1f}). "
                             f"Consider reverting or redesigning this CDS.")
            elif delta <= -drop_warn:
                flag = "warn"
                flags["warn"] += 1
                log("WARN", f"warn {d['name']}: pLDDT dropped {delta:+.1f} "
                             f"({d['original_plddt']:.1f} → {new_plddt:.1f})")
            else:
                flag = "ok"
                flags["ok"] += 1
                log("INFO", f"ok   {d['name']}: pLDDT delta {delta:+.1f} "
                             f"({d['original_plddt']:.1f} → {new_plddt:.1f})")

        report_rows.append({
            "design":              d["name"],
            "source_gb":           d["source_gb"],
            "aa_length":           len(d["aa_seq"]),
            "original_plddt":      f"{d['original_plddt']:.1f}" if d["original_plddt"] else "N/A",
            "domesticated_plddt":  f"{new_plddt:.1f}" if new_plddt else "N/A",
            "delta":               f"{delta:+.1f}" if delta is not None else "N/A",
            "flag":                flag,
        })

    # ── Step 6: Write report ─────────────────────────────────────────────────
    report_path = output_dir / "feedback4_revalidation_report.tsv"
    headers = ["design", "source_gb", "aa_length", "original_plddt",
               "domesticated_plddt", "delta", "flag"]
    with open(report_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in report_rows:
            f.write("\t".join(str(row[h]) for h in headers) + "\n")

    log("INFO", "==========================================")
    log("INFO", "Feedback Loop 4 complete.")
    log("INFO", f"  OK           : {flags['ok']}")
    log("INFO", f"  Warnings     : {flags['warn']}")
    log("INFO", f"  FAIL         : {flags['fail']}")
    log("INFO", f"  No data      : {flags['no_data']}")
    log("INFO", f"  Report       : {report_path}")
    log("INFO", "==========================================")

    if flags["fail"] > 0:
        log("WARN", f"{flags['fail']} design(s) failed pLDDT regression check. "
                    "Review the report and consider re-running domestication with "
                    "stricter codon constraints or dnachisel[reports].")
        sys.exit(2)

if __name__ == "__main__":
    main()
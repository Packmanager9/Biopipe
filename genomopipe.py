#!/usr/bin/env python3
# =============================================================================
# genomopipe.py — Master pipeline orchestrator
#
# Runs the full Genomopipe stack in order:
#
#   Phase 1 — genome_to_design.sh         (genome → proteins → structures)
#   Phase 2 — plasmid_design_moclo_v3.py  (Step 9: MoClo plasmid design)
#   Phase 3 — Feedback loops (in order):
#       Loop 6 — feedback6_blast_taxonomy_rerun.py
#                (BLAST taxonomy → BRAKER re-run with corrected OrthoDB partition)
#       [Future loops will be added here as they are implemented]
#
# USAGE:
#   python genomopipe.py config.yaml
#   python genomopipe.py config.json
#   python genomopipe.py config.txt
#
#   # Or drive entirely from CLI flags (every flag is optional):
#   python genomopipe.py \
#       --organism "Taraxacum officinale" \
#       --output_dir ./output \
#       --is_eukaryote true \
#       --auto_rnaseq \
#       --moclo_standard marillonnet \
#       --genes design_0_seq0.fasta design_1_seq3.fasta
#
#   # Mix: load a base config, then override specific keys on the CLI:
#   python genomopipe.py config.yaml --dry_run --skip_phase1
#
# RESUMABILITY:
#   Each phase writes a sentinel file in <output_dir>/latest/:
#     .genomopipe_phase1.done
#     .genomopipe_phase2.done
#     .genomopipe_feedback6.done
#   Re-running will skip phases that already completed. Use --force to
#   clear all sentinels and restart from scratch, or --skip_phase1 /
#   --skip_phase2 / --skip_feedback to selectively skip.
#
# CONFIG FILE FORMAT (YAML shown; JSON and plain-text .txt also accepted):
#
#   # ── Shared / global ─────────────────────────────────────────────────
#   output_dir:    ./output
#   scripts_dir:   .          # directory containing all .py and .sh scripts
#   email:         you@example.com
#
#   # ── Phase 1: genome_to_design.sh ────────────────────────────────────
#   organism:      "Taraxacum officinale"
#   is_eukaryote:  true
#   genemark_path: ~/genemark-etp-full/gmetp_linux_64/bin
#   bam:           ~          # path to existing BAM, or omit
#   auto_rnaseq:   false
#   force:         false      # --force: clear all genome_to_design checkpoints
#
#   # ── Phase 2: plasmid_design_moclo_v3.py ─────────────────────────────
#   assembly_method:       GoldenGate
#   moclo_standard:        marillonnet
#   enzyme_level0:         BsaI-HFv2
#   enzyme_level1:         BpiI
#   promoter:              pJ23119
#   rbs:                   RBS_strong
#   terminator:            rrnB_T1
#   ori:                   pMB1_ori
#   marker:                KanR
#   perform_domestication: true
#   output_prefix:         moclo_plasmid
#   genes: []              # list of .fasta paths; if empty, auto-discovered
#                          # from proteinmpnn_out/split_seqs in the latest run
#
#   # ── Phase 3: Feedback loops ──────────────────────────────────────────
#   fb6_min_hits:      5
#   fb6_evalue_cutoff: 1.0e-5
#   fb6_dry_run:       false
#
#   # ── Orchestrator control ─────────────────────────────────────────────
#   skip_phase1:    false
#   skip_phase2:    false
#   skip_feedback:  false   # skips ALL feedback loops
#   # Fine-grained feedback loop control:
#   run_fb6:        true
# =============================================================================

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

# ─────────────────────────────────────────────────────────────────────────────
# Defaults
# ─────────────────────────────────────────────────────────────────────────────

DEFAULTS: Dict[str, Any] = {
    # Shared
    "output_dir":    "./output",
    "scripts_dir":   ".",
    "email":         "PricklyPearEnterprises@gmail.com",

    # Phase 1 — genome_to_design.sh
    "organism":      "",
    "is_eukaryote":  False,
    "genemark_path": str(Path.home() / "genemark-etp-full/gmetp_linux_64/bin"),
    "bam":           None,
    "auto_rnaseq":   False,
    "force":         False,

    # Phase 2 — plasmid_design_moclo_v3.py
    "assembly_method":       "GoldenGate",
    "moclo_standard":        "marillonnet",
    "enzyme_level0":         "BsaI-HFv2",
    "enzyme_level1":         "BpiI",
    "promoter":              "pJ23119",
    "rbs":                   "RBS_strong",
    "terminator":            "rrnB_T1",
    "ori":                   "pMB1_ori",
    "marker":                "KanR",
    "backbone":              None,
    "perform_domestication": True,
    "output_prefix":         "moclo_plasmid",
    "genes":                 [],

    # Phase 3 — feedback loops
    "fb6_min_hits":      5,
    "fb6_evalue_cutoff": 1e-5,
    "fb6_dry_run":       False,

    # Orchestrator control
    "skip_phase1":   False,
    "skip_phase2":   False,
    "skip_feedback": False,
    "run_fb1": True,
    "run_fb2": True,
    "run_fb3": True,
    "run_fb4": True,
    "run_fb5": True,
    "run_fb6": True,
}

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────

_log_path: Optional[Path] = None


def log(level: str, msg: str) -> None:
    ts   = datetime.now().strftime("%F %T")
    line = f"[{ts}] [ORCHESTRATOR] [{level}] {msg}"
    print(line, flush=True)
    if _log_path:
        with open(_log_path, "a") as f:
            f.write(line + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# Config helpers
# ─────────────────────────────────────────────────────────────────────────────

def _to_bool(val: Any) -> bool:
    if isinstance(val, bool):
        return val
    if isinstance(val, int):
        return bool(val)
    s = str(val).strip().lower()
    if s in ("true", "1", "yes"):
        return True
    if s in ("false", "0", "no", "none", "null", "~", ""):
        return False
    raise ValueError(f"Cannot parse '{val}' as a boolean")


def load_config_file(path: Path) -> Dict[str, Any]:
    ext = path.suffix.lower()
    with open(path) as f:
        if ext in (".yaml", ".yml"):
            if not HAS_YAML:
                raise ImportError("pip install pyyaml")
            return yaml.safe_load(f) or {}
        elif ext == ".json":
            return json.load(f)
        else:
            # Plain-text key: value parser
            cfg: Dict[str, Any] = {}
            section: Optional[str] = None
            for raw in f:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if line.endswith(":") and " " not in line:
                    section = line[:-1].strip().lower()
                    cfg.setdefault(section, [])
                    continue
                if ":" in line:
                    k, v = line.split(":", 1)
                    cfg[k.strip()] = v.strip()
                    section = None
                elif section is not None:
                    cfg[section].append(line)
            return cfg


def build_config(args: argparse.Namespace) -> Dict[str, Any]:
    """Merge: DEFAULTS → config file → CLI flags (highest priority)."""
    cfg = dict(DEFAULTS)

    if args.config and Path(args.config).is_file():
        for k, v in load_config_file(Path(args.config)).items():
            cfg[k] = v

    # Map CLI flags into cfg — only when explicitly provided (not None)
    cli_map = {
        "organism":      args.organism,
        "output_dir":    args.output_dir,
        "scripts_dir":   args.scripts_dir,
        "is_eukaryote":  args.is_eukaryote,
        "genemark_path": args.genemark_path,
        "bam":           args.bam,
        "email":         args.email,
        "moclo_standard":        args.moclo_standard,
        "enzyme_level0":         args.enzyme_level0,
        "enzyme_level1":         args.enzyme_level1,
        "perform_domestication": args.perform_domestication,
        "output_prefix":         args.output_prefix,
        "fb6_min_hits":          args.fb6_min_hits,
        "fb6_evalue_cutoff":     args.fb6_evalue_cutoff,
    }
    for k, v in cli_map.items():
        if v is not None:
            cfg[k] = v

    # Boolean store_true flags
    if args.auto_rnaseq:    cfg["auto_rnaseq"]   = True
    if args.force:          cfg["force"]          = True
    if args.skip_phase1:    cfg["skip_phase1"]    = True
    if args.skip_phase2:    cfg["skip_phase2"]    = True
    if args.skip_feedback:  cfg["skip_feedback"]  = True
    if args.dry_run:        cfg["fb6_dry_run"]    = True
    if args.no_fb6:         cfg["run_fb6"]        = False
    if args.genes:          cfg["genes"]          = args.genes

    # Normalise booleans that may have come in as strings from a txt config
    for bool_key in ("is_eukaryote", "auto_rnaseq", "force",
                     "perform_domestication", "skip_phase1", "skip_phase2",
                     "skip_feedback", "fb6_dry_run", "run_fb6"):
        cfg[bool_key] = _to_bool(cfg[bool_key])

    cfg["output_dir"] = str(Path(cfg["output_dir"]).expanduser().resolve())
    cfg["scripts_dir"] = str(Path(cfg["scripts_dir"]).expanduser().resolve())

    return cfg


# ─────────────────────────────────────────────────────────────────────────────
# Sentinel / checkpoint helpers
# ─────────────────────────────────────────────────────────────────────────────

def _sentinel(run_dir: Path, name: str) -> Path:
    return run_dir / f".genomopipe_{name}.done"


def _is_done(run_dir: Path, name: str) -> bool:
    return _sentinel(run_dir, name).exists()


def _mark_done(run_dir: Path, name: str) -> None:
    _sentinel(run_dir, name).write_text(
        f"Completed {datetime.now().strftime('%F %T')}\n"
    )


def _clear_sentinels(run_dir: Path) -> None:
    for f in run_dir.glob(".genomopipe_*.done"):
        f.unlink()
    log("INFO", "Cleared all orchestrator sentinels")


def _resolve_latest_run(output_dir: Path) -> Optional[Path]:
    link = output_dir / "latest"
    if link.is_symlink() and link.exists():
        return link.resolve()
    runs = sorted(output_dir.glob("run_*"),
                  key=lambda p: p.stat().st_mtime, reverse=True)
    return runs[0] if runs else None


# ─────────────────────────────────────────────────────────────────────────────
# Script locator
# ─────────────────────────────────────────────────────────────────────────────

def _find_script(scripts_dir: str, name: str) -> Path:
    """Locate a script by searching scripts_dir then PATH."""
    candidates = [
        Path(scripts_dir) / name,
        Path(__file__).parent / name,
        Path(shutil.which(name) or ""),
    ]
    for p in candidates:
        if p.is_file():
            return p.resolve()
    raise FileNotFoundError(
        f"Cannot find '{name}'. Set scripts_dir in config or place it alongside genomopipe.py."
    )


# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 — genome_to_design.sh
# ─────────────────────────────────────────────────────────────────────────────

def run_phase1(cfg: Dict[str, Any], run_dir_ref: List[Optional[Path]]) -> Path:
    """
    Invoke genome_to_design.sh and return the resolved run directory.
    run_dir_ref is a one-element list used to communicate the run_dir
    back to the caller even if this phase is skipped.
    """
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    if cfg["skip_phase1"]:
        log("SKIP", "Phase 1 (genome_to_design.sh) — --skip_phase1 set")
        run_dir = _resolve_latest_run(output_dir)
        if run_dir is None:
            log("HALT", "skip_phase1 set but no existing run found in output_dir")
            sys.exit(1)
        run_dir_ref[0] = run_dir
        return run_dir

    run_dir = _resolve_latest_run(output_dir)

    # Check orchestrator-level sentinel (distinct from genome_to_design checkpoints)
    if run_dir and _is_done(run_dir, "phase1"):
        log("SKIP", f"Phase 1 already complete ({run_dir})")
        run_dir_ref[0] = run_dir
        return run_dir

    script = _find_script(cfg["scripts_dir"], "genome_to_design.sh")
    organism = cfg["organism"]
    if not organism:
        log("HALT", "No organism specified. Set 'organism' in config or pass --organism.")
        sys.exit(1)

    cmd = [
        "bash", str(script),
        organism,
        cfg["output_dir"],
        "true" if cfg["is_eukaryote"] else "false",
        f"--GENEMARK_PATH={cfg['genemark_path']}",
    ]
    if cfg.get("bam"):
        cmd.append(f"--bam={cfg['bam']}")
    if cfg["auto_rnaseq"]:
        cmd.append("--auto_rnaseq")
    if cfg["force"]:
        cmd.append("--force")

    log("START", f"Phase 1 — genome_to_design.sh")
    log("INFO",  f"  Command: {' '.join(cmd)}")

    phase1_log = output_dir / "genomopipe_phase1.log"
    with open(phase1_log, "w") as lf:
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        log("FAIL",
            f"Phase 1 exited {result.returncode}. "
            f"See {phase1_log} and {output_dir}/latest/logs/pipeline.log")
        sys.exit(result.returncode)

    run_dir = _resolve_latest_run(output_dir)
    if run_dir is None:
        log("HALT", "Phase 1 completed but no run directory found under output_dir")
        sys.exit(1)

    _mark_done(run_dir, "phase1")
    run_dir_ref[0] = run_dir
    log("DONE", f"Phase 1 complete — run dir: {run_dir}")
    return run_dir


# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — plasmid_design_moclo_v3.py
# ─────────────────────────────────────────────────────────────────────────────

def _build_moclo_config(cfg: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    """Construct the dict that plasmid_design_moclo_v3.py expects."""
    genes = cfg.get("genes") or []

    # Auto-discover from split_seqs if no genes specified
    if not genes:
        split_seqs = run_dir / "proteinmpnn_out" / "split_seqs"
        if split_seqs.is_dir():
            genes = sorted(str(p) for p in split_seqs.rglob("*.fasta"))
            log("INFO", f"Phase 2: auto-discovered {len(genes)} gene FASTA(s) "
                        f"from {split_seqs}")
        else:
            log("WARN", "Phase 2: no genes specified and split_seqs not found; "
                        "assembly will have no CDS parts")

    moclo: Dict[str, Any] = {
        "assembly_method":       cfg["assembly_method"],
        "moclo_standard":        cfg["moclo_standard"],
        "enzyme_level0":         cfg["enzyme_level0"],
        "enzyme_level1":         cfg["enzyme_level1"],
        "promoter":              cfg["promoter"],
        "rbs":                   cfg["rbs"],
        "terminator":            cfg["terminator"],
        "ori":                   cfg["ori"],
        "marker":                cfg["marker"],
        "perform_domestication": cfg["perform_domestication"],
        "output_prefix":         cfg["output_prefix"],
        "run_dir":               str(run_dir),
        "output_dir":            str(run_dir / "moclo_plasmids"),
        "genes":                 genes,
    }
    if cfg.get("backbone"):
        moclo["backbone"] = cfg["backbone"]

    return moclo


def run_phase2(cfg: Dict[str, Any], run_dir: Path) -> None:
    if cfg["skip_phase2"]:
        log("SKIP", "Phase 2 (plasmid_design_moclo_v3.py) — --skip_phase2 set")
        return

    if _is_done(run_dir, "phase2"):
        log("SKIP", "Phase 2 already complete")
        return

    script = _find_script(cfg["scripts_dir"], "plasmid_design_moclo_v3.py")

    moclo_cfg = _build_moclo_config(cfg, run_dir)

    # Write a temporary YAML config for the plasmid script so it receives
    # the full structured config (including nested lists for genes)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False, dir=run_dir
    ) as tmp:
        if HAS_YAML:
            yaml.dump(moclo_cfg, tmp, default_flow_style=False)
        else:
            json.dump(moclo_cfg, tmp, indent=2)
            tmp.name  # reassign suffix handled below
        tmp_path = tmp.name

    # Rename to .yaml if we used yaml, .json otherwise
    if not HAS_YAML:
        json_path = tmp_path.replace(".yaml", ".json")
        Path(tmp_path).rename(json_path)
        tmp_path = json_path

    log("START", "Phase 2 — plasmid_design_moclo_v3.py")
    log("INFO",  f"  Config: {tmp_path}")
    log("INFO",  f"  Genes:  {len(moclo_cfg['genes'])} file(s)")

    phase2_log = run_dir / "genomopipe_phase2.log"
    with open(phase2_log, "w") as lf:
        result = subprocess.run(
            [sys.executable, str(script), tmp_path],
            stdout=lf, stderr=subprocess.STDOUT
        )

    if result.returncode != 0:
        log("FAIL", f"Phase 2 exited {result.returncode}. See {phase2_log}")
        # Non-critical — continue to feedback loops
    else:
        _mark_done(run_dir, "phase2")
        log("DONE", "Phase 2 complete")


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 — Feedback loops
# ─────────────────────────────────────────────────────────────────────────────

def _build_fb6_config(cfg: Dict[str, Any], run_dir: Path) -> Dict[str, Any]:
    return {
        "run_dir":       str(run_dir),
        "min_hits":      cfg["fb6_min_hits"],
        "evalue_cutoff": cfg["fb6_evalue_cutoff"],
        "dry_run":       cfg["fb6_dry_run"],
        "genemark_path": cfg["genemark_path"],
        "is_eukaryote":  cfg["is_eukaryote"],
        "email":         cfg["email"],
        "output_dir":    str(run_dir / "feedback6_loop"),
    }


def run_feedback_loops(cfg: Dict[str, Any], run_dir: Path) -> None:
    if cfg["skip_feedback"]:
        log("SKIP", "All feedback loops — --skip_feedback set")
        return

    # ── Helper ────────────────────────────────────────────────────────────────
    def _run_loop(sentinel: str, label: str, script_name: str, lang: str,
                  extra_args: Optional[List[str]] = None) -> bool:
        """Run one feedback loop script. Returns True on success."""
        if _is_done(run_dir, sentinel):
            log("SKIP", f"{label} already complete")
            return True
        try:
            script = _find_script(cfg["scripts_dir"], script_name)
        except FileNotFoundError as e:
            log("SKIP", f"{label} — {e}")
            return False

        cmd = ([sys.executable, str(script)] if lang == "py"
               else ["bash", str(script)])
        cmd += ["--run_dir", str(run_dir)]
        if extra_args:
            cmd += extra_args

        loop_log = run_dir / f"genomopipe_{sentinel}.log"
        log("START", f"{label}")
        log("INFO",  f"  Command: {' '.join(cmd)}")
        with open(loop_log, "w") as lf:
            result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            log("FAIL", f"{label} exited {result.returncode}. See {loop_log}")
            return False
        _mark_done(run_dir, sentinel)
        log("DONE", f"{label} complete")
        return True

    # ── FB1: ColabFold quality → RFdiffusion resample ─────────────────────────
    if cfg.get("run_fb1", True):
        _run_loop("feedback1", "Feedback Loop 1 — ColabFold → RFdiffusion",
                  "feedback1_colabfold_to_rfdiffusion.sh", "sh")

    # ── FB2: pLDDT / ProteinMPNN resample ─────────────────────────────────────
    if cfg.get("run_fb2", True):
        _run_loop("feedback2", "Feedback Loop 2 — pLDDT / ProteinMPNN resample",
                  "feedback2_plddt_mpnn_resample.py", "py")

    # ── FB3: BLAST taxonomy → BRAKER re-annotation ────────────────────────────
    if cfg.get("run_fb3", True):
        _run_loop("feedback3", "Feedback Loop 3 — BLAST → BRAKER re-run",
                  "feedback3_blast_to_braker.sh", "sh")

    # ── FB4: Domesticated CDS revalidation ────────────────────────────────────
    if cfg.get("run_fb4", True):
        _run_loop("feedback4", "Feedback Loop 4 — domesticated CDS revalidation",
                  "feedback4_domesticated_cds_revalidate.py", "py")

    # ── FB5: Designed proteins → annotation ───────────────────────────────────
    if cfg.get("run_fb5", True):
        _run_loop("feedback5", "Feedback Loop 5 — designed proteins → annotation",
                  "feedback5_designed_proteins_to_annotation.sh", "sh")

    # ── FB6: BLAST taxonomy → OrthoDB partition fix ───────────────────────────
    if cfg.get("run_fb6", True):
        blast_file = run_dir / "blast_results.txt"
        if not blast_file.exists():
            log("SKIP", "Feedback Loop 6 — blast_results.txt not found. Skipping.")
        else:
            _run_loop(
                "feedback6", "Feedback Loop 6 — BLAST taxonomy re-run",
                "feedback6_blast_taxonomy_rerun.py", "py",
                extra_args=[
                    "--min_hits",      str(cfg["fb6_min_hits"]),
                    "--evalue_cutoff", str(cfg["fb6_evalue_cutoff"]),
                ] + (["--dry_run"] if cfg["fb6_dry_run"] else [])
            )
            audit = run_dir / "feedback6_loop" / "feedback6_taxonomy_audit.txt"
            if audit.exists():
                log("INFO", f"  Audit report: {audit}")


# ─────────────────────────────────────────────────────────────────────────────
# Summary report
# ─────────────────────────────────────────────────────────────────────────────

def write_summary(cfg: Dict[str, Any], run_dir: Path) -> None:
    lines = [
        "# Genomopipe Orchestrator — Run Summary",
        "",
        f"Timestamp  : {datetime.now().strftime('%F %T')}",
        f"Organism   : {cfg['organism']}",
        f"Output dir : {cfg['output_dir']}",
        f"Run dir    : {run_dir}",
        "",
        "## Phase completion",
    ]
    phases = [
        ("phase1",    "Phase 1 — genome_to_design.sh"),
        ("phase2",    "Phase 2 — plasmid_design_moclo_v3.py"),
        ("feedback6", "Feedback Loop 6 — BLAST taxonomy re-run"),
    ]
    for sentinel_name, label in phases:
        s = _sentinel(run_dir, sentinel_name)
        status = ("✓ done — " + s.read_text().strip()) if s.exists() else "— not run / failed"
        lines.append(f"  {label}: {status}")

    lines += [
        "",
        "## Key outputs",
        f"  Annotations    : {run_dir}/braker_out/",
        f"  Designs        : {run_dir}/designs/",
        f"  Structures     : {run_dir}/colabfold_out/",
        f"  Plasmid files  : {run_dir}/moclo_plasmids/",
        f"  BLAST results  : {run_dir}/blast_results.txt",
        f"  Feedback 6     : {run_dir}/feedback6_loop/",
        f"  Orchestrator log: {_log_path}",
    ]

    summary_path = run_dir / "genomopipe_summary.md"
    summary_path.write_text("\n".join(lines) + "\n")
    log("INFO", f"Summary written: {summary_path}")
    for line in lines:
        log("INFO", line)


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="genomopipe.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
            Genomopipe master orchestrator.
            Runs genome_to_design.sh → plasmid_design_moclo_v3.py → feedback loops.

            A config file (YAML / JSON / .txt) is recommended for complex runs.
            All config keys can also be passed as CLI flags shown below.
            CLI flags always override config file values.
        """),
    )

    p.add_argument("config", nargs="?",
                   help="Config file (.yaml / .json / .txt)")

    # ── Shared ────────────────────────────────────────────────────────────
    g = p.add_argument_group("Shared / global")
    g.add_argument("--output_dir",   metavar="PATH",
                   help="Root output directory (default: ./output)")
    g.add_argument("--scripts_dir",  metavar="PATH",
                   help="Directory containing pipeline scripts (default: .)")
    g.add_argument("--email",        metavar="EMAIL",
                   help="NCBI Entrez email")

    # ── Phase 1 ───────────────────────────────────────────────────────────
    g1 = p.add_argument_group("Phase 1 — genome_to_design.sh")
    g1.add_argument("--organism",      metavar="NAME",
                    help='Organism name or NCBI TaxID, e.g. "Taraxacum officinale"')
    g1.add_argument("--is_eukaryote",  metavar="BOOL",
                    help="true/false (default: false)")
    g1.add_argument("--genemark_path", metavar="PATH",
                    help="Path to GeneMark bin directory")
    g1.add_argument("--bam",           metavar="PATH",
                    help="Path to existing RNA-Seq BAM for BRAKER hints")
    g1.add_argument("--auto_rnaseq",   action="store_true",
                    help="Auto-download RNA-Seq from SRA for BRAKER hints")
    g1.add_argument("--force",         action="store_true",
                    help="Force full re-run of genome_to_design.sh (clears its checkpoints)")

    # ── Phase 2 ───────────────────────────────────────────────────────────
    g2 = p.add_argument_group("Phase 2 — plasmid_design_moclo_v3.py")
    g2.add_argument("--moclo_standard",        metavar="STD",
                    help="marillonnet | cidar | jump (default: marillonnet)")
    g2.add_argument("--enzyme_level0",         metavar="ENZ",
                    help="Level 0 domestication enzyme (default: BsaI-HFv2)")
    g2.add_argument("--enzyme_level1",         metavar="ENZ",
                    help="Level 1 assembly enzyme (default: BpiI)")
    g2.add_argument("--perform_domestication", metavar="BOOL",
                    help="true/false (default: true)")
    g2.add_argument("--output_prefix",         metavar="STR",
                    help="Output file stem (default: moclo_plasmid)")
    g2.add_argument("--genes",                 metavar="PATH", nargs="+",
                    help="One or more CDS FASTA paths; auto-discovered if omitted")

    # ── Phase 3 ───────────────────────────────────────────────────────────
    g3 = p.add_argument_group("Phase 3 — Feedback loops")
    g3.add_argument("--fb6_min_hits",      type=int,   metavar="N",
                    help="FB6: min BLAST hits to trigger partition switch (default: 5)")
    g3.add_argument("--fb6_evalue_cutoff", type=float, metavar="F",
                    help="FB6: max e-value for hits to count (default: 1e-5)")
    g3.add_argument("--dry_run",           action="store_true",
                    help="FB6: report partition mismatch but do NOT re-run BRAKER")
    g3.add_argument("--no_fb6",            action="store_true",
                    help="Skip Feedback Loop 6")

    # ── Orchestrator control ──────────────────────────────────────────────
    gc = p.add_argument_group("Orchestrator control")
    gc.add_argument("--skip_phase1",   action="store_true",
                    help="Skip Phase 1 (assumes a completed run already exists)")
    gc.add_argument("--skip_phase2",   action="store_true",
                    help="Skip Phase 2 (plasmid design)")
    gc.add_argument("--skip_feedback", action="store_true",
                    help="Skip all feedback loops")
    gc.add_argument("--reset",         action="store_true",
                    help="Clear all orchestrator sentinels before running "
                         "(does NOT clear genome_to_design.sh checkpoints; use --force for that)")

    return p


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    global _log_path

    parser = build_parser()
    args   = parser.parse_args()
    cfg    = build_config(args)

    # Set up logging before anything else
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    _log_path = output_dir / "genomopipe_orchestrator.log"

    log("INFO", "=" * 68)
    log("INFO", "Genomopipe orchestrator starting")
    log("INFO", f"  organism    : {cfg['organism'] or '(not set)'}")
    log("INFO", f"  output_dir  : {cfg['output_dir']}")
    log("INFO", f"  scripts_dir : {cfg['scripts_dir']}")
    log("INFO", f"  is_eukaryote: {cfg['is_eukaryote']}")
    log("INFO", f"  auto_rnaseq : {cfg['auto_rnaseq']}")
    log("INFO", f"  skip_phase1 : {cfg['skip_phase1']}")
    log("INFO", f"  skip_phase2 : {cfg['skip_phase2']}")
    log("INFO", f"  skip_feedback:{cfg['skip_feedback']}")
    log("INFO", f"  run_fb6     : {cfg['run_fb6']}")
    log("INFO", f"  fb6_dry_run : {cfg['fb6_dry_run']}")
    log("INFO", "=" * 68)

    # Validate that we can find the scripts before committing to a long run
    phase1_needed  = not cfg["skip_phase1"]
    phase2_needed  = not cfg["skip_phase2"]
    fb6_needed     = cfg["run_fb6"] and not cfg["skip_feedback"]

    missing = []
    if phase1_needed:
        try:
            _find_script(cfg["scripts_dir"], "genome_to_design.sh")
        except FileNotFoundError as e:
            missing.append(str(e))
    if phase2_needed:
        try:
            _find_script(cfg["scripts_dir"], "plasmid_design_moclo_v3.py")
        except FileNotFoundError as e:
            missing.append(str(e))
    if fb6_needed:
        try:
            _find_script(cfg["scripts_dir"], "feedback6_blast_taxonomy_rerun.py")
        except FileNotFoundError as e:
            missing.append(str(e))

    if missing:
        for m in missing:
            log("HALT", m)
        log("HALT", "Set --scripts_dir to the directory containing all pipeline scripts.")
        sys.exit(1)

    # ── Phase 1 ───────────────────────────────────────────────────────────
    run_dir_ref: List[Optional[Path]] = [None]
    run_dir = run_phase1(cfg, run_dir_ref)

    # --reset clears orchestrator sentinels AFTER phase1 resolves run_dir
    if args.reset:
        _clear_sentinels(run_dir)

    # ── Phase 2 ───────────────────────────────────────────────────────────
    run_phase2(cfg, run_dir)

    # ── Phase 3: Feedback loops ───────────────────────────────────────────
    run_feedback_loops(cfg, run_dir)

    # ── Summary ───────────────────────────────────────────────────────────
    write_summary(cfg, run_dir)

    log("INFO", "=" * 68)
    log("INFO", "Genomopipe orchestrator finished.")
    log("INFO", "=" * 68)


if __name__ == "__main__":
    main()

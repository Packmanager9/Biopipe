#!/usr/bin/env python3
# plasmid_design_moclo_v3.py  —  with YAML + JSON support + dynamic path resolution
# Usage: ./plasmid_design_moclo_v3.py /path/to/config.yaml
#                                     or .json or .txt

import os
import sys
import json
from pathlib import Path
import yaml
from typing import Dict, List, Optional, Tuple

from pydna.readers import read as pydna_read
from pydna.dsdna import Dseq, Dseqrecord
from pydna.assembly import Assembly
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.Data import CodonTable

# ─────────────────────────────────────────────
# Optional: superior domestication & optimization
# pip install dnachisel[reports]
# ─────────────────────────────────────────────
try:
    from dnachisel import DnaChiselProblem, AvoidPattern, CodonOptimize
    HAS_DNACHISEL = True
except ImportError:
    HAS_DNACHISEL = False

# ─────────────────────────────────────────────
# Dynamic path resolution
#
# Priority order (highest → lowest):
#   1. Config file keys: output_base, run_dir, parts_dir
#   2. Environment variables:
#        MOCLO_OUTPUT_BASE  – root of all pipeline runs
#        MOCLO_RUN_DIR      – explicit path to a specific run
#        MOCLO_PARTS_DIR    – explicit path to split_seqs
#   3. Auto-discovery: find the most-recently-modified run_* dir
#      under <output_base>/
#   4. Fallback: current working directory
# ─────────────────────────────────────────────

def _resolve_output_base() -> Path:
    """Return the pipeline output root, from env var or default ~/output."""
    env = os.environ.get("MOCLO_OUTPUT_BASE")
    if env:
        return Path(env).expanduser().resolve()
    return Path.home() / "output"


def _resolve_latest_run(output_base: Path) -> Optional[Path]:
    """
    Return the most-recently-modified run dir under output_base.
    Checks the 'latest' symlink first, then falls back to the
    newest run_* directory by mtime.
    """
    # Honour the symlink written by genome_to_design.sh
    latest_link = output_base / "latest"
    if latest_link.is_symlink() and latest_link.exists():
        return latest_link.resolve()

    # No symlink – pick newest run_* folder
    env_run = os.environ.get("MOCLO_RUN_DIR")
    if env_run:
        p = Path(env_run).expanduser().resolve()
        if p.is_dir():
            return p

    runs = sorted(
        output_base.glob("run_*"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return runs[0] if runs else None


def _resolve_parts_dir(run_dir: Optional[Path]) -> Path:
    """
    Locate the ProteinMPNN split_seqs directory.
    Checks env var first, then probes common sub-paths under run_dir,
    then falls back to cwd.
    """
    env = os.environ.get("MOCLO_PARTS_DIR")
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_dir():
            return p

    if run_dir and run_dir.is_dir():
        # Try the canonical sub-path produced by genome_to_design.sh
        candidate = run_dir / "proteinmpnn_out" / "split_seqs"
        if candidate.is_dir():
            return candidate

        # Broader search in case the user changed folder names
        matches = sorted(run_dir.rglob("split_seqs"), key=lambda p: p.stat().st_mtime, reverse=True)
        if matches:
            return matches[0]

        # Fall back to any proteinmpnn output directory
        mpnn_dirs = sorted(run_dir.rglob("proteinmpnn_out"), key=lambda p: p.stat().st_mtime, reverse=True)
        if mpnn_dirs:
            return mpnn_dirs[0]

    return Path.cwd()


# Resolve at import time so load_part() can use module-level PARTS_DIR,
# but allow design_plasmid() to override via config.
OUTPUT_BASE = _resolve_output_base()
LATEST_RUN  = _resolve_latest_run(OUTPUT_BASE)
PARTS_DIR   = _resolve_parts_dir(LATEST_RUN)

# ─────────────────────────────────────────────
# MoClo fusion site presets
# ─────────────────────────────────────────────
MOCLO_OVERHANGS = {
    "marillonnet": {
        "promoter":    ["GGAG", "TACT"],
        "rbs_5utr":    ["AATG", "GCTT"],
        "cds":         ["AATG", "GCTT"],
        "terminator":  ["TACT", "CTGC"],
    },
    "cidar": {
        "promoter":    ["GGAG", "TACT"],
        "rbs":         ["AATG", "GCTT"],
        "cds":         ["AATG", "GCTT"],
        "terminator":  ["TACT", "CTGC"],
    },
    "jump": {
        "promoter":    ["GGAG", "TACT"],
        "rbs":         ["AATG", "GCTT"],
        "cds":         ["AATG", "GCTT"],
        "terminator":  ["TACT", "CTGC"],
    },
}

STANDARD_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]

ECOLI_CODON_FREQ = {  # truncated example — expand or use dnachisel tables
    'ATG': 1.00, 'TGG': 1.00,
    'GCA': 0.59, 'GCC': 0.26, 'GCG': 0.36, 'GCT': 0.46,
    # ... add full 64 codons for best results
}

DEFAULT_ENZYME_PROPS = {
    "BsaI-HFv2": {"recog": "GGTCTC", "cut": "1/5", "overhang_len": 4, "opt_temp": 50, "inact_temp": 65, "buffer": "CutSmart",     "note": "preferred"},
    "BpiI":      {"recog": "GAAGAC", "cut": "2/6", "overhang_len": 4, "opt_temp": 37, "inact_temp": 65, "buffer": "NEBuffer 2.1", "note": "Level 1 common"},
}


# ─────────────────────────────────────────────
# Part loading — respects runtime parts_dir
# ─────────────────────────────────────────────

def load_part(key_or_path: str, parts_dir: Optional[Path] = None) -> Dseqrecord:
    """
    Load a part from an explicit path or by fuzzy-searching parts_dir.

    Search order:
      1. Treat key_or_path as a literal file path.
      2. Search parts_dir (or module-level PARTS_DIR) recursively for
         any FASTA whose name contains key_or_path.
      3. If parts_dir was passed explicitly (from config), also fall back
         to the module-level PARTS_DIR before giving up.
    """
    path = Path(key_or_path).expanduser()
    if path.is_file():
        return pydna_read(str(path))[0]

    search_dirs = []
    if parts_dir:
        search_dirs.append(parts_dir)
    if PARTS_DIR not in search_dirs:
        search_dirs.append(PARTS_DIR)

    for sdir in search_dirs:
        found = sorted(sdir.rglob(f"*{key_or_path}*.f*sta"))
        if found:
            return pydna_read(str(found[0]))[0]

    checked = ", ".join(str(d) for d in search_dirs)
    raise ValueError(f"Part '{key_or_path}' not found as a file or in: {checked}")


# ─────────────────────────────────────────────
# CDS domestication
# ─────────────────────────────────────────────

def domesticate_cds(
    seqrecord: Dseqrecord,
    forbidden_patterns: List[str],
    host_codon_freq: Dict[str, float] = ECOLI_CODON_FREQ,
    max_attempts_per_site: int = 8,
    max_total_iterations: int = 120
) -> Tuple[Dseqrecord, List[str]]:
    if HAS_DNACHISEL:
        try:
            problem = DnaChiselProblem(
                sequence=str(seqrecord.seq),
                constraints=[AvoidPattern(pat) for pat in forbidden_patterns],
                objectives=[CodonOptimize(species="e_coli")]
            )
            problem.resolve()
            new_seq = problem.sequence
            summary = problem.summary()
            new_record = Dseqrecord(new_seq, id=seqrecord.id + "_chisel_opt", name=seqrecord.name)
            return new_record, [f"DNA Chisel: {summary}"]
        except Exception as e:
            print(f"DNA Chisel failed: {e} — falling back")

    # Heuristic fallback
    seq = str(seqrecord.seq).upper()
    original = seq
    changes = []
    table = STANDARD_TABLE.forward_table
    iteration = 0

    while iteration < max_total_iterations:
        found_problem = False
        for pat in forbidden_patterns:
            pos = seq.find(pat)
            if pos == -1:
                continue
            found_problem = True

            start_scan = max(0, pos - 12) // 3 * 3
            end_scan = min(len(seq), pos + len(pat) + 12) // 3 * 3 + 3

            mutated = False
            for c_pos in range(start_scan, end_scan, 3):
                codon = seq[c_pos:c_pos+3]
                aa = table.get(codon)
                if not aa:
                    continue

                synonyms = [c for c in table.forward_table if table.forward_table[c] == aa]
                best_syn = None
                best_score = -1

                for syn in synonyms:
                    if syn == codon:
                        continue
                    test_seq = seq[:c_pos] + syn + seq[c_pos+3:]
                    creates_new = any(
                        test_seq.find(p) != -1 and test_seq.find(p) != pos
                        for p in forbidden_patterns
                    )
                    if creates_new:
                        continue
                    score = host_codon_freq.get(syn, 0.05)
                    if score > best_score:
                        best_score = score
                        best_syn = syn

                if best_syn:
                    seq = seq[:c_pos] + best_syn + seq[c_pos+3:]
                    changes.append(
                        f"Removed {pat} @ {pos} → codon {c_pos}:{codon} → {best_syn} (freq {best_score:.2f})"
                    )
                    mutated = True
                    break

            if not mutated:
                changes.append(f"Could not silently remove {pat} @ {pos}")

        if not found_problem:
            break
        iteration += 1

    if seq == original:
        changes.append("No changes needed or possible")

    new_record = Dseqrecord(seq, id=seqrecord.id + "_dom_opt", name=seqrecord.name)
    return new_record, changes


# ─────────────────────────────────────────────
# Main assembly function
# ─────────────────────────────────────────────

def design_plasmid(config: Dict):
    # ── Resolve paths from config, env, or auto-discovery ──────────────────
    # Allow config to override the output base and run dir
    cfg_output_base = config.get("output_base")
    cfg_run_dir     = config.get("run_dir")
    cfg_parts_dir   = config.get("parts_dir")

    output_base = Path(cfg_output_base).expanduser().resolve() if cfg_output_base else OUTPUT_BASE
    run_dir     = Path(cfg_run_dir).expanduser().resolve()     if cfg_run_dir     else _resolve_latest_run(output_base)
    parts_dir   = Path(cfg_parts_dir).expanduser().resolve()   if cfg_parts_dir   else _resolve_parts_dir(run_dir)

    print(f"Paths resolved:")
    print(f"  output_base : {output_base}")
    print(f"  run_dir     : {run_dir or '(not found)'}")
    print(f"  parts_dir   : {parts_dir}")

    # ── Assembly config ─────────────────────────────────────────────────────
    method    = config.get("assembly_method", "GoldenGate").lower()
    moclo_std = config.get("moclo_standard", None)
    prefix    = config.get("output_prefix", "moclo_plasmid")

    # Output dir: config > <run_dir>/moclo_plasmids > ./moclo_plasmids
    raw_out = config.get("output_dir")
    if raw_out:
        out_dir = Path(raw_out).expanduser().resolve()
    elif run_dir:
        out_dir = run_dir / "moclo_plasmids"
    else:
        out_dir = Path.cwd() / "moclo_plasmids"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Assembly: {method} | MoClo: {moclo_std or 'none'}")
    print(f"Output dir: {out_dir}")

    # ── Enzyme setup ────────────────────────────────────────────────────────
    enzyme_l0_name = config.get("enzyme_level0", "BsaI-HFv2")
    enzyme_l1_name = config.get("enzyme_level1", "BpiI")
    props_l0 = config.get("enzyme_props", {}).get(
        enzyme_l0_name, DEFAULT_ENZYME_PROPS.get(enzyme_l0_name, DEFAULT_ENZYME_PROPS["BsaI-HFv2"])
    )
    props_l1 = config.get("enzyme_props", {}).get(
        enzyme_l1_name, DEFAULT_ENZYME_PROPS.get(enzyme_l1_name, DEFAULT_ENZYME_PROPS["BpiI"])
    )

    if moclo_std and moclo_std in MOCLO_OVERHANGS:
        overhang_table = MOCLO_OVERHANGS[moclo_std]
    elif "overhang_table" in config:
        overhang_table = config["overhang_table"]
    else:
        overhang_table = None

    do_dom  = config.get("perform_domestication", False)
    forbidden = config.get(
        "domestication_enzymes",
        [enzyme_l0_name.split('-')[0], enzyme_l1_name.split('-')[0]],
    )

    # ── Load standard parts ─────────────────────────────────────────────────
    promoter   = load_part(config.get("promoter",   "pJ23119"),    parts_dir)
    rbs        = load_part(config.get("rbs",        "RBS_strong"), parts_dir)
    terminator = load_part(config.get("terminator", "rrnB_T1"),    parts_dir)
    backbone   = load_part(
        config.get("backbone", f"{config.get('ori', 'pMB1_ori')}+{config.get('marker', 'KanR')}"),
        parts_dir,
    )

    # ── Build transcription units ───────────────────────────────────────────
    tus   = []
    genes = config.get("genes", [])

    for i, gene_entry in enumerate(genes):
        gene_key = gene_entry if isinstance(gene_entry, str) else gene_entry.get("path", gene_entry)
        gene = load_part(gene_key, parts_dir)

        dom_changes = []
        if do_dom:
            gene, dom_changes = domesticate_cds(gene, forbidden)
            if dom_changes:
                print(f"Dom/Opt {gene.name}:")
                for ch in dom_changes:
                    print(f"  • {ch}")

        tu_name = f"Level1_TU_{i+1}_{gene.name}"

        # Work on local copies so we don't mutate the module-level objects
        prom_part = promoter
        rbs_part  = rbs
        term_part = terminator

        if overhang_table:
            prom_ov = overhang_table.get("promoter",   ["GGAG", "TACT"])
            rbs_ov  = overhang_table.get("rbs",        ["AATG", "GCTT"])
            cds_ov  = overhang_table.get("cds",        ["AATG", "GCTT"])
            term_ov = overhang_table.get("terminator", ["TACT", "CTGC"])

            prom_part = Dseqrecord(prom_ov[0] + str(prom_part.seq) + prom_ov[1], name=prom_part.name)
            rbs_part  = Dseqrecord(rbs_ov[0]  + str(rbs_part.seq)  + rbs_ov[1],  name=rbs_part.name)
            gene      = Dseqrecord(cds_ov[0]  + str(gene.seq)      + cds_ov[1],  name=gene.name)
            term_part = Dseqrecord(term_ov[0] + str(term_part.seq) + term_ov[1], name=term_part.name)

        tu = prom_part + rbs_part + gene + term_part
        tu.name = tu_name
        tus.append(tu)

    # ── In-silico assembly ──────────────────────────────────────────────────
    ass_limit  = 4 if overhang_table else 35
    ass        = Assembly([backbone] + tus, limit=ass_limit)
    candidates = ass.assemble_circular()
    final      = candidates[0] if candidates else backbone + Dseqrecord("".join(t.seq for t in tus))
    print(f"Assembly candidates: {len(candidates)}")

    # ── Protocol summary ────────────────────────────────────────────────────
    print("\nSuggested protocol:")
    print(f"  Level 0 domestication: remove internal {', '.join(forbidden)}")
    print(f"  Level 1 cycle: {props_l0['opt_temp']}°C 3-5 min → 16°C 4-5 min × 25–50")
    print(f"  Final digest: {props_l0['opt_temp']}°C 60 min")
    print(f"  Inactivate: {props_l0['inact_temp']}°C 20 min")

    # ── Write outputs ───────────────────────────────────────────────────────
    record = final.to_seqrecord(
        id=prefix,
        description=f"MoClo assembly - {moclo_std or 'custom'}",
    )
    pos = 0
    for part in [backbone] + tus:
        qualifiers = {"label": part.name, "note": "MoClo part"}
        feat = SeqFeature(FeatureLocation(pos, pos + len(part)), type="misc_feature", qualifiers=qualifiers)
        record.features.append(feat)
        pos += len(part)

    gb_path    = out_dir / f"{prefix}.gb"
    fasta_path = out_dir / f"{prefix}.fasta"

    SeqIO.write(record, gb_path, "genbank")
    SeqIO.write(record, fasta_path, "fasta")

    print(f"\nGenBank (SnapGene Viewer): {gb_path}")
    print(f"FASTA: {fasta_path}")
    print(f"Length: {len(final)} bp")


# ─────────────────────────────────────────────
# Config loading
# ─────────────────────────────────────────────

def load_config(path: Path) -> Dict:
    suffix = path.suffix.lower()
    try:
        if suffix in (".yaml", ".yml"):
            with open(path) as f:
                return yaml.safe_load(f) or {}
        elif suffix == ".json":
            with open(path) as f:
                return json.load(f)
        else:
            # Plain-text fallback
            config: Dict = {}
            section = None
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if line.endswith(":"):
                        section = line[:-1].strip().lower()
                        continue
                    if ":" in line:
                        k, v = [x.strip() for x in line.split(":", 1)]
                        if section in ("enzyme_props", "overhang_table"):
                            config.setdefault(section, {})[k] = v
                        else:
                            config[k] = v
                    elif section in ("genes", "overhangs", "domestication_enzymes"):
                        config.setdefault(section, []).append(line)
            return config
    except Exception as e:
        print(f"Error loading config file {path}: {e}")
        sys.exit(1)


def main():
    if len(sys.argv) != 2:
        print("Usage: plasmid_design_moclo_v3.py <config.yaml | config.json | instructions.txt>")
        print()
        print("Environment variable overrides (all optional):")
        print("  MOCLO_OUTPUT_BASE  – root directory of pipeline runs  (default: ~/output)")
        print("  MOCLO_RUN_DIR      – explicit path to a specific run")
        print("  MOCLO_PARTS_DIR    – explicit path to split_seqs")
        print()
        print("Config-file keys (take precedence over env vars):")
        print("  output_base, run_dir, parts_dir, output_dir")
        sys.exit(1)

    config = load_config(Path(sys.argv[1]))
    design_plasmid(config)


if __name__ == "__main__":
    main()
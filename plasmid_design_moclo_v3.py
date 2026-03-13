#!/usr/bin/env python3
# plasmid_design_moclo_v3.py  —  MoClo Golden Gate plasmid design
#                                 with codon optimization + CDS domestication
#
# Usage: ./plasmid_design_moclo_v3.py /path/to/config.yaml
#                                     or .json or .txt
#
# New in this version:
#   • Codon optimization step (before domestication)
#   • Automatic protein→DNA back-translation using host codon tables
#   • Supports E. coli, S. cerevisiae, H. sapiens/CHO, P. pastoris, B. subtilis
#   • Uses python-codon-tables if installed; falls back to built-in tables
#   • Uses DNA Chisel CodonOptimize if available (most powerful)
#   • New config keys: codon_optimize, expression_host, codon_optimize_method

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
# Optional: python-codon-tables
# pip install python-codon-tables
# ─────────────────────────────────────────────
try:
    from python_codon_tables import get_codon_table as _pct_get_table
    HAS_PYTHON_CODON_TABLES = True
except ImportError:
    HAS_PYTHON_CODON_TABLES = False

# ─────────────────────────────────────────────
# Dynamic path resolution
# ─────────────────────────────────────────────

def _resolve_output_base() -> Path:
    env = os.environ.get("MOCLO_OUTPUT_BASE")
    if env:
        return Path(env).expanduser().resolve()
    return Path.home() / "output"


def _resolve_latest_run(output_base: Path) -> Optional[Path]:
    latest_link = output_base / "latest"
    if latest_link.is_symlink() and latest_link.exists():
        return latest_link.resolve()
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
    env = os.environ.get("MOCLO_PARTS_DIR")
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_dir():
            return p
    if run_dir and run_dir.is_dir():
        candidate = run_dir / "proteinmpnn_out" / "split_seqs"
        if candidate.is_dir():
            return candidate
        matches = sorted(run_dir.rglob("split_seqs"), key=lambda p: p.stat().st_mtime, reverse=True)
        if matches:
            return matches[0]
        mpnn_dirs = sorted(run_dir.rglob("proteinmpnn_out"), key=lambda p: p.stat().st_mtime, reverse=True)
        if mpnn_dirs:
            return mpnn_dirs[0]
    return Path.cwd()


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

# ─────────────────────────────────────────────
# Built-in codon frequency tables
# {codon: (amino_acid, relative_frequency_0_to_1)}
# Frequencies are normalised so the most-used codon per AA = 1.0
# Sources: Kazusa codon usage database, NAR 2000 (RSCU-normalised)
# ─────────────────────────────────────────────

def _build_codon_table(raw: Dict[str, float]) -> Dict[str, Tuple[str, float]]:
    """
    Attach amino-acid labels and normalise per-AA so max codon = 1.0.
    raw: {codon: absolute_frequency_per_thousand}
    """
    fwd = STANDARD_TABLE.forward_table
    stop_codons = set(STANDARD_TABLE.stop_codons)

    # Group codons by amino acid
    aa_groups: Dict[str, Dict[str, float]] = {}
    for codon, freq in raw.items():
        if codon in stop_codons:
            aa = "*"
        elif codon in fwd:
            aa = fwd[codon]
        else:
            continue
        aa_groups.setdefault(aa, {})[codon] = freq

    # Normalise within each amino-acid group
    result: Dict[str, Tuple[str, float]] = {}
    for aa, codons in aa_groups.items():
        max_freq = max(codons.values()) or 1.0
        for codon, freq in codons.items():
            result[codon] = (aa, round(freq / max_freq, 4))
    return result


# ── E. coli K-12 ─────────────────────────────────────────────────────────────
# Frequencies per thousand from Kazusa (E. coli K-12)
_ECOLI_RAW: Dict[str, float] = {
    "TTT": 22.2, "TTC": 16.6,                                          # Phe
    "TTA":  4.0, "TTG": 13.1, "CTT": 11.1, "CTC": 11.0,              # Leu
    "CTA":  3.5, "CTG": 52.8,
    "ATT": 30.3, "ATC": 25.5, "ATA":  4.4,                            # Ile
    "ATG": 27.4,                                                        # Met
    "GTT": 18.1, "GTC": 15.3, "GTA": 11.0, "GTG": 26.0,              # Val
    "TCT":  8.7, "TCC": 8.7,  "TCA":  7.0, "TCG": 8.9,               # Ser
    "AGT":  8.8, "AGC": 16.1,
    "CCT":  7.3, "CCC":  5.5, "CCA":  8.5, "CCG": 23.2,              # Pro
    "ACT":  9.0, "ACC": 23.4, "ACA":  7.1, "ACG": 14.5,              # Thr
    "GCT": 15.3, "GCC": 26.1, "GCA": 21.2, "GCG": 33.6,              # Ala
    "TAT": 16.1, "TAC": 12.2,                                          # Tyr
    "CAT": 13.0, "CAC": 10.3,                                          # His
    "CAA": 15.0, "CAG": 34.2,                                          # Gln
    "AAT": 17.8, "AAC": 22.4,                                          # Asn
    "AAA": 36.5, "AAG": 10.8,                                          # Lys
    "GAT": 37.5, "GAC": 21.2,                                          # Asp
    "GAA": 39.4, "GAG": 18.3,                                          # Glu
    "TGT":  5.2, "TGC":  6.8,                                          # Cys
    "TGG": 15.3,                                                        # Trp
    "CGT": 21.0, "CGC": 22.6, "CGA":  3.6, "CGG": 5.4,               # Arg
    "AGA":  2.1, "AGG":  1.2,
    "GGT": 24.2, "GGC": 29.9, "GGA":  8.0, "GGG": 11.4,              # Gly
    "TAA":  2.0, "TAG":  0.3, "TGA":  1.1,                            # Stop
}

# ── Saccharomyces cerevisiae ──────────────────────────────────────────────────
_YEAST_RAW: Dict[str, float] = {
    "TTT": 18.0, "TTC": 21.1,
    "TTA": 11.7, "TTG": 27.2, "CTT": 12.3, "CTC":  5.4,
    "CTA":  6.9, "CTG": 11.0,
    "ATT": 30.1, "ATC": 17.2, "ATA": 11.4,
    "ATG": 20.9,
    "GTT": 18.8, "GTC": 11.8, "GTA":  7.4, "GTG": 10.8,
    "TCT": 23.5, "TCC": 14.2, "TCA": 18.7, "TCG":  8.6,
    "AGT": 14.2, "AGC":  9.8,
    "CCT": 13.5, "CCC":  6.4, "CCA": 41.8, "CCG":  5.3,
    "ACT": 20.3, "ACC": 13.0, "ACA": 17.8, "ACG":  8.0,
    "GCT": 20.9, "GCC": 12.7, "GCA": 16.4, "GCG":  6.2,
    "TAT": 18.7, "TAC": 15.0,
    "CAT": 13.6, "CAC":  7.7,
    "CAA": 30.1, "CAG": 12.1,
    "AAT": 36.6, "AAC": 25.0,
    "AAA": 41.9, "AAG": 30.8,
    "GAT": 37.6, "GAC": 20.2,
    "GAA": 45.5, "GAG": 19.2,
    "TGT":  8.1, "TGC":  4.8,
    "TGG": 10.4,
    "CGT":  6.4, "CGC":  2.6, "CGA":  2.7, "CGG":  1.7,
    "AGA": 21.3, "AGG": 9.2,
    "GGT": 24.0, "GGC":  9.8, "GGA": 10.9, "GGG":  6.0,
    "TAA":  1.1, "TAG":  0.5, "TGA":  0.6,
}

# ── Homo sapiens (also representative for CHO) ───────────────────────────────
_HUMAN_RAW: Dict[str, float] = {
    "TTT": 17.6, "TTC": 20.3,
    "TTA":  7.7, "TTG": 12.9, "CTT": 13.2, "CTC": 19.6,
    "CTA":  7.2, "CTG": 39.6,
    "ATT": 15.9, "ATC": 20.8, "ATA":  7.5,
    "ATG": 22.0,
    "GTT": 11.0, "GTC": 14.5, "GTA":  7.1, "GTG": 28.1,
    "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG":  4.4,
    "AGT": 15.1, "AGC": 19.5,
    "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG":  6.9,
    "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG":  6.1,
    "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG":  7.4,
    "TAT": 12.2, "TAC": 15.3,
    "CAT": 10.9, "CAC": 15.1,
    "CAA": 12.3, "CAG": 34.2,
    "AAT": 17.0, "AAC": 19.1,
    "AAA": 24.4, "AAG": 31.9,
    "GAT": 21.8, "GAC": 25.1,
    "GAA": 29.0, "GAG": 39.6,
    "TGT": 10.6, "TGC": 12.6,
    "TGG": 13.2,
    "CGT":  4.5, "CGC": 10.4, "CGA":  6.2, "CGG": 11.4,
    "AGA": 12.2, "AGG": 12.0,
    "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
    "TAA":  1.0, "TAG":  0.5, "TGA":  1.6,
}

# ── Pichia pastoris (Komagataella phaffii) ────────────────────────────────────
_PICHIA_RAW: Dict[str, float] = {
    "TTT": 24.4, "TTC": 19.5,
    "TTA": 16.0, "TTG": 28.3, "CTT": 14.1, "CTC":  7.4,
    "CTA":  8.3, "CTG": 10.1,
    "ATT": 30.6, "ATC": 20.4, "ATA":  8.0,
    "ATG": 21.3,
    "GTT": 21.9, "GTC": 13.5, "GTA": 10.3, "GTG": 13.5,
    "TCT": 21.5, "TCC": 13.1, "TCA": 17.2, "TCG":  7.1,
    "AGT": 16.3, "AGC": 11.7,
    "CCT": 14.8, "CCC":  7.2, "CCA": 39.1, "CCG":  4.6,
    "ACT": 22.1, "ACC": 14.0, "ACA": 18.0, "ACG":  7.3,
    "GCT": 22.5, "GCC": 13.8, "GCA": 17.8, "GCG":  6.8,
    "TAT": 21.5, "TAC": 16.2,
    "CAT": 15.1, "CAC":  8.4,
    "CAA": 31.8, "CAG": 13.5,
    "AAT": 38.2, "AAC": 26.1,
    "AAA": 43.2, "AAG": 32.1,
    "GAT": 40.1, "GAC": 22.8,
    "GAA": 47.6, "GAG": 20.3,
    "TGT":  9.2, "TGC":  5.6,
    "TGG": 10.9,
    "CGT":  7.0, "CGC":  3.0, "CGA":  3.1, "CGG":  2.1,
    "AGA": 22.8, "AGG": 10.5,
    "GGT": 25.0, "GGC": 10.8, "GGA": 12.1, "GGG":  6.4,
    "TAA":  1.3, "TAG":  0.6, "TGA":  0.7,
}

# ── Bacillus subtilis 168 ─────────────────────────────────────────────────────
_BSUBTILIS_RAW: Dict[str, float] = {
    "TTT": 25.3, "TTC": 13.9,
    "TTA": 11.4, "TTG": 14.0, "CTT": 13.2, "CTC":  7.4,
    "CTA":  7.2, "CTG": 10.4,
    "ATT": 34.2, "ATC": 21.3, "ATA":  7.0,
    "ATG": 26.4,
    "GTT": 26.2, "GTC": 11.2, "GTA": 13.7, "GTG": 22.5,
    "TCT": 15.0, "TCC":  8.0, "TCA": 11.4, "TCG":  6.5,
    "AGT": 16.0, "AGC": 12.9,
    "CCT": 12.0, "CCC":  6.1, "CCA": 11.9, "CCG": 14.2,
    "ACT": 18.2, "ACC": 12.9, "ACA": 10.6, "ACG": 13.9,
    "GCT": 26.0, "GCC": 13.5, "GCA": 19.8, "GCG": 20.1,
    "TAT": 21.4, "TAC": 12.4,
    "CAT": 16.2, "CAC":  8.0,
    "CAA": 21.7, "CAG": 16.9,
    "AAT": 31.2, "AAC": 21.8,
    "AAA": 42.8, "AAG": 14.0,
    "GAT": 42.5, "GAC": 18.3,
    "GAA": 45.0, "GAG": 16.2,
    "TGT":  7.2, "TGC":  5.8,
    "TGG": 14.4,
    "CGT": 20.0, "CGC": 11.4, "CGA":  4.0, "CGG":  4.3,
    "AGA":  5.1, "AGG":  2.6,
    "GGT": 29.0, "GGC": 16.5, "GGA": 14.0, "GGG":  8.8,
    "TAA":  2.4, "TAG":  0.7, "TGA":  1.1,
}

# Compiled tables (built on first access via a module-level dict)
_BUILTIN_TABLES: Dict[str, Dict[str, Tuple[str, float]]] = {
    "e_coli":      _build_codon_table(_ECOLI_RAW),
    "s_cerevisiae": _build_codon_table(_YEAST_RAW),
    "h_sapiens":   _build_codon_table(_HUMAN_RAW),
    "p_pastoris":  _build_codon_table(_PICHIA_RAW),
    "b_subtilis":  _build_codon_table(_BSUBTILIS_RAW),
}

# Aliases: map common user-supplied names to internal keys
_HOST_ALIASES: Dict[str, str] = {
    # E. coli variants
    "ecoli": "e_coli", "e.coli": "e_coli", "escherichia coli": "e_coli",
    "e_coli k12": "e_coli", "e_coli b": "e_coli", "bl21": "e_coli",
    # Yeast S. cerevisiae
    "yeast": "s_cerevisiae", "saccharomyces cerevisiae": "s_cerevisiae",
    "s.cerevisiae": "s_cerevisiae", "s_cerevisiae": "s_cerevisiae",
    # Human / CHO
    "human": "h_sapiens", "homo sapiens": "h_sapiens", "h.sapiens": "h_sapiens",
    "cho": "h_sapiens", "mammalian": "h_sapiens", "hek293": "h_sapiens",
    # Pichia / Komagataella
    "pichia": "p_pastoris", "pichia pastoris": "p_pastoris",
    "komagataella phaffii": "p_pastoris", "p.pastoris": "p_pastoris",
    # B. subtilis
    "bacillus subtilis": "b_subtilis", "b.subtilis": "b_subtilis",
    "b_subtilis 168": "b_subtilis",
}

# dnachisel species strings (used when HAS_DNACHISEL=True)
_DNACHISEL_SPECIES: Dict[str, str] = {
    "e_coli":       "e_coli",
    "s_cerevisiae": "s_cerevisiae",
    "h_sapiens":    "h_sapiens",
    "p_pastoris":   "p_pastoris",
    "b_subtilis":   "b_subtilis",
}

# python-codon-tables organism strings
_PCT_SPECIES: Dict[str, str] = {
    "e_coli":       "Escherichia coli",
    "s_cerevisiae": "Saccharomyces cerevisiae",
    "h_sapiens":    "Homo sapiens",
    "p_pastoris":   "Komagataella phaffii",
    "b_subtilis":   "Bacillus subtilis",
}


def _resolve_host(host_str: str) -> str:
    """Normalise an expression host string to an internal key."""
    key = host_str.lower().strip()
    if key in _BUILTIN_TABLES:
        return key
    if key in _HOST_ALIASES:
        return _HOST_ALIASES[key]
    # Partial match
    for alias, canonical in _HOST_ALIASES.items():
        if alias in key or key in alias:
            return canonical
    print(f"[WARN] Unknown expression host '{host_str}'. "
          f"Falling back to 'e_coli'. "
          f"Supported: {', '.join(sorted(_BUILTIN_TABLES.keys()))}")
    return "e_coli"


def _get_codon_table_for_host(host_key: str) -> Dict[str, Tuple[str, float]]:
    """
    Return {codon: (aa, norm_freq)} for host_key.
    Prefers python-codon-tables if available; falls back to built-in tables.
    """
    if HAS_PYTHON_CODON_TABLES:
        pct_name = _PCT_SPECIES.get(host_key, host_key)
        try:
            raw_pct = _pct_get_table(pct_name)  # {aa: {codon: freq}}
            flat: Dict[str, float] = {}
            for aa_codons in raw_pct.values():
                flat.update(aa_codons)
            return _build_codon_table(flat)
        except Exception as exc:
            print(f"[WARN] python-codon-tables failed for '{pct_name}': {exc}. "
                  "Using built-in table.")
    return _BUILTIN_TABLES.get(host_key, _BUILTIN_TABLES["e_coli"])


# ─────────────────────────────────────────────────────────────────────────────
# Sequence-type detection
# ─────────────────────────────────────────────────────────────────────────────

_PROTEIN_ONLY_CHARS = set("FHILMPQRWYfhilmpqrwy")


def _is_protein_sequence(seq: str) -> bool:
    """Return True if seq contains amino-acid characters absent from DNA."""
    return bool(_PROTEIN_ONLY_CHARS.intersection(seq))


# ─────────────────────────────────────────────────────────────────────────────
# Codon optimization
# ─────────────────────────────────────────────────────────────────────────────

def back_translate(aa_seq: str, codon_table: Dict[str, Tuple[str, float]]) -> str:
    """
    Convert an amino-acid sequence to a DNA CDS using the most-frequent codon
    for each amino acid in the given host codon table.

    Uses max-frequency codon selection (CAI maximisation). The ATG start codon
    is forced at position 0 if the protein sequence starts with M.
    """
    # Build aa → best_codon lookup
    aa_to_best: Dict[str, str] = {}
    aa_to_all: Dict[str, List[Tuple[float, str]]] = {}
    for codon, (aa, freq) in codon_table.items():
        if aa == "*":
            continue
        aa_to_all.setdefault(aa, []).append((freq, codon))
    for aa, options in aa_to_all.items():
        aa_to_best[aa] = max(options)[1]  # highest freq

    # Standard table as fallback for any AA not in codon_table
    fwd = STANDARD_TABLE.forward_table
    back_table: Dict[str, List[str]] = {}
    for codon, aa in fwd.items():
        back_table.setdefault(aa, []).append(codon)

    dna_codons: List[str] = []
    for i, aa in enumerate(aa_seq.upper().replace("*", "")):
        if aa == "M" and i == 0:
            dna_codons.append("ATG")
        elif aa in aa_to_best:
            dna_codons.append(aa_to_best[aa])
        elif aa in back_table:
            dna_codons.append(back_table[aa][0])
        else:
            # Unknown AA — encode as NNN placeholder
            print(f"[WARN] Unknown amino acid '{aa}' at position {i}; encoding as NNN.")
            dna_codons.append("NNN")

    return "".join(dna_codons)


def codon_optimize_dna(
    dna_seq: str,
    codon_table: Dict[str, Tuple[str, float]],
) -> str:
    """
    Re-encode an existing CDS using the most-frequent synonymous codons from
    codon_table. Preserves the reading frame (assumes seq starts at codon 0).
    Ignores partial trailing codons.
    """
    fwd = STANDARD_TABLE.forward_table
    stop_codons = set(STANDARD_TABLE.stop_codons)

    # Build aa → best_codon
    aa_to_best: Dict[str, str] = {}
    aa_to_all: Dict[str, List[Tuple[float, str]]] = {}
    for codon, (aa, freq) in codon_table.items():
        if aa == "*":
            continue
        aa_to_all.setdefault(aa, []).append((freq, codon))
    for aa, options in aa_to_all.items():
        aa_to_best[aa] = max(options)[1]

    dna_out: List[str] = []
    seq = dna_seq.upper()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in stop_codons:
            dna_out.append(codon)
            break
        aa = fwd.get(codon)
        if aa and aa in aa_to_best:
            dna_out.append(aa_to_best[aa])
        else:
            dna_out.append(codon)  # keep original if unresolvable
    return "".join(dna_out)


def optimize_and_back_translate(
    seqrecord: Dseqrecord,
    host_key: str,
    method: str = "auto",
) -> Tuple[Dseqrecord, List[str]]:
    """
    Main entry point for codon optimization.

    method: "auto"         — use DNA Chisel if available, then python-codon-tables, then built-in
            "dnachisel"    — require DNA Chisel (error if not installed)
            "max_frequency"— use max-frequency lookup from host table (no dnachisel)
            "harmonize"    — placeholder; falls back to max_frequency

    Handles both protein and DNA input:
      - Protein input: always back-translates to optimized DNA
      - DNA input:     re-encodes using synonymous substitutions

    Returns (optimized_Dseqrecord, [change_log_strings])
    """
    seq_str = str(seqrecord.seq).strip()
    is_protein = _is_protein_sequence(seq_str)
    changes: List[str] = []
    codon_table = _get_codon_table_for_host(host_key)
    pct_name = _PCT_SPECIES.get(host_key, host_key)
    dc_species = _DNACHISEL_SPECIES.get(host_key, "e_coli")

    # ── DNA Chisel path ───────────────────────────────────────────────────────
    if (method in ("auto", "dnachisel")) and HAS_DNACHISEL:
        try:
            if is_protein:
                # Back-translate first using max-freq, then chisel will optimise further
                dna_seq = back_translate(seq_str, codon_table)
                changes.append(f"Back-translated protein → DNA ({len(dna_seq)} bp) "
                                f"using {host_key} max-frequency table")
            else:
                dna_seq = seq_str

            problem = DnaChiselProblem(
                sequence=dna_seq,
                constraints=[],
                objectives=[CodonOptimize(species=dc_species)],
            )
            problem.resolve()
            optimized_seq = problem.sequence
            changes.append(
                f"DNA Chisel CodonOptimize ({dc_species}): "
                f"{problem.objectives_text_summary()}"
            )
            new_rec = Dseqrecord(
                optimized_seq,
                id=seqrecord.id + "_codon_opt",
                name=seqrecord.name,
            )
            return new_rec, changes
        except Exception as exc:
            print(f"[WARN] DNA Chisel codon optimization failed: {exc}. Falling back.")

    if method == "dnachisel":
        raise RuntimeError("dnachisel method requested but DNA Chisel is not available "
                           "or failed. pip install dnachisel[reports]")

    # ── Max-frequency / python-codon-tables path ──────────────────────────────
    if is_protein:
        dna_seq = back_translate(seq_str, codon_table)
        src = "python-codon-tables" if HAS_PYTHON_CODON_TABLES else "built-in table"
        changes.append(
            f"Back-translated protein ({len(seq_str)} aa) → DNA ({len(dna_seq)} bp) "
            f"using {host_key} max-frequency table ({src})"
        )
    else:
        original_dna = seq_str
        dna_seq = codon_optimize_dna(seq_str, codon_table)
        n_changes = sum(
            1 for i in range(0, min(len(original_dna), len(dna_seq)) - 2, 3)
            if original_dna[i:i+3] != dna_seq[i:i+3]
        )
        src = "python-codon-tables" if HAS_PYTHON_CODON_TABLES else "built-in table"
        changes.append(
            f"Re-encoded {n_changes} codons for {host_key} max-frequency ({src})"
        )

    new_rec = Dseqrecord(
        dna_seq,
        id=seqrecord.id + "_codon_opt",
        name=seqrecord.name,
    )
    return new_rec, changes


# ─────────────────────────────────────────────
# Part loading
# ─────────────────────────────────────────────

def load_part(key_or_path: str, parts_dir: Optional[Path] = None) -> Dseqrecord:
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
# CDS domestication (restriction-site removal)
# ─────────────────────────────────────────────

ECOLI_CODON_FREQ: Dict[str, float] = {
    c: freq for c, (aa, freq) in _BUILTIN_TABLES["e_coli"].items()
}


def domesticate_cds(
    seqrecord: Dseqrecord,
    forbidden_patterns: List[str],
    host_codon_freq: Dict[str, float] = ECOLI_CODON_FREQ,
    max_attempts_per_site: int = 8,
    max_total_iterations: int = 120,
) -> Tuple[Dseqrecord, List[str]]:
    """
    Remove internal restriction enzyme recognition sites via silent (synonymous)
    substitutions. Respects host codon frequency when choosing replacements.

    Uses DNA Chisel if available; otherwise falls back to a heuristic scan.
    """
    if HAS_DNACHISEL:
        try:
            problem = DnaChiselProblem(
                sequence=str(seqrecord.seq),
                constraints=[AvoidPattern(pat) for pat in forbidden_patterns],
                objectives=[],  # codon optimization already done upstream
            )
            problem.resolve()
            new_seq = problem.sequence
            summary = problem.summary()
            new_record = Dseqrecord(
                new_seq,
                id=seqrecord.id + "_dom",
                name=seqrecord.name,
            )
            return new_record, [f"DNA Chisel AvoidPattern: {summary}"]
        except Exception as e:
            print(f"[WARN] DNA Chisel domestication failed: {e} — falling back to heuristic")

    # Heuristic fallback
    seq = str(seqrecord.seq).upper()
    original = seq
    changes: List[str] = []
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

                synonyms = [c for c in table if table[c] == aa]
                best_syn, best_score = None, -1

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
                        f"Removed {pat} @ {pos}: codon {c_pos} {codon}→{best_syn} "
                        f"(freq {best_score:.2f})"
                    )
                    mutated = True
                    break

            if not mutated:
                changes.append(f"Could not silently remove {pat} @ {pos}")

        if not found_problem:
            break
        iteration += 1

    if seq == original:
        changes.append("No restriction-site changes needed or possible")

    new_record = Dseqrecord(seq, id=seqrecord.id + "_dom", name=seqrecord.name)
    return new_record, changes


# ─────────────────────────────────────────────
# Main assembly function
# ─────────────────────────────────────────────

def design_plasmid(config: Dict):
    # ── Resolve paths ──────────────────────────────────────────────────────
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

    # ── Assembly config ────────────────────────────────────────────────────
    method    = config.get("assembly_method", "GoldenGate").lower()
    moclo_std = config.get("moclo_standard", None)
    prefix    = config.get("output_prefix", "moclo_plasmid")

    raw_out = config.get("output_dir")
    if raw_out:
        out_dir = Path(raw_out).expanduser().resolve()
    elif run_dir:
        out_dir = run_dir / "moclo_plasmids"
    else:
        out_dir = Path.cwd() / "moclo_plasmids"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Codon optimization config ──────────────────────────────────────────
    do_codon_opt   = config.get("codon_optimize", False)
    host_str       = config.get("expression_host", "e_coli")
    opt_method     = config.get("codon_optimize_method", "auto")
    host_key       = _resolve_host(host_str) if do_codon_opt else "e_coli"

    print(f"\nAssembly : {method} | MoClo: {moclo_std or 'none'}")
    print(f"Output   : {out_dir}")
    if do_codon_opt:
        print(f"Codon opt: ENABLED — host={host_key}, method={opt_method}")
        print(f"  DNA Chisel    : {'available' if HAS_DNACHISEL else 'NOT installed'}")
        print(f"  python-codon-tables: {'available' if HAS_PYTHON_CODON_TABLES else 'NOT installed'}")
    else:
        print("Codon opt: disabled (set codon_optimize: true to enable)")

    # ── Enzyme setup ───────────────────────────────────────────────────────
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

    do_dom   = config.get("perform_domestication", False)
    forbidden = config.get(
        "domestication_enzymes",
        [enzyme_l0_name.split("-")[0], enzyme_l1_name.split("-")[0]],
    )

    # ── Build host codon frequency dict for domestication heuristic ────────
    host_codon_freq: Dict[str, float] = {
        c: freq for c, (aa, freq) in _get_codon_table_for_host(host_key).items()
    }

    # ── Load standard parts ────────────────────────────────────────────────
    promoter   = load_part(config.get("promoter",   "pJ23119"),    parts_dir)
    rbs        = load_part(config.get("rbs",        "RBS_strong"), parts_dir)
    terminator = load_part(config.get("terminator", "rrnB_T1"),    parts_dir)
    backbone   = load_part(
        config.get("backbone", f"{config.get('ori', 'pMB1_ori')}+{config.get('marker', 'KanR')}"),
        parts_dir,
    )

    # ── Build transcription units ──────────────────────────────────────────
    tus   = []
    genes = config.get("genes", [])

    for i, gene_entry in enumerate(genes):
        gene_key = gene_entry if isinstance(gene_entry, str) else gene_entry.get("path", gene_entry)
        gene = load_part(gene_key, parts_dir)

        gene_name = gene.name or f"gene_{i+1}"
        print(f"\n── Gene {i+1}: {gene_name} ──")

        # ── Step A: Codon optimization / back-translation ──────────────────
        if do_codon_opt:
            gene, opt_changes = optimize_and_back_translate(gene, host_key, method=opt_method)
            print(f"  Codon optimization:")
            for ch in opt_changes:
                print(f"    • {ch}")
        else:
            # If input is protein sequence but no codon_optimize flag, we still
            # need to back-translate to DNA — use E. coli max-frequency as a
            # neutral default and warn the user.
            if _is_protein_sequence(str(gene.seq)):
                print(f"  [WARN] '{gene_name}' appears to be a protein sequence "
                      "but codon_optimize is False.")
                print("         Back-translating with E. coli defaults. "
                      "Set codon_optimize: true and expression_host to optimize properly.")
                gene, opt_changes = optimize_and_back_translate(
                    gene, "e_coli", method="max_frequency"
                )

        # ── Step B: CDS domestication (restriction site removal) ───────────
        dom_changes: List[str] = []
        if do_dom:
            gene, dom_changes = domesticate_cds(gene, forbidden, host_codon_freq)
            if dom_changes:
                print(f"  Domestication:")
                for ch in dom_changes:
                    print(f"    • {ch}")

        # ── Step C: Assemble transcription unit ────────────────────────────
        tu_name = f"Level1_TU_{i+1}_{gene_name}"

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

    # ── In-silico assembly ─────────────────────────────────────────────────
    ass_limit  = 4 if overhang_table else 35
    ass        = Assembly([backbone] + tus, limit=ass_limit)
    candidates = ass.assemble_circular()
    final      = candidates[0] if candidates else backbone + Dseqrecord("".join(t.seq for t in tus))
    print(f"\nAssembly candidates: {len(candidates)}")

    # ── Protocol summary ───────────────────────────────────────────────────
    print("\nSuggested protocol:")
    print(f"  Codon optimization  : {'Yes — ' + host_key if do_codon_opt else 'No'}")
    print(f"  Level 0 domestication: remove internal {', '.join(forbidden)}")
    print(f"  Level 1 cycle: {props_l0['opt_temp']}°C 3-5 min → 16°C 4-5 min × 25–50")
    print(f"  Final digest: {props_l0['opt_temp']}°C 60 min")
    print(f"  Inactivate: {props_l0['inact_temp']}°C 20 min")

    # ── Write outputs ──────────────────────────────────────────────────────
    record = final.to_seqrecord(
        id=prefix,
        description=f"MoClo assembly - {moclo_std or 'custom'} - "
                    f"{'codon_opt:' + host_key if do_codon_opt else 'no_codon_opt'}",
    )
    pos = 0
    for part in [backbone] + tus:
        qualifiers = {"label": part.name, "note": "MoClo part"}
        feat = SeqFeature(
            FeatureLocation(pos, pos + len(part)),
            type="misc_feature",
            qualifiers=qualifiers,
        )
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
# Enzyme defaults
# ─────────────────────────────────────────────

DEFAULT_ENZYME_PROPS = {
    "BsaI-HFv2": {
        "recog": "GGTCTC", "cut": "1/5", "overhang_len": 4,
        "opt_temp": 50, "inact_temp": 65, "buffer": "CutSmart",
        "note": "preferred",
    },
    "BpiI": {
        "recog": "GAAGAC", "cut": "2/6", "overhang_len": 4,
        "opt_temp": 37, "inact_temp": 65, "buffer": "NEBuffer 2.1",
        "note": "Level 1 common",
    },
}


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
        print("Usage: plasmid_design_moclo_v3.py <config.yaml | config.json | config.txt>")
        print()
        print("Codon optimization keys (new):")
        print("  codon_optimize: true/false      — enable codon optimization (default: false)")
        print("  expression_host: <host>         — target expression host (default: e_coli)")
        print("  codon_optimize_method: <method> — auto | max_frequency | dnachisel (default: auto)")
        print()
        print("Supported expression hosts:")
        for k in sorted(_BUILTIN_TABLES.keys()):
            print(f"  {k}")
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

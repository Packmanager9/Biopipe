#!/usr/bin/env python3
# =============================================================================
# feedback6_blast_taxonomy_rerun.py
#
# FEEDBACK LOOP 6: BLAST results → taxonomy verification → BRAKER re-run
#                  with corrected OrthoDB partition
#
# WORKFLOW:
#   1. Parse blast_results.txt from Step 7
#   2. For each unique subject accession, look up its taxonomic lineage via
#      NCBI Entrez and assign it to the SINGLE most-specific OrthoDB partition
#      keyword (e.g. a vertebrate hit counts as "vertebrata", not also
#      "metazoa" or "eukaryota")
#   3. Compare the dominant inferred lineage against the partition that was
#      actually used in Step 3 (read from the BRAKER step log)
#   4. If >= min_hits accessions map to a different, more specific partition,
#      download it and re-run BRAKER with the corrected OrthoDB hints
#   5. Run gffread (+ GTF repair) and protein extraction on the new annotation
#   6. Write a taxonomy audit report: original partition, suggested partition,
#      per-partition hit counts, top BLAST subjects, and whether a re-run
#      was triggered
#
# WHY THIS HELPS:
#   The initial OrthoDB group is determined by the organism's own Entrez
#   taxonomy.  If the designed proteins (Steps 5-5b) turn out to be
#   homologous to sequences from a DIFFERENT lineage, adding that lineage's
#   OrthoDB partition as protein hints may recover gene models that were
#   missed or fragmented in the original BRAKER run.
#
# USAGE:
#   python feedback6_blast_taxonomy_rerun.py config.yaml
#   python feedback6_blast_taxonomy_rerun.py config.json
#   python feedback6_blast_taxonomy_rerun.py config.txt
#   python feedback6_blast_taxonomy_rerun.py \
#       --run_dir ./output/latest \
#       --min_hits 5 \
#       --evalue_cutoff 1e-5 \
#       --dry_run
#
# CONFIG KEYS (yaml/json/txt, all optional):
#   run_dir          Pipeline run directory (default: ./output/latest)
#   min_hits         Minimum number of BLAST hits from a given partition to
#                    trigger a switch (default: 5)
#   evalue_cutoff    Maximum e-value for a hit to be counted (default: 1e-5)
#   dry_run          Report only; do NOT download or re-run BRAKER (default: false)
#   output_dir       Where to write outputs (default: <run_dir>/feedback6_loop)
#   genemark_path    Path to GeneMark bin dir
#   is_eukaryote     true/false (default: true)
#   email            NCBI Entrez email (required for taxonomy lookups)
# =============================================================================

import os
import re
import csv
import sys
import json
import time
import shutil
import argparse
import subprocess
import textwrap
from collections import Counter
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

try:
    from Bio import Entrez, SeqIO
    HAS_BIO = True
except ImportError:
    HAS_BIO = False

# ─────────────────────────────────────────────
# OrthoDB partition map
#
# Maps lineage keywords (as they appear lowercased in Entrez taxonomy strings)
# to the OrthoDB download slug at bioinf.uni-greifswald.de.
#
# SPECIFICITY_ORDER defines most→least specific for per-record assignment:
# each BLAST-hit record is counted exactly once, under its most specific
# matching keyword.  This prevents "eukaryota" from being inflated by every
# vertebrate or plant hit.
# ─────────────────────────────────────────────
ORTHODB_MAP: Dict[str, str] = {
    "mammalia":      "Vertebrata",   # no separate mammal partition
    "vertebrata":    "Vertebrata",
    "arthropoda":    "Arthropoda",
    "fungi":         "Fungi",
    "viridiplantae": "Viridiplantae",
    "alveolata":     "Alveolata",
    "stramenopiles": "Stramenopiles",
    "amoebozoa":     "Amoebozoa",
    "euglenozoa":    "Euglenozoa",
    "metazoa":       "Metazoa",
    "eukaryota":     "Eukaryota",
}

# Ordered from most to least specific — first match wins for per-record counting
SPECIFICITY_ORDER: List[str] = [
    "mammalia", "vertebrata", "arthropoda", "fungi", "viridiplantae",
    "alveolata", "stramenopiles", "amoebozoa", "euglenozoa",
    "metazoa", "eukaryota",
]

ORTHODB_BASE = "https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12"

# ─────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────
_log_file: Optional[str] = None

def log(level: str, msg: str) -> None:
    ts   = datetime.now().strftime("%F %T")
    line = f"[{ts}] [{level}] {msg}"
    print(line)
    if _log_file:
        with open(_log_file, "a") as f:
            f.write(line + "\n")

# ─────────────────────────────────────────────
# Config helpers
# ─────────────────────────────────────────────
def _str_to_bool(val) -> bool:
    """Convert bool, int, or string 'true'/'false'/'1'/'0' to bool.
    Raises ValueError on unrecognised strings."""
    if isinstance(val, bool):
        return val
    if isinstance(val, int):
        return bool(val)
    s = str(val).strip().lower()
    if s in ("true", "1", "yes"):
        return True
    if s in ("false", "0", "no"):
        return False
    raise ValueError(f"Cannot interpret '{val}' as a boolean")


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
            cfg: dict = {}
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    k, v = line.split(":", 1)
                    cfg[k.strip()] = v.strip()
            return cfg


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Feedback Loop 6: BLAST taxonomy → BRAKER partition correction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("config", nargs="?", help="YAML / JSON / txt config file")
    p.add_argument("--run_dir",       default=None)
    p.add_argument("--min_hits",      type=int,   default=None)
    p.add_argument("--evalue_cutoff", type=float, default=None)
    p.add_argument("--dry_run",       action="store_true", default=False)
    p.add_argument("--output_dir",    default=None)
    p.add_argument("--genemark_path", default=None)
    p.add_argument("--is_eukaryote",  default=None)
    p.add_argument("--email",         default=None)
    return p.parse_args()


def resolve_config(args: argparse.Namespace) -> dict:
    """Merge defaults → config file → CLI flags (highest priority)."""
    cfg: dict = {
        "run_dir":       str(Path("./output/latest").resolve()),
        "min_hits":      5,
        "evalue_cutoff": 1e-5,
        "dry_run":       False,
        "output_dir":    None,
        "genemark_path": str(Path.home() / "genemark-etp-full/gmetp_linux_64/bin"),
        "is_eukaryote":  "true",
        "email":         "PricklyPearEnterprises@gmail.com",
    }

    if args.config and Path(args.config).is_file():
        for k, v in load_config_file(args.config).items():
            if k in cfg:
                cfg[k] = v

    # CLI flags override config file (only when explicitly provided)
    for k in ("run_dir", "min_hits", "evalue_cutoff", "output_dir",
              "genemark_path", "is_eukaryote", "email"):
        val = getattr(args, k, None)
        if val is not None:
            cfg[k] = val
    if args.dry_run:                  # argparse store_true → always safe bool
        cfg["dry_run"] = True

    # Normalise types now so downstream code never has to guess
    cfg["min_hits"]     = int(cfg["min_hits"])
    cfg["evalue_cutoff"] = float(cfg["evalue_cutoff"])
    cfg["dry_run"]      = _str_to_bool(cfg["dry_run"])   # FIX: was bool("false") == True
    cfg["is_eukaryote"] = _str_to_bool(cfg["is_eukaryote"])

    return cfg

# ─────────────────────────────────────────────
# BLAST result parsing
# ─────────────────────────────────────────────
def parse_blast_subjects(blast_file: Path,
                         evalue_cutoff: float) -> List[str]:
    """Return a de-duplicated list of subject accessions passing evalue filter."""
    seen: set = set()
    accs: List[str] = []
    with open(blast_file) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 11:
                continue
            try:
                ev = float(row[10])
            except ValueError:
                continue
            if ev <= evalue_cutoff:
                acc = row[1].split(".")[0]   # strip version number
                if acc not in seen:
                    seen.add(acc)
                    accs.append(acc)
    return accs

# ─────────────────────────────────────────────
# Entrez taxonomy lookup
# ─────────────────────────────────────────────
def fetch_lineages(accs: List[str],
                   email: str,
                   batch_size: int = 50) -> Counter:
    """
    For each subject accession, determine its single most-specific OrthoDB
    partition keyword and count how many accessions fall into each partition.

    Each accession is counted EXACTLY ONCE — under its most specific
    matching keyword in SPECIFICITY_ORDER.  This prevents generic terms
    like "eukaryota" from being inflated by every vertebrate hit.

    Returns Counter of {partition_keyword: accession_count}.
    """
    if not HAS_BIO:
        raise ImportError("Biopython required: pip install biopython")

    Entrez.email = email
    partition_counts: Counter = Counter()
    total = len(accs)

    for i in range(0, total, batch_size):
        batch = accs[i : i + batch_size]
        log("INFO", f"Fetching taxonomy batch {i // batch_size + 1} / "
                    f"{(total + batch_size - 1) // batch_size} "
                    f"({len(batch)} accessions)")
        try:
            handle  = Entrez.efetch(db="protein", id=",".join(batch),
                                    rettype="gb", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            for rec in records:
                # GBSeq_taxonomy is a semicolon-separated lineage string,
                # e.g. "Eukaryota; Viridiplantae; Streptophyta; ..."
                lineage_str = rec.get("GBSeq_taxonomy", "").lower()

                # FIX: assign each record to its SINGLE most-specific keyword
                for kw in SPECIFICITY_ORDER:
                    if kw in lineage_str:
                        partition_counts[kw] += 1
                        break   # stop at first (most specific) match

        except Exception as e:
            log("WARN", f"Entrez batch failed: {e}")

        time.sleep(0.4)   # respect NCBI rate limit (3 req/s without API key)

    return partition_counts

# ─────────────────────────────────────────────
# Detect original OrthoDB partition from step log
# ─────────────────────────────────────────────
def detect_original_partition(run_dir: Path) -> Optional[str]:
    """
    Determine which OrthoDB partition was used in Step 3 by searching
    the step3_annotate log for the orthodb filename.

    Returns a SPECIFICITY_ORDER keyword (e.g. 'viridiplantae') or None.
    """
    # The per-step log written by run_step() is the most reliable source
    candidate_logs = [
        run_dir / "logs" / "step3_annotate.log",
        run_dir / "braker_out" / "braker.log",
    ]

    for log_path in candidate_logs:
        if not log_path.exists():
            continue
        try:
            text = log_path.read_text(errors="replace").lower()
            # Look for the orthodb filename pattern: orthodb_<group>.fasta
            m = re.search(r"orthodb_([a-z]+)(?:_corrected)?\.fasta", text)
            if m:
                kw = m.group(1)
                # Map slug back to keyword (e.g. "vertebrata" → "vertebrata")
                if kw in ORTHODB_MAP:
                    return kw
                # Handle aliased slugs: slug "vertebrata" matches keyword "vertebrata"
                # and "mammalia" both map to the same slug — return the slug keyword
                for k, slug in ORTHODB_MAP.items():
                    if slug.lower() == kw:
                        return k
        except Exception as e:
            log("WARN", f"Could not read {log_path}: {e}")

    return None

# ─────────────────────────────────────────────
# OrthoDB download + header sanitization
# ─────────────────────────────────────────────
def download_orthodb(group_slug: str, dest_path: Path) -> bool:
    """Download, decompress, and sanitize an OrthoDB partition FASTA."""
    if not HAS_BIO:
        raise ImportError("Biopython required: pip install biopython")

    url    = f"{ORTHODB_BASE}/{group_slug}.fa.gz"
    gz     = dest_path.parent / (dest_path.stem + ".fa.gz")

    log("INFO", f"Downloading {url}")
    rc = subprocess.run(["wget", "-q", url, "-O", str(gz)]).returncode
    if rc != 0 or not gz.exists():
        log("WARN", f"wget failed for {url}")
        return False

    rc = subprocess.run(["gunzip", "-f", str(gz)]).returncode
    raw = gz.with_suffix("")                    # .fa.gz → .fa
    if raw.exists():
        shutil.move(str(raw), str(dest_path))

    if not dest_path.exists():
        log("WARN", f"Decompressed file not found at {dest_path}")
        return False

    # Sanitize headers (spaces and pipes → underscores, clear description)
    try:
        records = []
        for rec in SeqIO.parse(str(dest_path), "fasta"):
            rec.id          = rec.id.replace(" ", "_").replace("|", "_")
            rec.description = ""
            records.append(rec)
        SeqIO.write(records, str(dest_path), "fasta")
        log("INFO", f"Sanitized {len(records)} OrthoDB headers → {dest_path}")
    except Exception as e:
        log("WARN", f"Header sanitization failed: {e}")

    return dest_path.exists() and dest_path.stat().st_size > 0

# ─────────────────────────────────────────────
# BRAKER re-run
# ─────────────────────────────────────────────
def run_braker(genome: Path,
               prot_seq: Path,
               out_dir: Path,
               anno_prefix: str,
               genemark_path: str,
               bam: Optional[str],
               num_threads: int,
               step_log: Path) -> bool:
    """
    Invoke BRAKER inside braker_env via `conda run`, with the environment
    variables that BRAKER requires explicitly forwarded.

    FIX: the original version passed no AUGUSTUS_* env vars, so BRAKER
    always failed to locate Augustus binaries and configs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve conda prefix for braker_env so we can build the required paths
    result = subprocess.run(
        ["conda", "run", "-n", "braker_env", "bash", "-c", "echo $CONDA_PREFIX"],
        capture_output=True, text=True
    )
    conda_prefix = result.stdout.strip() or ""
    aug_config   = str(Path(conda_prefix) / "config") if conda_prefix else ""
    aug_bin      = str(Path(conda_prefix) / "bin")    if conda_prefix else ""

    # Write-check the Augustus config dir; copy if not writable (same logic
    # as genome_to_design.sh step 3)
    aug_config_writable = aug_config
    if aug_config and not os.access(aug_config, os.W_OK):
        writable = str(Path.home() / "augustus_config_writable")
        if not Path(writable).exists():
            shutil.copytree(aug_config, writable)
        aug_config_writable = writable
        log("INFO", f"Augustus config not writable — using {writable}")

    # Remove stale species dir so BRAKER trains fresh
    species_dir = Path(aug_config_writable) / "species" / anno_prefix
    if species_dir.exists():
        shutil.rmtree(str(species_dir))
        log("INFO", f"Removed stale Augustus species dir: {species_dir}")

    env_extras = {
        "AUGUSTUS_CONFIG_PATH":   aug_config_writable,
        "AUGUSTUS_BIN_PATH":      aug_bin,
        "AUGUSTUS_SCRIPTS_PATH":  aug_bin,
        "GENEMARK_PATH":          genemark_path,
        "PROTHINT_PATH":          str(Path(genemark_path) / "ProtHint" / "bin"),
        "BAMTOOLS_PATH":          str(Path(aug_bin) / "bamtools") if aug_bin else "",
    }

    braker_cmd = [
        "conda", "run", "--no-capture-output", "-n", "braker_env",
        "braker.pl",
        f"--AUGUSTUS_CONFIG_PATH={aug_config_writable}",
        f"--AUGUSTUS_BIN_PATH={aug_bin}",
        f"--AUGUSTUS_SCRIPTS_PATH={aug_bin}",
        f"--GENEMARK_PATH={genemark_path}",
        f"--PROTHINT_PATH={genemark_path}/ProtHint/bin",
        f"--genome={genome}",
        f"--prot_seq={prot_seq}",
        f"--species={anno_prefix}",
        f"--workingdir={out_dir}",
        f"--threads={num_threads}",
        "--softmasking",
        "--nocleanup",
        "--augustus_args=--min_intron_len=10 --max_intron_len=25000",
    ]
    if bam and Path(bam).exists():
        braker_cmd.append(f"--bam={bam}")

    log("INFO", f"BRAKER: genome={genome.name}  prot_seq={prot_seq.name}")
    log("INFO", f"BRAKER: workingdir={out_dir}")

    merged_env = {**os.environ, **env_extras}
    with open(step_log, "w") as lf:
        rc = subprocess.run(braker_cmd,
                            stdout=lf,
                            stderr=subprocess.STDOUT,
                            env=merged_env).returncode

    if rc != 0:
        log("WARN", f"BRAKER exited {rc} — check {step_log}")
        # Surface the last 30 lines for quick diagnosis
        try:
            lines = step_log.read_text().splitlines()[-30:]
            for l in lines:
                log("WARN", f"  {l}")
        except Exception:
            pass

    return rc == 0

# ─────────────────────────────────────────────
# GTF repair (mirrors genome_to_design.sh step 3)
# ─────────────────────────────────────────────
def repair_gtf(gtf_in: Path, gtf_out: Path) -> bool:
    """
    Inject missing transcript_id fields.  GeneMark sometimes emits records
    without transcript_id; gffread requires it.  Mirrors the repair logic
    in genome_to_design.sh.
    """
    ok = repaired = 0
    try:
        with open(gtf_in) as fin, open(gtf_out, "w") as fout:
            for line in fin:
                if line.startswith("#") or not line.strip():
                    fout.write(line)
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    fout.write(line)
                    continue
                attrs = parts[8]
                if "transcript_id" in attrs:
                    ok += 1
                    fout.write(line)
                    continue
                m = re.search(r'gene_id\s+"?([^";]+)"?', attrs)
                if m:
                    tid = m.group(1).strip() + ".t1"
                    parts[8] = attrs.rstrip("; ") + f'; transcript_id "{tid}";'
                    repaired += 1
                fout.write("\t".join(parts) + "\n")
        log("INFO", f"GTF repair: {ok} OK, {repaired} transcript_id fields injected → {gtf_out}")
        return True
    except Exception as e:
        log("WARN", f"GTF repair failed: {e}")
        return False

# ─────────────────────────────────────────────
# gffread protein extraction (mirrors genome_to_design.sh step 3)
# ─────────────────────────────────────────────
def run_gffread(gtf: Path, genome: Path, out_faa: Path,
                conda_prefix: str) -> bool:
    """
    FIX: the original script assumed BRAKER produces a .faa directly.
    It does not — gffread must convert GTF + genome.fa → protein FASTA.
    """
    gffread = Path(conda_prefix) / "bin" / "gffread" if conda_prefix \
              else Path(shutil.which("gffread") or "gffread")

    if not gffread.exists():
        log("WARN", f"gffread not found at {gffread}. "
                    "Install with: mamba install -n braker_env -c bioconda gffread")
        return False

    log("INFO", f"Running gffread: {gtf.name} + {genome.name} → {out_faa.name}")
    rc = subprocess.run(
        [str(gffread), str(gtf), "-g", str(genome), "-y", str(out_faa)],
        capture_output=True, text=True,
    )
    if rc.returncode != 0:
        log("WARN", f"gffread failed: {rc.stderr.strip()}")
    return rc.returncode == 0 and out_faa.exists()

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main() -> None:
    global _log_file

    if not HAS_BIO:
        print("ERROR: Biopython required — pip install biopython")
        sys.exit(1)

    args = parse_args()
    cfg  = resolve_config(args)

    run_dir    = Path(cfg["run_dir"]).expanduser().resolve()
    min_hits   = cfg["min_hits"]
    evalue_cut = cfg["evalue_cutoff"]
    dry_run    = cfg["dry_run"]
    genemark   = cfg["genemark_path"]
    is_euk     = cfg["is_eukaryote"]
    email      = cfg["email"]
    output_dir = (Path(cfg["output_dir"]).expanduser().resolve()
                  if cfg["output_dir"] else run_dir / "feedback6_loop")

    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = output_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    _log_file = str(log_dir / "feedback6.log")

    log("INFO", "=" * 60)
    log("INFO", "Feedback Loop 6: BLAST taxonomy → BRAKER partition correction")
    log("INFO", f"run_dir        : {run_dir}")
    log("INFO", f"min_hits       : {min_hits}")
    log("INFO", f"evalue_cutoff  : {evalue_cut}")
    log("INFO", f"dry_run        : {dry_run}")
    log("INFO", f"output_dir     : {output_dir}")
    log("INFO", "=" * 60)

    blast_file = run_dir / "blast_results.txt"
    if not blast_file.exists():
        log("HALT", f"blast_results.txt not found in {run_dir}")
        sys.exit(1)

    num_threads = os.cpu_count() or 4

    # ── Step 1: Parse BLAST subjects ─────────────────────────────────────────
    log("INFO", "Parsing BLAST results...")
    accs = parse_blast_subjects(blast_file, evalue_cut)
    if not accs:
        log("HALT", "No BLAST hits passed the evalue filter")
        sys.exit(1)
    log("INFO", f"Unique subject accessions (e-val ≤ {evalue_cut}): {len(accs)}")

    # ── Step 2: Per-accession most-specific lineage assignment ────────────────
    log("INFO", "Querying NCBI taxonomy (cap: 200 accessions)...")
    partition_counts = fetch_lineages(accs[:200], email)

    if not partition_counts:
        log("WARN", "No lineage data returned — check network connection and email")

    log("INFO", "Partition counts (one count per unique accession):")
    for kw, count in partition_counts.most_common():
        log("INFO", f"  {kw:<20} {count:>5}")

    # ── Step 3: Detect original partition ────────────────────────────────────
    original_partition = detect_original_partition(run_dir)
    log("INFO", f"Original BRAKER partition: {original_partition or 'unknown'}")

    # ── Step 4: Determine suggested partition ────────────────────────────────
    # Pick the most frequent partition that:
    #   (a) has >= min_hits accessions
    #   (b) maps to a DIFFERENT slug than the original partition
    #       (handles mammalia vs vertebrata aliasing)
    original_slug = ORTHODB_MAP.get(original_partition, "") if original_partition else ""

    suggested: Optional[str] = None
    for kw, count in partition_counts.most_common():
        if count < min_hits:
            break
        if ORTHODB_MAP.get(kw, "") != original_slug:
            suggested = kw
            break

    # ── Step 5: Write audit report ────────────────────────────────────────────
    report_lines = [
        "Taxonomy Audit Report",
        "=" * 60,
        f"Original OrthoDB partition : {original_partition or 'unknown'}"
        + (f" (slug: {original_slug})" if original_slug else ""),
        f"Suggested partition        : {suggested or 'no change needed'}"
        + (f" (slug: {ORTHODB_MAP[suggested]})" if suggested else ""),
        "",
        "Partition hit counts (one count per unique BLAST subject accession):",
    ]
    for kw, count in partition_counts.most_common():
        marker = " ← original" if kw == original_partition else \
                 " ← SUGGESTED" if kw == suggested else ""
        report_lines.append(f"  {kw:<20} {count:>5}{marker}")
    report_lines.append("=" * 60)

    if suggested is None:
        report_lines.append(
            "CONCLUSION: BLAST hits are consistent with the original partition. "
            "No re-run triggered."
        )
        log("INFO", "No partition switch needed.")
    else:
        report_lines.append(
            f"CONCLUSION: '{suggested}' accounts for the majority of high-confidence "
            f"BLAST hits. Switching from '{original_partition or 'unknown'}' "
            f"→ '{suggested}' (slug: {ORTHODB_MAP[suggested]})."
        )
        if dry_run:
            report_lines.append("DRY RUN: BRAKER re-run skipped.")
            log("INFO", f"DRY RUN: would switch to '{suggested}' partition.")

    report_path = output_dir / "feedback6_taxonomy_audit.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    log("INFO", f"Audit report: {report_path}")

    if suggested is None or dry_run:
        log("INFO", "Done.")
        sys.exit(0)

    # ── Step 6: Download corrected OrthoDB partition ──────────────────────────
    slug        = ORTHODB_MAP[suggested]
    new_orthodb = Path.home() / f"orthodb_{suggested}_corrected.fasta"

    if new_orthodb.exists() and new_orthodb.stat().st_size > 1000:
        log("INFO", f"Reusing cached {new_orthodb}")
    else:
        log("INFO", f"Downloading OrthoDB {slug} partition...")
        if not download_orthodb(slug, new_orthodb):
            log("HALT", f"Failed to download OrthoDB {slug} partition")
            sys.exit(1)

    # ── Step 7: Resolve genome ───────────────────────────────────────────────
    genome: Optional[Path] = None
    for g_name in ("masked.fna", "clean.fna"):
        g = run_dir / g_name
        if g.exists():
            genome = g
            break
    if genome is None:
        log("HALT", "No masked.fna or clean.fna found in run_dir")
        sys.exit(1)

    bam_path: Optional[str] = None
    merged_bam = run_dir / "rnaseq_merged.bam"
    if merged_bam.exists():
        bam_path = str(merged_bam)
        log("INFO", f"Carrying forward RNA-Seq BAM: {bam_path}")

    # ── Step 8: BRAKER re-run ────────────────────────────────────────────────
    anno_prefix = f"corrected_{suggested}"
    braker_out  = output_dir / "braker_corrected"
    braker_log  = log_dir / "braker_corrected.log"

    log("INFO", f"Re-running BRAKER with '{suggested}' OrthoDB partition...")
    braker_ok = run_braker(genome, new_orthodb, braker_out, anno_prefix,
                           genemark, bam_path, num_threads, braker_log)
    if not braker_ok:
        log("WARN", "BRAKER re-run reported errors — attempting protein extraction anyway")

    # ── Step 9: GTF repair + gffread → protein FASTA ─────────────────────────
    # FIX: BRAKER does not write a .faa — must run gffread on the output GTF.
    # Also repair missing transcript_id fields as in genome_to_design.sh.
    gtf_raw   = braker_out / "braker.gtf"
    gtf_fixed = braker_out / "braker_fixed.gtf"

    # Prefer genome.fa written by BRAKER (--nocleanup); fall back to input genome
    braker_genome = braker_out / "genome.fa"
    if not braker_genome.exists():
        log("WARN", "BRAKER genome.fa not found — falling back to input genome")
        braker_genome = genome

    src_faa   = braker_out / f"{anno_prefix}.faa"   # gffread output location
    out_faa   = output_dir / "proteins_corrected.faa"

    gffread_ok = False
    if gtf_raw.exists():
        if repair_gtf(gtf_raw, gtf_fixed):
            # Resolve conda prefix once for gffread path
            res = subprocess.run(
                ["conda", "run", "-n", "braker_env", "bash", "-c", "echo $CONDA_PREFIX"],
                capture_output=True, text=True,
            )
            conda_prefix = res.stdout.strip()
            gffread_ok = run_gffread(gtf_fixed, braker_genome, src_faa, conda_prefix)
    else:
        log("WARN", f"braker.gtf not found in {braker_out} — skipping gffread")

    # ── Step 10: Filter proteins ──────────────────────────────────────────────
    if gffread_ok and src_faa.exists():
        log("INFO", "Filtering proteins (longest isoform, min-length cutoff)...")
        all_records = list(SeqIO.parse(str(src_faa), "fasta"))
        by_gene: Dict[str, object] = {}
        for rec in all_records:
            m      = re.match(r"^(.+)\.t\d+$", rec.id)
            gid    = m.group(1) if m else rec.id
            if gid not in by_gene or len(rec.seq) > len(by_gene[gid].seq):  # type: ignore[union-attr]
                by_gene[gid] = rec

        min_len  = 200 if is_euk else 100
        proteins = sorted(
            [r for r in by_gene.values() if len(r.seq) > min_len],  # type: ignore[union-attr]
            key=lambda r: r.id,
        )
        SeqIO.write(proteins, str(out_faa), "fasta")
        log("INFO",
            f"Input isoforms : {len(all_records)}\n"
            f"Unique genes   : {len(by_gene)}\n"
            f"Written        : {len(proteins)} proteins (>{min_len} aa) → {out_faa}")
    else:
        log("WARN", "Protein extraction skipped — gffread did not produce a .faa")

    # ── Final summary ─────────────────────────────────────────────────────────
    log("INFO", "=" * 60)
    log("INFO", "Feedback Loop 6 complete.")
    log("INFO", f"Original partition  : {original_partition or 'unknown'}")
    log("INFO", f"Corrected partition : {suggested}")
    log("INFO", f"New OrthoDB hints   : {new_orthodb}")
    log("INFO", f"Corrected BRAKER    : {braker_out}")
    log("INFO", f"Corrected proteins  : {out_faa}")
    log("INFO", f"Audit report        : {report_path}")
    log("INFO", "")
    log("INFO", "To continue the design pipeline from the corrected proteins:")
    log("INFO", f"  cp {out_faa} {run_dir}/proteins.faa")
    log("INFO", "  rm <run_dir>/.step5_rfdiffusion.done   # and later checkpoints")
    log("INFO", "  Then re-run genome_to_design.sh (it will resume from Step 5)")
    log("INFO", "=" * 60)


if __name__ == "__main__":
    main()

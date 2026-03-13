# Genomopipe

An end-to-end automated pipeline that takes an organism name and produces computationally designed protein structures and lab-ready plasmid assembly files - from raw genome fetching through annotation, backbone design, sequence design, structure prediction, MoClo Golden Gate cloning design, and a suite of iterative feedback loops that feed design results back into earlier stages to improve quality. A native desktop GUI (Genomopipe App) provides a full graphical interface for configuring, launching, monitoring, and exploring results.

---

## Table of Contents

1. [Overview](#overview)
2. [Scripts](#scripts)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [The Desktop GUI (Genomopipe App)](#the-desktop-gui-genomopipe-app)
6. [Usage](#usage)
7. [The Orchestrator (genomopipe.py)](#the-orchestrator-genomopipepy)
8. [Running Scripts Individually](#running-scripts-individually)
9. [Pipeline Steps In Detail](#pipeline-steps-in-detail)
10. [Codon Optimization](#codon-optimization)
11. [Feedback Loops](#feedback-loops)
12. [Resumability and Checkpoints](#resumability-and-checkpoints)
13. [Config File Reference](#config-file-reference)
14. [Path Resolution](#path-resolution)
15. [Output Structure](#output-structure)
16. [Tips and Notes](#tips-and-notes)
17. [Troubleshooting](#troubleshooting)
18. [Citations](#citations)

---

## Overview

```
                        ┌─────────────────────────┐
                        │     genomopipe.py        │
                        │   (master orchestrator)  │
                        └────────────┬────────────┘
                                     │
              ┌──────────────────────┼──────────────────────┐
              │                      │                      │
           Phase 1                Phase 2               Phase 3
    genome_to_design.sh   plasmid_design_moclo_v3.py   Feedback loops
              │
    ┌─────────▼──────────────────────────────────────────────────┐
    │  [Step 1]  genome_fetch.sh    → NCBI genome download        │
    │  [Step 1b] Sanitize FASTA     → Normalize headers           │
    │  [Step 2]  BBTools QC         → Adapter trim / QC           │
    │  [Step 2b] RepeatMasker       → Mask repeats (eukaryotes)   │
    │  [Step 3]  BRAKER / Prokka    → Gene annotation             │
    │  [Step 4]  Protein Extract    → Filter longest isoforms     │
    │  [Step 5]  RFdiffusion        → Novel backbone design       │
    │  [Step 5b] ProteinMPNN        → Sequence design             │
    │  [Step 6]  ColabFold          → Structure prediction        │
    │  [Step 7]  BLAST              → Validate + assign function  │
    │  [Step 8]  Rename             → Name outputs by function    │
    └─────────────────────────────────────────────────────────────┘
              │
    ┌─────────▼──────────────────────────────────────────────────┐
    │  [Step 9]  MoClo Plasmid Design                             │
    │            plasmid_design_moclo_v3.py                       │
    │            Domesticate CDSs → add fusion sites → assemble   │
    │            Outputs: .gb (SnapGene) + .fasta (synthesis)     │
    └─────────────────────────────────────────────────────────────┘
              │
    ┌─────────▼──────────────────────────────────────────────────┐
    │  Feedback Loops (Phase 3)                                   │
    │                                                             │
    │  [FB1]  ColabFold → RFdiffusion motif re-scaffolding       │
    │  [FB2]  pLDDT filter → ProteinMPNN resampling              │
    │  [FB3]  BLAST hits → BRAKER re-annotation (enriched hints) │
    │  [FB4]  Domesticated CDS → ColabFold re-validation         │
    │  [FB5]  Validated designs → annotation hints               │
    │  [FB6]  BLAST taxonomy → BRAKER with corrected OrthoDB     │
    └─────────────────────────────────────────────────────────────┘
```

---

## Scripts

| Script | Language | Role |
|---|---|---|
| `genomopipe.py` | Python | **Master orchestrator** - single entry point for the full pipeline |
| `genome_to_design.sh` | Bash | Phase 1 - genome fetch through structure prediction (Steps 1–8) |
| `genome_fetch.sh` | Bash | Called by genome_to_design.sh; NCBI genome download with fallback chain |
| `plasmid_design_moclo_v3.py` | Python | Phase 2 (Step 9) - MoClo Golden Gate plasmid design |
| `feedback1_colabfold_to_rfdiffusion.sh` | Bash | FB1 - exploit best ColabFold hits as RFdiffusion motifs |
| `feedback2_plddt_mpnn_resample.py` | Python | FB2 - resample ProteinMPNN for low-pLDDT designs |
| `feedback3_blast_to_braker.sh` | Bash | FB3 - enrich BRAKER hints with BLAST hit sequences |
| `feedback4_domesticated_cds_revalidate.py` | Python | FB4 - re-validate domesticated CDSs with ColabFold |
| `feedback5_designed_proteins_to_annotation.sh` | Bash | FB5 - use validated designs as annotation hints |
| `feedback6_blast_taxonomy_rerun.py` | Python | FB6 - correct OrthoDB partition based on BLAST taxonomy |

FB6 is integrated into the orchestrator. FB1–FB5 are run as standalone scripts; they share the same config format and run directory conventions.

---

## Requirements

### Hardware
- Multi-core CPU (pipeline auto-detects and uses all available threads via `nproc`)
- GPU strongly recommended for ColabFold and RFdiffusion (falls back to CPU)
- 50–200 GB disk space depending on organism genome size
- 32+ GB RAM recommended for large eukaryotic genomes (RepeatModeler on genomes >500 Mbp can peak at 40 GB)

### Python dependencies

Required for the orchestrator and all Python feedback scripts:

```bash
pip install pyyaml biopython pydna
# Optional - substantially improves CDS domestication quality in Step 9:
pip install dnachisel[reports]
# Optional - provides authoritative NCBI codon usage tables for all organisms:
pip install python-codon-tables
```

All Python scripts require Python 3.7+ and run outside any Conda environment; they invoke Conda environments internally via `conda run`.

### Conda environments

| Environment | Used for |
|---|---|
| `braker_env` | Genome annotation, QC, protein extraction, BLAST, gffread, utility steps |
| `SE3nv` | RFdiffusion backbone design and ProteinMPNN sequence design |
| `colabfold` | ColabFold structure prediction |

#### `braker_env`
```bash
conda create -n braker_env
conda activate braker_env
mamba install -c bioconda -c conda-forge \
    braker augustus repeatmodeler repeatmasker bbduk prokka \
    blast gffread biopython star samtools sra-tools \
    python perl
```

#### `SE3nv`
```bash
conda env create -f ~/RFdiffusion/env/SE3nv.yml
conda activate SE3nv
pip install ~/RFdiffusion/env/SE3Transformer
pip install -e ~/RFdiffusion --no-deps
```

#### `colabfold`
```bash
conda create -n colabfold
conda activate colabfold
pip install colabfold[alphafold]
```

#### `bioenv` (used by `genome_fetch.sh`)
```bash
conda activate bioenv
mamba install -c bioconda ncbi-genome-download entrez-direct cd-hit biopython
```

### External tools

| Tool | Default path | Notes |
|---|---|---|
| GeneMark-ETP | `~/genemark-etp-full/gmetp_linux_64/bin` | License required from Georgia Tech |
| RFdiffusion | `~/RFdiffusion/` | Clone from GitHub |
| RFdiffusion models | `~/models/Base_ckpt.pt` | Download `Base_ckpt.pt` separately |
| ProteinMPNN | `~/ProteinMPNN/` | Clone from GitHub |

---

## Installation

```bash
# 1. Place all scripts in the same directory
chmod +x genome_to_design.sh \
         feedback1_colabfold_to_rfdiffusion.sh \
         feedback3_blast_to_braker.sh \
         feedback5_designed_proteins_to_annotation.sh \
         genome_fetch.sh \
         genomopipe.py \
         plasmid_design_moclo_v3.py \
         feedback2_plddt_mpnn_resample.py \
         feedback4_domesticated_cds_revalidate.py \
         feedback6_blast_taxonomy_rerun.py

# 2. Install Miniconda if not already present
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/etc/profile.d/conda.sh

# 3. Set up Conda environments (see Requirements above)

# 4. Install Python dependencies
pip install pyyaml biopython pydna

# 5. Download GeneMark license key - place .gm_key in $HOME
# https://genemark.bme.gatech.edu/license_download.cgi

# 6. Download RFdiffusion model weights
mkdir -p ~/models
# See: https://github.com/RosettaCommons/RFdiffusion

# 7. Install and launch the desktop app
cd bioforge_app
npm install
npm start
```

---

## The Desktop GUI (Genomopipe App)

The Genomopipe desktop app is an Electron application providing a full graphical interface for configuring, launching, monitoring, and exploring pipeline results. It is designed to be the primary interface for working with the pipeline.

### Launching

```bash
cd bioforge_app
npm start
```

### Interface layout

The app has a vertical icon sidebar on the left for tab navigation, a persistent run-directory file tree in a collapsible left panel, and a main content area on the right. The title bar shows current pipeline status (IDLE / RUNNING / ERROR) and a scrolling ticker at the bottom displays the most recent log line. Settings persist to `~/.config/Genomopipe/bioforge-settings.json`.

### Tabs

#### ⚙ Configure

The main launch panel. Set organism name, output directory, and optional config file, then choose annotation options (eukaryote mode, auto RNA-Seq download, BAM file) and pipeline control flags (force re-run). Three action buttons:

- **RUN PIPELINE** - launches the pipeline in the background; switches to the Run tab automatically
- **LOAD RESULTS** - loads an existing completed or in-progress run from the output directory
- **BROWSE RESULTS** - opens an OS folder picker to load results from any path

#### ▶ Run

Live pipeline execution view. Shows a progress bar, per-step status pills (colour-coded: active cyan / done green / failed red), and a full scrolling log with syntax highlighting (`[START]` cyan, `[DONE]` green, `[FAIL]` red, `[WARN]` amber, `[SKIP]` grey). Auto-scroll can be toggled. A kill button sends SIGTERM to the running process.

The sidebar auto-refreshes every 4 seconds while running, showing new checkpoint sentinel files as they appear. The app also detects pipelines launched externally from the terminal - hit ↻ on the sidebar or load the output directory and it begins live-refreshing automatically. While the pipeline is actively running, `.failed` sentinel files are shown in grey rather than red to indicate they may be stale leftovers from a previous `--force` run.

#### 🔬 Structures

Loads all PDB files from `designs/` (RFdiffusion) and `colabfold_out/` (ColabFold) into a list. Clicking a structure loads it into the embedded 3DMol viewer.

**Viewer controls:**

| Input | Action |
|---|---|
| W / S | Rotate X axis |
| A / D | Rotate Y axis |
| Q / E | Rotate Z axis |
| Scroll | Zoom |
| Drag | Pan |
| Shift + Drag | Rotate |
| R | Reset view |
| F | Toggle spin |
| X | Export PNG |
| Escape or ← Back | Return to previous tab |

Style (Cartoon / Stick / Sphere / Surface) and colour scheme (Spectrum / Chain / Secondary structure / pLDDT B-factor) are selectable from dropdowns. Toolbar buttons open the current structure in PyMol, VMD, or ChimeraX, or reveal the file in the system file manager.

> **Note:** 3DMol requires `3Dmol-min.js` in the `bioforge_app/` directory. Download from https://3dmol.org or bundle via npm. Without it, a placeholder is shown but all other panels function normally.

#### 🧬 Sequences

Loads `proteins.faa` and all split FASTA files from `proteinmpnn_out/split_seqs/`. Sequences are rendered with per-residue colour coding:

| Class | Colour | Residues |
|---|---|---|
| Positive | Blue | H K R |
| Negative | Red | D E |
| Hydrophobic | Amber | A V L I M F W Y |
| Polar | Green | S T N Q |
| Special | Purple | P G |
| Cysteine | Gold | C |

A search box filters by sequence ID or sequence content. Copy All and Open in JalView buttons are in the toolbar.

#### 🧫 Plasmids

Loads all GenBank files from `moclo_plasmids/`. Clicking a file renders a linear feature map canvas and a feature table below it. Feature types are colour-coded (CDS amber, gene cyan, promoter green, rep_origin purple, terminator orange, etc.). Open in SnapGene Viewer and Reveal in file manager buttons are in the toolbar.

#### 💥 BLAST

Parses `blast_results.txt` into a filterable, sortable table. Columns: Query, Subject, %ID, Length, E-value, Bits, Function. Filter by any text across query / subject / function. Sort by E-value, %ID, or query name.

#### 📁 Files

A full file browser with address bar navigation.

- **Address bar** - always shows the current directory; type any path and press Enter or Go to jump directly to it
- **← Back** - browser-style history navigation through all directories visited this session
- **↑ Up** - navigate to parent directory
- **Browse…** - OS folder picker
- **TREE / FLAT toggle** - switch between flat navigator mode (click a folder to enter it, address bar updates) and tree mode (folders expand inline with ▶ / ▼ chevrons). The button glows cyan while tree mode is active.

Clicking a file previews it in the right pane. PDB, FASTA, and GenBank files route to their specialised viewers (Structures / Sequences / Plasmids tabs). TSV and CSV files render as a sortable table. All other text files display as plain text.

#### 🔄 Feedback

Cards for all six feedback loops, each showing the loop name, a description, configurable parameter fields (pre-filled with sensible defaults), and Run / View Output buttons. A shared log at the bottom shows live output from whichever loop is running. A Kill button stops the active loop. Cards are automatically marked DONE if their sentinel file is detected in the loaded run directory.

#### ⚙ Settings

Persistent settings for all paths and external tool commands.

| Setting | Default |
|---|---|
| Conda profile.sh | `~/miniconda3/etc/profile.d/conda.sh` |
| Scripts directory | App directory |
| Default output dir | `~/output` |
| GeneMark bin dir | `~/genemark-etp-full/gmetp_linux_64/bin` |
| JalView command | `jalview` |
| SnapGene Viewer | `snapgene-viewer` |
| PyMol command | `pymol` |
| VMD command | `vmd` |
| ChimeraX command | `chimerax` |

---

## Usage

### Recommended: run everything with the orchestrator

`genomopipe.py` is the single entry point. It runs Phases 1–3 in order, handles resumability via per-phase sentinel files, and automatically passes config values to each downstream script. All keys can be supplied in a config file (YAML, JSON, or plain text), as CLI flags, or both - CLI flags always win.

```bash
# Full run from a config file
python genomopipe.py genomopipe_config.yaml

# Minimal CLI invocation - prokaryote
python genomopipe.py --organism "Escherichia coli" --output_dir ./output

# Eukaryote with automatic RNA-Seq download
python genomopipe.py \
    --organism "Taraxacum officinale" \
    --output_dir ./output \
    --is_eukaryote true \
    --auto_rnaseq

# Load a base config, override individual keys on the CLI
python genomopipe.py genomopipe_config.yaml \
    --organism "Arabidopsis thaliana" \
    --fb6_dry_run

# Resume an interrupted run - completed phases are skipped automatically
python genomopipe.py genomopipe_config.yaml

# Skip Phase 1 (use an existing completed run)
python genomopipe.py genomopipe_config.yaml --skip_phase1

# Run only feedback loops on an existing run
python genomopipe.py genomopipe_config.yaml --skip_phase1 --skip_phase2

# Force Phase 1 to restart from scratch (clears genome_to_design.sh checkpoints)
python genomopipe.py genomopipe_config.yaml --force

# Clear only orchestrator sentinels, keep genome_to_design.sh checkpoints
python genomopipe.py genomopipe_config.yaml --reset
```

### Orchestrator CLI reference

| Flag | Type | Default | Description |
|---|---|---|---|
| `config` | positional | - | YAML / JSON / .txt config file |
| `--organism` | string | - | Organism name (quoted) or NCBI TaxID |
| `--output_dir` | path | `./output` | Root output directory |
| `--scripts_dir` | path | `.` | Directory containing all pipeline scripts |
| `--email` | string | - | NCBI Entrez email (required for FB6 taxonomy lookups) |
| `--is_eukaryote` | bool | `false` | `true` for eukaryotes |
| `--genemark_path` | path | `~/genemark-etp-full/.../bin` | GeneMark bin dir |
| `--bam` | path | - | Existing RNA-Seq BAM for BRAKER hints |
| `--auto_rnaseq` | flag | false | Auto-download RNA-Seq from SRA |
| `--force` | flag | false | Clear genome_to_design.sh checkpoints and restart Phase 1 |
| `--moclo_standard` | string | `marillonnet` | `marillonnet` \| `cidar` \| `jump` |
| `--enzyme_level0` | string | `BsaI-HFv2` | Level 0 domestication enzyme |
| `--enzyme_level1` | string | `BpiI` | Level 1 assembly enzyme |
| `--perform_domestication` | bool | `true` | Run CDS domestication in Phase 2 |
| `--output_prefix` | string | `moclo_plasmid` | Output file stem for Phase 2 |
| `--genes` | paths | - | CDS FASTA paths; auto-discovered from split_seqs if omitted |
| `--codon_optimize` | flag | false | Enable codon optimization / back-translation before domestication |
| `--expression_host` | string | `e_coli` | Target expression host: `e_coli`, `s_cerevisiae`, `h_sapiens`, `p_pastoris`, `b_subtilis` |
| `--codon_optimize_method` | string | `auto` | `auto` \| `max_frequency` \| `dnachisel` |
| `--fb6_min_hits` | int | `5` | FB6: min BLAST hits needed to trigger partition switch |
| `--fb6_evalue_cutoff` | float | `1e-5` | FB6: max e-value for a hit to count |
| `--dry_run` | flag | false | FB6: audit only, do not re-run BRAKER |
| `--no_fb6` | flag | false | Skip Feedback Loop 6 |
| `--skip_phase1` | flag | false | Skip genome_to_design.sh (requires existing run) |
| `--skip_phase2` | flag | false | Skip plasmid design |
| `--skip_feedback` | flag | false | Skip all feedback loops |
| `--reset` | flag | false | Clear orchestrator sentinels before running |

---

## Running Scripts Individually

Each script can also be run standalone when iterating on a single phase.

**Phase 1:**
```bash
./genome_to_design.sh "Organism name" /path/to/output [true|false] [options]

./genome_to_design.sh "Escherichia coli" ./output
./genome_to_design.sh "Taraxacum officinale" ./output true --auto_rnaseq
./genome_to_design.sh "Arabidopsis thaliana" ./output true --bam=/data/rnaseq.bam
./genome_to_design.sh "Homo sapiens" ./output true \
    --GENEMARK_PATH=/opt/genemark/bin --auto_rnaseq
./genome_to_design.sh "Mus musculus" ./output true --force
```

**Phase 2:**
```bash
python plasmid_design_moclo_v3.py config.yaml   # or .json or .txt
```

**Feedback loops (all accept config file or CLI flags):**
```bash
# FB1 - motif re-scaffolding
./feedback1_colabfold_to_rfdiffusion.sh config.yaml
./feedback1_colabfold_to_rfdiffusion.sh \
    --run_dir=./output/latest --top_n=3 --iterations=2

# FB2 - pLDDT-gated resampling
python feedback2_plddt_mpnn_resample.py config.yaml
python feedback2_plddt_mpnn_resample.py \
    --run_dir=./output/latest --plddt_pass=75 --max_iterations=3

# FB3 - enriched BRAKER hints
./feedback3_blast_to_braker.sh config.yaml
./feedback3_blast_to_braker.sh \
    --run_dir=./output/latest --evalue_cutoff=1e-10 --pident_min=40

# FB4 - domestication re-validation
python feedback4_domesticated_cds_revalidate.py config.yaml
python feedback4_domesticated_cds_revalidate.py \
    --run_dir=./output/latest --plddt_drop_warn=5.0 --plddt_drop_fail=15.0

# FB5 - designs as annotation hints
./feedback5_designed_proteins_to_annotation.sh config.yaml
./feedback5_designed_proteins_to_annotation.sh \
    --run_dir=./output/latest --plddt_min=75

# FB6 - taxonomy partition correction (also run by orchestrator)
python feedback6_blast_taxonomy_rerun.py config.yaml
python feedback6_blast_taxonomy_rerun.py \
    --run_dir=./output/latest --min_hits=5 --dry_run
```

---

## Pipeline Steps In Detail

### Phase 1: genome_to_design.sh (Steps 1–8)

#### Step 1. Genome Fetch (`genome_fetch.sh`)
Queries NCBI for the organism by name or TaxID, determines the relevant taxonomic group, expands to all descendant TaxIDs for broader coverage, and downloads the best available genome assembly in priority order:
1. RefSeq reference assemblies
2. Any RefSeq assembly
3. Any GenBank assembly
4. Falls back to fetching individual nucleotide sequences and clustering with `cd-hit-est`

#### Step 1b. FASTA Header Sanitization
Replaces spaces and pipe characters in sequence headers with underscores to prevent downstream tool failures.

#### Step 2. QC with BBTools
Runs `bbduk.sh` to trim adapter sequences. Falls back to the raw sanitized genome if this step fails.

#### Step 2b. Repeat Masking
Runs `RepeatModeler` to build a de novo repeat library, then `RepeatMasker` to soft-mask the genome. Skipped automatically for genomes under 1 Mbp. Critical for eukaryotes (failure halts the pipeline when `is_eukaryote=true`).

For large genomes (>500 Mbp) this step takes 12–30 hours. Progress is logged to `logs/pipeline.log` via `[REPEATMASK]` markers at each phase (BuildDatabase, RepeatModeler start, round headers, library size, RepeatMasker start, masked percentage). Full unfiltered RepeatModeler output - including the per-second countdown - goes to `logs/repeatmodeler_detail.log` so nothing is lost.

#### Step 3. Gene Annotation
Eukaryotes use BRAKER (Augustus + GeneMark-ETP) with OrthoDB protein hints auto-selected by taxonomic lineage (fungi, viridiplantae, vertebrata, arthropoda, etc.) and optionally supplemented with RNA-Seq evidence. The `--auto_rnaseq` mode queries SRA via Entrez, downloads with `prefetch` + `fasterq-dump`, aligns with STAR, and merges BAMs with `samtools`. Falls back to proteins-only if no RNA-Seq is found. If BRAKER fails due to insufficient intron evidence, the pipeline automatically retries in `--esmode`. Prokaryotes use Prokka.

#### Step 4. Protein Extraction and Filtering
Keeps only the longest isoform per gene; removes proteins shorter than 100 aa (200 aa for eukaryotes). Outputs `proteins.faa`.

#### Step 5. RFdiffusion Backbone Design
Generates novel protein backbones. Uses a pre-computed ColabFold structure from `colabfold_pre/` as a fixed motif if present; otherwise runs unconditional de novo design. Uses GNU `parallel` if available.

#### Step 5b. ProteinMPNN Sequence Design
Designs 8 amino acid sequence candidates per backbone at sampling temperature 0.1. Splits outputs into individual FASTA files for ColabFold input.

#### Step 6. ColabFold Structure Prediction
Predicts and relaxes structures with 3 recycle passes using the AlphaFold2 backend.

#### Step 7. BLAST Validation
Remote `blastp` against NCBI `nr` for each designed sequence. Combines all results into `blast_results.txt`.

#### Step 8. Rename Outputs
Renames ColabFold output files using the putative function name from BLAST results.

---

### Phase 2: MoClo Plasmid Design (`plasmid_design_moclo_v3.py`) - Step 9

Takes designed sequences from `proteinmpnn_out/split_seqs/` and produces lab-ready Golden Gate assembly files.

1. **Protein detection** - automatically detects whether input sequences are amino-acid (protein FASTA from ProteinMPNN) or DNA. Protein inputs are always back-translated to DNA before downstream steps.
2. **Codon optimization** *(new — if `codon_optimize: true`)* - re-encodes CDSs for the target expression host before restriction-site removal. For protein inputs this is the back-translation step itself; for DNA inputs it applies synonymous substitutions genome-wide. See [Codon Optimization](#codon-optimization) for full details.
3. **CDS domestication** - silently removes internal restriction enzyme recognition sites using synonymous codon substitution scored by host codon frequency. Uses DNA Chisel if installed, otherwise a built-in heuristic fallback.
4. **Fusion-site addition** - wraps each part with the correct 4-bp MoClo overhangs for the chosen standard (Marillonnet, CIDAR, or JUMP).
5. **In-silico assembly** - assembles promoter → RBS → CDS → terminator → backbone and verifies a valid circular product.
6. **Protocol summary** - prints suggested thermocycler conditions for the enzyme pair.
7. **Outputs** - `.gb` GenBank file (SnapGene-ready) and `.fasta` for synthesis ordering, written to `<run_dir>/moclo_plasmids/`.

When run via the orchestrator, `genes` are auto-discovered from `proteinmpnn_out/split_seqs/` if not explicitly listed in the config.

---

## Codon Optimization

Codon optimization recodes each designed CDS using synonymous codons preferred by the target expression host. This is distinct from CDS domestication (which only removes restriction sites): codon optimization considers the entire sequence, maximising expression potential by selecting codons that match the host's translational machinery.

### When to use it

Enable codon optimization whenever the designed protein will be expressed in a host other than the organism it was originally designed against — or any time you want maximum expression yield. It is particularly important when:

- Expressing a designed protein in *E. coli* after designing from a eukaryotic genome
- Targeting *P. pastoris* or CHO cells for secreted protein production
- Minimising the use of rare codons that slow ribosomes or cause frameshifts

### How it works

The codon optimization step runs **before** CDS domestication in Phase 2. The pipeline:

1. **Detects sequence type** — protein FASTA from ProteinMPNN is back-translated directly to optimized DNA; DNA inputs are re-encoded via synonymous substitution.
2. **Selects codons** from the host's codon frequency table using one of three methods (see below).
3. **Passes the optimized DNA** to the existing domestication step, which then removes any restriction sites the new codon choices may have introduced.

### Methods

| Method | When used | Description |
|---|---|---|
| `auto` | default | Tries DNA Chisel first, then python-codon-tables, then built-in tables |
| `max_frequency` | explicit | Always uses the most-frequent codon per amino acid from the host table |
| `dnachisel` | explicit | Forces DNA Chisel's constraint-satisfaction solver (error if not installed) |

**DNA Chisel** (`pip install dnachisel[reports]`) gives the best results because it simultaneously satisfies codon-optimization objectives and restriction-site avoidance constraints in a single pass, trading slightly lower average codon frequency for globally fewer conflicts.

**python-codon-tables** (`pip install python-codon-tables`) uses the most recent NCBI codon usage database tables and supports all organisms with sufficient NCBI entries. Recommended over the built-in tables when your expression host is not one of the five pre-loaded organisms.

**Built-in tables** (always available, no extra install) cover five common expression hosts derived from Kazusa codon usage database values.

### Supported expression hosts

| Key | Organism | Use case |
|---|---|---|
| `e_coli` | *Escherichia coli* K-12 | Bacterial expression (default) |
| `s_cerevisiae` | *Saccharomyces cerevisiae* | Yeast expression, secretion |
| `h_sapiens` | *Homo sapiens* | Human / CHO mammalian expression |
| `p_pastoris` | *Komagataella phaffii* (Pichia) | High-yield secreted protein |
| `b_subtilis` | *Bacillus subtilis* 168 | Gram-positive / spore-display |

Common aliases are also accepted (e.g. `ecoli`, `yeast`, `human`, `cho`, `pichia`).

For any other organism, install `python-codon-tables` and pass the NCBI species name as `expression_host` — the library will fetch the correct table automatically.

### Configuration

```yaml
codon_optimize:        true         # false by default
expression_host:       e_coli       # see table above
codon_optimize_method: auto         # auto | max_frequency | dnachisel
```

```bash
# CLI equivalents
python genomopipe.py config.yaml \
    --codon_optimize \
    --expression_host s_cerevisiae \
    --codon_optimize_method max_frequency
```

### Combining with domestication

Both steps are applied in sequence and are independently controllable:

```yaml
codon_optimize:        true    # re-encode for host first
perform_domestication: true    # then remove restriction sites
expression_host:       e_coli
```

When both are enabled and DNA Chisel is installed, the `auto` method runs a single unified pass that satisfies codon-optimization objectives and restriction-site avoidance constraints simultaneously, which is more efficient and typically produces fewer synonymous changes overall.

### Feedback Loop 4 interaction

FB4 (`feedback4_domesticated_cds_revalidate.py`) re-validates domesticated sequences with ColabFold to catch any pLDDT regressions caused by synonymous changes. This is equally relevant after codon optimization — run FB4 after Phase 2 whenever `codon_optimize: true` to confirm that re-encoding did not disrupt co-translational folding.

---

## Feedback Loops

Each loop reads outputs from earlier phases, feeds information back upstream, and writes results into a dedicated `feedback<N>_loop/` subdirectory. Each loop is independently resumable. FB6 runs automatically as part of the orchestrator; FB1–FB5 are run as standalone scripts or from the GUI Feedback tab.

### FB1 - ColabFold → RFdiffusion motif re-scaffolding (`feedback1_colabfold_to_rfdiffusion.sh`)

Ranks all ColabFold predictions by mean pLDDT and selects the top-N structures. Each is fed back to RFdiffusion as a fixed motif, with variable-length flanking regions scaffolded around it. ProteinMPNN designs sequences for the new backbones and ColabFold re-evaluates them. Repeats for a configurable number of iterations. The "exploit the best hits" strategy: rather than treating all RFdiffusion designs equally, it focuses redesign effort on structures that already showed high confidence, converging toward better scaffolds with each cycle.

| Key | Default | Description |
|---|---|---|
| `top_n` | 3 | Number of top ColabFold structures to use as motifs |
| `plddt_min` | 70 | Minimum mean pLDDT to qualify a structure |
| `motif_residues` | `A1-50` | Fixed residue range passed to RFdiffusion |
| `flank_min` / `flank_max` | 10 / 40 | Scaffolded flanking region lengths |
| `iterations` | 2 | Number of ColabFold → RFdiffusion cycles |
| `num_designs` | 4 | RFdiffusion designs per motif per iteration |
| `num_seqs` | 8 | ProteinMPNN sequences per backbone |

**Outputs:** `feedback1_loop/iteration_<N>/` with new designs, sequences, and structures; `fb1_summary.tsv` comparing original vs. redesigned pLDDT.

### FB2 - pLDDT filter → ProteinMPNN resampling (`feedback2_plddt_mpnn_resample.py`)

Splits all ColabFold predictions into PASS (mean pLDDT ≥ threshold) and RETRY (below threshold). For RETRY designs, the corresponding RFdiffusion PDB is sent back to ProteinMPNN at a slightly higher sampling temperature to diversify the sequence search. ColabFold re-evaluates the new candidates. Repeats up to `max_iterations`; designs that never pass are flagged as unconverged. The key insight is that a poor ColabFold prediction reflects a bad sequence, not a bad backbone - more diverse sequence sampling often finds a sequence that folds cleanly onto the same scaffold.

| Key | Default | Description |
|---|---|---|
| `plddt_pass` | 75 | Mean pLDDT threshold to accept a design |
| `plddt_warn` | 60 | Below this, a warning is printed even if passing |
| `max_iterations` | 3 | Max resample rounds before marking unconverged |
| `resample_temp` | 0.2 | ProteinMPNN sampling temperature on retry |
| `resample_n` | 16 | New sequences per retry backbone |

**Outputs:** `feedback2_loop/iteration_<N>/`; `fb2_summary.tsv` with design ID, iteration reached, final pLDDT, pass/fail.

### FB3 - BLAST hits → BRAKER re-annotation (`feedback3_blast_to_braker.sh`)

Parses `blast_results.txt` for high-confidence hits and fetches the full subject sequences from NCBI. De-duplicates and combines these with the original OrthoDB hint file, then re-runs BRAKER with the enriched protein hint set. Compares gene counts and mean protein length between the original and enriched annotations and writes a diff report. The initial OrthoDB partition covers a broad taxonomic group; BLAST hits from designed sequences reveal which specific protein families are actually present - families that may be under-represented in the generic partition. Adding them as direct hints can recover missed exons and fix gene merges and splits in exactly the protein classes being designed.

| Key | Default | Description |
|---|---|---|
| `evalue_cutoff` | 1e-10 | Max e-value for a BLAST hit to be included as a hint |
| `pident_min` | 40 | Min percent identity for inclusion |
| `max_hits_per_query` | 5 | Max subject sequences fetched per query |
| `bam_file` | - | RNA-Seq BAM to carry forward into re-annotation |

**Outputs:** `feedback3_loop/blast_subjects.faa`, `enriched_hints.faa`, `braker_enriched/`, `proteins_enriched.faa`, `fb3_annotation_diff.txt`.

**To continue the design pipeline from enriched proteins:**
```bash
cp ./output/latest/feedback3_loop/proteins_enriched.faa \
   ./output/latest/proteins.faa
rm ./output/latest/.step5_rfdiffusion.done   # and any later checkpoints
./genome_to_design.sh "Organism" ./output true   # resumes from Step 5
```

### FB4 - Domesticated CDS → ColabFold re-validation (`feedback4_domesticated_cds_revalidate.py`)

Reads the `.gb` file(s) produced by Step 9 and extracts CDSs that were domesticated. Calculates the number and positions of synonymous changes relative to the pre-domestication source. Runs ColabFold on the domesticated amino acid sequences and compares pLDDT scores against the original Step 6 prediction. Flags sequences where pLDDT dropped more than the configured threshold. Synonymous codon substitutions can alter local mRNA structure and translation rate, occasionally affecting co-translational folding. This loop catches those regressions before a sequence reaches synthesis.

| Key | Default | Description |
|---|---|---|
| `plasmid_dir` | `<run_dir>/moclo_plasmids` | Where to read `.gb` files |
| `plddt_drop_warn` | 5.0 | pLDDT drop that triggers a WARNING |
| `plddt_drop_fail` | 15.0 | pLDDT drop that triggers a FAIL flag |

**Outputs:** `feedback4_loop/domesticated_seqs/`, `colabfold_dom/`, `fb4_report.tsv` (num_changes, original_pLDDT, dom_pLDDT, delta, flag).

### FB5 - Validated designs → annotation hints (`feedback5_designed_proteins_to_annotation.sh`)

Filters ColabFold-predicted designed proteins by pLDDT threshold and combines the passing sequences with the original OrthoDB hint file to produce `designed_hints.faa`. This enriched hint set can be used two ways: (a) re-run BRAKER on the *same* organism to recover genes in exactly the protein families represented by the designs; (b) pass as `--prot_seq` when launching a fresh run on a *related* organism, bootstrapping annotation quality with curated high-confidence sequences rather than generic OrthoDB proteins. Optionally triggers the target-organism run automatically via `target_organism`.

| Key | Default | Description |
|---|---|---|
| `plddt_min` | 75 | Min pLDDT to include a design as a hint |
| `target_organism` | - | If set, automatically launch genome_to_design.sh on a new species |
| `target_output` | - | Output directory for the target organism run |
| `target_eukaryote` | true | Eukaryote flag for the target organism |
| `include_blast_func` | true | Annotate hint sequences with BLAST function labels |

**Outputs:** `feedback5_loop/designed_hints.faa`, `fb5_manifest.tsv`.

### FB6 - BLAST taxonomy → BRAKER re-run with corrected OrthoDB partition (`feedback6_blast_taxonomy_rerun.py`)

Parses `blast_results.txt` and fetches the taxonomic lineage of each subject accession. Each accession is assigned to its single most-specific OrthoDB partition keyword using a fixed specificity order. Compares the dominant partition against the one actually used in Step 3 by reading the step log. Partition aliases (e.g. `mammalia` and `vertebrata` both map to slug `Vertebrata`) are resolved before comparison. If a different, more specific partition accumulates ≥ `min_hits` accessions, downloads it, re-runs BRAKER with the corrected hints, runs GTF repair and gffread, and filters proteins. Writes a taxonomy audit report listing per-partition hit counts, original and suggested partitions, and whether a re-run was triggered.

| Key | Default | Description |
|---|---|---|
| `min_hits` | 5 | Min accessions in a partition to trigger a switch |
| `evalue_cutoff` | 1e-5 | Max BLAST e-value to count an accession |
| `dry_run` | false | Audit only; do not download or re-run BRAKER |

**Outputs:** `feedback6_loop/feedback6_taxonomy_audit.txt`, `braker_corrected/`, `proteins_corrected.faa`.

**To continue the design pipeline from corrected proteins:**
```bash
cp ./output/latest/feedback6_loop/proteins_corrected.faa \
   ./output/latest/proteins.faa
rm ./output/latest/.step5_rfdiffusion.done   # and any later checkpoints
python genomopipe.py genomopipe_config.yaml --skip_phase1   # resumes from Step 5
```

---

## Resumability and Checkpoints

### Orchestrator sentinels

The orchestrator writes a `.done` sentinel in the run directory after each phase completes:

| Sentinel | Phase |
|---|---|
| `.genomopipe_phase1.done` | genome_to_design.sh |
| `.genomopipe_phase2.done` | plasmid_design_moclo_v3.py |
| `.genomopipe_feedback6.done` | Feedback Loop 6 |

Re-running `genomopipe.py` with the same config skips any phase whose sentinel exists. To rerun a specific phase, delete its sentinel:

```bash
rm ./output/latest/.genomopipe_phase2.done
python genomopipe.py genomopipe_config.yaml
```

To clear all orchestrator sentinels without touching step checkpoints:
```bash
python genomopipe.py genomopipe_config.yaml --reset
```

### genome_to_design.sh step checkpoints

Within Phase 1, each step has its own `.step*.done` / `.step*.failed` checkpoint. Re-running picks up from the last successful step.

```bash
# Resume an interrupted run
./genome_to_design.sh "Taraxacum officinale" ./output true

# Re-run a single step
rm ./output/latest/.step4_protein_extract.done
./genome_to_design.sh "Taraxacum officinale" ./output true

# Restart Phase 1 entirely
./genome_to_design.sh "Taraxacum officinale" ./output true --force
# or via orchestrator:
python genomopipe.py genomopipe_config.yaml --force
```

### Sentinel reference

| File | Meaning |
|---|---|
| `.step1_genome_fetch.done` | Step completed (contains timestamp + elapsed time) |
| `.step3_annotate.failed` | Step failed (contains exit code + timestamp) |
| `.genomopipe_phase1.done` | Steps 1–8 complete (orchestrator) |
| `.genomopipe_phase2.done` | Step 9 complete (orchestrator) |
| `.genomopipe_feedback6.done` | FB6 complete |
| `.pipeline_complete.done` | Full pipeline including Step 9 complete |

`--force` clears both `.done` and `.failed` sentinels so stale failure markers from previous runs do not persist.

### Feedback loop checkpoints

Each feedback loop script uses the same `run_step` mechanism as `genome_to_design.sh`, writing `.step*.done` sentinels inside its own `feedback<N>_loop/` directory. Re-running a feedback script resumes from the last completed step within that loop.

---

## Config File Reference

All scripts accept YAML, JSON, or plain-text `.txt` config files. The orchestrator config covers all phases. All keys are optional - omitted keys fall back to the defaults shown in the tables above. Three ready-to-use example files are included: `genomopipe_config.yaml`, `genomopipe_config.json`, and `genomopipe_config.txt`.

### YAML (recommended)

```yaml
# ── Shared / global ─────────────────────────────────────────────────────────
output_dir:   ./output
scripts_dir:  .               # directory containing all pipeline scripts
email:        you@example.com # required for FB6 taxonomy lookups

# ── Phase 1: genome_to_design.sh ────────────────────────────────────────────
organism:      "Taraxacum officinale"
is_eukaryote:  true
genemark_path: ~/genemark-etp-full/gmetp_linux_64/bin
# bam: /data/rnaseq_merged.bam
auto_rnaseq:   false
force:         false

# ── Phase 2: plasmid_design_moclo_v3.py ─────────────────────────────────────
assembly_method:       GoldenGate
moclo_standard:        marillonnet   # marillonnet | cidar | jump
enzyme_level0:         BsaI-HFv2
enzyme_level1:         BpiI
promoter:              pJ23119
rbs:                   RBS_strong
terminator:            rrnB_T1
ori:                   pMB1_ori
marker:                KanR
# backbone: pUC19_prepared.fasta     # overrides ori + marker
perform_domestication: true
output_prefix:         moclo_plasmid
genes: []   # leave empty to auto-discover from proteinmpnn_out/split_seqs
            # or list explicit paths:
            # genes:
            #   - proteinmpnn_out/split_seqs/design_0_seq0.fasta

# ── Codon optimization (Phase 2, before domestication) ────────────────────────
codon_optimize:        false   # true → re-encode CDSs for expression_host
expression_host:       e_coli  # e_coli | s_cerevisiae | h_sapiens | p_pastoris | b_subtilis
codon_optimize_method: auto    # auto | max_frequency | dnachisel

# ── Phase 3: Feedback Loop 6 (orchestrator-integrated) ──────────────────────
fb6_min_hits:      5
fb6_evalue_cutoff: 1.0e-5
fb6_dry_run:       false   # true → audit only, do not re-run BRAKER

# ── Orchestrator control ─────────────────────────────────────────────────────
skip_phase1:   false
skip_phase2:   false
skip_feedback: false
run_fb6:       true
```

### JSON

```json
{
  "output_dir":   "./output",
  "scripts_dir":  ".",
  "email":        "you@example.com",
  "organism":     "Taraxacum officinale",
  "is_eukaryote": true,
  "auto_rnaseq":  false,
  "force":        false,
  "assembly_method":       "GoldenGate",
  "moclo_standard":        "marillonnet",
  "enzyme_level0":         "BsaI-HFv2",
  "enzyme_level1":         "BpiI",
  "promoter":              "pJ23119",
  "rbs":                   "RBS_strong",
  "terminator":            "rrnB_T1",
  "ori":                   "pMB1_ori",
  "marker":                "KanR",
  "perform_domestication": true,
  "output_prefix":         "moclo_plasmid",
  "genes":                 [],
  "codon_optimize":        false,
  "expression_host":       "e_coli",
  "codon_optimize_method": "auto",
  "fb6_min_hits":      5,
  "fb6_evalue_cutoff": 1e-5,
  "fb6_dry_run":       false,
  "skip_phase1":   false,
  "skip_phase2":   false,
  "skip_feedback": false,
  "run_fb6":       true
}
```

### Plain text (`.txt`)

```
# key: value pairs. Lists use bare values under a section header ending in ':'.
output_dir: ./output
scripts_dir: .
email: you@example.com

organism: Taraxacum officinale
is_eukaryote: true
auto_rnaseq: false
force: false

assembly_method: GoldenGate
moclo_standard: marillonnet
enzyme_level0: BsaI-HFv2
enzyme_level1: BpiI
promoter: pJ23119
rbs: RBS_strong
terminator: rrnB_T1
ori: pMB1_ori
marker: KanR
perform_domestication: true
output_prefix: moclo_plasmid

# Leave the section empty to auto-discover genes from split_seqs
genes:
proteinmpnn_out/split_seqs/design_0_seq0.fasta
proteinmpnn_out/split_seqs/design_1_seq3.fasta

fb6_min_hits: 5
fb6_evalue_cutoff: 1e-5
fb6_dry_run: false

skip_phase1: false
skip_phase2: false
skip_feedback: false
run_fb6: true
```

> For complex configs (nested `enzyme_props`, custom `overhang_table`) YAML or JSON are preferred over plain text.

### MoClo plasmid config (`plasmid_design_moclo_v3.py` standalone)

Extended YAML example for standalone use:

```yaml
assembly_method: GoldenGate
moclo_standard: marillonnet
output_prefix: beet_custom_gg_v3
enzyme_level0: BsaI-HFv2
enzyme_level1: BpiI

enzyme_props:           # optional enzyme property overrides
  BsaI-HFv2:
    recog: GGTCTC
    cut: "1/5"
    overhang_len: 4
    opt_temp: 50
    inact_temp: 65
    buffer: CutSmart

promoter:   pJ23119
rbs:        RBS_B0034
terminator: rrnB_T1
marker:     KanR
ori:        p15A

overhang_table:         # optional: override moclo_standard preset
  promoter:   [GGAG, TACT]
  rbs:        [AATG, GCTT]
  cds:        [AATG, GCTT]
  terminator: [TACT, CTGC]

genes:
  - proteinmpnn_out/split_seqs/design_0_seq0.fasta
  - proteinmpnn_out/split_seqs/design_1_seq3.fasta

perform_domestication: true
domestication_enzymes: [BsaI, BpiI]
output_dir: ./custom_plasmids
```

### Full config key reference

| Key | Default | Description |
|---|---|---|
| `output_dir` | `./output` | Root directory for all run subdirectories |
| `scripts_dir` | `.` | Directory containing all pipeline scripts |
| `email` | - | NCBI Entrez email (required for FB6) |
| `organism` | - | Organism name (quoted) or NCBI TaxID |
| `is_eukaryote` | `false` | `true` for eukaryotes |
| `genemark_path` | `~/genemark-etp-full/.../bin` | GeneMark bin directory |
| `bam` | - | Existing RNA-Seq BAM path |
| `auto_rnaseq` | `false` | Auto-download RNA-Seq from SRA |
| `force` | `false` | Clear genome_to_design.sh checkpoints and restart |
| `assembly_method` | `GoldenGate` | Assembly strategy |
| `moclo_standard` | `marillonnet` | Overhang preset: `marillonnet`, `cidar`, or `jump` |
| `enzyme_level0` | `BsaI-HFv2` | Level 0 domestication enzyme |
| `enzyme_level1` | `BpiI` | Level 1 assembly enzyme |
| `promoter` | `pJ23119` | Promoter part name or path |
| `rbs` | `RBS_strong` | RBS part name or path |
| `terminator` | `rrnB_T1` | Terminator part name or path |
| `ori` | `pMB1_ori` | Origin of replication |
| `marker` | `KanR` | Antibiotic resistance marker |
| `backbone` | `<ori>+<marker>` | Full backbone path (overrides ori + marker) |
| `perform_domestication` | `true` | Run CDS domestication |
| `output_prefix` | `moclo_plasmid` | Output file stem |
| `genes` | auto-discovered | List of CDS FASTA paths |
| `codon_optimize` | `false` | Enable codon optimization before domestication |
| `expression_host` | `e_coli` | Target host: `e_coli`, `s_cerevisiae`, `h_sapiens`, `p_pastoris`, `b_subtilis` |
| `codon_optimize_method` | `auto` | `auto` \| `max_frequency` \| `dnachisel` |
| `fb6_min_hits` | `5` | Min accessions in a partition to trigger a switch |
| `fb6_evalue_cutoff` | `1e-5` | Max BLAST e-value to count |
| `fb6_dry_run` | `false` | Audit only; do not re-run BRAKER |
| `skip_phase1` | `false` | Skip genome_to_design.sh |
| `skip_phase2` | `false` | Skip plasmid design |
| `skip_feedback` | `false` | Skip all feedback loops |
| `run_fb6` | `true` | Run Feedback Loop 6 via orchestrator |

---

## Path Resolution

`plasmid_design_moclo_v3.py` resolves all paths dynamically so it works on any machine without editing the script. Priority order (highest to lowest):

| Priority | Source | Example |
|---|---|---|
| 1 | Config file key | `parts_dir: /data/split_seqs` |
| 2 | Environment variable | `export MOCLO_PARTS_DIR=/data/split_seqs` |
| 3 | Auto-discovery | Follows `latest/` symlink, then newest `run_*` by mtime |
| 4 | Current working directory | Wherever the script is invoked from |

| Environment variable | Controls | Default |
|---|---|---|
| `MOCLO_OUTPUT_BASE` | Root of all pipeline runs | `~/output` |
| `MOCLO_RUN_DIR` | Explicit run directory | Newest `run_*` under output_base |
| `MOCLO_PARTS_DIR` | Explicit split_seqs path | Auto-searched under run_dir |

```bash
# Default - works automatically after genome_to_design.sh has run once
python plasmid_design_moclo_v3.py config.yaml

# Pipeline outputs are on a different disk
export MOCLO_OUTPUT_BASE=/mnt/storage/genome_runs
python plasmid_design_moclo_v3.py config.yaml

# Point to a specific older run
export MOCLO_RUN_DIR=/mnt/storage/genome_runs/run_20240101_0900
python plasmid_design_moclo_v3.py config.yaml
```

---

## Output Structure

```
output/
├── genomopipe_orchestrator.log        # Master orchestrator log
├── genomopipe_phase1.log              # Phase 1 stdout/stderr
├── genomopipe_phase2.log              # Phase 2 stdout/stderr
├── genomopipe_feedback6.log           # FB6 stdout/stderr
├── latest -> run_20240315_1430/       # Symlink to most recent run
└── run_20240315_1430/
    │
    ├── .genomopipe_phase1.done        # Orchestrator sentinels
    ├── .genomopipe_phase2.done
    ├── .genomopipe_feedback6.done
    ├── .step*.done / .step*.failed    # genome_to_design.sh step checkpoints
    │
    ├── Taraxacum_officinale/          # Raw genome files (.fna)
    ├── clean.fna                      # QC-cleaned genome
    ├── masked.fna                     # Repeat-masked genome (eukaryotes)
    │
    ├── braker_out/
    │   ├── braker.gtf
    │   ├── braker_fixed.gtf           # GTF with repaired transcript_id fields
    │   ├── Taraxacum_officinale.faa   # All annotated proteins (all isoforms)
    │   └── genome.fa                  # BRAKER's header-sanitized genome copy
    ├── proteins.faa                   # Filtered proteins (longest isoform, >100aa)
    │
    ├── designs/                       # RFdiffusion backbone PDBs
    ├── proteinmpnn_out/
    │   ├── seqs/                      # Multi-sequence .fa per backbone
    │   └── split_seqs/                # Single-sequence .fasta files for ColabFold
    ├── colabfold_out/                 # Predicted/relaxed structures (.pdb, .json)
    │
    ├── blast_results/                 # Per-design BLAST outputs
    ├── blast_results.txt              # Combined BLAST summary (tabular)
    │
    ├── moclo_plasmids/
    │   ├── moclo_plasmid.gb           # GenBank file (SnapGene-ready)
    │   └── moclo_plasmid.fasta        # FASTA for synthesis ordering
    │
    ├── RM_*/                          # RepeatModeler working directory
    │   ├── round-1/ … round-N/
    │   └── consensi.fa.classified     # Final classified repeat library
    ├── repeat_out/                    # RepeatMasker output and statistics
    │
    ├── feedback1_loop/
    │   ├── iteration_1/
    │   │   ├── designs/
    │   │   ├── proteinmpnn_out/
    │   │   └── colabfold_out/
    │   ├── fb1_summary.tsv
    │   └── logs/
    │
    ├── feedback2_loop/
    │   ├── iteration_1/
    │   │   ├── proteinmpnn_retry/
    │   │   └── colabfold_retry/
    │   ├── fb2_summary.tsv
    │   └── logs/
    │
    ├── feedback3_loop/
    │   ├── blast_subjects.faa
    │   ├── enriched_hints.faa
    │   ├── braker_enriched/
    │   ├── proteins_enriched.faa
    │   ├── fb3_annotation_diff.txt
    │   └── logs/
    │
    ├── feedback4_loop/
    │   ├── domesticated_seqs/
    │   ├── colabfold_dom/
    │   ├── fb4_report.tsv
    │   └── logs/
    │
    ├── feedback5_loop/
    │   ├── designed_hints.faa
    │   ├── fb5_manifest.tsv
    │   └── logs/
    │
    ├── feedback6_loop/
    │   ├── feedback6_taxonomy_audit.txt
    │   ├── braker_corrected/
    │   ├── proteins_corrected.faa
    │   └── logs/
    │
    ├── logs/
    │   ├── pipeline.log               # Master log (all steps mirrored here)
    │   ├── repeatmodeler_detail.log   # Full unfiltered RepeatModeler output
    │   ├── step1_genome_fetch.log
    │   ├── step2b_repeatmask.log
    │   ├── step3_annotate.log
    │   └── ...
    ├── README.md                      # Auto-generated run report
    └── genomopipe_summary.md          # Orchestrator phase summary
```

---

## Tips and Notes

**Orchestrator vs. standalone scripts:** The orchestrator is the right tool for full runs and automated workflows. Run scripts directly when iterating on a single phase during development or debugging.

**Config layering:** The priority chain is always `defaults → config file → CLI flags`. Keep a base config for an organism and override specific keys per experiment on the command line without editing the file.

**Feedback loop ordering:** Run FB2 before FB1 - converging sequences first makes better motifs. Run FB3 and FB6 independently (they address different aspects of annotation quality). Run FB4 after Phase 2 plasmid design. Run FB5 last, once designs are fully validated.

**FB6 dry run first:** Use `fb6_dry_run: true` on the first run with a new organism to read the taxonomy audit before committing to a full BRAKER re-run, which takes as long as the original annotation.

**FB5 cross-species bootstrapping:** If running Genomopipe on multiple related organisms, use FB5's `target_organism` key to automatically seed the second organism's annotation with high-confidence protein designs from the first run.

**RepeatModeler runtime:** On genomes >500 Mbp, RepeatModeler typically runs 15–30 hours across 5–6 rounds before producing a repeat library. This is expected and normal. The pipeline does not hang during this period - check `logs/pipeline.log` for `[REPEATMASK]` progress markers or `logs/repeatmodeler_detail.log` for round-level detail. Use the GUI's live sidebar refresh to watch sentinel files appearing.

**Motif scaffolding from the start:** Pre-run ColabFold on `proteins.faa` and place results in `colabfold_pre/` before launching the pipeline to use a known structure as a motif from the very first RFdiffusion run. FB1 automates this iteratively for subsequent cycles.

**RNA-Seq for annotation quality:** `--auto_rnaseq` or a pre-existing BAM substantially improves BRAKER annotation for eukaryotes. For well-studied organisms, RNA-Seq is almost always available in SRA.

**OrthoDB proteins:** Downloaded automatically to `~/orthodb_<group>.fasta` on first use and reused on subsequent runs. Manual download: https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/

**Parallelism:** All parallelizable steps use GNU `parallel` if available, falling back to sequential loops. `conda install -c conda-forge parallel` gives significant speedups on multi-core systems.

**DNA Chisel:** `pip install dnachisel[reports]` substantially improves CDS domestication quality in Step 9. Without it the built-in heuristic fallback is used.

**Codon optimization host selection:** Use `e_coli` for bacterial expression, `s_cerevisiae` or `p_pastoris` for yeast secretion, `h_sapiens` for CHO/HEK293 mammalian expression. When in doubt, use `auto` method — DNA Chisel will simultaneously optimise codons and remove restriction sites in a single pass if it is installed.

**python-codon-tables for exotic hosts:** If your expression organism is not one of the five built-in hosts, install `pip install python-codon-tables` and set `expression_host` to the full NCBI species name (e.g. `Corynebacterium glutamicum`). The library fetches the correct Kazusa-derived codon table automatically.

**Codon optimization and FB4:** Always run Feedback Loop 4 after Phase 2 when `codon_optimize: true`. FB4 re-folds the optimized sequences with ColabFold and flags any designs where synonymous changes caused a pLDDT drop, giving you an early warning before synthesis.

**Detecting externally-launched runs in the GUI:** The app detects pipelines launched from the terminal. Load the output directory via Browse Results or hit ↻ on the run sidebar - the app detects the active run and begins live-refreshing automatically.

---

## Troubleshooting

**Script not found on startup** the orchestrator validates all script paths before starting. Set `--scripts_dir` to the directory containing all `.sh` and `.py` pipeline files, or place everything alongside `genomopipe.py`.

**`gffread not found`**
```bash
mamba install -n braker_env -c bioconda gffread
```

**BRAKER fails with "less than 1000 introns"** the pipeline automatically retries in `--esmode`. For fragmented genomes, supplement with `--auto_rnaseq`.

**RepeatModeler fails on small genomes** expected - masking is skipped automatically for genomes under 1 Mbp.

**`consensi.fa` is 0 bytes / RepeatMasker never ran** if the pipeline is still running, this is normal - `consensi.fa` is only written after all RepeatModeler rounds complete, which can take many hours. Do not interrupt. Check `logs/repeatmodeler_detail.log` for round-level progress. If the pipeline has already exited and `consensi.fa` is still empty, RepeatModeler failed silently; check `logs/step2b_repeatmask.log` for error detail.

**ColabFold GPU errors** remove `--use-gpu-relax` from the `colabfold_batch` call in Step 6 for CPU-only runs.

**GeneMark license errors** ensure `.gm_key` is present in `$HOME`. Download from: https://genemark.bme.gatech.edu/license_download.cgi

**SRA download failures** `prefetch` requires network access and may fail behind some firewalls. Provide a BAM manually with `--bam`.

**Part not found in Phase 2** `plasmid_design_moclo_v3.py` reports exactly which directories were searched. Set `MOCLO_PARTS_DIR`, add a `parts_dir` key to the config, or pass a full path in the `genes` list.

**Domestication loop does not converge** install `dnachisel[reports]` - its constraint-satisfaction solver resolves most cases the heuristic fallback cannot handle.

**FB2 designs always unconverged** lower `plddt_pass`, increase `resample_n`, or raise `resample_temp`. Some backbones are genuinely difficult to thread - consider running FB1 to redesign those scaffolds before retrying FB2.

**FB3 / FB6 BRAKER re-run fails** both loops require the same environment variables as `genome_to_design.sh`. Verify `genemark_path` is correct and that the Augustus config directory is writable. Both loops copy it to `~/augustus_config_writable` automatically as a fallback.

**FB4 flags many domesticated sequences** a large number of FAIL flags suggests the domestication heuristic is clustering changes in a small region. Installing `dnachisel[reports]` distributes synonymous substitutions more evenly across the CDS.

**FB5 target organism run fails immediately** verify that `genome_to_design.sh` is on PATH or in `scripts_dir` and that all Conda environments are set up. The target run is a full Phase 1 invocation with all the same dependencies.

**FB6 always suggests the same partition** if `original_partition` is unknown (BRAKER log not found or not parseable), every partition keyword passes the mismatch check. Verify that `logs/step3_annotate.log` exists and contains the OrthoDB filename pattern. Supply `run_dir` explicitly with `--run_dir` if the log lives elsewhere.

**Phase 2 fails but pipeline continues** Phase 2 is non-critical in the orchestrator - a failure is logged and feedback loops still run. Check `genomopipe_phase2.log` in the output directory for the error.

**3DMol not loading in the GUI** the structure viewer requires `3Dmol-min.js` in the `bioforge_app/` directory. Download from https://3dmol.org or install via npm. Without it, a placeholder is shown but all other panels function normally.

**GUI sidebar shows stale `.step*.failed` file during a `--force` run** `--force` clears `.failed` sentinels at pipeline startup. If the sidebar has not yet refreshed, it may briefly show the old state. Hit ↻ or wait up to 4 seconds for the auto-refresh to fire.

**Prokka fails: `Can't locate XML/Simple.pm`** the `XML::Simple` Perl module is missing from `braker_env`:
```bash
conda install -n braker_env -c bioconda perl-xml-simple
```

**Prokka fails: `needs blastp 2.2 or higher`** Prokka 1.13 has a version-parsing bug - it strips the last digit from the version string, so `2.15` becomes `2.1` which fails the `>= 2.2` check. Two fixes are needed:

1. Ensure `braker_env`'s blastp symlinks to the base environment's copy (which Prokka parses correctly):
```bash
ln -sf ~/miniconda3/bin/blastp ~/miniconda3/envs/braker_env/bin/blastp
```

2. Lower Prokka's minimum required blastp version to `2.1` so the parsed value passes:
```bash
sudo sed -i 's/MINVER  => "2\.2"/MINVER  => "2.1"/g' ~/miniconda3/bin/prokka
```

Verify both changes took effect:
```bash
grep -n "MINVER" ~/miniconda3/bin/prokka | grep -A1 "159:\|165:"
ls -la ~/miniconda3/envs/braker_env/bin/blastp
```

---

## Citations

If you use Genomopipe in published work, please cite:

- **Genomopipe**: Ditzler, C.J. (2026). Genomopipe: An automated pipeline from genome annotation to diffusion-based protein structure design (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18934145

```bibtex
@software{ditzler2026genomopipe,
  author  = {Ditzler, Cole James},
  title   = {Genomopipe: An automated pipeline from genome annotation
             to diffusion-based protein structure design},
  year    = {2026},
  version = {v1.0.0},
  doi     = {10.5281/zenodo.18934145},
  url     = {https://doi.org/10.5281/zenodo.18934145}
}
```

Please also cite the underlying tools:

- **BRAKER**: Brůna et al. (2021). *NAR Genomics and Bioinformatics*, 3(1). https://doi.org/10.1093/nargab/lqaa108
- **Augustus**: Stanke et al. (2008). *Nucleic Acids Research*, 36(Web Server issue). https://doi.org/10.1093/nar/gkn220
- **GeneMark-ETP**: Brůna et al. (2023). *bioRxiv*. https://doi.org/10.1101/2023.01.13.524024
- **Prokka**: Seemann (2014). *Bioinformatics*, 30(14). https://doi.org/10.1093/bioinformatics/btu153
- **RepeatModeler / RepeatMasker**: Smit & Hubley. RepeatMasker Open-4.0. http://www.repeatmasker.org
- **RFdiffusion**: Watson et al. (2023). *Nature*, 620. https://doi.org/10.1038/s41586-023-06415-8
- **ProteinMPNN**: Dauparas et al. (2022). *Science*, 378(6615). https://doi.org/10.1126/science.add2187
- **ColabFold**: Mirdita et al. (2022). *Nature Methods*, 19. https://doi.org/10.1038/s41592-022-01488-1
- **BLAST**: Camacho et al. (2009). *BMC Bioinformatics*, 10. https://doi.org/10.1186/1471-2105-10-421
- **STAR**: Dobin et al. (2013). *Bioinformatics*, 29(1). https://doi.org/10.1093/bioinformatics/bts635
- **Biopython**: Cock et al. (2009). *Bioinformatics*, 25(11). https://doi.org/10.1093/bioinformatics/btp163
- **DNA Chisel**: Zulkower & Rosser (2022). *Bioinformatics*, 38(16). https://doi.org/10.1093/bioinformatics/btac296
- **3Dmol.js**: Rego & Koes (2015). *Bioinformatics*, 31(8). https://doi.org/10.1093/bioinformatics/btu829
- **Electron**: https://www.electronjs.org

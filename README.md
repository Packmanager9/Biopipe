# Genomopipe

An end-to-end automated pipeline that takes an organism name and produces computationally designed protein structures taken from raw genome fetching through annotation, backbone design, sequence design, and structure prediction.

---

## Overview

```
Organism Name
     │
     ▼
 [Step 1]  genome_fetch     → Downloads reference genome from NCBI (RefSeq / GenBank)
     │
     ▼
 [Step 1b] Sanitize FASTA   → Normalizes sequence headers
     │
     ▼
 [Step 2]  BBTools QC       → Adapter trimming and quality control
     │
     ▼
 [Step 2b] RepeatMasker     → Masks repetitive elements (eukaryotes)
     │
     ▼
 [Step 3]  BRAKER / Prokka  → Gene annotation (eukaryote or prokaryote)
     │
     ▼
 [Step 4]  Protein Extract  → Filters longest isoforms, removes short sequences
     │
     ▼
 [Step 5]  RFdiffusion      → De novo protein backbone design
     │
     ▼
 [Step 5b] ProteinMPNN      → Sequence design for RFdiffusion backbones
     │
     ▼
 [Step 6]  ColabFold        → Structure prediction and relaxation
     │
     ▼
 [Step 7]  BLAST            → Remote validation of designed sequences
     │
     ▼
 [Step 8]  Rename           → Annotate outputs with putative function names
```

---

## Requirements

### Hardware
- Multi-core CPU (pipeline auto-detects and uses all available threads)
- GPU recommended for ColabFold and RFdiffusion (falls back to CPU)
- ~50–200 GB disk space depending on organism genome size

### Software and Conda Environments

Three Conda environments are required. They must be created before running the pipeline.

| Environment | Used for |
|---|---|
| `braker_env` | Genome annotation, QC, protein extraction, BLAST, utility steps |
| `SE3nv` | RFdiffusion backbone design and ProteinMPNN sequence design |
| `colabfold` | ColabFold structure prediction |

#### `braker_env` has key packages
```bash
conda create -n braker_env
conda activate braker_env
mamba install -c bioconda -c conda-forge \
    braker augustus repeatmodeler repeatmasker bbduk prokka \
    blast gffread biopython star samtools sra-tools \
    python perl
```

#### `SE3nv` possible RFdiffusion environment
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

### External Tools

| Tool | Default Path | Notes |
|---|---|---|
| GeneMark-ETP | `~/genemark-etp-full/gmetp_linux_64/bin` | License required from Georgia Tech |
| RFdiffusion | `~/RFdiffusion/` | Clone from GitHub |
| RFdiffusion models | `~/models/` | Download `Base_ckpt.pt` separately |
| ProteinMPNN | `~/ProteinMPNN/` | Clone from GitHub |

### `genome_fetch` dependencies
```bash
conda activate bioenv
mamba install -c bioconda ncbi-genome-download entrez-direct cd-hit biopython
```

---

## Installation

```bash
# 1. Clone or copy the scripts to your home directory
cp genome_fetch.sh ~/Desktop/favbooks/genome_fetch.sh
chmod +x ~/Desktop/favbooks/genome_fetch.sh
chmod +x genome_to_design.sh

# 2. Install Miniconda (if not already present)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/etc/profile.d/conda.sh

# 3. Set up Conda environments (see Requirements above)

# 4. Download GeneMark license key and place gmes binary at the default path
# https://genemark.bme.gatech.edu/license_download.cgi

# 5. Download RFdiffusion model weights
mkdir -p ~/models
# See: https://github.com/RosettaCommons/RFdiffusion for download links
```

---

## Usage

### Basic syntax

```bash
./genome_to_design.sh "Organism name" /path/to/output [eukaryote_flag] [options]
```

### Arguments

| Argument | Position | Required | Description |
|---|---|---|---|
| `organism` | 1 | Yes | Organism name (quoted) or NCBI TaxID |
| `output_dir` | 2 | Yes | Path to output directory (created if absent) |
| `is_eukaryote` | 3 | No | `true` for eukaryotes, `false` for prokaryotes (default: `false`) |
| `--GENEMARK_PATH=` | flag | No | Path to GeneMark bin dir (default: `~/genemark-etp-full/gmetp_linux_64/bin`) |
| `--bam=` | flag | No | Path to existing RNA-Seq BAM file, or comma-separated list, for BRAKER hints |
| `--auto_rnaseq` | flag | No | Auto-query and download RNA-Seq from SRA for BRAKER hints |
| `--force` | flag | No | Force a full re-run from scratch, ignoring all checkpoints |

### Examples

```bash
# Prokaryote call: minimal
./genome_to_design.sh "Escherichia coli" ./output

# Eukaryote will attempt to run with automatic RNA-Seq download for annotation hints
./genome_to_design.sh "Taraxacum officinale" ./output true --auto_rnaseq

# Eukaryote run with a pre-existing BAM file
./genome_to_design.sh "Arabidopsis thaliana" ./output true --bam=/data/rnaseq.bam

# Custom GeneMark path
./genome_to_design.sh "Homo sapiens" ./output true \
    --GENEMARK_PATH=/opt/genemark/bin --auto_rnaseq

# Force full re-run
./genome_to_design.sh "Mus musculus" ./output true --force
```

---

## Resumability

The pipeline uses checkpoint files (`.step*.done` / `.step*.failed`) to support resuming interrupted runs. Re-running the same command will skip all steps that already completed successfully.

```bash
# Resume an interrupted run, just re-run the original command:
./genome_to_design.sh "Taraxacum officinale" ./output true --auto_rnaseq

# Re-run a single specific step:
rm ./output/latest/.step4_protein_extract.done
./genome_to_design.sh "Taraxacum officinale" ./output true --auto_rnaseq

# Start completely fresh:
./genome_to_design.sh "Taraxacum officinale" ./output true --auto_rnaseq --force
```

Each run is stored in a timestamped subdirectory (`run_YYYYMMDD_HHMM/`) under the output directory, with a `latest/` symlink pointing to the most recent run.

---

## Output Structure

```
output/
├── latest -> run_20240315_1430/      # Symlink to most recent run
└── run_20240315_1430/
    ├── Taraxacum_officinale/         # Raw genome files (.fna)
    ├── clean.fna                     # QC-cleaned genome
    ├── masked.fna                    # Repeat-masked genome (eukaryotes)
    ├── braker_out/                   # BRAKER annotation outputs
    │   ├── braker.gtf
    │   ├── braker_fixed.gtf
    │   ├── Taraxacum_officinale.faa  # Annotated proteins (all isoforms)
    │   └── genome.fa
    ├── proteins.faa                  # Filtered proteins (longest isoform, >100aa)
    ├── designs/                      # RFdiffusion backbone PDBs
    ├── proteinmpnn_out/
    │   ├── seqs/                     # Multi-sequence .fa files from ProteinMPNN
    │   └── split_seqs/               # Single-sequence .fasta files for ColabFold
    ├── colabfold_out/                # Predicted/relaxed structures (.pdb, .json)
    ├── blast_results/                # Per-design BLAST outputs
    ├── blast_results.txt             # Combined BLAST summary
    ├── logs/
    │   ├── pipeline.log              # Master log (all steps)
    │   ├── step1_genome_fetch.log
    │   ├── step3_annotate.log
    │   └── ...                       # Per-step logs
    ├── README.md                     # Auto-generated run report
    └── .step*.done / .step*.failed   # Checkpoint sentinels
```

---

## Pipeline Steps In Detail

### Step 1. Genome Fetch (`genome_fetch.sh`)
Queries NCBI for the organism by name or TaxID, determines the relevant taxonomic group, and attempts to download the best available genome assembly in priority order:
1. RefSeq reference assemblies
2. Any RefSeq assembly
3. Any GenBank assembly
4. Falls back to fetching individual nucleotide sequences and clustering with `cd-hit-est`

### Step 1b. FASTA Header Sanitization
Replaces spaces and pipe characters in sequence headers with underscores to prevent downstream tool failures.

### Step 2. QC with BBTools
Runs `bbduk.sh` to trim adapter sequences. If this step fails or produces no output, the pipeline falls back to the raw (sanitized) genome.

### Step 2b. Repeat Masking
Runs `RepeatModeler` to build a de novo repeat library, then `RepeatMasker` to soft-mask the genome. Automatically skipped for genomes under 1 MB. Required for eukaryotes (`repeat_critical=true` when `is_eukaryote=true`).

### Step 3. Gene Annotation
- **Eukaryotes**: BRAKER (Augustus + GeneMark-ETP), optionally with RNA-Seq BAM hints (`--bam`) or automatically downloaded RNA-Seq (`--auto_rnaseq`). OrthoDB protein hints are automatically downloaded for the correct taxonomic group (fungi, viridiplantae, vertebrata, etc.).
- **Prokaryotes**: Prokka.

The `--auto_rnaseq` mode queries SRA via Entrez, downloads runs with `prefetch` + `fasterq-dump`, aligns with STAR, and merges BAMs with `samtools`.

### Step 4. Protein Extraction and Filtering
Parses the annotated `.faa`, keeps only the longest isoform per gene, and filters out proteins shorter than 100 aa (200 aa for eukaryotes). Outputs `proteins.faa`.

### Step 5. RFdiffusion Backbone Design
Generates novel protein backbone structures. If a `colabfold_pre/` directory exists with pre-computed structures, the top-ranked structure is used as a motif for scaffolding. Otherwise runs unconditional de novo design.

### Step 5b. ProteinMPNN Sequence Design
Designs amino acid sequences for each RFdiffusion backbone, producing 8 candidate sequences per backbone. Sequences are split into individual FASTA files for ColabFold input.

### Step 6. ColabFold Structure Prediction
Predicts and relaxes the structure of each ProteinMPNN-designed sequence using ColabFold (AlphaFold2 backend), with 3 recycle passes.

### Step 7. BLAST Validation
Runs remote `blastp` against the NCBI `nr` database for each designed sequence to identify putative function by homology.

### Step 8. Rename Outputs
Renames ColabFold output files to include the putative function name derived from BLAST results, replacing the generic design ID with a human-readable label.

---

## Tips and Notes

**Motif scaffolding with RFdiffusion:** To use a known structure as a scaffold motif, pre-run ColabFold on `proteins.faa` and place the results in `colabfold_pre/` before launching the pipeline. The script will automatically detect and use the top-ranked structure.

**RNA-Seq for annotation quality:** Using RNA-Seq evidence (either `--bam` or `--auto_rnaseq`) substantially improves BRAKER annotation quality for eukaryotes. For well-studied organisms, RNA-Seq data is usually available in SRA.

**OrthoDB proteins:** Downloaded automatically to `~/orthodb_<group>.fasta` on first use and reused on subsequent runs. Manual download is also possible from [https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/).

**Parallelism:** All parallelizable steps use GNU `parallel` if available, falling back to sequential loops. Install with `conda install -c conda-forge parallel` for significant speedups on multi-core systems.

---

## Citations

If you use Genomopipe in published work, please cite:

- **Genomopipe**: Ditzler, C.J. (2026). Genomopipe: An automated pipeline from genome annotation to diffusion-based protein structure design (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

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

Genomopipe also depends on the following tools so please cite them accordingly:

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

---

## Troubleshooting

**`gffread not found`**  to install into `braker_env`:
```bash
mamba install -n braker_env -c bioconda gffread
```

**BRAKER fails with "less than 1000 introns"** if this happens the pipeline automatically retries in `--esmode`. For very fragmented or small genomes, consider supplementing with RNA-Seq via `--auto_rnaseq`.

**RepeatModeler fails on small genomes** this is expected behavior; the pipeline skips masking automatically for genomes under 1 MB.

**ColabFold GPU errors** to fix: remove `--use-gpu-relax` from the `colabfold_batch` call in Step 6 if running on CPU only.

**GeneMark license errors** (key is a DNA sequence, and may not look normal) ensure the `.gm_key` license file is present in `$HOME`. Download from: https://genemark.bme.gatech.edu/license_download.cgi

**SRA download failures** using `prefetch` requires network access and may fail behind certain firewalls. Manually provide a BAM file via `--bam` as an alternative.

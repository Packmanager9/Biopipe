#!/bin/sh
# Usage: genome_fetch "organism_name_or_taxid"
# Example: genome_fetch "Escherichia coli"
# genome_fetch 562
# Assumes tools installed in Conda 'bioenv' (ncbi-genome-download, cd-hit-est, entrez-direct via conda install -c bioconda ncbi-genome-download cd-hit entrez-direct)
# Also assumes biopython for group determination: conda install -c conda-forge biopython
if [ -z "$1" ]; then
  echo "Error: Provide organism name or taxid as argument."
  exit 1
fi
organism_input="$1"
organism_query=$(echo "$organism_input" | sed 's/_/ /g')
output_dir=$(echo "$organism_input" | sed 's/ /_/g')
mkdir -p "$output_dir"
cd "$output_dir" || exit 1
log_file="$(pwd)/genome_fetch.log"
log() {
  echo "[$(date '+%F %T')] $*" | tee -a "$log_file"
}
log "Starting genome fetch for '$organism_input'"
# Load Conda (adjust path if Miniconda is in /opt)
. ~/miniconda3/etc/profile.d/conda.sh
conda activate bioenv || { log "Error: Conda env activation failed"; exit 1; }
log "Conda environment activated: bioenv"
# Get taxonomy ID if name provided
if expr "$organism_input" : '^[0-9]\+$' >/dev/null; then
  taxid="$organism_input"
else
  log "Searching for taxonomy ID for '$organism_query'"
  taxid=$(esearch -db taxonomy -query "$organism_query" | efetch -format uid)
  if [ -z "$taxid" ]; then
    log "Error: No taxonomy ID found for '$organism_query'."
    exit 1
  fi
fi
log "Fetching data for organism (TaxID: $taxid)..."
# Determine NCBI group based on taxonomy
log "Determining NCBI group based on taxonomy"
group=$(python3 - <<EOF
from Bio import Entrez
Entrez.email = "PricklyPearEnterprises@gmail.com"
try:
    handle = Entrez.efetch(db="taxonomy", id="$taxid", retmode="xml")
    record = Entrez.read(handle)
    lineage = [tax['ScientificName'].lower() for tax in record[0]['LineageEx']]
    if 'bacteria' in lineage:
        print('bacteria')
    elif 'archaea' in lineage:
        print('archaea')
    elif 'fungi' in lineage:
        print('fungi')
    elif 'viridiplantae' in lineage:
        print('plant')
    elif 'metazoa' in lineage:
        if 'vertebrata' in lineage:
            if 'mammalia' in lineage:
                print('vertebrate_mammalian')
            else:
                print('vertebrate_other')
        else:
            print('invertebrate')
    elif 'protozoa' in lineage:
        print('protozoa')
    elif 'viruses' in lineage:
        print('viral')
    else:
        print('all')
except Exception as e:
    print('all') # Fallback
EOF
)
log "Determined NCBI group: $group"
# Get descendant taxids (includes subspecies) to handle ambiguity
log "Fetching descendant taxids"
desc_taxids=$(esearch -db taxonomy -query "txid$taxid[Subtree]" | efetch -format uid | tr '\n' ',' | sed 's/,$//')
if [ -z "$desc_taxids" ]; then
  desc_taxids="$taxid"
fi
log "Including descendant taxids for broader search: $desc_taxids"
# Common options for ncbi-genome-download
common_opts="--formats fasta --output-folder . --verbose --retries 3 --flat-output --assembly-levels chromosome,complete"
# Construct taxid option
taxid_opt="--taxids $desc_taxids"
# Step 1: Try RefSeq reference genomes (any level)
log "Starting Step 1: RefSeq reference genomes (any level)..."
ncbi-genome-download $group $common_opts $taxid_opt --refseq-categories reference --section refseq 2>&1 | tee -a "$log_file"
if find . -name "*.fna.gz" -print -quit | grep -q .; then
  log "Unzipping downloaded files (recursive)"
  find . -name "*.fna.gz" -exec gunzip -f {} \;
  log "Reference genome (RefSeq) downloaded successfully."
  conda deactivate
  exit 0
fi
# Step 2: Fall back to any RefSeq genomes (any level)
log "Starting Step 2: Any RefSeq genomes..."
ncbi-genome-download $group $common_opts $taxid_opt --section refseq 2>&1 | tee -a "$log_file"
if find . -name "*.fna.gz" -print -quit | grep -q .; then
  log "Unzipping downloaded files (recursive)"
  find . -name "*.fna.gz" -exec gunzip -f {} \;
  log "Genome (RefSeq) downloaded successfully."
  conda deactivate
  exit 0
fi
# Step 3: Fall back to GenBank genomes (any level)
log "Starting Step 3: GenBank genomes..."
ncbi-genome-download $group $common_opts $taxid_opt --section genbank 2>&1 | tee -a "$log_file"
if find . -name "*.fna.gz" -print -quit | grep -q .; then
  log "Unzipping downloaded files (recursive)"
  find . -name "*.fna.gz" -exec gunzip -f {} \;
  log "Genome (GenBank) downloaded successfully."
  conda deactivate
  exit 0
fi
# Step 4: No assemblies - gather individual sequences from subtree and cluster
log "No assemblies found. Gathering available genomic sequences from subtree..."
# Construct query with OR for each txid[Organism]
txid_query=""
for dt in $(echo $desc_taxids | tr ',' ' '); do
  if [ -n "$txid_query" ]; then
    txid_query="$txid_query OR "
  fi
  txid_query="${txid_query}txid$dt[Organism]"
done
query="($txid_query) AND biomol_genomic[PROP] NOT partial[title]"
log "Constructed nucleotide query: $query"
esearch -db nucleotide -query "$query" | efetch -format fasta > all_sequences.fasta 2>&1 | tee -a "$log_file"
if [ ! -s all_sequences.fasta ]; then
  log "Error: No genomic sequences found for this organism or descendants."
  conda deactivate
  exit 1
fi
log "Clustering sequences with cd-hit-est"
cd-hit-est -i all_sequences.fasta -o representative_sequences.fasta -c 0.99 -n 10 -M 0 -T 0 2>&1 | tee -a "$log_file"
rm all_sequences.fasta
log "Parsimonious representative sequences gathered and clustered successfully in representative_sequences.fasta."
conda deactivate
log "Genome fetch completed"
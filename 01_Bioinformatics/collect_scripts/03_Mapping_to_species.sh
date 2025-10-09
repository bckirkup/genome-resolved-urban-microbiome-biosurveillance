#!/usr/bin/env bash
## Last update 03-09-2025 by Suleiman & AbdulAziz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err


# Author: Suleiman & Assistant
# Purpose: Merge ARG, VFDB, and GTDB abundance outputs with species info for global biosurveillance analysis
# Date: Updated 2025-09-12


# ==============================
#        CONFIGURATION
# ==============================
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Environmental_metagenomics_global_derep"
MERGED_DIR="$BASE_DIR/merged_outputs"
THREADS=56

# Input files
GTDB_IN="$MERGED_DIR/gtdb_combined.tsv"
ARG_IN="$MERGED_DIR/arg_combined.tsv"
VFDB_IN="$MERGED_DIR/vfdb_combined.tsv"
GTDB_MAP_BAC="$MERGED_DIR/gtdb_genome_species.bac120_map.tsv"
GTDB_MAP_ARCH="$MERGED_DIR/gtdb_genome_species.ar53_map.tsv"

# Outputs
FINAL_OUT="$MERGED_DIR/final_outputs"
mkdir -p "$FINAL_OUT"

# ==============================
# Step 1: Build unified genome-to-species map
# ==============================
echo "[Step 1] Merging GTDB species maps..."
awk 'FNR==1 && NR!=1 {next} {print}' "$GTDB_MAP_BAC" "$GTDB_MAP_ARCH" > "$FINAL_OUT/genome_to_species.tsv"
echo "   Genome-to-species map created."

# ==============================
# Step 2: Add Species to GTDB Abundance
# ==============================
echo "[Step 2] Annotating GTDB with species info..."
awk -v map="$FINAL_OUT/genome_to_species.tsv" '
BEGIN {
  FS=OFS="\t"
  while ((getline < map) > 0) {
    species[$1] = $2
  }
}
NR==1 { print $0, "Species"; next }
{
  sid = $4  # Genome column
  print $0, (sid in species ? species[sid] : "Unmapped")
}' "$GTDB_IN" > "$FINAL_OUT/gtdb_with_species.tsv"
echo "   GTDB abundance annotated."

# ==============================
# Step 3: Extract Genome ID from ARG/VFDB headers
# ==============================
extract_genome_from_arg_vfdb() {
  local infile=$1
  local outfile=$2
  local type=$3

  awk -v type="$type" 'BEGIN{FS=OFS="\t"}
  NR==1 { print $0, "Genome"; next }
  {
    genome = "Unmapped"
    if (type == "arg") {
      match($3, /sample=([^|]+)/, m)
      genome = m[1]
    } else if (type == "vfdb") {
      split($3, a, "|")
      genome = a[1]
    }
    print $0, genome
  }' "$infile" > "$outfile"
}

echo "[Step 3] Extracting genome ID from ARG and VFDB..."
extract_genome_from_arg_vfdb "$ARG_IN" "$FINAL_OUT/arg_with_genome.tsv" "arg"
extract_genome_from_arg_vfdb "$VFDB_IN" "$FINAL_OUT/vfdb_with_genome.tsv" "vfdb"
echo "   Genome extraction done."

# ==============================
# Step 4: Annotate ARG and VFDB with species info
# ==============================
annotate_species() {
  local infile=$1
  local outfile=$2

  awk -v map="$FINAL_OUT/genome_to_species.tsv" 'BEGIN {
    FS=OFS="\t"
    while ((getline < map) > 0) {
      species[$1] = $2
    }
  }
  NR==1 { print $0, "Species"; next }
  {
    print $0, ($NF in species ? species[$NF] : "Unmapped")
  }' "$infile" > "$outfile"
}

echo "[Step 4] Annotating ARG and VFDB..."
annotate_species "$FINAL_OUT/arg_with_genome.tsv" "$FINAL_OUT/arg_with_species.tsv"
annotate_species "$FINAL_OUT/vfdb_with_genome.tsv" "$FINAL_OUT/vfdb_with_species.tsv"
echo "   Annotation complete."

# ==============================
# Done
# ==============================
echo "
 All outputs ready in: $FINAL_OUT"
echo "  - GTDB: gtdb_with_species.tsv"
echo "  - ARG : arg_with_species.tsv"
echo "  - VFDB: vfdb_with_species.tsv"

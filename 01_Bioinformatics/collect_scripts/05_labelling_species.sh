#!/usr/bin/env bash

# Author: Suleiman & Assistant
# Purpose: Add enhanced taxonomic labels, skipping empty species/genus/etc.
# Date: 2025-09-12

# ==============================
#        CONFIGURATION
# ==============================
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Environmental_metagenomics_global_derep"
MERGED_DIR="$BASE_DIR/merged_outputs/final_outputs"
INPUT_FILE="$MERGED_DIR/gtdb_species_ARG_VFDB.tsv"
OUTPUT_FILE="$MERGED_DIR/gtdb_species_ARG_VFDB_labeled.tsv"
DELIM=$'\t'

# ==============================
#         PROCESSING
# ==============================
awk -F"$DELIM" -v OFS="$DELIM" '
BEGIN {
  print "Genome", "Sample", "Read", "Bin_ID", "RPKM", "TPM", "Covered Fraction", "Relative Abundance (%)", "Trimmed Mean", "Mean", "Variance", "Taxon", "Enhanced_Taxon"
}
NR > 1 {
  tax = $12
  enhanced = "Unknown"

  n = split(tax, parts, ";")

  for (i = n; i >= 1; i--) {
    if (parts[i] ~ /^s__/) {
      name = parts[i]; gsub(/^s__/, "", name)
      if (name != "") { enhanced = name; break }
    } else if (parts[i] ~ /^g__/) {
      name = parts[i]; gsub(/^g__/, "", name)
      if (name != "") { enhanced = name " spp."; break }
    } else if (parts[i] ~ /^f__/) {
      name = parts[i]; gsub(/^f__/, "", name)
      if (name != "") { enhanced = name; break }
    } else if (parts[i] ~ /^o__/) {
      name = parts[i]; gsub(/^o__/, "", name)
      if (name != "") { enhanced = name; break }
    } else if (parts[i] ~ /^c__/) {
      name = parts[i]; gsub(/^c__/, "", name)
      if (name != "") { enhanced = name; break }
    } else if (parts[i] ~ /^p__/) {
      name = parts[i]; gsub(/^p__/, "", name)
      if (name != "") { enhanced = name; break }
    } else if (parts[i] ~ /^d__/) {
      name = parts[i]; gsub(/^d__/, "", name)
      if (name != "") { enhanced = name; break }
    }
  }

  print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, enhanced
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo " Enhanced taxon labels written to: $OUTPUT_FILE"

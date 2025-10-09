#!/bin/bash

# ==============================
#        CONFIGURATION
# ==============================
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Environmental_metagenomics_global_derep"
MERGED_DIR="$BASE_DIR/merged_outputs/final_outputs"
THREADS=56

# ==== Input files ====
ARG_GENOME="$MERGED_DIR/arg_with_genome.tsv"
VF_GENOME="$MERGED_DIR/vfdb_with_genome.tsv"
GTDB_FILE="$MERGED_DIR/gtdb_with_species.tsv"

# ==== Output files ====
ALL_GENOMES="$MERGED_DIR/all_ARG_VFDB_genomes.txt"
GTDB_MATCHED="$MERGED_DIR/gtdb_species_ARG_VFDB.tsv"

# ==== Step 1: Extract genome column (column 7) ====
cut -f7 "$ARG_GENOME" | grep -v "^Genome" > "$MERGED_DIR/tmp_arg_genomes.txt"
cut -f7 "$VF_GENOME"  | grep -v "^Genome" > "$MERGED_DIR/tmp_vfdb_genomes.txt"

# ==== Step 2: Combine and deduplicate ====
cat "$MERGED_DIR"/tmp_*_genomes.txt | sort -u > "$ALL_GENOMES"
rm "$MERGED_DIR"/tmp_*_genomes.txt  # cleanup

# ==== Step 3: Filter GTDB (match on 4th column) ====
awk -v list="$ALL_GENOMES" '
BEGIN {
    FS = OFS = "\t"
    while ((getline < list) > 0) {
        ids[$1]
    }
}
NR == 1 || ($4 in ids)
' "$GTDB_FILE" > "$GTDB_MATCHED"

echo " Subset of GTDB with both ARG and VFDB genomes saved to:"
echo "   $GTDB_MATCHED"

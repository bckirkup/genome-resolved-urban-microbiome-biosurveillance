#!/usr/bin/env bash
## Last update 03-09-2025 by Suleiman & AbdulAziz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err


# ---- CONFIGURATION ----
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Environmental_metagenomics_global_derep"
THREADS=54

# Output paths
MERGED_DIR="$BASE_DIR/merged_outputs"
mkdir -p "$MERGED_DIR"

# Subdirectories containing cleaned abundance data
declare -A INPUT_PATHS=(
  [gtdb]="$BASE_DIR/taxonomy/coverm1"
  [arg]="$BASE_DIR/arg_abundance/abundance1"
  [vfdb]="$BASE_DIR/VFDB_OUT/abundance1"
)

# Subfolder names for each environment
ENV_LIST=(ambulance hosp_env hosp_sewage public_transport)

# 1. Merge cleaned abundance files for each type
for TYPE in gtdb arg vfdb; do
  echo "Merging cleaned $TYPE abundance files..."
  OUTFILE="$MERGED_DIR/${TYPE}_combined.tsv"
  HEADER_WRITTEN=false

  for ENV in "${ENV_LIST[@]}"; do
    CLEAN_DIR="${INPUT_PATHS[$TYPE]}/${ENV}_cleaned"
    for f in "$CLEAN_DIR"/*.tsv; do
      if [[ "$HEADER_WRITTEN" == false ]]; then
        head -n1 "$f" > "$OUTFILE"
        HEADER_WRITTEN=true
      fi
      tail -n +2 "$f" >> "$OUTFILE"
    done
  done

  echo "   Combined file created: $OUTFILE"
done

# 2. Create mapping between GTDB Genome IDs and Sample (for downstream join)
# Assumes gtdb summary is available
GTDB_SUMMARY="$BASE_DIR/taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"
GTDB_MAPPING="$MERGED_DIR/gtdb_genome_species.ar53_map.tsv"

awk -F '\t' 'NR>1 {print $1"\t"$2}' "$GTDB_SUMMARY" > "$GTDB_MAPPING"
echo "   Genome-to-species mapping written to: $GTDB_MAPPING"


GTDB_SUMMARY1="$BASE_DIR/taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
GTDB_MAPPING="$MERGED_DIR/gtdb_genome_species.bac120_map.tsv"

awk -F '\t' 'NR>1 {print $1"\t"$2}' "$GTDB_SUMMARY1" > "$GTDB_MAPPING"
echo "   Genome-to-species mapping written to: $GTDB_MAPPING"


echo -e "\n< All merged abundance tables and species mapping ready at: $MERGED_DIR"

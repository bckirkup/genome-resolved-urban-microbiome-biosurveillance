#!/usr/bin/env bash
## Last update 03-09-2025 by Suleiman & AbdulAziz
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err


# ---- Environment Setup ----
eval "$(conda shell.bash hook)" || true

# ---- Input Directories (update paths as needed) ----


# Base path
BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/Environmental_metagenomics_global_derep"

# Tools and environments
TOOLS=("taxonomy/coverm1" "arg_abundance/abundance1" "VFDB_OUT/abundance1")
TOOL_LABELS=("GTDB" "ARG" "VFDB")
ENVIRONMENTS=("ambulance" "hosp_env" "hosp_sewage" "public_transport")

# Process each tool directory
for i in "${!TOOLS[@]}"; do
  TOOL="${TOOLS[$i]}"
  LABEL="${TOOL_LABELS[$i]}"

  echo "[Cleaning $LABEL abundance]"
  for ENV in "${ENVIRONMENTS[@]}"; do
    INPUT_DIR="$BASE_DIR/$TOOL/$ENV"
    OUTPUT_DIR="$BASE_DIR/$TOOL/${ENV}_cleaned"
    mkdir -p "$OUTPUT_DIR"

    for f in "$INPUT_DIR"/*.tsv; do
      [[ -e "$f" ]] || continue  # skip if no files
      sample=$(basename "$f" .tsv | sed 's/_\(arg\|VFDB\)_abundance//')
      outfile="$OUTPUT_DIR/${sample}_${LABEL}_cleaned.tsv"

      # Keep header
      head -n 1 "$f" > "$outfile"

      # Filter rows: keep if any numeric column is nonzero
      awk -v sample="$sample" -v env="$ENV" 'NR>1 {
        keep = 0
        for (i = 2; i <= NF; i++) {
          if ($i ~ /^[0-9.]+$/ && $i+0 != 0) {
            keep = 1
            break
          }
        }
        if (keep) {
          print sample "\t" env "\t" $0
        }
      }' "$f" >> "$outfile"
    done
  done
done

echo " All abundances cleaned and saved under *_cleaned folders."

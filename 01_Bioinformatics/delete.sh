#!/bin/bash

set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <Project> [--run]"
    exit 1
fi

PROJECT="$1"
RUN=0
if [[ "${2:-}" == "--run" ]]; then
    RUN=1
fi

BASE_DIR="/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal/DATA/$PROJECT"
SRA_IDS="$BASE_DIR/SRA_ids"

if [ ! -f "$SRA_IDS" ]; then
    echo "SRA_ids file not found at $SRA_IDS"
    exit 1
fi
## Write the files you want to keep here. Any other thing is deleted
KEEP=(
    "assembly_qc"
    "bins"
    "checkm2_\$id"
    "fastqscreen_\$id"
    "megahit"
    "nonpareil_\$id"
    "quast"
    "\$id_R1.decontam.paired.fq.gz"
    "\$id_1.fastq.gz"
    "\$id_2.fastq.gz"
    "\$id_R1.qc.nophix.tagged_screen.txt"
    "\$id_R2.decontam.paired.fq.gz"
    "\$id_R2.qc.nophix.tagged_screen.txt"
    "\$id_singletons.decontam.fq.gz"
    "bbduk_phix_\$id.txt"
    "\$id_R1.qc.nophix.tagged_screen.html"
    "\$id_R2.qc.nophix.tagged_screen.html"
    "fastqc_posthost"
    "fastqc_prehost"
    "fastqscreen_\$id"
    "nonpareil_\$id"
)

echo "BASE: $BASE_DIR"
if [ "$RUN" -eq 1 ]; then
    echo "!!! ACTUAL DELETION ENABLED !!!"
else
    echo "Running in DRY RUN mode (no deletion will happen)"
fi

while read -r id; do
    [[ -z "$id" ]] && continue   # skip blank lines
    SAMPLE_DIR="$BASE_DIR/$id"
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "Warning: Directory $SAMPLE_DIR not found, skipping"
        continue
    fi

    # Prepare keep list for this sample
    KEEP_PATTERNS=()
    for pattern in "${KEEP[@]}"; do
        KEEP_PATTERNS+=( "${pattern//\$id/$id}" )
    done

    # Only operate inside sample dir!
    for item in "$SAMPLE_DIR"/*; do
        [ -e "$item" ] || continue
        basename_item="$(basename "$item")"
        keep=0
        for kp in "${KEEP_PATTERNS[@]}"; do
            if [[ "$basename_item" == "$kp" ]]; then
                keep=1
                break
            fi
        done
        if [ "$keep" -eq 0 ]; then
            echo "To delete: $item"
            if [ "$RUN" -eq 1 ]; then
                if [ -d "$item" ]; then
                    rm -rf "$item"
                else
                    rm -f "$item"
                fi
            fi
        fi
    done

done < "$SRA_IDS"

if [ "$RUN" -eq 0 ]; then
    echo
    echo "Dry run complete. To actually delete, run: bash $0 $PROJECT --run"
fi

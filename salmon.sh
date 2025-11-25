#!/bin/bash

# Set number of threads
THREADS=64

# Base directory
BASE_DIR="/DKHB/users/timtd/transcriptomes"
TRANSCRIPTOME="${BASE_DIR}/transcripts.fasta"
SALMON_INDEX="${BASE_DIR}/salmon_index"

# Build index if it doesn’t already exist
if [ ! -d "$SALMON_INDEX" ]; then
    echo "Building Salmon index..."
    salmon index -t "$TRANSCRIPTOME" -i "$SALMON_INDEX" -p $THREADS
    echo "Index built at $SALMON_INDEX"
fi

# Sample reads (interleaved)
declare -A READ_FILES
READ_FILES=( 
    ["DIAVIR1"]="${BASE_DIR}/DIAVIR1_ZKRN250004453_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR2"]="${BASE_DIR}/DIAVIR2_ZKRN250004454_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR3"]="${BASE_DIR}/DIAVIR3_ZKRN250004455_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR4"]="${BASE_DIR}/DIAVIR4_ZKRN250004456_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR5"]="${BASE_DIR}/DIAVIR5_ZKRN250004457_1A_22MWC3LT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR6"]="${BASE_DIR}/DIAVIR6_ZKRN250004458_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR7"]="${BASE_DIR}/DIAVIR7_ZKRN250004459_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR8"]="${BASE_DIR}/DIAVIR8_ZKRN250004460_1A_22MWFGLT4_L3_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR9"]="${BASE_DIR}/DIAVIR9_ZKRN250004461_1A_22MWFGLT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR10"]="${BASE_DIR}/DIAVIR10_ZKRN250004462_1A_22MWFGLT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR11"]="${BASE_DIR}/DIAVIR11_ZKRN250004463_1A_22MWFGLT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR12"]="${BASE_DIR}/DIAVIR12_ZKRN250004464_1A_22MWFGLT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
    ["DIAVIR13"]="${BASE_DIR}/DIAVIR13_ZKRN250004465_1A_22MWFGLT4_L1_1__paired__trimmed_pairs_R1_fastq_gz_R1.fastq.gz"
)

# Quantify expression
for SAMPLE in "${!READ_FILES[@]}"; do
    READ_FILE="${READ_FILES[$SAMPLE]}"
    OUT_DIR="${BASE_DIR}/salmon_quant/${SAMPLE}"

    echo "Processing ${SAMPLE}..."

    if [ ! -f "$READ_FILE" ]; then
        echo "Warning: Read file not found: $READ_FILE. Skipping..."
        continue
    fi

    salmon quant -i "$SALMON_INDEX" \
                 -l A \
                 -r "$READ_FILE" \
                 -p $THREADS \
                 --validateMappings \
                 -o "$OUT_DIR"

    echo "Finished quantifying ${SAMPLE}, results in ${OUT_DIR}/quant.sf"
done

echo "All quantifications complete."

#!/bin/bash

DB_PATH="kraken-standard-db"
READ_DIR="transcriptomes"
OUT_DIR="kraken_out"
THREADS=128

mkdir -p $OUT_DIR

# Loop through interleaved fastq.gz files for DIAVIR samples only
for f in $READ_DIR/DIAVIR*.fastq.gz; do
    sample=$(basename "$f" | cut -d'_' -f1)  # Extract sample ID like DIAVIR1
    echo "Classifying $sample"

    kraken2 --db $DB_PATH \
            --threads $THREADS \
            --output $OUT_DIR/${sample}_output.txt \
            --report $OUT_DIR/${sample}_report.txt \
            --use-names \
            --gzip-compressed \
            $f
done


Rscript Kraken_process.R

#!/bin/bash
#ssh timtd@gretel.nib.si

cd /scratch/timtd/trimmed/assembled/

LINEAGE="stramenopiles_odb10"
MODE="transcriptome"

for samp in DIAVIR*_assembled; do
    fasta="$samp/transcripts.fasta"

    # If the uncompressed file doesn't exist, fall back to the .gz
    if [ ! -f "$fasta" ]; then
        fasta="$samp/transcripts.min200.fasta.gz"
    fi

    echo "Running BUSCO on $samp using $fasta"

    busco -i "$fasta" \
          -l $LINEAGE \
          -m $MODE \
          -o busco_$samp \
          -f \
		  -c 30
done
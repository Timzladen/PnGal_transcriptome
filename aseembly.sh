#!/bin/bash


# Loop to create directories for DIAVIR1 to DIAVIR13
for i in {1..13}; do
    # Create the directory with the desired format
    mkdir -p "DIAVIR${i}_assembled"
    echo "Directory DIAVIR${i}_assembled created."
done

echo "Starting processing DIAVIR10"
spades.py --rna --pe-12 'DIAVIR10_ZKRN250004462-1A_22MWFGLT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR10_assembled/ --threads 64 --memory 250  
echo "Starting processing DIAVIR12"
spades.py --rna --pe-12 'DIAVIR12_ZKRN250004464-1A_22MWFGLT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR12_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR13"
spades.py --rna --pe-12 'DIAVIR13_ZKRN250004465-1A_22MWFGLT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR183_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR1"
spades.py --rna --pe-12 'DIAVIR1_ZKRN250004453-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR1_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR2"
spades.py --rna --pe-12 'DIAVIR2_ZKRN250004454-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR2_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR3"
spades.py --rna --pe-12 'DIAVIR3_ZKRN250004455-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR3_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR4"
spades.py --rna --pe-12 'DIAVIR4_ZKRN250004456-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR4_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR5"
spades.py --rna --pe-12 'DIAVIR5_ZKRN250004457-1A_22MWC3LT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR5_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR6"
spades.py --rna --pe-12 'DIAVIR6_ZKRN250004458-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR6_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR7"
spades.py --rna --pe-12 'DIAVIR7_ZKRN250004459-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR7_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR8"
spades.py --rna --pe-12 'DIAVIR8_ZKRN250004460-1A_22MWFGLT4_L3_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR8_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR9"
spades.py --rna --pe-12 'DIAVIR9_ZKRN250004461-1A_22MWFGLT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR9_assembled/ --threads 64 --memory 250
echo "Starting processing DIAVIR11"
spades.py --rna --pe-12 'DIAVIR11_ZKRN250004463-1A_22MWFGLT4_L1_1 (paired, trimmed pairs)R1.fastq.gz' -o assembled/DIAVIR11_assembled/ --threads 64 --memory 250




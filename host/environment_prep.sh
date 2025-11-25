#!/bin/bash
#deseq2 script Gretel

#Number of MAPPED reads tot he MMTESP transcriptomes in viral treatments is very low!


echo -e "Sample\tFragments\tApprox_Reads"

for d in salmon_quant/DIAVIR*; do
    sample=$(basename "$d")
    quant="$d/quant.sf"
    if [[ -f "$quant" ]]; then
        fragments=$(awk 'NR > 1 {s += $5} END {print s}' "$quant")
        reads=$(awk 'NR > 1 {s += $5} END {print s * 2}' "$quant")
        echo -e "$sample\t$fragments\t$reads"
    fi
done

#Sample  Fragments       Approx_Reads
#DIAVIR1 331510  663020
#DIAVIR10        4.26448e+06     8.52895e+06
#DIAVIR11        1.55807e+07     3.11614e+07
#DIAVIR12        3.41053e+06     6.82105e+06
#DIAVIR13        101904  203808
#DIAVIR2 349685  699370
#DIAVIR3 295997  591994
#DIAVIR4 1.31744e+07     2.63487e+07
#DIAVIR5 352115  704230
#DIAVIR6 1.16009e+07     2.32018e+07
#DIAVIR7 3.73756e+06     7.47512e+06
#DIAVIR8 267545  535090
#DIAVIR9 8.22175e+06     1.64435e+07


conda create -n deseq2_env -c conda-forge -c bioconda \
  r-base=4.2 \
  bioconductor-deseq2=1.38.0 \
  r-readr \
  mamba
  
 conda activate deseq2
 
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximport")

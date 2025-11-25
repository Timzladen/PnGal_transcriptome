#!/bin/bash

#run the assembly on CLC Genomics QCed and trimmed fastq files
bash assembly.sh

#cluster contigs
bash cluster.sh

#run salmon read quantification
bash salmon.sh

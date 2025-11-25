#!/bin/bash
#CLUSTERING
cd-hit -i merged_transcripts.fasta.transdecoder.pep -o clustered.pep -c 0.98 -n 5 -d 0 -T 8 -M 16000

#In R parse. clstr to create a transcript to cluster map
parse_cd_hit_clusters <- function(clstr_file) {
  lines <- readLines(clstr_file)
  
  cluster_id <- NULL
  results <- list()
  
  for (line in lines) {
    if (startsWith(line, ">Cluster")) {
      cluster_id <- sub(">Cluster ", "cluster_", line)
    } else {
      # Extract the transcript ID from quotes (after >)
      match <- regmatches(line, regexpr("(?<=>)[^ ]+", line, perl = TRUE))
      transcript_id <- match[1]
      results[[length(results) + 1]] <- data.frame(
        transcript_id = transcript_id,
        cluster = cluster_id,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine into a data frame
  bind_rows(results)
}

#### extract cluster representatives
# Extract representative (first) transcript per cluster
awk '/^>Cluster/{cluster=$0; next} !seen[cluster]++ {match($0, />[^ ]+/); print substr($0, RSTART+1, RLENGTH-1)}' clustered.pep.clstr > cluster_representatives.txt
# Remove the trailing '...' from the IDs
sed 's/\.p[0-9]\+\.*$//' cluster_representatives.txt > clean_cluster_representatives.txt

seqkit grep -f clean_cluster_representatives.txt merged_transcripts.fasta > clustered_transcripts.fasta


salmon index -t clustered_transcripts.fasta -i ../salmon_index_clust/clustered_index --threads 40
salmon quant -i clustered_index \
  -l A \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o sample_quant




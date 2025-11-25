library(tximport)
library(readr)
library(DESeq2)

# Load sample metadata (CSV with columns: sample,condition)
samples <- read.csv("sample_data.csv", stringsAsFactors = TRUE)
#rownames(samples) <- samples$sample
#samples$sample <- NULL  # optional: remove now-redundant colu

# Build file paths based on sample names
base_path <- "/users/timtd/transcriptomes/MMETSP/fasta/salmon_quant_renamed"
files <- file.path(base_path, samples$sample, "quant.sf")

txi <- tximport(files, type = "salmon", txOut = TRUE)

dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)


# Subset to samples of interest BEFORE running DESeq()
dds_subset <- dds[, samples$condition %in% c("virus5", "control5")]
dds_subset$condition <- droplevels(dds_subset$condition)

# Add biomass covariate
dds_subset$biomass <- log10(samples$biomass[samples$condition %in% c("virus5", "control5")])
design(dds_subset) <- ~ biomass + condition

# Now run DESeq only on the subset
dds_subset <- DESeq(dds_subset)
##NORMALIZE DATA PER DIATOM Cell

#add cell numbers as a covariate in the dds model
colData(dds)$diatom_scaled <- scale(log10(colData(dds)$diatom_cells))
design(dds) <- ~ diatom_scaled + condition
dds <- DESeq(dds)
F <- lfcShrink(dds, coef = "condition_virus5_vs_control5", type = "apeglm")



# Filter out low-count transcripts (min 10 reads in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]


# Run DESeq2
dds <- DESeq(dds)

# Get results (log fold changes, p-values)
res <- results(dds)

#Subset to Virus5 and Control5
keep_samples <- dds$condition %in% c("virus5", "control5")
dds_subset <- dds[, keep_samples]

# Drop unused levels
dds_subset$condition <- droplevels(dds_subset$condition)

# Relevel to make control5 the reference
dds_subset$condition <- relevel(dds_subset$condition, ref = "control5")

# Rerun DESeq2
dds_subset <- DESeq(dds_subset)

# Now shrink
resLFC <- lfcShrink(dds_subset, coef = "condition_virus5_vs_control5", type = "apeglm")


# Filter significant results
de_transcripts <- subset(resLFC, padj < 0.05 & abs(log2FoldChange) > 1)

# Save results to CSV
write.csv(as.data.frame(de_transcripts), "DESeq2_results_C5-V5.csv")

# Optional: print how many DE transcripts you got
cat("Found", nrow(de_transcripts), "significantly DE transcripts.\n")

#To get to FASTA
writeLines(rownames(de_transcripts), "de_transcripts_ids.txt")
quit()

#adds sequence data to the DE_transcript file
seqtk subseq ../assembled/merged_transcripts.fasta de_transcripts_ids.txt > DE_transcripts.fasta




library(DESeq2)
de <- read.csv("DESeq2_results_C5-V5.csv", row.names = 1)
colnames(de)<- c("X.query", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
# Load eggNOG results (skip comments)
eggnog <- read.delim("/scratch/timtd/transcriptomes/trimmed/salmon_quant/eggnog_full_named_clean.tsv", header = TRUE, sep = "\t", quote = "")

# Merge by transcript ID
de$X.query <- rownames(de)
merged <- merge(de, eggnog, by = "X.query")
colnames(merged)[1] <- "transcript_id"

# Save merged table
write.csv(merged, "DESeq2_with_annotations.csv")

#Quick barplot for top represented DE GO terms
cogs <- merged$COG_category

# Split multiple entries, flatten, and count
cog_split <- strsplit(as.character(cogs), "")
cog_table <- table(unlist(cog_split))
cog_df <- as.data.frame(cog_table)
colnames(cog_df) <- c("COG", "Count")

# Optional: Sort by count
cog_df <- cog_df[order(-cog_df$Count), ]

#Add COG descriptions
cog_desc <- c(
  A = "RNA processing & modification",
  B = "Chromatin structure & dynamics",
  C = "Energy production & conversion",
  D = "Cell cycle control & mitosis",
  E = "Amino acid metabolism & transport",
  F = "Nucleotide metabolism & transport",
  G = "Carbohydrate metabolism & transport",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  J = "Translation & ribosome structure",
  K = "Transcription",
  L = "DNA replication & repair",
  M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility",
  O = "Post-translational modification",
  P = "Inorganic ion transport & metabolism",
  Q = "Secondary metabolite biosynthesis",
  R = "General function prediction",
  S = "Function unknown",
  T = "Signal transduction",
  U = "Intracellular trafficking",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Y = "Nuclear structure",
  Z = "Cytoskeleton"
)

# Add description to the data frame
cog_df$Description <- cog_desc[cog_df$COG]

ggplot(head(cog_df, 15), aes(x = reorder(Description, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top COG Functional Categories (DE Transcripts)",
    x = "Function",
    y = "Number of DE Transcripts"
  )
  
#See what organisms the DE terms are associated with
# Extract NCBI taxon ID from seed_ortholog column
merged$taxid <- sub("\\..*", "", merged$seed_ortholog)

library(DBI)
library(RSQLite)

# Connect to the taxa database
db <- dbConnect(SQLite(), "/users/timtd/transcriptomes/trimmed/salmon_quant/eggnogDB/eggnog.taxa.db")
unique_taxids <- unique(merged$taxid)
# Build SQL-compatible string
taxid_str <- paste(unique_taxids, collapse = ",")


query <- paste0("SELECT taxid, spname FROM species WHERE taxid IN (", taxid_str, ")")
taxa_map <- dbGetQuery(db, query)

merged_named <- merge(merged, taxa_map, by = "taxid", all.x = TRUE)

###Only diatoms
# Match known diatom genera
diatom_hits <- grepl("Thalassiosira|Pseudo-nitzschia|Fragilariopsis|Skeletonema|Cyclotella|Chaetoceros", 
                     merged_named$spname, ignore.case = TRUE)

# Filter DE transcripts with diatom-like annotations
diatom_transcripts <- merged_named[diatom_hits, ]

# KEGG pathway column may have multiple entries separated by commas
kegg_terms <- strsplit(as.character(diatom_transcripts$KEGG_Pathway), ",")
kegg_flat <- unlist(kegg_terms)
kegg_counts <- sort(table(kegg_flat), decreasing = TRUE)

# Quick look at top KEGG pathways
head(kegg_counts, 20)

kegg_df <- as.data.frame(head(kegg_counts, 15))
colnames(kegg_df) <- c("KEGG_Pathway", "Count")

library(ggplot2)
ggplot(kegg_df, aes(x = reorder(KEGG_Pathway, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "darkcyan") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top KEGG Pathways (Diatom DE Transcripts)",
       x = "Pathway", y = "Number of Transcripts")

#Name KEGG pathways
pathways <- unique(unlist(strsplit(diatom_transcripts$KEGG_Pathway, ",")))
# Clean pathway list
pathways_clean <- pathways[grepl("^ko\\d{5}$", pathways)]  # Only keep properly formatted KEGG IDs
filtered_ids <- pathways_clean[!grepl("^ko01", pathways_clean)] #only keep pathways that are meaningful not global or general once that are not downloadable
#fetch kegg names
kegg_names <- sapply(filtered_ids, function(pid) {
  tryCatch({
    lines <- readLines(paste0("https://rest.kegg.jp/get/", pid), warn = FALSE)
    name_line <- grep("^NAME", lines, value = TRUE)
    sub("^NAME\\s+", "", name_line)
  }, error = function(e) {
    warning(paste("Skipping", pid, ":", conditionMessage(e)))
    NA
  })
})

# Make data frame for merging
kegg_name_df <- data.frame(
  KEGG_Pathway = names(kegg_names),
  KEGG_Name = unname(kegg_names)
)

kegg_counts_df <- as.data.frame(kegg_counts)
colnames(kegg_counts_df) <- c("KEGG_Pathway", "Count")
kegg_annotated <- merge(kegg_counts_df, kegg_name_df, by = "KEGG_Pathway", all.x = TRUE)

ggplot(head(kegg_annotated[order(-kegg_annotated$Count), ], 15),
       aes(x = reorder(KEGG_Name, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top KEGG Pathways (Diatom DE Transcripts)",
       x = "Pathway", y = "Transcript Count")
  
###Up-Down regulation 
diatom_transcripts$direction <- ifelse(diatom_transcripts$log2FoldChange > 0,
                                       "Up_in_virus", "Down_in_virus")
									   
library(tidyr)
library(dplyr)

# Expand multiple KEGG pathway values into separate rows
kegg_long <- diatom_transcripts %>%
  filter(!is.na(KEGG_Pathway)) %>%
  mutate(KEGG_Pathway = strsplit(as.character(KEGG_Pathway), ",")) %>%
  unnest(KEGG_Pathway)									   
kegg_long <- kegg_long %>% filter(grepl("^ko\\d{5}$", KEGG_Pathway))
kegg_counts <- kegg_long %>%
  group_by(KEGG_Pathway, direction) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
  
 kegg_annotated <- merge(kegg_counts, kegg_name_df, by = "KEGG_Pathway", all.x = TRUE)
 
 library(ggplot2)

top_kegg <- kegg_annotated %>%
  group_by(KEGG_Name) %>%
  summarise(total = sum(n)) %>%
  top_n(15, total) %>%
  pull(KEGG_Name)

ggplot(kegg_annotated %>% filter(KEGG_Name %in% top_kegg),
       aes(x = reorder(KEGG_Name, n), y = n, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top KEGG Pathways in Diatom DE Transcripts",
       x = "Pathway", y = "Number of Transcripts", fill = "Regulation") +
  scale_fill_manual(values = c("Up_in_virus" = "firebrick", "Down_in_virus" = "steelblue"))

###Annotate upregulated NA transcripts
unannotated_up <- diatom_transcripts[
  diatom_transcripts$Description == "-" &
  diatom_transcripts$direction == "Up_in_virus", ]

writeLines(unannotated_up$transcript_id, "unannotated_upregulated_ids.txt")

conda activate transdecoder
TransDecoder.LongOrfs -t unannotated_upregulated.fasta
TransDecoder.Predict -t unannotated_upregulated.fasta
conda activate diamond
diamond blastp -q unannotated_upregulated.fasta.transdecoder.pep -d /biodbs/diamond/db_nr_20241121/nr.dmnd -o UPV5_diamond_hits.tsv -e 1e-5 -k 1 -p 64

#In R
subject_ids <- unique(diamond_hits$subject_id)
writeLines(subject_ids, "subject_ids.txt")

conda activate blast
blastdbcmd -db /biodbs/blastdb/120.ncbi/nr -entry_batch subject_ids.txt -outfmt "%i %t" -out subject_descriptions.txt


##Summarize annotated transcripts
#read in diamond hits of unannotated transcripts
library(dplyr)
library(tidyr)

#source("Annotate_domains.R") #this runs the script to annotate "other" fields by the domains of the taxids if this was not done before

diamond_hits <- read.delim("diamond_top1000.tsv", header = FALSE,
                           col.names = c("transcript_id", "hit_id", "pident", "length", "evalue", "bitscore", "hit_description"))
diamond_hits <- diamond_hits %>%
  mutate(organism_type = case_when(
    grepl("RNA virus|picorna", hit_description, ignore.case = TRUE) ~ "PnGalRNAV",
    grepl("bacter|microbium|fluvi", hit_description, ignore.case = TRUE)     ~ "Bacteria",
    grepl("diatom|pseudo-nitzschia|thalassiosira|fragilariopsis|cyclotella|chaetoceros", hit_description, ignore.case = TRUE) ~ "Diatom",
	grepl("phage", hit_description, ignore.case = TRUE) ~ "Virus_diamond",
    TRUE ~ "Other"
  ))
diamond_labels <- setNames(diamond_hits$organism_type, diamond_hits$transcript_id)
# Join annotation

tpm_annotated <- tpm_long %>%
  left_join(eggnog_full_named %>% select(transcript_id, organism_type), by = "transcript_id")
  
tpm_annotated <- tpm_annotated %>%
  mutate(organism_type = ifelse(transcript_id %in% names(diamond_labels),
                                diamond_labels[transcript_id],
                                organism_type))


tpm_summary <- tpm_annotated %>%
  mutate(organism_type = ifelse(is.na(organism_type), "Unannotated", organism_type)) %>%
  group_by(sample, organism_type) %>%
  summarise(group_TPM = sum(TPM), .groups = "drop")

tpm_summary_wide <- tpm_summary %>%
  pivot_wider(names_from = organism_type, values_from = group_TPM, values_fill = 0)

tpm_summary_final <- tpm_summary_wide %>%
  mutate(total_TPM = rowSums(across(where(is.numeric))),
         annotated_TPM = total_TPM - Unannotated)
write.csv(tpm_summary_final, "annotated_reads_summary.csv")
tpm_summary_long <- tpm_summary %>%
  mutate(organism_type = factor(organism_type, levels = c("Diatom", "Bacteria", "Virus", "PnGalRNAV", "Other", "Unannotated", "Archaea")))

png("Transcript_comp_groups.png", width=6, height=5, res=400, units="in")
ggplot(tpm_summary_long, aes(x = sample, y = group_TPM, fill = organism_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Transcriptomic Composition by Supergroup",
       y = "Summed TPM", x = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


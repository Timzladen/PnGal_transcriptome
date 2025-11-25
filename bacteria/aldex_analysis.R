####Bacteria ALDEX
library(tximport)
library(ALDEx2)
library(tidyverse)

base_path <- "/scratch/timtd/transcriptomes/trimmed/salmon_quant_clst"
# Load sample metadata (CSV with columns: sample,condition)
samples <- read.csv("/scratch/timtd/transcriptomes/trimmed/salmon_quant/sample_data.csv", stringsAsFactors = TRUE)
files <- file.path(base_path, samples$sample, "quant.sf")

#rownames(samples) <- samples$sample
eggnog_full_named <- read.delim("/scratch/timtd/transcriptomes/trimmed/salmon_quant/eggnog_full_named_clean.tsv", header = TRUE, sep = "\t", quote = "")
# Extract genus from species names
library(stringr)
eggnog_full_named$genus <- word(eggnog_full_named$spname, 1)


txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "no")
raw_counts <- txi$counts
raw_counts_int <- round(raw_counts)
conditions <- samples$condition  # e.g., "control5" vs "virus5"

subset_samples <- which(conditions %in% c("control5", "virus5"))
raw_counts_sub <- raw_counts_int[, subset_samples]

# Subset the condition vector to match
conditions_sub <- droplevels(conditions[subset_samples])

# Run ALDEx2 on the subset
aldex_res <- aldex(raw_counts_sub, conditions_sub, test = "t", effect = TRUE, denom = "all")

# ---- 1. Merge ALDEx2 results with annotation ----
# Ensure `aldex_res` has transcript IDs as rownames
aldex_df <- as.data.frame(aldex_res)
aldex_df$transcript_id <- rownames(aldex_df)

# Merge with eggnog annotations (assuming it's named `eggnog_full_named`)
aldex_annotated <- aldex_df %>%
  left_join(eggnog_full_named, by = "transcript_id")

# ---- 2. Filter for bacterial transcripts ----
aldex_bacteria <- aldex_annotated %>%
  filter(organism_type == "Bacteria")

# ---- 3. Identify significantly differentially expressed transcripts ----
sig_bacteria <- aldex_bacteria %>%
  filter(we.eBH < 0.1) %>%  # FDR cutoff
  mutate(direction = ifelse(diff.btw > 0, "Up_in_Virus", "Down_in_Virus"))

# ---- 4. Quick summary: COG enrichment ----
cog_summary <- sig_bacteria %>%
  filter(!is.na(COG_category)) %>%
  count(COG_category, direction) %>%
  group_by(direction) %>%
  arrange(desc(n)) %>%
  slice_head(n = 10)  # Top 10 per direction
  
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
cog_summary$Description <- cog_desc[cog_summary$COG_category]

# ---- 5. Plot COG enrichment ----
ggplot(cog_summary, aes(x = reorder(Description, n), y = n, fill = direction)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Top COG Categories (Bacterial DE Genes)",
       x = "COG Category", y = "Number of DE Genes") +
  theme_minimal() +
  scale_fill_manual(values = c("Up_in_Virus" = "firebrick", "Down_in_Virus" = "steelblue"))

# ---- 6. KEGG enrichment ----
pathways <- unique(unlist(strsplit(sig_bacteria$KEGG_Pathway, ",")))
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
kegg_long <- sig_bacteria %>%
  filter(!is.na(KEGG_Pathway)) %>%
  mutate(KEGG_Pathway = strsplit(as.character(KEGG_Pathway), ",")) %>%
  unnest(KEGG_Pathway)									   
kegg_long <- kegg_long %>% filter(grepl("^ko\\d{5}$", KEGG_Pathway))
kegg_counts <- kegg_long %>%
  group_by(KEGG_Pathway, direction, genus) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
  
 kegg_annotated <- merge(kegg_counts, kegg_name_df, by = "KEGG_Pathway", all.x = TRUE)
 

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
  scale_fill_manual(values = c("Up_in_Virus" = "firebrick", "Down_in_Virus" = "steelblue"))


###Link function to bacteria

ggplot(kegg_annotated, aes(x = KEGG_Pathway, y = genus, size = n, color = direction)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up_in_Virus" = "firebrick", "Down_in_Virus" = "steelblue")) +
  theme_minimal() +
  labs(title = "Functionally Annotated Bacterial DE Genes by Genus",
       x = "KEGG Pathway", y = "Genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Filter to top KEGG names first
top_kegg <- kegg_annotated %>%
  group_by(KEGG_Name) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  top_n(30, total) %>%
  pull(KEGG_Name)

# Filter main annotated data
kegg_top_annotated <- kegg_annotated %>%
  filter(KEGG_Name %in% top_kegg)

# Plot using KEGG_Name for x-axis
ggplot(kegg_top_annotated, aes(x = KEGG_Name, y = genus, size = n, color = direction)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Up_in_Virus" = "firebrick", "Down_in_Virus" = "steelblue")) +
  theme_minimal() +
  labs(title = "Functionally Annotated Bacterial DE Genes by Genus (Top 30 Pathways)",
       x = "KEGG Pathway", y = "Genus", size = "Gene Count", color = "Regulation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
# Filter unannotated DE bacterial genes
unannotated_candidates <- sig_bacteria %>%
  filter(is.na(Description) | Description == "" | Description == "-") %>%
  filter(effect > 1, we.eBH < 0.1)

# View and save
write.csv(unannotated_candidates, "unannotated_candidate_transcripts.csv", row.names = FALSE)


##AFTER ALPHAFOLD2 and FOLDSEEK refine anotations
# Load libraries
library(tidyverse)
library(ggplot2)
library(readr)

new_anotation<- read.csv("sig_bacteria_clst.csv")

# Load refined annotations (if not already loaded)
# sig_bacteria <- read.csv("bacteria_sig_reano.csv")
# kegg_name_df <- read.csv("kegg_name_df.csv")  # Optional: if you already have KEGG ID → Name

library(tidyverse)
library(readr)


# Step 1: Expand KEGG_Pathway strings into individual IDs
df_long <- new_anotation %>%
  filter(!is.na(KEGG_Pathway)) %>%
  mutate(transcript_row = row_number()) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(KEGG_Pathway = str_trim(KEGG_Pathway))

# Step 2: Merge KEGG names
df_long <- df_long %>%
  left_join(kegg_annotated, by = "KEGG_Pathway")

# Step 3: Collapse back to one row per transcript with merged KEGG names
df_merged <- df_long %>%
  group_by(transcript_row) %>%
  summarise(
    across(-c(KEGG_Pathway, KEGG_Name), first),
    KEGG_Pathway = paste(unique(KEGG_Pathway), collapse = ","),
    KEGG_Name = if (any(!is.na(KEGG_Name))) paste(unique(na.omit(KEGG_Name)), collapse = ",") else NA_character_,
    .groups = "drop"
  )

library(tidyverse)
library(ggplot2)

# Step 5: Prioritize annotation: KEGG_Name > Description > refined_description > "Unannotated"
df_merged <- df_merged %>%
  filter(domain != "PnGalRNAV")
df_annotated <- df_merged %>%
  mutate(
    KEGG_Name_clean = na_if(KEGG_Name, ""),
    KEGG_Name_clean = na_if(KEGG_Name_clean, "-"),
    Description_clean = na_if(Description, ""),
    Description_clean = na_if(Description_clean, "-"),
    #Refined_clean = na_if(refined_description, ""),
    #Refined_clean = na_if(Refined_clean, "-")
  ) %>%
  mutate(
    # Use only the first KEGG name if multiple exist
    KEGG_Main = str_split(KEGG_Name_clean, ",") %>% map_chr(~ if (length(.x) > 0) .x[1] else NA_character_),
    
    # Prioritized annotation logic
    Full_Annotation = case_when(
      !is.na(KEGG_Main) ~ KEGG_Main,
      !is.na(Description_clean) ~ Description_clean,
      #!is.na(Refined_clean) ~ Refined_clean,
      TRUE ~ "Unannotated"
    ),
    
    # Truncate for plotting
    Annotation = str_trunc(Full_Annotation, 30)
  )
  
#Select top 30 most up- and down-regulated transcripts (by absolute effect)
df_top100 <- df_annotated %>%
  mutate(abs_effect = abs(effect)) %>%
  arrange(desc(abs_effect)) %>%
  slice_head(n = 128)

# Summarize effect size for top 30 most DE genes
annotation_summary <- df_top100 %>%
  group_by(genus.x, Annotation, direction.x) %>%
  summarise(
    Mean_Effect = mean(effect, na.rm = TRUE),
    Abs_Effect = mean(abs(effect), na.rm = TRUE),
    .groups = "drop"
  )
  
# Reorder Annotation levels alphabetically (case-insensitive)
library(forcats)

annotation_summary <- annotation_summary %>%
  mutate(
    Annotation = fct_relevel(Annotation, unique(Annotation)[order(tolower(unique(Annotation)))])
  )

# Step 7: Bubble plot: X = annotation, Y = genus, bubble size = number of genes
pdf("Top80_DE_bac_genus.pdf", width= 16, height=9)
ggplot(annotation_summary, aes(x = Annotation, y = genus.x, size = Abs_Effect, fill = direction.x)) +
  geom_point(alpha = 0.8, shape = 21, color = "black") +
  scale_size(range = c(3, 15)) +
  theme_minimal() +
  labs(
    title = "Mean Effect Size of Top 80 DE Genes by Genus and Annotation",
    x = "Annotation (truncated)",
    y = "Genus",
    size = "Abs. Effect",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )
dev.off()


###### ---------------------------------
####           CAZY processing
###### ---------------------------------
# Load libraries
library(tidyverse)

# ---------------------------
# 1. Load and parse overview.txt
# ---------------------------
overview <- read_tsv("/scratch/timtd/transcriptomes/trimmed/assembled/dbcan_output/overview.txt", show_col_types = FALSE)

overview_parsed <- overview %>%
  mutate(
    Raw_CAZy = coalesce(HMMER, DIAMOND, dbCAN_sub),
    CAZy_Family = str_extract(Raw_CAZy, "^[A-Z]+[0-9]+(_[0-9]+)?"),
    EC_Number = if_else(`EC#` == "-", NA_character_, `EC#`)
  ) %>%
  select(Gene_ID = `Gene ID`, CAZy_Family, EC_Number) %>%
  filter(!is.na(CAZy_Family))
# ---------------------------
# 2. Load eggNOG annotations (must contain: transcript ID, CAZy/KEGG/COG fields, taxonomy etc.)
# ---------------------------
# Clean transcript ID suffixes

overview_parsed <- overview_parsed %>%
  mutate(transcript_id = str_replace(Gene_ID, "\\.p[0-9]+$", ""))

eggnog_full_named <- eggnog_full_named %>%
  mutate(transcript_id = str_replace(transcript_id, "\\.p[0-9]+$", ""))

# Make sure Gene_ID matches (e.g. transcript ID or protein ID used in dbCAN)
cazy_annotated_full <- eggnog_full_named %>%
  left_join(overview_parsed, by = "transcript_id")
# ---------------------------
# 2. Load substrate mapping and join
# ---------------------------
substrate_map <- read_tsv("/users/timtd/dbCAN3_db/fam-substrate-mapping.tsv", show_col_types=FALSE)

#keep only the top-priority substrate per family
substrate_map_clean <- substrate_map %>%
  group_by(Family) %>%
  slice(1) %>%
  ungroup()
  
cazy_annotated_full <- cazy_annotated_full %>%
  left_join(substrate_map_clean, by = c("CAZy_Family" = "Family"))
#UPDATE for unannotated alphafolds
cazy_annotated_full <- cazy_annotated_full %>%
  mutate(CAZy = if_else(transcript_id == "NODE_10_length_16375_cov_4049.107042_g4_i1", "GH149", CAZy))  
write.csv(cazy_annotated_full, "CAZYnalysis/cazy_annotated_full.csv")

# ---------------------------
# 3. Add TPM values from txi object (assuming txi already loaded)
# ---------------------------
colnames(txi$abundance) <- samples$sample
# Melt txi TPM matrix into long format
library(tidyverse)
txi_tpm_long <- txi$abundance %>%
  as.data.frame() %>%
  rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = "Sample", values_to = "TPM")
  
  
txi_tpm_long <- txi_tpm_long %>%
  mutate(transcript_id = str_replace(transcript_id, "\\.p[0-9]+$", ""))  

cazy_annotated_full <- cazy_annotated_full %>%
  mutate(transcript_id_clean = str_replace(transcript_id, "\\.p[0-9]+$", ""))

# Join TPM info into annotation
cazy_annotated_full_expr <- cazy_annotated_full %>%
  left_join(txi_tpm_long, by = c("transcript_id_clean" = "transcript_id"))
  


# ---------------------------
# 5. Save full annotated table
# ---------------------------
write_csv(cazy_annotated_full_expr, "annotated_cazy_transcripts.csv")

# ---------------------------
# 6. Summarize by substrate
# ---------------------------
substrate_summary <- cazy_annotated_full_expr %>%
  group_by(Substrate_high_level) %>%
  summarise(
    N_transcripts = n(),
    N_families = n_distinct(CAZy_Family),
    N_EC = n_distinct(EC_Number.x, na.rm = TRUE)
  ) %>%
  arrange(desc(N_transcripts))

write_csv(substrate_summary, "substrate_summary.csv")

###Plot
# Optional: filter for bacterial transcripts
cazy_bacteria <- cazy_annotated_full_expr %>%
  filter(domain == "Bacteria")  # adjust column name if different

###OPTIONAL: NORMALIZE
# Compute total TPM per sample (for normalization)
total_tpm_per_sample <- cazy_bacteria %>%
  group_by(Sample) %>%
  summarise(total_TPM_sample = sum(TPM, na.rm = TRUE))

# Join and normalize
cazy_summary_by_sample_norm <- cazy_summary_by_sample %>%
  left_join(total_tpm_per_sample, by = "Sample") %>%
  mutate(relative_TPM = total_TPM / total_TPM_sample)
  
# Ensure NA genera are consistently labeled
cazy_summary_by_sample_norm <- cazy_summary_by_sample_norm %>%
  mutate(genus = if_else(is.na(genus), "Unclassified", as.character(genus)))
###----------------------------------------------------------
cazy_summary_by_sample <- cazy_bacteria %>%
  filter(!is.na(Substrate_high_level)) %>%
  group_by(Sample, Substrate_high_level, genus) %>%
  summarise(
    total_TPM = sum(TPM, na.rm = TRUE),
    n_transcripts = n(),
    .groups = "drop"
  )
  
# Define algal polysaccharide substrates of interest
target_substrates <- c("beta-glucan", "alginate", "chitin", "xylan", "beta-mannan", "exo-polysaccharide", "host glycan")
desired_sample_order <- c("DIAVIR9", "DIAVIR11", "DIAVIR4", "DIAVIR6", "DIAVIR8",
                          "DIAVIR5", "DIAVIR3", "DIAVIR12", "DIAVIR7", "DIAVIR10",
                          "DIAVIR13", "DIAVIR2", "DIAVIR1")

cazy_summary_by_sample_norm <- cazy_summary_by_sample_norm %>%
  mutate(Sample = factor(Sample, levels = desired_sample_order))
top15_genera <- cazy_summary_by_sample %>%
  group_by(genus) %>%
  summarise(total_TPM = sum(total_TPM, na.rm = TRUE)) %>%
  arrange(desc(total_TPM)) %>%
  slice_head(n = 15) %>%
  pull(genus)
  
cazy_summary_by_sample_norm <- cazy_summary_by_sample_norm %>%
  mutate(genus_grouped = if_else(genus %in% top15_genera, genus, "Other"))
 
library(ggsci)
pdf("cazy_summary_by_sample.pdf", width=9, height=16)
ggplot(cazy_summary_by_sample_norm %>% filter(Substrate_high_level %in% target_substrates),
       aes(x = Sample, y = relative_TPM, fill = genus_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Substrate_high_level, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_d3("category20")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
         strip.text.y = element_text(angle = 0, hjust = 0),
    strip.placement = "outside",
    panel.spacing.y = unit(1, "lines"),
	strip.text = element_text(size = 24),  # doubled from default 12 to 24
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  labs(title = "CAZyme TPMs per substrate and bacterial genus",
       y = "Total TPM", x = "Sample", fill = "Genus")
dev.off()


# ---------------------------
# 7. Analysis and exploration
# ---------------------------
library(tidyverse)
library(pheatmap)
 
  
# ---- Define CAZy marker families per pathway ----
pul_markers <- list(
  Laminarin = c("GH16", "GH30", "PL7", "SusD", "TonB"),
  Alginate = c("PL6", "PL7", "PL14", "SusD", "TonB"),
  Chitin = c("GH18", "GH19", "GH20", "CBM5", "CBM12"),
  Xylan = c("GH10", "GH11", "GH30"),
  Starch = c("GH13", "GH77", "SusD", "CBM20"),
  `Beta-glucan` = c("GH5", "GH16", "GH17")
)


# ---- Use your annotated_expr data frame ----
# Ensure Gene_ID, CAZy_Family, Sample, TPM are available

# ---- Filter out very low-expression (optional) ----
filtered_expr <- cazy_annotated_full_expr %>%
  filter(TPM > 0.1, !is.na(CAZy_Family))
  
# ---- Create presence table ----
presence <- filtered_expr %>%
  mutate(Present = 1) %>%
  distinct(Sample, CAZy_Family, Present) %>%
  pivot_wider(names_from = CAZy_Family, values_from = Present, values_fill = 0)

# ---- Compute completeness per pathway ----
completeness_df <- presence %>%
  # Exclude Sample column from conversion
  mutate(across(-Sample, ~ as.numeric(.x))) %>%
  # Replace any NAs with 0
  mutate(across(-Sample, ~ replace_na(.x, 0))) %>%
  # Now summarise to count presence across samples
  summarise(across(-Sample, sum)) %>%
  pivot_longer(everything(), names_to = "CAZy_Family", values_to = "Count")
  
# Function to compute completeness
calc_pathway_completeness <- function(df, markers) {
  sapply(markers, function(marker_vec) {
    rowMeans(df[, intersect(colnames(df), marker_vec)], na.rm = TRUE)
  }) %>% as.data.frame()
}

completeness <- calc_pathway_completeness(presence, pul_markers)
completeness$Sample <- presence$Sample

# ---- Plot pathway completeness ----
completeness_long <- completeness %>%
  pivot_longer(-Sample, names_to = "Pathway", values_to = "Completeness")
desired_sample_order <- c("DIAVIR9", "DIAVIR11", "DIAVIR4", "DIAVIR6", "DIAVIR8",
                          "DIAVIR5", "DIAVIR3", "DIAVIR12", "DIAVIR7", "DIAVIR10",
                          "DIAVIR13", "DIAVIR2", "DIAVIR1")

# Factor the Sample column with specified order
completeness_long$Sample <- factor(completeness_long$Sample, levels = desired_sample_order)
# Add Condition variable
completeness_long <- completeness_long %>%
  mutate(Condition = case_when(
    Sample == "DIAVIR9" ~ "Initial",
    Sample %in% c("DIAVIR11", "DIAVIR4", "DIAVIR6") ~ "control5",
    Sample %in% c("DIAVIR8", "DIAVIR5", "DIAVIR3") ~ "virus5",
    Sample %in% c("DIAVIR12", "DIAVIR7", "DIAVIR10") ~ "control9",
    Sample %in% c("DIAVIR13", "DIAVIR2", "DIAVIR1") ~ "virus9",
    TRUE ~ "Other"
  ))
# Keep Sample ordering fixed
completeness_long$Sample <- factor(completeness_long$Sample, levels = desired_sample_order)  
library(ggnewscale)
library(ggsci)

pdf("CAZYnalysis/cazy_completness.pdf", width=9, height= 12)
ggplot(completeness_long, aes(x = Sample, y = Completeness, fill = Pathway)) +
  geom_col(position = "dodge") +
  #coord_flip() +
  facet_grid(rows = vars(Condition), scales = "free_y", space = "free_y") +
  labs(
    title = "Pathway completeness per sample",
    y = "Completeness (0–1)",
    x = "Sample"
  ) +
  theme_minimal(base_size = 13) +
  scale_fill_npg() +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.placement = "outside",
    panel.spacing.y = unit(1, "lines"),
	strip.text = element_text(size = 24),  # doubled from default 12 to 24
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
dev.off() 
  
completeness_df <- presence %>%
  # Exclude Sample column from conversion
  mutate(across(-Sample, ~ as.numeric(.x))) %>%
  # Replace any NAs with 0
  mutate(across(-Sample, ~ replace_na(.x, 0))) %>%
  # Now summarise to count presence across samples
  summarise(across(-Sample, sum)) %>%
  pivot_longer(everything(), names_to = "CAZy_Family", values_to = "Count")


completeness_df %>%
  ggplot(aes(x = reorder(CAZy_Family, Count), y = Count)) +
  geom_col() +
  coord_flip() +
  labs(x = "CAZy Family", y = "Number of Samples with Expression",
       title = "CAZy Family Presence Across Samples") +
  theme_minimal()
  
marker_genes <- c("susd", "tonb", "gh16", "gh30", "pl7", "pl6", "pl14",
                  "gh18", "gh19", "gh20", "cbm5", "cbm12", "gh10", "gh11",
                  "gh13", "gh77", "cbm20", "gh5", "gh17", "gh1", "gh3")

heatmap_markers <- cazy_annotated_full_expr %>%
  filter(
    str_detect(tolower(Description), paste(marker_genes, collapse = "|")) |
    str_detect(tolower(PFAMs), paste(marker_genes, collapse = "|")) |
    str_detect(tolower(Preferred_name), paste(marker_genes, collapse = "|")) |
    tolower(CAZy_Family) %in% tolower(marker_genes)
  ) %>%
  mutate(Marker = case_when(
    str_detect(tolower(Description), "susd") ~ "SusD",
    str_detect(tolower(Description), "tonb") ~ "TonB",
    !is.na(CAZy_Family) ~ CAZy_Family,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Marker))
  
heatmap_data <- heatmap_markers %>%
  group_by(Sample, Marker) %>%
  summarise(Mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Mean_TPM, values_fill = 0) %>%
  column_to_rownames("Marker")  
				  
# Heatmap
pheatmap(as.matrix(heatmap_data),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Mean TPM of marker CAZymes across samples")


		 
#### -------------- Compute completness of pathways per genera using the dbCAN substrate mapping ----------------------


library(tidyverse)
library(patchwork)


# ---------------------------------------
# 2. Prepare completeness per sample
# ---------------------------------------
# Ensure NA genera are consistently labeled
cazy_summary_by_sample_norm <- cazy_summary_by_sample_norm %>%
  mutate(genus = if_else(is.na(genus), "Unclassified", as.character(genus)))

# Calculate top 15 genera INCLUDING Unclassified
top15_genera <- cazy_summary_by_sample_norm %>%
  filter(Substrate_high_level %in% target_substrates) %>%
  group_by(genus) %>%
  summarise(total = sum(total_TPM, na.rm = TRUE)) %>%
  arrange(desc(total)) %>%
  slice_max(order_by = total, n = 15) %>%
  pull(genus)

# Build rel_expr_matrix
rel_expr_matrix <- cazy_summary_by_sample_norm %>%
  filter(genus %in% top15_genera | genus == "Unclassified") %>%   # ✅ include "Unclassified"
  select(Sample, genus, Substrate_high_level, relative_TPM)

# Ensure Substrate-level completeness reference
ref_families_per_substrate <- cazy_annotated_full_expr %>%
  filter(domain == "Bacteria",
         TPM > 0,
         !is.na(Substrate_high_level),
         !is.na(CAZy_Family),
         Substrate_high_level %in% target_substrates) %>%
  distinct(Substrate_high_level, CAZy_Family)

# Recompute completeness per sample
completeness_per_sample <- cazy_annotated_full_expr %>%
  filter(domain == "Bacteria",
         TPM > 0,
         !is.na(Substrate_high_level),
         !is.na(CAZy_Family),
         Substrate_high_level %in% target_substrates) %>%
  distinct(Sample, Substrate_high_level, CAZy_Family) %>%
  group_by(Sample, Substrate_high_level) %>%
  summarise(detected = n_distinct(CAZy_Family), .groups = "drop") %>%
  left_join(
    ref_families_per_substrate %>%
      dplyr::count(Substrate_high_level, name = "total_ref_fams"),
    by = "Substrate_high_level"
  ) %>%
  mutate(completeness = detected / total_ref_fams)

# Substrates that have both expression and completeness data
valid_subs <- intersect(
  rel_expr_matrix %>% filter(relative_TPM > 0) %>% pull(Substrate_high_level) %>% unique(),
  completeness_per_sample %>% filter(completeness > 0) %>% pull(Substrate_high_level) %>% unique()
)

condition_map <- function(df) {
  df %>% mutate(Condition = case_when(
    Sample == "DIAVIR9" ~ "Initial-Day1", 
    Sample %in% c("DIAVIR11", "DIAVIR4", "DIAVIR6") ~ "control5",
    Sample %in% c("DIAVIR8", "DIAVIR5", "DIAVIR3") ~ "virus5",
    Sample %in% c("DIAVIR12", "DIAVIR7", "DIAVIR10") ~ "control9",
    Sample %in% c("DIAVIR13", "DIAVIR2", "DIAVIR1") ~ "virus9",
    TRUE ~ "Other"
  ))
}

# Extract red and green from the NPG palette
npg_colors <- pal_npg("nrc")(10)
npg_red   <- npg_colors[3]  # usually "#E64B35FF"
npg_green <- npg_colors[2]  # usually "#00A087FF"

plots <- list()


for (i in seq_along(valid_subs)) {
  sub <- valid_subs[i]

  sub_expr <- rel_expr_matrix %>% filter(Substrate_high_level == sub)
  # Label NA genus as 'Unclassified'
sub_expr <- sub_expr %>%
  mutate(genus = if_else(is.na(genus), "Unclassified", as.character(genus)))


total_tpm_by_genus_sample <- cazy_summary_by_sample_norm %>%
  filter(Substrate_high_level == sub) %>%
  mutate(genus = if_else(is.na(genus), "Unclassified", as.character(genus))) %>%
  filter(genus %in% c(top15_genera, "Unclassified")) %>%
  condition_map()

# Rebuild genus_levels with Unclassified included
genus_levels <- total_tpm_by_genus_sample %>%
  group_by(genus) %>%
  summarise(total = sum(total_TPM), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(genus)


  # Ensure both plots use same genus order
  sub_expr <- sub_expr %>%
    complete(Sample, genus = genus_levels, fill = list(relative_TPM = 0)) %>%
    mutate(genus = factor(genus, levels = genus_levels))

  total_tpm_by_genus_sample <- total_tpm_by_genus_sample %>%
    complete(genus = genus_levels, Sample, fill = list(total_TPM = 0, Condition = "Other")) %>%
    mutate(genus = factor(genus, levels = genus_levels))
	
tpm_by_genus_condition <- total_tpm_by_genus_sample %>%
  group_by(genus, Condition) %>%
  summarise(mean_TPM = mean(total_TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(genus = factor(genus, levels = genus_levels))
  
  bar_plot_data <- completeness_per_sample %>%
  filter(Substrate_high_level == sub) %>%
  mutate(Sample = factor(Sample, levels = desired_sample_order))

  # Expression heatmap
expr_plot <- ggplot(sub_expr,
                    aes(x = Sample,
                        y = genus,
                        fill = log10(relative_TPM + 1e-7))) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick")) +
  labs(title = sub, fill = "log10(RelTPM)", y = "Genus", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    legend.key.height = unit(0.3, "cm")  # ✅ smaller colorbar without shrinking text
  )

  # Add x-axis labels ONLY to the bottom-most heatmap
  if (i == length(valid_subs)) {
    expr_plot <- expr_plot +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_line()
      )
  }

  # Completeness barplot
bar_plot <-  ggplot(bar_plot_data,aes(x = Sample, y = completeness, fill = completeness)) +
  geom_col() +
  scale_fill_gradient2(
    low = "#D73027",     # red
    mid = "#FEE08B",     # yellow
    high = "#1A9850",    # green
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  ylim(0, 1) +
  theme_minimal() +
  labs(y = NULL, x = NULL, fill = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 3.5),
    legend.position = "none"
  )

  # Total TPM side barplot
side_bar <- ggplot(tpm_by_genus_condition,
                   aes(y = genus, x = mean_TPM, fill = Condition)) +
  geom_col(position = "stack", width = 0.6) +
  theme_minimal() +
  scale_fill_npg() +
  labs(x = "Mean total TPM", y = NULL, fill = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "lines")
  )
  # Define layout using patchwork::area
  layout <- c(
    area(t = 1, l = 1, b = 1, r = 1),  # bar_plot
    area(t = 2, l = 1, b = 2, r = 1),  # expr_plot
    area(t = 2, l = 2, b = 2, r = 2)   # side_bar
  )

  # Combine into one plot
  plots[[sub]] <- bar_plot + expr_plot + side_bar +
    plot_layout(design = layout, widths = c(4, 2), heights = c(1, 5))
}

pdf("Completness_relativeTPM_cazy.pdf", width=10, height=15)
wrap_plots(plots, ncol = 1) +
  plot_annotation(title = "CAZy substrate expression, pathway completeness, and total TPM by genus (per substrate)")
dev.off()



###Extract cazy completntess
# Step 1: Filter only CAZy families observed in your data
observed_families <- cazy_annotated_full_expr %>%
  filter(domain == "Bacteria",
         TPM > 0,
         !is.na(Substrate_high_level),
         !is.na(CAZy_Family),
         Substrate_high_level %in% target_substrates) %>%
  distinct(Sample, Substrate_high_level, CAZy_Family) %>%
  mutate(present = 1)

# Step 2: Get all combinations of samples × observed CAZy families per substrate
ref_families <- observed_families %>%
  distinct(Substrate_high_level, CAZy_Family)

all_combinations <- ref_families %>%
  crossing(Sample = unique(cazy_annotated_full_expr$Sample)) %>%
  left_join(observed_families, by = c("Sample", "Substrate_high_level", "CAZy_Family")) %>%
  mutate(present = replace_na(present, 0))

# Step 3: Pivot into wide matrix
presence_matrix <- all_combinations %>%
  pivot_wider(names_from = Sample, values_from = present)

# Step 4: Save to file
write_csv(presence_matrix, "data_derived_cazy_family_presence.csv")


###          ANNOTATE ALDEX RESULTS WITH CAZY FAMILIES         ######
library(tidyverse)
library(ggplot2)
library(forcats)
library(cowplot)
library(viridis)	
library(stringr)	 
# Step 1: Merge ALDEx2 results with CAZy annotations
aldex_annotated <- aldex_df %>%
  left_join(cazy_annotated_full, by = "transcript_id") %>%
  filter(!is.na(CAZy))  # keep only CAZyme-annotated transcripts
  
#Reintegrate transcripts with eggnog descriptions but no dbCAn annotation
# STEP 1: Identify susD and tonB hits in sig_bacteria
susd_tonb_rescued <- new_anotation %>%
  filter(
    str_detect(tolower(Description), "susd|ragb") |
    str_detect(tolower(Description), "tonb")
  ) %>%
  mutate(
    CAZy_Family = case_when(
      str_detect(tolower(Description), "susd|ragb") ~ "SusD",
      str_detect(tolower(Description), "tonb") ~ "TonB",
      TRUE ~ NA_character_
    )
  )

# STEP 2: Remove these from aldex_annotated if they somehow leaked in
aldex_cleaned <- aldex_annotated %>%
  filter(!transcript_id %in% susd_tonb_rescued$transcript_id)

# STEP 3: Bind the rescued susD/tonB transcripts to aldex_annotated
aldex_fixed <- bind_rows(aldex_cleaned, susd_tonb_rescued)

# Check presence
table(aldex_fixed$CAZy_Family)

# Step 2: Optional: highlight some marker CAZy families
marker_genes <- c("GH16", "GH30", "PL7", "SusD", "TonB", "PL6", "PL14",
                  "GH18", "GH19", "GH20", "CBM5", "CBM12", "GH10", "GH11",
                  "GH13", "GH77", "CBM20", "GH5", "GH17", "GH149", "GH149*")

aldex_markers <- aldex_fixed %>%
  filter(str_detect(CAZy, str_c(marker_genes, collapse = "|")))

aldex_fixed <- aldex_fixed %>%
  mutate(
    CAZy_Label = case_when(
      str_detect(tolower(Description), "susd") ~ "SusD",
      str_detect(tolower(Description), "tonb") ~ "TonB",
      !is.na(CAZy) & CAZy != "" & CAZy != "-" ~ CAZy,
      !is.na(CAZy_Family) & CAZy_Family != "" & CAZy_Family != "-" ~ CAZy_Family,
      TRUE ~ NA_character_
    )
  )
  
aldex_fixed <- aldex_fixed %>%
  mutate(
    CAZy_Label_clean = str_replace(CAZy_Label, "_\\d+$", ""),
    CAZy_Label_full  = CAZy_Label  # optional, if you want to retain full names
  )  
# Define significance threshold
sig_threshold <- 0.1

# Prepare data for plotting with susd/tonB transcripts

library(ggplot2)
library(ggsci)

p <- ggplot(aldex_plot_df, aes(x = effect, y = transcript_id, fill = CAZy_Label_clean)) +
  geom_col() +

  # Asterisks for significant genes
  geom_text(
    aes(label = asterisk),
    hjust = ifelse(aldex_plot_df$effect > 0, -0.3, 1.3),
    size = 4   # ~ 10–11 pt equivalent in pdf
  ) +

  # Right-side label: truncated description + genus
  geom_text(
    aes(label = label_with_genus),
    x = max(abs(aldex_plot_df$effect), na.rm = TRUE) * 1.25,
    hjust = 0,
    size = 3.5,          # readable on A4
    color = "gray30"
  ) +

  coord_cartesian(clip = "off") +
  scale_fill_d3(palette = "category20") +
  theme_minimal(base_size = 14) +   # increase base font scaling
  theme(
    axis.text.y  = element_text(color = "black", size = 10),  # ≥10 pt
    axis.text.x  = element_text(color = "black", size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.title.x = element_text(size = 12),
    plot.margin  = margin(10, 380, 10, 10, "pt"),

    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text  = element_text(size = 12),   # ≥12 pt
    legend.title = element_text(size = 12),

    plot.title    = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    x = "Effect Size (positive = up in virus5)",
    y = "Transcript ID",
    fill = "CAZy Family",
    title = "Top CAZymes by ALDEx2 Effect Size",
    subtitle = "* indicates FDR < 0.05"
  )

# --- Landscape A4 ---
# A4 dimensions (inches): 11.69 × 8.27 → using slightly extended height (9 in allowed)
pdf("cazy_susC.pdf", width = 11.7, height = 9)
p
dev.off()


###if we want to show tonB/susD only fro the exact significantly expressed transcripts
aldex_plot_df <- aldex_fixed %>%
  filter(if_any(c(CAZy, CAZy_Family), ~ str_detect(.x, str_c(marker_genes, collapse = "|")))) %>%
  mutate(
    is_significant = we.eBH < 0.1,
    short_desc = if_else(!is.na(Description), str_trunc(Description, 30), "No description")
  ) %>%
  arrange(desc(abs(effect))) %>%
  slice_head(n = 50) %>%
  mutate(
    transcript_id = fct_reorder(transcript_id, effect),
    asterisk = if_else(is_significant, "*", "")
  )  

aldex_plot_df$label_with_genus <- paste0(str_trunc(aldex_plot_df$Description, 80), " (", aldex_plot_df$genus, ")")
# Plot with manual second axis
library(ggsci)
p1 <- ggplot(aldex_plot_df, aes(x = effect, y = transcript_id, fill = CAZy)) +
  geom_col() +

  # Asterisks for significant genes
  geom_text(
    aes(label = asterisk),
    hjust = ifelse(aldex_plot_df$effect > 0, -0.3, 1.3),
    size = 5
  ) +

  # Right-side label: truncated description + genus
  geom_text(
    aes(label = label_with_genus),
    x = max(abs(aldex_plot_df$effect), na.rm = TRUE) * 1.25,
    hjust = 0,
    size = 2.9,
    color = "gray30"
  ) +

  coord_cartesian(clip = "off") +
  scale_fill_d3(palette = "category20") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(color = "black", size=8),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(5.5, 380, 5.5, 5.5, "pt"),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),       # smaller legend boxes
    legend.text = element_text(size = 8),    # smaller legend text
    legend.title = element_text(size = 9)    # optional: slightly smaller title
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  labs(
    x = "Effect Size (positive = up in virus5)",
    y = "Transcript ID",
    fill = "CAZy Family",
    title = "Top CAZymes by ALDEx2 Effect Size",
    subtitle = "* indicates FDR < 0.05"
  )


####use cazy substrate mappings
# Load the family-to-substrate map
sub_map <- read.delim("fam-substrate-mapping-08252022.tsv", sep = "\t", header = TRUE)

# Example structure: columns should include "Family" and "Substrate"
# head(sub_map)

# Load your dbCAN overview (must include a column like 'HMMER' or 'Family')
dbcan <- read.delim("overview.tsv", sep = "\t", header = TRUE)

# Extract family name (e.g., from 'GH16(24-327)')
dbcan$Family <- sub("\\(.*", "", dbcan$HMMER)

# Join to substrate map
dbcan_annotated <- merge(dbcan, sub_map, by.x = "Family", by.y = "Family")

# Count CAZymes by substrate category
library(dplyr)
summary <- dbcan_annotated %>%
  group_by(Substrate) %>%
  summarise(n_genes = n())

# Optional: merge with TPMs for expression analysis


#Aminopeptidases
aminopeptidase_hits <- aldex_fixed %>%
  filter(str_detect(tolower(Description), "aminopeptidase") |
         str_detect(tolower(Preferred_name), "aminopeptidase") |
         str_detect(tolower(PFAMs), "peptidase") |
         str_detect(tolower(Description), "exopeptidase") |
         str_detect(tolower(Preferred_name), "exopeptidase") |
         str_detect(tolower(Description), "exoprotease") |
         str_detect(tolower(Preferred_name), "exoprotease")) %>%
  filter(organism_type == "Bacteria")
		 
aminopeptidase_sig <- aminopeptidase_hits %>%
  filter(we.eBH < 0.1)  # or adjusted p-value / FDR threshold

# Replace with actual IDs
aminopeptidase_ids <- aminopeptidase_sig$transcript_id

# Filter TPMs for those genes
ap_tpm <- txi_tpm_long %>%
  filter(transcript_id %in% aminopeptidase_ids)  
  

# Wide format TPM matrix
amino_mat <- ap_tpm %>%
  pivot_wider(names_from = Sample, values_from = TPM)

# Add shortened description
amino_mat <- amino_mat %>%
  left_join(aldex_fixed %>% select(transcript_id, Description), by = "transcript_id") %>%
  mutate(short_desc = ifelse(is.na(Description), transcript_id, str_trunc(Description, 40)))

# Set description as rownames
amino_mat_final <- amino_mat %>%
  column_to_rownames("short_desc") %>%
  select(-transcript_id, -Description)
  
library(pheatmap)

# Optional log transformation
log_amino_mat <- log10(as.matrix(amino_mat_final) + 1)

pheatmap(log_amino_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize_row = 9,
         fontsize_col = 10,
         color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
         main = "Expression of Aminopeptidase/Exoprotease Genes (Bacterial)")




# Adjust the column name to match your data (e.g., EC or EC_Number)
peptidase_transcripts <- aldex_fixed %>%
  filter(str_detect(EC, "^3\\.4\\.11|^3\\.4\\.16")) %>%
  filter(organism_type=="Bacteria") %>%
  distinct(transcript_id, Description, EC) 
  
# === Step 2: Join with TPM matrix ===

# Ensure txi_tpm_long is in long format: transcript_id | Sample | TPM
peptidase_expr <- peptidase_transcripts %>%
  left_join(txi_tpm_long, by = "transcript_id") %>%
  mutate(short_desc = str_trunc(Description, 40))  # For heatmap labeling

# === Step 3: Pivot for heatmap ===

heatmap_data <- peptidase_expr %>%
  select(transcript_id, Sample, TPM, short_desc, EC) %>%
  pivot_wider(names_from = Sample, values_from = TPM, values_fill = 0)

# Matrix prep
rownames_mat <- paste(heatmap_data$short_desc, heatmap_data$EC, sep = " | ")
heatmap_matrix <- as.matrix(heatmap_data %>% select(-transcript_id, -short_desc, -EC))
rownames(heatmap_matrix) <- rownames_mat

# === Step 4: Plot ===
pdf("aminopep.pdf", width=9, height=12)
pheatmap(heatmap_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",  # Standardize expression across samples
         fontsize_row = 3,
         fontsize_col = 10,
         color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
         main = "Peptidase-related Gene Expression (TPM)")
dev.off()		 
#### ZOOM in on virus upregulated transcripts
# Assuming: `peptidase_expr` from the previous step
# And `sample_metadata`: data.frame(Sample, Condition)

# Step 1: Join condition info
peptidase_expr <- peptidase_expr %>%
  left_join(samples, by = c("Sample"= "sample"))
  peptidase_expr <- peptidase_expr %>%


# Step 2: Average expression per transcript per condition
mean_expr <- peptidase_expr %>%
  group_by(transcript_id, short_desc, EC, condition) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_TPM, values_fill = 0)

# Step 3: Define upregulated in virus
upregulated_ids <- mean_expr %>%
  filter(virus5 > control5 * 2 | virus9 > control9 *2 | virus5 > control9 * 2) %>%  # Can adjust threshold
  pull(transcript_id)


# Step 4: Subset original matrix
zoomed_expr <- peptidase_expr %>%
  filter(transcript_id %in% upregulated_ids) %>%
  select(transcript_id, Sample, TPM, short_desc, EC) %>%
  pivot_wider(names_from = Sample, values_from = TPM, values_fill = 0)

# Step 5: Prepare matrix for heatmap
heatmap_matrix_zoom <- as.matrix(zoomed_expr %>% select(-transcript_id, -short_desc, -EC))
rownames(heatmap_matrix_zoom) <- paste(zoomed_expr$short_desc, zoomed_expr$EC, sep = " | ")

# Step 6: Plot

pheatmap(heatmap_matrix_zoom,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         fontsize_row = 5,
         color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100),
         main = "Upregulated Peptidase Genes in Virus Samples")
		 
		 
#COG plots (Tinta 2023) 
# 1. Clean and filter
cog_data <- cazy_annotated_full_expr %>%
  filter(!is.na(COG_category), !is.na(genus)) %>%
  filter(domain=="Bacteria")%>%
  rename(COG_function = COG_category) %>%
  mutate(
    genus = str_trim(genus),
    Sample = as.character(Sample)  # Ensure Sample stays character
  ) %>%
  filter(str_detect(COG_function, "^[A-Z]$"))

# 2. Top 10 genera
top_genera <- cog_data %>%
  dplyr::count(genus, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  pull(genus)

# 3. Collapse low abundance genera
cog_data <- cog_data %>%
  mutate(genus_grouped = if_else(genus %in% top_genera, genus, "Other"))

# 4. Summarise TPMs per Sample × COG × Genus
summary_data <- cog_data %>%
  group_by(Sample, COG_function, genus_grouped) %>%
  summarise(total_tpm = sum(TPM), .groups = "drop")

# ✅ CHECK
stopifnot(!any(is.na(summary_data$Sample)))

# 5. Normalize within each Sample × COG
summary_data <- summary_data %>%
  group_by(Sample, COG_function) %>%
  mutate(pie_fraction = total_tpm / sum(total_tpm)) %>%
  ungroup()

# ✅ 3. Pivot only now: one row per Sample × COG with all genus columns
pie_data <- summary_data %>%
  select(Sample, COG_function, genus_grouped, pie_fraction) %>%
  pivot_wider(
    names_from = genus_grouped,
    values_from = pie_fraction,
    values_fill = 0
  )



cog_order <- sort(unique(pie_data$COG_function))

pie_data <- pie_data %>%
   mutate(Sample = factor(Sample, levels = desired_sample_order),
         x = as.numeric(Sample)) %>%
  mutate(        y = as.numeric(factor(COG_function, levels = cog_order))
  )
  
genus_cols <- setdiff(names(pie_data), c("Sample", "COG_function", "x", "y"))

#for scaling per total TPM
abund_data <- cog_data %>%
  group_by(Sample, COG_function, genus_grouped) %>%
  summarise(genus_tpm = sum(TPM), .groups = "drop")

# Compute total per Sample + COG
radius_data <- abund_data %>%
  group_by(Sample, COG_function) %>%
  summarise(total_tpm = sum(genus_tpm), .groups = "drop")
  
pie_data <- pie_data %>%
  left_join(radius_data, by = c("Sample", "COG_function")) %>%
  mutate(scaled_radius = sqrt(total_tpm / max(total_tpm, na.rm = TRUE)) * 0.9)
  

  
 #for scaling per rleativer  
pie_data <- pie_data %>%
  rowwise() %>%
  mutate(total_tpm = sum(c_across(all_of(genus_cols)), na.rm = TRUE)) %>%
  ungroup()

max_tpm <- max(pie_data$total_tpm, na.rm = TRUE)

pie_data <- pie_data %>%
  mutate(scaled_radius = sqrt(total_tpm / max_tpm) * 0.4)  # try 0.5 or 0.6
  
pdf("COG_pies_TPM_absolute.pdf", width = 16, height = 9)

ggplot() +
  geom_scatterpie(
    data = pie_data,
    aes(x = x, y = y, r = scaled_radius),
    cols = genus_cols,
    color = NA
  ) +
  scale_fill_d3(palette="category20") +
  scale_x_continuous(
    breaks = seq_along(desired_sample_order),
    labels = desired_sample_order,
    expand = c(0, 0)  # removes extra padding
  ) +
  scale_y_continuous(
    breaks = seq_along(cog_order),
    labels = cog_desc[cog_order],
    expand = c(0, 0)
  ) +
  coord_equal(
  xlim = c(0.5, length(sample_order) + 0.5),
  ylim = c(0.5, length(cog_order) + 0.5),
  expand = FALSE, clip="off"
) +
#for the absolute plot
geom_scatterpie_legend(
  x = max(pie_data$x) + 1.5,  # or adjust to move legend to right
  y = max(pie_data$y),        # near top of plot
  n = 3,                      # number of circles to show
  radius = c(0.3, 0.6, 0.9),  # legend radii to display
  labeller = function(r) round((r / 0.9)^2 * max(pie_data$total_tpm), 0)
)+
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  ) +
  labs(
    title = "Genus Composition of One-letter COG Functions per Sample",
    fill = "Genus"
  )

dev.off()


#DIATOM SCGs detected:
## --- SCG coverage for Bacillariophyta (2836|Bacillariophyta)
## --- using full-transcriptome quantifications (salmon_quant_clst) ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)

# --- Inputs from your ALDEX clustered quantification environment ----
stopifnot(exists("txi"), exists("eggnog_full_named"), exists("samples"))

# Assign TPM matrix from full quantification
tpm_mat <- txi$abundance
colnames(tpm_mat) <- samples$sample
tpm_annot <- eggnog_full_named

# --- Filter to diatom transcripts (Bacillariophyta) based on max_annot_lvl ----
if ("max_annot_lvl" %in% names(tpm_annot)) {
  diatom_ids <- tpm_annot$transcript_id[
    grepl("2836\\|Bacillariophyta", tpm_annot$max_annot_lvl, ignore.case = TRUE)
  ]
  message("Found ", length(diatom_ids), " transcripts annotated as Bacillariophyta.")
} else {
  stop("Column 'max_annot_lvl' not found in tpm_annot.")
}

# --- Define SCG COGs ----
scg_cogs <- c("COG0012","COG0016","COG0018","COG0172","COG0215",
              "COG0495","COG0525","COG0533","COG0541","COG0552")

# --- Identify SCG-annotated transcripts ----
ann_scg <- tpm_annot %>%
  filter(transcript_id %in% diatom_ids) %>%
  select(transcript_id, eggNOG_OGs) %>%
  mutate(eggNOG_OGs = as.character(eggNOG_OGs)) %>%
  separate_rows(eggNOG_OGs, sep = "[,;]") %>%
  mutate(COG = str_extract(eggNOG_OGs, "^[A-Z]+[0-9]+"),
         COG = str_replace(COG, "^KOG", "COG")) %>%   # <-- convert KOG→COG
  filter(COG %in% scg_cogs) %>%
  distinct(transcript_id, COG)

message("Identified ", nrow(ann_scg), " SCG transcripts (Bacillariophyta only).")

# --- Check overlap with TPM matrix ----
ann_scg <- ann_scg %>% filter(transcript_id %in% rownames(tpm_mat))
if (nrow(ann_scg) == 0) stop("No SCG transcripts found in TPM matrix after filtering.")

# --- TPM extraction and reshaping ----
scg_tpm_long <- tpm_mat[rownames(tpm_mat) %in% ann_scg$transcript_id, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("transcript_id") %>%
  pivot_longer(-transcript_id, names_to = "sample", values_to = "TPM") %>%
  left_join(ann_scg, by = "transcript_id") %>%
  left_join(samples %>% select(sample, condition), by = "sample")

# --- Collapse duplicates per COG per sample ----
scg_tpm_collapsed <- scg_tpm_long %>%
  group_by(sample, condition, COG) %>%
  summarise(
    TPM_max = max(TPM, na.rm = TRUE),
    TPM_sum = sum(TPM, na.rm = TRUE),
    n_tx    = n_distinct(transcript_id),
    .groups = "drop"
  )

# --- Compute coverage (presence) ----
tpm_threshold <- 1
scg_presence <- scg_tpm_collapsed %>%
  mutate(present = as.integer(TPM_max > tpm_threshold))

scg_coverage <- scg_presence %>%
  group_by(sample, condition) %>%
  summarise(
    detected_cogs = sum(present),
    total_cogs    = length(unique(COG)),
    coverage_percent = 100 * detected_cogs / total_cogs,
    .groups = "drop"
  ) %>%
  arrange(condition, sample)

print(scg_coverage)

# --- Plots and outputs ----
dir.create("diatoM_SCGs", showWarnings = FALSE)

pdf("diatoM_SCGs/scg_coverage_Bacillariophyta_fullTranscriptome.pdf", width = 5, height = 5)
ggplot(scg_coverage, aes(sample, coverage_percent, fill = condition)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(detected_cogs, "/", total_cogs)),
            vjust = -0.3, size = 3.5) +
  coord_cartesian(ylim = c(0, 110)) +
  theme_minimal(base_size = 13) +
  labs(title = "Bacillariophyta SCG coverage per sample (whole transcriptome)",
       x = "Sample", y = "SCG coverage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

write.csv(scg_coverage, "diatoM_SCGs/scg_coverage_per_sample_Bacillariophyta_fullTranscriptome.csv", row.names = FALSE)
write.csv(scg_tpm_collapsed, "diatoM_SCGs/scg_tpm_collapsed_Bacillariophyta_fullTranscriptome.csv", row.names = FALSE)

p_tile <- scg_presence %>%
  mutate(COG = factor(COG, levels = scg_cogs)) %>%
  ggplot(aes(x = sample, y = COG, fill = factor(present))) +
  geom_tile(color = "grey60") +
  scale_fill_manual(values = c("0"="white","1"="black"),
                    labels = c("Absent","Present"), name = "Detection") +
  theme_minimal(base_size = 12) +
  labs(title = paste0("Bacillariophyta SCG detection (TPM_max >", tpm_threshold, ")"),
       x = "Sample", y = "COG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf("diatoM_SCGs/scg_presence.pdf", width=5, height=5)
p_tile
dev.off()
p_cond <- scg_coverage %>%
  group_by(condition) %>%
  summarize(mean_cov = mean(coverage_percent),
            se_cov = sd(coverage_percent)/sqrt(dplyr::n()),
            .groups = "drop") %>%
  ggplot(aes(condition, mean_cov, fill = condition)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean_cov - se_cov, ymax = mean_cov + se_cov), width = 0.2) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal(base_size = 13) +
  labs(title = "Bacillariophyta SCG coverage by condition",
       x = "Condition", y = "SCG coverage (%)") +
  theme(legend.position = "none")
pdf("diatoM_SCGs/scg_condition.pdf", width=5, height=5)
p_cond
dev.off()


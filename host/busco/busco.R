## --- BUSCO-based transcriptome completeness per sample -----------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)

# --- 1. Load BUSCO full table (adjust path as needed) ----
busco_file <- "busco_stramenopiles/run_stramenopiles_odb10/full_table.tsv"

busco_tbl <- read.delim(busco_file, comment.char = "#", header = FALSE, sep = "\t",
                        col.names = c("BuscoID", "Status", "Sequence", "Score", "Length"))

# Keep only "Complete" BUSCOs
busco_complete <- busco_tbl %>%
  filter(Status %in% c("Complete", "Duplicated")) %>%
  distinct(BuscoID)

n_busco <- nrow(busco_complete)
message("Detected ", n_busco, " complete BUSCO orthologs.")

# --- 2. Ensure tpm_mat and samples are in your environment ----
stopifnot(exists("tpm_mat"), exists("samples"))
stopifnot(!is.null(colnames(tpm_mat)))

# --- 3. Extract transcript IDs that correspond to BUSCO hits ----
# BUSCO adds IDs from the transcriptome FASTA — match to your TPM IDs by overlap
busco_ids <- busco_complete$BuscoID
tpm_ids <- rownames(tpm_mat)

# Some BUSCO runs produce gene IDs like "g1234.t1" while your TPM rows may match
# only partially. If needed, simplify both sides:
clean_ids <- function(x) gsub("\\..*", "", x)
tpm_ids_clean <- clean_ids(tpm_ids)
busco_ids_clean <- clean_ids(busco_ids)

# --- 4. TPM detection per sample ----
# Subset TPM matrix to BUSCO transcripts
tpm_busco <- tpm_mat[tpm_ids[tpm_ids_clean %in% busco_ids_clean], , drop = FALSE]

# Define detection threshold
tpm_threshold <- 1

# Binary presence matrix (1 = expressed)
presence_mat <- tpm_busco > tpm_threshold
coverage <- colSums(presence_mat) / nrow(tpm_busco) * 100

# --- 5. Combine with sample metadata ----
busco_coverage <- data.frame(
  sample = colnames(tpm_busco),
  busco_detected = colSums(presence_mat),
  total_busco = nrow(tpm_busco),
  coverage_percent = coverage
) %>%
  left_join(samples %>% select(sample, condition), by = "sample") %>%
  arrange(condition, sample)

print(busco_coverage)

# --- 6. Plot BUSCO completeness per sample ----
p_busco <- ggplot(busco_coverage, aes(x = sample, y = coverage_percent, fill = condition)) +
  geom_col(color = "black") +
  geom_text(aes(label = paste0(busco_detected, "/", total_busco)),
            vjust = -0.3, size = 3.5) +
  coord_cartesian(ylim = c(0, 110)) +
  theme_minimal(base_size = 13) +
  labs(title = "BUSCO completeness per sample (TPM > 1)",
       x = "Sample", y = "Detected BUSCOs (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
print(p_busco)

# Optional export
write.csv(busco_coverage, "BUSCO_completeness_per_sample.csv", row.names = FALSE)
ggsave("BUSCO_completeness_per_sample.pdf", p_busco, width = 6, height = 4)

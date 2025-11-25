#KRAKEN processing
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(forcats)

# Set path to Kraken report files
report_files <- list.files("kraken_out", pattern = "_report.txt", full.names = TRUE)


parse_kraken_report <- function(file) {
  df <- read_tsv(file, col_names = FALSE, quote = "", show_col_types = FALSE)
  colnames(df) <- c("percent", "reads_clade", "reads_direct", "rank_code", "ncbi_taxid", "name")
  df <- df %>%
    mutate(name = str_trim(name),
           sample = str_extract(basename(file), "DIAVIR[0-9]+"))
  return(df)
}

kraken_all <- bind_rows(lapply(report_files, parse_kraken_report))

 kraken_genus <- kraken_all %>%
  filter(rank_code %in% c("F", "G")) %>%
  select(sample, rank_code, name, percent)

kraken_wide <- kraken_genus %>%
  pivot_wider(names_from = name, values_from = percent, values_fill = 0)

# Filter to only genus-level entries
genus_abundance <- kraken_genus %>%
  filter(rank_code == "G")


# Filter top 15 genera for clarity
top_genera <- kraken_genus %>%
  filter(rank_code == "G") %>%
  group_by(name) %>%
  summarise(mean_abundance = mean(percent, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:15) %>%
  pull(name)


# Collapse low-abundance genera into "Other"  
genus_abundance <- genus_abundance %>%
  mutate(name = ifelse(name %in% top_genera, name, "Other")) %>%
  group_by(sample, name) %>%
  summarise(percent = sum(percent, na.rm = TRUE), .groups = "drop")
  
 # Normalize so that each sample sums to 100%
genus_abundance <- genus_abundance %>%
  group_by(sample) %>%
  mutate(rel_abund = 100 * percent / sum(percent)) %>%
  ungroup()


# Create a named color palette
n_colors <- length(unique(genus_abundance$name))
colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_colors)  

desired_sample_order <- c("DIAVIR9", "DIAVIR11", "DIAVIR4", "DIAVIR6", "DIAVIR8",
                          "DIAVIR5", "DIAVIR3", "DIAVIR12", "DIAVIR7", "DIAVIR10",
                          "DIAVIR13", "DIAVIR2", "DIAVIR1")

genus_abundance$sample <- factor(genus_abundance$sample, levels = desired_sample_order)  
  
ggplot(genus_abundance, aes(x = sample, y = rel_abund, fill = fct_reorder(name, rel_abund, .fun = sum, .desc = TRUE))) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "Genus",
       title = "Genus-Level Relative Abundance (Top 20 Genera)") +
	   scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
# Filter top 15 families for clarity
# Filter to only genus-level entries
family_abundance <- kraken_genus %>%
  filter(rank_code == "F")


top_family <- kraken_genus %>%
  filter(rank_code == "F") %>%
  group_by(name) %>%
  summarise(mean_abundance = mean(percent, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:15) %>%
  pull(name)


# Collapse low-abundance genera into "Other"  
family_abundance <- family_abundance %>%
  mutate(name = ifelse(name %in% top_family, name, "Other")) %>%
  group_by(sample, name) %>%
  summarise(percent = sum(percent, na.rm = TRUE), .groups = "drop")
  
 # Normalize so that each sample sums to 100%
family_abundance <- family_abundance %>%
  group_by(sample) %>%
  mutate(rel_abund = 100 * percent / sum(percent)) %>%
  ungroup()


# Create a named color palette
n_colors <- length(unique(family_abundance$name))
colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_colors)  

desired_sample_order <- c("DIAVIR9", "DIAVIR11", "DIAVIR4", "DIAVIR6", "DIAVIR8",
                          "DIAVIR5", "DIAVIR3", "DIAVIR12", "DIAVIR7", "DIAVIR10",
                          "DIAVIR13", "DIAVIR2", "DIAVIR1")

family_abundance$sample <- factor(family_abundance$sample, levels = desired_sample_order)  
  
ggplot(family_abundance, aes(x = sample, y = rel_abund, fill = fct_reorder(name, rel_abund, .fun = sum, .desc = TRUE))) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Relative Abundance (%)", fill = "family",
       title = "Family-Level Relative Abundance (Top 15 Families)") +
	   scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

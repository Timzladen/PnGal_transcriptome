# ---- Setup ----
library(WGCNA)
library(tidyverse)

allowWGCNAThreads()        # use multithreading if available
options(stringsAsFactors = FALSE)
library(tximport)
library(ALDEx2)
library(dplyr)
library(ggplot2)
library(tidyr)
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
# Inputs you already have:
#   txi      : tximport list
#   samples  : data.frame with rownames = sample IDs, includes 'condition' and 'biomass'
# If you already computed vsd_mat, start from there. Otherwise:

dds_full <- DESeq2::DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)

# Filter: keep features expressed in >= 2 samples (adjust as needed)
keep <- rowSums(DESeq2::counts(dds_full) > 0) >= 2
dds_full <- dds_full[keep, ]

# VST (good for WGCNA)
vsd <- DESeq2::vst(dds_full, blind = TRUE)
vsd_mat <- SummarizedExperiment::assay(vsd)   # genes x samples

# Optional: variance filter (keeps top 75% most variable genes)
keep_var <- apply(vsd_mat, 1, var, na.rm = TRUE)
cut_var  <- quantile(keep_var, 0.25, na.rm = TRUE)
vsd_mat  <- vsd_mat[keep_var > cut_var, ]

# WGCNA expects samples in rows, genes in columns
datExpr <- t(vsd_mat)

# Basic sanity check
gsg <- goodSamplesGenes(datExpr, verbose = 3)
stopifnot(gsg$allOK)

# Cluster samples (optional outlier check)
sampleTree <- hclust(dist(datExpr), method = "average")
pdf("WGCNA/sample_clustering_new.pdf", 7, 5)
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex = 0.7)
dev.off()

# ---- Traits ----
# Binary one-hot for condition + continuous biomass
# Relevel as needed so "control5" is the reference in downstream modeling
samples$condition <- factor(samples$condition)
traits_onehot <- model.matrix(~ 0 + condition, data = samples)
colnames(traits_onehot) <- gsub("^condition", "", colnames(traits_onehot))
traits <- cbind(traits_onehot, biomass = as.numeric(samples$biomass))
stopifnot(identical(rownames(datExpr), rownames(traits)))

# ---- Soft-threshold (signed) ----
powers <- 1:20
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)

pdf("WGCNA/soft_threshold.pdf", 7, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "blue", lty = 2)

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")
dev.off()



# ---- blockwiseModules ----
# Key parameters to tune: minModuleSize, deepSplit, mergeCutHeight, maxBlockSize
net <- blockwiseModules(
  datExpr,
  power = 14,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  corType = "pearson",
  verbose = 3
)

moduleColors <- net$colors
MEs          <- net$MEs
geneTree     <- net$dendrograms[[1]]

# Plots: dendrogram + module colors
pdf("WGCNA/gene_dendro_module_colors.pdf", 9, 6)
plotDendroAndColors(geneTree, moduleColors, "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Order eigengenes and correlate with traits
MEs       <- orderMEs(MEs)
ME_names  <- colnames(MEs)

modTraitCor  <- cor(MEs, traits, use = "pairwise.complete.obs")
modTraitP    <- corPvalueStudent(modTraitCor, nSamples = nrow(traits))

# BH-adjust p-values per trait
modTraitQ <- apply(modTraitP, 2, p.adjust, method = "BH")
colnames(modTraitQ) <- colnames(traits)
rownames(modTraitQ) <- rownames(modTraitP)

# Heatmap with p (or q) overlay
textMat <- paste0(signif(modTraitCor, 2), "\n(",
                  ifelse(modTraitQ < 0.1, paste0("q=", signif(modTraitQ, 2)), paste0("p=", signif(modTraitP, 2))), ")")
dim(textMat) <- dim(modTraitCor)

pdf("WGCNA/module_trait_heatmap.pdf", 8, 8)
labeledHeatmap(
  Matrix = modTraitCor,
  xLabels = colnames(traits),
  yLabels = rownames(modTraitCor),
  ySymbols = rownames(modTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMat,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1,1),
  main = "Module–trait correlations (signed network)"
)
dev.off()


# Module assignments
module_df <- tibble(
  transcript_id = colnames(datExpr),
  module        = moduleColors
)
readr::write_tsv(module_df, "WGCNA/wgcna_module_assignments_blockwise.tsv")

# Eigengenes
MEs_out <- MEs %>%
  as.data.frame() %>%
  rownames_to_column("sample")
readr::write_tsv(MEs_out, "WGCNA/wgcna_module_eigengenes.tsv")

# Correlations / p / q tables
readr::write_tsv(
  as_tibble(modTraitCor, rownames = "module"),
  "WGCNA/module_trait_correlations.tsv"
)
readr::write_tsv(
  as_tibble(modTraitP, rownames = "module"),
  "WGCNA/module_trait_pvalues.tsv"
)
readr::write_tsv(
  as_tibble(modTraitQ, rownames = "module"),
  "WGCNA/module_trait_qvalues.tsv"
)

# ---- (Optional) merge very similar modules post hoc ----
# If you want to be more aggressive merging:
# ME_diss   <- 1 - cor(MEs)
# ME_tree   <- hclust(as.dist(ME_diss), method = "average")
# pdf("WGCNA/ME_clustering.pdf", 7, 5); plot(ME_tree); dev.off()
# mergeCutHeight2 <- 0.20  # stricter merge (cor > 0.80)
# merge <- mergeCloseModules(datExpr, moduleColors, cutHeight = mergeCutHeight2, verbose = 3)
# moduleColors <- merge$colors
# MEs <- orderMEs(merge$newMEs)
# (then re-run the module–trait correlation block & re-save)

message("Done: blockwiseModules pipeline complete.")


##ASSIGN AND INSPECT
black_genes     <- module_df %>% filter(module == "black")     %>% pull(transcript_id)
turquoise_genes <- module_df %>% filter(module == "turquoise") %>% pull(transcript_id)
brown_genes     <- module_df %>% filter(module == "brown")     %>% pull(transcript_id)

run_kegg_enrichment <- function(gene_set, res_annot) {
  term2gene <- res_annot %>%
    filter(!is.na(KEGG_ko)) %>%
    separate_rows(KEGG_ko, sep = ",") %>%
    select(term = KEGG_ko, gene = transcript_id) %>%
    distinct()

  enricher(
    gene          = intersect(gene_set, term2gene$gene),
    universe      = term2gene$gene,
    TERM2GENE     = term2gene,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.1
  )
}

enrich_black <- run_kegg_enrichment(black_genes, res_annot)
enrich_turqouise <- run_kegg_enrichment(turquoise_genes, res_annot)
enrich_brown <- run_kegg_enrichment(brown_genes, res_annot)

geneModuleMembership <- cor(datExpr, MEs, use = "p")
geneTraitSignificance <- cor(datExpr, traits[,"virus5", drop = FALSE], use = "p")

hub_black <- tibble(
  transcript_id = colnames(datExpr),
  module        = moduleColors,
  kME           = geneModuleMembership[, "MEblack"],
  GS_virus5     = geneTraitSignificance[, 1]
) %>%
  filter(module == "black") %>%
  arrange(desc(kME))

head(hub_black, 20)


#---- Function definition ----#
run_module_kegg_enrichment <- function(module_name, module_df, res_annot, out_dir = "WGCNA") {
  message("Running KEGG enrichment for module: ", module_name)

  # 1) Map KOs to transcripts (TERM2GENE)
  term2gene_full <- res_annot %>%
    dplyr::filter(!is.na(KEGG_ko)) %>%
    tidyr::separate_rows(KEGG_ko, sep = ",") %>%
    dplyr::select(term = KEGG_ko, gene = transcript_id, max_annot_lvl) %>%
    dplyr::distinct()

  term2gene <- term2gene_full %>%
    dplyr::select(term, gene) %>%
    dplyr::distinct()

  # 2) Module genes (transcript IDs)
  mod_genes <- module_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(transcript_id) %>%
    unique()

  # 3) Background = all genes with KOs
  background_genes <- unique(term2gene$gene)

  # 4) Enrichment
  enrich_results <- clusterProfiler::enricher(
    gene          = intersect(mod_genes, background_genes),
    universe      = background_genes,
    TERM2GENE     = term2gene,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.1
  )

  # If no terms enriched, still write an empty table with columns and return
  if (is.null(enrich_results) || nrow(enrich_results@result) == 0) {
    out_file_empty <- file.path(out_dir, paste0(module_name, "_module_enrichment.csv"))
    readr::write_csv(tibble::as_tibble(enrich_results@result), out_file_empty)
    message("→ Written (empty): ", out_file_empty)
    return(enrich_results)
  }

  # 5) TERM2NAME (KO → Description)
  term2name <- res_annot %>%
    dplyr::filter(!is.na(KEGG_ko)) %>%
    tidyr::separate_rows(KEGG_ko, sep = ",") %>%
    dplyr::select(KEGG_ko, Description) %>%
    dplyr::distinct() %>%
    dplyr::filter(stringr::str_detect(KEGG_ko, "^ko:")) %>%
    dplyr::rename(term = KEGG_ko, name = Description) %>%
    dplyr::group_by(term) %>% dplyr::slice(1) %>% dplyr::ungroup()

  # 6) Attach KO names
  res_tbl <- enrich_results@result %>%
    dplyr::left_join(term2name, by = c("ID" = "term")) %>%
    dplyr::mutate(Description = ifelse(!is.na(name), name, Description)) %>%
    dplyr::select(-name)

  # 7) Add per-KO summary of max_annot_lvl for genes in THIS module (from geneID column)
  #    geneID is a slash-separated list of transcript IDs contributing to that KO from 'gene'
  maxlvl_summary <- res_tbl %>%
    dplyr::select(ID, geneID) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(gene = geneID) %>%
    dplyr::left_join(term2gene_full %>% dplyr::select(term, gene, max_annot_lvl),
                     by = c("ID" = "term", "gene" = "gene")) %>%
    dplyr::filter(!is.na(max_annot_lvl)) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      module_max_annot_lvl_list = paste(unique(max_annot_lvl), collapse = ";"),
      module_max_annot_lvl_n    = dplyr::n(),
      .groups = "drop"
    )

  res_tbl <- res_tbl %>%
    dplyr::left_join(maxlvl_summary, by = "ID")

  # 8) Export final table
  out_file <- file.path(out_dir, paste0(module_name, "_module_enrichment.csv"))
  readr::write_csv(res_tbl, out_file)
  message("→ Written: ", out_file)

  # Put back into the enrichResult object for downstream plotting if desired
  enrich_results@result <- res_tbl
  return(enrich_results)
}

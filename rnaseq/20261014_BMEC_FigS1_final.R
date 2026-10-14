# Supplementary Figure S1 -- RNA-seq of hiPSC-derived EECM-BMECs
#
# Reproduces all RNA-seq panels of Fig. S1 from the gene count matrices
# Samples: EECM-BMECs from 3 healthy-control (HC), 2 MS donor lines, one primary human BMEC sample (GSE97575, GSM2572306) as reference.
# Cell-line names: HC1 = 1B, HC2 = 2B, HC3 = 3A, MS3 = 9B, MS4 = 10A.

# Paths (relative to the script folder)
data_dir <- "data"
prefix <- paste0(format(Sys.Date(), "%Y%m%d"), "_")
save_dir <- file.path("results", paste0(prefix, "BMEC"))
dir.create(save_dir, recursive = TRUE)
save_dir <- paste0(save_dir, "/")

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(dendextend)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(stringr)
library(openxlsx)
library(magick)

heatmap_colors <- colorRampPalette(c("#84a98c", "#FAFAF7", "#ff8a6a"))(100)

# donor -> cell line
line_names <- c(HC1 = "1B", HC2 = "2B", HC3 = "3A", MS3 = "9B", MS4 = "10A")

# save each panel as raw + LZW tiff
save_tiff_pair <- function(draw_fun, name, width, height) {
  f_raw <- paste0(save_dir, prefix, name, "_nocompression.tiff")
  tiff(f_raw, width = width, height = height, units = "in", res = 600)
  draw_fun()
  dev.off()
  image_write(image_read(f_raw), paste0(save_dir, prefix, name, ".tiff"), format = "tiff", compression = "LZW")
}
bold_heatmap_names <- function(hm, which = c("col_names", "row_names")) {
  for (w in which) {
    i <- which(hm$gtable$layout$name == w)
    if (length(i)) hm$gtable$grobs[[i]]$gp$fontface <- "bold"
  }
  hm
}

# INPUT
counts_ipsc <- read.delim(file.path(data_dir, "counts_EECM_BMECs.tsv"), row.names = 1, check.names = FALSE)
counts_primary <- read.delim(file.path(data_dir, "counts_primary_BMEC_line.tsv"), row.names = 1, check.names = FALSE)

# donor ID taken from the column name (e.g. "HC1_1_N_Merged" -> "HC1")
donor_ipsc <- sub("_.*$", "", colnames(counts_ipsc))
stopifnot(setequal(donor_ipsc, names(line_names)))

# PANELS A-B -- marker heatmaps (primary BMEC + 5 iPSC lines)
common_genes <- intersect(rownames(counts_ipsc), rownames(counts_primary))
counts_combined <- cbind(counts_ipsc[common_genes, ], counts_primary[common_genes, "SRR5441866", drop = FALSE])

info_bmec6 <- data.frame(
  row.names = colnames(counts_combined),
  donor = c(donor_ipsc, "Primary"),
  group = c(ifelse(grepl("^HC", donor_ipsc), "HC", "MS"), "Primary")
)
info_bmec6$label <- ifelse(info_bmec6$group == "Primary", "Primary\nBMEC",
                           paste0(line_names[info_bmec6$donor], "\n(", info_bmec6$group, ")"))

dds_bmec6 <- DESeqDataSetFromMatrix(countData = round(counts_combined),
                                    colData = info_bmec6, design = ~ group)
dds_bmec6 <- dds_bmec6[rowSums(counts(dds_bmec6)) >= 10, ]
vst_bmec6 <- vst(dds_bmec6, blind = TRUE)
vst_mat_bmec6 <- assay(vst_bmec6)
rownames(vst_mat_bmec6) <- sub("\\..*$", "", rownames(vst_mat_bmec6))
colnames(vst_mat_bmec6) <- info_bmec6$label

# column order: Primary, 1B, 2B, 3A, 9B, 10A
col_order6 <- c("Primary\nBMEC", paste0(line_names, "\n(", c("HC", "HC", "HC", "MS", "MS"), ")"))

group_factor_bmec6 <- factor(info_bmec6$group, levels = c("Primary", "HC", "MS"),
                             labels = c("Primary BMECs", "HC derived BMECs", "MS derived BMECs"))
annotation_bmec6 <- data.frame(group = group_factor_bmec6, row.names = info_bmec6$label)[col_order6, , drop = FALSE]
group_colors_3 <- list(group = c("Primary BMECs" = "#F7DAE5", "HC derived BMECs" = "#86A2E2", "MS derived BMECs" = "#CDCD8F"))

draw_bmec_marker_heatmap <- function(genes, title, name) {
  ids <- mapIds(org.Hs.eg.db, keys = genes, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  mat <- vst_mat_bmec6[ids, col_order6]
  rownames(mat) <- genes
  hm <- pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE, scale = "none",
                 color = heatmap_colors,
                 annotation_col = annotation_bmec6, annotation_colors = group_colors_3,
                 cellwidth = 35, cellheight = 35, fontsize_col = 9, fontsize_row = 13, angle_col = 0,
                 main = title, silent = TRUE)
  hm <- bold_heatmap_names(hm)
  save_tiff_pair(function() { grid::grid.newpage(); grid::grid.draw(hm$gtable) }, name, 10, 7)
  invisible(hm)
}

## Panel A
hm_A <- draw_bmec_marker_heatmap(c("PECAM1", "CDH5", "CLDN5", "OCLN", "TJP1", "SLC2A1", "ABCB1", "VWF"),
                                 "BMEC identity markers", "BMEC_A_identity_markers")
## Panel B
hm_B <- draw_bmec_marker_heatmap(c("ICAM1", "ICAM2", "SELE", "SELP", "VCAM1", "ALCAM", "MCAM", "CD99"),
                                 "Adhesion molecules", "BMEC_B_adhesion_molecules")

# DIFFERENTIAL EXPRESSION -- MS vs HC (1B, 2B, 3A vs 9B, 10A)
info_ipsc <- data.frame(
  row.names = colnames(counts_ipsc),
  donor = donor_ipsc,
  line = line_names[donor_ipsc],
  group = factor(ifelse(grepl("^HC", donor_ipsc), "HC", "MS"), levels = c("HC", "MS"))
)
info_ipsc$label <- paste0(info_ipsc$line, "\n(", info_ipsc$group, ")")
print(info_ipsc[, c("donor", "line", "group")])

dds_bmec <- DESeqDataSetFromMatrix(countData = round(counts_ipsc), colData = info_ipsc, design = ~ group)
dds_bmec <- dds_bmec[rowSums(counts(dds_bmec)) >= 10, ]
dds_bmec <- DESeq(dds_bmec)
vst_bmec <- vst(dds_bmec, blind = FALSE)

padj_cutoff <- 0.05
lfc_cutoff <- 1

## unshrunk results, used for the DEG list
res_unshrunk <- results(dds_bmec, contrast = c("group", "MS", "HC"))
res_unshrunk_df <- as.data.frame(res_unshrunk)
res_unshrunk_df$ensembl_id <- rownames(res_unshrunk_df)
res_unshrunk_df$ensembl_base <- sub("\\..*$", "", res_unshrunk_df$ensembl_id)
res_unshrunk_df$symbol <- unname(mapIds(org.Hs.eg.db, keys = res_unshrunk_df$ensembl_base,
                                        keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))

## shrunk LFC, only for the flag and the volcano
res_shrunk <- lfcShrink(dds_bmec, contrast = c("group", "MS", "HC"), res = res_unshrunk, type = "normal")
res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df$ensembl_base <- sub("\\..*$", "", rownames(res_shrunk_df))

## 2B is the only male donor, so drop the 4 Y-linked genes (KDM5D-DT, KDNBY, FAM224B, BCORP1). PAR1 genes like CRLF2 stay.
y_linked_ids <- c("ENSG00000288049", "ENSG00000260197", "ENSG00000230663", "ENSG00000215580")
res_unshrunk_df <- res_unshrunk_df[!res_unshrunk_df$ensembl_base %in% y_linked_ids, ]
res_shrunk_df <- res_shrunk_df[!res_shrunk_df$ensembl_base %in% y_linked_ids, ]

## DEGs: padj < 0.05 and |LFC| > 1 (unshrunk)
sig_genes <- res_unshrunk_df[!is.na(res_unshrunk_df$padj) & res_unshrunk_df$padj < padj_cutoff &
                               abs(res_unshrunk_df$log2FoldChange) > lfc_cutoff, ]

## flag DEGs where shrunk |LFC| is still > 1 (others grey in the table)
shrunk_lookup <- setNames(res_shrunk_df$log2FoldChange, res_shrunk_df$ensembl_base)
sig_genes$shrunk_log2FC <- shrunk_lookup[sig_genes$ensembl_base]
sig_genes$survives_shrinkage_filter <- abs(sig_genes$shrunk_log2FC) > lfc_cutoff

cat("BMEC DEGs (Y-linked removed):", nrow(sig_genes),
    "| Up:", sum(sig_genes$log2FoldChange > 0), "| Down:", sum(sig_genes$log2FoldChange < 0),
    "| Survive shrinkage:", sum(sig_genes$survives_shrinkage_filter), "\n")

## fix NA/placeholder symbols by hand (HGNC, GRCh38), else use the Ensembl ID
manual_symbol_overrides <- c(
  ENSG00000218537 = "MIF-AS1", 
  ENSG00000244280 = "ECEL1P2",
  ENSG00000220848 = "RPS18P9", 
  ENSG00000154608 = "CEP170P1",
  ENSG00000286619 = "PM20D1-AS1", 
  ENSG00000230630 = "DNM3OS",
  ENSG00000245025 = "RHOBTB2-AS1", 
  ENSG00000260910 = "SWINGN",
  ENSG00000283050 = "GTF2IP12", 
  ENSG00000231770 = "TMEM44-AS1",
  ENSG00000223715 = "LINC01208", 
  ENSG00000259436 = "LNROP"
)

placeholder <- grepl("^LOC[0-9]+$|^(RP[0-9]|AP0|AC[0-9]|CTD-|CTB-)", sig_genes$symbol)
sig_genes$symbol_final <- ifelse(
  !is.na(sig_genes$symbol) & sig_genes$symbol != "" & !placeholder,
  sig_genes$symbol,
  ifelse(sig_genes$ensembl_base %in% names(manual_symbol_overrides),
         manual_symbol_overrides[sig_genes$ensembl_base], sig_genes$ensembl_base))

# DEG TABLE (Excel, colour-coded)
# grey rows = shrunk |LFC| <= 1
deg_table <- sig_genes %>%
  arrange(padj) %>%
  select(ensembl_base, symbol = symbol_final, log2FoldChange,
         `log2FoldChange_moderated (Love et al., 2014)` = shrunk_log2FC,
         pvalue, padj, survives_shrinkage_filter)
survived <- deg_table$survives_shrinkage_filter
deg_table$survives_shrinkage_filter <- NULL

wb <- createWorkbook()
addWorksheet(wb, "BMEC_DEGs")
writeData(wb, "BMEC_DEGs", deg_table)
addStyle(wb, "BMEC_DEGs", createStyle(fontColour = "#000000"), rows = which(survived) + 1,
         cols = 1:ncol(deg_table), gridExpand = TRUE)
addStyle(wb, "BMEC_DEGs", createStyle(fontColour = "#A6A6A6"), rows = which(!survived) + 1,
         cols = 1:ncol(deg_table), gridExpand = TRUE)
saveWorkbook(wb, paste0(save_dir, prefix, "BMEC_DEG_table_unshrunk_colorcoded.xlsx"))

# PANEL D -- top-50 DEG heatmap (shrinkage survivors, z-scored)
sig_survivors <- sig_genes[sig_genes$survives_shrinkage_filter, ]
top50 <- sig_survivors %>% arrange(padj) %>% slice_head(n = 50)

mat <- assay(vst_bmec)[top50$ensembl_id, ]
rownames(mat) <- top50$symbol_final
colnames(mat) <- info_ipsc[colnames(mat), "label"]
mat_z <- t(scale(t(mat)))

ann_col <- data.frame(Group = info_ipsc$group, row.names = info_ipsc$label)
col_hc <- hclust(dist(t(mat_z)))
col_hc <- dendextend::rotate(col_hc, order = paste0(line_names, "\n(", c("HC", "HC", "HC", "MS", "MS"), ")"))
breaks_z <- seq(-max(abs(mat_z)), max(abs(mat_z)), length.out = 101)

hm_D <- pheatmap(mat_z, color = heatmap_colors, breaks = breaks_z,
                 cluster_rows = TRUE, cluster_cols = col_hc,
                 annotation_col = ann_col,
                 annotation_colors = list(Group = c(HC = "#86A2E2", MS = "#CDCD8F")),
                 angle_col = 0, fontsize_row = 9, fontsize_col = 11, cellwidth = 25, silent = TRUE)
hm_D <- bold_heatmap_names(hm_D, "col_names")
save_tiff_pair(function() { grid::grid.newpage(); grid::grid.draw(hm_D$gtable) }, "BMEC_D_top50_DEG_heatmap", 10, 7)

# PANEL C -- volcano (shrunk LFC, top 30 by padj labelled)
res_shrunk_df$sig <- "NS"
res_shrunk_df$sig[!is.na(res_shrunk_df$padj) & res_shrunk_df$padj < padj_cutoff &
                    res_shrunk_df$log2FoldChange > lfc_cutoff] <- "Up in MS"
res_shrunk_df$sig[!is.na(res_shrunk_df$padj) & res_shrunk_df$padj < padj_cutoff &
                    res_shrunk_df$log2FoldChange < -lfc_cutoff] <- "Down in MS"

top30 <- sig_genes %>% arrange(padj) %>% slice_head(n = 30)
top30$label <- top30$symbol_final
top30$is_ensembl <- grepl("^ENSG", top30$label)
top30_shrink <- merge(top30[, c("ensembl_base", "label", "is_ensembl")],
                      res_shrunk_df[, c("ensembl_base", "log2FoldChange", "padj")], by = "ensembl_base")

no_line_genes <- c("ERAP2", "SPAG4", "SIM1", "MYH14")
volcano_colors <- scale_color_manual(
  values = c("Up in MS" = "#ff8a6a", "Down in MS" = "#84a98c", "NS" = "grey80"),
  breaks = c("Up in MS", "Down in MS", "NS"),
  labels = c("Upregulated in MS", "Downregulated in MS", "Not significant"))
x_breaks <- c(-5, -2.5, 0, 2.5, 5)
xlim <- c(-6, 6)

erap2 <- top30_shrink[top30_shrink$label == "ERAP2", ]
rest <- top30_shrink[top30_shrink$label != "ERAP2", ]
rest_noline <- rest[rest$label %in% no_line_genes, ]
rest_line <- rest[!rest$label %in% no_line_genes, ]
combined_line <- rbind(
  transform(rest_line[!rest_line$is_ensembl, ], fontface_val = "bold", size_val = 3.6),
  transform(rest_line[rest_line$is_ensembl, ], fontface_val = "plain", size_val = 2.9))

top_panel <- ggplot(res_shrunk_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.6) + volcano_colors +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = erap2, aes(label = label), size = 3.6, color = "black", fontface = "bold",
                  segment.size = 0.3, min.segment.length = Inf, show.legend = FALSE) +
  coord_cartesian(xlim = xlim, ylim = c(90, 95)) +
  scale_x_continuous(breaks = x_breaks) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", plot.margin = margin(b = 0),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18)) +
  labs(y = NULL)

bottom_panel <- ggplot(res_shrunk_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.6) + volcano_colors +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = combined_line, aes(label = label, fontface = fontface_val, size = size_val),
                  color = "black", max.overlaps = 30, segment.size = 0.3, min.segment.length = 0, show.legend = FALSE) +
  geom_text_repel(data = rest_noline, aes(label = label), size = 3.6, color = "black", fontface = "bold",
                  max.overlaps = 30, segment.size = 0.3, min.segment.length = Inf, show.legend = FALSE) +
  scale_size_identity() +
  coord_cartesian(xlim = xlim, ylim = c(0, 15)) +
  scale_x_continuous(breaks = x_breaks) +
  theme_minimal(base_size = 9) +
  theme(plot.margin = margin(t = 2),
        axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.2, "cm"),
        # legend inside the plot, top-left, black box (as in Fig. S1C)
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.6)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold change (MS / HC)", y = "-log10 adjusted p-value", color = NULL)

volcano <- (top_panel / bottom_panel + plot_layout(heights = c(1, 4))) &
  theme(axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(color = "black"))
volcano[[2]] <- volcano[[2]] +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.text.x = element_text(color = "black"))

set.seed(42)
save_tiff_pair(function() print(volcano), "BMEC_C_volcano", width = 12, height = 8)

# PANEL E -- GO enrichment (BP + CC + MF, shrinkage survivors)
#    universe = all genes tested after the >=10-count filter (Y-free)
universe <- res_unshrunk_df$ensembl_base
go_input <- sig_survivors$ensembl_base

simplify_by_ont <- function(ont, cutoff_val) {
  ego <- enrichGO(gene = go_input, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = ont, pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (nrow(as.data.frame(ego)) == 0) return(NULL)
  simplify(ego, cutoff = cutoff_val, by = "p.adjust", select_fun = min)
}
# redundancy cutoffs: BP 0.7, CC 0.7, MF 0.4
ego_bp <- simplify_by_ont("BP", 0.7)
ego_cc <- simplify_by_ont("CC", 0.7)
ego_mf <- simplify_by_ont("MF", 0.4)

go_df <- do.call(rbind, list(
  if (!is.null(ego_bp)) cbind(as.data.frame(ego_bp), ONTOLOGY = "BP"),
  if (!is.null(ego_cc)) cbind(as.data.frame(ego_cc), ONTOLOGY = "CC"),
  if (!is.null(ego_mf)) cbind(as.data.frame(ego_mf), ONTOLOGY = "MF")))
go_df$GeneRatioNum <- sapply(strsplit(go_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
go_df <- go_df[order(go_df$ONTOLOGY, go_df$p.adjust), ]
go_df$Description <- factor(go_df$Description, levels = rev(unique(go_df$Description)))
print(table(go_df$ONTOLOGY))

go_dot <- ggplot(go_df, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_size(range = c(1, 5)) +
  scale_color_gradient(low = "#3B5B92", high = "#C3A8D1", trans = "log10", name = "p.adjust",
                       limits = c(min(go_df$p.adjust), 0.05), breaks = c(0.005, 0.01, 0.05)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "GeneRatio", y = NULL, title = "GO enrichment: BMEC DEGs (MS vs HC)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9, face = "bold"),
        strip.text.y = element_text(size = 10, face = "bold"),
        panel.spacing = unit(0.15, "lines"),
        plot.margin = margin(t = 2, b = 2, l = 2, r = 2),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.9, "cm")) +
  guides(color = guide_colorbar(reverse = TRUE, barheight = unit(3, "cm"), barwidth = unit(0.5, "cm")),
         size = guide_legend(override.aes = list(color = "#3B5B92")))

save_tiff_pair(function() print(go_dot), "BMEC_E_GO_dotplot", width = 7, height = max(4, nrow(go_df) * 0.13))


# Record software versions
writeLines(capture.output(sessionInfo()), paste0(save_dir, "sessionInfo.txt"))

cat("done:", save_dir, "\n")
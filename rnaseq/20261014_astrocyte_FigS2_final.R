# Supplementary Figure S2 -- RNA-seq of hiPSC-derived astrocytes
#
# Reproduces all RNA-seq panels of Fig. S2 from the gene count matrix
# Samples: hiPSC-derived astrocytes from 3 healthy-control (HC), 2 MS donor lines, 3 primary fetal astrocyte samples as reference.
# Cell-line names: HC1 = 1B, HC2 = 2B, HC3 = 3A, MS3 = 9B, MS4 = 10A.

# Paths (relative to the script folder)
data_dir <- "data"
prefix <- paste0(format(Sys.Date(), "%Y%m%d"), "_")
save_dir <- file.path("results", paste0(prefix, "astrocyte"))
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

# INPUT -- sample metadata (matched by SRA run ID)
counts <- read.delim(file.path(data_dir, "counts_astrocytes.tsv"), row.names = 1, check.names = FALSE)

srr_to_donor <- c(SRR7901104 = "HC1", SRR7901110 = "HC2", SRR7901116 = "MS4",
                  SRR7901122 = "MS2", SRR7901128 = "HC3", SRR7901134 = "MS1",
                  SRR7901140 = "MS3", SRR7901146 = "Fetal1", SRR7901148 = "Fetal2",
                  SRR7901150 = "Fetal3")
srr <- sub("^.*(SRR[0-9]+).*$", "\\1", colnames(counts))
stopifnot(all(srr %in% names(srr_to_donor)))

sample_info <- data.frame(row.names = colnames(counts), srr = srr, donor = srr_to_donor[srr])
sample_info$group <- sub("[0-9]+$", "", sample_info$donor) # HC / MS / Fetal
sample_info$label <- ifelse(sample_info$group == "Fetal",
                            paste0("Primary\n", sample_info$donor),
                            paste0(line_names[sample_info$donor], "\n(", sample_info$group, ")"))
print(sample_info[, c("srr", "donor", "group")])

# PANELS A-C -- marker heatmaps (Fetal1-3 + 1B, 2B, 3A, 9B, 10A)
p1_order <- c("Fetal1", "Fetal2", "Fetal3", "HC1", "HC2", "HC3", "MS3", "MS4")
p1_samples <- rownames(sample_info)[match(p1_order, sample_info$donor)]
info_p1 <- sample_info[p1_samples, ]

dds_p1 <- DESeqDataSetFromMatrix(countData = round(counts[, p1_samples]), colData = info_p1, design = ~ group)
dds_p1 <- dds_p1[rowSums(counts(dds_p1)) >= 10, ]
vst_p1 <- vst(dds_p1, blind = TRUE)
vst_mat_p1 <- assay(vst_p1)

marker_ids <- c(
  ALDH1L1 = "ENSG00000144908",
  S100B = "ENSG00000160307",
  AQP4 = "ENSG00000171885",
  GJA1 = "ENSG00000152661",
  GFAP = "ENSG00000131095",
  SLC1A3 = "ENSG00000079215",
  SLC1A2 = "ENSG00000110436",
  GLUL = "ENSG00000135821",
  VIM = "ENSG00000026025",
  STAT3 = "ENSG00000168610",
  C3 = "ENSG00000125730",
  CP = "ENSG00000047457",
  SERPINA3 = "ENSG00000196136",
  S100A10 = "ENSG00000197747",
  NES = "ENSG00000132688",
  SOX2 = "ENSG00000181449",
  FABP7 = "ENSG00000164434",
  TOP2A = "ENSG00000131747"
)

marker_expr <- vst_mat_p1[match(marker_ids, sub("\\..*$", "", rownames(vst_mat_p1))), ]
rownames(marker_expr) <- names(marker_ids)
colnames(marker_expr) <- info_p1$label

identity_genes <- c("ALDH1L1", "S100B", "AQP4", "GJA1", "GFAP", "SLC1A3", "SLC1A2", "GLUL")
reactive_genes <- c("C3", "SERPINA3", "CP", "STAT3", "S100A10", "GFAP", "VIM")
precursor_genes <- c("NES", "SOX2", "FABP7", "TOP2A")

annotation_p1 <- data.frame(
  group = factor(info_p1$group, levels = c("Fetal", "HC", "MS"),
                 labels = c("Primary fetal astrocytes", "HC iPSC-derived astrocytes", "MS iPSC-derived astrocytes")),
  row.names = info_p1$label)
group_colors_p1 <- list(group = c("Primary fetal astrocytes" = "#F7DAE5",
                                  "HC iPSC-derived astrocytes" = "#86A2E2",
                                  "MS iPSC-derived astrocytes" = "#CDCD8F"))

draw_marker_heatmap <- function(genes, title, name) {
  mat <- marker_expr[genes, ]
  breaks <- seq(min(mat), max(mat), length.out = 101) # colour scale per panel
  hm <- pheatmap(mat, annotation_col = annotation_p1, annotation_colors = group_colors_p1,
                 color = colorRampPalette(c("#84a98c", "#FAFAF7", "#ff8a6a"))(100), breaks = breaks,
                 cluster_rows = TRUE, cluster_cols = FALSE,
                 angle_col = 0, fontsize_col = 9, fontsize_row = 13,
                 cellwidth = 35, cellheight = 35, main = title, silent = TRUE)
  hm <- bold_heatmap_names(hm)
  save_tiff_pair(function() { grid::grid.newpage(); grid::grid.draw(hm$gtable) }, name, 10, 7)
  invisible(hm)
}

hm_A <- draw_marker_heatmap(identity_genes, "Core astrocyte identity markers", "astro_A_identity_markers")
hm_B <- draw_marker_heatmap(reactive_genes, "Reactive astrocyte markers", "astro_B_reactive_markers")
hm_C <- draw_marker_heatmap(precursor_genes, "Progenitor markers", "astro_C_progenitor_markers")

# DIFFERENTIAL EXPRESSION -- MS vs HC (1B, 2B, 3A vs 9B, 10A)
p2_samples <- rownames(sample_info)[sample_info$donor %in% names(line_names)]
info_p2 <- sample_info[p2_samples, ]
info_p2$group <- factor(info_p2$group, levels = c("HC", "MS"))

dds_p2 <- DESeqDataSetFromMatrix(countData = round(counts[, p2_samples]), colData = info_p2, design = ~ group)
dds_p2 <- dds_p2[rowSums(counts(dds_p2)) >= 10, ]
dds_p2 <- DESeq(dds_p2)
vst_p2 <- vst(dds_p2, blind = FALSE)

padj_cutoff <- 0.05
lfc_cutoff <- 1

## unshrunk results, used for the DEG list
res_p2 <- results(dds_p2, contrast = c("group", "MS", "HC"))
res_df <- as.data.frame(res_p2)
res_df$ensembl_id <- rownames(res_df)

sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) > lfc_cutoff, ]
sig_genes$ensembl_base <- sub("\\..*$", "", sig_genes$ensembl_id)
sig_genes$symbol <- unname(mapIds(org.Hs.eg.db, keys = sig_genes$ensembl_base,
                                  keytype = "ENSEMBL", column = "SYMBOL"))
sig_genes <- sig_genes[order(sig_genes$padj), ]

## shrunk LFC, only for the flag and the volcano
resLFC <- lfcShrink(dds_p2, contrast = c("group", "MS", "HC"), res = res_p2, type = "normal")
res_df_shrink <- as.data.frame(resLFC)
res_df_shrink$ensembl_id <- rownames(res_df_shrink)

sig_genes$shrunk_log2FC <- setNames(res_df_shrink$log2FoldChange, res_df_shrink$ensembl_id)[sig_genes$ensembl_id]
sig_genes$survives_shrinkage_filter <- abs(sig_genes$shrunk_log2FC) > lfc_cutoff

cat("Astrocyte DEGs:", nrow(sig_genes),
    "| Up:", sum(sig_genes$log2FoldChange > 0), "| Down:", sum(sig_genes$log2FoldChange < 0),
    "| Survive shrinkage:", sum(sig_genes$survives_shrinkage_filter), "\n")

## Donor-sex check (2B is the only male donor): DEGs annotated to chrY.
## Only pseudoautosomal (PAR1) genes are found (CD99, ZBED1) -> none removed.
gene_chr <- AnnotationDbi::select(org.Hs.eg.db, keys = sig_genes$ensembl_base,
                                  keytype = "ENSEMBL", columns = c("SYMBOL", "CHR"))
cat("chrY-annotated DEGs:",
    paste(unique(gene_chr$SYMBOL[!is.na(gene_chr$CHR) & gene_chr$CHR == "Y"]), collapse = ", "), "\n")

## fix NA/placeholder symbols by hand (HGNC, GRCh38), else use the Ensembl ID
manual_symbol_overrides <- c(
  ENSG00000177822 = "TENM3-AS1", ENSG00000251273 = "LINC02228",
  ENSG00000231764 = "DLX6-AS1", ENSG00000286214 = "COPG2IT1",
  ENSG00000180229 = "HERC2P3", ENSG00000246100 = "LINC00900",
  ENSG00000237552 = "LINC02567", ENSG00000285744 = "MITA1",
  ENSG00000211829 = "TRDC", ENSG00000230910 = "ADGRB3-DT",
  ENSG00000234814 = "SVIL2P", ENSG00000242797 = "GLYCTK-AS1",
  ENSG00000230202 = "RPL29P4", ENSG00000286619 = "PM20D1-AS1",
  ENSG00000254416 = "LINC02732", ENSG00000179277 = "MEIS3P1",
  ENSG00000185607 = "ACTBP7", ENSG00000226702 = "MIR217HG",
  ENSG00000223749 = "MIR503HG", ENSG00000232044 = "SILC1",
  ENSG00000224078 = "SNHG14", ENSG00000235505 = "CASP4LP",
  ENSG00000175170 = "FAM182B", ENSG00000189212 = "DPY19L2P1",
  ENSG00000164393 = "ADGRF2P", ENSG00000270607 = "LEISA1",
  ENSG00000244227 = "LRRC77P", ENSG00000237949 = "LINC00844",
  ENSG00000259867 = "DYNLRB2-AS1")

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
addWorksheet(wb, "Astrocyte_DEGs")
writeData(wb, "Astrocyte_DEGs", deg_table)
addStyle(wb, "Astrocyte_DEGs", createStyle(fontColour = "#000000"), rows = which(survived) + 1,
         cols = 1:ncol(deg_table), gridExpand = TRUE)
addStyle(wb, "Astrocyte_DEGs", createStyle(fontColour = "#A6A6A6"), rows = which(!survived) + 1,
         cols = 1:ncol(deg_table), gridExpand = TRUE)
saveWorkbook(wb, paste0(save_dir, prefix, "astro_DEG_table_unshrunk_colorcoded.xlsx"))

# PANEL E -- top-50 DEG heatmap (shrinkage survivors, z-scored)
sig_survivors <- sig_genes[sig_genes$survives_shrinkage_filter, ]
top50 <- sig_survivors %>% arrange(padj) %>% slice_head(n = 50)

mat <- assay(vst_p2)[top50$ensembl_id, ]
rownames(mat) <- top50$symbol_final
colnames(mat) <- info_p2[colnames(mat), "label"]
mat_z <- t(scale(t(mat)))

ann_col <- data.frame(Group = info_p2$group, row.names = info_p2$label)
col_hc <- hclust(dist(t(mat_z)))
col_hc <- dendextend::rotate(col_hc, order = c("3A\n(HC)", "1B\n(HC)", "2B\n(HC)", "10A\n(MS)", "9B\n(MS)"))
breaks_z <- seq(-max(abs(mat_z)), max(abs(mat_z)), length.out = 101)

hm_E <- pheatmap(mat_z, color = heatmap_colors, breaks = breaks_z,
                 cluster_rows = TRUE, cluster_cols = col_hc,
                 annotation_col = ann_col,
                 annotation_colors = list(Group = c(HC = "#86A2E2", MS = "#CDCD8F")),
                 angle_col = 0, fontsize_row = 6, fontsize_col = 8, cellwidth = 25, silent = TRUE)
hm_E <- bold_heatmap_names(hm_E, "col_names")
save_tiff_pair(function() { grid::grid.newpage(); grid::grid.draw(hm_E$gtable) }, "astro_E_top50_DEG_heatmap", 10, 7)

# PANEL D -- volcano (shrunk LFC, top 30 by padj labelled)
res_df_shrink$sig <- "NS"
res_df_shrink$sig[!is.na(res_df_shrink$padj) & res_df_shrink$padj < padj_cutoff &
                    res_df_shrink$log2FoldChange > lfc_cutoff] <- "Up in MS"
res_df_shrink$sig[!is.na(res_df_shrink$padj) & res_df_shrink$padj < padj_cutoff &
                    res_df_shrink$log2FoldChange < -lfc_cutoff] <- "Down in MS"

top30 <- sig_genes %>% arrange(padj) %>% slice_head(n = 30)
top30$label <- top30$symbol
top30_shrink <- merge(top30[, c("ensembl_id", "label")],
                      res_df_shrink[, c("ensembl_id", "log2FoldChange", "padj")], by = "ensembl_id")

top_labels <- top30_shrink[top30_shrink$label %in% c("CD74", "GDF10"), ] # top strip
bottom_labels <- top30_shrink[!top30_shrink$label %in% c("CD74", "GDF10"), ]
no_line <- c("RBM47", "NR2E1", "HLA-DRA")
bottom_noline <- bottom_labels[bottom_labels$label %in% no_line, ]
bottom_line <- bottom_labels[!bottom_labels$label %in% no_line, ]

volcano_colors <- scale_color_manual(
  values = c("Up in MS" = "#ff8a6a", "Down in MS" = "#84a98c", "NS" = "grey80"),
  breaks = c("Up in MS", "Down in MS", "NS"),
  labels = c("Upregulated in MS", "Downregulated in MS", "Not significant"))
xlim <- c(-10, 10); top_ylim <- c(50, 82); bottom_ylim <- c(0, 40)

top_panel <- ggplot(res_df_shrink, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.6) + volcano_colors +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_labels, aes(label = label), size = 3.6, color = "black", fontface = "bold",
                  segment.size = 0.3, min.segment.length = Inf, show.legend = FALSE) +
  coord_cartesian(xlim = xlim, ylim = top_ylim) +
  scale_y_continuous(breaks = seq(top_ylim[1], top_ylim[2], 10)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", plot.margin = margin(b = 0),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18)) +
  labs(y = NULL)

bottom_panel <- ggplot(res_df_shrink, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.6) + volcano_colors +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = bottom_noline, aes(label = label), size = 3.6, color = "black", fontface = "bold",
                  max.overlaps = 30, segment.size = 0.3, min.segment.length = Inf, show.legend = FALSE) +
  geom_text_repel(data = bottom_line, aes(label = label), size = 3.6, color = "black", fontface = "bold",
                  max.overlaps = 30, segment.size = 0.3, min.segment.length = 0, show.legend = FALSE) +
  coord_cartesian(xlim = xlim, ylim = bottom_ylim) +
  scale_y_continuous(breaks = seq(bottom_ylim[1], bottom_ylim[2], 10)) +
  theme_minimal(base_size = 9) +
  theme(plot.margin = margin(t = 2),
        axis.text.x = element_text(size = 16), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 16), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.2, "cm"),
        # legend inside the plot, top-right, black box (as in Fig. S2D)
        legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.6)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(x = "log2 fold change (MS / HC)", y = "-log10 adjusted p-value", color = NULL)

volcano <- (top_panel / bottom_panel + plot_layout(heights = c(1, 3))) &
  theme(axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.text.y = element_text(color = "black"))
volcano[[2]] <- volcano[[2]] +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5),
        axis.text.x = element_text(color = "black"))

set.seed(42)
save_tiff_pair(function() print(volcano), "astro_D_volcano", width = 12, height = 8)

# PANEL F -- GO Biological Process (pooled up+down, shrinkage survivors)
#    universe = all genes tested after the >=10-count filter
universe <- sub("\\..*", "", rownames(res_p2))
go_input <- sig_survivors$ensembl_base

go_bp <- enrichGO(gene = go_input, universe = universe, OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2)
cat("Significant BP terms before simplify:", nrow(as.data.frame(go_bp)), "\n")
go_bp_simpl <- simplify(go_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

go_df <- as.data.frame(go_bp_simpl)[1:10, ]
go_df$GeneRatioNum <- sapply(strsplit(go_df$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
go_df$Description <- factor(go_df$Description, levels = rev(go_df$Description))

go_dot <- ggplot(go_df, aes(x = GeneRatioNum, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_size(range = c(1, 5)) +
  scale_color_gradient(low = "#3B5B92", high = "#C3A8D1", trans = "log10", name = "p.adjust") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "GeneRatio", y = NULL, title = "GO Biological Process: MS vs HC astrocytes") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 15, face = "bold"),
        plot.margin = margin(t = 2, b = 2, l = 2, r = 2),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.9, "cm")) +
  guides(color = guide_colorbar(reverse = TRUE, barheight = unit(3, "cm"), barwidth = unit(0.5, "cm")),
         size = guide_legend(override.aes = list(color = "#3B5B92")))

save_tiff_pair(function() print(go_dot), "astro_F_GO_BP_top10", width = 9.5, height = 7)


# Record software versions
writeLines(capture.output(sessionInfo()), paste0(save_dir, "sessionInfo.txt"))

cat("done:", save_dir, "\n")
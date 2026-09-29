# ============================================================================== 
# PAIP1 普通转录组分析
# ============================================================================== 

# 1. 环境准备
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(enrichplot)
library(ggsignif)
library(stringr)

setwd("D:/LK/普通转录组/")

color_ctrl <- "#f933ff"
color_cko  <- "#2591a8"
my_colors <- c("Control" = color_ctrl, "Paip1-cKO" = color_cko)
my_labels <- c(
  "Control" = "Control",
  "Paip1-cKO" = expression(italic("Paip1") * "-cKO")
)

# 2. 数据读取与样本信息
counts_file <- "quant/PAIP1_counts.txt"

raw_counts <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

count_matrix <- raw_counts[, 7:ncol(raw_counts)]
colnames(count_matrix) <- c("cKO_1", "cKO_2", "cKO_3", "Ctrl_1", "Ctrl_2", "Ctrl_3")

clean_ensembl <- gsub("\\..*", "", raw_counts$Geneid)
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = clean_ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
gene_symbols[is.na(gene_symbols)] <- clean_ensembl[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)
rownames(count_matrix) <- gene_symbols

colData <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c(rep("Paip1-cKO", 3), rep("Control", 3)),
    levels = c("Control", "Paip1-cKO")
  )
)

# 3. DESeq2 差异分析
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = colData,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Paip1-cKO", "Control"))
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>%
  drop_na(padj)

padj_threshold <- 0.05
log2FC_threshold <- 0.58

res_df <- res_df %>%
  mutate(
    Status = case_when(
      padj < padj_threshold & log2FoldChange > log2FC_threshold ~ "Up",
      padj < padj_threshold & log2FoldChange < -log2FC_threshold ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(res_df, "PAIP1_Optimized_DESeq2_Results.csv", row.names = FALSE)

# 4. 靶向火山图
num_up <- sum(res_df$Status == "Up")
num_down <- sum(res_df$Status == "Down")

key_target_genes <- c(
  "Paip1", "Ybx2", "Rad51", "Dmc1", "Rpa2",
  "Stra8", "Dazl", "Sycp1", "Sycp3", "Mlh1"
)

label_genes <- res_df %>%
  filter(gene_id %in% key_target_genes & Status != "NS")

pdf("Fig1_Targeted_Volcano_Plot.pdf", width = 8, height = 7)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Status)) +
  geom_point(alpha = 0.6, size = 3.5) +
  scale_color_manual(
    values = c(
      "Up" = color_cko,
      "Down" = color_ctrl,
      "NS" = "#DFDFDF"
    )
  ) +
  geom_vline(
    xintercept = c(-log2FC_threshold, log2FC_threshold),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  geom_hline(
    yintercept = -log10(padj_threshold),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_id),
    size = 5.5,
    fontface = "bold.italic",
    box.padding = 0.8,
    point.padding = 0.3,
    color = "black",
    segment.color = "black",
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = paste0("Down: ", num_down, " genes"),
    vjust = 2,
    hjust = -0.2,
    color = color_ctrl,
    size = 6,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Up: ", num_up, " genes"),
    vjust = 2,
    hjust = 1.2,
    color = color_cko,
    size = 6,
    fontface = "bold"
  ) +
  coord_cartesian(xlim = c(-8, 4.5)) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  ) +
  labs(
    title = expression("Volcano Plot (" * italic("Paip1") * "-cKO vs Control)"),
    x = expression(Log[2] ~ Fold ~ Change),
    y = expression(-Log[10] ~ P[adj])
  )

dev.off()

# 5. PCA
vsd <- vst(dds, blind = FALSE)

pdf("QC_1_PCA_Plot.pdf", width = 7, height = 5.5)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 6, alpha = 0.9) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = my_colors, labels = my_labels) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  ggtitle("PCA (Corrected Conditions)")

dev.off()

# 6. 靶向基因热图
key_target_genes <- c(
  "Paip1", "Sycp1", "Sycp3", "Hormad1", "Mlh1", "Spo11",
  "Atr", "H2ax", "Rpa2", "Rad51", "Dmc1"
)

vst_mat <- assay(vsd)
target_pathway_genes <- intersect(key_target_genes, rownames(vst_mat))
targeted_matrix <- vst_mat[target_pathway_genes, ]

ordered_cols <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "cKO_1", "cKO_2", "cKO_3")
targeted_matrix_ordered <- targeted_matrix[, ordered_cols]
ann_colors <- list(condition = my_colors)

pdf("Fig4_Targeted_Heatmap.pdf", width = 6.5, height = 7)

pheatmap(
  targeted_matrix_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = 3,
  scale = "row",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = colData,
  annotation_colors = ann_colors,
  fontsize = 14,
  main = "Key Genes Expression"
)

dev.off()

# 7. 靶向基因箱线图
norm_counts <- counts(dds, normalized = TRUE)

genes_for_boxplot <- c(
  "Paip1", "Ybx2", "Rad51", "Dmc1", "Rpa2",
  "Stra8", "Dazl", "Sycp1", "Sycp3", "Mlh1"
)

dir.create("Boxplots_Individual", showWarnings = FALSE)

for (gene in genes_for_boxplot) {
  if (gene %in% rownames(norm_counts)) {
    gene_data <- data.frame(
      Sample = colnames(norm_counts),
      Expression = norm_counts[gene, ],
      Group = colData$condition
    )

    gene_res <- res_df[res_df$gene_id == gene, ]
    sig_star <- "ns"

    if (nrow(gene_res) > 0 && !is.na(gene_res$pvalue)) {
      raw_p <- gene_res$pvalue
      if (raw_p < 0.001) {
        sig_star <- "***"
      } else if (raw_p < 0.01) {
        sig_star <- "**"
      } else if (raw_p < 0.05) {
        sig_star <- "*"
      }
    }

    y_max <- max(gene_data$Expression)
    y_pos <- y_max + (y_max * 0.1)

    p_box <- ggplot(gene_data, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5) +
      geom_jitter(width = 0.15, size = 4, color = "black", alpha = 0.8) +
      scale_fill_manual(values = my_colors, labels = my_labels) +
      scale_x_discrete(labels = my_labels) +
      geom_signif(
        comparisons = list(c("Control", "Paip1-cKO")),
        annotations = sig_star,
        y_position = y_pos,
        tip_length = 0.02,
        textsize = 8,
        vjust = 0.2,
        color = "black"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
      theme_classic(base_size = 20) +
      labs(
        title = bquote(italic(.(gene))),
        y = "Normalized Counts",
        x = ""
      ) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 26),
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 18),
        axis.title.y = element_text(face = "bold", size = 20)
      )

    file_name <- paste0("Boxplots_Individual/Fig5_Boxplot_", gene, ".pdf")
    ggsave(file_name, plot = p_box, width = 4.5, height = 5.5)
  }
}

# 8. GSEA
res_gsea <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

gene_list_gsea <- res_gsea$log2FoldChange
names(gene_list_gsea) <- res_gsea$gene_id

id_trans <- bitr(
  names(gene_list_gsea),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

res_gsea_entrez <- res_gsea %>%
  inner_join(id_trans, by = c("gene_id" = "SYMBOL"))

gene_list_entrez <- sort(
  setNames(res_gsea_entrez$log2FoldChange, res_gsea_entrez$ENTREZID),
  decreasing = TRUE
)

gse_bp <- gseGO(
  geneList = gene_list_entrez,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

core_pathway_ids <- c("GO:0140013", "GO:0061982", "GO:0006302")
valid_ids <- intersect(core_pathway_ids, as.data.frame(gse_bp)$ID)

if (length(valid_ids) > 0) {
  pdf("Fig6_Targeted_GSEA_Meiosis_DSB_Triple.pdf", width = 8.5, height = 7)

  p_gsea <- gseaplot2(
    gse_bp,
    geneSetID = valid_ids,
    color = c("#E64B35", "#3C5488", "#00A087"),
    title = "Meiosis Arrest & DSB Repair Dysfunction",
    pvalue_table = TRUE,
    base_size = 18
  )

  print(p_gsea)
  dev.off()
}

# 9. 下调基因 GO 富集
down_genes <- res_df %>%
  filter(Status == "Down") %>%
  pull(gene_id)

down_entrez <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

pdf("Fig7_GO_Enrichment_Dotplot_Down.pdf", width = 8, height = 7)

p_go <- dotplot(
  ego_down,
  split = "ONTOLOGY",
  showCategory = 8
) +
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  scale_color_viridis_c(option = "C", direction = -1) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(color = "black", size = 12)
  ) +
  labs(
    title = expression(
      "GO Enrichment of Down-regulated Genes in " * italic("Paip1") * "-cKO"
    )
  )

print(p_go)
dev.off()

# 10. 下调基因 KEGG 富集
ekegg_down <- enrichKEGG(
  gene = down_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

ekegg_down <- setReadable(
  ekegg_down,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID"
)

pdf("Fig8_KEGG_Enrichment_Barplot_Down.pdf", width = 8, height = 6)

p_kegg <- barplot(ekegg_down, showCategory = 15) +
  scale_fill_viridis_c(option = "C", direction = -1) +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_text(color = "black", size = 12)) +
  labs(
    title = expression(
      "KEGG Pathways of Down-regulated Genes in " * italic("Paip1") * "-cKO"
    )
  )

print(p_kegg)
dev.off()

# 11. GO 基因-通路网络图
ego_bp_down <- enrichGO(
  gene = down_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)

foldchanges <- res_df$log2FoldChange
names(foldchanges) <- res_df$gene_id

pdf("Fig9_GO_Cnetplot_Network.pdf", width = 12, height = 9)

p_cnet <- cnetplot(
  ego_bp_down,
  categorySize = "pvalue",
  foldChange = foldchanges,
  showCategory = 5,
  colorEdge = TRUE,
  circular = FALSE
) +
  scale_color_gradient2(
    name = "log2FC",
    low = "navy",
    high = "firebrick3",
    mid = "white"
  ) +
  labs(title = "Gene-Concept Network (Top Biological Processes)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

print(p_cnet)
dev.off()

# 12. Top 100 差异基因热图
top_degs <- res_df %>%
  filter(Status != "NS") %>%
  arrange(padj) %>%
  head(100) %>%
  pull(gene_id)

deg_matrix <- vst_mat[top_degs, ]
deg_matrix_ordered <- deg_matrix[, ordered_cols]

pdf("Fig10_Global_DEG_Heatmap.pdf", width = 7, height = 9)

pheatmap(
  deg_matrix_ordered,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = 3,
  scale = "row",
  show_rownames = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  annotation_col = colData,
  annotation_colors = ann_colors,
  main = "Top 100 Differentially Expressed Genes"
)

dev.off()

# 13. GO 富集结果导出
write.csv(
  as.data.frame(ego_down),
  file = "Table_GO_Enrichment_Down_All.csv",
  row.names = FALSE
)

up_genes <- res_df %>%
  filter(Status == "Up") %>%
  pull(gene_id)

if (length(up_genes) > 0) {
  up_entrez <- bitr(
    up_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  ego_up <- enrichGO(
    gene = up_entrez$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  write.csv(
    as.data.frame(ego_up),
    file = "Table_GO_Enrichment_Up_All.csv",
    row.names = FALSE
  )
}

# 14. 靶向 GO 双向棒棒糖图
target_down_ids <- c(
  "GO:0051321",
  "GO:0007127",
  "GO:0140013",
  "GO:0045132",
  "GO:0007286",
  "GO:0048515",
  "GO:0120316",
  "GO:0035082",
  "GO:0003341"
)

target_up_ids <- c(
  "GO:1901317",
  "GO:2000241",
  "GO:0000979",
  "GO:0032102"
)

if (file.exists("Table_GO_Enrichment_Up_All.csv")) {
  df_down <- read.csv("Table_GO_Enrichment_Down_All.csv")
  df_up <- read.csv("Table_GO_Enrichment_Up_All.csv")

  plot_down <- df_down %>%
    filter(ID %in% target_down_ids) %>%
    mutate(
      Direction = "Down-regulated",
      LogP = log10(p.adjust)
    )

  plot_up <- df_up %>%
    filter(ID %in% target_up_ids) %>%
    mutate(
      Direction = "Up-regulated",
      LogP = -log10(p.adjust)
    )

  plot_df <- bind_rows(plot_down, plot_up) %>%
    mutate(
      Description = case_when(
        Description == "RNA polymerase II core promoter sequence-specific DNA binding" ~ "RNAP II promoter DNA binding",
        Description == "negative regulation of response to external stimulus" ~ "neg. regulation of response to stimulus",
        Description == "regulation of flagellated sperm motility" ~ "regulation of sperm motility",
        Description == "regulation of reproductive process" ~ "regulation of reproduction",
        TRUE ~ Description
      )
    )

  ordered_descriptions <- c(
    "meiotic cell cycle",
    "meiosis I",
    "meiotic nuclear division",
    "meiotic chromosome segregation",
    "spermatid development",
    "spermatid differentiation",
    "sperm flagellum assembly",
    "axoneme assembly",
    "cilium movement",
    "regulation of sperm motility",
    "regulation of reproduction",
    "RNAP II promoter DNA binding",
    "neg. regulation of response to stimulus"
  )

  plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(ordered_descriptions)
  )

  pdf("Fig11_Targeted_GO_Diverging_Lollipop.pdf", width = 10.5, height = 7.5)

  p_lollipop <- ggplot(
    plot_df,
    aes(x = LogP, y = Description, color = Direction)
  ) +
    geom_segment(
      aes(x = 0, xend = LogP, y = Description, yend = Description),
      linewidth = 1.2,
      alpha = 0.8
    ) +
    geom_point(aes(size = Count), alpha = 1) +
    scale_color_manual(
      values = c(
        "Down-regulated" = "#f933ff",
        "Up-regulated" = "#2591a8"
      )
    ) +
    scale_size_continuous(range = c(4, 10), name = "Gene Count") +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.8
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold", color = "black", size = 13),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title.x = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
    ) +
    labs(
      title = expression("Targeted Biological Processes in " * italic("Paip1") * "-cKO Testis"),
      x = expression(Signed ~ -Log[10](italic(P)[adj])),
      y = ""
    ) +
    scale_x_continuous(labels = function(x) abs(x)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 35))

  print(p_lollipop)
  dev.off()
}

cat("Ordinary RNA-seq analysis completed.\n")

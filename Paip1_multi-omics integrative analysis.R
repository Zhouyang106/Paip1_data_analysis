# ==============================================================================
# PAIP1 多组学联合分析
# ==============================================================================

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(pheatmap)
library(ggVennDiagram)
library(igraph)
library(ggraph)
library(stringr)
library(tibble)

setwd("D:/LK/Proteome/")

color_ctrl <- "#f933ff"
color_cko <- "#2591a8"

fc_threshold <- 0.58
p_threshold <- 0.05

# ------------------------------------------------------------------------------
# 1. 读取转录组和蛋白组差异结果
# ------------------------------------------------------------------------------

prot_res <- read.csv("Table_Proteome_DEPs.csv")

rna_res <- read.csv(
  "D:/LK/普通转录组/PAIP1_Optimized_DESeq2_Results.csv"
)

if ("gene_id" %in% colnames(rna_res)) {
  rna_res <- rna_res %>% rename(Gene = gene_id)
} else if ("X" %in% colnames(rna_res)) {
  rna_res <- rna_res %>% rename(Gene = X)
} else if ("row.names" %in% colnames(rna_res)) {
  rna_res <- rna_res %>% rename(Gene = row.names)
}

rna_res <- rna_res %>%
  rename_with(
    ~ "log2FC_RNA",
    .cols = any_of("log2FoldChange")
  ) %>%
  rename_with(
    ~ "pval_RNA",
    .cols = any_of("pvalue")
  ) %>%
  rename_with(
    ~ "padj_RNA",
    .cols = any_of("padj")
  )

if (!"Status_RNA" %in% colnames(rna_res)) {
  rna_res <- rna_res %>%
    mutate(
      Status_RNA = case_when(
        !is.na(padj_RNA) &
          padj_RNA < p_threshold &
          log2FC_RNA > fc_threshold ~ "Up",
        !is.na(padj_RNA) &
          padj_RNA < p_threshold &
          log2FC_RNA < -fc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )
}

multiomics <- inner_join(
  rna_res,
  prot_res,
  by = "Gene"
)

write.csv(
  multiomics,
  "Table_Multiomics_All.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 2. 九象限联合分析
# ------------------------------------------------------------------------------

multiomics <- multiomics %>%
  mutate(
    Sig_RNA = case_when(
      !is.na(pval_RNA) &
        pval_RNA < p_threshold &
        abs(log2FC_RNA) > fc_threshold &
        log2FC_RNA > 0 ~ "Up",
      !is.na(pval_RNA) &
        pval_RNA < p_threshold &
        abs(log2FC_RNA) > fc_threshold &
        log2FC_RNA < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    Sig_Prot = case_when(
      !is.na(pval_Prot) &
        pval_Prot < p_threshold &
        abs(log2FC_Prot) > fc_threshold &
        log2FC_Prot > 0 ~ "Up",
      !is.na(pval_Prot) &
        pval_Prot < p_threshold &
        abs(log2FC_Prot) > fc_threshold &
        log2FC_Prot < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    Quadrant = case_when(
      Sig_RNA == "Up" & Sig_Prot == "Up" ~ "Q1",
      Sig_RNA == "NS" & Sig_Prot == "Up" ~ "Q2",
      Sig_RNA == "Down" & Sig_Prot == "Up" ~ "Q3",
      Sig_RNA == "Up" & Sig_Prot == "NS" ~ "Q4",
      Sig_RNA == "NS" & Sig_Prot == "NS" ~ "Q5",
      Sig_RNA == "Down" & Sig_Prot == "NS" ~ "Q6",
      Sig_RNA == "Up" & Sig_Prot == "Down" ~ "Q7",
      Sig_RNA == "NS" & Sig_Prot == "Down" ~ "Q8",
      Sig_RNA == "Down" & Sig_Prot == "Down" ~ "Q9",
      TRUE ~ "Other"
    ),
    Label = ifelse(
      Gene %in% c("Ankrd31", "Spdya", "Ercc8", "Esco1", "Myo19"),
      Gene,
      NA
    )
  )

write.csv(
  multiomics,
  "Table_Multiomics_Quadrants.csv",
  row.names = FALSE
)

uncoupling_targets <- multiomics %>%
  filter(Sig_RNA %in% c("NS", "Up"), Sig_Prot == "Down") %>%
  arrange(pval_Prot) %>%
  slice_head(n = 15)

pdf(
  "Fig_Multiomics_Nine_Quadrant_Plot.pdf",
  width = 9,
  height = 8
)

ggplot(
  multiomics,
  aes(log2FC_RNA, log2FC_Prot, color = Quadrant)
) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(
    values = c(
      "Q1" = "#e64b35",
      "Q2" = "#f39b7f",
      "Q3" = "#8491b4",
      "Q4" = "#f39b7f",
      "Q5" = "#cccccc",
      "Q6" = "#91d1c2",
      "Q7" = "#DC0000",
      "Q8" = "#DC0000",
      "Q9" = "#4dbbd5",
      "Other" = "#cccccc"
    )
  ) +
  geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "black"
  ) +
  geom_text_repel(
    data = uncoupling_targets,
    aes(label = Gene),
    size = 5,
    fontface = "bold.italic",
    color = "black",
    box.padding = 0.8
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    )
  ) +
  labs(
    title = "Transcriptome-Proteome Uncoupling",
    x = expression(Log[2] * FC ~ "(RNA-seq)"),
    y = expression(Log[2] * FC ~ "(Proteomics)")
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5))

dev.off()

# ------------------------------------------------------------------------------
# 3. 蛋白下调分层分析
# ------------------------------------------------------------------------------

protein_down_targets <- multiomics %>%
  filter(Sig_Prot == "Down") %>%
  mutate(
    Target_Group = "Proteome Down"
  )

write.csv(
  protein_down_targets,
  "Table_Multiomics_Proteome_Down_Targets.csv",
  row.names = FALSE
)

pdf(
  "Fig_Combined_Down_Scatter.pdf",
  width = 8,
  height = 7
)

ggplot(multiomics, aes(log2FC_RNA, log2FC_Prot)) +
  geom_point(
    color = "grey85",
    size = 2.5,
    alpha = 0.6
  ) +
  geom_point(
    data = protein_down_targets,
    color = color_cko,
    size = 3.5,
    alpha = 0.85
  ) +
  geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_hline(
    yintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "grey30"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 22
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Proteome Down-regulation Targets",
    x = expression(Log[2] * FC ~ "(RNA-seq)"),
    y = expression(Log[2] * FC ~ "(Proteomics)")
  )

dev.off()

# ------------------------------------------------------------------------------
# 4. 蛋白下调靶点 GO 富集
# ------------------------------------------------------------------------------

target_genes <- unique(protein_down_targets$Gene)
universe_genes <- unique(multiomics$Gene)

target_entrez <- bitr(
  target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

universe_entrez <- bitr(
  universe_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

if (length(target_entrez) > 0 && length(universe_entrez) > 0) {
  ego_down_multi <- enrichGO(
    gene = unique(target_entrez),
    universe = unique(universe_entrez),
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "none",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 2,
    readable = TRUE
  )

  down_multi_df <- as.data.frame(ego_down_multi)

  if (nrow(down_multi_df) > 0) {
    target_multi_terms <- down_multi_df %>%
      filter(
        grepl(
          "DNA|repair|chromosome|meio|repro|cell cycle|strand|telomere",
          Description,
          ignore.case = TRUE
        )
      ) %>%
      arrange(pvalue) %>%
      slice_head(n = 10)

    write.csv(
      down_multi_df,
      "Table_Multiomics_Proteome_Down_GO_All.csv",
      row.names = FALSE
    )

    write.csv(
      target_multi_terms,
      "Table_Combined_Down_Targeted_DSB.csv",
      row.names = FALSE
    )

    if (nrow(target_multi_terms) > 0) {
      target_multi_terms <- target_multi_terms %>%
        mutate(log_P = -log10(pvalue))

      target_multi_terms$Description <- factor(
        target_multi_terms$Description,
        levels = target_multi_terms$Description
      )

      pdf(
        "Fig_Combined_Down_GO_Bubble.pdf",
        width = 8,
        height = 7.5
      )

      ggplot(
        target_multi_terms,
        aes(x = "", y = Description)
      ) +
        geom_point(
          aes(size = Count, fill = log_P),
          shape = 21,
          color = "black",
          stroke = 1.2
        ) +
        scale_fill_gradient(
          low = "#2591a8",
          high = "#f933ff",
          name = expression(bold(-Log[10] * P))
        ) +
        scale_size_continuous(
          range = c(5, 12),
          name = "Gene Count"
        ) +
        scale_y_discrete(
          labels = function(x) stringr::str_wrap(x, width = 45)
        ) +
        theme_bw(base_size = 16) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_text(
            size = 15,
            color = "black",
            face = "bold",
            lineheight = 1.1
          ),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 20
          )
        ) +
        labs(
          title = "Targeted Pathway Enrichment (Proteome Down)"
        )

      dev.off()
    }
  }
}

# ------------------------------------------------------------------------------
# 5. 翻译解偶联靶点
# ------------------------------------------------------------------------------

uncoupled_targets <- multiomics %>%
  filter(
    Sig_Prot == "Down",
    Sig_RNA %in% c("NS", "Up")
  )

write.csv(
  uncoupled_targets,
  "Table_Multiomics_Uncoupled_Targets.csv",
  row.names = FALSE
)

uncoupled_entrez <- bitr(
  unique(uncoupled_targets$Gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

if (
  length(uncoupled_entrez) > 0 &&
  length(universe_entrez) > 0
) {
  ego_unc <- enrichGO(
    gene = unique(uncoupled_entrez),
    universe = unique(universe_entrez),
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "none",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 2,
    readable = TRUE
  )

  unc_df <- as.data.frame(ego_unc)

  if (nrow(unc_df) > 0) {
    target_unc_terms <- unc_df %>%
      filter(
        grepl(
          "DNA|repair|chromosome|meio|repro|cell cycle|strand|telomere",
          Description,
          ignore.case = TRUE
        )
      ) %>%
      arrange(pvalue) %>%
      slice_head(n = 8)

    write.csv(
      unc_df,
      "Table_Multiomics_Uncoupled_GO_All.csv",
      row.names = FALSE
    )

    write.csv(
      target_unc_terms,
      "Table_Multiomics_Targeted_DSB.csv",
      row.names = FALSE
    )

    if (nrow(target_unc_terms) > 0) {
      target_unc_terms <- target_unc_terms %>%
        mutate(LogP = -log10(pvalue)) %>%
        mutate(Description = stringr::str_to_sentence(Description))

      target_unc_terms$Description <- factor(
        target_unc_terms$Description,
        levels = rev(target_unc_terms$Description)
      )

      pdf(
        "Fig_Multiomics_DSB_Enrichment.pdf",
        width = 9.5,
        height = 6
      )

      ggplot(
        target_unc_terms,
        aes(FoldEnrichment, Description)
      ) +
        geom_segment(
          aes(
            x = 0,
            xend = FoldEnrichment,
            y = Description,
            yend = Description
          ),
          color = "grey70",
          linewidth = 1.2
        ) +
        geom_point(
          aes(size = Count, fill = LogP),
          shape = 21,
          color = "black",
          stroke = 1.2
        ) +
        scale_size_continuous(
          range = c(6, 12),
          name = "Gene Count"
        ) +
        scale_fill_gradient(
          low = "#4DBBD5",
          high = "#E64B35",
          name = expression(-Log[10](italic(P)-value))
        ) +
        theme_bw(base_size = 14) +
        theme(
          axis.text.y = element_text(
            color = "black",
            face = "bold",
            size = 11
          ),
          axis.title = element_text(face = "bold"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(
            hjust = 0.5,
            face = "bold"
          )
        ) +
        labs(
          title = "Functional Enrichment of Uncoupled Targets",
          x = "Fold Enrichment",
          y = ""
        )

      dev.off()
    }
  }
}

# ------------------------------------------------------------------------------
# 6. 机制韦恩图
# ------------------------------------------------------------------------------

rna_down_list <- multiomics %>%
  filter(Sig_RNA == "Down") %>%
  pull(Gene)

prot_down_list <- multiomics %>%
  filter(Sig_Prot == "Down") %>%
  pull(Gene)

pdf("Fig_Venn_Co_Downregulated.pdf", width = 6, height = 5)

p1 <- ggVennDiagram(
  list(
    "RNA Down" = unique(rna_down_list),
    "Protein Down" = unique(prot_down_list)
  ),
  label_alpha = 0,
  edge_size = 1,
  set_color = "black"
) +
  scale_fill_gradient(
    low = "#F4F6F6",
    high = "#5DADE2"
  ) +
  theme(legend.position = "none") +
  labs(title = "Co-downregulated Targets") +
  theme(plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 16
  ))

print(p1)
dev.off()

rna_normal_up_list <- multiomics %>%
  filter(Sig_RNA %in% c("NS", "Up")) %>%
  pull(Gene)

pdf("Fig_Venn_Uncoupled.pdf", width = 6, height = 5)

p2 <- ggVennDiagram(
  list(
    "RNA Unchanged/Up" = unique(rna_normal_up_list),
    "Protein Down" = unique(prot_down_list)
  ),
  label_alpha = 0,
  edge_size = 1,
  set_color = "black"
) +
  scale_fill_gradient(
    low = "#FDEDEC",
    high = "#E74C3C"
  ) +
  theme(legend.position = "none") +
  labs(title = "Translationally Uncoupled Targets") +
  theme(plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 16
  ))

print(p2)
dev.off()

# ------------------------------------------------------------------------------
# 7. 多组学机制网络图
# ------------------------------------------------------------------------------

custom_links <- data.frame(
  Term = c(
    rep("Meiosis Arrest & DNA Repair", 4),
    rep("Flagella Assembly & Cytoskeleton", 4),
    rep("Glycolysis & Energy Metabolism", 4)
  ),
  Gene = c(
    "Ankrd31", "Spdya", "Esco1", "Ercc8",
    "Myh7", "Actn3", "Myo19", "Cep83",
    "Eno3", "Pgam2", "Prkaa2", "Atp2a1"
  )
)

node_fc <- custom_links %>%
  left_join(
    prot_res %>% select(Gene, log2FC_Prot),
    by = "Gene"
  )

terms <- unique(custom_links$Term)

term_nodes <- data.frame(
  name = terms,
  log2FC_Prot = NA_real_,
  node_type = "Pathway"
)

gene_nodes <- node_fc %>%
  transmute(
    name = Gene,
    log2FC_Prot = log2FC_Prot,
    node_type = "Gene"
  )

all_nodes <- bind_rows(
  term_nodes,
  gene_nodes
)

net <- graph_from_data_frame(
  d = custom_links,
  vertices = all_nodes,
  directed = FALSE
)

set.seed(42)

pdf(
  "Fig_Publication_Mechanism_Network.pdf",
  width = 12,
  height = 9
)

ggraph(net, layout = "kk") +
  geom_edge_link(
    alpha = 0.5,
    color = "grey60",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_node_point(
    data = function(x) x %>% filter(node_type == "Gene"),
    aes(size = 8, color = log2FC_Prot)
  ) +
  scale_color_gradient2(
    low = "#E64B35",
    mid = "#4DBBD5",
    high = "grey80",
    midpoint = -1.5,
    name = expression(Log[2] * FC)
  ) +
  geom_node_point(
    data = function(x) x %>% filter(node_type == "Pathway"),
    aes(size = 15),
    color = "#3C5488",
    shape = 18
  ) +
  geom_node_text(
    data = function(x) x %>% filter(node_type == "Gene"),
    aes(label = name),
    repel = TRUE,
    size = 5,
    fontface = "bold.italic",
    color = "black",
    box.padding = 0.8,
    point.padding = 0.5,
    max.overlaps = Inf
  ) +
  geom_node_label(
    data = function(x) x %>% filter(node_type == "Pathway"),
    aes(label = name),
    repel = TRUE,
    size = 5.5,
    fontface = "bold",
    color = "white",
    fill = "#3C5488",
    box.padding = 1,
    point.padding = 0.5,
    label.padding = unit(0.4, "lines"),
    label.r = unit(0.3, "lines")
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    legend.position = "right"
  ) +
  labs(title = "Multi-omics Integrated Mechanism Network") +
  scale_size_identity()

dev.off()

# ------------------------------------------------------------------------------
# 8. 靶向蛋白表达图
# ------------------------------------------------------------------------------

target_protein_genes <- data.frame(
  Gene = c(
    "Ankrd31", "Spdya", "Esco1", "Ercc8",
    "Myh7", "Actn3", "Myo19", "Cep83",
    "Eno3", "Pgam2", "Prkaa2", "Atp2a1"
  ),
  Category = c(
    rep("Meiosis & DNA Repair", 4),
    rep("Flagella & Cytoskeleton", 4),
    rep("Energy Metabolism", 4)
  )
)

target_protein_plot <- target_protein_genes %>%
  left_join(prot_res, by = "Gene") %>%
  mutate(
    sig_label = case_when(
      pval_Prot < 0.001 ~ "***",
      pval_Prot < 0.01 ~ "**",
      pval_Prot < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(Category, log2FC_Prot)

target_protein_plot$Gene <- factor(
  target_protein_plot$Gene,
  levels = target_protein_plot$Gene
)

pdf(
  "Fig_Target_Proteins_Lollipop.pdf",
  width = 9,
  height = 6
)

ggplot(
  target_protein_plot,
  aes(Gene, log2FC_Prot, color = Category)
) +
  geom_segment(
    aes(
      x = Gene,
      xend = Gene,
      y = 0,
      yend = log2FC_Prot
    ),
    linewidth = 1.2
  ) +
  geom_point(size = 5) +
  geom_text(
    aes(
      label = sig_label,
      y = log2FC_Prot - 0.15
    ),
    color = "black",
    size = 5,
    fontface = "bold"
  ) +
  scale_color_manual(
    values = c(
      "Meiosis & DNA Repair" = "#E74C3C",
      "Flagella & Cytoskeleton" = "#3498DB",
      "Energy Metabolism" = "#2ECC71"
    )
  ) +
  facet_grid(
    ~ Category,
    scales = "free_x",
    space = "free_x"
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.8
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold.italic",
      color = "black"
    ),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.05))
  ) +
  labs(
    title = "Protein Changes of Core Targets",
    x = "",
    y = expression(Log[2] * FC)
  )

dev.off()

# ------------------------------------------------------------------------------
# 9. RNA-seq 与蛋白组样本表达矩阵重建
# ------------------------------------------------------------------------------

rna_counts_file <- "D:/LK/普通转录组/quant/PAIP1_counts.txt"

raw_counts <- read.table(
  rna_counts_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

rna_count_matrix <- raw_counts[, 7:ncol(raw_counts)]
colnames(rna_count_matrix) <- c(
  "cKO_1", "cKO_2", "cKO_3",
  "Ctrl_1", "Ctrl_2", "Ctrl_3"
)

clean_ensembl <- gsub(
  "\\..*",
  "",
  raw_counts$Geneid
)

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = clean_ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols[is.na(gene_symbols)] <- clean_ensembl[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)

rownames(rna_count_matrix) <- gene_symbols

colData_rna <- data.frame(
  row.names = colnames(rna_count_matrix),
  condition = factor(
    c(rep("cKO", 3), rep("Control", 3)),
    levels = c("Control", "cKO")
  )
)

dds_rna <- DESeqDataSetFromMatrix(
  countData = rna_count_matrix,
  colData = colData_rna,
  design = ~ condition
)

keep_rna <- rowSums(counts(dds_rna) >= 10) >= 3
dds_rna <- dds_rna[keep_rna, ]
dds_rna <- DESeq(dds_rna)

vsd_rna <- vst(dds_rna, blind = FALSE)
rna_vst <- assay(vsd_rna)

# ------------------------------------------------------------------------------
# 10. RNA-Protein 靶基因热图与相关性分析
# ------------------------------------------------------------------------------

key_target_genes <- c(
  "Paip1", "Sycp1", "Sycp3", "Hormad1", "Mlh1",
  "Spo11", "Atr", "H2ax", "Rpa2", "Rad51", "Dmc1"
)

prot_mat <- log2_expr_filtered
prot_mat <- prot_mat[, sample_names, drop = FALSE]

target_rna_genes <- intersect(
  key_target_genes,
  rownames(rna_vst)
)

target_prot_genes <- intersect(
  key_target_genes,
  rownames(prot_mat)
)

common_genes <- intersect(
  target_rna_genes,
  target_prot_genes
)

if (length(common_genes) > 0) {
  rna_target_matrix <- rna_vst[
    common_genes,
    sample_names,
    drop = FALSE
  ]

  prot_target_matrix <- prot_mat[
    common_genes,
    sample_names,
    drop = FALSE
  ]

  pdf(
    "Fig5_Targeted_Proteome_Heatmap.pdf",
    width = 6.5,
    height = 7
  )

  pheatmap(
    prot_target_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col = 3,
    scale = "row",
    color = colorRampPalette(
      c("navy", "white", "firebrick3")
    )(50),
    annotation_col = data.frame(
      condition = factor(
        c(
          "Control", "Control", "Control",
          "cKO", "cKO", "cKO"
        ),
        levels = c("Control", "cKO")
      ),
      row.names = sample_names
    ),
    annotation_colors = list(
      condition = c(
        "Control" = color_ctrl,
        "cKO" = color_cko
      )
    ),
    fontsize = 14,
    main = "Key Proteins Expression"
  )

  dev.off()

  cor_matrix <- cor(
    t(rna_target_matrix),
    t(prot_target_matrix),
    method = "pearson"
  )

  p_matrix <- matrix(
    NA_real_,
    nrow = length(common_genes),
    ncol = length(common_genes),
    dimnames = list(common_genes, common_genes)
  )

  for (i in seq_along(common_genes)) {
    for (j in seq_along(common_genes)) {
      test_res <- cor.test(
        as.numeric(rna_target_matrix[i, ]),
        as.numeric(prot_target_matrix[j, ]),
        method = "pearson"
      )
      p_matrix[i, j] <- test_res$p.value
    }
  }

  star_matrix <- matrix(
    "",
    nrow = nrow(p_matrix),
    ncol = ncol(p_matrix),
    dimnames = dimnames(p_matrix)
  )

  star_matrix[p_matrix < 0.05] <- "*"
  star_matrix[p_matrix < 0.01] <- "**"
  star_matrix[p_matrix < 0.001] <- "***"

  write.csv(
    cor_matrix,
    "Table_RNA_Protein_Correlation.csv"
  )

  write.csv(
    p_matrix,
    "Table_RNA_Protein_Correlation_Pvalue.csv"
  )

  pdf(
    "Fig6_RNA_Protein_Correlation_Stars.pdf",
    width = 7.5,
    height = 7.5
  )

  pheatmap(
    cor_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = star_matrix,
    fontsize_number = 18,
    number_color = "black",
    breaks = seq(-1, 1, length.out = 51),
    color = colorRampPalette(
      c("#3498DB", "white", "#E74C3C")
    )(50),
    fontsize = 14,
    main = "Cross-omics Correlation"
  )

  dev.off()
}

cat("\nMulti-omics analysis completed.\n")

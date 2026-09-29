# ==============================================================================
# PAIP1-cKO 定量蛋白组分析
# ==============================================================================

library(tidyverse)
library(limma)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(enrichplot)
library(stringr)
library(tibble)

setwd("D:/LK/Proteome/")

color_ctrl <- "#f933ff"
color_cko <- "#2591a8"

expr_cols <- c("WT-1", "WT-2", "WT-3", "CKO-1", "CKO-2", "CKO-3")
sample_names <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "cKO_1", "cKO_2", "cKO_3")
group_levels <- c("Control", "cKO")

# ------------------------------------------------------------------------------
# 1. 数据读取与预处理
# ------------------------------------------------------------------------------

prot_raw <- read.delim(
  "data/all.proteins.normalized.xls",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

prot_raw <- prot_raw %>%
  filter(!is.na(Gene), Gene != "") %>%
  mutate(
    mean_expr = rowMeans(
      dplyr::select(., all_of(expr_cols)),
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(mean_expr)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  dplyr::select(-mean_expr)

expr_matrix <- as.matrix(prot_raw[, expr_cols])
rownames(expr_matrix) <- prot_raw$Gene

log2_expr <- log2(expr_matrix + 1)
keep <- rowSums(log2_expr > 1) >= 3
log2_expr_filtered <- log2_expr[keep, , drop = FALSE]

colnames(log2_expr_filtered) <- sample_names

colData_prot <- data.frame(
  condition = factor(
    c(rep("Control", 3), rep("cKO", 3)),
    levels = group_levels
)
)
rownames(colData_prot) <- sample_names

# ------------------------------------------------------------------------------
# 2. PCA
# ------------------------------------------------------------------------------

pca_res <- prcomp(t(log2_expr_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
pca_df$Group <- factor(
  c(rep("Control", 3), rep("cKO", 3)),
  levels = group_levels
)
pca_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pdf("Fig1_Proteome_PCA.pdf", width = 7.5, height = 6)

ggplot(pca_df, aes(PC1, PC2, color = Group, fill = Group)) +
  stat_ellipse(
    geom = "polygon",
    level = 0.95,
    type = "norm",
    alpha = 0.15,
    linetype = "dashed",
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_point(
    size = 4,
    alpha = 0.9,
    shape = 21,
    color = "black",
    stroke = 1
  ) +
  geom_text_repel(
    aes(label = Sample),
    size = 4.5,
    fontface = "bold",
    color = "black",
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Control" = color_ctrl, "cKO" = color_cko)) +
  scale_fill_manual(values = c("Control" = color_ctrl, "cKO" = color_cko)) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  scale_y_continuous(expand = expansion(mult = 0.2)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", linewidth = 1)
  ) +
  labs(
    title = "Proteome PCA",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  )

dev.off()

# ------------------------------------------------------------------------------
# 3. Limma 差异蛋白分析
# ------------------------------------------------------------------------------

group_list <- factor(
  c(rep("Control", 3), rep("cKO", 3)),
  levels = group_levels
)

design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)

fit <- lmFit(log2_expr_filtered, design)
contrast_matrix <- makeContrasts(cKO - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

prot_res <- topTable(fit2, coef = 1, number = Inf) %>%
  tibble::rownames_to_column("Gene") %>%
  rename(
    log2FC_Prot = logFC,
    pval_Prot = P.Value,
    padj_Prot = adj.P.Val
  )

log2fc_cutoff_prot <- 0.58
pval_cutoff_prot <- 0.05

prot_res <- prot_res %>%
  mutate(
    Status_Prot = case_when(
      pval_Prot < pval_cutoff_prot & log2FC_Prot > log2fc_cutoff_prot ~ "Up",
      pval_Prot < pval_cutoff_prot & log2FC_Prot < -log2fc_cutoff_prot ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(
  prot_res,
  "Table_Proteome_DEPs.csv",
  row.names = FALSE
)

write.csv(
  prot_res,
  "PAIP1_Proteome_limma_Results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. 差异蛋白火山图
# ------------------------------------------------------------------------------

count_up <- sum(prot_res$Status_Prot == "Up")
count_down <- sum(prot_res$Status_Prot == "Down")

star_genes <- c(
  "Paip1", "Dazl", "Piwil1", "Piwil2", "Sycp3",
  "Boule", "Prm1", "Prm2", "Tnp1", "Tnp2",
  "Pgk2", "Acr", "Spag1", "Dmc1", "Rad51",
  "Stra8", "Mlh1"
)

label_prot <- prot_res %>%
  filter(Gene %in% star_genes)

pdf("Fig2_Proteome_Volcano_Targeted.pdf", width = 8, height = 7)

ggplot(prot_res, aes(log2FC_Prot, -log10(pval_Prot), color = Status_Prot)) +
  geom_point(alpha = 0.55, size = 3) +
  scale_color_manual(
    values = c(
      "Up" = color_cko,
      "Down" = color_ctrl,
      "NS" = "#E5E5E5"
    )
  ) +
  geom_vline(
    xintercept = c(-log2fc_cutoff_prot, log2fc_cutoff_prot),
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_hline(
    yintercept = -log10(pval_cutoff_prot),
    linetype = "dashed",
    color = "grey30"
  ) +
  geom_text_repel(
    data = label_prot,
    aes(label = Gene),
    size = 5,
    fontface = "bold.italic",
    color = "black",
    box.padding = 0.7,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = paste0("Down: ", count_down),
    hjust = -0.1,
    vjust = 1.6,
    color = color_ctrl,
    size = 5.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Up: ", count_up),
    hjust = 1.1,
    vjust = 1.6,
    color = color_cko,
    size = 5.5,
    fontface = "bold"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = expression("Proteome (" * italic("Paip1") * "-cKO vs Control)"),
    x = expression(Log[2] * FC),
    y = expression(-Log[10] * P)
  )

dev.off()

# ------------------------------------------------------------------------------
# 5. 下调蛋白 GO 富集
# ------------------------------------------------------------------------------

down_prots <- prot_res %>%
  filter(Status_Prot == "Down") %>%
  pull(Gene)

if (length(down_prots) > 0) {
  entrez_down <- bitr(
    unique(down_prots),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  if (nrow(entrez_down) > 0) {
    ego_prot_down <- enrichGO(
      gene = unique(entrez_down$ENTREZID),
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable = TRUE
    )

    go_down_df <- as.data.frame(ego_prot_down)

    if (nrow(go_down_df) > 0) {
      write.csv(
        go_down_df,
        "Table_Proteome_Down_GO.csv",
        row.names = FALSE
      )

      pdf("Fig3_Proteome_Down_GO_Barplot.pdf", width = 8, height = 6)
      print(
        barplot(
          ego_prot_down,
          showCategory = 15,
          title = "Down-regulated Proteins"
        ) +
          theme_classic(base_size = 14) +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      )
      dev.off()
    }
  }
}

# ------------------------------------------------------------------------------
# 6. 下调蛋白 GO 客观分类展示
# ------------------------------------------------------------------------------

if (file.exists("Table_Proteome_Down_GO.csv")) {
  go_res <- read.csv("Table_Proteome_Down_GO.csv")

  selected_terms <- data.frame(
    Description = c(
      "multicellular organismal movement",
      "striated muscle contraction",
      "sarcomere organization",
      "vesicle docking",
      "glucose catabolic process",
      "glycolytic process",
      "pyruvate metabolic process",
      "ATP metabolic process",
      "purine ribonucleotide metabolic process",
      "ribose phosphate metabolic process",
      "transcription initiation at RNA polymerase II promoter"
    ),
    Category = c(
      rep("Cytoskeleton & Motility", 4),
      rep("Energy Metabolism", 4),
      rep("Nucleotide & Transcription", 3)
    )
  )

  plot_df <- inner_join(
    selected_terms,
    go_res,
    by = "Description"
  ) %>%
    mutate(LogP = -log10(p.adjust)) %>%
    arrange(Category, LogP)

  if (nrow(plot_df) > 0) {
    plot_df$Description <- factor(
      plot_df$Description,
      levels = plot_df$Description
    )

    plot_df$Category <- factor(
      plot_df$Category,
      levels = c(
        "Cytoskeleton & Motility",
        "Energy Metabolism",
        "Nucleotide & Transcription"
      )
    )

    cat_colors <- c(
      "Cytoskeleton & Motility" = "#E64B35",
      "Energy Metabolism" = "#4DBBD5",
      "Nucleotide & Transcription" = "#00A087"
    )

    pdf(
      "Fig5_Proteome_Objective_GO_Barplot.pdf",
      width = 12,
      height = 10
    )

    ggplot(plot_df, aes(LogP, Description, fill = Category)) +
      geom_col(
        color = "black",
        linewidth = 0.5,
        width = 0.7,
        alpha = 0.85
      ) +
      scale_fill_manual(values = cat_colors) +
      facet_grid(
        Category ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      geom_vline(
        xintercept = -log10(0.05),
        linetype = "dashed",
        color = "red",
        linewidth = 0.8
      ) +
      theme_bw(base_size = 15) +
      theme(
        strip.text.y = element_text(angle = 270, face = "bold", size = 13),
        strip.background = element_rect(fill = "grey90", color = "black"),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(face = "bold", size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        title = "Key Functional Modules of Down-regulated Proteins",
        x = expression(-Log[10](italic(P)[adj])),
        y = ""
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

    dev.off()
  }
}

# ------------------------------------------------------------------------------
# 7. 靶向基因蛋白表达箱线图
# ------------------------------------------------------------------------------

genes_for_boxplot <- c(
  "Paip1", "Ybx2", "Rad51", "Dmc1", "Rpa2", "Rec8",
  "Stra8", "Sycp1", "Sycp3", "Mlh1", "Hormad1",
  "Spo11", "Atr", "H2ax"
)

prot_sub <- as.data.frame(log2_expr_filtered) %>%
  tibble::rownames_to_column("Gene") %>%
  filter(tolower(Gene) %in% tolower(genes_for_boxplot))

detected_genes <- prot_sub$Gene
missing_genes <- genes_for_boxplot[
  !tolower(genes_for_boxplot) %in% tolower(detected_genes)
]

cat(
  "Detected target proteins:",
  length(detected_genes), "/", length(genes_for_boxplot), "\n"
)

if (length(missing_genes) > 0) {
  cat("Missing target proteins:", paste(missing_genes, collapse = ", "), "\n")
}

prot_long <- prot_sub %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Log2_Expression"
  ) %>%
  mutate(
    Group = case_when(
      str_detect(Sample, "^Ctrl") ~ "Control",
      str_detect(Sample, "^cKO") ~ "Paip1-cKO",
      TRUE ~ "Other"
    ),
    Group = factor(Group, levels = c("Control", "Paip1-cKO")),
    Gene = factor(
      Gene,
      levels = genes_for_boxplot[
        tolower(genes_for_boxplot) %in% tolower(Gene)
      ]
    )
  ) %>%
  filter(Group != "Other")

pdf(
  "Fig_Proteome_Target_Genes_Boxplot.pdf",
  width = 12,
  height = 9
)

ggplot(prot_long, aes(Group, Log2_Expression, fill = Group)) +
  geom_boxplot(
    width = 0.5,
    color = "black",
    linewidth = 0.8,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_point(
    position = position_jitter(width = 0.15, seed = 123),
    size = 2.5,
    color = "black",
    alpha = 0.85
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.x.npc = "center",
    label.y.npc = 0.92,
    size = 5
  ) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c(
      "Control" = color_ctrl,
      "Paip1-cKO" = color_cko
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 15, face = "bold.italic"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", angle = 15, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.7, color = "black")
  ) +
  labs(
    y = expression("Log"[2] * " (Protein Expression)"),
    title = "Target Protein Abundance"
  )

dev.off()

# ------------------------------------------------------------------------------
# 8. 差异蛋白全局热图
# ------------------------------------------------------------------------------

sig_genes <- prot_res %>%
  filter(Status_Prot != "NS") %>%
  pull(Gene)

sig_genes <- intersect(sig_genes, rownames(log2_expr_filtered))

if (length(sig_genes) > 0) {
  sig_matrix <- log2_expr_filtered[sig_genes, , drop = FALSE]

  z_matrix <- t(apply(sig_matrix, 1, function(x) {
    s <- sd(x)
    if (is.na(s) || s == 0) {
      rep(0, length(x))
    } else {
      (x - mean(x)) / s
    }
  }))

  pdf("FigB_Proteome_DEPs_Heatmap.pdf", width = 6, height = 8)

  pheatmap(
    z_matrix,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    show_rownames = FALSE,
    cluster_cols = FALSE,
    main = "Differentially Expressed Proteins",
    fontsize = 12,
    annotation_col = colData_prot,
    annotation_colors = list(
      condition = c("Control" = color_ctrl, "cKO" = color_cko)
    )
  )

  dev.off()
}

# ------------------------------------------------------------------------------
# 9. 重点蛋白展示
# ------------------------------------------------------------------------------

focused_genes <- c("Paip1", "H2ax", "Rad51", "Dmc1", "Sycp3")

focus_df <- prot_res %>%
  filter(Gene %in% focused_genes)

if (nrow(focus_df) > 0) {
  pdf("FigC_Focused_Genes_Barplot.pdf", width = 7, height = 5)

  ggplot(focus_df, aes(Gene, log2FC_Prot, fill = Status_Prot)) +
    geom_col(color = "black", width = 0.6) +
    scale_fill_manual(
      values = c(
        "Up" = color_cko,
        "Down" = color_ctrl,
        "NS" = "grey80"
      )
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        face = "bold.italic",
        angle = 45,
        hjust = 1,
        color = "black"
      ),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = "Key Protein Changes",
      x = "",
      y = "Log2 Fold Change"
    )

  dev.off()
}

# ------------------------------------------------------------------------------
# 10. 差异蛋白定向 GO/DSB 探索性分析
# ------------------------------------------------------------------------------

universe_genes <- unique(prot_res$Gene)
deg_genes <- prot_res %>%
  filter(Status_Prot != "NS") %>%
  pull(Gene)

universe_entrez <- bitr(
  universe_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

deg_entrez <- bitr(
  unique(deg_genes),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)$ENTREZID

if (length(deg_entrez) > 0 && length(universe_entrez) > 0) {
  ego_dys <- enrichGO(
    gene = unique(deg_entrez),
    universe = unique(universe_entrez),
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "none",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 3,
    readable = TRUE
  )

  target_terms <- as.data.frame(ego_dys) %>%
    filter(
      grepl(
        "DNA|repair|chromosome|meio|repro|cell cycle|strand",
        Description,
        ignore.case = TRUE
      )
    ) %>%
    arrange(pvalue)

  write.csv(
    target_terms,
    "Table_Targeted_DSB_Pathways.csv",
    row.names = FALSE
  )

  if (nrow(target_terms) > 0) {
    plot_df <- target_terms %>%
      head(8) %>%
      mutate(
        LogP = -log10(pvalue),
        Description = stringr::str_to_sentence(Description)
      )

    plot_df$Description <- factor(
      plot_df$Description,
      levels = rev(plot_df$Description)
    )

    pdf(
      "Fig_Proteome_Targeted_DSB_Pathways.pdf",
      width = 10,
      height = 6
    )

    ggplot(plot_df, aes(LogP, Description, fill = LogP)) +
      geom_col(
        color = "black",
        linewidth = 0.6,
        width = 0.7
      ) +
      scale_fill_gradient(
        low = "#2591a8",
        high = "#f933ff",
        name = expression(-Log[10](italic(P)-value))
      ) +
      geom_vline(
        xintercept = -log10(0.05),
        linetype = "dashed",
        color = "grey30"
      ) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.y = element_text(
          color = "black",
          face = "bold",
          size = 11
        ),
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      ) +
      labs(
        title = "Targeted DSB- and Meiosis-related Pathways",
        x = expression(-Log[10](italic(P)-value)),
        y = ""
      )

    dev.off()

    pdf(
      "Fig_Proteome_Targeted_DSB_Pathways_Bubble.pdf",
      width = 8,
      height = 6.5
    )

    ggplot(plot_df, aes(x = "", y = Description)) +
      geom_point(
        aes(size = Count, fill = LogP),
        shape = 21,
        color = "black",
        stroke = 1.2
      ) +
      scale_fill_gradient(
        low = "#2591a8",
        high = "#f933ff",
        name = expression(bold(-Log[10] * italic(P)))
      ) +
      scale_size_continuous(range = c(5, 12), name = "Gene Count") +
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
          size = 14,
          color = "black",
          face = "bold"
        ),
        panel.grid = element_blank(),
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      ) +
      labs(title = "Targeted DSB- and Meiosis-related Pathways")

    dev.off()
  }
}

# ------------------------------------------------------------------------------
# 11. 核心靶蛋白棒棒糖图
# ------------------------------------------------------------------------------

target_genes <- data.frame(
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

target_plot_df <- target_genes %>%
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

target_plot_df$Gene <- factor(
  target_plot_df$Gene,
  levels = target_plot_df$Gene
)

target_plot_df$Category <- factor(
  target_plot_df$Category,
  levels = c(
    "Meiosis & DNA Repair",
    "Flagella & Cytoskeleton",
    "Energy Metabolism"
  )
)

pdf(
  "Fig_Target_Proteins_Lollipop.pdf",
  width = 9,
  height = 6
)

ggplot(
  target_plot_df,
  aes(Gene, log2FC_Prot, color = Category)
) +
  geom_segment(
    aes(
      x = Gene,
      xend = Gene,
      y = 0,
      yend = log2FC_Prot
    ),
    linewidth = 1.2,
    alpha = 0.8
  ) +
  geom_point(size = 5) +
  geom_text(
    aes(
      label = sig_label,
      y = log2FC_Prot - 0.15
    ),
    color = "black",
    size = 5,
    vjust = 0.5,
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
    axis.text.y = element_text(color = "black"),
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

cat("\nProtein analysis completed.\n")

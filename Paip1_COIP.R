# ==============================================================================
# PAIP1 Co-IP/MS 分析
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("D:/LK/COIP/")

# ------------------------------------------------------------------------------
# 1. 原始数据读取与组别校正
# ------------------------------------------------------------------------------

raw_data <- readr::read_tsv("report.pg_matrix.tsv")

clean_data <- raw_data %>%
  dplyr::select(
    Protein.Group,
    Protein.Names,
    Genes,
    First.Protein.Description,
    True_PAIP1_IP = IgG,
    True_IgG_Ctrl = PAIP1
  ) %>%
  mutate(
    True_PAIP1_IP_filled = replace_na(True_PAIP1_IP, 1000),
    True_IgG_Ctrl_filled = replace_na(True_IgG_Ctrl, 1000),
    Fold_Change = True_PAIP1_IP_filled / True_IgG_Ctrl_filled
  ) %>%
  filter(
    !str_detect(
      coalesce(Genes, ""),
      "(?i)^Krt|^Igk|^Igh|^Igl|^Alb"
    )
  )

# ------------------------------------------------------------------------------
# 2. 高可信度互作蛋白筛选
# ------------------------------------------------------------------------------

golden_targets <- clean_data %>%
  filter(
    True_PAIP1_IP_filled > 10000,
    Fold_Change >= 10
  ) %>%
  arrange(desc(Fold_Change))

check_genes <- clean_data %>%
  filter(str_detect(
    coalesce(Genes, ""),
    "(?i)Paip1|Ybx2"
  )) %>%
  select(
    Genes,
    True_IgG_Ctrl,
    True_PAIP1_IP,
    Fold_Change
  )

print(check_genes)

write_csv(
  golden_targets,
  "PAIP1_True_Interactors_High_Confidence.csv"
)

# ------------------------------------------------------------------------------
# 3. CRAPome 风格污染清洗
# ------------------------------------------------------------------------------

interactome_data <- read_csv(
  "PAIP1_True_Interactors_High_Confidence.csv"
)

crapome_pattern <- paste0(
  "(?i)^(Act|Myl|Myh|Tpm|Vim|Calm|S100|",
  "Hba|Hbb|Des$|Tagln|Tmod|H1-|H2[ab]|H3|H4|Hist|",
  "Gapdh|Hsp|Cct|Atp[0-9]|Vdac|Ywha|Ig[kgl]|Krt|Alb)"
)

final_clean_data <- interactome_data %>%
  filter(
    !str_detect(
      coalesce(Genes, ""),
      crapome_pattern
    )
  ) %>%
  arrange(desc(Fold_Change))

check_core <- final_clean_data %>%
  filter(str_detect(
    coalesce(Genes, ""),
    "(?i)Paip1|Ybx2|Piwil|Eif|Rpl|Rps"
  )) %>%
  select(
    Genes,
    True_PAIP1_IP_filled,
    Fold_Change
  ) %>%
  arrange(desc(Fold_Change)) %>%
  slice_head(n = 15)

print(check_core)

write_csv(
  final_clean_data,
  "PAIP1_Clean_Interactome_Final.csv"
)

# ------------------------------------------------------------------------------
# 4. IP-MS 丰度-富集度散点图
# ------------------------------------------------------------------------------

plot_data <- clean_data %>%
  filter(True_PAIP1_IP_filled > 1000) %>%
  mutate(
    log10_Abundance = log10(True_PAIP1_IP_filled),
    log2_FC = log2(Fold_Change),
    Significance = case_when(
      str_detect(
        coalesce(Genes, ""),
        "(?i)^Paip1$|^Ybx2$"
      ) ~ "Bait & Core Targets",
      Fold_Change >= 10 &
        True_PAIP1_IP_filled > 10000 ~ "Significant Interactors",
      TRUE ~ "Background"
    )
  )

top_genes <- plot_data %>%
  filter(Significance != "Background") %>%
  arrange(desc(log2_FC)) %>%
  slice_head(n = 15)

core_genes <- plot_data %>%
  filter(Significance == "Bait & Core Targets")

label_data <- bind_rows(
  top_genes,
  core_genes
) %>%
  distinct(Genes, .keep_all = TRUE)

ip_plot <- ggplot(
  plot_data,
  aes(log10_Abundance, log2_FC)
) +
  geom_point(
    data = filter(plot_data, Significance == "Background"),
    color = "#E0E0E0",
    alpha = 0.6,
    size = 1.5
  ) +
  geom_point(
    data = filter(plot_data, Significance == "Significant Interactors"),
    color = "#4C84C3",
    alpha = 0.8,
    size = 2.5
  ) +
  geom_point(
    data = filter(plot_data, Significance == "Bait & Core Targets"),
    color = "#D3423E",
    alpha = 1,
    size = 3.5
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = Genes),
    size = 4.5,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = 50,
    fontface = "italic"
  ) +
  geom_hline(
    yintercept = log2(10),
    linetype = "dashed",
    color = "grey50"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "PAIP1 Interactome Profiling",
    x = expression(Log[10] * "(PAIP1 IP Abundance)"),
    y = expression(Log[2] * "(Fold Change: IP / IgG)")
  )

ggsave(
  "PAIP1_Interactome_ScatterPlot.pdf",
  ip_plot,
  width = 8,
  height = 6
)

ggsave(
  "PAIP1_Interactome_ScatterPlot.png",
  ip_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------------------------
# 5. 客观分类散点图
# ------------------------------------------------------------------------------

plot_df_objective <- clean_data %>%
  mutate(
    Log2_IgG = log2(True_IgG_Ctrl_filled),
    Log2_PAIP1 = log2(True_PAIP1_IP_filled),
    Log2_FC = log2(Fold_Change),
    Rank = rank(-Fold_Change, ties.method = "first"),
    Point_Category = case_when(
      str_detect(
        coalesce(Genes, ""),
        "(?i)^Paip1$|^Ybx2$"
      ) ~ "Paip1 & Ybx2",
      True_PAIP1_IP_filled > 10000 &
        Fold_Change >= 10 &
        !str_detect(coalesce(Genes, ""), crapome_pattern) ~ "Enriched (FC >= 10)",
      TRUE ~ "Non-enriched"
    ),
    Point_Category = factor(
      Point_Category,
      levels = c(
        "Non-enriched",
        "Enriched (FC >= 10)",
        "Paip1 & Ybx2"
      )
    )
  )

top_n_labels <- 15

top_interactors <- plot_df_objective %>%
  filter(Point_Category == "Enriched (FC >= 10)") %>%
  arrange(desc(Fold_Change)) %>%
  slice_head(n = top_n_labels) %>%
  pull(Genes)

plot_df_objective <- plot_df_objective %>%
  mutate(
    Label = case_when(
      str_detect(
        coalesce(Genes, ""),
        "(?i)^Paip1$|^Ybx2$"
      ) ~ Genes,
      Genes %in% top_interactors ~ Genes,
      TRUE ~ ""
    )
  )

pdf(
  "FigS_IPMS_Scatter_Objective.pdf",
  width = 7.5,
  height = 7
)

ggplot(
  plot_df_objective,
  aes(Log2_IgG, Log2_PAIP1, color = Point_Category)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey50",
    linewidth = 1
  ) +
  geom_point(
    aes(
      size = Point_Category,
      alpha = Point_Category
    )
  ) +
  geom_text_repel(
    aes(label = Label),
    size = 5,
    fontface = "italic",
    box.padding = 0.8,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Non-enriched" = "grey80",
      "Enriched (FC >= 10)" = "#2591a8",
      "Paip1 & Ybx2" = "#E74C3C"
    )
  ) +
  scale_size_manual(
    values = c(
      "Non-enriched" = 1.5,
      "Enriched (FC >= 10)" = 3.5,
      "Paip1 & Ybx2" = 4.5
    )
  ) +
  scale_alpha_manual(
    values = c(
      "Non-enriched" = 0.4,
      "Enriched (FC >= 10)" = 0.9,
      "Paip1 & Ybx2" = 1
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = c(0.8, 0.2),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(
    title = "Protein Abundance in IP-MS",
    x = expression("Log"[2] * " (IgG Intensity)"),
    y = expression("Log"[2] * " (PAIP1 IP Intensity)")
  )

dev.off()

# ------------------------------------------------------------------------------
# 6. 蛋白富集排序图
# ------------------------------------------------------------------------------

waterfall_df <- plot_df_objective %>%
  filter(True_PAIP1_IP_filled > 1000)

if (nrow(waterfall_df) > 0) {
  pdf(
    "FigS_IPMS_Waterfall_Objective.pdf",
    width = 8,
    height = 6
  )

  ggplot(
    waterfall_df,
    aes(Rank, Log2_FC, fill = Point_Category)
  ) +
    geom_col(width = 1) +
    geom_text_repel(
      aes(label = Label, color = Point_Category),
      size = 4.5,
      fontface = "italic",
      direction = "y",
      nudge_x = max(waterfall_df$Rank) * 0.15,
      box.padding = 0.6,
      max.overlaps = Inf,
      show.legend = FALSE,
      segment.color = "grey50"
    ) +
    scale_fill_manual(
      values = c(
        "Non-enriched" = "grey85",
        "Enriched (FC >= 10)" = "#2591a8",
        "Paip1 & Ybx2" = "#E74C3C"
      )
    ) +
    scale_color_manual(
      values = c(
        "Non-enriched" = "grey",
        "Enriched (FC >= 10)" = "#1e7284",
        "Paip1 & Ybx2" = "#c0392b"
      )
    ) +
    geom_hline(
      yintercept = log2(10),
      linetype = "dashed",
      color = "black"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = c(0.7, 0.8),
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18)
    ) +
    labs(
      title = "Protein Enrichment Ranking",
      x = "Rank",
      y = expression("Log"[2] * " (Enrichment Fold Change)")
    )

  dev.off()
}

# ------------------------------------------------------------------------------
# 7. Top 150 互作蛋白 GO 富集
# ------------------------------------------------------------------------------

ip_data <- read_csv(
  "PAIP1_Clean_Interactome_Final.csv"
)

top_interactors <- ip_data %>%
  arrange(desc(Fold_Change)) %>%
  slice_head(n = 150) %>%
  pull(Genes) %>%
  unique()

entrez_ip <- bitr(
  top_interactors,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

if (nrow(entrez_ip) > 0) {
  ego_ip_top150 <- enrichGO(
    gene = unique(entrez_ip$ENTREZID),
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  top150_go_df <- as.data.frame(ego_ip_top150)

  write.csv(
    top150_go_df,
    "Table_IPMS_Top150_GO_Enrichment.csv",
    row.names = FALSE
  )

  if (nrow(top150_go_df) > 0) {
    pdf(
      "Fig_IPMS_GO_Enrichment.pdf",
      width = 8,
      height = 7
    )

    print(
      dotplot(
        ego_ip_top150,
        split = "ONTOLOGY",
        showCategory = 8
      ) +
        facet_grid(ONTOLOGY ~ ., scales = "free") +
        theme_bw(base_size = 14) +
        labs(title = "GO Enrichment of PAIP1 Interactome")
    )

    dev.off()
  }
}

writeLines(
  top_interactors,
  "Genes_for_STRING_Input.txt"
)

# ------------------------------------------------------------------------------
# 8. 核心互作网络基因 GO 富集
# ------------------------------------------------------------------------------

refined_ip_genes <- c(
  "Paip1", "Ybx2", "Dazl", "Piwil2",
  "Eif4a1", "Eif2s2", "Eef1a1",
  "Hnrnpa1", "Hnrnpc", "Hnrnpa2b1", "Alyref",
  "Thrap3", "Bclaf1", "Erh", "Ncl", "Npm1",
  "Rpl30", "Rpl10a", "Rpl13", "Rpl28", "Rpl12",
  "Rpl23a", "Rpl27", "Rplp0", "Rplp2", "Rpl34",
  "Rpl35", "Rps19", "Rps14", "Rps7", "Rps4x",
  "Rps29", "Rps16", "Rps3", "Rps25", "Rps11",
  "Rps20", "Rps13", "Rps18", "Mrps27", "Ppp1cb",
  "Ppp1r12b"
)

entrez_core <- bitr(
  refined_ip_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

if (nrow(entrez_core) > 0) {
  ego_ip_core <- enrichGO(
    gene = unique(entrez_core$ENTREZID),
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )

  df_core <- as.data.frame(ego_ip_core)

  write.csv(
    df_core,
    "Table_IPMS_Core_GO_Enrichment.csv",
    row.names = FALSE
  )

  target_ip_terms <- data.frame(
    Description = c(
      "RNA splicing",
      "mRNA stabilization",
      "translational initiation",
      "regulation of translation",
      "ribosome biogenesis",
      "cytoplasmic translation",
      "ribonucleoprotein complex biogenesis"
    ),
    Module = c(
      "Upstream mRNA Processing",
      "Upstream mRNA Processing",
      "Translation Initiation",
      "Translation Initiation",
      "Ribosomal Machinery",
      "Ribosomal Machinery",
      "mRNP Granules & Germ Cell"
    )
  )

  plot_ip_df <- df_core %>%
    filter(Description %in% target_ip_terms$Description) %>%
    left_join(target_ip_terms, by = "Description") %>%
    mutate(LogP = -log10(p.adjust))

  if (nrow(plot_ip_df) > 0) {
    biological_order <- c(
      "RNA splicing",
      "mRNA stabilization",
      "translational initiation",
      "regulation of translation",
      "ribosome biogenesis",
      "cytoplasmic translation",
      "ribonucleoprotein complex biogenesis"
    )

    plot_ip_df$Description <- factor(
      plot_ip_df$Description,
      levels = rev(biological_order)
    )

    module_colors <- c(
      "Upstream mRNA Processing" = "#c8b6e2",
      "Translation Initiation" = "#ffd966",
      "Ribosomal Machinery" = "#7bc5f1",
      "mRNP Granules & Germ Cell" = "#2153a3"
    )

    pdf(
      "Fig12_IPMS_Core_Module_Barplot.pdf",
      width = 10,
      height = 6.5
    )

    ggplot(
      plot_ip_df,
      aes(LogP, Description, fill = Module)
    ) +
      geom_col(
        color = "black",
        width = 0.7,
        alpha = 0.9
      ) +
      scale_fill_manual(values = module_colors) +
      geom_vline(
        xintercept = -log10(0.05),
        linetype = "dashed",
        color = "red",
        linewidth = 1
      ) +
      theme_classic(base_size = 16) +
      theme(
        axis.text.y = element_text(
          face = "bold",
          color = "black",
          size = 14
        ),
        axis.text.x = element_text(
          color = "black",
          size = 12
        ),
        axis.title.x = element_text(
          face = "bold",
          size = 15
        ),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 18
        )
      ) +
      labs(
        title = expression(
          "Functional Enrichment of " *
          italic("Paip1") *
          " Core Interactome"
        ),
        x = expression(-Log[10](italic(P)[adj])),
        y = ""
      )

    dev.off()
  }
}

cat("\nCo-IP/MS analysis completed.\n")

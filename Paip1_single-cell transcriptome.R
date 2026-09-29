# Paip1 单细胞转录组分析

# 1. 环境与参数
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggsci)
  library(pheatmap)
  library(stringr)
  library(ggnewscale)
  library(ggforce)
  library(ggrepel)
  library(tidyr)
  library(RColorBrewer)
  library(limma)
  library(pcaMethods)
  library(enrichplot)
  library(Seurat)
  library(harmony)
  library(monocle)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

setwd("D:/LK/")

cols_all <- c(
  "#D6A2D6", "#6BC4B8", "#FFB77A", "#D15F76",
  "#9B8BC6", "#89D3C7", "#3E8E9F", "#F4E76E", "#C97B5F"
)

cols_germ <- c(
  "#9B8BC6", "#B8F0E0", "#5DD4C1", "#2A6A7D",
  "#4A96D7", "#F4E76E", "#C97B5F"
)

cols_monocle <- c("#B8F0E0", "#5DD4C1", "#2A6A7D", "#4A96D7")

# 2. QC、标准化、降维与聚类
cko_data <- Read10X(data.dir = "cKO/")
wt_data <- Read10X(data.dir = "WT/")


CKO <- CreateSeuratObject(cko_data, project = "cko", min.cells = 3, min.features = 200)
WT <- CreateSeuratObject(wt_data, project = "wt", min.cells = 3, min.features = 200)


sc_obj <- merge(WT, y = c(CKO), add.cell.ids = c("WT", "cKO"))
head(colnames(sc_obj))


sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")


VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000)


sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(sc_obj), 10)
plot1 <- VariableFeaturePlot(sc_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("QC_VariableFeatures.png", plot = plot1 + plot2, width = 12, height = 5, dpi = 300)


all.genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all.genes)
sc_obj <- RunPCA(sc_obj, npcs = 30, verbose = FALSE)


ElbowPlot(sc_obj, ndims = 30, reduction = "pca")


sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:15)
sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.4)


sc_obj$orig.ident <- factor(sc_obj$orig.ident, levels = c("wt", "cko"))
p1 <- DimPlot(sc_obj, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sc_obj, reduction = "umap", label = TRUE)
ggsave("Global_UMAP_Group.pdf", plot = plot_grid(p1, p2), width = 12, height = 5, dpi = 300)


saveRDS(sc_obj, file = "sc_obj_step1_clustered.rds")


# 3. Marker 鉴定与细胞类型注释
sc_obj <- readRDS("sc_obj_step1_clustered.rds")


sc_obj <- JoinLayers(sc_obj)
markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "global_markers.csv")


marker_list <- list(
  Sertoli     = c("Ctsl", "Clu", "Cst12", "Cst9", "Defb19", "Ldhb"),
  Leydig      = c("Cyp17a1", "Hsd3b1", "Star", "Hsd3b6", "Fabp3", "Akr1c1"),
  Myoid       = c("Dcn", "Igfbp7", "Gsn", "Myl9", "Mgp", "Rbbp7"),
  Macrophages = c("Lyz2", "C1qb", "C1qc", "C1qa", "Pf4", "Apoe"),
  SPG         = c("Stra8", "Ptma", "Uchl1", "Dazl", "Crabp1", "Nmt2"),
  SPC         = c("Tex12", "Tex101", "Sycp3", "Insl6", "Spag6l", "Tbpl1"),
  Spermatid   = c("Tex29", "Tex36", "Tbc1d23", "Cst13", "Prm2", "Tfam")
)


fp1 <- FeaturePlot(sc_obj, features = c(marker_list$Sertoli, marker_list$Leydig),
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Sertoli_Leydig.png", plot = fp1, width = 13, height = 10, dpi = 300)


fp2 <- FeaturePlot(sc_obj, features = c(marker_list$Myoid, marker_list$Macrophages),
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Myoid_Macrophage.png", plot = fp2, width = 13, height = 10, dpi = 300)


fp3 <- FeaturePlot(sc_obj, features = c(marker_list$SPG, marker_list$SPC),
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Germ1.png", plot = fp3, width = 13, height = 10, dpi = 300)


fp4 <- FeaturePlot(sc_obj, features = marker_list$Spermatid,
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 3)
ggsave("Marker_FeaturePlot_Germ2.png", plot = fp4, width = 9, height = 10, dpi = 300)


lapply(names(marker_list), function(cell_type) {
  p <- VlnPlot(sc_obj, features = marker_list[[cell_type]], pt.size = 0, ncol = 3)
  ggsave(paste0("Marker_VlnPlot_", cell_type, ".png"), plot = p, width = 15, height = 10, dpi = 300)
})


all_markers <- unlist(marker_list)


dot1 <- DotPlot(sc_obj, features = all_markers, group.by = "seurat_clusters") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

ggsave("Marker_DotPlot_All.pdf", plot = dot1, width = 14, height = 6, dpi = 300)


new_cluster_ids <- c(
  '0' = "Leydig cell",
  '1' = "Late SPC",
  '2' = "Round Spermatid",
  '3' = "Late SPC",
  '4' = "Elongating Spermatid",
  '5' = "Early SPC",
  '6' = "Myoid cell",
  '7' = "Round Spermatid",
  '8' = "Sertoli cell",
  '9' = "Round Spermatid",
  '10' = "Early SPC",
  '11' = "Late SPC",
  '12' = "Elongating Spermatid",
  '13' = "Macrophages cell",
  '14' = "Late SPC",
  '15' = "SPG"
)

sc_obj <- RenameIdents(sc_obj, new_cluster_ids)
sc_obj$cell_type <- Idents(sc_obj)


Idents(sc_obj) <- factor(Idents(sc_obj), levels = c("Sertoli cell","Leydig cell","Myoid cell","Macrophages cell",
                                                    "SPG", "Early SPC", "Late SPC","Round Spermatid","Elongating Spermatid"))


saveRDS(sc_obj, file = "sc_obj_step2_annotated.rds")


# 4. 全局细胞类型整理

sc_obj$cell_type_merged <- as.character(sc_obj$cell_type)
sc_obj$cell_type_merged[sc_obj$cell_type_merged %in% c("Round Spermatid", "Elongating Spermatid")] <- "Spermatids"
sc_obj$cell_type_merged <- factor(
  sc_obj$cell_type_merged,
  levels = c(
    "Sertoli cell", "Leydig cell", "Myoid cell", "Macrophages cell",
    "SPG", "Early SPC", "Late SPC", "Spermatids"
  )
)

sc_obj$cell_type_main <- as.character(sc_obj$cell_type_merged)
sc_obj$cell_type_main[sc_obj$cell_type_main %in% c(
  "Sertoli cell", "Leydig cell", "Myoid cell", "Macrophages cell"
)] <- "Somatic cells"
sc_obj$cell_type_main[sc_obj$cell_type_main %in% c("Early SPC", "Late SPC")] <- "SPC"
sc_obj$cell_type_main <- factor(
  sc_obj$cell_type_main,
  levels = c("Somatic cells", "SPG", "SPC", "Spermatids")
)

cols_celltype_merged <- c(
  "Sertoli cell" = "#D6A2D6",
  "Leydig cell" = "#6BC4B8",
  "Myoid cell" = "#FFB77A",
  "Macrophages cell" = "#D15F76",
  "SPG" = "#9B8BC6",
  "Early SPC" = "#89D3C7",
  "Late SPC" = "#3E8E9F",
  "Spermatids" = "#F4E76E"
)

cols_celltype_main <- c(
  "Somatic cells" = "#D3D3D3",
  "SPG" = "#9B8BC6",
  "SPC" = "#89D3C7",
  "Spermatids" = "#F4E76E"
)

Idents(sc_obj) <- sc_obj$cell_type_merged

p_umap_merged <- DimPlot(
  sc_obj,
  group.by = "cell_type_merged",
  split.by = "orig.ident",
  cols = cols_celltype_merged,
  label = TRUE
)

ggsave(
  "UMAP_CellType_Merged.pdf",
  plot = p_umap_merged,
  width = 12,
  height = 5,
  dpi = 300
)

cellnum <- as.data.frame(table(sc_obj$cell_type_merged, sc_obj$orig.ident))
colnames(cellnum) <- c("Celltype", "Group", "Value")

group_totals <- cellnum %>%
  group_by(Group) %>%
  summarise(Total = sum(Value), .groups = "drop")

cellnum <- cellnum %>%
  left_join(group_totals, by = "Group") %>%
  mutate(Proportion = Value / Total * 100)

p_celltype_proportion <- ggplot(
  cellnum,
  aes(x = Group, y = Proportion, fill = Celltype)
) +
  geom_bar(stat = "identity", width = 0.4, color = "white") +
  scale_fill_manual(values = cols_celltype_merged) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  labs(x = NULL, y = "Cell proportion (%)") +
  theme_classic()

ggsave(
  "CellType_Proportion.pdf",
  plot = p_celltype_proportion,
  width = 7,
  height = 6,
  dpi = 300
)

germ_types <- c("SPG", "Early SPC", "Late SPC", "Spermatids")

germ_count <- cellnum %>%
  filter(Celltype %in% germ_types)

p_germ_count <- ggplot(
  germ_count,
  aes(x = Celltype, y = Value, fill = Group)
) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  labs(x = "Cell type", y = "Number of cells", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("wt" = "#f933ff", "cko" = "#2591a8"))

ggsave(
  "GermCell_Count.pdf",
  plot = p_germ_count,
  width = 7,
  height = 5,
  dpi = 300
)

somatic_cells <- c(
  "Sertoli cell", "Leydig cell", "Myoid cell", "Macrophages cell"
)

sertoli_counts <- cellnum %>%
  filter(Celltype == "Sertoli cell") %>%
  select(Group, Sertoli_Value = Value)

somatic_data <- cellnum %>%
  filter(Celltype %in% somatic_cells) %>%
  left_join(sertoli_counts, by = "Group") %>%
  mutate(Relative_Abundance = Value / Sertoli_Value)

somatic_data$Celltype <- factor(somatic_data$Celltype, levels = somatic_cells)

p_somatic <- ggplot(
  somatic_data,
  aes(x = Celltype, y = Relative_Abundance, fill = Group)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("wt" = "#f933ff", "cko" = "#2591a8")) +
  labs(
    x = NULL,
    y = "Relative cell abundance\n(Normalized to Sertoli cells)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85, 0.85)
  )

ggsave(
  "Somatic_Cell_Relative_Abundance.pdf",
  plot = p_somatic,
  width = 6,
  height = 5.5,
  dpi = 300
)

saveRDS(sc_obj, file = "sc_obj_celltype_merged.rds")

Idents(sc_obj) <- sc_obj$cell_type

# 5. 生殖细胞再聚类与注释
sc_obj <- readRDS("sc_obj_step2_annotated.rds")


germ_cells <- subset(sc_obj, idents = c("SPG", "Early SPC", "Late SPC", "Round Spermatid", "Elongating Spermatid"))


germ_sce <- CreateSeuratObject(counts = GetAssayData(germ_cells, assay = "RNA", layer = 'counts'),
                               meta.data = germ_cells@meta.data)

germ_sce <- NormalizeData(germ_sce) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)


germ_sce <- RunHarmony(germ_sce, group.by.var = "orig.ident")


germ_sce <- FindNeighbors(germ_sce, reduction = "harmony", dims = 1:20)
germ_sce <- FindClusters(germ_sce, resolution = 0.5)
germ_sce <- RunUMAP(germ_sce, reduction = "harmony", dims = 1:20)


current_ids <- 0:15

new_germ_ids <- c(
  '3' = "Sertoli cell", '9' = "Sertoli cell",
  '13' = "SPG",
  '14' = "Leptotene",
  '4' = "Zygotene",
  '5' = "Pachytene", '12' = "Pachytene",
  '1' = "Diplotene", '7' = "Diplotene", '6' = "Diplotene",
  '2' = "Round Std", '8' = "Round Std", '10' = "Round Std",
  '0' = "Elong Std", '11' = "Elong Std", '15' = "Elong Std"
)

all_clusters <- levels(germ_sce)
names(new_germ_ids) <- as.character(names(new_germ_ids))


germ_sce <- RenameIdents(germ_sce, new_germ_ids)
germ_sce$celltype <- Idents(germ_sce)


germ_sce$celltype <- factor(germ_sce$celltype, levels = c("Sertoli cell","SPG","Leptotene","Zygotene","Pachytene","Diplotene","Round Std","Elong Std"))
Idents(germ_sce) <- germ_sce$celltype


save(germ_sce, file = "germ_sce_step3_annotated.Rdata")


meta_df <- germ_sce@meta.data[, c("celltype", "orig.ident")]


cell_stats <- meta_df %>%
  group_by(orig.ident, celltype) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(Total = sum(Count)) %>%
  mutate(Proportion = Count / Total * 100)


target_cells <- c("Leptotene", "Zygotene", "Pachytene", "Diplotene")
plot_data <- cell_stats %>%
  filter(celltype %in% target_cells) %>%
  mutate(celltype = factor(celltype, levels = target_cells))


cols_meiosis <- c("#B8F0E0", "#5DD4C1", "#2A6A7D", "#4A96D7")
names(cols_meiosis) <- target_cells


p_bar_all <- ggplot(plot_data, aes(x = orig.ident, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.5, color = "white") +
  scale_fill_manual(values = cols_meiosis) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0, 0), limits = c(0, 105)) +
  scale_x_discrete(labels = toupper) +
  labs(x = NULL, y = "Percentage of Cell Type (%)", title = "Cell Proportion") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  )

ggsave("Proportion_Barplot_All.pdf", plot = p_bar_all, width = 6, height = 5)


p_bar_facet <- ggplot(plot_data, aes(x = Proportion, y = orig.ident, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(celltype ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = cols_meiosis, guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)),
            hjust = -0.2, size = 3.5) +
  labs(x = "Proportion (%)", y = NULL) +
  theme_classic() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)))

ggsave("Proportion_Barplot_Faceted.pdf", plot = p_bar_facet, width = 7, height = 6)


sunburst_data <- plot_data %>%
  mutate(orig.ident = as.factor(orig.ident)) %>%
  arrange(orig.ident, celltype) %>%
  group_by(orig.ident) %>%
  mutate(

    ymax = cumsum(Proportion),
    ymin = lag(ymax, default = 0),

    label_pos = (ymin + ymax) / 2
  )
sunburst_data$r0 <- as.numeric(sunburst_data$orig.ident) - 0.4
sunburst_data$r  <- as.numeric(sunburst_data$orig.ident) + 0.4


p_sunburst <- ggplot(sunburst_data) +
  geom_arc_bar(aes(x0 = 0, y0 = 0,
                   r0 = r0, r = r,
                   start = ymin * 2 * pi / 100,
                   end = ymax * 2 * pi / 100,
                   fill = celltype),
               color = "white", size = 0.5) +

  geom_text(data = subset(sunburst_data, celltype == "Pachytene" | Proportion > 5),
            aes(x = (r0 + 0.45) * sin(label_pos * 2 * pi / 100),
                y = (r0 + 0.45) * cos(label_pos * 2 * pi / 100),
                label = sprintf("%.1f%%", Proportion)),
            size = 3, fontface = "bold") +

  annotate("text", x = 0, y = 0, label = "Group", fontface = "bold") +
  scale_fill_manual(values = cols_meiosis) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "right", legend.title = element_blank())

ggsave("Proportion_Sunburst.pdf", plot = p_sunburst, width = 8, height = 6)


# 6. PAIP1 表达模式

sc_obj$orig.ident <- factor(sc_obj$orig.ident, levels = c("wt", "cko"))
germ_sce$orig.ident <- factor(germ_sce$orig.ident, levels = c("wt", "cko"))

wt_global <- subset(sc_obj, orig.ident == "wt")
wt_germ <- subset(germ_sce, orig.ident == "wt")

wt_global_8 <- wt_global
Idents(wt_global_8) <- wt_global_8$cell_type_merged

wt_global_4 <- wt_global
Idents(wt_global_4) <- wt_global_4$cell_type_main

wt_germ$celltype_paip1 <- as.character(wt_germ$celltype)
wt_germ$celltype_paip1[wt_germ$celltype_paip1 %in% c("Round Std", "Elong Std")] <- "Spermatids"
wt_germ$celltype_paip1 <- factor(
  wt_germ$celltype_paip1,
  levels = c(
    "Sertoli cell", "SPG", "Leptotene", "Zygotene",
    "Pachytene", "Diplotene", "Spermatids"
  )
)
Idents(wt_germ) <- wt_germ$celltype_paip1

p_dot_global <- DotPlot(
  wt_global_4,
  features = "Paip1"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#f933ff") +
  labs(
    title = "Paip1 Expression Pattern in WT",
    x = NULL,
    y = NULL
  )

p_dot_germ <- DotPlot(
  wt_germ,
  features = "Paip1"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "#f933ff") +
  labs(
    title = "Paip1 Expression Pattern in WT Germ Cells",
    x = NULL,
    y = NULL
  )

pdf("Fig_scRNA_Paip1_WT_DotPlots.pdf", width = 12, height = 5)
cowplot::plot_grid(p_dot_global, p_dot_germ, ncol = 2)
dev.off()

p_feat_germ <- FeaturePlot(
  wt_germ,
  features = "Paip1",
  pt.size = 0.5,
  order = TRUE,
  cols = c("lightgrey", "#f933ff")
) +
  ggtitle("Paip1 Spatial Distribution (WT Germ Cells)")

pdf("Fig_scRNA_Paip1_WT_FeaturePlot.pdf", width = 6, height = 5)
print(p_feat_germ)
dev.off()

paip1_germ_colors <- c(
  "Sertoli cell" = "#D3D3D3",
  "SPG" = "#9B8BC6",
  "Leptotene" = "#B8F0E0",
  "Zygotene" = "#5DD4C1",
  "Pachytene" = "#2A6A7D",
  "Diplotene" = "#4A96D7",
  "Spermatids" = "#F4E76E"
)

p_vln_germ <- VlnPlot(
  wt_germ,
  features = "Paip1",
  pt.size = 0.1,
  cols = paip1_germ_colors
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Endogenous Paip1 Expression (WT)",
    y = "Expression Level"
  )

pdf(
  "Fig_scRNA_Paip1_WT_VlnPlot_Spermatids_Merged.pdf",
  width = 8,
  height = 5
)
print(p_vln_germ)
dev.off()

global_colors_merged <- c(
  "Sertoli cell" = "#D6A2D6",
  "Leydig cell" = "#6BC4B8",
  "Myoid cell" = "#FFB77A",
  "Macrophages cell" = "#D15F76",
  "SPG" = "#9B8BC6",
  "Early SPC" = "#89D3C7",
  "Late SPC" = "#3E8E9F",
  "Spermatids" = "#F4E76E"
)

p_vln_global <- VlnPlot(
  wt_global_8,
  features = "Paip1",
  pt.size = 0,
  cols = global_colors_merged
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Endogenous Paip1 Expression (WT Global Cells)",
    y = "Expression Level"
  )

pdf(
  "Fig_scRNA_Paip1_WT_VlnPlot_Global_Merged.pdf",
  width = 9,
  height = 5
)
print(p_vln_global)
dev.off()

p_vln_global_main <- VlnPlot(
  wt_global_4,
  features = "Paip1",
  pt.size = 0.1,
  cols = cols_celltype_main
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Endogenous Paip1 Expression (WT Main Lineages)",
    y = "Expression Level"
  )

pdf(
  "Fig_scRNA_Paip1_WT_VlnPlot_4Groups.pdf",
  width = 6.5,
  height = 5
)
print(p_vln_global_main)
dev.off()

target_meiotic_cells <- c(
  "Leptotene", "Zygotene", "Pachytene", "Diplotene"
)

wt_meiotic <- subset(
  wt_germ,
  celltype %in% target_meiotic_cells
)
wt_meiotic$celltype <- factor(
  wt_meiotic$celltype,
  levels = target_meiotic_cells
)
Idents(wt_meiotic) <- wt_meiotic$celltype

meiotic_colors <- c(
  "Leptotene" = "#B8F0E0",
  "Zygotene" = "#5DD4C1",
  "Pachytene" = "#2A6A7D",
  "Diplotene" = "#4A96D7"
)

p_vln_meiotic <- VlnPlot(
  wt_meiotic,
  features = "Paip1",
  pt.size = 0.1,
  cols = meiotic_colors
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Endogenous Paip1 Expression (Meiotic Prophase I)",
    y = "Expression Level"
  )

pdf(
  "Fig_scRNA_Paip1_WT_VlnPlot_Meiotic_Prophase.pdf",
  width = 6,
  height = 5
)
print(p_vln_meiotic)
dev.off()

Idents(sc_obj) <- sc_obj$cell_type
Idents(germ_sce) <- germ_sce$celltype

# 7. Monocle 2 轨迹分析
load("germ_sce_step3_annotated.Rdata")


target_cells <- c("Leptotene", "Zygotene", "Pachytene", "Diplotene")
seurat_sub <- subset(germ_sce, idents = target_cells)


data <- GetAssayData(seurat_sub, assay = "RNA", layer = 'counts')
pd <- new("AnnotatedDataFrame", data = seurat_sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~celltype", cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)


cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"Leptotene"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return(1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))


p1 <- plot_cell_trajectory(cds, color_by = "celltype") + scale_color_manual(values = cols_monocle)
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave("Monocle_Trajectory_Combined.pdf", plot = p1 + p2, width = 16, height = 8)


save(cds, file = "monocle_cds_final.RData")


# 8. 差异表达与 GO 富集
Meio <- subset(germ_sce, idents = "Pachytene")
Idents(Meio) <- "orig.ident"

DE_results <- FindMarkers(
  Meio,
  ident.1 = "cko",
  ident.2 = "wt",
  logfc.threshold = 0.5,
  min.pct = 0.1,
  test.use = "wilcox"
)


DE_results_fil <- DE_results %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

write.csv(DE_results_fil, "DE_results_Pachytene.csv")


gene_list <- rownames(DE_results_fil)
eg <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
eg <- na.omit(eg)


go_res <- enrichGO(
  gene = eg$ENTREZID,
  OrgDb = "org.Mm.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


if (!exists("go_res")) stop("Error: 'go_res' object not found. Please run Step 5.2 first.")
go_df <- as.data.frame(go_res)


process_go_data <- function(df, ontology_type, top_n = 10) {
  subset_df <- df[df$ONTOLOGY == ontology_type, ]
  if (nrow(subset_df) == 0) return(NULL)


  subset_df <- head(subset_df, top_n)


  subset_df$Description_short <- str_trunc(subset_df$Description, width = 50, side = "right")


  subset_df$Description_short <- factor(subset_df$Description_short, levels = rev(subset_df$Description_short))

  return(subset_df)
}

data_BP <- process_go_data(go_df, "BP")
data_CC <- process_go_data(go_df, "CC")
data_MF <- process_go_data(go_df, "MF")


mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11, color = "black"),
  plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  panel.grid.minor = element_blank()
)


plot_GO_bar <- function(data, title_prefix, color_palette = "Blues") {
  if (is.null(data)) {
    warning(paste("No data for", title_prefix))
    return(NULL)
  }

  ggplot(data, aes(x = Count, y = Description_short, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +

    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +

    scale_fill_distiller(palette = color_palette, direction = 1, name = "-log10(padj)") +
    labs(x = "Gene Number", y = NULL,
         title = paste0(title_prefix, " Enrichment")) +
    theme_bw() +
    mytheme
}


p_BP <- plot_GO_bar(data_BP, "Biological Process", "Oranges")
if (!is.null(p_BP)) {
  print(p_BP)
  ggsave("GO_Enrichment_BP.pdf", plot = p_BP, width = 9, height = 7)
}


p_CC <- plot_GO_bar(data_CC, "Cellular Component", "Reds")
if (!is.null(p_CC)) {
  print(p_CC)
  ggsave("GO_Enrichment_CC.pdf", plot = p_CC, width = 9, height = 7)
}


p_MF <- plot_GO_bar(data_MF, "Molecular Function", "Blues")
if (!is.null(p_MF)) {
  print(p_MF)
  ggsave("GO_Enrichment_MF.pdf", plot = p_MF, width = 9, height = 7)
}


if (!is.null(data_BP)) {
  p_custom <- ggplot(data_BP, aes(x = Count, y = Description_short, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +

    scale_fill_gradient(low = "#077878", high = "#f933ff", name = "-log10(padj)") +
    labs(x = "Gene Number", y = NULL, title = "BP Enrichment (Custom Color)") +
    theme_bw() +
    mytheme

  ggsave("GO_Enrichment_BP_CustomColor.pdf", plot = p_custom, width = 9, height = 7)
}


# 9. 跨物种分析
create_seurat_expression_matrix <- function(seurat_objs, genes) {
  expression_list <- list()

  for(species in names(seurat_objs)) {
    obj <- seurat_objs[[species]]


    genes_present <- intersect(genes, rownames(obj))
    missing_genes <- setdiff(genes, rownames(obj))
    if(length(missing_genes) > 0) {
      warning(paste(species, "missing genes:", paste(missing_genes, collapse=", ")))
    }


    data_matrix <- LayerData(obj, assay = "RNA", layer = "data")


    avg_expr <- rowMeans(data_matrix[genes_present, , drop=FALSE], na.rm = TRUE)


    full_expr <- rep(0, length(genes))
    names(full_expr) <- genes
    full_expr[names(avg_expr)] <- avg_expr

    expression_list[[species]] <- full_expr
  }

  expression_matrix <- do.call(cbind, expression_list)
  rownames(expression_matrix) <- genes
  colnames(expression_matrix) <- names(seurat_objs)

  return(expression_matrix)
}


clean_seurat_factors <- function(obj) {
  obj@meta.data$group <- droplevels(as.factor(obj@meta.data$group))
  obj@meta.data$celltype <- droplevels(as.factor(obj@meta.data$celltype))
  return(obj)
}


to_mouse_symbol <- function(genes) {
  return(str_to_title(genes))
}


calculate_expression_diff <- function(obj_spc, obj_std, gene_list, label) {

  avg_spc <- AverageExpression(obj_spc, features = gene_list, verbose = FALSE)$RNA
  avg_spc <- rowMeans(avg_spc)


  avg_std <- AverageExpression(obj_std, features = gene_list, verbose = FALSE)$RNA
  avg_std <- rowMeans(avg_std)


  res <- data.frame(
    Gene = names(avg_spc),
    Mean_SPC = as.numeric(avg_spc),
    Mean_Std = as.numeric(avg_std[names(avg_spc)])
  )


  res$Difference <- res$Mean_SPC - res$Mean_Std
  res$Species <- label


  res <- res %>% arrange(desc(abs(Difference)))

  return(res)
}


setwd("D:/sheep/")
load("sce2_dim20_reso0.3_zhushi.Rdata")
Tsheep <- sce
Idents(Tsheep) <- "group"
rm(sce)


clean_genes_sheep <- gsub("^gene-", "", rownames(Tsheep))
rownames(Tsheep) <- clean_genes_sheep
rownames(Tsheep@assays$RNA) <- clean_genes_sheep


setwd("D:/Y/mouse/")
load("Tmouse.Rdata")
Tmouse <- data3
Tmouse$group <- Tmouse$orig.ident
Idents(Tmouse) <- "group"
Tmouse[["RNA"]] <- JoinLayers(Tmouse[["RNA"]])
rm(data3)


clean_genes_mouse <- toupper(rownames(Tmouse))
rownames(Tmouse) <- clean_genes_mouse
rownames(Tmouse@assays$RNA) <- clean_genes_mouse


setwd("D:/Y/human/")
load("Thuman.Rdata")
Thuman <- data3
Idents(Thuman) <- "group"
Thuman[["RNA"]] <- JoinLayers(Thuman[["RNA"]])
rm(data3)


setwd("D:/LK/Cattle/")
load("Tcattle.Rdata")
Tcattle <- immune.combined
Tcattle[["RNA"]] <- JoinLayers(Tcattle[["RNA"]])
rm(immune.combined)


setwd("D:/LK/Cross/")


seurat_objects <- list(
  Human = Thuman,
  Mouse = Tmouse,
  Sheep = Tsheep
)


new_gene <- Three_species_share_genes$gene
selected_genes <- new_gene


expression_matrix <- create_seurat_expression_matrix(seurat_objects, selected_genes)

cat("Expression matrix generated.\n")


expression_matrix_log <- log2(expression_matrix + 1)


custom_colors <- colorRampPalette(c("#313695", "#abd9e9",
                                    "#fee090", "#fdae61", "#f46d43", "#d73027"))(100)

real_data_heatmap_new <- pheatmap(expression_matrix_log,
                                  scale = "row",
                                  cluster_rows = TRUE,
                                  cluster_cols = FALSE,
                                  border_color = "grey60",
                                  color = custom_colors,
                                  show_rownames = TRUE,
                                  show_colnames = TRUE,
                                  fontsize_row = 10,
                                  fontsize_col = 12,
                                  angle_col = 45,
                                  cellwidth = 30,
                                  cellheight = 20,
                                  annotation_legend = FALSE,
                                  main = "New Gene Set Expression Patterns Across Species\n(Row Z-score Scaled)",
                                  silent = FALSE)

ggsave("new_gene_heatmap.pdf", plot = real_data_heatmap_new, width = 9, height = 12, dpi = 300)


setwd("D:/sheep/")
load("g_sce_subset_dim20_reso0.6_4_zhushi.Rdata")
sheep <- g_sce_subset

clean_genes <- gsub("^gene-", "", rownames(sheep))
rownames(sheep) <- clean_genes
rownames(sheep@assays$RNA) <- clean_genes
rm(g_sce_subset)


setwd("D:/Y/mouse/")
load("g_sce_dim20_reso1.2_4_zhushi.Rdata")
mouse <- g_sce
Idents(mouse) <- "celltype"
mouse@meta.data$group <- mouse@meta.data$orig.ident
rm(g_sce)


setwd("D:/Y/human/")
load("g_sce_dim20_reso1.2_4_zhushi.Rdata")
human <- g_sce
Idents(human) <- "celltype"
rm(g_sce)


setwd("D:/LK/Cross/")


Idents(human) <- "group"
human_DE_results <- FindMarkers(
  object = human,
  ident.1 = "11yo",
  ident.2 = "7yo",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  slot = "data",
  min.diff.pct = 0.05,
  only.pos = FALSE
)
write.csv(as.data.frame(human_DE_results), file = "human_DE_results.csv")


Idents(mouse) <- "orig.ident"
mouse_DE_results <- FindMarkers(
  object = mouse,
  ident.1 = "PND8",
  ident.2 = "PND5",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  slot = "data",
  min.diff.pct = 0.05,
  only.pos = FALSE
)
write.csv(as.data.frame(mouse_DE_results), file = "mouse_DE_results.csv")


Idents(sheep) <- "group"
sheep_DE_results <- FindMarkers(
  object = sheep,
  ident.1 = "3m",
  ident.2 = "1m",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  slot = "data",
  min.diff.pct = 0.05,
  only.pos = FALSE
)
write.csv(as.data.frame(sheep_DE_results), file = "sheep_DE_results.csv")


TF <- Post_tran_TF$gene


target_spc_types1 <- c("Preleptotene", "Leptotene/Zygotene", "Pachytene", "Diplotene")
target_spc_types2 <- c("Round Std", "Elong Std")


target_group_human1 <- c("11yo", "13yo")
target_group_human2 <- c("14yo", "25yo")
human_SPC <- subset(human, subset = group %in% target_group_human1 & celltype %in% target_spc_types1)
human_Std <- subset(human, subset = group %in% target_group_human2 & celltype %in% target_spc_types2)


target_group_mouse1 <- c("PND8", "PND10", "PND14")
target_group_mouse2 <- c("PND30", "PND35")
mouse_SPC <- subset(mouse, subset = group %in% target_group_mouse1 & celltype %in% target_spc_types1)
mouse_Std <- subset(mouse, subset = group %in% target_group_mouse2 & celltype %in% target_spc_types2)


target_group_sheep1 <- c("3m", "5m")
target_group_sheep2 <- c("6m", "12m")
sheep_SPC <- subset(sheep, subset = group %in% target_group_sheep1 & celltype %in% target_spc_types1)
sheep_Std <- subset(sheep, subset = group %in% target_group_sheep2 & celltype %in% target_spc_types2)


human_SPC <- clean_seurat_factors(human_SPC)
human_Std <- clean_seurat_factors(human_Std)
mouse_SPC <- clean_seurat_factors(mouse_SPC)
mouse_Std <- clean_seurat_factors(mouse_Std)
sheep_SPC <- clean_seurat_factors(sheep_SPC)
sheep_Std <- clean_seurat_factors(sheep_Std)


intersect_human <- intersect(TF, rownames(human_SPC))
intersect_human <- intersect(intersect_human, rownames(human_Std))


intersect_sheep <- intersect(TF, rownames(sheep_SPC))
intersect_sheep <- intersect(intersect_sheep, rownames(sheep_Std))


TF_mouse_format <- to_mouse_symbol(TF)
intersect_mouse <- intersect(TF_mouse_format, rownames(mouse_SPC))
intersect_mouse <- intersect(intersect_mouse, rownames(mouse_Std))

cat("Intersections found: Human:", length(intersect_human), "| Sheep:", length(intersect_sheep), "| Mouse:", length(intersect_mouse), "\n")


diff_human <- calculate_expression_diff(human_SPC, human_Std, intersect_human, "Human")
diff_mouse <- calculate_expression_diff(mouse_SPC, mouse_Std, intersect_mouse, "Mouse")
diff_sheep <- calculate_expression_diff(sheep_SPC, sheep_Std, intersect_sheep, "Sheep")


write.csv(diff_human, "TF_Diff_Human.csv")
write.csv(diff_mouse, "TF_Diff_Mouse.csv")
write.csv(diff_sheep, "TF_Diff_Sheep.csv")


diff_human$Gene_Std_Name <- toupper(diff_human$Gene)
diff_mouse$Gene_Std_Name <- toupper(diff_mouse$Gene)
diff_sheep$Gene_Std_Name <- toupper(diff_sheep$Gene)

diff_combined <- bind_rows(diff_human, diff_mouse, diff_sheep)
diff_combined$Species <- factor(diff_combined$Species, levels = c("Human", "Mouse", "Sheep"))
diff_combined$Direction <- ifelse(diff_combined$Difference > 0, "Higher in SPC", "Higher in Std")
diff_combined$Direction <- factor(diff_combined$Direction, levels = c("Higher in Std", "Higher in SPC"))


gene_order <- diff_combined %>%
  group_by(Gene_Std_Name) %>%
  summarise(Total_Abs_Diff = sum(abs(Difference))) %>%
  arrange(Total_Abs_Diff) %>%
  pull(Gene_Std_Name)
diff_combined$Gene_Std_Name <- factor(diff_combined$Gene_Std_Name, levels = gene_order)

p1 <- ggplot(diff_combined, aes(x = Gene_Std_Name, y = Difference, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Species, scales = "free_x") +
  coord_flip() +
  scale_fill_manual(values = c("Higher in Std" = "#4575b4", "Higher in SPC" = "#d73027")) +
  theme_bw() +
  labs(title = "Expression Intensity Difference (SPC - Std) of TF Genes",
       y = "Difference (Mean Expression: SPC - Std)", x = "Gene Symbol (Standardized)") +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 8))

ggsave("Diff_TF.pdf", plot = p1, width = 12, height = 5, dpi = 300)


plot_data <- diff_combined %>%
  dplyr::select(Gene_Std_Name, Species, Mean_SPC, Mean_Std, Direction) %>%
  pivot_longer(cols = c("Mean_SPC", "Mean_Std"), names_to = "Condition", values_to = "Expression")

plot_data$Condition <- gsub("Mean_", "", plot_data$Condition)
plot_data$Condition <- factor(plot_data$Condition, levels = c("SPC", "Std"))

p3 <- ggplot(plot_data, aes(x = Condition, y = Expression, group = Gene_Std_Name)) +
  geom_line(aes(color = Direction), size = 0.8, alpha = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~Species, scales = "free_y") +
  scale_color_manual(values = c("Higher in Std" = "#4575b4", "Higher in SPC" = "#d73027")) +
  theme_bw() +
  labs(title = "Gene Expression Intensity Changes: SPC vs Std",
       y = "Mean Expression Intensity", x = "Condition", color = "Trend") +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  geom_text_repel(aes(label = Gene_Std_Name), size = 3, direction = "y",
                  min.segment.length = 0, max.overlaps = 15)

ggsave("SlopeChart_General.pdf", plot = p3, width = 16, height = 7, dpi = 300)


problem_genes <- c("DDX5", "DAZL", "PABPC1", "PABPC1", "PAIP2", "RBMXL2", "PABPC2", "PAIP1", "MOV10L1", "DDX4", "YBX2")

p3_high <- ggplot(plot_data, aes(x = Condition, y = Expression, group = Gene_Std_Name)) +
  geom_line(aes(color = Direction), size = 0.8, alpha = 0.1) +
  geom_line(data = plot_data %>% filter(Gene_Std_Name %in% problem_genes),
            aes(color = Direction), size = 1.2, alpha = 1) +
  geom_point(size = 2, alpha = 0.2) +
  geom_point(data = plot_data %>% filter(Gene_Std_Name %in% problem_genes),
             size = 3, color = "black") +
  geom_text_repel(data = plot_data %>% filter(Gene_Std_Name %in% problem_genes),
                  aes(label = Gene_Std_Name),
                  size = 4, box.padding = 0.5, point.padding = 0.5, min.segment.length = 0) +
  facet_wrap(~Species, scales = "free_y") +
  scale_color_manual(values = c("Higher in Std" = "#4575b4", "Higher in SPC" = "#d73027")) +
  theme_bw() +
  labs(title = "Gene Expression Intensity Changes: SPC vs Std (Highlighted)",
       y = "Mean Expression Intensity", x = "Condition", color = "Trend") +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  theme(strip.text = element_text(face = "bold", size = 12))

ggsave("SlopeChart_Highlighted.pdf", plot = p3_high, width = 16, height = 7, dpi = 300)


process_and_plot_go <- function(file_path, suffix) {


  tryCatch({
    df <- readxl::read_excel(file_path)
  }, error = function(e) {
    stop(paste("Cannot read file:", file_path))
  })


  df_plot <- df %>%
    mutate(p.adjust = as.numeric(`p.adjust`)) %>%
    mutate(logP = -log10(p.adjust)) %>%
    mutate(Ontology = gsub("GO:", "", source))


  df_separated <- df_plot %>%
    mutate(
      Main_Ontology = case_when(
        Ontology == "MF" ~ "Molecular Function (MF)",
        Ontology %in% c("BP") ~ "Biological Process (BP)",
        TRUE ~ Ontology
      )
    ) %>%
    arrange(Main_Ontology, logP) %>%
    group_by(Main_Ontology) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    ungroup()


  p_dotplot_final <- ggplot(df_separated, aes(x = logP, y = Description)) +
    geom_point(aes(color = logP), size = 4) +
    facet_wrap(~ Main_Ontology, scales = "free_y", ncol = 2) +
    scale_color_gradient(low = "#fee0d2", high = "#de2d26", name = expression("-log"[10]*"(P.adj)")) +
    labs(x = expression("-log"[10]*"(Adjusted P-value)"), y = NULL,
         title = paste("GO Enrichment Dot Plot -", suffix)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold")
    )

  ggsave(paste0("GO_", suffix, "_Combined.pdf"), plot = p_dotplot_final, width = 12, height = 7, dpi = 300)


  ontologies_to_plot <- c("Molecular Function (MF)", "Biological Process (BP)")

  for (current_ontology in ontologies_to_plot) {
    df_subset <- df_separated %>%
      filter(Main_Ontology == current_ontology) %>%
      group_by(Main_Ontology) %>%
      mutate(Description = factor(Description, levels = rev(Description))) %>%
      ungroup()

    if (nrow(df_subset) == 0) {
      message(paste("No data for", current_ontology, "in", suffix))
      next
    }

    p_dotplot <- ggplot(df_subset, aes(x = logP, y = Description)) +
      geom_point(aes(color = logP), size = 4) +
      scale_color_gradient(low = "#fee0d2", high = "#de2d26", name = expression("-log"[10]*"(P.adj)")) +
      labs(x = expression("-log"[10]*"(Adjusted P-value)"), y = NULL,
           title = paste(current_ontology, "Enrichment -", suffix)) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.position = "bottom"
      )

    file_name <- gsub("\\s+|\\(|\\)", "_", current_ontology)
    ggsave(paste0(file_name, "_", suffix, ".pdf"), plot = p_dotplot, width = 10, height = 6, dpi = 300)
    message(paste("Saved separate plot for", current_ontology, "in", suffix))
  }
}


setwd("D:/LK/Cross/")


if(file.exists("term_human.xlsx")) {
  process_and_plot_go("term_human.xlsx", "human")
}


if(file.exists("term_mouse.xlsx")) {
  process_and_plot_go("term_mouse.xlsx", "mouse")
}


if(file.exists("term_sheep.xlsx")) {
  process_and_plot_go("term_sheep.xlsx", "sheep")
}


# 10. 蛋白质组分析
setwd("D:/蛋白质组/danbaizhifenxi/")


my_colors <- c("Control" = "#0f99b2", "cko" = "#f74ed6", "CKO" = "#f74ed6")
volcano_colors <- c("Up" = "#f933ff", "Down" = "#077878", "NS" = "grey")


raw_data <- read_excel("data1.xlsx")
protein_ids <- raw_data[[1]]
expr_data <- as.matrix(raw_data[-1])
rownames(expr_data) <- protein_ids


expr_data[is.na(expr_data)] <- min(expr_data, na.rm = TRUE) * 0.5


norm_data <- normalizeBetweenArrays(expr_data, method = "quantile")


gene_symbols <- mapIds(org.Mm.eg.db, keys = rownames(norm_data), column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")

rownames(norm_data) <- ifelse(is.na(gene_symbols), rownames(norm_data), gene_symbols)

head(norm_data)


group_list <- factor(c(rep("Control", 3), rep("cko", 3)), levels = c("Control", "cko"))


pca_res <- pca(t(norm_data), method = "svd", nPcs = 2, scale = "uv")
pca_scores <- as.data.frame(scores(pca_res))
pca_scores$group <- group_list
pca_scores$sample <- colnames(norm_data)


p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 3) +
  ggforce::geom_mark_ellipse(alpha = 0.2, show.legend = FALSE) +
  geom_text_repel(aes(label = sample), show.legend = FALSE) +
  labs(title = "PCA Analysis",
       x = paste0("PC1 (", round(pca_res@R2[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(pca_res@R2[2] * 100, 1), "%)")) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme_bw()

print(p_pca)
ggsave("PCA_plot.pdf", p_pca, width = 6, height = 5)


design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)


fit <- lmFit(norm_data, design)


contrast.matrix <- makeContrasts(cko - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


deg_results <- topTable(fit2, adjust = "BH", sort.by = "P", number = Inf)
deg_results$Protein <- rownames(deg_results)


logFC_cutoff <- 1
pval_cutoff <- 0.05

deg_results$Significance <- case_when(
  deg_results$adj.P.Val < pval_cutoff & deg_results$logFC > logFC_cutoff ~ "Up",
  deg_results$adj.P.Val < pval_cutoff & deg_results$logFC < -logFC_cutoff ~ "Down",
  TRUE ~ "NS"
)


write.csv(deg_results, "DE_results_Limma.csv")


up_n <- sum(deg_results$Significance == "Up")
down_n <- sum(deg_results$Significance == "Down")


deg_results$Label <- NA
top_genes <- head(order(deg_results$adj.P.Val), 20)
deg_results$Label[top_genes] <- deg_results$Protein[top_genes]

p_vol <- ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey30") +
  geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20) +
  annotate("text", x = max(deg_results$logFC), y = max(-log10(deg_results$adj.P.Val)),
           label = paste("Up:", up_n), color = volcano_colors["Up"], hjust = 1) +
  annotate("text", x = min(deg_results$logFC), y = max(-log10(deg_results$adj.P.Val)),
           label = paste("Down:", down_n), color = volcano_colors["Down"], hjust = 0) +
  labs(title = "Volcano Plot", x = expression(log[2]("Fold Change")), y = expression(-log[10]("Adj P-value"))) +
  theme_bw()

print(p_vol)
ggsave("Volcano_plot.pdf", p_vol, width = 7, height = 6)


if(exists("YBX2_PAIP1_bindingRNA")) {
  gene_list <- YBX2_PAIP1_bindingRNA$Symbol


  target_expr <- norm_data[rownames(norm_data) %in% gene_list, , drop=FALSE]

  if(nrow(target_expr) > 0) {

    df_long <- as.data.frame(target_expr) %>%
      mutate(Protein = rownames(.)) %>%
      pivot_longer(-Protein, names_to = "Sample", values_to = "Expression") %>%
      mutate(Group = ifelse(grepl("Control", Sample), "Control", "cko"))


    p_vln <- ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white") +
      labs(title = "Expression of Target Gene Set") +
      scale_fill_manual(values = my_colors) +
      theme_bw()

    ggsave("TargetSet_Violin.pdf", p_vln, width = 6, height = 5)
  }
}


sig_genes <- subset(deg_results, Significance != "NS")$Protein


gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

if(nrow(gene_ids) > 0) {

  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)


  write.csv(as.data.frame(ego), "GO_Enrichment_All.csv")


  kk <- enrichKEGG(gene = gene_ids$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05)


  if(!is.null(ego)) {
    go_df <- as.data.frame(ego)
    top_go <- go_df %>%
      group_by(ONTOLOGY) %>%
      slice_min(p.adjust, n = 5) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))

    p_go <- ggplot(top_go, aes(x = ONTOLOGY, y = Description, size = Count, color = -log10(p.adjust))) +
      geom_point() +
      scale_color_gradient(low = "#077878", high = "#f933ff") +
      theme_bw() +
      labs(title = "GO Enrichment (Top 5)")

    print(p_go)
    ggsave("GO_Bubble_plot.pdf", p_go, width = 10, height = 8)
  }
}


log_data <- log2(norm_data + 1)
z_score_data <- t(scale(t(log_data)))


top_vars <- head(order(apply(z_score_data, 1, var), decreasing = TRUE), 50)
mat_top50 <- z_score_data[top_vars, ]


ann_col <- data.frame(Group = group_list)
rownames(ann_col) <- colnames(mat_top50)
ann_colors <- list(Group = c("Control" = "#0f99b2", "cko" = "#f74ed6"))


pheatmap(mat_top50,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         main = "Top 50 Variable Proteins",
         color = colorRampPalette(c("#077878", "white", "#f933ff"))(100),
         filename = "Heatmap_Top50.pdf", width = 8, height = 10)


group_means <- data.frame(
  Control = rowMeans(log_data[, group_list == "Control"]),
  cko = rowMeans(log_data[, group_list == "cko"])
)


mat_means <- t(scale(t(group_means)))


pheatmap(mat_means,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Group Mean Expression (All Proteins)",
         color = colorRampPalette(c("#31BD8B", "white", "#CF65F0"))(100),
         filename = "Heatmap_GroupMeans.pdf", width = 5, height = 8)

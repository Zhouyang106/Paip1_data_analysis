# 1. Load Libraries
# Data Manipulation & Visualization
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

# Single Cell Analysis
library(Seurat)      # Recommended V5
library(harmony)     # Batch correction
library(monocle)     # Trajectory analysis (Monocle 2)
library(clusterProfiler) # Enrichment analysis
library(org.Mm.eg.db)    # Mouse annotation database

# 2. Set Global Parameters
# Set working directory (Modified for portability)
work_dir <- "D:/LK/"
setwd(work_dir)

# Define Color Palettes (Global)
cols_all <- c("#D6A2D6","#6BC4B8","#FFB77A","#D15F76","#9B8BC6","#89D3C7","#3E8E9F","#F4E76E","#C97B5F")
cols_germ <- c("#9B8BC6","#B8F0E0","#5DD4C1","#2A6A7D","#4A96D7","#F4E76E","#C97B5F")
cols_monocle <- c("#B8F0E0","#5DD4C1","#2A6A7D","#4A96D7") # For Leptotene to Diplotene
# ==============================================================================
# Step 1: QC, Normalization, Integration and Clustering
# ==============================================================================

# 1.1 Data Loading
cko_data <- Read10X(data.dir = "cKO/")
wt_data <- Read10X(data.dir = "WT/")

# Create Seurat Objects
CKO <- CreateSeuratObject(cko_data, project = "cko", min.cells = 3, min.features = 200)
WT <- CreateSeuratObject(wt_data, project = "wt", min.cells = 3, min.features = 200)

# Merge Objects
sc_obj <- merge(WT, y = c(CKO), add.cell.ids = c("WT", "cKO"))
head(colnames(sc_obj))

# 1.2 Quality Control (QC)
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")

# Visualization before filtering
VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filtering
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000)

# 1.3 Normalization & Feature Selection
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)

# Visualization of Variable Features
top10 <- head(VariableFeatures(sc_obj), 10)
plot1 <- VariableFeaturePlot(sc_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave("QC_VariableFeatures.png", plot = plot1 + plot2, width = 12, height = 5, dpi = 300)

# 1.4 Scaling & PCA
all.genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all.genes)
sc_obj <- RunPCA(sc_obj, npcs = 30, verbose = FALSE)

# Determine Dimensionality
ElbowPlot(sc_obj, ndims = 30, reduction = "pca")

# 1.5 UMAP & Clustering
sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:15)
sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.4)

# 1.6 Visualization
sc_obj$orig.ident <- factor(sc_obj$orig.ident, levels = c("wt", "cko"))
p1 <- DimPlot(sc_obj, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(sc_obj, reduction = "umap", label = TRUE)
ggsave("Global_UMAP_Group.pdf", plot = plot_grid(p1, p2), width = 12, height = 5, dpi = 300)

# Save Checkpoint 1
saveRDS(sc_obj, file = "sc_obj_step1_clustered.rds")
# ==============================================================================
# Step 2: Marker Detection & Cell Type Annotation
# ==============================================================================
sc_obj <- readRDS("sc_obj_step1_clustered.rds")

# 2.1 Find Markers (Optional)
sc_obj <- JoinLayers(sc_obj)
markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "global_markers.csv")

# 2.2 Visualization of Canonical Markers

marker_list <- list(
  Sertoli     = c("Ctsl", "Clu", "Cst12", "Cst9", "Defb19", "Ldhb"),
  Leydig      = c("Cyp17a1", "Hsd3b1", "Star", "Hsd3b6", "Fabp3", "Akr1c1"), 
  Myoid       = c("Dcn", "Igfbp7", "Gsn", "Myl9", "Mgp", "Rbbp7"),
  Macrophages = c("Lyz2", "C1qb", "C1qc", "C1qa", "Pf4", "Apoe"),
  SPG         = c("Stra8", "Ptma", "Uchl1", "Dazl", "Crabp1", "Nmt2"),
  SPC         = c("Tex12", "Tex101", "Sycp3", "Insl6", "Spag6l", "Tbpl1"),
  Spermatid   = c("Tex29", "Tex36", "Tbc1d23", "Cst13", "Prm2", "Tfam")
)

# A. FeaturePlots 
# 组合1: Sertoli + Leydig
fp1 <- FeaturePlot(sc_obj, features = c(marker_list$Sertoli, marker_list$Leydig), 
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Sertoli_Leydig.png", plot = fp1, width = 13, height = 10, dpi = 300)

# 组合2: Myoid + Macrophages
fp2 <- FeaturePlot(sc_obj, features = c(marker_list$Myoid, marker_list$Macrophages), 
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Myoid_Macrophage.png", plot = fp2, width = 13, height = 10, dpi = 300)

# 组合3: SPG + SPC
fp3 <- FeaturePlot(sc_obj, features = c(marker_list$SPG, marker_list$SPC), 
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 4)
ggsave("Marker_FeaturePlot_Germ1.png", plot = fp3, width = 13, height = 10, dpi = 300)

# 组合4: Spermatid
fp4 <- FeaturePlot(sc_obj, features = marker_list$Spermatid, 
                   pt.size = 0.5, min.cutoff = "q9", cols = c("gray", "red"), ncol = 3)
ggsave("Marker_FeaturePlot_Germ2.png", plot = fp4, width = 9, height = 10, dpi = 300)

# B. Violin Plots
# 批量绘制保存，使用 lapply 循环 (比写7遍 VlnPlot 更简洁)
# 注意：此时 Idents(sc_obj) 最好是数字编号 cluster，如果是已注释的名字也可以
lapply(names(marker_list), function(cell_type) {
  p <- VlnPlot(sc_obj, features = marker_list[[cell_type]], pt.size = 0, ncol = 3) # pt.size=0 去除噪点
  ggsave(paste0("Marker_VlnPlot_", cell_type, ".png"), plot = p, width = 15, height = 10, dpi = 300)
})

# C. DotPlot (气泡图 - 汇总所有基因)
# 展平所有基因为一个向量
all_markers <- unlist(marker_list)

# 气泡图绘制
dot1 <- DotPlot(sc_obj, features = all_markers, group.by = "seurat_clusters") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), # 优化X轴标签
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "right"
  ) + 
  scale_color_gradientn(values = seq(0, 1, 0.2), 
                        colours = c('#330066', '#336699', '#66CC66', '#FFCC33')) # 你原先的配色

ggsave("Marker_DotPlot_All.pdf", plot = dot1, width = 14, height = 6, dpi = 300)

# 2.3 Assign Cell Types
# Fixed spelling: Leyding -> Leydig
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
sc_obj$cell_type <- Idents(sc_obj) # Save annotation to metadata

# Reorder levels for plotting
Idents(sc_obj) <- factor(Idents(sc_obj), levels = c("Sertoli cell","Leydig cell","Myoid cell","Macrophages cell", 
                                                    "SPG", "Early SPC", "Late SPC","Round Spermatid","Elongating Spermatid"))

# Save Checkpoint 2
saveRDS(sc_obj, file = "sc_obj_step2_annotated.rds")
# ==============================================================================
# Step 3: Germ Cell Sub-clustering and Harmony Integration
# ==============================================================================
sc_obj <- readRDS("sc_obj_step2_annotated.rds")

# 3.1 Subset Germ Cells
germ_cells <- subset(sc_obj, idents = c("SPG", "Early SPC", "Late SPC", "Round Spermatid", "Elongating Spermatid"))

# 3.2 Re-process Germ Cells
# Create new object to remove old slots/caches
germ_sce <- CreateSeuratObject(counts = GetAssayData(germ_cells, assay = "RNA", layer = 'counts'), 
                               meta.data = germ_cells@meta.data)

germ_sce <- NormalizeData(germ_sce) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

# 3.3 Harmony Integration (Batch Correction)
germ_sce <- RunHarmony(germ_sce, group.by.var = "orig.ident")

# 3.4 Clustering based on Harmony
germ_sce <- FindNeighbors(germ_sce, reduction = "harmony", dims = 1:20)
germ_sce <- FindClusters(germ_sce, resolution = 0.5)
germ_sce <- RunUMAP(germ_sce, reduction = "harmony", dims = 1:20)

# 3.5 Fine-grained Annotation of Germ Cells
# Logic optimized from loop to direct assignment
current_ids <- 0:15
# Define mapping based on user logic (Cluster ID -> Cell Type)
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
# Check for unassigned clusters
all_clusters <- levels(germ_sce)
names(new_germ_ids) <- as.character(names(new_germ_ids))

# Apply annotation carefully
germ_sce <- RenameIdents(germ_sce, new_germ_ids)
germ_sce$celltype <- Idents(germ_sce)

# Reorder
germ_sce$celltype <- factor(germ_sce$celltype, levels = c("Sertoli cell","SPG","Leptotene","Zygotene","Pachytene","Diplotene","Round Std","Elong Std"))
Idents(germ_sce) <- germ_sce$celltype

# Save Checkpoint 3
save(germ_sce, file = "germ_sce_step3_annotated.Rdata")

# 3.6 Visualization (Sunburst & Barplots)
library(ggforce) # Required for Sunburst plot

# 1. Data Preparation: Calculate Proportions
# 提取细胞类型和分组信息
meta_df <- germ_sce@meta.data[, c("celltype", "orig.ident")]

# 统计各组各细胞类型的数量
cell_stats <- meta_df %>%
  group_by(orig.ident, celltype) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(Total = sum(Count)) %>%
  mutate(Proportion = Count / Total * 100) # 计算百分比

# 筛选感兴趣的减数分裂相关细胞 (用于绘图)
target_cells <- c("Leptotene", "Zygotene", "Pachytene", "Diplotene")
plot_data <- cell_stats %>%
  filter(celltype %in% target_cells) %>%
  mutate(celltype = factor(celltype, levels = target_cells)) # 固定因子顺序

# 定义配色 (与 Monocle 轨迹图保持一致)
cols_meiosis <- c("#B8F0E0", "#5DD4C1", "#2A6A7D", "#4A96D7")
names(cols_meiosis) <- target_cells

# 2. Stacked Bar Plot (总体比例)
p_bar_all <- ggplot(plot_data, aes(x = orig.ident, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.5, color = "white") +
  scale_fill_manual(values = cols_meiosis) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), expand = c(0, 0), limits = c(0, 105)) +
  scale_x_discrete(labels = toupper) + # X轴标签大写
  labs(x = NULL, y = "Percentage of Cell Type (%)", title = "Cell Proportion") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right"
  )

ggsave("Proportion_Barplot_All.pdf", plot = p_bar_all, width = 6, height = 5)

# 3. Faceted Bar Plot (分面展示，突出显示 Pachytene)
# 准备数据：标记是否为 Pachytene 以便加粗显示（可选），此处保留原逻辑
p_bar_facet <- ggplot(plot_data, aes(x = Proportion, y = orig.ident, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(celltype ~ ., scales = "free_y", space = "free") + # 纵向分面
  scale_fill_manual(values = cols_meiosis, guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
            hjust = -0.2, size = 3.5) + # 添加百分比标签
  labs(x = "Proportion (%)", y = NULL) +
  theme_classic() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) # 留出右侧空间给标签

ggsave("Proportion_Barplot_Faceted.pdf", plot = p_bar_facet, width = 7, height = 6)

# 4. Sunburst Plot (旭日图)
# 准备旭日图所需的坐标数据
sunburst_data <- plot_data %>%
  mutate(orig.ident = as.factor(orig.ident)) %>%
  arrange(orig.ident, celltype) %>% # 必须先排序
  group_by(orig.ident) %>%
  mutate(
    # 计算累积百分比用于画圆弧
    ymax = cumsum(Proportion),
    ymin = lag(ymax, default = 0),
    # 计算标签位置 (中间点)
    label_pos = (ymin + ymax) / 2
  )
sunburst_data$r0 <- as.numeric(sunburst_data$orig.ident) - 0.4
sunburst_data$r  <- as.numeric(sunburst_data$orig.ident) + 0.4

# 绘图
p_sunburst <- ggplot(sunburst_data) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, 
                   r0 = r0, r = r, 
                   start = ymin * 2 * pi / 100, # 转换为弧度
                   end = ymax * 2 * pi / 100, 
                   fill = celltype),
               color = "white", size = 0.5) +
  # 添加文字标签 (仅针对 Pachytene 或比例大于5%的)
  geom_text(data = subset(sunburst_data, celltype == "Pachytene" | Proportion > 5),
            aes(x = (r0 + 0.45) * sin(label_pos * 2 * pi / 100),
                y = (r0 + 0.45) * cos(label_pos * 2 * pi / 100),
                label = sprintf("%.1f%%", Proportion)),
            size = 3, fontface = "bold") +
  # 中心添加组别标签
  annotate("text", x = 0, y = 0, label = "Group", fontface = "bold") +
  scale_fill_manual(values = cols_meiosis) +
  coord_equal() +
  theme_void() + # 移除坐标轴和背景
  theme(legend.position = "right", legend.title = element_blank())

ggsave("Proportion_Sunburst.pdf", plot = p_sunburst, width = 8, height = 6)
# ==============================================================================
# Step 4: Trajectory Analysis with Monocle 2
# ==============================================================================
load("germ_sce_step3_annotated.Rdata")

# 4.1 Setup Monocle Object
target_cells <- c("Leptotene", "Zygotene", "Pachytene", "Diplotene")
seurat_sub <- subset(germ_sce, idents = target_cells)

# Extract Data (Seurat V5 compatible)
data <- GetAssayData(seurat_sub, assay = "RNA", layer = 'counts')
pd <- new("AnnotatedDataFrame", data = seurat_sub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# 4.2 Preprocessing
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# 4.3 Feature Selection (Ordering Genes)
# Using Differential Gene Test over Cell Types
diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~celltype", cores = 4) # Adjust cores
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)

# 4.4 Dimensionality Reduction & Ordering
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)

# Define Root State (Automated)
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"Leptotene"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return(1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

# 4.5 Visualization
p1 <- plot_cell_trajectory(cds, color_by = "celltype") + scale_color_manual(values = cols_monocle)
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave("Monocle_Trajectory_Combined.pdf", plot = p1 + p2, width = 16, height = 8)

# Save Monocle Object
save(cds, file = "monocle_cds_final.RData")

# ==============================================================================
# Step 5: Differential Expression & Enrichment Analysis (GO)
# ==============================================================================
# Focus: Pachytene Stage (WT vs cKO)

# 5.1 Find Markers
Meio <- subset(germ_sce, idents = "Pachytene")
Idents(Meio) <- "orig.ident" # Switch to Group ID

DE_results <- FindMarkers(
  Meio,
  ident.1 = "cko",
  ident.2 = "wt",
  logfc.threshold = 0.5,
  min.pct = 0.1,
  test.use = "wilcox"
)

# Filter Significant Genes
DE_results_fil <- DE_results %>% 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

write.csv(DE_results_fil, "DE_results_Pachytene.csv")

# 5.2 GO Enrichment
gene_list <- rownames(DE_results_fil)
eg <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
eg <- na.omit(eg)

# Run Enrichment
go_res <- enrichGO(
  gene = eg$ENTREZID,
  OrgDb = "org.Mm.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 5.3 Visualization (Barplot function)
library(stringr)

# 1. Process Data
# 将 enrichResult 对象转换为数据框
if (!exists("go_res")) stop("Error: 'go_res' object not found. Please run Step 5.2 first.")
go_df <- as.data.frame(go_res)

# 分别提取 BP, CC, MF 前 10 个条目 (按 qvalue/p.adjust 排序)
# 辅助函数：提取并处理描述文本
process_go_data <- function(df, ontology_type, top_n = 10) {
  subset_df <- df[df$ONTOLOGY == ontology_type, ]
  if (nrow(subset_df) == 0) return(NULL)
  
  # 取前 N 个
  subset_df <- head(subset_df, top_n)
  
  # 截断过长的描述文本 (保留右侧或中间可能更合适，这里保留原逻辑 side="right")
  subset_df$Description_short <- str_trunc(subset_df$Description, width = 50, side = "right")
  
  # 设置因子水平，保证画图时顺序是固定的 (倒序，让 barplot 从上到下显示)
  subset_df$Description_short <- factor(subset_df$Description_short, levels = rev(subset_df$Description_short))
  
  return(subset_df)
}

data_BP <- process_go_data(go_df, "BP")
data_CC <- process_go_data(go_df, "CC")
data_MF <- process_go_data(go_df, "MF")

# 2. Define Plotting Theme & Function
# 自定义通用主题
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11, color = "black"),
  plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  panel.grid.minor = element_blank()
)

# 绘图函数
plot_GO_bar <- function(data, title_prefix, color_palette = "Blues") {
  if (is.null(data)) {
    warning(paste("No data for", title_prefix))
    return(NULL)
  }
  
  ggplot(data, aes(x = Count, y = Description_short, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    # 自动换行，防止某些标签依然过长
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) + 
    # 颜色渐变
    scale_fill_distiller(palette = color_palette, direction = 1, name = "-log10(padj)") +
    labs(x = "Gene Number", y = NULL,
         title = paste0(title_prefix, " Enrichment")) +
    theme_bw() +
    mytheme
}

# 3. Generate and Save Plots
# 绘制 BP (生物学过程) - 使用橙色系
p_BP <- plot_GO_bar(data_BP, "Biological Process", "Oranges")
if (!is.null(p_BP)) {
  print(p_BP)
  ggsave("GO_Enrichment_BP.pdf", plot = p_BP, width = 9, height = 7)
}

# 绘制 CC (细胞组分) - 使用红色系
p_CC <- plot_GO_bar(data_CC, "Cellular Component", "Reds")
if (!is.null(p_CC)) {
  print(p_CC)
  ggsave("GO_Enrichment_CC.pdf", plot = p_CC, width = 9, height = 7)
}

# 绘制 MF (分子功能) - 使用蓝色系
p_MF <- plot_GO_bar(data_MF, "Molecular Function", "Blues")
if (!is.null(p_MF)) {
  print(p_MF)
  ggsave("GO_Enrichment_MF.pdf", plot = p_MF, width = 9, height = 7)
}

# 4. (Optional) Custom Color Gradient Version
if (!is.null(data_BP)) {
  p_custom <- ggplot(data_BP, aes(x = Count, y = Description_short, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +
    # 自定义渐变色
    scale_fill_gradient(low = "#077878", high = "#f933ff", name = "-log10(padj)") + 
    labs(x = "Gene Number", y = NULL, title = "BP Enrichment (Custom Color)") +
    theme_bw() +
    mytheme
  
  ggsave("GO_Enrichment_BP_CustomColor.pdf", plot = p_custom, width = 9, height = 7)
}


##cross
# ==============================================================================
# Step 0: Environment Setup, Libraries, and Helper Functions
# ==============================================================================

# 1. Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(readxl)
library(stringr)

# 2. Define Helper Functions (To be used later)

# Function: Create Expression Matrix for Cross-species Heatmap
create_seurat_expression_matrix <- function(seurat_objs, genes) {
  expression_list <- list()
  
  for(species in names(seurat_objs)) {
    obj <- seurat_objs[[species]]
    
    # Check gene existence
    genes_present <- intersect(genes, rownames(obj))
    missing_genes <- setdiff(genes, rownames(obj))
    if(length(missing_genes) > 0) {
      warning(paste(species, "missing genes:", paste(missing_genes, collapse=", ")))
    }
    
    # Get Data (Compatible with Seurat V5)
    # Using LayerData to safely extract normalized data
    data_matrix <- LayerData(obj, assay = "RNA", layer = "data")
    
    # Calculate Mean Expression
    avg_expr <- rowMeans(data_matrix[genes_present, , drop=FALSE], na.rm = TRUE)
    
    # Fill missing genes with 0
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

# Function: Clean Factors in Seurat Object
clean_seurat_factors <- function(obj) {
  obj@meta.data$group <- droplevels(as.factor(obj@meta.data$group))
  obj@meta.data$celltype <- droplevels(as.factor(obj@meta.data$celltype))
  return(obj)
}

# Function: Convert Human/TF symbols (All Caps) to Mouse symbols (Title Case)
to_mouse_symbol <- function(genes) {
  return(str_to_title(genes)) 
}

# Function: Calculate Expression Difference (SPC vs Std)
calculate_expression_diff <- function(obj_spc, obj_std, gene_list, label) {
  # Calculate Mean SPC
  avg_spc <- AverageExpression(obj_spc, features = gene_list, verbose = FALSE)$RNA
  avg_spc <- rowMeans(avg_spc)
  
  # Calculate Mean Std
  avg_std <- AverageExpression(obj_std, features = gene_list, verbose = FALSE)$RNA
  avg_std <- rowMeans(avg_std)
  
  # Merge
  res <- data.frame(
    Gene = names(avg_spc),
    Mean_SPC = as.numeric(avg_spc),
    Mean_Std = as.numeric(avg_std[names(avg_spc)])
  )
  
  # Calculate Difference
  res$Difference <- res$Mean_SPC - res$Mean_Std
  res$Species <- label
  
  # Sort
  res <- res %>% arrange(desc(abs(Difference)))
  
  return(res)
}
# ==============================================================================
# Step 1: Data Loading & Cross-Species Heatmap
# ==============================================================================

# 1.1 Load Sheep Data
setwd("D:/Y/") 
load("sce2_dim20_reso0.3_zhushi.Rdata")
Tsheep <- sce
Idents(Tsheep) <- "group"
rm(sce)

# Clean Sheep Genes (remove "gene-" prefix)
clean_genes_sheep <- gsub("^gene-", "", rownames(Tsheep))
rownames(Tsheep) <- clean_genes_sheep
rownames(Tsheep@assays$RNA) <- clean_genes_sheep

# 1.2 Load Mouse Data
setwd("D:/Y/mouse/") 
load("Tmouse.Rdata")
Tmouse <- data3
Tmouse$group <- Tmouse$orig.ident
Idents(Tmouse) <- "group"
Tmouse[["RNA"]] <- JoinLayers(Tmouse[["RNA"]]) # Seurat V5 Fix
rm(data3)

# Standardize Mouse Genes to Uppercase (for matching)
clean_genes_mouse <- toupper(rownames(Tmouse))
rownames(Tmouse) <- clean_genes_mouse
rownames(Tmouse@assays$RNA) <- clean_genes_mouse

# 1.3 Load Human Data
setwd("D:/Y/human/") 
load("Thuman.Rdata")
Thuman <- data3
Idents(Thuman) <- "group"
Thuman[["RNA"]] <- JoinLayers(Thuman[["RNA"]])
rm(data3)

# 1.4 Load Cattle Data (Optional, loaded but not used in list below based on user code)
setwd("D:/LK/Cattle/") 
load("Tcattle.Rdata")
Tcattle <- immune.combined
Tcattle[["RNA"]] <- JoinLayers(Tcattle[["RNA"]])
rm(immune.combined)

# ------------------------------------------------------------------------------
# 1.5 Generate Heatmap
# ------------------------------------------------------------------------------
setwd("D:/LK/Cross/") 

# Define Object List
seurat_objects <- list(
  Human = Thuman,
  Mouse = Tmouse,
  Sheep = Tsheep
)

# Load Gene Set (Ensure 'Three_species_share_genes' exists in env)
# if(!exists("Three_species_share_genes")) stop("Please load Three_species_share_genes object first!")
new_gene <- Three_species_share_genes$gene
selected_genes <- new_gene

# Generate Matrix
expression_matrix <- create_seurat_expression_matrix(seurat_objects, selected_genes)

cat("Expression matrix generated.\n")

# Log Transformation
expression_matrix_log <- log2(expression_matrix + 1)

# Plot Heatmap
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
# ==============================================================================
# Step 2: Sub-clustering Data Load & Differential Expression (DE)
# ==============================================================================

# 2.1 Load Sub-clustered Data (4 stages)
# Sheep
setwd("D:/Y/") 
load("g_sce_subset_dim20_reso0.6_4_zhushi.Rdata")
sheep <- g_sce_subset
# Fix Sheep Names
clean_genes <- gsub("^gene-", "", rownames(sheep))
rownames(sheep) <- clean_genes
rownames(sheep@assays$RNA) <- clean_genes
rm(g_sce_subset)

# Mouse
setwd("D:/Y/mouse/") 
load("g_sce_dim20_reso1.2_4_zhushi.Rdata")
mouse <- g_sce
Idents(mouse) <- "celltype"
mouse@meta.data$group <- mouse@meta.data$orig.ident
rm(g_sce)

# Human
setwd("D:/Y/human/") 
load("g_sce_dim20_reso1.2_4_zhushi.Rdata")
human <- g_sce
Idents(human) <- "celltype"
rm(g_sce)

# 2.2 Run FindMarkers
setwd("D:/LK/Cross/") 

# Human DE (11yo vs 7yo)
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
write.csv(as.data.frame(human_DE_results), file = "human_DE_results.csv") # Fixed typo 'huamn'

# Mouse DE (PND8 vs PND5)
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

# Sheep DE (3m vs 1m)
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
# ==============================================================================
# Step 3: Transcription Factor (TF) Intensity Analysis
# ==============================================================================

# Ensure 'Post_tran_TF' exists
# if(!exists("Post_tran_TF")) stop("Please load Post_tran_TF object first!")
TF <- Post_tran_TF$gene

# 3.1 Define Targets and Subset Data
target_spc_types1 <- c("Preleptotene", "Leptotene/Zygotene", "Pachytene", "Diplotene")
target_spc_types2 <- c("Round Std", "Elong Std")

# Human Subset
target_group_human1 <- c("11yo", "13yo")
target_group_human2 <- c("14yo", "25yo")
human_SPC <- subset(human, subset = group %in% target_group_human1 & celltype %in% target_spc_types1)
human_Std <- subset(human, subset = group %in% target_group_human2 & celltype %in% target_spc_types2)

# Mouse Subset
target_group_mouse1 <- c("PND8", "PND10", "PND14")
target_group_mouse2 <- c("PND30", "PND35")
mouse_SPC <- subset(mouse, subset = group %in% target_group_mouse1 & celltype %in% target_spc_types1)
mouse_Std <- subset(mouse, subset = group %in% target_group_mouse2 & celltype %in% target_spc_types2)

# Sheep Subset
target_group_sheep1 <- c("3m", "5m")
target_group_sheep2 <- c("6m", "12m")
sheep_SPC <- subset(sheep, subset = group %in% target_group_sheep1 & celltype %in% target_spc_types1)
sheep_Std <- subset(sheep, subset = group %in% target_group_sheep2 & celltype %in% target_spc_types2)

# Clean factors
human_SPC <- clean_seurat_factors(human_SPC)
human_Std <- clean_seurat_factors(human_Std)
mouse_SPC <- clean_seurat_factors(mouse_SPC)
mouse_Std <- clean_seurat_factors(mouse_Std)
sheep_SPC <- clean_seurat_factors(sheep_SPC)
sheep_Std <- clean_seurat_factors(sheep_Std)

# 3.2 Find TF Intersections
# Human
intersect_human <- intersect(TF, rownames(human_SPC))
intersect_human <- intersect(intersect_human, rownames(human_Std))

# Sheep
intersect_sheep <- intersect(TF, rownames(sheep_SPC))
intersect_sheep <- intersect(intersect_sheep, rownames(sheep_Std))

# Mouse (Handle Case Sensitivity)
TF_mouse_format <- to_mouse_symbol(TF)
intersect_mouse <- intersect(TF_mouse_format, rownames(mouse_SPC))
intersect_mouse <- intersect(intersect_mouse, rownames(mouse_Std))

cat("Intersections found: Human:", length(intersect_human), "| Sheep:", length(intersect_sheep), "| Mouse:", length(intersect_mouse), "\n")

# 3.3 Calculate Differences
diff_human <- calculate_expression_diff(human_SPC, human_Std, intersect_human, "Human")
diff_mouse <- calculate_expression_diff(mouse_SPC, mouse_Std, intersect_mouse, "Mouse")
diff_sheep <- calculate_expression_diff(sheep_SPC, sheep_Std, intersect_sheep, "Sheep")

# Save tables
write.csv(diff_human, "TF_Diff_Human.csv")
write.csv(diff_mouse, "TF_Diff_Mouse.csv")
write.csv(diff_sheep, "TF_Diff_Sheep.csv")

# 3.4 Visualization: Bar Plot (TF Difference)
# Normalize Gene Names for Plotting
diff_human$Gene_Std_Name <- toupper(diff_human$Gene)
diff_mouse$Gene_Std_Name <- toupper(diff_mouse$Gene)
diff_sheep$Gene_Std_Name <- toupper(diff_sheep$Gene)

diff_combined <- bind_rows(diff_human, diff_mouse, diff_sheep)
diff_combined$Species <- factor(diff_combined$Species, levels = c("Human", "Mouse", "Sheep"))
diff_combined$Direction <- ifelse(diff_combined$Difference > 0, "Higher in SPC", "Higher in Std")
diff_combined$Direction <- factor(diff_combined$Direction, levels = c("Higher in Std", "Higher in SPC"))

# Order genes based on total absolute difference
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

# 3.5 Visualization: Slope Chart (General)
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

ggsave("SlopeChart_General.pdf", plot = p3, width = 16, height = 7, dpi = 300) # Replaced Chinese filename

# 3.6 Visualization: Slope Chart (Highlighted Genes)
problem_genes <- c("DDX5", "DAZL", "PABPC1", "PABPC1", "PAIP2", "RBMXL2", "PABPC2", "PAIP1", "MOV10L1", "DDX4", "YBX2")

p3_high <- ggplot(plot_data, aes(x = Condition, y = Expression, group = Gene_Std_Name)) +
  geom_line(aes(color = Direction), size = 0.8, alpha = 0.1) + # Background lines
  geom_line(data = plot_data %>% filter(Gene_Std_Name %in% problem_genes), 
            aes(color = Direction), size = 1.2, alpha = 1) + # Highlight lines
  geom_point(size = 2, alpha = 0.2) +
  geom_point(data = plot_data %>% filter(Gene_Std_Name %in% problem_genes), 
             size = 3, color = "black") + # Highlight points
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

ggsave("SlopeChart_Highlighted.pdf", plot = p3_high, width = 16, height = 7, dpi = 300) # Replaced Chinese filename
# ==============================================================================
# Step 4: GO Enrichment Visualization (Human, Mouse, Sheep)
# ==============================================================================

# Helper function to process and plot GO data
process_and_plot_go <- function(file_path, suffix) {
  
  # 1. Read Data
  tryCatch({
    df <- readxl::read_excel(file_path)
  }, error = function(e) {
    stop(paste("Cannot read file:", file_path))
  })
  
  # 2. Process Data
  df_plot <- df %>%
    mutate(p.adjust = as.numeric(`p.adjust`)) %>%
    mutate(logP = -log10(p.adjust)) %>%
    mutate(Ontology = gsub("GO:", "", source))
  
  # 3. Separation Logic
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
    mutate(Description = factor(Description, levels = Description)) %>% # For Faceted Plot
    ungroup()
  
  # 4. Faceted Plot (Combined)
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
  
  # 5. Separate Plots Loop (BP / MF)
  ontologies_to_plot <- c("Molecular Function (MF)", "Biological Process (BP)")
  
  for (current_ontology in ontologies_to_plot) {
    df_subset <- df_separated %>%
      filter(Main_Ontology == current_ontology) %>%
      group_by(Main_Ontology) %>%
      mutate(Description = factor(Description, levels = rev(Description))) %>% # Rev for single plot
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

# ------------------------------------------------------------------------------
# Run GO Plotting for each Species
# ------------------------------------------------------------------------------
setwd("D:/LK/Cross/") 

# Run for Human
if(file.exists("term_human.xlsx")) {
  process_and_plot_go("term_human.xlsx", "human")
}

# Run for Mouse
if(file.exists("term_mouse.xlsx")) {
  process_and_plot_go("term_mouse.xlsx", "mouse")
}

# Run for Sheep (Assuming filename follows pattern, if exists)
if(file.exists("term_sheep.xlsx")) {
  process_and_plot_go("term_sheep.xlsx", "sheep")
}

##Proteome
# ------------------------------------------------------------------------------
# 1. 环境配置与库加载
# ------------------------------------------------------------------------------
# 加载库
library(readxl)          # 读取Excel
library(dplyr)           # 数据处理
library(tidyr)           # 数据转换
library(stringr)         # 字符串处理
library(ggplot2)         # 绘图核心
library(ggrepel)         # 标签防重叠
library(limma)           # 差异分析
library(pcaMethods)      # PCA分析
library(pheatmap)        # 热图
library(clusterProfiler) # 富集分析
library(org.Mm.eg.db)    # 小鼠基因注释
library(enrichplot)      # 富集绘图
library(RColorBrewer)    # 调色板

# 设置工作目录
setwd("D:/蛋白质组/danbaizhifenxi/")

# 定义统一配色
my_colors <- c("Control" = "#0f99b2", "cko" = "#f74ed6", "CKO" = "#f74ed6") # 注意cko大小写统一
volcano_colors <- c("Up" = "#f933ff", "Down" = "#077878", "NS" = "grey")

# ------------------------------------------------------------------------------
# 2. 数据导入与预处理
# ------------------------------------------------------------------------------
# 读取数据 (假设第一列是Protein ID，后续是表达量)
raw_data <- read_excel("data1.xlsx")
protein_ids <- raw_data[[1]]
expr_data <- as.matrix(raw_data[-1])
rownames(expr_data) <- protein_ids

# 缺失值填充 (使用最小值的一半填充NA)
expr_data[is.na(expr_data)] <- min(expr_data, na.rm = TRUE) * 0.5

# 归一化 (Quantile normalization)
norm_data <- normalizeBetweenArrays(expr_data, method = "quantile")

# 转换ID：UniProt -> Gene Symbol (便于后续分析)
gene_symbols <- mapIds(org.Mm.eg.db, keys = rownames(norm_data), column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
# 如果找不到Symbol保留原ID
rownames(norm_data) <- ifelse(is.na(gene_symbols), rownames(norm_data), gene_symbols)

head(norm_data)

# ------------------------------------------------------------------------------
# 3. 质控分析：PCA
# ------------------------------------------------------------------------------
# 准备分组信息 (假设前3个是Control，后3个是cko)
group_list <- factor(c(rep("Control", 3), rep("cko", 3)), levels = c("Control", "cko"))

# 运行PCA
pca_res <- pca(t(norm_data), method = "svd", nPcs = 2, scale = "uv")
pca_scores <- as.data.frame(scores(pca_res))
pca_scores$group <- group_list
pca_scores$sample <- colnames(norm_data)

# 绘制PCA图
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

# ------------------------------------------------------------------------------
# 4. 差异分析 (Limma)
# ------------------------------------------------------------------------------
# 创建设计矩阵
design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list)

# 线性模型拟合
fit <- lmFit(norm_data, design)

# 创建对比 (cko - Control)
contrast.matrix <- makeContrasts(cko - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 获取结果
deg_results <- topTable(fit2, adjust = "BH", sort.by = "P", number = Inf)
deg_results$Protein <- rownames(deg_results)

# 标记显著性 (FC > 1.2 或 1.5, P < 0.05，根据需求调整)
logFC_cutoff <- 1  # log2(2) = 1
pval_cutoff <- 0.05

deg_results$Significance <- case_when(
  deg_results$adj.P.Val < pval_cutoff & deg_results$logFC > logFC_cutoff ~ "Up",
  deg_results$adj.P.Val < pval_cutoff & deg_results$logFC < -logFC_cutoff ~ "Down",
  TRUE ~ "NS"
)

# 保存差异表格
write.csv(deg_results, "DE_results_Limma.csv")

# ------------------------------------------------------------------------------
# 5. 可视化：火山图
# ------------------------------------------------------------------------------
# 统计数量
up_n <- sum(deg_results$Significance == "Up")
down_n <- sum(deg_results$Significance == "Down")

# 添加标签 (Top显著的20个)
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

# ------------------------------------------------------------------------------
# 6. 特定基因集分析 (Violin Plot)
# ------------------------------------------------------------------------------
# 注意：请确保 YBX2_PAIP1_bindingRNA 对象已加载，或在此处定义 gene_list
# gene_list <- c("GeneA", "GeneB", ...) 
if(exists("YBX2_PAIP1_bindingRNA")) {
  gene_list <- YBX2_PAIP1_bindingRNA$Symbol
  
  # 提取特定基因表达量
  target_expr <- norm_data[rownames(norm_data) %in% gene_list, , drop=FALSE]
  
  if(nrow(target_expr) > 0) {
    # 转换为长格式
    df_long <- as.data.frame(target_expr) %>%
      mutate(Protein = rownames(.)) %>%
      pivot_longer(-Protein, names_to = "Sample", values_to = "Expression") %>%
      mutate(Group = ifelse(grepl("Control", Sample), "Control", "cko"))
    
    # 总体小提琴图
    p_vln <- ggplot(df_long, aes(x = Group, y = Expression, fill = Group)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white") +
      labs(title = "Expression of Target Gene Set") +
      scale_fill_manual(values = my_colors) +
      theme_bw()
    
    ggsave("TargetSet_Violin.pdf", p_vln, width = 6, height = 5)
  }
}

# ------------------------------------------------------------------------------
# 7. 富集分析 (GO & KEGG)
# ------------------------------------------------------------------------------
# 提取显著差异蛋白并ID转换
sig_genes <- subset(deg_results, Significance != "NS")$Protein

# 转换ID为 ENTREZID
gene_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

if(nrow(gene_ids) > 0) {
  # GO 分析 (一次性运行所有 Ontology)
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Mm.eg.db,
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
  
  # 保存GO结果
  write.csv(as.data.frame(ego), "GO_Enrichment_All.csv")
  
  # KEGG 分析
  kk <- enrichKEGG(gene = gene_ids$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05)
  
  # 气泡图绘制 (GO Top 5 per category)
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

# ------------------------------------------------------------------------------
# 8. 全局热图 (Heatmap)
# ------------------------------------------------------------------------------
# 8.1 样本层面的 Top 50 变异蛋白热图
# Z-score 标准化
log_data <- log2(norm_data + 1)
z_score_data <- t(scale(t(log_data)))

# 选取方差最大的前50个蛋白
top_vars <- head(order(apply(z_score_data, 1, var), decreasing = TRUE), 50)
mat_top50 <- z_score_data[top_vars, ]

# 样本注释
ann_col <- data.frame(Group = group_list)
rownames(ann_col) <- colnames(mat_top50)
ann_colors <- list(Group = c("Control" = "#0f99b2", "cko" = "#f74ed6"))

# 绘制
pheatmap(mat_top50,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         main = "Top 50 Variable Proteins",
         color = colorRampPalette(c("#077878", "white", "#f933ff"))(100),
         filename = "Heatmap_Top50.pdf", width = 8, height = 10)

# 8.2 组均值热图 (所有蛋白)
# 计算组均值
group_means <- data.frame(
  Control = rowMeans(log_data[, group_list == "Control"]),
  cko = rowMeans(log_data[, group_list == "cko"])
)

# Z-score 均值数据
mat_means <- t(scale(t(group_means)))

# 绘制
pheatmap(mat_means,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Group Mean Expression (All Proteins)",
         color = colorRampPalette(c("#31BD8B", "white", "#CF65F0"))(100),
         filename = "Heatmap_GroupMeans.pdf", width = 5, height = 8)




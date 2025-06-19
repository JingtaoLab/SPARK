library(ggrepel)
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(UCell)
library(SingleR)
library(limma)
library(stringr)
library(jsonlite) # 
library(org.Hs.eg.db)
library(patchwork)
library(presto)
library(scRNAtoolVis)
library(corrplot)
library(ggsci)
library(BiocParallel)

### 方法内判断输入样本的格式
readData <- function(dirPath, fileList, project) {
  filesCount <- length(fileList)
  if (filesCount == 0) {
    print("样本路径读取为空")
    return(NULL)
  }else if (filesCount == 1) {
    obj.counts <- NULL
    if (endsWith(fileList, ".csv")) {
      obj.counts <- read.csv(fileList, header = T, row.names = 1)
    }else if (endsWith(fileList, ".txt")) {
      obj.counts <- read.table(fileList, header = T)
    }
    print(fileList)
    obj <- CreateSeuratObject(obj.counts, project = project, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^MT\\.")
    return(obj)
  }else if (filesCount == 3) {
    print(dirPath)
    obj.counts <- Read10X(dirPath)
    obj <- CreateSeuratObject(obj.counts, project = project, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^MT\\.")
    return(obj)
  }
  return(NULL)
}

### 判断输入样本的格式
dir <- inputPath
fileList <- list.files(dir, include.dirs = F, full.names = TRUE, recursive = F)
pbmc_expr <- readData(dir, fileList, "expr")

allCellCount <- jsonlite::toJSON(data.frame(gene = dim(pbmc_expr)[1], cell = dim(pbmc_expr)[2]), pretty = F)
write_json(allCellCount, paste0(outputPath, "cellCount.json", seq = ""))

# 控制变量组
dir <- controlPath
fileList <- list.files(dir, include.dirs = F, full.names = TRUE, recursive = F)
pbmc_ctrl <- readData(dir, fileList, "ctrl")

pbmc <- merge(pbmc_expr, pbmc_ctrl, add.cell.ids = c("expr", "ctrl"))

## 添加分组信息
pbmc$sample = stringr::str_split_fixed(colnames(pbmc), "_", n = 2)[, 1]

if (isSaveData) {
  # 临时保存pbmc数据
  save(pbmc, pbmc_ctrl, pbmc_expr, file = paste0(outputPath, "pbmc.rds", seq = ""));
}
rm(pbmc_ctrl, pbmc_expr)
gc()

### 质量控制 QC
nFeature_RNA_value <- round(as.matrix(quantile(pbmc$nFeature_RNA, 96 / 100))[1], 2)
nCount_RNA_value <- round(as.matrix(quantile(pbmc$nCount_RNA, 96 / 100))[1], 2)
percent_mt_value <- round(as.matrix(quantile(pbmc$percent.mt, 90 / 100))[1], 2)

p1 <- VlnPlot(pbmc, features = "percent.mt") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = percent_mt_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = percent_mt_value, label = percent_mt_value)
p2 <- VlnPlot(pbmc, features = "nCount_RNA") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = nCount_RNA_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = nCount_RNA_value, label = nCount_RNA_value)
p3 <- VlnPlot(pbmc, features = "nFeature_RNA") &
  geom_hline(linetype = 'dotdash', col = 'red', yintercept = nFeature_RNA_value, size = 1) &
  NoLegend() &
  annotate(geom = "label", x = 2, y = nFeature_RNA_value, label = nFeature_RNA_value)

pQC <- wrap_plots(p1, p2, p3, ncol = 3)

#### filter
aging <- subset(pbmc, nCount_RNA >= 800 &
                  nCount_RNA < nCount_RNA_value &
                  percent.mt <= percent_mt_value &
                  nFeature_RNA < nFeature_RNA_value &
                  nFeature_RNA > 500)

##data processing
aging <- NormalizeData(aging, normalization.method = "LogNormalize", scale.factor = 10000) ####默认参数
# PCA
aging <- FindVariableFeatures(aging, selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA()  
ElbowPlot(aging, ndims = 40)


if (isReduction) {
  # 去批次
  library(harmony)
  aging <- RunHarmony(aging, group.by.vars = "sample")
  #降维聚类
  aging <- FindNeighbors(aging, reduction = "harmony", dims = 1:dims) %>% 
    FindClusters(resolution = 0.2)
  aging <- RunUMAP(aging, reduction = "harmony", dims = 1:dims, label = T) %>% 
    RunTSNE(reduction = "harmony", dims = 1:dims, label = T)
}else {
  aging <- FindNeighbors(aging, dims = 1:20) %>% RunUMAP(dims = 1:20)
  aging <- FindClusters(aging, resolution = 0.2) 
}

##plot
DimPlot(aging, reduction = "umap", group.by = c("sample")) 
aging.markers=FindAllMarkers(aging,only.pos = T,assay = "RNA",logfc.threshold = 0.25)

# read nreference
young1 <- readRDS(referenceFilePath)

test_label <- t(FetchData(young1, vars = c("ident")))

test_data <- as.data.frame(GetAssayData(young1, slot = "data"))
rm(young1)
gc()

test_ref_list <- list(count = test_data, label = test_label)

aging_for_SingleR <- GetAssayData(aging, slot = "data") ##获取标准化矩阵

aging.hesc <- SingleR(test = aging_for_SingleR, ref = test_ref_list$count, labels = test_ref_list$label)
aging@meta.data$labels <- aging.hesc$labels

## annotation
Idents(aging) <- aging$labels
aging$celltype <- aging@active.ident

plot_celltytpe <- DimPlot(aging, group.by = c("seurat_clusters", "labels"), reduction = "umap", label = T)

## 
Idents(aging) <- factor(Idents(aging), levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes", "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm", "Macrophage", "Lymphocyte", "Myoid/Interstitial", "Adult Leydig", "Sertoli", "Smooth muscle", "Endothelial"))
## 
pumap_nolabel <- DimPlot(aging, reduction = "umap") +
  theme(legend.text = element_text(size = 15),
        plot.title = element_text(size = 20)) +
  scale_color_d3("category20")


### Marker expression
genes_to_check <- c("DAZL", "DDX4", "MAGEA4", "UTF1", "FCGR3A", "KIT", 
                    "DMRT1", "DMRTB1", "STRA8", "SYCP3", "SPO11", "MLH3", 
                    "ZPBP", "ID4", "PIWIL4", "UCHL1", "TNP1", "TNP2", "PRM2", 
                    "SOX9", "WT1", "AMH", "PRND", "FATE1", "VWF", "PECAM1", 
                    "CDH5", "DLK1", "IGF1", "CYP11A1", "STAR", "NOTCH3", 
                    "ACTA2", "MYH11", "CYP26B1", "WFDC1", "CD14", "CD163", 
                    "C1QA", "C1QC", "CD8A", "CD8B", "PTPRC")
p_all_markers <- DotPlot(aging, features = genes_to_check, assay = 'RNA', group.by = 'celltype') + coord_flip()

dot_data <- p_all_markers$data

colnames(dot_data) <- c("AverageExpression_unscaled", "Precent Expressed", "Features", "celltype", "Average Expression")

####plot
p_marker_dotplot = ggplot(dot_data, aes(celltype, Features, size = ``Precent Expressed``)) +
  geom_point(shape = 21, aes(fill = ``Average Expression``), position = position_dodge(0)) +
  theme_minimal() +
  xlab(NULL) +
  ylab(NULL) +
  scale_size_continuous(range = c(1, 10)) +
  theme_bw() +
  scale_fill_gradient(low = "grey", high = "#E54924") +
  theme(legend.position = "right", legend.box = "vertical", #图例位置
        legend.margin = margin(t = 0, unit = 'cm'),
        legend.spacing = unit(0, "in"),
        axis.text.x = element_text(color = "black", size = 30, angle = 30,
                                   hjust = 1), #x轴
        axis.text.y = element_text(color = "black", size = 20), #y轴
        legend.text = element_text(size = 30, color = "black"), #图例
        legend.title = element_text(size = 30, color = "black"), #图例
        axis.title.y = element_text(vjust = 1, size = 30),
        text = element_text(size = 30)
  ) +
  labs(x = " ", y = "Features")


## 特定marker基因的Featureplot展示
p_markers_features <- FeaturePlot(aging, features = c("DDX4", "DAZL", "VIM", "PTPRC", "CD68", "CD3D"), ncol = 3)

## 细胞比例计算
## 样本间的比例变化
celltype_ratio <- prop.table(table(Idents(aging), aging$sample), margin = 2)
celltype_ratio <- as.data.frame(celltype_ratio)
colnames(celltype_ratio)[1] <- "celltype"
colnames(celltype_ratio)[2] <- "sample"
colnames(celltype_ratio)[3] <- "freq"

p_ratio <- ggplot(celltype_ratio) +
  geom_bar(aes(x = sample, y = freq, fill = celltype), stat = "identity", colour = '#222222') +
  theme_classic() +
  labs(x = 'Sample', y = 'Ratio') +
  coord_flip() +
  theme( #legend.position = "none",
    axis.text.x = element_text(color = "black", size = 20,
                               hjust = 1), #x轴
    axis.text.y = element_text(color = "black", size = 20),
    text = element_text(size = 20, face = "bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black')) +
  scale_fill_d3("category20") +                  
  scale_y_continuous(expand = c(0, 0))  

write.csv(celltype_ratio, file = paste0(outputPath, "AllCellRatio_data.csv"), row.names = FALSE)

## 
## germ cells
germ_cell <- c("SSCs", "Differentiating S'gonia", "Early primary S'cytes", "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm")
germ <- subset(aging, celltype %in% germ_cell)
Idents(germ) <- factor(Idents(germ), levels = c("SSCs", "Differentiating S'gonia", "Early primary S'cytes", "Late primary S'cytes", "Round S'tids", "Elongated S'tids", "Sperm"))
germ$celltype <- germ@active.ident

##
sperm_genes = c("TNP1", "PRM1", "PRM2")
#去除特定基因-删除特定基因
germ <- CreateSeuratObject(germ@assays$RNA@counts[-which(rownames(germ) %in% as.character(sperm_genes)),],
                           meta.data = germ@meta.data)

## 体细胞
somatic <- subset(aging,celltype %in% germ_cell,invert = TRUE)
sperm_genes = c("TNP1", "TNP2", "TNP3", "PRM1", "PRM2", "PRM3")
somatic <- CreateSeuratObject(somatic@assays$RNA@counts[-which(rownames(somatic) %in% as.character(sperm_genes)),],
                              meta.data = somatic@meta.data)
if (isSaveData) {
  # 保存注释后的aging数据
  save(somatic,germ, file = paste0(outputPath, "germ_AND_somatic.rds", seq = ""))
}

## 生殖细胞UMAP

# cellOBj <- somatic
# celltype_prefix <- "somatic_"

## 定义运行的函数
UmapCellratioFun <- function(cellOBj, celltype_prefix) {
  ### 归一化后pca降维，寻找合适的维度拐点
  # 归一化
  cellOBj <- NormalizeData(cellOBj, normalization.method = "LogNormalize", scale.factor = 10000) ####默认参数
  # 寻找高变异基因->scale归一化->跑pca降维(主成分分析)
  cellOBj <- FindVariableFeatures(cellOBj, selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA()  ###nfeatures一般选2000-5000，对结果影响较大，需要手动选择
  # 从拐点图选择合适的维度值
  pca_num <- ElbowPlot(cellOBj, ndims = 40)
  
  # 去批次
  library(harmony)
  cellOBj <- RunHarmony(cellOBj, group.by.vars = "sample")
  #降维聚类
  cellOBj <- FindNeighbors(cellOBj, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)
  cellOBj <- RunUMAP(cellOBj, reduction = "harmony", dims = 1:20, label = T) %>% RunTSNE(reduction = "harmony", dims = 1:20, label = T)
  
  ## 用整体注释的细胞类型重新定义生殖细胞
  Idents(cellOBj) <- cellOBj$celltype
  
  ## >>>> 出图
  pumap <- DimPlot(cellOBj, reduction = "umap") +
    theme(legend.text = element_text(size = 15),
          plot.title = element_text(size = 20)) +
    scale_color_d3("category20")
  
  ggsave(
    filename = paste0(outputPath, celltype_prefix, "UMAP.png", seq = ""), # 保存的文件名称。通过后缀来决定生成什么格式的图片
    width = 2000,             # 宽
    height = 1200,            # 高
    units = "px",          # 单位
    dpi = 300,              # 分辨率DPI
    plot = pumap,
    limitsize = FALSE
  )
  
  ##  这个地方数据需要调 整，注释数据集中diff很少
  
  ## 统计细胞数量
  count <- as.data.frame(table(Idents(cellOBj), cellOBj$sample))
  
  write.csv(count, file = paste0(outputPath, celltype_prefix, "cell_count.csv"), row.names = FALSE)
  
  ## 样本间比较图--副图
  
  pumap <- DimPlot(cellOBj, reduction = "umap", group.by = "sample") +
    theme(legend.text = element_text(size = 15),
          plot.title = element_text(size = 20)) +
    scale_color_d3("category20")
  
  ggsave(
    filename = paste0(outputPath, celltype_prefix, "UMAP_Sample.png", seq = ""), #保存的文件名称。通过后缀来决定生成什么格式的图片
    width = 1800,             # 宽
    height = 1200,            # 高
    units = "px",          # 单位
    dpi = 300,              # 分辨率DPI
    plot = pumap,
    limitsize = FALSE
  )
  
  
  ## 细胞比例计算
  ## 样本间的比例变化
  celltype_ratio <- prop.table(table(Idents(cellOBj), cellOBj$sample), margin = 2)
  celltype_ratio <- as.data.frame(celltype_ratio)
  colnames(celltype_ratio)[1] <- "celltype"
  colnames(celltype_ratio)[2] <- "sample"
  colnames(celltype_ratio)[3] <- "freq"
  ## celltype_ratio是个数据框，可以提取细胞比例信息
  
  ## 保存细胞比例数据
  write.csv(celltype_ratio, file = paste0(outputPath, celltype_prefix, "CellRatio_data.csv"), row.names = FALSE)
  
  p_ratio <- ggplot(celltype_ratio) +
    geom_bar(aes(x = sample, y = freq, fill = celltype), stat = "identity", colour = '#222222') +
    theme_classic() +
    labs(x = 'Sample', y = 'Ratio') +
    #coord_flip()+
    theme( #legend.position = "none",
      axis.text.x = element_text(color = "black", size = 20,
                                 hjust = 1), #x轴
      axis.text.y = element_text(color = "black", size = 20),
      text = element_text(size = 20, face = "bold"),
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.text.y = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_blank(),
      axis.line = element_line(color = 'black')) +
    scale_fill_d3("category20") +                   ## 细胞类型对应颜色
    scale_y_continuous(expand = c(0, 0)) ## 去掉与y轴之间的距离
  
  
  ## 使用scRNAtoolVis的方法进行比例计算
  p_ratio_v2 <- cellRatioPlot(object = cellOBj, sample.name = "sample", celltype.name = "celltype", col.width = 0.5, flow.alpha = 0.3) +
    scale_fill_d3("category20") +
    theme( #legend.position = "none",
      axis.text.x = element_text(color = "black", size = 20, hjust = 0.5, vjust = 0.5), #x轴
      axis.text.y = element_text(color = "black", size = 20),
      text = element_text(size = 15, face = "bold"))
  
  
  ## 差异分析
  
  # 基因差异分析
  DEG_all <- rbind.data.frame()
  
  for (item in rev(unique(cellOBj$celltype))) {
    tmp <- subset(cellOBj, idents = item)
    errCheck = tryCatch({
      tmp.markers <- FindMarkers(tmp, group.by = "sample", ident.1 = "expr", ident.2 = "ctrl", min.pct = 0.1, logfc.threshold = 0.25)
      tmp.markers$gene <- rownames(tmp.markers)
      0
    }, error = function(e) {
      # print(e)
      2
    })
    if (errCheck == 2) {
      # print("遇到错误，结束循环")
      next
    }
    tmp.markers$cluster <- item
    DEG_all <- rbind.data.frame(DEG_all, tmp.markers)
  }
  
  write.csv(DEG_all, file = paste0(outputPath, celltype_prefix, "different_expression_gene.csv", seq = ""))
  
  ###去掉线粒体基因，差异基因的火山图可视化
  mt.gene = unique(grep("^MT\\.|^MT-", DEG_all$gene, value = T))
  library(dplyr)
  DEG_all2 = subset(DEG_all, !(gene %in% mt.gene))
  p_deg_volcano <- jjVolcano(diffData = DEG_all2, legend.position = c(0.94, 0.97), topGeneN = 8) +
    scale_fill_d3("category20") +
    theme(
      axis.text.y = element_text(color = "black", size = 25),
      text = element_text(size = 25))
  
  
  ####差异基因数统计
  count <- subset(DEG_all, p_val < 0.05 & abs(avg_log2FC) > 0.25)
  count <- data.frame(table(count$cluster, count$avg_log2FC > 0.25))
  colnames(count) <- c("cluster", "log2FC", "Freq")
  count$Freq2 <- ifelse(count$log2FC == "TRUE", count$Freq, 0 - count$Freq)
  count$fill <- ifelse(count$log2FC == "TRUE", "Up", "Down")
  count$fill <- factor(count$fill, levels = c("Up", "Down"))
  deg_bar_plot <- ggplot(count, aes(x = cluster, y = Freq2, fill = fill)) +
    geom_bar(stat = 'identity', position = 'stack') +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 20, angle = 30, hjust = 1, vjust = 1),
          axis.text.y = element_text(color = "black", size = 25),
          axis.title.x = element_text(color = "black", size = 25),
          axis.title.y = element_text(color = "black", size = 25),
          plot.title = element_text(size = 25, hjust = 0.5),
          legend.text = element_text(size = 25)) +
    guides(fill = guide_legend(title = NULL)) +
    geom_text(label = count$Freq, nudge_x = 0, nudge_y = 1, size = 6) +
    xlab("Clusters") +  #x轴标签
    ylab("DEG counts") +  #y轴标签
    labs(title = "Differention Expression Gene Counts")  #设置标题
  
  ## 保存差异基因数目
  write.csv(count, file = paste0(outputPath, celltype_prefix, "DEG_gene_count.csv", seq = ""))
  
  library(clusterProfiler)
  
  #### 针对差异基因进行通路富集分析，区分上调基因和下调基因
  go.enrich = function(gene) {
    eg = bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"), OrgDb = "org.Hs.eg.db")
    ego <- enrichGO(gene = eg[, 2], OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.4,
                    qvalueCutoff = 0.2, readable = T)
    if (is.null(ego) ||
        is.null(ego@result) ||
        length(rownames(ego@result)) == 0) {
      return(NULL);
    }
    go = data.frame(ego@result)
    go$GeneRatio2 <- sapply(go$GeneRatio, function(x) eval(parse(text = x)))
    return(go)
  }
  
  ### 细胞类型的富集分析,针对基因的
  gene.enrich = function(data) {
    go <- rbind.data.frame()
    for (i in unique(data$cluster)) {
      test = subset(data, cluster == i)
      result = go.enrich(test$gene)
      if (is.null(result)) {
        next;
      }
      result$cluster = i
      go = rbind(go, result)
    }
    return(go)
  }
  
  ####去掉线粒体的基因
  DEG_all=subset(DEG_all, !(gene%in%c(grep("MT-|MT\\.",DEG_all$gene,value = T))) )
  ###定义通路点图画图函数
  GO.plot=function(data,ontology,deg){
    #点图#
    ####拆分BP和MF分开画图，美观性较好
    p1 = ggplot(subset(data, ONTOLOGY == ontology ), aes(x = cluster, y = reorder(Description, -pvalue), size = Count, color = -log10(pvalue))) +
      geom_point() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 13, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(color = "black", size = 13),
            axis.title.x = element_text(color = "black", size = 15),
            axis.title.y = element_text(color = "black", size = 15),
            plot.title = element_text(size = 25, hjust = 0.5)
      ) + scale_color_gradient(low = "lightgrey", high = "red") + xlab("Clusters") +  #x轴标签
      ylab("Pathway") +  #y轴标签
      labs(title = paste(deg," Regulate GO Terms Enrichment (",ontology,")",sep=""))   #设置标题
    return(p1)
  }
  ### 上调基因
  up = subset(DEG_all, p_val < 0.05 & avg_log2FC > 0.25)
  ### 读取蛋白编码基因的文件
  genes = read.csv(proteinPath, header = F)
  colnames(genes)[1] = "Gene"
  up = subset(up, gene %in% genes$Gene)
  ####当有报错的时候跳过报错
  tmp = try({up.go = gene.enrich(up)}, silent = TRUE)
  if (inherits(tmp, "try-error")) {
    print("An error occurred.")
  } else{
    print("No error occurred.") 
  }
  
  #####筛选通路
  
  # 使用grepl()函数进行筛选
  up.go <- up.go[!grepl(paste(keywords, collapse = "|"), up.go$Description),]
  if(dim(up.go)[1]!=0){
    up.go <- subset(up.go, Count > 2)
    write.csv(up.go, file = paste0(outputPath, celltype_prefix, "upGo.csv", seq = ""))
    upGoTopN <- up.go %>% group_by(cluster) %>% arrange(p.adjust) %>% slice_head(n = 15) %>% arrange(desc(Count))
    
  }
  
  ### 下调基因
  down = subset(DEG_all, p_val < 0.05 & avg_log2FC < -0.25)
  ####只筛选蛋白编码基因的文件
  down = subset(down, gene %in% genes$Gene)
  ####跳过报错
  tmp = try({down.go = gene.enrich(down)}, silent = TRUE)
  if (inherits(tmp, "try-error")) {
    print("An error occurred.")
  } else{
    print("No error occurred.") 
  }
  
  #####筛选通路
  # 使用grepl()函数进行筛选
  down.go <- down.go[!grepl(paste(keywords, collapse = "|"), down.go$Description),]
  if(dim(down.go)[1]!=0){
    down.go <- subset(down.go, Count > 2)
    write.csv(down.go, file = paste0(outputPath, celltype_prefix, "downGo.csv", seq = ""))
    downGoTopN <- down.go %>% arrange(p.adjust) %>% slice_head(n = 15) %>% arrange(desc(Count))
    
  }    
}

####添加一步判断，如果没有生殖细胞则跳过生殖细胞的分析，同理体细胞
if (dim(table(germ$celltype,germ$sample))[2]!=1){
  ## 运行生殖细胞
  UmapCellratioFun(germ, "germ_")
}else{
  print('There is no germ cells')
}
if (dim(table(somatic$celltype,somatic$sample))[2]!=1){
  ## 运行体细胞
  UmapCellratioFun(somatic, "somatic_")
}else{
  print('There is no somatic cells')
}


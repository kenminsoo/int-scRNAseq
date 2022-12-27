library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(DoubletFinder)

#some parameters

cmf = 'filtered_feature_bc_matrix'
project = 'nb'

######

adj.matrix <- Read10X(cmf)
srat <- CreateSeuratObject(adj.matrix,project = project)

#erase adj.matrix
adj.matrix <- NULL

#get meta data
meta <- srat@meta.data

#add percentage of mitochrodrial and ribosomal genes
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")

#makedir for figures
system('mkdir figures')
system('mkdir figures/prelim_analysis')
system('touch figures/prelim_analysis/info.txt')
system('echo "unnormalized raw data figures" > figures/prelim_analysis/info.txt')
system('mkdir figures/norm_analysis')
system('touch figures/norm_analysis/info.txt')
system('echo "normalized data figures" > figures/norm_analysis/info.txt')

#raw plots showing basically metadata features of the data
vln_metadata <- VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

ggsave('figures/prelim_analysis/vln_metadata.pdf', plot=vln_metadata, device = pdf)

scatter_mitovsUMI <- FeatureScatter(srat, feature1 = 'nCount_RNA', feature2 = 'percent.mt') + ggtitle('Percent Mitochondrial Genes vs. # UMI Counted') 
  
ggsave('figures/prelim_analysis/scatter_mitoVSumi.pdf', plot=scatter_mitovsUMI, device = pdf)

#filter out cells with bad QC metrics
#metrics
#Features > 300
#UMI > 500
#mt % < 20
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 300, 'Low_nFeature', 'Pass')
srat[['QC']] <- ifelse(srat@meta.data$nCount_RNA < 500 & srat@meta.data$QC != 'Pass', paste('Low_nCount',', Low_nCount'), srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
table(srat[['QC']])

#filtered violin plot
filt_vln_metadata <- VlnPlot(subset(srat, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

ggsave('figures/prelim_analysis/filt_vln_metadata.pdf', plot=filt_vln_metadata, device = pdf)

#normalization

norm_srat <- NormalizeData(srat)

#genes with high variation
norm_srat <- FindVariableFeatures(norm_srat, selection.method = "vst", nfeatures = 2000)

var_ft <- VariableFeaturePlot(norm_srat)

#find top 10 var genes
top10 <- head(VariableFeatures(norm_srat), 10)
var_ft_lab <- LabelPoints(var_ft, top10, repel = TRUE, xnudge = 0, ynudge = 0) + ggtitle('Variance vs. Average Expression')
ggsave('figures/norm_analysis/var_ft_lab.pdf', plot=var_ft_lab, device = pdf)

#PCA
all.genes <- rownames(norm_srat)
norm_srat <- ScaleData(norm_srat, features = all.genes)
norm_srat <- RunPCA(norm_srat, features = VariableFeatures(object = norm_srat))

DimPlot(norm_srat, reduction = "pca")

ElbowPlot(norm_srat)

norm_srat <- FindNeighbors(norm_srat, dims = 1:15)

norm_srat <- FindClusters(norm_srat, resolution = 0.5)

norm_srat <- RunUMAP(norm_srat, dims = 1:15, verbose = F)

umap_unlabled_unfilt <- DimPlot(norm_srat,label.size = 4,repel = T,label = T) + ggtitle('UMAP Visualization (Unfiltered)')
ggsave('figures/norm_analysis/umap_unlabeled_unfilt.pdf', plot=umap_unlabled_unfilt, device = pdf)

#Let's try to identify NCC,

FeaturePlot(norm_srat, features = c('CD14','LYZ','CD52','CD7'))
FeaturePlot(norm_srat, features = c("PHOX2B", "SYP", "CKIT", "PTRH2"))

#well that failed
FeaturePlot(norm_srat, features = "percent.mt")


#take only ones that passed QC
norm_srat <- subset(norm_srat, subset = QC == 'Pass')
umap_unlabeled_filt <- DimPlot(norm_srat, label.size = 4, repel = T, label = T) + ggtitle("UMAP Visualization (Filtered)")

ggsave('figures/norm_analysis/umap_unlabeled_filt.pdf', plot=umap_unlabeled_filt, device = pdf)

FeaturePlot(norm_srat, features = 'percent.mt')

#rerun UMAP
all.genes <- rownames(norm_srat)
norm_srat <- ScaleData(norm_srat, features = all.genes)
norm_srat <- RunPCA(norm_srat, features = VariableFeatures(object = norm_srat))

DimPlot(norm_srat, reduction = "pca")

ElbowPlot(norm_srat)

norm_srat <- FindNeighbors(norm_srat, dims = 1:15)

norm_srat <- FindClusters(norm_srat, resolution = 0.5)

norm_srat <- RunUMAP(norm_srat, dims = 1:15, verbose = F)

umap_unlabeled_filt_rerun <- DimPlot(norm_srat, label.size = 4, repel = T, label = T) + ggtitle("UMAP Visualization (Filtered Rerun)")

ggsave('figures/norm_analysis/umap_unlabeled_filt_rerun.pdf', plot=umap_unlabeled_filt_rerun, device = pdf)

#cell cycle
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
norm_srat <- CellCycleScoring(norm_srat, s.features = s.genes, g2m.features = g2m.genes)
FeaturePlot(norm_srat, features = c("S.Score","G2M.Score"))

#doublet cell removal
norm_srat <- RunTSNE(norm_srat)

#pK identification
sweep_res_list_nb <- paramSweep_v3(norm_srat, PCs = 1:15, sct = FALSE)
sweep_stats_nb <- summarizeSweep(sweep_res_list_nb, GT = FALSE)
bcmvn_nb <- find.pK(sweep_stats_nb)

#homotypic doublet proportion estimate
cluster_res <- norm_srat@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(cluster_res)
nExp_poi <- round(0.075*nrow(norm_srat@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))

norm_srat <- doubletFinder_v3(norm_srat, PCs = 1:15, pN = 0.25, pK = 0.1, nExp = nExp_poi_adj, reuse.pANN = FALSE)



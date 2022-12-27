library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(DoubletFinder)
library(sctransform)
library(glmGamPoi)
library(writexl)

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
system('mkdir figures/QC')
system('touch figures/QC/info.txt')
system('echo "QC1 includes mito and gene quantity. QC2 includes doublets" > figures/QC/info.txt')

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
write.table(table(srat[['QC']]), file='figures/QC/QC1.txt')

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

#Pick number of components
elbow <- ElbowPlot(norm_srat)
ggsave('figures/norm_analysis/elbow.pdf', plot=elbow, device = pdf)

print("Please look at produced elbowplot in /figures/norm_analysis")
print("Pick number of components to use")
num_components <- readline()

num_components <- as.integer(num_components)


norm_srat <- FindNeighbors(norm_srat, dims = 1:num_components)

norm_srat <- FindClusters(norm_srat, resolution = 0.5)

norm_srat <- RunUMAP(norm_srat, dims = 1:num_components, verbose = F)

umap_unlabled_unfilt <- DimPlot(norm_srat,label.size = 4,repel = T,label = T) + ggtitle('UMAP Visualization (Unfiltered)')
ggsave('figures/norm_analysis/umap_unlabeled_unfilt.pdf', plot=umap_unlabled_unfilt, device = pdf)

#plot pre-filtering clusters showing different metrics in our metadata
preQC1_fp_mt <- FeaturePlot(norm_srat, features = "percent.mt", label.size = 4,repel = T,label = T) + ggtitle("Percent Mitochondrial Genes (pre-filter)")
preQC1_fp_rib <- FeaturePlot(norm_srat, features = "percent.rb", label.size = 4,repel = T,label = T) + ggtitle("Percent Ribosomal Genes (pre-filter)")
preQC1_fp_UMI <- FeaturePlot(norm_srat, features = "nFeature_RNA", label.size = 4,repel = T,label = T) + ggtitle("Number of UMIs (pre-filter)")

system('mkdir figures/QC/QC1')

ggsave('figures/QC/QC1/preQC1_fp_mt.pdf', plot=preQC1_fp_mt, device = pdf)
ggsave('figures/QC/QC1/preQC1_fp_rib.pdf', plot=preQC1_fp_rib, device = pdf)
ggsave('figures/QC/QC1/preQC1_fp_UMI.pdf', plot=preQC1_fp_UMI, device = pdf)

#take only ones that passed QC
norm_srat <- subset(norm_srat, subset = QC == 'Pass')
umap_unlabeled_filt <- DimPlot(norm_srat, label.size = 4, repel = T, label = T) + ggtitle("UMAP Visualization (Filtered)")

ggsave('figures/norm_analysis/umap_unlabeled_filt.pdf', plot=umap_unlabeled_filt, device = pdf)

write.table(table(norm_srat[['QC']]), 'figures/QC/QC1/post_filt_QC1.txt')

#doublet cell removal
norm_srat <- RunTSNE(norm_srat)

#pK identification
sweep_res_list_nb <- paramSweep_v3(norm_srat, PCs = 1:num_components, sct = FALSE)
sweep_stats_nb <- summarizeSweep(sweep_res_list_nb, GT = FALSE)
bcmvn_nb <- find.pK(sweep_stats_nb)

#homotypic doublet proportion estimate
cluster_res <- norm_srat@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(cluster_res)
nExp_poi <- round(0.075*nrow(norm_srat@meta.data))
nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))

norm_srat <- doubletFinder_v3(norm_srat, PCs = 1:num_components, pN = 0.25, pK = 0.1, nExp = nExp_poi_adj, reuse.pANN = FALSE)

norm_srat@meta.data$is_doublet <- norm_srat@meta.data[10:10]
norm_srat@meta.data$dub_score <- norm_srat@meta.data[9:9]
QC2_doub_score <- FeaturePlot(norm_srat, features = "dub_score", label.size = 4, repel = T, label = T) + ggtitle("Doublet Score")

system('mkdir figures/QC/QC2')
ggsave('figures/QC/QC2/dub_score_feat.pdf', plot = QC2_doub_score, device = pdf)


#we can add another doubletfinder command... it does something, but did not read into it. add in future
norm_srat[['QC']] <- ifelse(norm_srat@meta.data$is_doublet == 'Doublet' & norm_srat@meta.data$QC == 'Pass', 'Doublet', norm_srat@meta.data$QC)

write.table(table(norm_srat[['QC']]), 'figures/QC/QC2/pre_filt_QC2.txt')

norm_srat <- subset(norm_srat, subset = QC == 'Pass')

write.table(table(norm_srat[['QC']]), 'figures/QC/QC2/post_filt_QC2.txt')

#cell cycle
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
norm_srat <- CellCycleScoring(norm_srat, s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
norm_srat <- RunPCA(norm_srat, features = c(s.genes, g2m.genes))
preReg_cc_dim <- DimPlot(norm_srat, reduction = 'pca', group.by="Phase") + ggtitle("Cell Cycle, Unregressed")
ggsave('figures/QC/preReg_cc_dim.pdf', plot=preReg_cc_dim, device = pdf)

#prompt user to regress out cell cycle
print("Regress normally? N will take difference between S score and G2M score to preserve differentiating cells?")
print("Y/N/M/D - M is mito only, D is don't regress")
regress_out_cc <- readline()

print("SCT or RNA normalization method?")
print("SCT/RNA")
norm_meth <- readline()

if (norm_meth == "SCT"){
if (regress_out_cc == "Y"){
norm_srat <- SCTransform(norm_srat,method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), verbose = FALSE)
#checks if the regression worked
norm_srat <- RunPCA(norm_srat, features = c(s.genes, g2m.genes))
postReg_cc_dim <- DimPlot(norm_srat) + ggtitle("Cell Cycle, Regressed/Controlled")
ggsave('figures/QC/postReg_cc_dim.pdf', plot=postReg_cc_dim, device = pdf)

} else if (regress_out_cc == "N") {
norm_srat$CC.Difference <- norm_srat$S.Score - norm_srat$G2M.Score
norm_srat <- SCTransform(norm_srat, method = "glmGamPoi", vars.to.regress = c("CC.Difference", "percent.mt"), verbose = FALSE)

norm_srat <- RunPCA(norm_srat, features = c(s.genes, g2m.genes))
postReg_cc_dim <- DimPlot(norm_srat) + ggtitle("Cell Cycle, Regressed with CC Difference and % MT")
ggsave('figures/QC/postReg_cc_dim.pdf', plot=postReg_cc_dim, device = pdf)
} else if (regress_out_cc == "M") {
norm_srat <- SCTransform(norm_srat, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
} else {
norm_srat <- SCTransform(norm_srat, method = "glmGamPoi", verbose = FALSE)
}
} else if (norm_meth == "RNA") {
if (regress_out_cc == "Y"){
norm_srat <- ScaleData(norm_srat, vars.to.regress = c("S.Score", "G2M.Score","percent.mt"), features = rownames(norm_srat))
norm_srat <- RunPCA(norm_srat, features = c(s.genes, g2m.genes))
postReg_cc_dim <- DimPlot(norm_srat) + ggtitle("Cell Cycle, Regressed/Controlled")
ggsave('figures/QC/postReg_cc_dim.pdf', plot=postReg_cc_dim, device = pdf)

} else if (regress_out_cc == "N") {
norm_srat$CC.Difference <- norm_srat$S.Score - norm_srat$G2M.Score
norm_srat <- ScaleData(norm_srat, vars.to.regress = c("CC.Difference", "percent.mt"), features = rownames(norm_srat))

norm_srat <- RunPCA(norm_srat, features = c(s.genes, g2m.genes))
postReg_cc_dim <- DimPlot(norm_srat) + ggtitle("Cell Cycle, Regressed with CC Difference and % MT")
ggsave('figures/QC/postReg_cc_dim.pdf', plot=postReg_cc_dim, device = pdf)
} else if (regress_out_cc == "M") {
norm_srat <- ScaleData(norm_srat, vars.to.regress = "percent.mt", features = rownames(norm_srat))

} else {
norm_srat <- ScaleData(norm_srat, features = rownames(norm_srat))

}
}

print("do graphs look ok?, Y/N")
readline()
#rerun the variable features
#now rerun PCA clustering etc with new normalized data
#uses 2x the number of components as originally selected if SCT used

if (norm_meth == "RNA"){
num_components <- num_components/2
norm_srat <- FindVariableFeatures(norm_srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(norm_srat), 10)
var_ft_postfilt <- VariableFeaturePlot(norm_srat)
varft_postfilt <- LabelPoints(plot = var_ft_postfilt, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave('figures/final_analysis/var_ft_postfilt.pdf', plot = varft_postfilt, device = pdf)
} else if (norm_meth == "SCT"){
top10 <- head(VariableFeatures(norm_srat), 10)
var_ft_postfilt <- VariableFeaturePlot(norm_srat, assay = "SCT")
varft_postfilt <- LabelPoints(plot = var_ft_postfilt, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave('figures/final_analysis/var_ft_postfilt.pdf', plot = varft_postfilt, device = pdf)
}

system('mkdir figures/final_analysis')

norm_srat <- RunPCA(norm_srat, verbose = F)

elbow2 <- ElbowPlot(norm_srat)
ggsave('figures/final_analysis/elbow2.pdf', plot=elbow2, device = pdf)

norm_srat <- RunUMAP(norm_srat, dims= 1:num_components*2, verbose = F)
norm_srat <- FindNeighbors(norm_srat, dims = 1:num_components*2, verbose = F)
norm_srat <- FindClusters(norm_srat, verbose = F)

write.table(table(norm_srat[[]]$seurat_clusters), 'figures/final_analysis/cluster_info.txt')

final_unlabeled_dim <- DimPlot(norm_srat, reduction = "umap"  , label = T)
ggsave('figures/final_analysis/final_unlabeled_dim.pdf', plot=final_unlabeled_dim, device = pdf)

print("Now, please look at figures/final_analysis/final_unlabeled_dim.pdf, and check if it is okay")
print("Would you like to add/remove components? (NOTE: your choice will be multiplied by a factor of two)")
print("Y/N")
reduce_compo_choice <- readline()
while (reduce_compo_choice == 'Y'){
print("Add/Remove X: (choose an integer)")
reducing_factor <- readline()
reducing_factor <- as.integer(reducing_factor)

num_components <- num_components + reducing_factor

norm_srat <- RunPCA(norm_srat, verbose = FALSE)
norm_srat <- RunUMAP(norm_srat, dims = 1:num_components*2, verbose = F)
norm_srat <- FindNeighbors(norm_srat, dims = 1:num_components*2, verbose = F)
norm_srat <- FindClusters(norm_srat, verbose = F)

write.table(table(norm_srat[[]]$seurat_clusters), 'figures/final_analysis/cluster_info.txt')

final_unlabeled_dim <- DimPlot(norm_srat, reduction = "umap"  , label = T)
ggsave('figures/final_analysis/final_unlabeled_dim.pdf', plot=final_unlabeled_dim, device = pdf)


print("Do you need to reduce again?")
print("Y/N")
reduce_compo_choice <- readline()
}

#we will now run FeaturePlots for a series of genes that may be of interest
#test what happens if invalid gene is put in
#we will now run FeaturePlots for a series of genes that may be of interest

print("Do you want to run featureplots for a gene of interest? Y/N")
inSilico_hist <- readline()

system('mkdir figures/final_analysis/gene_features')


#test what happens if invalid gene is put in
while (inSilico_hist == 'Y') {
print("Enter gene of interest")
GOI <- readline()

GOI <- toupper(GOI)

if (GOI %in% all.genes){

GOI_plot <- FeaturePlot(norm_srat, GOI) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle(GOI)

ggsave(paste('figures/final_analysis/gene_features/',GOI, '.pdf', sep = ""), plot = GOI_plot, device = pdf)

print("Continue for different gene? Y/N")
inSilico_hist <- readline()
}else{
print("Gene is invalid. Re-enter genes")
}
}

#this is an attempt to find an optimal resolution that makes sense biologically
print("This will try to help you find an optimal resoluton that makes sense for your experiment.")
print("Please choose a step size for resolution (which goes from 0.2-2")
res_step <- readline()
res_step <- as.numeric(res_step)
res <- 0.2

system('mkdir figures/final_analysis/resolution')

DefaultAssay(norm_srat) <- "RNA"
norm_srat <- NormalizeData(norm_srat)
norm_srat <- FindVariableFeatures(norm_srat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(norm_srat)
norm_srat <- ScaleData(norm_srat, features = all.genes)

while (res < 2) {
if (norm_meth == "SCT"){
DefaultAssay(norm_srat) <- "SCT"
} else {
print("skipping SCT transfer since in RNA mode")
}

system(paste('mkdir figures/final_analysis/resolution/',res, sep=""))

norm_srat <- FindClusters(norm_srat, verbose = F, resolution = res)
write.table(table(norm_srat[[]]$seurat_clusters), paste('figures/final_analysis/resolution/',res,'/groups.txt', sep="")) 
res_dim <- DimPlot(norm_srat, label = T) + ggtitle(paste("DimPlot with Resolution", res))
ggsave(paste('figures/final_analysis/resolution/',res,'/dimplot',res,'.pdf', sep = ""), plot = res_dim, device = pdf)

DefaultAssay(norm_srat) <- "RNA"
all_markers <- FindAllMarkers(norm_srat, only.pos = T)
top_10fc_markers <- as.data.frame(all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10))
top_10p_val_markers <- as.data.frame(all_markers %>% group_by(cluster) %>% slice_min(p_val_adj,n = 10))
found_markers <- as.data.frame(all_markers %>% arrange(cluster))

write_xlsx(top_10fc_markers, paste('figures/final_analysis/resolution/',res,'/',res,'_top_10fc_markers.xlsx', sep=""))
write_xlsx(top_10p_val_markers, paste('figures/final_analysis/resolution/',res,'/',res,'_top_10p_val_markers.xlsx', sep=""))
write_xlsx(found_markers, paste('figures/final_analysis/resolution/',res,'/',res,'_all_markers.xlsx', sep=""))

#heat map of top 5 log2FC genes in each cluster
heatmap_genes <- all_markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = F)
#norm_srat <- ScaleData(norm_srat, features = as.character(unique(heatmap_genes$gene)), assay = 'RNA')
norm_srat <- SetIdent(norm_srat, value = norm_srat@meta.data$seurat_clusters)
heatmap <- DoHeatmap(norm_srat, features = as.character(unique(heatmap_genes$gene)),assay = "RNA") + ggtitle(paste("Top Log2 FC Heatmap with Resolution", res)) 
ggsave(paste('figures/final_analysis/resolution/',res,'/heatmap_FC',res,'.pdf', sep=""), plot = heatmap, device = pdf)

#heat map of top 5 p val gene in each cluster
heatmap_genes_p <- all_markers %>% group_by(cluster) %>% slice_min(p_val_adj, n = 5, with_ties = F)
#norm_srat <- ScaleData(norm_srat, features = as.character(unique(heatmap_genes$gene)), assay = 'RNA')
norm_srat <- SetIdent(norm_srat, value = norm_srat@meta.data$seurat_clusters)
heatmap <- DoHeatmap(norm_srat, features = as.character(unique(heatmap_genes_p$gene)),assay = "RNA") + ggtitle(paste("Top adj. P Value Heatmap with Resolution", res))
ggsave(paste('figures/final_analysis/resolution/',res,'/heatmap_pval',res,'.pdf', sep=""), plot = heatmap, device = pdf)

#violin plot of malignant genes
vln_plasma <- VlnPlot(norm_srat, features = c('MZB1', 'BRSK1', 'AC026202.3', 'JSRP1', 'LINC00582', 'PARM1'), ncol = 3, assay = "RNA", pt.size = 1)
ggsave(paste('figures/final_analysis/resolution/',res,'/plasma_vln',res,'.pdf', sep=""), plot = vln_plasma, device = pdf)
#violin plot of top genes (2x)
heatmap_genes %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 2, with_ties = F) -> top2_FC
heatmap_genes_p %>% group_by(cluster) %>% slice_min(p_val_adj, n = 2, with_ties = F) -> top2_pval

vln_top_fc <- VlnPlot(norm_srat, features = as.character(unique(top2_FC$gene)), ncol = 5, assay = "RNA", pt.size = 0)
ggsave(paste('figures/final_analysis/resolution/',res,'/vln_top_fc',res,'.pdf', sep=""), plot = vln_top_fc, device = pdf)

vln_top_pval <- VlnPlot(norm_srat, features = as.character(unique(top2_pval$gene)), ncol = 5, assay = "RNA", pt.size = 0)
ggsave(paste('figures/final_analysis/resolution/',res,'/vln_top_pval',res,'.pdf', sep=""), plot = vln_top_pval, device = pdf)

res <- res + res_step

}

print("please look at the resolutions file and choose an optimal resolution")
opt_res <- readline()
opt_res <- as.numeric(opt_res)

if (norm_meth == "SCT"){
DefaultAssay(norm_srat) <- "SCT"
} else {
print("skipping SCT transfer since in RNA mode")
}

norm_srat <- FindClusters(norm_srat, verbose = F, resolution = opt_res)
norm_srat <- SetIdent(norm_srat, value = norm_srat@meta.data$seurat_clusters)
write.table(table(norm_srat[[]]$seurat_clusters), paste('figures/final_analysis/best_groups.txt', sep=""))
best_res_dim <- DimPlot(norm_srat, label = T) + ggtitle(paste("DimPlot with Resolution", res))
ggsave(paste('figures/final_analysis//best_res_dimplot.pdf', sep = ""), plot = best_res_dim, device = pdf)

print("Do you want to run featureplots for a gene of interest with optimal resolution? Y/N")
inSilico_hist <- readline()

#test what happens if invalid gene is put in
system('mkdir figures/final_analysis/gene_features_optres')

while (inSilico_hist == 'Y') {
print("Enter gene of interest")
GOI <- readline()

GOI <- toupper(GOI)

if (GOI %in% all.genes){

GOI_plot <- FeaturePlot(norm_srat, GOI, label = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle(GOI)

ggsave(paste('figures/final_analysis/gene_features_optres/',GOI, '.pdf', sep = ""), plot = GOI_plot, device = pdf)

print("Continue for different gene? Y/N")
inSilico_hist <- readline()
}

else{
print("Gene is invalid. Re-enter genes")
}

}

#now annotate the cells
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
sce <- as.SingleCellExperiment(DietSeurat(norm_srat))
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)

system('mkdir figures/final_analysis/annotated')
write.table(table(hpca.main$pruned.labels), 'figures/final_analysis/annotated/main_labels_hpca.txt')
write.table(table(hpca.fine$pruned.labels), 'figures/final_analysis/annotated/fine_labels_hpca.txt')

norm_srat@meta.data$hpca.main   <- hpca.main$pruned.labels
norm_srat@meta.data$hpca.fine   <- hpca.fine$pruned.labels
norm_srat <- SetIdent(norm_srat, value = "hpca.fine")
fine_labeled_dimplot <- DimPlot(norm_srat, label = T , repel = T, label.size = 3) + NoLegend()
ggsave('figures/final_analysis/annotated/fine_labeled_dimplot.pdf', plot = fine_labeled_dimplot, device = pdf)
norm_srat <- SetIdent(norm_srat, value = "hpca.main")
main_labeled_dimplot <- DimPlot(norm_srat, label = T , repel = T, label.size = 3) + NoLegend()
ggsave('figures/final_analysis/annotated/main_labeled_dimplot.pdf', plot = main_labeled_dimplot, device = pdf)

print("done. continue with more downstream analysis after looking at this data")


#allow to run featureplots with labels
print("Do you want to run featureplots for a gene of interest with cell labels? Y/N")
inSilico_hist <- readline()

#test what happens if invalid gene is put in
system('mkdir figures/final_analysis/gene_features_cell_label')

while (inSilico_hist == 'Y') {
print("Enter gene of interest")
GOI <- readline()

GOI <- toupper(GOI)

if (GOI %in% all.genes){

GOI_plot <- FeaturePlot(norm_srat, GOI, label = T) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle(GOI)

ggsave(paste('figures/final_analysis/gene_features_cell_label/',GOI, '.pdf', sep = ""), plot = GOI_plot, device = pdf)

print("Continue for different gene? Y/N")
inSilico_hist <- readline()
}

else{
print("Gene is invalid. Re-enter genes")
}

}

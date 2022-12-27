#uses packagess from seurat processing. 
#must run seurat processing.R before running this file.
#based upon your selection, this script will look at a certain type of cell found in analysis
#this will find:
#Marker genes for each sub-celltpe cluster
#Malignancy scores for each
#RNA-velocity for specific celltype
#RNA-velocity, if chosen, for all cells found
#Allows researcher to zoom in on specific areas of the original clustering. 


system('mkdir figures/subset')
print("Choose what cell type you'd like to subset on")
act_subset <- readline()

srat_sub <- subset(norm_srat, idents = act_subset)

#Recluster 
srat_sub <- RunPCA(srat_sub, verbose = F)

elbow_sub <- ElbowPlot(srat_sub, ndims = 40)
ggsave('figures/subset/elbow_sub.pdf', plot=elbow_sub, device = pdf)


srat_sub <- RunUMAP(srat_sub, dims = 1:num_components, verbose = F)
srat_sub <- FindNeighbors(srat_sub, dims = 1:num_components, verbose = F)
srat_sub <- FindClusters(srat_sub, verbose = F, resolution = opt_res)
dimplot_sub <- DimPlot(srat_sub, label = T) + ggtitle(paste(act_subset, "Subset"))
ggsave('figures/subset/dimplot_sub.pdf', plot=dimplot_sub, device = pdf)

#Find variable genes in clusters with RNA assay
DefaultAssay(srat_sub) <- "RNA"
all_markers_sub <- FindAllMarkers(srat_sub, only.pos = T)
#It will probably recalculate here

top_10fc_markers_sub <- as.data.frame(all_markers_sub %>% group_by(cluster) %>% slice_max(avg_log2FC,n = 10))
top_10p_val_markers_sub <- as.data.frame(all_markers_sub %>% group_by(cluster) %>% slice_min(p_val_adj,n = 10))
found_markers_sub <- as.data.frame(all_markers_sub %>% arrange(cluster))

write_xlsx(top_10fc_markers_sub, "figures/subset/top10_fc_markers_subset.xlsx")
write_xlsx(top_10p_val_markers_sub, "figures/subset/top10_pval_markers_subset.xlsx")
write_xlsx(found_markers_sub, "figures/subset/all_markers.xlsx")

heatmap_genes_sub <- all_markers_sub %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5, with_ties = F)

heatmap_genes_sub_p <- all_markers_sub %>% group_by(cluster) %>% slice_min(p_val_adj, n = 5, with_ties = F)

srat_sub <- SetIdent(srat_sub, value = srat_sub@meta.data$seurat_clusters)
heatmap_sub <- DoHeatmap(srat_sub, features = as.character(unique(heatmap_genes_sub$gene)), assay = "RNA") + ggtitle(paste("Top Log2 FC for", act_subset))
ggsave('figures/subset/heatmap_LOG2_sub.pdf', plot = heatmap_sub, width = 24, height = 24, device = pdf)

heatmap_p_sub <- DoHeatmap(srat_sub, features = as.character(unique(heatmap_genes_sub_p$gene)), assay = "RNA") + ggtitle(paste("Top p value for", act_subset))
ggsave('figures/subset/heatmap_p_val_sub.pdf', plot = heatmap_p_sub, width = 24, height = 24, device = pdf)

heatmap_genes_sub %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 2, with_ties = F) -> top2_FC_sub

heatmap_genes_sub_p %>% group_by(cluster) %>% slice_min(p_val_adj, n = 2, with_ties = F) -> top2_pval_sub

vln_top_fc_sub <- VlnPlot(srat_sub, features = as.character(unique(top2_FC_sub$gene)), ncol = 5, assay = "RNA", pt.size = 0)
ggsave('figures/subset/vlnPlot_sub_FC.pdf', plot = vln_top_fc_sub, width = 24, height = 24, device = pdf)

vln_top_pval_sub <- VlnPlot(srat_sub, features = as.character(unique(top2_pval_sub$gene)), ncol = 5 ,assay = "RNA", pt.size = 0)

ggsave('figures/subset/vlnPlot_sub_Pval.pdf', plot = vln_top_pval_sub, width = 24, height = 24, device = pdf)

#now let's start RNA Velocity

#documenting RNA velocity workflow

##non python commands (commandline)
###Move to a shell script that runs, then runs this workflow.

#kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno Homo_sapiens.GRCh38.cdna.all.fa.gz /home/kenmn/bioinfo_tools/refdata-gex-GRCh38-2020-A/genes/genes.gtf
#kb count -i index.idx -g t2g.txt -x 10xv2 --workflow lamanno --loom -c1 cdna_t2c.txt -c2 intron_t2c.txt /home/kenmn/nfs_scratch/nb12799275/fastq/SRR10156295_S1_L001_R1_001.fastq /home/kenmn/nfs_scratch/nb12799275/fastq/SRR10156295_S1_L001_R2_001.fastq

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample = anndata.read_loom("velocity/possorted_genome_bam_5TKBT.loom")

sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters.csv")

#fix names, in format: possorted_genome_bam_5TKBT:AAACCTGAGACCCACCx
sample_obs['x'] = 'possorted_genome_bam_5TKBT:' + sample_obs['x'].astype(str)
sample_obs['x'] = sample_obs['x'].str[:-2]
sample_obs['x'] = sample_obs['x'].astype(str) +'x'

umap_cord['Unnamed: 0'] = 'possorted_genome_bam_5TKBT:' + umap_cord['Unnamed: 0'].astype(str)
umap_cord['Unnamed: 0'] = umap_cord['Unnamed: 0'].str[:-2]
umap_cord['Unnamed: 0'] = umap_cord['Unnamed: 0'].astype(str) + 'x'

cell_clusters['Unnamed: 0'] = umap_cord['Unnamed: 0']

#filter
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]

umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'CellID'})

sample_index = pd.DataFrame(sample.obs.index)
umap_ordered = sample_index.merge(umap_cord, on = "CellID")
umap_ordered = umap_ordered.iloc[:,1:]

sample.obsm['X_umap'] = umap_ordered.values

cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'CellID'})
cell_clusters_reordered = sample_index.merge(cell_clusters, on = "CellID")
cell_clusters_reordered = cell_clusters_reordered.iloc[:,1:]

sample.obs['ClusterName'] = cell_clusters_reordered.values
sample.obs['ClusterName'] = pd.Categorical(sample.obs.ClusterName)

#export from seurat
ident_colors = ['#F8766D','#E18A00','#BE9C00','#8CAB00','#24B700','#00BE70','#00C1AB','#00BBDA','#00ACFC','#8B93FF','#D575FE','#F962DD','#FF65AC', '#3E1C14']

sample.uns['ClusterName_colors'] = ident_colors

scv.pl.proportions(sample, save = 'prop.pdf')

scv.pp.filter_and_normalize(sample)
scv.pp.moments(sample)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
scv.pl.velocity_embedding(sample, basis = 'umap', title = "Overall RNA Velocity of NB12799275", save = 'velo_rna.pdf', color = 'ClusterName')

scv.tl.velocity_clusters(sample)

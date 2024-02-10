#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:21:52 2024

@author: liza
"""

import numpy as np
from scipy.sparse import csc_matrix
import scanpy as sc
import pandas as pd
import decoupler as dc

# Plotting options, change to your liking
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(5, 6))

#%% FROM SEURAT FILES GC

# Read the .mtx file along with the barcodes and gene names
adata_GC = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/GC.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GCgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GCbarcodes.tsv", header=None, sep='\t')[0]

adata_GC = adata_GC.transpose()

# Check dimensions and assign names
if adata_GC.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_GC.var_names = genes
    adata_GC.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_GC.var_names = genes
adata_GC.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_GC.var_names_make_unique(join="-")
adata_GC.obs_names_make_unique(join="-")

#%% FROM SEURAT FILES GP

# Read the .mtx file along with the barcodes and gene names
adata_GP = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/GP.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GPgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GPbarcodes.tsv", header=None, sep='\t')[0]

adata_GP = adata_GP.transpose()

# Check dimensions and assign names
if adata_GP.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_GP.var_names = genes
    adata_GP.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_GP.var_names = genes
adata_GP.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_GP.var_names_make_unique(join="-")
adata_GP.obs_names_make_unique(join="-")
#%%FROM SEURAT FILES WP

# Read the .mtx file along with the barcodes and gene names
adata_WP = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/WP.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WPgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WPbarcodes.tsv", header=None, sep='\t')[0]

adata_WP = adata_WP.transpose()

# Check dimensions and assign names
if adata_WP.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_WP.var_names = genes
    adata_WP.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_WP.var_names = genes
adata_WP.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_WP.var_names_make_unique(join="-")
adata_WP.obs_names_make_unique(join="-")


#%% FROM SEURAT FILES WC

# Read the .mtx file along with the barcodes and gene names
adata_WC = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/WC.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WCgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WCbarcodes.tsv", header=None, sep='\t')[0]

adata_WC = adata_WC.transpose()

# Check dimensions and assign names
if adata_WC.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_WC.var_names = genes
    adata_WC.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_WC.var_names = genes
adata_WC.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_WC.var_names_make_unique(join="-")
adata_WC.obs_names_make_unique(join="-")

#%% Put all the data together

sc.pp.filter_cells(adata_GC, min_genes=200)
sc.pp.filter_genes(adata_GC, min_cells=10)

sc.pp.filter_cells(adata_GP, min_genes=200)
sc.pp.filter_genes(adata_GP, min_cells=10)

sc.pp.filter_cells(adata_WC, min_genes=200)
sc.pp.filter_genes(adata_WC, min_cells=10)

sc.pp.filter_cells(adata_WP, min_genes=200)
sc.pp.filter_genes(adata_WP, min_cells=10)

adata_GC.obs['sample'] = 'dblGATA_Control'
adata_GP.obs['sample'] = 'dblGATA_PyMT'
adata_WC.obs['sample'] = 'WT_Control'
adata_WP.obs['sample'] = 'WT_PyMT'

adata = sc.concat([adata_GC, adata_GP, adata_WC, adata_WP])

# Subset the adata object for dblGATA samples
adata_dblGATA = adata[adata.obs['sample'].isin(['dblGATA_Control', 'dblGATA_PyMT'])]


# Subset the adata object for WT samples
adata_WT = adata[adata.obs['sample'].isin(['WT_Control', 'WT_PyMT'])]

#%% Filtering the data for dead cells

# Identifying mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('mt')  # Adjust the prefix if needed

# Identifying ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))  # Adjust these prefixes as needed

# Calculating quality control metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Plots for QC
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# Filtering cells based on gene count thresholds
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)
adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]

# Filtering cells based on mitochondrial gene percentage
adata = adata[adata.obs.pct_counts_mt < 13, :]

# Removing mitochondrial and ribosomal genes from the dataset
adata = adata[:, ~adata.var['mt'] & ~adata.var['ribo']]

#%% Check for ribosomal genes and mitochondrial genes

ribo_genes = adata.var_names.str.startswith(('Rps', 'Rpl'))
if any(ribo_genes):
    print("Ribosomal genes starting with 'Rps' or 'Rpl' found:")
    print(adata.var_names[ribo_genes])
else:
    print("No ribosomal genes starting with 'Rps' or 'Rpl' found.")


mt_genes = adata.var_names.str.startswith('mt-')
if any(mt_genes):
    print("Mitochondrial genes starting with 'mt-' found:")
    print(adata.var_names[mt_genes])
else:
    print("No mitochondrial genes starting with 'mt-' found.")
    
#%% Check if the values are raw counts and MAKE A COPY of RAW counts

from scipy.sparse import issparse

# Check if the data is stored as a sparse matrix
if issparse(adata.X):
    # Convert to a dense format for viewing
    dense_X = adata.X.toarray()
    print("Adata is sparse-matrix")
else:
    # If it's already a dense format, just assign it
    dense_X = adata.X

print(dense_X[:5, :5])

adata.raw = adata

#%% Data normalization

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%% Scaling data

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca(adata, color = 'sample')

sc.pl.pca_variance_ratio(adata, log=True)

#%% Plotting UMAP

print(adata.X[:5, :5])

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

sc.tl.umap(adata)

sc.pl.umap(adata)

sc.tl.leiden(adata, resolution = 0.3)
sc.pl.umap(adata, color=['leiden'], frameon = False, legend_loc = 'on data')
sc.pl.umap(adata, color=['sample'], frameon = False)

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

#%% Subplots for 4 samples

import matplotlib.pyplot as plt


# Create a figure with 4 subplots (2x2 layout)
fig, axs = plt.subplots(2, 2, figsize=(12, 12))

# Define samples and subplot titles
samples = ['WT_Control', 'dblGATA_Control','WT_PyMT', 'dblGATA_PyMT']
subplot_titles = ['WT Control', 'dblGATA Control', 'WT PyMT', 'dblGATA PyMT']

# Plot UMAP for each sample in a different subplot
for i, sample in enumerate(samples):
    ax = axs[i // 2, i % 2]
    sc.pl.umap(adata[adata.obs['sample'] == sample], ax=ax, color='leiden', 
               title=subplot_titles[i], show=False, frameon=False, legend_loc = 'on data')

# Adjust the layout and display the plot
plt.tight_layout()
plt.show()

#%% Save TOP ranked genes

import matplotlib.pyplot as plt
import seaborn as sns


top_n_genes = 25

TOP_genes_cluster = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(top_n_genes)

TOP_genes_cluster.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/TOPMarkers.csv', sep=',', encoding='utf-8', header='true', index=False)

#%% UMAP colors for genes

import matplotlib.colors as mcolors

# Create a custom color map from grey to violet
colors = ['#d7d4d9', '#380357']
cmap = mcolors.LinearSegmentedColormap.from_list('grey_to_purple', colors, N=10)

#sc.pl.umap(adata, color=['leiden'], frameon = False, legend_loc = 'on data')

#%% B-cells - 0

sc.pl.umap(adata, color = ['Cd79a', 'Ms4a1', 'Cd79b', 'Ighd', 'Cd19'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Neutrophils - 1

sc.pl.umap(adata, color = ['G0s2', 'Clec4d', 'Cd14'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Alveolar macrophage - 4

sc.pl.umap(adata, color=['Ear2', 'Csf1r', 'Cx3cr1'], frameon=False, legend_loc='on data', cmap=cmap)

#%% CD8+ T-cell - 6

sc.pl.umap(adata, color = ['Cd3d','Cd8a', 'Cd8b1', 'Il7r', 'Ccr7', 'Klrc1', 'Eomes', 'Cxcr3'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Plasmacytoid dendritic cell - 13

sc.pl.umap(adata, color = ['Ms4a6c', 'Plac8', 'Bst2', 'Irf7', 'Siglech'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Dendritic -  5

sc.pl.umap(adata, color = ['Naaa', 'Irf8', 'Itgae'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% NK Cell - 3

sc.pl.umap(adata, color = ['Nkg7', 'Klra8', 'Klra4', 'Klrb1c'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Interstitial macrophage - 8

sc.pl.umap(adata, color = ['C1qc', 'C1qa', 'Pf4', 'Mertk', 'Vcan'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% CD4+ T-cells - 2

sc.pl.umap(adata, color = ['Cd4', 'Tnfrsf4', 'Il7r', 'Ccr7', 'Icos'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%%  T-reg - 

sc.pl.umap(adata, color = ['Foxp3', 'Ctla4', 'Il2ra', 'Tnfrsf18'], frameon=False, legend_loc = 'on data', cmap=cmap)

#%% Th1

sc.pl.umap(adata, color = ['Tbx21', 'Ifng'], frameon=False, legend_loc = 'on data', cmap=cmap)


#%% ANNOTATION CELL CLUSTERS

cell_type = {
    '0': "B",
    '1': "Neut",
    '2': "CD4+",
    '3': "NK",
    '4': "Alv Macro",
    '5': "Dend",
    '6': "CD8+",
    '7': "T",
    '8': "Inter Macro",
    '9': "Non_imune",
    '10': "Th2",
    '11': "Non_imune",
    '12': "Non_imune",
    '13': "Pcdc"
    }

adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type)

sc.set_figure_params(figsize=(6, 6))


sc.pl.umap(adata, color='cell_type', legend_loc = 'on data', title='Cell types', frameon=False)

#%% T-cells cluster

adata_T = adata[(adata.obs['cell_type'] == 'CD8+') | (adata.obs['cell_type'] == 'CD4+') | (adata.obs['cell_type'] == 'T')| (adata.obs['cell_type'] == 'Th2')].copy()


#%% T cluster subcluster

sc.pp.neighbors(adata_T, n_neighbors=4, n_pcs=20)

sc.tl.leiden(adata_T, resolution = 0.2)

sc.tl.umap(adata_T)

sc.pl.umap(adata_T)

sc.tl.rank_genes_groups(adata_T, 'leiden', method='wilcoxon')

sc.pl.umap(adata_T, color=['leiden'], frameon = False, legend_loc = 'on data')
sc.pl.umap(adata_T, color=['sample'], frameon = False)

#%% T cells marker genes

marker_genes = pd.DataFrame(adata_T.uns['rank_genes_groups']['names'])

marker_genes.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/TOPMarkers_T.csv', sep=',', encoding='utf-8', header='true', index=False)

#%% Shai annotation, Dot-plots for clusters

samples_to_plot = ['WT_PyMT', 'dblGATA_PyMT']

sc.set_figure_params(figsize=(10, 20))


# Create a dictionary to map cluster labels to T cell states


# List of marker genes for different T cell states
marker_genes = {
    'Naive CD4+': ['Cd4', 'Tcf7', 'Ccr7', 'Lef1'],
    'Th1': ['Ifng', 'Tbx21', 'Tnf', 'Il2ra'],
    'Tr1': ['Il10', 'Ifng', 'Lag3', 'Maf'],
    'Th2': ['Gata3', 'Cxcr3'],
    'T-regs': ['Foxp3', 'Ctla4', 'Il2ra', 'Il10'],
    'Naive CD8+': ['Cd8a', 'Tcf7', 'Ccr7', 'Lef1'],
    'Effector CD8+': ['Gzma', 'Gzmk', 'Gzmb', 'Prf1', 'Ifng', 'Cd27'],
    'Memory CD8+': ['Il7r', 'Klrg1', 'Bcl2'],
    'Exhausted CD8+': ['Pdcd1', 'Lag3', 'Ctla4']
}

# Iterate through the samples and create a dot plot for each
for sample in samples_to_plot:
    # Subset the data for the current sample
    adata_sample = adata_T[adata_T.obs['sample'] == sample].copy()
    
    # Perform dot plot
    sc.pl.dotplot(adata_sample, marker_genes, groupby='leiden', swap_axes=True, vmax=5)
    
    # Annotate clusters with their corresponding T cell states
    cluster_to_state = {
        '0': 'Naive CD4+ T',
        '1': 'Th1',
        '2': 'Tr1',
        '3': 'Th2',
        '4': 'Tregs',
        '5': 'Naive CD8+ T',
        '6': 'Effector CD8+ T',
        '7': 'Memory CD8+ T',
        '8': 'Exhausted CD8+ T'
    }   
    # Show the plot
    plt.show()


#%%Look on the marker genes

sc.set_figure_params(figsize=(5, 6))


marker_genes = {
    'Naive CD4+': ['Cd4', 'Tcf7', 'Ccr7', 'Lef1'],
    'Th1': ['Ifng', 'Tbx21', 'Tnf', 'Il2ra'],
    'Tr1': ['Il10', 'Ifng', 'Lag3', 'Maf'],
    'Th2': ['Gata3', 'Cxcr3'],
    'T-regs': ['Foxp3', 'Ctla4', 'Il2ra', 'Il10'],
    'Naive CD8+': ['Cd8a', 'Tcf7', 'Ccr7', 'Lef1'],
    'Effector CD8+': ['Gzma', 'Gzmk', 'Gzmb', 'Prf1', 'Ifng', 'Cd27'],
    'Memory CD8+': ['Il7r', 'Klrg1', 'Bcl2'],
    'Exhausted CD8+': ['Pdcd1', 'Lag3', 'Ctla4']
}

sample_names = ['WT_PyMT', 'dblGATA_PyMT']

sc.pl.umap(adata_T[adata_T.obs['sample'] == 'WT_PyMT'], color = ['Stat2'], frameon=False, legend_loc = 'on data', cmap=cmap)


#%% ANNOTATION CELL CLUSTERS in T-cells sub

cell_type_T = {
    '0': "Exh_CD4&Treg",
    '1': "Effector_CD8",
    '2': "Mix",
    '3': "Naive_CD4&CD8",
    '4': "Naive_CD4&CD8",
    '5': "Th2",
    '6': "Exh_CD4",
    '7': "Exh_CD4&CD8",
    '8': "B_cont"
    }

adata_T.obs['cell_type'] = adata_T.obs['leiden'].map(cell_type_T)

sc.set_figure_params(figsize=(8, 8))

sc.pl.umap(adata_T[adata_T.obs['sample'] == 'WT_PyMT'], color='cell_type', legend_loc = 'on data', title='WT T-cells Cell types', frameon=False)
sc.pl.umap(adata_T[adata_T.obs['sample'] == 'dblGATA_PyMT'], color='cell_type', legend_loc = 'on data', title='dblGATA T-cells Cell types', frameon=False)

adata_T = adata_T[adata_T.obs['cell_type'] != 'B_cont'].copy()

#%% Subplots for 4 samples

import matplotlib.pyplot as plt


# Create a figure with 4 subplots (2x2 layout)
fig, axs = plt.subplots(2, 2, figsize=(12, 12))

# Define samples and subplot titles
samples = ['WT_Control', 'dblGATA_Control','WT_PyMT', 'dblGATA_PyMT']
subplot_titles = ['WT Control', 'dblGATA Control', 'WT PyMT', 'dblGATA PyMT']

# Plot UMAP for each sample in a different subplot
for i, sample in enumerate(samples):
    ax = axs[i // 2, i % 2]
    sc.pl.umap(adata_T[adata_T.obs['sample'] == sample], ax=ax, color='cell_type', 
               title=subplot_titles[i], show=False, frameon=False, legend_loc = 'on data')

# Adjust the layout and display the plot
plt.tight_layout()
plt.show()


#%% CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional targets compiled from 12 different resources. 

net = dc.get_collectri(organism='mouse', split_complexes=False)

#%% To infer TF enrichment scores we will run the univariate linear model (ulm) method. 
#For each cell in our dataset (adata) and each TF in our network (net), it fits a linear model that predicts the observed gene expression 
#based solely on the TF’s TF-Gene interaction weights. Once fitted, the obtained t-value of the slope is the score. If it is positive, we interpret that the TF is active and if it is negative we interpret that it is inactive.

dc.run_ulm(
    mat=adata_T,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True
)

#%% activity infered for TF across cells WT

acts_WT = dc.get_acts(adata_T[adata_T.obs['sample'] == 'WT_PyMT'], obsm_key='ulm_estimate')
acts = dc.get_acts(adata_T, obsm_key='ulm_estimate')

#%% activity infered for TF across cells dblGATA

acts_dblGATA = dc.get_acts(adata_T[adata_T.obs['sample'] == 'dblGATA_PyMT'], obsm_key='ulm_estimate')

#%%

sc.pl.umap(acts_WT, color=['Stat2', 'cell_type'], cmap='RdBu_r', vcenter=0)
sc.pl.umap(acts_dblGATA, color=['Stat2', 'cell_type'], cmap='RdBu_r', vcenter=0)

sc.pl.violin(acts_WT, keys=['Stat2'], groupby='cell_type', rotation=90)
sc.pl.violin(acts_dblGATA, keys=['Stat2'], groupby='cell_type', rotation=90)

sc.pl.umap(acts, color=['Ddit3', 'cell_type'], cmap='RdBu_r', vcenter=0)
sc.pl.violin(acts, keys=['Ddit3'], groupby='cell_type', rotation=90)


#Here we observe the activity infered for PAX5 across cells, which it is particulary active in B cells. 
#Interestingly, PAX5 is a known TF crucial for B cell identity and function. 
#The inference of activities from “foot-prints” of target genes is more informative than just looking at the molecular readouts of a given TF, as an example here is the gene expression of PAX5, which is not very informative by itself since it is just expressed in few cells:


#%% top TF per cell type

df_WT = dc.rank_sources_groups(acts_WT, groupby='cell_type', reference='rest', method='t-test_overestim_var')
df_dblGATA = dc.rank_sources_groups(acts_dblGATA, groupby='cell_type', reference='rest', method='t-test_overestim_var')

df = dc.rank_sources_groups(acts, groupby='cell_type', reference='rest', method='t-test_overestim_var')

#%% extract the top 10 markers per cell type:

n_markers = 5
source_markers_WT = df_WT.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
source_markers_dblGATA = df_dblGATA.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

#%% plot the obtained markers:

sc.pl.matrixplot(acts_WT, source_markers_WT, 'cell_type', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

sc.pl.matrixplot(acts_dblGATA, source_markers_dblGATA, 'cell_type', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

sc.pl.matrixplot(acts, source_markers, 'sample', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

#%% individual TFs by plotting their distributions

sc.pl.violin(acts_dblGATA, keys=['Irf7'], groupby='cell_type', rotation=90)

#%% Helicopter view of TFs per cell_type

# Filter rows where 'pvals' and 'pvals_adj' are < 0.05
filtered_df_WT = df_WT[(df_WT['pvals'] < 0.05) & (df_WT['pvals_adj'] < 0.05)]
filtered_df_dblGATA = df_dblGATA[(df_dblGATA['pvals'] < 0.05) & (df_dblGATA['pvals_adj'] < 0.05)]

# Combine the filtered dataframes and create a new 'dataset' column
filtered_df_WT['dataset'] = 'WT'
filtered_df_dblGATA['dataset'] = 'dblGATA'
combined_df = pd.concat([filtered_df_WT, filtered_df_dblGATA])

# Set the same y-axis limits for both plots
y_limit = max(combined_df['statistic'].max(), abs(combined_df['meanchange'].min()))

# Create a figure and plot violins for both datasets with the violins on the x-axis
plt.figure(figsize=(10, 6))
sns.violinplot(x='statistic', y='dataset', hue='group', data=combined_df, inner='quart', palette="Set3", cut=0, bw=0.2)
plt.xlim(-y_limit, y_limit)  # Set x-axis limits

# Add a single legend
plt.legend(loc='upper right', title='Group')

# Set labels and title
plt.xlabel('TFs footprint score')
plt.ylabel('Sample')
plt.title('Rankings by score (pvals & pvals_adj < 0.05)')

# Show the plot
plt.tight_layout()
plt.show()

#%% Functional enrichment of biological terms (The Molecular Signatures Database (MSigDB) )

msigdb = dc.get_resource('MSigDB')
#%% Convert to mouse

mouse_msigdb = dc.translate_net(msigdb, target_organism = 'mouse', unique_by = ('geneset', 'genesymbol'))

#%% # Filter by msigdb['collection'].unique()

msigdb_c = msigdb[msigdb['collection']=='go_molecular_function']

# Remove duplicated entries
msigdb_c = msigdb_c[~msigdb_c.duplicated(['geneset', 'genesymbol'])]

#%% Read database with mice genes and change human to mice

genes = pd.read_csv("/home/liza/Documents/PhD/human_mouse_1to1_orthologs.csv", sep=',')
genes = genes.rename(columns={"human": "genesymbol"})

merged_df = pd.merge(msigdb_c, genes, on='genesymbol', how='left')

# Select the desired columns
final_df = merged_df[['mouse', 'collection', 'geneset']]

# Rename the 'mouse' column to 'genesymbol'
final_df.rename(columns={'mouse': 'genesymbol'}, inplace=True)

# Remove duplicated entries
final_df = final_df[~final_df.duplicated(['geneset', 'genesymbol'])]

#%% Run Enrichment with Over Representation Analysis (ORA)

dc.run_ora(
    mat=adata_T,
    net=final_df,
    source='geneset',
    target='genesymbol',
    verbose=True
)

#%% ORA
acts = dc.get_acts(adata_T, obsm_key='ora_estimate')

#%% ORA

# We need to remove inf and set them to the maximum value observed
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e


#df = dc.rank_sources_groups(acts, groupby='cell_type', reference='rest', method='t-test_overestim_var')

def ora_analysis_for_cell_type_df(adata, cell_type, sample1, sample2):
    # Subset data for the current cell type
    adata_subset = adata[(adata.obs['cell_type'] == cell_type) & (adata.obs['sample'].isin([sample1, sample2]))]

    # Perform the differential expression analysis
    df = dc.rank_sources_groups(adata_subset, groupby='sample', reference=sample1, method='t-test_overestim_var')

    # Filter by pvalue and padj
    filtered_df = df[(df['pvals'] <= 0.05) & (df['pvals_adj'] <= 0.05)]
    return filtered_df


df_CD8_effector = ora_analysis_for_cell_type_df(acts, 'Effector_CD8', 'dblGATA_PyMT', 'WT_PyMT')
#%% ORA Naive in WT pver dblGATA

df_Naive_CD4_CD8 = ora_analysis_for_cell_type_df(acts, 'Naive_CD4&CD8', 'dblGATA_PyMT', 'WT_PyMT')

#%% ORA Exh_CD4&CD8&Treg

df_Exh_CD4_Treg = ora_analysis_for_cell_type_df(acts, 'Exh_CD4&Treg', 'dblGATA_PyMT', 'WT_PyMT')
df_Exh_CD4_CD8 = ora_analysis_for_cell_type_df(acts, 'Exh_CD4&CD8', 'dblGATA_PyMT', 'WT_PyMT')


#%% Barplot for ORA

# Set the aesthetic style of the plots
sns.set(style="whitegrid")


# Iterate through the groups and create a plot for each
for i, (group_name, group_data) in enumerate(df_CD8_effector.groupby('group')):
    # Separate positive and negative statistics
    positive_stats = group_data[group_data['statistic'] > 0].sort_values(by='statistic', ascending=False)
    negative_stats = group_data[group_data['statistic'] < 0].sort_values(by='statistic', ascending=True)
    
    # Concatenate top 10 positive and negative statistics
    top_positive = positive_stats.head(10)
    top_negative = negative_stats.head(10)
    top_stats = pd.concat([top_negative, top_positive], axis=0)
    
    # Create a new figure for each group
    plt.figure(figsize=(20, 8))
    
    # Create a horizontal bar plot
    ax = sns.barplot(
        x='statistic', 
        y='names', 
        data=top_stats, 
        orient='h', 
        color='lightblue'  # Default color, will be overwritten
    )
    
    # Color bars based on the sign of the 'statistic' values
    for bar, statistic in zip(ax.patches, top_stats['statistic']):
        if statistic < 0:
            bar.set_facecolor('#d7d4d9')
        else:
            bar.set_facecolor('#380357')
    
    # Set title for the plot and increase title size
    ax.set_title("<- dblGATA | ORA CD8+ effector cells | WT ->", fontsize=22)
    # Set xlabel and remove ylabel
    ax.set_xlabel("Enrichment Score", fontsize=18)
    ax.set_ylabel('')  # Remove y-axis label

    # Increase label size
    ax.xaxis.label.set_size(16)
    ax.yaxis.label.set_size(16)

    # Increase tick size
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

    # Optional: If you want to display the actual statistic value on the bars
    for index, value in enumerate(top_stats['statistic']):
        ax.text(value, index, f'{value:.2f}', color='black', va='center', fontsize=14)

    # Adjust layout
    plt.tight_layout()
    
    # Show the plot
    plt.show()

#%% Calculate DE genes

def de_analysis_for_cell_type_df(adata, cell_type, sample1, sample2, remove_ribo_mito=True):
    # Subset data for the current cell type
    adata_subset = adata[(adata.obs['cell_type'] == cell_type) & (adata.obs['sample'].isin([sample1, sample2]))]

    # Perform the differential expression analysis
    sc.tl.rank_genes_groups(adata_subset, groupby='sample', reference=sample1)

    # Extract DE genes results
    de_genes = adata_subset.uns['rank_genes_groups']

    # Convert to DataFrame
    results_df = pd.DataFrame({
        'genes': de_genes['names'][sample2],
        'logfoldchanges': de_genes['logfoldchanges'][sample2],
        'pvals': de_genes['pvals'][sample2],
        'pvals_adj': de_genes['pvals_adj'][sample2]
    })
    
    
    # Filter out rows where logfoldchanges is NaN or infinite
    results_df = results_df.replace([np.inf, -np.inf], np.nan).dropna(subset=['logfoldchanges'])

    # Filter by pvalue and padj
    filtered_df = results_df[(results_df['pvals'] <= 0.05) & (results_df['pvals_adj'] <= 0.05)]
    return filtered_df
    print(f"Saved DE genes to {filtered_df}")

df_CD4 = de_analysis_for_cell_type_df(adata, 'CD4+', 'dblGATA_PyMT', 'WT_PyMT')

df_CD8_effector_T = de_analysis_for_cell_type_df(adata_T, 'Effector_CD8', 'dblGATA_PyMT', 'WT_PyMT')

df_Exh_CD4_CD8_T = de_analysis_for_cell_type_df(adata_T, 'Exh_CD4&CD8', 'dblGATA_PyMT', 'WT_PyMT')

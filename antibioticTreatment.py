#!/usr/bin/env python
# coding: utf-8

import diffxpy.api as de
import pandas as pd
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import umap
import umap.plot
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
from gff3 import Gff3
import gffutils
import scanpy as sc

#-----CIP-----------------
CIP_0h = pd.read_csv('TC-T0.selected.counts.tsv', sep='\t')
CIP_0h.columns = ["gene"] + Persister_06275_0h.columns[1:].tolist()
CIP_0h.set_index('gene', inplace=True)

CIP_1h = pd.read_csv('CPR-T1.selected.counts.tsv', sep='\t')
CIP_1h.columns = ["gene"] + Persister_06275_5h.columns[1:].tolist()
CIP_1h.set_index('gene', inplace=True)


CIP_2h = pd.read_csv('CPR-T2.selected.counts.tsv', sep='\t')
CIP_2h.columns = ["gene"] + Persister_06275_6h.columns[1:].tolist()
CIP_2h.set_index('gene', inplace=True)

CIP_4h = pd.read_csv('CPR-T4.selected.counts.tsv', sep='\t')
CIP_4h.columns = ["gene"] + Persister_06275_7h.columns[1:].tolist()
CIP_4h.set_index('gene', inplace=True)

filename = "Escherichia_coli_bw25113.ASM75055v1.46.gff3"
database_filename = "Escherichia_coli_bw25113-128"
db = gffutils.create_db(filename, database_filename, merge_strategy="warning")

genes = [i for i in [f.id for f in db.all_features()] if 'gene' in i]

gene_names = []
for gene in genes:
    gene_names.append(db[gene][8]['Name'][0])

gene_dict_rev = dict(zip(gene_names, genes))
gene_dict = dict(zip(genes, gene_names))

all_df = pd.DataFrame()
all_df['gene'] = genes

frames = [CIP_0h, CIP_1h, CIP_2h, CIP_4h]
len(frames[0].T), len(frames[1].T),len(frames[2].T), len(frames[3].T) 

times = []
t = ['0h','CIP_1h','CIP_2h','CIP_4h']
for i,frame in enumerate(frames):
    frames[i].columns = [name +'_'+ t[i] for name in frame.columns]
    times = times+ (len(frame.T)*[t[i]])

all_df = pd.merge(all_df, CIP_0h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, CIP_1h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, CIP_2h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, CIP_4h, on='gene', how='outer').fillna(0, downcast='infer')

all_df.gene = all_df.gene.apply(lambda x: gene_dict[x])
all_df.set_index('gene', inplace=True)
all_df = all_df.T
all_df['time_point'] = times
all_df.loc[:,'Total UMIs'] = all_df.drop(['time_point'],axis = 1).sum(axis=1)

sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)

all_df = all_df[(all_df['Total UMIs'] <3000) & (all_df['Total UMIs'] >100)]
all_df=all_df[-((all_df['Total UMIs'] >1500) & (all_df['time_point'] =='CIP_2h'))] 
all_df=all_df[-((all_df['Total UMIs'] >2500) & (all_df['time_point'] =='CIP_1h'))] 
all_df=all_df[-((all_df['Total UMIs'] >1000) & (all_df['time_point'] =='CIP_4h'))] 
sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)

all_df['gene_count'] = all_df.drop(['time_point','Total UMIs'],axis = 1).astype(bool).sum(axis=1)
sns.violinplot(data=all_df, x='time_point', y='gene_count', cut=0)

all_df.iloc[:, :-3].to_csv('CIP.csv')
data = sc.read('CIP.csv')
data.obs['time'] = all_df['time_point']
data.var_names_make_unique()

rRNA = pd.read_table('rRNA22.csv',sep=',')
rRNA.head()
pattern = '|'.join(rRNA['gene'].to_list())
rRNAbool = data.var_names.str.contains(pattern)
data = data[:,~rRNAbool]

sc.pl.highest_expr_genes(data, n_top=20)

sc.pp.filter_cells(data, min_genes=10)
sc.pp.filter_genes(data, min_cells=3)

sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
adata= data
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,n_bins=20)
sc.pl.highly_variable_genes(adata)

adata.obs['n_counts']=adata.X.sum(axis=1)
sc.pp.regress_out(adata, ['n_counts'])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=4)
sc.tl.umap(adata)
sc.pl.umap(adata, color='time',size=8)
sc.tl.leiden(adata, resolution =1)
sc.pl.umap(adata, color='leiden',size=8,legend_loc='on data',legend_fontsize =10)

marker_genes_dict=['ompF', 'tsx','lamB',]
sc.pl.violin(adata, marker, groupby='time', standard_scale='var',dendrogram=False)

sc.tl.rank_genes_groups(adata, 'leiden', n_genes=5)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,standard_scale ='var')

marker_genes={'SOS response':[  'sulA',  'recN','uvrB','uvrA','yebG', 'ruvA','dinI','recF',  'dinG', 'uvrD', 'lexA' ],
                   'ROS degradation':['sodB', 'sodA']}
sc.pl.dotplot(adata, marker_genes, groupby='time', standard_scale='var',dendrogram=False)

marker_genes={'TCA':['mdh', 'gltA', 'acnB', 'acnA','icd','sucA','sucB','lpd','sucD','sucC','sdhA','sdhB','sdhC'],
                   'glycine metabolism':['glyA'],
                   'glycolysis':['gapA'],
                   'glyoxylate cycle':['mdh','gltA','acnB','acnA','aceA','aceB']}
sc.pl.dotplot(adata, marker_genes, groupby='time', standard_scale='var',dendrogram=False)

#---AMP----
Persister_06275_0h = pd.read_csv('Persister_06275-0h.cell.counts.tsv', sep='\t')
Persister_06275_0h.set_index('gene', inplace=True)
Persister_06275_1h = pd.read_csv('Persister_06275-1h.cell.counts.tsv', sep='\t')
Persister_06275_1h.set_index('gene', inplace=True)
Persister_06275_2h = pd.read_csv('Persister_06275-2h.cell.counts.tsv', sep='\t')
Persister_06275_2h.set_index('gene', inplace=True)
Persister_06275_4h = pd.read_csv('Persister_06275-4h.cell.counts.tsv', sep='\t')
Persister_06275_4h.set_index('gene', inplace=True)

filename = "Escherichia_coli_bw25113.ASM75055v1.46.gff3"
database_filename = "Escherichia_coli_bw25113-35"
db = gffutils.create_db(filename, database_filename, merge_strategy="warning")

genes = [i for i in [f.id for f in db.all_features()] if 'gene' in i]

gene_names = []
for gene in genes:
    gene_names.append(db[gene][8]['Name'][0])

gene_dict_rev = dict(zip(gene_names, genes))
gene_dict = dict(zip(genes, gene_names))

all_df = pd.DataFrame()
all_df['gene'] = genes

gene=pd.DataFrame(list(gene_dict.items()))
gene.to_csv('gene.csv')

all_df.head()

frames = [Persister_06275_0h, Persister_06275_1h, Persister_06275_2h, Persister_06275_4h]
len(frames[0].T), len(frames[1].T),len(frames[2].T), len(frames[3].T) 

times = []
t = ['0h','1h','2h','4h']
for i,frame in enumerate(frames):
    frames[i].columns = [name +'_'+ t[i] for name in frame.columns]
    times = times+ (len(frame.T)*[t[i]])

all_df = pd.merge(all_df, Persister_06275_0h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, Persister_06275_1h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, Persister_06275_2h, on='gene', how='outer').fillna(0, downcast='infer')
all_df = pd.merge(all_df, Persister_06275_4h, on='gene', how='outer').fillna(0, downcast='infer')

all_df.gene = all_df.gene.apply(lambda x: gene_dict[x])

all_df.set_index('gene', inplace=True)

all_df = all_df.T
all_df['time_point'] = times
all_df.loc[:,'Total UMIs'] = all_df.drop(['time_point'],axis = 1).sum(axis=1)

sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)

all_df = all_df[(all_df['Total UMIs'] <3000) & (all_df['Total UMIs'] >300)]

sns.violinplot(data=all_df, x='time_point', y='Total UMIs', cut=0)

all_df['gene_count'] = all_df.drop(['time_point','Total UMIs'],axis = 1).astype(bool).sum(axis=1)
sns.violinplot(data=all_df, x='time_point', y='gene_count', cut=0)

all_df.iloc[:, :-3].to_csv('AMP.csv')

data = sc.read('AMP.csv')
data.obs['time'] = all_df['time_point']
data.var_names_make_unique()

rRNA = pd.read_table('rRNA22.csv',sep=',')
rRNA.head()
pattern = '|'.join(rRNA['gene'].to_list())
rRNAbool = data.var_names.str.contains(pattern)
data.obs['percent_rRNA'] = np.sum(data[:,rRNAbool].X,axis=1) / np.sum(data.X,axis=1)
sc.pl.violin(data, keys='percent_rRNA',groupby = 'time')

data = data[:,~rRNAbool]
sc.pl.highest_expr_genes(data, n_top=20)

sc.pp.filter_cells(data, min_genes=10)
sc.pp.filter_genes(data, min_cells=3)

sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
adata= data
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

adata.obs['n_counts']=adata.X.sum(axis=1)

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=17)
sc.tl.umap(adata)
sc.pl.umap(adata, color='time',size=15)
sc.tl.leiden(adata, resolution =1)
sc.pl.umap(adata, color='leiden',size=15,legend_loc='on data',legend_fontsize =10)

sc.tl.dendrogram(adata, 'leiden')
sc.tl.rank_genes_groups(adata, 'leiden', n_genes=10)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,standard_scale = 'var',groups=('6','9','10'))

sc.tl.rank_genes_groups(adata, 'leiden', n_genes=20)
result = data.uns['rank_genes_groups']
groups = result['names'].dtype.names
degs=pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(20)
degs
degs.to_csv('AMP-degs.csv')

sc.pl.correlation_matrix(data, groupby='leiden')

sc.tl.rank_genes_groups(data, 'time', n_genes=5)
sc.pl.rank_genes_groups_dotplot(data, n_genes=5,standard_scale ='var')

sc.tl.rank_genes_groups(data, 'leiden', n_genes=5)
sc.pl.rank_genes_groups_dotplot(data, n_genes=5,standard_scale ='var')

marker_genes=['glpK','aceA','osmY','fusA','ssrA','fhuA','secY','recA','sulA','ftsN']
sc.pl.umap(data, color=['recA'],size=30,legend_loc='on data',legend_fontsize =10,color_map='coolwarm')

marker_genes={'SOS response':[ 'uvrB','recF', 'uvrA', 'dinG', 'uvrD', 'lexA', 'tisB', 'cho', 'dinI', 'dinB', 'yebG', 'umuD','umuC', 'sulA', 'ruvB', 'ruvA', 'recN','recA'],
                   'ROS degradation':['sodB', 'sodA', 'katG', 'katE']}
sc.pl.dotplot(data, marker_genes, groupby='leiden', standard_scale='var',dendrogram=True)
sc.pl.dotplot(data, marker_genes, groupby='time', standard_scale='var',dendrogram=True)

marker_genes={'TCA':['mdh', 'gltA', 'acnB', 'acnA','icd','sucA','sucB','lpd','sucD','sucC','sdhA','sdhB','sdhC'],
                   'glycine metabolism':['glyA'],
                   'glycolysis':['gapA'],
                   'glyoxylate cycle':['mdh','gltA','acnB','acnA','aceA','aceB']}
sc.pl.dotplot(data, marker_genes, groupby='leiden', standard_scale='var',dendrogram=True)






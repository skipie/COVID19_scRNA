import scanpy as sc
import numpy,pandas
import matplotlib.pyplot as plt
import seaborn
genes = ['CCR7', 'CD3D', 'CD8A', 'IL7R','GZMK', 'LEF1', 'TCF7', 'AQP3', 'PRDM1','NKG7','PRF1','CD69','GNLY','GZMB',
'FCER1A','CD34','NKG7','FCGR3A','PPBP','CD79A','MS4A1','CD19','CD27','IGHD','CD14','CD68','S100A12','HBB']

#cells_df = pandas.read_csv('E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/meta.data.csv', sep=',', header=0, index_col=0)
cells_df = pandas.read_csv('/share/pub/qiuf/COVID/03-procesure/Gsemerged_data_processure/meta.data.csv', sep=',', header=0,low_memory=False)
counts_df = pandas.read_csv('/share/pub/qiuf/COVID/03-procesure/Gsemerged_data_processure/Rpoly_counts.csv', sep=',', index_col=0,
                            engine='c', na_filter=False, low_memory=False)
#counts_df = pandas.read_csv('E:/COVID/03-procesure/Gsemerged_data_processure/recluster_file/Rpoly_counts.csv', sep=',', index_col=0,header=0,
#                            engine='c', na_filter=False, low_memory=False , nrows=5)
#, nrows=5

clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

xmax={'CCR7':5, 'CD3D':5, 'CD8A':5, 'IL7R':5,'GZMK':5, 'LEF1':5, 'TCF7':5, 'AQP3':5, 'PRDM1':5,'NKG7':5,'PRF1':5,'CD69':5,'GNLY':5,'GZMB':5,
'FCER1A':5,'CD34':5,'NKG7':5,'FCGR3A':5,'PPBP':5,'CD79A':5,'MS4A1':5,'CD19':5,'CD27':5,'IGHD':5,'CD14':5,'CD68':5,'S100A12':5,'HBB':5}
cluster_color = {0:'lightskyblue', 1:'coral', 2:'dodgerblue', 3:'deepskyblue', 4:'limegreen',
                 5:'royalblue', 6:'m', 7:'mediumblue', 8:'crimson', 9:'darkgoldenrod',
                 10:'green'}
#, 11:'cornflowerblue', 12:'hotpink'
seaborn.set_context('talk')
fig, axes = plt.subplots(1, len(genes), figsize=(35,13), facecolor='w')
fig.subplots_adjust(hspace=0, wspace=0)
for igene,gene in enumerate(genes):
    gene_df = counts_df.loc[[gene]].T
    x= numpy.arange(len(gene_df.index.values))
    gene_df.index=x
    gene_df['cluster'] = cells_df.loc[x, 'RPoly']

    gene_df = gene_df.loc[gene_df['cluster'].isin(clusters)]
    seaborn.violinplot(data=gene_df, y='cluster', x=gene, orient='h', order=clusters, linewidth=0,
                       palette=cluster_color, inner=None, scale='width', cut=0, ax=axes[igene])
    axes[igene].set_xticks([])
    axes[igene].set_xlim(0, xmax[gene])
    if igene!=0:
        axes[igene].set_yticks([])
    else:
        axes[igene].set_yticklabels(list(map(str, numpy.array(clusters)+1)))
    axes[igene].set_ylabel('')
plt.savefig('/share/pub/qiuf/COVID/01-data/CELL/plot/temp_violin.pdf')
plt.show()
plt.close()



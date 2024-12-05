import pandas as pd
from numpy import trunc
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
import seaborn as sns

def extract_clustered_table(res, data):
    if res.dendrogram_row is None:
        print("Apparently, rows were not clustered.")
        return -1    
    if res.dendrogram_col is not None:
        new_cols = data.columns[res.dendrogram_col.reordered_ind]
        new_ind = data.index[res.dendrogram_row.reordered_ind]        
        return data.loc[new_ind, new_cols]    
    else:
        new_ind = data.index[res.dendrogram_row.reordered_ind]
        return data.loc[new_ind,:]

input_median = snakemake.input.median[0]
input_sd = snakemake.input.sd[1]
output_clustermap = snakemake.output.euclidean
output_sd = snakemake.output.sd
output_clustermap_cosine = snakemake.output.cosine
analysis = snakemake.wildcards.analyses
cluster_x = snakemake.params.cluster_x[analysis]
#input_median = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_36\data_by_neuron_median.csv'
#input_sd = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_36\data_by_neuron_sd.csv'

median = pd.read_csv(input_median, dtype={'neuron': str, 'n_subjects': str})
median['neuron'] = median['neuron'] + ' | n=' + median['n_subjects']
median = median.drop(columns='n_subjects')
median = median.set_index(['neuron'])
clusters = list(median.columns[median.columns.str.contains('clusters')])
clusters = median[clusters]
median = median.drop(columns=clusters)

sd = pd.read_csv(input_sd, dtype={'neuron': str, 'n_subjects': str})
sd['neuron'] = sd['neuron'] + ' | n=' + sd['n_subjects']
sd = sd.drop(columns='n_subjects')
sd = sd.set_index(['neuron'])
clusters = list(sd.columns[sd.columns.str.contains('clusters')])
clusters = sd[clusters]
sd = sd.drop(columns=clusters)

#c = int(clusters.columns.str.extract('(\d+)').max())
#cmap = dict(zip(range(0,c), sns.color_palette('rocket', c)))
#row_colors = clusters.stack().map(cmap).unstack()

x=3
cbl=3.5
if cluster_x:
    x=1.0
    cbl=2.5

w = round(len(median.columns)/x)
sns.set_theme(color_codes=True)
sns.set(font_scale=2.0)

palette = sns.color_palette('vlag', 7)
levels = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
cmap, norm = from_levels_and_colors(levels, palette, extend='both')

g=sns.clustermap(median, cmap='vlag', tree_kws=dict(linewidths=1.8),
                 metric='euclidean', method='ward', 
                 xticklabels=cluster_x, yticklabels=True, 
                 center=0, vmin=-cbl, vmax=cbl,
                 col_cluster=cluster_x, figsize=(w,22))
plt.savefig(output_clustermap, dpi=300, bbox_inches='tight')
plt.close()

plt.figure(figsize=(w+12,22))
sd_ordered = extract_clustered_table(g, sd)
m = sd_ordered.stack().median()
min_h = 0
max_h = 2
g=sns.heatmap(sd_ordered, cmap='vlag', center=m, vmin=min_h, vmax=max_h)
plt.savefig(output_sd, dpi=300, bbox_inches='tight')
plt.close()

g=sns.clustermap(median, cmap='vlag', tree_kws=dict(linewidths=1.8),
                 metric='cosine', 
                 xticklabels=cluster_x, yticklabels=True, 
                 center=0, vmin=-cbl, vmax=cbl,
                 col_cluster=cluster_x, figsize=(w,22))
plt.savefig(output_clustermap_cosine, dpi=300, bbox_inches='tight')
plt.close()
#row_colors=row_colors, 

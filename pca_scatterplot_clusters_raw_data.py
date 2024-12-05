import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_pc = snakemake.input.pc[0]
input_clusters = snakemake.input.clusters[1]
output_file = snakemake.output[0]
#input_pc = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\raw_data\pca_principal_components.csv'
#input_clusters = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\raw_data\data_raw_clustered.csv'

data = pd.read_csv(input_pc, dtype={'subject_id': str})
clusters = pd.read_csv(input_clusters, dtype={'subject_id': str})

clusters = pd.merge(clusters, data, how='inner', on='subject_id')

cols = list(clusters.columns[((clusters.columns.str.contains('clusters')) | 
        (clusters.columns.str.contains('PC')) | 
        (clusters.columns.str.contains('subject_id')))])
clusters = clusters[cols]
clusters[clusters.columns[clusters.columns.str.contains('clusters')]] = \
    clusters[clusters.columns[clusters.columns.str.contains('clusters')]]. \
    astype(int).astype(str)

other_cols = list(clusters.columns[~(clusters.columns.str.contains('clusters'))])
clusters = pd.melt(clusters, id_vars=other_cols, 
                   var_name='n_clusters', value_name='cluster')

sum_clusters = clusters.groupby(['n_clusters', 'cluster'])['subject_id'].count(). \
               reset_index(name='n')

clusters = pd.merge(clusters, sum_clusters, how='inner', on=['n_clusters', 'cluster'])
clusters['n'] = clusters['n'].astype(int).astype(str)
clusters['cluster'] = clusters['cluster'] + ' | n=' + clusters['n']
clusters = clusters.drop(columns='n')
clusters = clusters.sort_values(by=['n_clusters', 'cluster'])

#ax = sns.relplot(x='PC_1', y='PC_2', hue='cluster', col='n_clusters', 
#                 data=clusters, kind='scatter', height=5, aspect=1.0)
clust = list(clusters['n_clusters'].unique())
n_clusters = len(clust)
x = 1

f = plt.figure(figsize=(20,5))
for i, c in enumerate(clust):
    data_c = clusters[clusters['n_clusters']==c]
    f.add_subplot(x, n_clusters/x, i+1)
    sns.scatterplot(x='PC_1', y='PC_2', hue='cluster', data=data_c)

plt.savefig(output_file, bbox_inches='tight', dpi=150)
plt.close()

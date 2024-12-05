import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_file = snakemake.input.data[1]
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\raw_data\data_raw_clustered.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})

data[data.columns[data.columns.str.contains('clusters')]] = \
    data[data.columns[data.columns.str.contains('clusters')]]. \
    astype(int).astype(str)

clusters = list(data.columns[data.columns.str.contains('clusters')])
features = list(data.columns[~(data.columns.str.contains('clusters'))])

data = pd.melt(data, id_vars=features, value_vars=clusters, 
                 var_name='n_clusters', value_name='cluster')

sum_clusters = data.groupby(['n_clusters', 'cluster'])['subject_id'].count(). \
               reset_index(name='n')

data = pd.merge(data, sum_clusters, on=['n_clusters', 'cluster'], how='inner')
data['n'] = data['n'].astype(int).astype(str)
data['cluster'] = data['cluster'] + ' | n=' + data['n']
data = data.drop(columns='n')

data = pd.melt(data, id_vars=['subject_id', 'n_clusters', 'cluster'],
                 var_name='feature', value_name='value')
data = data.sort_values(by=['n_clusters', 'feature', 'cluster'])

#ax = sns.catplot(x='feature', y='value', hue='cluster', row='n_clusters', 
#                 data=data, kind='box', height=5, aspect=5.0)
clusters = list(data['n_clusters'].unique())
n_clusters = len(clusters)
x = 1
f = plt.figure(figsize=(15,10))
for i, c in enumerate(clusters):
    data_c = data[data['n_clusters']==c]
    f.add_subplot(n_clusters/x, x, i+1)    
    sns.boxplot(x='feature', y='value', hue='cluster', data=data_c)
    plt.legend(bbox_to_anchor=(1.01, 1)) #, borderaxespad=0

plt.savefig(output_file, bbox_inches='tight', dpi=150)
plt.close()

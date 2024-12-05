import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_data = snakemake.input.data[1]
output_file = snakemake.output[0]
analysis = snakemake.wildcards.analyses
cluster_x = snakemake.params.cluster_x[analysis]
#input_data = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\longitudinal\raw_data\data_clustered.csv'

data = pd.read_csv(input_data, dtype={'subject_id': str})

data = data.set_index(['subject_id'])
clusters = list(data.columns[data.columns.str.contains('clusters')])
clusters = data[clusters]
data = data.drop(columns=clusters)

c = int(clusters.columns.str.extract('(\d+)').max())
cmap = dict(zip(range(0,c), sns.color_palette('rocket', c)))
row_colors = clusters.stack().map(cmap).unstack()

x=3
if cluster_x:
    x=1
w = round(len(data.columns)/x)
sns.set_theme(color_codes=True)
sns.set(font_scale=1.5)
g=sns.clustermap(data, row_colors=row_colors, cmap='vlag', 
                 metric='euclidean', method='ward',
                 xticklabels=True, yticklabels=True, center=0,
                 col_cluster=cluster_x, figsize=(w,22))
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

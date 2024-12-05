import pandas as pd

input_som = snakemake.input.som[0]
input_som_clustered = snakemake.input.som_clustered[0]
input_data_clustered = snakemake.input.data[1]
output_file = snakemake.output[0]

#input_som = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\16_neurons\som_id_groups.csv'
#input_som_clustered = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\16_neurons\som_clusters_clustered.csv'
#input_data_clustered = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\raw_data\data_raw_clustered.csv'

som = pd.read_csv(input_som, dtype={'subject_id': str, 'neuron': str})
som = som.drop(['x', 'y'], axis=1)

som_clustered = pd.read_csv(input_som_clustered, dtype={'neuron': str})

cols_clust = list(som_clustered.columns[som_clustered.columns.str.contains('clusters')])
cols_filter = ['neuron'] + cols_clust
som_clustered = som_clustered.filter(cols_filter)

cols = som_clustered.columns[som_clustered.columns.str.contains('clusters')]
cols_new = cols + ' som'
som_clustered = som_clustered.rename(columns=dict(zip(cols, cols_new)))

data_clustered = pd.read_csv(input_data_clustered, dtype={'subject_id': str})
#data_clustered = data_clustered.filter(['subject_id', '2 clusters', '3 clusters', '5 clusters'])
cols = data_clustered.columns[data_clustered.columns.str.contains('clusters')]
cols_new = cols + ' raw'
data_clustered = data_clustered.rename(columns=dict(zip(cols, cols_new)))

merged = pd.merge(som, som_clustered, on='neuron', how='inner')
merged = pd.merge(merged, data_clustered, on='subject_id', how='inner')

merged.to_csv(output_file, header=True, index=False)

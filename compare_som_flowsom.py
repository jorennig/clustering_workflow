import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_som = snakemake.input.som[0]
input_flow_som = snakemake.input.flowsom[0]
output_neurons = snakemake.output.neurons
#input_som = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som\clusters_merged.csv'
#input_flow_som = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som\flow_som_results.csv'

som = pd.read_csv(input_som, dtype={'subject_id': str})
flowsom = pd.read_csv(input_flow_som, dtype={'subject_id': str})

som_all = pd.merge(som, flowsom, how='inner', on='subject_id')
som_all = som_all.set_index('subject_id')

som_all = som_all.filter(['neuron', '3 clusters som', 'som_id', 'meta_cluster_id'])
rename = {'neuron': 'neuron_som', '3 clusters som': 'clusters_som', 'som_id': 'neuron_flowsom', 'meta_cluster_id': 'clusters_flowsom'}
som_all = som_all.rename(columns=rename)
som_all['clusters_som'] = som_all['clusters_som'] + 1

table_neuron = pd.crosstab(som_all['neuron_som'], som_all['neuron_flowsom'])
table_neuron = table_neuron.sort_values(by=list(table_neuron.columns), 
                                        ascending=False)
plt.figure(figsize = (30,25))
sns.heatmap(table_neuron, cmap='viridis', annot=False, cbar=True)
plt.savefig(output_neurons, bbox_inches='tight', dpi=150)
plt.close()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

input_file = snakemake.input.data
input_flow_som = snakemake.input.flowsom
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_25\clusters_merged.csv'
#input_flow_som = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_25\flow_som_results.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str, 'neuron': str})
clusters = list(data.columns[data.columns.str.contains('clusters')]) + ['distance']
data = data.drop(columns=clusters)

flowsom = pd.read_csv(input_flow_som, dtype={'subject_id': str, 'som_id': str}). \
                      rename(columns={'som_id': 'neuron'})
flowsom = flowsom.drop(columns=['meta_cluster_id'])
neuron = ['24', '25']
flowsom = flowsom[flowsom['neuron'].isin(neuron)]

data = data[data['subject_id'].isin(flowsom['subject_id'])]
data = data.drop(columns=['neuron'])
data = pd.merge(data, flowsom, on='subject_id', how='inner')

data = pd.melt(data, id_vars=['subject_id', 'neuron'], 
               var_name='feature_name', value_name='numeric_value'). \
               sort_values(by=['neuron', 'subject_id', 'feature_name'])
data[['feature', 'period']] = data['feature_name'].str.split('_', 1, expand=True)
data = data.drop(columns='feature_name')
data['size'] = 1

data_sum = data.groupby(['neuron', 'feature', 'period']).median().reset_index()
data_sum['subject_id'] = 'median_neuron_' + data_sum['neuron']
data_sum['size'] = 1.5

data = data.append(data_sum)

sns.relplot(data=data, x='period', y='numeric_value', col='feature', 
            hue='subject_id', style='size', kind='line', 
            size='size', row='neuron')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

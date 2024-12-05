import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

input_file = snakemake.input.data
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_25\clusters_merged.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str, 'neuron': str})
clusters = list(data.columns[data.columns.str.contains('clusters')]) + ['distance']
data = data.drop(columns=clusters)

neuron = ['22', '23', '03', '18', '12', '08', '09', '13', '17']

data = data[data['neuron'].isin(neuron)]

data = pd.melt(data, id_vars=['subject_id', 'neuron'], 
               var_name='feature_name', value_name='numeric_value'). \
               sort_values(by=['neuron', 'subject_id', 'feature_name'])
data[['feature', 'period']] = data['feature_name'].str.split('_', 1, expand=True)
data = data.drop(columns='feature_name')

data_sum = data.groupby(['neuron', 'feature', 'period']).median().reset_index()
data_sum['subject_id'] = 'median_neuron_' + data_sum['neuron']
data_sum['neuron'] = 'median'

#data = data.append(data_sum)

sns.relplot(data=data_sum, x='period', y='numeric_value', col='feature', 
            hue='subject_id', kind='line')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

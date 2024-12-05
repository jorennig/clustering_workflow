import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_file = snakemake.input.data[0]
output_py_flow = snakemake.output.pf
output_change_slope = snakemake.output.cs
#input_file = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/som_results.csv'

dtype={'subject_id': str, 'neuron': str}
data = pd.read_csv(input_file, dtype=dtype)

check_pairs = {'pysom_baseline': 'flowsom_baseline', 
               'pysom_change': 'flowsom_change',
               'pysom_slopes': 'flowsom_slopes'}
n_com = len(check_pairs.keys())

f = plt.figure(figsize=(20,4))
for i, (key, value) in enumerate(check_pairs.items()):    
    data_c = data[(data['analysis']==key) | (data['analysis']==value)]
    data_t = pd.pivot_table(data_c, index='subject_id', values='cluster', columns='analysis').dropna()
        
    table = pd.crosstab(data_t[key], data_t[value])
    cols = list(table.columns)
    table = table.sort_values(by=cols[0], ascending=False)
    
    f.add_subplot(1, n_com, i+1)    
    sns.heatmap(table, cmap='viridis', annot=True, cbar=False)

plt.savefig(output_py_flow, bbox_inches='tight', dpi=150)
plt.close()


check_pairs = {'pysom_change': 'pysom_slopes',
               'flowsom_change': 'flowsom_slopes'}
n_com = len(check_pairs.keys())
f = plt.figure(figsize=(15,4))
for i, (key, value) in enumerate(check_pairs.items()):    
    data_c = data[(data['analysis']==key) | (data['analysis']==value)]
    data_t = pd.pivot_table(data_c, index='subject_id', values='cluster', columns='analysis').dropna()
    
    table = pd.crosstab(data_t[key], data_t[value])
    cols = list(table.columns)
    table = table.sort_values(by=cols[0], ascending=False)
    
    f.add_subplot(1, n_com, i+1)    
    sns.heatmap(table, cmap='viridis', annot=True, cbar=False)

plt.savefig(output_change_slope, bbox_inches='tight', dpi=150)
plt.close()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

input_file = snakemake.input.data[0]
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\16_neurons\clusters_merged.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})

clusters = list(data.columns[data.columns.str.contains('cluster')])
n_clusters = len(clusters)
x = 2

f = plt.figure(figsize=(10,6))
for i, c in enumerate(clusters):
    table = round(pd.crosstab(data['disease group'], data[c], normalize='index')*100).T
    #table = pd.crosstab(data[c], data['disease group'])
    table = table.sort_values(by=['control', 'pre-manifest hd', 'stage 1 hd', 
                                  'stage 2 hd'], ascending=False)
    f.add_subplot(x, n_clusters/x, i+1)    
    sns.heatmap(table, cmap='viridis', annot=True, cbar=False)

plt.savefig(output_file, bbox_inches='tight', dpi=150)
plt.close()

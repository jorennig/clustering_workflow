import pandas as pd
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

input_file = snakemake.input.data[0]
output_heatmap = snakemake.output.bar

input_file = r'//pstore.bas.roche.com/data/rpmda/HD242932/jira/RPMDA-6476-som-implementation/results/9_neurons/som_clusters_clustered.csv'

data = pd.read_csv(input_file) #, dtype={'subject_id': str}

clusters = list(data.columns[data.columns.str.contains('cluster')])
data[clusters] = data[clusters].astype(int).astype(str)

for c in clusters:
    table = pd.crosstab(data['neuron'], data[c])
    
    ax = sns.heatmap(table, cmap='viridis', annot=True, cbar=False)
    bottom, top = ax.get_ylim()
    plt.ylabel('neuron')
#    output_file = path_output + 'heatmap' + c + '_raw_data.png'
#    plt.savefig(output_file, dpi=150, bbox_inches='tight')
#    plt.close()

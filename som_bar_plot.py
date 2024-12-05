import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def percent_group(chunk, col_name):    
    name = col_name + ' percent'
    chunk[name] = chunk.groupby(col_name) \
                       [col_name].transform('count')
    chunk[name] = chunk[name]/len(chunk)*100
    return chunk

input_file = snakemake.input.data[0]
output_file = snakemake.output.bar
#input_file = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\16_neurons\som_id_groups.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.drop(columns=['subject_id', 'distance']) #'x', 'y', 
#data['disease code'] = data['disease code']
#data['neuron'] = data['neuron']

data = data.groupby(['neuron']).apply(percent_group, ('disease code'))
#data = data.groupby(['disease code']).apply(percent_group, ('neuron'))
data = data.drop_duplicates(subset=['neuron', 'disease code'], keep='first')

data_d = pd.pivot_table(data, values='disease code percent', 
                        index=['neuron'], 
                        columns=['disease group'])
data_d = data_d.fillna(0)

colors = ['#EDB233', '#90C3EC', '#C02942', '#79BD9A']
c_map = sns.set_palette(sns.color_palette(colors))

data_d = data_d.sort_values(by=['control', 'pre-manifest hd', 
                                'stage 1 hd', 'stage 2 hd'], ascending=False)
ax = data_d.plot(kind='bar', stacked=True, colormap=c_map, legend=True)
ax.legend(bbox_to_anchor=(1.0, 1.0))
plt.savefig(output_file, bbox_inches='tight', dpi=150)
plt.close()

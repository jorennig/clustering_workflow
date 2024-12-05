import pandas as pd
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

def percent_group(chunk, col_name):
    name = col_name + ' percent'
    chunk[name] = chunk.groupby(col_name) \
                       [col_name].transform('count')
    chunk[name] = chunk[name]/len(chunk)*100
    return chunk

input_data = snakemake.input.data[1]
input_demo = snakemake.input.demo

#input_data = r'//pstore.bas.roche.com/data/rpmda/HD242932/jira/RPMDA-6476-som-implementation/results/raw_data/data_raw_clustered.csv'
#input_demo = r'//pstore/data/rpmda/HD242932/clinical_data/standardised/results/demographics.csv'

data = pd.read_csv(input_data, dtype={'subject_id': str})
clusters = list(data.columns[data.columns.str.contains('cluster')])
#data[clusters] = data[clusters].astype(int).astype(str)

demo = pd.read_csv(input_demo, dtype={'subject_id': str})
groups = demo.filter(['subject_id', 'disease group'])
replace = {'control': 0, 'pre-manifest hd': 1, 'stage 1 hd': 2, 'stage 2 hd': 3}
groups['disease code'] = groups['disease group'].replace(replace)

data = pd.merge(data, groups, on='subject_id', how='inner')

replace = {0:1, 1:0}
data['2 clusters'] = data['2 clusters'].replace(replace)

colors = ['#EDB233', '#90C3EC', '#C02942', '#79BD9A']
c_map = sns.set_palette(sns.color_palette(colors))

for c in clusters:
    data = data.groupby([c]).apply(percent_group, ('disease code'))
    
    data_d = data.drop_duplicates(subset=[c, 'disease code'], keep='first')
    
    data_d = pd.pivot_table(data_d, values='disease code percent', 
                            index=[c], 
                            columns=['disease group'])
    data_d = data_d.fillna(0)
        
    data_d = data_d.sort_values(by=['control', 'pre-manifest hd', 
                                'stage 1 hd', 'stage 2 hd'], ascending=False)
    ax = data_d.plot(kind='bar', stacked=True, colormap=c_map, legend=True)
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    
    study = 'HD242932'
    ticket = 'RPMDA-6476-som-implementation'
    output_bar = r'//pstore/data/rpmda/' + study + '/jira/' + ticket + '/results/raw_data/groups_per_cluster_' + c + '.png'
    
    plt.savefig(output_bar, bbox_inches='tight', dpi = 150)
    plt.close()

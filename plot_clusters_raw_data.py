import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

input_data = snakemake.input.data[1]
#input_data = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\raw_data\data_raw_clustered.csv'

study = 'HD242932'
ticket = 'RPMDA-6476-som-implementation'
path_output = r'//pstore/data/rpmda/' + study + '/jira/' + ticket + '/results/raw_data/'

data = pd.read_csv(input_data, dtype={'subject_id': str})
clusters = list(data.columns[data.columns.str.contains('cluster')])
data[clusters] = data[clusters].astype(int).astype(str)

for c in clusters:
    output_file = path_output + 'pairplot_' + c + '_raw_data.png'
    sns.set_theme(style='ticks')
    sns.pairplot(data, hue=c, corner=True, plot_kws={'alpha':0.8})    
    plt.savefig(output_file)
    plt.close()

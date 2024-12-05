import pandas as pd

input_data = snakemake.input.data[1]
input_dm = snakemake.input.dm[0]
output_file = snakemake.output[0]
#input_data = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\36_neurons\som_id.csv'
#input_dm = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6476-som-implementation\results\36_neurons\neuron_distance.csv'

data = pd.read_csv(input_data, dtype={'subject_id': str, 'neuron': str})
dist = pd.read_csv(input_dm)

data = pd.merge(data, dist, on=['x', 'y'], how='inner'). \
                            sort_values(by='subject_id')

data.to_csv(output_file, index=False)

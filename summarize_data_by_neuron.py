import pandas as pd

input_file = snakemake.input.data
input_som = snakemake.input.som[0]
output_median = snakemake.output.median
output_sd = snakemake.output.sd
#input_file = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6837-feature-extraction-genhd-1\data/data_longitudinal.csv'
#input_som = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal\som_36\som_id_distance.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})

som = pd.read_csv(input_som, dtype={'subject_id': str, 'neuron': str})
som = som.drop(['x', 'y'], axis=1)

data = pd.merge(data, som, on='subject_id', how='inner')

data = pd.melt(data, id_vars=['subject_id', 'neuron', 'distance'], 
                                var_name='feature_name', 
                                value_name='numeric_value')

data_sum = data.groupby(['feature_name', 'neuron']). \
            agg(numeric_value=('numeric_value', 'median'),
                sd=('numeric_value', 'std'),
                n_subjects=('numeric_value', 'count')).reset_index(). \
                sort_values(by=['feature_name', 'neuron'])
data_sum['sd'] = data_sum['sd'].fillna(0)

data_median = pd.pivot_table(data_sum, values='numeric_value', 
                      index=['neuron', 'n_subjects'], 
                      columns='feature_name').reset_index()

data_sd = pd.pivot_table(data_sum, values='sd', 
                      index=['neuron', 'n_subjects'], 
                      columns='feature_name').reset_index()

data_median.to_csv(output_median, header=True, index=False)
data_sd.to_csv(output_sd, header=True, index=False)

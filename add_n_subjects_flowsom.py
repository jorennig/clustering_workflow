import pandas as pd

input_som_id = snakemake.input.som_id[0]
input_median = snakemake.input.median[1]
input_sd = snakemake.input.sd[2]
output_median = snakemake.output.median
output_sd = snakemake.output.sd
#input_som_id = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_144\flow_som_results.csv'
#input_median = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_144\flow_som_feature_medians_nodes.csv'
#input_sd = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\longitudinal_norm\som_144\flow_som_feature_std_nodes.csv'

data = pd.read_csv(input_som_id, dtype={'subject_id': str, 'som_id': str}).rename(columns={'som_id':'neuron'})
data['neuron'] = data['neuron'].str.zfill(2)

data_sum = data.groupby(['neuron']). \
            agg(n_subjects=('subject_id', 'count')).reset_index(). \
            sort_values(by=['neuron'])

median = pd.read_csv(input_median, dtype={'Unnamed: 0': str}).rename(columns={'Unnamed: 0':'neuron'})
median['neuron'] = median['neuron'].str.zfill(2)
median.columns = median.columns.str.replace('X','')
median = pd.merge(median, data_sum, on='neuron', how='inner')
median = median.dropna()
median.to_csv(output_median, header=True, index=False)

sd = pd.read_csv(input_sd, dtype={'Unnamed: 0': str}).rename(columns={'Unnamed: 0':'neuron'})
sd['neuron'] = sd['neuron'].str.zfill(2)
sd.columns = sd.columns.str.replace('X','')
sd = pd.merge(sd, data_sum, on='neuron', how='inner')
sd = sd.fillna(0)
sd.to_csv(output_sd, header=True, index=False)

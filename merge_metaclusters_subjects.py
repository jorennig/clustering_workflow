import pandas as pd

output_file = snakemake.output[0]

input_meta_clusters = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/meta_clusters.csv'
input_pysom_baseline = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/baseline/som_36/som_id.csv'
input_pysom_change = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/change_from_baseline/som_25/som_id.csv'
input_pysom_slopes = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/slopes/som_25/som_id.csv'
input_flowsom_baseline = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/baseline/som_36/flow_som_results.csv'
input_flowsom_change = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/change_from_baseline/som_25/flow_som_results.csv'
input_flowsom_slopes = r'//pstore/data/rpmda/BN40423-v2/jira/RPMDA-6838-clustering-genhd-1/results/slopes/som_25/flow_som_results.csv'

dtype={'subject_id': str, 'neuron': str, 'som_id': str}
meta_clusters = pd.read_csv(input_meta_clusters, dtype=dtype)

pysom_baseline = pd.read_csv(input_pysom_baseline, dtype=dtype).drop(['x', 'y'], axis=1)
pysom_baseline['analysis'] = 'pysom_baseline'

pysom_change = pd.read_csv(input_pysom_change, dtype=dtype).drop(['x', 'y'], axis=1)
pysom_change['analysis'] = 'pysom_change'

pysom_slopes = pd.read_csv(input_pysom_slopes, dtype=dtype).drop(['x', 'y'], axis=1)
pysom_slopes['analysis'] = 'pysom_slopes'

flowsom_baseline = pd.read_csv(input_flowsom_baseline, dtype=dtype). \
        drop(['meta_cluster_id'], axis=1).rename(columns={'som_id': 'neuron'})
flowsom_baseline['analysis'] = 'flowsom_baseline'
flowsom_baseline['neuron'] = flowsom_baseline['neuron'].str.zfill(2)

flowsom_change = pd.read_csv(input_flowsom_change, dtype=dtype). \
        drop(['meta_cluster_id'], axis=1).rename(columns={'som_id': 'neuron'})
flowsom_change['analysis'] = 'flowsom_change'
flowsom_change['neuron'] = flowsom_change['neuron'].str.zfill(2)

flowsom_slopes = pd.read_csv(input_flowsom_slopes, dtype=dtype). \
        drop(['meta_cluster_id'], axis=1).rename(columns={'som_id': 'neuron'})
flowsom_slopes['analysis'] = 'flowsom_slopes'
flowsom_slopes['neuron'] = flowsom_slopes['neuron'].str.zfill(2)

som_results = pd.concat([pysom_baseline, pysom_change, pysom_slopes, 
                           flowsom_baseline, flowsom_change, flowsom_slopes])

som_results = pd.merge(som_results, meta_clusters, on=['neuron', 'analysis'], 
                       how='inner').sort_values(by=['analysis', 'cluster']). \
                       reset_index(drop=True)

som_results.to_csv(output_file, header=True, index=False)

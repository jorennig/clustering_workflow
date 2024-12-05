import pandas as pd

input_file = snakemake.input.data[0]
input_demo = snakemake.input.demo
output_file = snakemake.output[0]

data = pd.read_csv(input_file, dtype={'subject_id': str})

demo = pd.read_csv(input_demo, dtype={'subject_id': str})
groups = demo.filter(['subject_id', 'disease group'])
replace = {'control': 0, 'pre-manifest hd': 1, 'stage 1 hd': 2, 'stage 2 hd': 3}
groups['disease code'] = groups['disease group'].replace(replace)

data = pd.merge(data, groups, on='subject_id', how='inner'). \
                                sort_values(by='subject_id')

data.to_csv(output_file, index=False)

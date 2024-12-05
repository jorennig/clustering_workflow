import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_file = snakemake.input.data
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6837-feature-extraction-genhd-1\data\data_baseline.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.set_index(['subject_id'])

r = int(np.floor(len(list(data.columns))/2))
c = int(np.floor(len(list(data.columns))/2))

if len(list(data.columns)) > 10:
    cols = pd.DataFrame(data.columns)
    cols = cols[0].str.split(pat='_', expand=True).rename(columns={0:'f', 1:'p'})
    r = cols['p'].nunique()
    c = cols['f'].nunique()

ax = data.hist(bins=10,) # layout=(c, r), figsize=(r*2, c*2)
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

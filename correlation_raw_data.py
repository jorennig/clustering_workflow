import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

input_file = snakemake.input.data
output_file = snakemake.output[0]
#input_file = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6400-clustering-extract-normalize-features\data\data_baseline.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.set_index(['subject_id'])

cols = list(data.columns)
r, p = stats.spearmanr(data)
r = pd.DataFrame(r)
r.columns = cols
r.index = cols
r = r.round(3)

mask = np.triu(np.ones_like(r, dtype=bool))
n_unique = len(pd.unique(r.round(1).values.ravel('K')))-1
cmap = sns.color_palette('viridis', n_unique)
ax = sns.heatmap(r, cmap=cmap, mask=mask)
plt.yticks(rotation=0)
plt.savefig(output_file, dpi=150, bbox_inches='tight')
plt.close()

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

input_file = snakemake.input.data
output_pc = snakemake.output.pc
output_loadings = snakemake.output.loadings
output_components = snakemake.output.components
output_cumulative = snakemake.output.cumulative

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.set_index(['subject_id'])
cols = list(data.columns)

vals = data.values
pca = PCA()
principal_components = pca.fit_transform(vals)
variance_explained = pca.explained_variance_ratio_
variance_explained = pd.DataFrame(variance_explained)*100

pc_cols = 'PC_' + pd.DataFrame(range(1, len(variance_explained)+1)).astype(str)
variance_explained.index = list(pc_cols[0])
variance_explained = variance_explained[0:6]

variance_cumulative = pd.DataFrame(np.cumsum(pca.explained_variance_ratio_))

principal_components =  pd.DataFrame(principal_components, 
                                    columns=list(pc_cols[0]))
principal_components['subject_id'] = data.index
principal_components.to_csv(output_pc, header=True, index=False)

loadings = (pca.components_.T * np.sqrt(pca.explained_variance_)).T
loadings = pd.DataFrame(loadings, index=list(pc_cols[0]), columns=cols)
loadings.to_csv(output_loadings, header=True, index=True)

ax = variance_explained.plot.bar(legend=False)
ax.set_ylabel('percent variance explained')
plt.savefig(output_components, bbox_inches='tight', dpi=150)
plt.close()

ax = variance_cumulative.plot(legend=False)
ax.set_xlabel('number of components')
ax.set_ylabel('cumulative explained variance')
plt.savefig(output_cumulative, bbox_inches='tight', dpi=150)
plt.close()

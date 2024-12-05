import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering

input_file = snakemake.input.data
n_clust = snakemake.params.n_clust
output_dendro = snakemake.output.dendro
output_file = snakemake.output.cluster
#input_file = r'//pstore/data/rpmda/HD242932/jira/RPMDA-6400-clustering-extract-normalize-features/data/data_baseline.csv'

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.set_index(['subject_id'])

plt.figure(figsize=(10, 7))
dend = shc.dendrogram(shc.linkage(data, method='ward', metric='euclidean'))
plt.savefig(output_dendro, bbox_inches='tight', dpi=150)
plt.close()

clusters = np.empty([len(data), len(n_clust)])
for i,c in enumerate(n_clust):
    cluster = AgglomerativeClustering(n_clusters=c, affinity='euclidean', 
                                      linkage='ward')
    clusters[:,i] = cluster.fit_predict(data)

names = list(map(str, n_clust))
clusters = pd.DataFrame(clusters, columns=names)
clusters.columns = clusters.columns + ' clusters'

data = pd.concat([data.reset_index(), clusters], axis=1)

data.to_csv(output_file, header=True, index=False)

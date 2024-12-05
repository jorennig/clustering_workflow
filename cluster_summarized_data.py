import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as shc
from sklearn.cluster import AgglomerativeClustering

input_data = snakemake.input.data[0]
n_clust = snakemake.params.n_clust
output_dendro = snakemake.output.dendro
output_file = snakemake.output.cluster
#input_data = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\data_by_neuron.csv'
#output_dendro = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\clustering_som_dendrogram.png'
#output_file = r'\\pstore.bas.roche.com\data\rpmda\HD242932\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\som_neurons_clustered.csv'

data = pd.read_csv(input_data, dtype={'neuron': str, 'n_subjects': str})

index_cols = data[['neuron', 'n_subjects']]
data['neuron'] = data['neuron'] + ' | n=' + data['n_subjects']
data = data.drop(columns='n_subjects')
data = data.set_index(['neuron'])

plt.figure(figsize=(10, 7))
dend = shc.dendrogram(shc.linkage(data, method='ward', metric='euclidean'))
plt.savefig(output_dendro, bbox_inches='tight', dpi = 150)
plt.close()

#n_clust = [2, 3, 4]
clusters = np.empty([len(data), len(n_clust)])
for i,c in enumerate(n_clust):
    cluster = AgglomerativeClustering(n_clusters=c, affinity='euclidean', 
                                      linkage='ward')
    clusters[:,i] = cluster.fit_predict(data)

names = list(map(str, n_clust))
clusters = pd.DataFrame(clusters, columns=names)
clusters.columns = clusters.columns + ' clusters'

data = pd.concat([data.reset_index(drop=True), clusters], axis=1)
data['neuron'] = index_cols['neuron']
data['n_subjects'] = index_cols['n_subjects']

data.to_csv(output_file, header=True, index=False)

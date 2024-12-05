import pandas as pd
import numpy as np
from som import SOM
from scipy.spatial.distance import cdist

def winner(som_map, vector):
    matDist = np.zeros(np.shape(som_map)[:2])
    for x in range(matDist.shape[1]):
        for y in range(matDist.shape[0]):
            matDist[y, x] = cdist([som_map[y, x, :].ravel()], [vector.ravel()], 
                                  metric='euclidean')[0][0]
    idx = np.argmin(matDist)
    return idx

input_file = snakemake.input.data
output_som = snakemake.output.som
output_cluster = snakemake.output.cluster
n_neurons = int(snakemake.wildcards.n)

data = pd.read_csv(input_file, dtype={'subject_id': str})
data = data.set_index(['subject_id'])
array = data.to_numpy()

nr = int(np.sqrt(n_neurons))
som = SOM(nr, nr)
som.fit(array, 30000, save_e=True, interval=500)

df_som = pd.DataFrame(data=np.array([winner(som.map, vec) for vec in array]),
                      columns=['neuron'], index=data.index).reset_index() # double check
winner_neurons = pd.DataFrame(som.winner_neurons(array), 
                              columns = ['x','y'])
df_som = pd.concat([df_som, winner_neurons], axis=1)
df_som['subject_id'] = df_som['subject_id'].astype(str).str.zfill(4)
df_som['neuron'] = df_som['neuron'].astype(str).str.zfill(2)

som.save(output_som)
df_som.to_csv(output_cluster, header=True, index=False)

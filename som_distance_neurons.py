import pandas as pd
import numpy as np
from som import SOM

input_som = snakemake.input.som[0]
n_neurons = int(snakemake.wildcards.n)
output_dist = snakemake.output[0]

nr = int(np.sqrt(n_neurons))
som = SOM(nr,nr)
som.load(input_som)

som.distance_map(metric='euclidean')
dm = pd.DataFrame(np.array(som.distmap))
dm['y'] = dm.index
dm['y'] = dm['y'].values[::-1].astype(str)

dm = pd.melt(dm, id_vars='y', var_name='x', value_name='distance'). \
             sort_values(by=['x', 'y']).reset_index()
dm = dm.filter(['x', 'y', 'distance'])

dm.to_csv(output_dist, header=True, index=False)

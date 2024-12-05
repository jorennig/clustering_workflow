import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from som import SOM

input_data = snakemake.input.data
input_som = snakemake.input.som[0]
n_neurons = int(snakemake.wildcards.n)

output_error = snakemake.output.error
output_distance = snakemake.output.distance
output_point_map = snakemake.output.point
#input_data = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6837-feature-extraction-genhd-1\data\data_baseline.csv'

data = pd.read_csv(input_data, dtype={'subject_id': str}). \
                                sort_values(by=['subject_id'])
data = data.set_index(['subject_id'])
array = data.to_numpy()

nr = int(np.sqrt(n_neurons))
som = SOM(nr,nr)
som.load(input_som)

som.plot_error_history()
plt.savefig(output_error, bbox_inches='tight')
plt.close()

som.plot_distance_map()
plt.savefig(output_distance, bbox_inches='tight')
plt.close()

group = np.zeros(len(array), dtype=int)
som.plot_point_map(array, group, 'n')
plt.savefig(output_point_map, bbox_inches='tight')
plt.close()

import pandas as pd

output_file = snakemake.output[0]

col_names = ['neuron','cluster']

pysom_baseline = {12:1, 16:1, 9:1, 10:1, 15:1, 22:1, 17:1, 6:1, 11:1, 4:1, 5:1, 7:1, 8:1,
                  18:2, 24:2, 29:2, 13:2, 23:2, 20:2, 26:2, 19:2, 25:2, 31:2,
                  32:3, 2:3, 1:3, 30:3, 28:3, 0:3, 35:3, 27:3, 3:3, 33:3, 34:3, 14:3, 21:3}
pysom_baseline = pd.DataFrame(list(pysom_baseline.items()), columns=col_names) 
pysom_baseline['analysis'] = 'pysom_baseline'

pysom_change = {22:1, 23:1, 3:1, 18:1, 12:1, 8:1, 9:1, 13:1, 17:1, 
                1:2, 5:2, 0:2, 6:2,
                7:3, 21:3, 16:3, 10:3, 11:3, 2:3, 19:3, 24:3, 4:3, 20:3,
                14:4, 15:4}
pysom_change = pd.DataFrame(list(pysom_change.items()), columns=col_names) 
pysom_change['analysis'] = 'pysom_change'

pysom_slopes = {7:1, 3:1, 1:1, 9:1, 0:1, 10:1, 5:1, 6:1, 
                11:2, 13:2, 
                4:3, 20:3, 8:3, 23:3, 
                22:4, 2:4, 21:4, 16:4, 12:4, 17:4,
                19:5, 14:5, 18:5, 15:5, 24:5}
pysom_slopes = pd.DataFrame(list(pysom_slopes.items()), columns=col_names) 
pysom_slopes['analysis'] = 'pysom_slopes'

flowsom_baseline = {6:1, 4:1, 5:1, 10:1, 11:1,
                    3:2, 8:2, 9:2, 1:2, 2:2, 13:2, 7:2, 14:2,
                    31:3, 32:3, 33:3, 12:3,
                    19:4, 25:4, 15:4, 16:4, 27:4, 28:4, 20:4, 21:4, 26:4,
                    18:5, 17:5, 22:5, 23:5,
                    29:6, 24:6, 30:6, 36:6, 34:6, 35:6}
flowsom_baseline = pd.DataFrame(list(flowsom_baseline.items()), columns=col_names) 
flowsom_baseline['analysis'] = 'flowsom_baseline'

flowsom_change = {1:1, 21:1, 22:1, 6:1, 11:1, 7:1, 13:1, 16:1, 12:1, 17:1, 
                  20:2, 25:2, 
                  19:3, 18:3, 23:3, 24:3, 
                  3:4, 2:4, 8:4, 9:4, 
                  15:5, 10:5, 14:5, 4:5, 5:5}
flowsom_change = pd.DataFrame(list(flowsom_change.items()), columns=col_names) 
flowsom_change['analysis'] = 'flowsom_change'

flowsom_slopes = {7:1, 3:1, 1:1, 9:1, 0:1, 10:1, 5:1, 6:1, 
                 11:2, 13:2, 
                 4:3, 20:3, 8:3, 23:3, 
                 22:4, 2:4, 21:4, 16:4, 12:4, 17:4,
                 19:5, 14:5, 18:5, 15:5, 24:5}
flowsom_slopes = pd.DataFrame(list(flowsom_slopes.items()), columns=col_names) 
flowsom_slopes['analysis'] = 'flowsom_slopes'

meta_clusters = pd.concat([pysom_baseline, pysom_change, pysom_slopes,
                           flowsom_baseline, flowsom_change, flowsom_slopes])
meta_clusters['neuron'] = meta_clusters['neuron'].astype(str).str.zfill(2)

meta_clusters.to_csv(output_file, header=True, index=False)

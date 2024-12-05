import pandas as pd
import matplotlib.pyplot as plt 

input_coord = snakemake.input.coord[4]
input_sizes = snakemake.input.sizes[5]
input_graph = snakemake.input.graph[6]
output_file = snakemake.output[0]
#input_coord = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\mst_coordinates.csv'
#input_sizes = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\mst_node_sizes.csv'
#input_graph = r'\\pstore.bas.roche.com\data\rpmda\BN40423-v2\jira\RPMDA-6838-clustering-genhd-1\results\baseline\som_36\mst_graph.txt'

coord = pd.read_csv(input_coord).rename(columns={'Unnamed: 0':'som_id'})
size = pd.read_csv(input_sizes).rename(columns={'Unnamed: 0':'som_id', 'x':'size'})
graph = pd.read_csv(input_graph, sep=' ')
graph = graph.rename(columns={graph.columns[0]: 's', graph.columns[1]: 'e'})
graph = graph + 1
df = pd.DataFrame([[1, 2]], columns=list('se'))
graph = graph.append(df)

mst = pd.merge(coord, size, on='som_id', how='inner')

fig = plt.figure(figsize=(5, 4))
fig.set(dpi=300)
plt.scatter(mst['V1'], mst['V2'], mst['size']**2)
for i, txt in enumerate(mst['som_id']):
    plt.annotate(txt, (mst['V1'][i], mst['V2'][i]))
       
for i, row in graph.iterrows():
    print(row.T)
    s = row['s']
    e = row['e']
    coord_s = coord[coord['som_id'] == s]
    coord_e = coord[coord['som_id'] == e]
    coord_a = pd.concat([coord_s, coord_e], axis=0)
    x = coord_a['V1'][0:2] 
    y = coord_a['V2'][0:2]
    plt.plot(x, y, color='blue')

plt.tight_layout()
plt.savefig(output_file, bbox_inches="tight")

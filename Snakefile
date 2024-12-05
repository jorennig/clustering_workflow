configfile: 'config.yaml'

rule histogram_raw_data:
    input: script='scripts/histogram_raw_data.py',
           data=lambda w: config['input_files'][w.analyses]
    output: 'results/{analyses}/raw_data/histogram.png'
    script: 'scripts/histogram_raw_data.py'

rule correlation_raw_data:
    input: script='scripts/correlation_raw_data.py',
           data=lambda w: config['input_files'][w.analyses]
    output: 'results/{analyses}/raw_data/heatmap_correlations.png'
    script: 'scripts/correlation_raw_data.py'

rule pca_raw_data:
    input: script='scripts/pca_raw_data.py',
           data=lambda w: config['input_files'][w.analyses]
    output: pc='results/{analyses}/raw_data/pca_principal_components.csv',
            loadings='results/{analyses}/raw_data/pca_factor_loadings.csv',
            components='results/{analyses}/raw_data/pca_components.png',
            cumulative='results/{analyses}/raw_data/pca_cumulative.png',
    script: 'scripts/pca_raw_data.py'

rule cluster_raw_data:
    input: script='scripts/cluster_raw_data.py',
           data=lambda w: config['input_files'][w.analyses]
    params: n_clust=config['n_clust']
    output: dendro='results/{analyses}/raw_data/dendrogram.png',
            cluster='results/{analyses}/raw_data/data_clustered.csv'
    script: 'scripts/cluster_raw_data.py'

rule clustermap_raw_data:
    input: script='scripts/clustermap_raw_data.py',
           data=rules.cluster_raw_data.output
    params: cluster_x=config['cluster_x']
    output: 'results/{analyses}/raw_data/clustermap.png'
    script: 'scripts/clustermap_raw_data.py'

rule box_plot_clusters_raw:
    input: script='scripts/box_plot_clusters_raw.py',
           data=rules.cluster_raw_data.output
    output: 'results/{analyses}/raw_data/box_plot_clusters.png'
    script: 'scripts/box_plot_clusters_raw.py'

rule pca_scatterplot_clusters_raw_data:
    input: script='scripts/pca_scatterplot_clusters_raw_data.py',
           pc=rules.pca_raw_data.output,
           clusters=rules.cluster_raw_data.output
    output: 'results/{analyses}/raw_data/pca_scatterplot_clusters.png'
    script: 'scripts/pca_scatterplot_clusters_raw_data.py'


rule som_clustering:
    input: script='scripts/som_clustering.py',
           data=lambda w: config['input_files'][w.analyses]
    output: som='results/{analyses}/som_{n}/som.pickle',
            cluster='results/{analyses}/som_{n}/som_id.csv'
    script: 'scripts/som_clustering.py'

rule som_distance_neurons:
    input: script='scripts/som_distance_neurons.py',
           som=rules.som_clustering.output
    output: 'results/{analyses}/som_{n}/neuron_distance.csv'
    script: 'scripts/som_distance_neurons.py'

rule merge_distance:
    input: script='scripts/merge_distance.py',
           data=rules.som_clustering.output,
           dm=rules.som_distance_neurons.output
    output: 'results/{analyses}/som_{n}/som_id_distance.csv'
    script: 'scripts/merge_distance.py'

rule som_plot:
    input: script='scripts/som_plot.py',
           data=lambda w: config['input_files'][w.analyses],
           som=rules.som_clustering.output
    output: error='results/{analyses}/som_{n}/SOM_training_error.png',
            distance='results/{analyses}/som_{n}/SOM_distance_map.png',
            point='results/{analyses}/som_{n}/SOM_point_map.png'
    script: 'scripts/som_plot.py'

rule summarize_data_by_neuron:
    input: script='scripts/summarize_data_by_neuron.py',
           data=lambda w: config['input_files'][w.analyses],
           som=rules.merge_distance.output
    output: median='results/{analyses}/som_{n}/data_by_neuron_median.csv',
            sd='results/{analyses}/som_{n}/data_by_neuron_sd.csv',
    script: 'scripts/summarize_data_by_neuron.py'

rule cluster_summarized_data:
    input: script='scripts/cluster_summarized_data.py',
           data=rules.summarize_data_by_neuron.output
    params: n_clust=config['n_clust']
    output: cluster='results/{analyses}/som_{n}/data_by_neuron_clustered.csv',
            dendro='results/{analyses}/som_{n}/clustering_som_dendrogram.png'
    script: 'scripts/cluster_summarized_data.py'

rule box_plot_clusters_som:
    input: script='scripts/box_plot_clusters_som.py',
           data=rules.cluster_summarized_data.output
    output: 'results/{analyses}/som_{n}/box_plot_clusters_som.png'
    script: 'scripts/box_plot_clusters_som.py'

rule clustermap_data_som:
    input: script='scripts/clustermap_data_som.py',
           median=rules.summarize_data_by_neuron.output,
           sd=rules.summarize_data_by_neuron.output
    params: cluster_x=config['cluster_x']
    output: euclidean='results/{analyses}/som_{n}/clustermap_data_som.png',
            sd='results/{analyses}/som_{n}/clustermap_data_som_sd.png',
            cosine='results/{analyses}/som_{n}/clustermap_data_som_cosine.png' 
    script: 'scripts/clustermap_data_som.py'

rule merge_clusters:
    input: script='scripts/merge_clusters.py',
           som=rules.merge_distance.output,
           som_clustered=rules.cluster_summarized_data.output,
           data=rules.cluster_raw_data.output
    output: 'results/{analyses}/som_{n}/clusters_merged.csv'
    script: 'scripts/merge_clusters.py'

rule pca_scatterplot_clusters_som:
    input: script='scripts/pca_scatterplot_clusters_som.py',
           pc=rules.pca_raw_data.output,
           clusters=rules.merge_clusters.output
    output: 'results/{analyses}/som_{n}/pca_scatterplot_clusters_som.png'
    script: 'scripts/pca_scatterplot_clusters_som.py'

rule flow_som:
    input: script='scripts/flow_som.R',
           data=lambda w: config['input_files'][w.analyses]
    output: results='results/{analyses}/som_{n}/flow_som_results.csv',
            median='results/{analyses}/som_{n}/flow_som_feature_medians_nodes.csv',
            sd='results/{analyses}/som_{n}/flow_som_feature_std_nodes.csv',
            plot='results/{analyses}/som_{n}/flow_som_mst.png',
            coord='results/{analyses}/som_{n}/mst_coordinates.csv',
            size='results/{analyses}/som_{n}/mst_node_sizes.csv',
            graph='results/{analyses}/som_{n}/mst_graph.txt'
    script: 'scripts/flow_som.R'

rule compare_som_flowsom:
    input: script='scripts/compare_som_flowsom.py',
           som=rules.merge_clusters.output,
           flowsom=rules.flow_som.output
    output: neurons='results/{analyses}/som_{n}/confusion_matrix_neurons.png'
    script: 'scripts/compare_som_flowsom.py'

rule mst_visualization:
    input: script='scripts/mst_visualization.py',
           coord=rules.flow_som.output,
           sizes=rules.flow_som.output,
           graph=rules.flow_som.output,
    output: 'results/{analyses}/som_{n}/mst_visualization.png'
    script: 'scripts/mst_visualization.py'

rule add_n_subjects_flowsom:
    input: script='scripts/add_n_subjects_flowsom.py',
           som_id=rules.flow_som.output,
           median=rules.flow_som.output,
           sd=rules.flow_som.output
    output: median='results/{analyses}/som_{n}/flow_som_feature_medians_nodes_n_subjects.csv',
            sd='results/{analyses}/som_{n}/flow_som_feature_sd_nodes_n_subjects.csv'
    script: 'scripts/add_n_subjects_flowsom.py'

rule clustermap_data_flowsom:
    input: script='scripts/clustermap_data_som.py',
           median=rules.add_n_subjects_flowsom.output,
           sd=rules.add_n_subjects_flowsom.output,
    params: cluster_x=config['cluster_x']
    output: euclidean='results/{analyses}/som_{n}/clustermap_data_flowsom.png',
            sd='results/{analyses}/som_{n}/clustermap_data_flowsom_sd.png',
            cosine='results/{analyses}/som_{n}/clustermap_data_flowsom_cosine.png' 
    script: 'scripts/clustermap_data_som.py'

rule plot_outliers_change_from_baseline:
    input: script='scripts/plot_outliers_change_from_baseline.py',
           data='results/change_from_baseline/som_25/clusters_merged.csv'
    output: 'results/change_from_baseline/som_25/plot_outliers_change_from_baseline.png'
    script: 'scripts/plot_outliers_change_from_baseline.py'

rule plot_decline_change_from_baseline:
    input: script='scripts/plot_decline_change_from_baseline.py',
           data='results/change_from_baseline/som_25/clusters_merged.csv'
    output: 'results/change_from_baseline/som_25/plot_decline_change_from_baseline.png'
    script: 'scripts/plot_decline_change_from_baseline.py'

rule plot_improvement_change_from_baseline:
    input: script='scripts/plot_improvement_change_from_baseline.py',
           data='results/change_from_baseline/som_25/clusters_merged.csv'
    output: 'results/change_from_baseline/som_25/plot_improvement_change_from_baseline.png'
    script: 'scripts/plot_improvement_change_from_baseline.py'

rule plot_flowsom_nodes_containing_som_outliers:
    input: script='scripts/plot_flowsom_nodes_containing_som_outliers.py',
           data='results/change_from_baseline/som_25/clusters_merged.csv',
           flowsom='results/change_from_baseline/som_25/flow_som_results.csv'
    output: 'results/change_from_baseline/som_25/plot_flowsom_nodes_containing_som_outliers.png'
    script: 'scripts/plot_flowsom_nodes_containing_som_outliers.py'

rule define_meta_clusters:
    input: script='scripts/define_meta_clusters.py',
    output: 'results/meta_clusters.csv'
    script: 'scripts/define_meta_clusters.py'

rule merge_metaclusters_subjects:
    input: script='scripts/merge_metaclusters_subjects.py',
    output: 'results/som_results.csv'
    script: 'scripts/merge_metaclusters_subjects.py'

rule confusion_matrix_meta_clusters:
    input: script='scripts/confusion_matrix_meta_clusters.py',
           data=rules.merge_metaclusters_subjects.output,
    output: pf='results/confusion_matrix_pysom_flowsom.png',
            cs='results/confusion_matrix_change_slopes.png',
    script: 'scripts/confusion_matrix_meta_clusters.py'


rule all:
    input: expand('results/{analyses}/raw_data/histogram.png', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/heatmap_correlations.png', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/pca_principal_components.csv', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/dendrogram.png', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/clustermap.png', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/box_plot_clusters.png', analyses=config['analyses']),
           expand('results/{analyses}/raw_data/pca_scatterplot_clusters.png', analyses=config['analyses']),
           expand('results/{analyses}/som_{n}/som_id.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/neuron_distance.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/som_id_distance.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/SOM_point_map.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/data_by_neuron_median.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/clustering_som_dendrogram.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/box_plot_clusters_som.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/clustermap_data_som.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/clusters_merged.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/pca_scatterplot_clusters_som.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/flow_som_results.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/confusion_matrix_neurons.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/mst_visualization.png', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/flow_som_feature_medians_nodes_n_subjects.csv', analyses=config['analyses'], n=config['n_neurons']),
           expand('results/{analyses}/som_{n}/clustermap_data_flowsom.png', analyses=config['analyses'], n=config['n_neurons']),
           rules.plot_outliers_change_from_baseline.output,
           rules.plot_decline_change_from_baseline.output,
           rules.plot_improvement_change_from_baseline.output,
           rules.plot_flowsom_nodes_containing_som_outliers.output,
           rules.define_meta_clusters.output,
           rules.merge_metaclusters_subjects.output,
           rules.confusion_matrix_meta_clusters.output,

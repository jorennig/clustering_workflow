# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("FlowSOM")

vignette("FlowSOM")

library(flowCore)
library(FlowSOM)
library(ggplot2)

input_data <- snakemake@input[["data"]]
output_results <- snakemake@output[["results"]]
output_median <- snakemake@output[["median"]]
output_sd <- snakemake@output[["sd"]]
output_plot <- snakemake@output[["plot"]]
output_coord <- snakemake@output[["coord"]]
output_size <- snakemake@output[["size"]]
output_graph <- snakemake@output[["graph"]]
n_neurons = strtoi(snakemake@wildcards[["n"]], base=0L)

x_nodes <- sqrt(n_neurons)
y_nodes <- sqrt(n_neurons)
max_meta <- 8

# rlen : the number of times it loops over training data
# distf : Distance function (1=manhattan, 2=euclidean, 3=chebyshev, 4=cosine)
# init: default FALSE = initialized random
# initf: if init=TRUE, set init function: Initialize_PCA or Initialize_KWSP

rawData <- read.csv(input_data, header=T)

# convert to numeric matrix, col1 = subject_id
matData <- data.matrix(rawData[c(2:ncol(rawData))])

# set seed for reproducibility
s = 2610
set.seed(s)

# create flowFrame object (required input format for FlowSOM)
data_flowFrame <- flowCore::flowFrame(matData)

# no pre-processing necessary
som_obj <- ReadInput(data_flowFrame, compensate = FALSE, 
                     transform = FALSE, scale = FALSE)

som_obj <- BuildSOM(som_obj, 
                    xdim=x_nodes, ydim=y_nodes,
                    rlen=10,
                    distf=2,
                    init = TRUE,
                    initf = Initialize_PCA) 
som_obj <- BuildMST(som_obj)

# Apply metaclustering
metacl <- MetaClustering(som_obj$map$codes,
                         "metaClustering_consensus", 
                         max = max_meta,
                         seed = s)

# get som id and meta cluster id
results <- data.frame(subject_id=rawData[1],
                      som_id=som_obj$map$mapping[, 1],
                      meta_cluster_id=metacl[som_obj$map$mapping[, 1]])

write.csv(results, output_results, col.names=T, row.names=F, sep=",")
write.csv(som_obj$map$medianValues, output_median, col.names=T, row.names=T, sep=",")
write.csv(som_obj$map$sdValues, output_sd, col.names=T, row.names=T, sep=",")

# get x and y coordinates of nodes
write.csv(som_obj$MST$l, output_coord, col.names=T, row.names=T)
# get sizes of nodes
write.csv(som_obj$MST$size, output_size, col.names=T, row.names=T)
# contains iGraph attributes --> the connectors for the MST
write_graph(som_obj$MST$graph, output_graph, "edgelist")

png(output_plot, width=5000, height=4000)
PlotStars(som_obj, view='MST')
dev.off()

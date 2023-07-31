# Supplementary File

## 1. Installation

Install the latest version from GitHub
`````{R}
devtools::install_github("behnam-yousefi/PLEXI/PLEXI")
library(PLEXI)
`````
PLEXI will also install Keras for R, which need to be activated by installation of Miniconda. For this, according to the installation guidline of Keras ([here](https://cran.r-project.org/web/packages/keras/vignettes/index.html)):
`````{R}
library(keras)
install_keras()
`````
This is required only once for the installation.

## 2. Apply on simulated networks

To test the ```PLEXI``` package, a toy example of multilayer network can be generated using the ```network_gen()``` function:
`````{R}
myNet = network_gen(n.nodes = 100, n.var.nodes = 5, n.var.nei = 90, noise.sd = .01)
`````
The process is described as the following:
1. two fully connected networks with ```n.nodes``` number of matching nodes and uniform random edge weights is generated.
2. a number of ```n.var.nodes``` nodes are randomly selected to have different edge weights with ```n.var.nei``` number of nodes (```n.var.nei``` $<$ ```n.nodes``` $- 1$) between the two networks.
3. a set of random *Gaussian* noise with zero mean and sd = ```noise.sd``` is generated and added to all of the edge weights.

The generated multiplex network and the set of the randomly selected nodes are accessible by the following lines, respectively.
`````{R}
graph_data = myNet$data_graph
var_nodes = myNet$var_nodes
`````
We then feed ```graph_data``` to the PLEXI pipeline specialized for a two-layer multiplex network (Scenario *I*), which is composed of two commands:
`````{R}
embeddingSpaceList = plexi_embedding_2layer(graph_data, train.rep = 50)
plexi_output = plexi_node_detection_2layer(embeddingSpaceList)
print(plexi_output$high_var_nodes_index)
`````
the ```plexi_embedding_2layer()``` function represents all the nodes in a common embedding space (Section 2.1); and the ```plexi_node_detection_2layer()``` function calculates the node-pair distances and asignes a p-value to each node-pair (Section 2.2). This process is repeated ```train.rep``` times to improve the robustness. The source code available at [Usage_Examples/network_generation_ex.R](Usage_Examples/network_generation_ex.R).

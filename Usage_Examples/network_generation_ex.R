rm(list=ls())

library(PLEXI)

myNet = network_gen(n.nodes = 100, n.var.nodes = 5, n.var.nei = 90, noise.sd = .01)

graph_data = myNet$data_graph
var_nodes = myNet$var_nodes

embeddingSpaceList = plexi_embedding_2layer(graph_data, train.rep = 50)
plexi_output = plexi_node_detection_2layer(embeddingSpaceList)

print(plexi_output$high_var_nodes_index)
print(var_nodes)
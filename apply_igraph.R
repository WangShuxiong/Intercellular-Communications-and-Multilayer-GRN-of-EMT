library(igraph)
# input the cluster-cluster interactions
adjmatrix = matrix(c(0.301375268249383,0.367203967487594,0.501198574002041,0.337190490687900,0.597227765214592,0.734136424178725,1,0.674353519816110,0.501396446305173,0.609047946619900,0.831535121802832,0.559384718842893,0.569614810342009,0.701419574895473,0.955048571411523,0.644385517821623),nrow=4,ncol=4,byrow = FALSE)
#Create graphs from adjacency matrices
g=graph_from_adjacency_matrix(adjmatrix, mode = "directed", weighted = TRUE,diag = TRUE, add.colnames = NULL, add.rownames = NA)
strength(g,mode = "in")
strength(g,mode = "out")
closeness(g, vids = V(g), mode = "in", normalized = TRUE)
closeness(g, vids = V(g), mode = "out", normalized = TRUE)
page_rank(g, algo = "prpack",vids = V(g), directed = TRUE, damping = 0.85)$vector
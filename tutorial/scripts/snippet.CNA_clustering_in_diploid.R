#// always use very fine resolution

n_neighbors = 10
k_superclones <- 10; k_subclones <- 10

findClusters_embedding <- 'umap'
findClusters_method <- 'hdbscan'
findClusters_ncomponents <- 2

message('n_neighbors = ', n_neighbors)
message('k_superclones = ', k_superclones)
message('k_subclones = ', k_subclones)

objd <- copykit::runPca(objd)
# umap
objd <- runUmap(objd, n_neighbors=n_neighbors, n_threads=20, min_dist = 0)

objd  <- findClusters(
  objd, 
  k_superclones = k_superclones,
  k_subclones = k_subclones,
  embedding = findClusters_embedding, 
  method = findClusters_method, 
  ncomponents = findClusters_ncomponents)
# consensu phylo
objd <- calcConsensus(objd, consensus_by = 'subclones')
objd <- runConsensusPhylo(objd)
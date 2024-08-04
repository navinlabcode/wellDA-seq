message('Snippet: Clustering CNAs of aneuploid cells')
print(objd)

if (ncol(objd) < 200) {
  n_neighbors = 3
  k_superclones <- 3; k_subclones <- 3
} else if (ncol(objd) < 500) {
  n_neighbors = 5
  k_superclones <- 5; k_subclones <- 5
} else if (ncol(objd) < 1000) {
  n_neighbors = 10
  k_superclones <- 10; k_subclones <- 10
} else {
  n_neighbors = 20
  k_superclones <- 20; k_subclones <- 20
}

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

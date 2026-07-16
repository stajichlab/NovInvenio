# Cluster candidates before TBLASTN validation

Candidates are clustered with mmseqs2 (`easy-cluster --min-seq-id 0.3 -c 0.8`) before genomic validation. Only cluster representatives are searched against outgroup genomes via TBLASTN; hits are expanded to all cluster members.

## Trade-off

Clustering simplifies the overall pattern of genes and builds consensus across similar proteins (including paralogs). It reduces the number of TBLASTN searches from one-per-candidate to one-per-cluster.

The cost: a non-representative cluster member may have an outgroup genome hit that the representative does not, leaving it undetected. At 30% identity / 80% coverage, members of the same cluster are similar but not identical — the representative is a consensus proxy, not a guarantee.

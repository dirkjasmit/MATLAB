```markdown
# MATLAB
MATLAB tools and functions

This repository contains some useful functions for MATLAB

```

```markdown
# 1. clusterPermutationTest3D_par.m is a script to perform cluster permutation test for between or within subject t-test designs. Support file GetConnectMask provides an adjacency matrix for EEGLAB structs. Cluster permutation according to Maris and Oostenveld (2007).

function [p_values, observedClusters, observedStats, permutedStats, clusterMask, stats, scores] = clusterPermutationTest3D_between(D, G, adjacencyMatrix, varargin)

D            : per subject ERPs (channel x time x subject) matrix, either values or difference scores between conditions.
G            : either
- a grouping variable. Only two values allowed for now (0 and 1).
- a data struct with the same dimensions as D that holds values of the second level of the within-subject condition
- empty [], to be used when D already holds the difference scores between the conditions to-be-tested.
adjacencyMatrix : matrix describing which channels are connected (numChannels x numChannels). Take a look at GetConnectMask(), a handy tool to get an adjacency mask for EEGLAB data.

**optional input (format "key", value )**

- **"covariates"** : covariates. UNDER CONSTRUCTION alpha version is ready: t-statiostics are calulcated 
- **"numPermutations"** : how many permutations to perform
- **"pThreshold"** : threshold at which clusters are defined (def=0.01)
- **"direction"** : 'pos' or 'neg': positive or negative effects (def = 'pos')
- **"numClusters"** : number of clusters to extract (def = 1)

**return values:**

- p_values: one sided pvalues for the clusters
- observedClusters: indices of the clusters
- observedStats: summed-t statistic for the cluster
- permutedStats: summed-t statistics of the randmomizations
- clusterMask: 2D matrix (channel x time) of the clusters (1, 2, ...)
- stats : further statistical output, including the 2D (channel x time) t statistics



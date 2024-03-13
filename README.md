# COMSE: Feature selection

*Mar 13, 2024 with COMSE version 1.0.0*

**Abstract**

COMSE is an unsupervised feature selection framework using community detection to capture informative genes from scRNA-seq data. COMSE identified homogenous cell substates with high resolution, as demonstrated by distinguishing different cell cycle stages. Evaluations based on real and simulated scRNA-seq datasets showed COMSE outperformed methods even with high dropout rates in cell clustering assignment. We also demonstrate that by identifying communities of genes associated with batch effects, COMSE parses signals reflecting biological difference from noise arising due to differences in sequencing protocols, thereby enabling integrated analysis of scRNA-seq datasets generated differently.

**Package**
* R:  igraph, Rdimtools


## Contents

[toc]



## 1 Introduction to COMSE

COMSE first partitions all genes into different communities in latent space inferred by Principle Component Analysis (PCA) using the Louvain algorithm[1]. Within each community, we apply a denoising procedure to remove noise introduced during sequencing or other procedures. It then selects highly informative genes from each community based on the Laplacian score [2] (Fig. 1). For more information on COMSE, we recommend the user to check the following article:

> COMSE: Analysis of Single-Cell RNA-seq Data Using Community Detection Based Feature Selection (https://doi.org/10.1101/2023.06.03.543526)

Please cite this article if you use SIGNET in your research. 

## 2 Requirements

### 2.1 Package installation

For R platform, the core dependent packages are `igraph` and `Rdimtools`. 

``` R
install.packages("igraph")
install.packages("Rdimtools")
```


### 2.2 Input: expression matrix

The input of COMSE is the **expression matrix** of scRNA-seq:

* Each column represents a cell sample and each row represents a gene. 
* The row name of the matrix should be the gene-symbol of gene ID.
* Expression units: The preferred expression values are raw values. Since we will use the binarized data matrix for MLP training, the expression units (raw, TPM, FPKM/RPKM) have almost no effect on the binarization matrix.

## 3 COMSE Running demo
``` R
library(igraph)
library(Rdimtools)
source("Functions_COMSE_Denoise.R")
## Loading data
data = readRDS("demo_data.RDS")
counts = data[[1]]
meta = data[[2]]
result = COMSE(counts)
```



## Reference

[1] Blondel,V.D., Guillaume,J.-L., Lambiotte,R. and Lefebvre,E. (2008) Fast unfolding of communities in large networks. Journal of Statistical Mechanics: Theory and Experiment, 10.1088/1742-5468/2008/10/P10008.

[2] He,X., Cai,D. and Niyogi,P. (2005) Laplacian Score for Feature Selection. In 18th International Conference on Neural Information Processing Systems.pp. 507–514.

[3] Qinhuan Luo, Yao Chen, and Xun Lan. Comse: Analysis of single-cell rna-seq data using community detection based feature selection. bioRxiv, 2023.

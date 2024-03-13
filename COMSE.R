library(igraph)
library(Rdimtools)

# Test data                  
data = readRDS("~/Desktop/demo_data.RDS")
counts = data[[1]]
meta = data[[2]]
result = COMSE(counts,feature_gene_num = 2000)

## selected gene : result[[1]]





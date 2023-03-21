library(igraph)
library(Rdimtools)
# Function
BuildSimKNN <- function(data, k){
  ## data_normalization, intra-samples
  dist_matrix = as.matrix(dist(data))
  neighbours <- apply(dist_matrix,1,function(x) sort(x,index.return=TRUE)$ix[2:(k+1)])  %>% as.integer
  adj <- as.matrix(Matrix::sparseMatrix(i=rep(1:nrow(data),each=k),
                                        j=neighbours,
                                        x=1,
                                        dims=c(nrow(data),nrow(data))))
  adj <- (adj | t(adj)) * 1
  #sigma = mean(dist_matrix[adj!=0])
  #print(sigma)
  #dist_matrix = exp(-dist_matrix*dist_matrix/(2*sigma*sigma))
  #dist_matrix[adj!=1] = 0
  return(adj)
}
RMLaplacianScore <- function(data,k=30,t,slices = 100,sample_size=50){
  if(length(which(colSums(data)==0))){
    data = data[,-which(colSums(data)==0)]
  }
  sample_num = min(sample_size,dim(data)[2])
  if(sample_num == dim(data)[2]){
    return(do.lscore(t(data),type = c("knn",k), t= t)$lscore)
  }
  score_total = matrix(0,nrow = slices, ncol = dim(data)[1])
  for(j in 1:slices){
    sample_index = sample(1:dim(data)[2],sample_num)
    score_total[j,] = do.lscore(t(data[,sample_index]),type = c("knn",k), t= t)$lscore
  }
  score_total[is.nan(score_total)]<-10
  score = apply(score_total,2,mean)
  return(score)
}
sub_gene_calculation<-function(group,num_gene){
  gene_num = as.numeric(group)
  total_gene_num = sum(gene_num)
  small = (gene_num*num_gene/total_gene_num)-floor(gene_num*num_gene/total_gene_num)
  small = small>=median(small)
  gene_subnum = floor(gene_num*num_gene/total_gene_num) + small
  return(gene_subnum)
}
feature_gene_shrink <- function(group,target_num,feature_gene,weight=1){
  ratio = group/length(feature_gene) * weight
  ratio = ratio/sum(ratio)
  gene_subnum = round(ratio * target_num)
  gene_shrinked = feature_gene[1:gene_subnum[1]]
  for(i in 2:length(gene_subnum)){
    if(gene_subnum[i]>0){
      gene_shrinked = c(gene_shrinked,feature_gene[(sum(group[1:(i-1)])+1):(sum(group[1:(i-1)])+gene_subnum[i])])
    }
  }
  return(gene_shrinked)
}
sort_feature <- function(genes, relation,threshold){
  genes_sort = character()
  genes_sort[1] = genes[1]
  for(i in 2:(length(genes)-1)){
    #genes_candidate = rowMeans(as.data.frame(relation[-which(rownames(relation)%in%genes_sort),genes_sort]))
    genes_candidate = apply(as.data.frame(relation[-which(rownames(relation)%in%genes_sort),genes_sort]),1,max)
    genes_candidate2 = names(genes_candidate)[which(genes_candidate<=threshold)]
    if(length(genes_candidate2)){
      genes_add = genes[min(which(genes%in%genes_candidate2))]
      genes_sort = c(genes_sort,genes_add)
    }
    else{
      genes_add = genes[min(which(genes%in%names(genes_candidate)))]
      genes_sort = c(genes_sort,genes_add)
    }
  }
  genes_sort[length(genes)] = genes[which(!genes%in%genes_sort)]
  return(genes_sort)
}
Partition <- function(counts,k_neighbors = 30, min.cells = 3){
  data = counts[which(rowSums(counts!=0)>=min.cells),]
  data = log1p(sweep(data,1,rowSums(data),FUN = "/")*1e4)
  pca = prcomp(data, scale=TRUE)
  print("PCA Completed")
  pca.var <- pca$sdev^2 
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
  data_pca = pca$x[,1:5]
  graph = BuildSimKNN(data_pca,k = k_neighbors)
  graph[graph!=0] = 1
  G = graph.adjacency(graph,mode = "undirected")
  print("Graph Construction")
  fc <-multilevel.community(G,weights = E(G)$weight)
  mem = fc$membership
  names(mem) = rownames(data_pca)
  return(mem)
}
COMSE<-function(counts,group_min = 30 , k_neighbors = 30, min.cells = 3, feature_gene_num = 2000, t= 200, q= 0.95,ssize = 100, sslices = 100){
  mem = Partition(counts)
  names(mem) = rownames(counts)
  group = table(mem)[which(table(mem)>=group_min)]
  data = log1p(sweep(counts,2,colSums(counts),FUN = "/")*1e4)
  gene_feature = character()
  #weight = numeric()
  for(i in 1:length(group)){
    print(i)
    data_sub = data[which(mem==names(group)[i]),]
    score = RMLaplacianScore(data_sub,t=t,slices=sslices,sample_size = ssize)
    names(score) = rownames(data)[which(mem==names(group)[i])]
    genes = names(sort(score))
    system_uncertainty = cor_matrix[genes,genes]
    system_total = system_uncertainty[upper.tri(system_uncertainty)]
    threshold = max(0.5,quantile(system_total, q))
    #weight[i] = quantile(system_total, q)
    print(threshold)
    genes_sorted = sort_feature(genes,system_uncertainty,threshold)
    #which(genes_sorted%in%gene1)
    gene_feature = c(gene_feature,genes_sorted)
  }
  selected_gene =  gsub("_","-",feature_gene_shrink(table(mem)[which(table(mem)>group_min)],feature_gene_num,gene_feature))
  result = list(HIG = selected_gene, membership = mem, sorted_gene = gene_feature)
  return(result)
}


## demo data
load("~/Desktop/Projects/project_GRN/data_brain/sceMouseBrain.RData")
counts<-sceMouseBrain@assays$data$counts
result = COMSE(counts)




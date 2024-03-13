## Core functions
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
GenePartition <- function(data,k_neighbors = 30){
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
BuildSimKNN2 <- function(data, k){
  ## data_normalization, intra-samples
  dist_matrix = as.matrix(dist(data))
  neighbours <- apply(dist_matrix,1,function(x) sort(x,index.return=TRUE)$ix[2:(k+1)])  %>% as.integer
  adj <- as.matrix(Matrix::sparseMatrix(i=rep(1:nrow(data),each=k),
                                        j=neighbours,
                                        x=1,
                                        dims=c(nrow(data),nrow(data))))
  #sigma = mean(dist_matrix[adj!=0])
  #print(sigma)
  #dist_matrix = exp(-dist_matrix*dist_matrix/(2*sigma*sigma))
  #dist_matrix[adj!=1] = 0
  return(adj)
}
CellPartition <- function(data,k_neighbors = 30){
   pca = prcomp(t(data), scale=TRUE)
   print("PCA Completed")
   pca.var <- pca$sdev^2 
   pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
   data_pca = pca$x[,1:10]
   #data_umap=umap(data_pca)
   graph = BuildSimKNN2(data_pca,k = max((k_neighbors-1),1))
   diag(graph) = 1
   mean_exp = as.matrix(data)%*%t(graph)/k_neighbors
   res = data - mean_exp
   return(res)
 }
DenoiseRandom <-function(counts,meta = NA,group_min = 30 ,k_neighbors_gene = 30,k_neighbors_cell = 15,CellMethod="KNN", cor_threshold = 0.5,feature_gene_num = 2000, t= 200, q= 0.95,ssize = 100, sslices = 100,min.cells = 3, Denoise = FALSE,is.log =FALSE){
  if(is.log){
    data = counts[which(rowSums(counts!=0)>=min.cells),]
    logdata = sweep(data,2,colSums(data),FUN = "/") + 1
    logdata2 = sweep(data,1,rowSums(data),FUN = "/") + 1
  }else{
    data = counts[which(rowSums(counts!=0)>=min.cells),]
    logdata = log1p(sweep(data,2,colSums(data),FUN = "/")*1e4)
    logdata2 = log1p(sweep(data,1,rowSums(data),FUN = "/")*1e4)
  }
  mem = GenePartition(logdata,k_neighbors_gene)
  group = table(mem)[which(table(mem)>=group_min)]
  print(group)
  gamma = data.frame(matrix(0,nrow = length(group),ncol = dim(data)[2]))
  data_iter = logdata2
  #res2 = data.frame(matrix(0,nrow(data),ncol(data)))
  if(Denoise){
    res1 = data.frame(matrix(0,nrow(data),ncol(data)))
    if(!is.na(meta)){
      print(meta)
      for(i in 1:dim(data)[1]){
        fit = bayesglm(as.numeric(data_iter[i,]) ~  meta,prior.mean = 0, prior.scale = 1, prior.df = Inf)
        res1[i,] = fit$residuals
        counter(i,1000)
      }
    }
    else{
      if(CellMethod =="KNN"){
        res1 = CellPartition(logdata2,k_neighbors_cell)
      }else{
        alpha = 50/k_neighbors_cell
        seed =  1234
        lda = LDA(t(data),control = list(alpha = alpha, seed = seed), k = k_neighbors_cell)
        group_info = lda@gamma
        for(i in 1:dim(data)[1]){
          fit = lm(as.numeric(data_iter[i,]) ~  group_info)
          res1[i,] = fit$residuals
          counter(i,1000)
        }
      }
    }
    print("step2 started!")
    for(k in 1:length(group)){
      data_inter = as.data.frame(res1[which(mem==names(group)[k]),])
      colnames(data_inter) = colnames(data)
      data_inter$gene = rownames(data_inter)
      data_inter2 = melt(data_inter)
      data_inter = data_inter[,1:(dim(data_inter)[2]-1)]
      fit <- lme(value~1,random = ~1|variable,data=data_inter2,control = lmeControl(opt = "optim"))
      gamma_eval = fit$coefficients$random$variable[,1] + fit$coefficients$fixed
      data_iter[which(mem==names(group)[k]),] = t(t(data_iter[which(mem==names(group)[k]),]) - gamma_eval)
      gamma[k,] =  gamma_eval
      print(k)
    }
    
  }
  print("Feature Gene Selection Started!")
  cor_matrix_list = list()
  for(i in 1:length(group)){
    data_sub = logdata[which(mem==names(group)[i]),]
    cor_matrix_list[[i]] = cor(t(data_sub))
  }
  print("Similarity Matrix Calculated!")
  gene_feature = character()
  #weight = numeric()
  print("Gene Rank start calculation!")
  
  for(i in 1:length(group)){
    print(i)
    data_sub = logdata2[which(mem==names(group)[i]),]
    score = RMLaplacianScore(data_sub,t=t,sample_size = ssize,slices = sslices, k =10)
    names(score) = rownames(data)[which(mem==names(group)[i])]
    genes = names(sort(score))
    system_uncertainty = cor_matrix_list[[i]]
    system_total = system_uncertainty[upper.tri(system_uncertainty)]
    threshold = max(cor_threshold,quantile(system_total, q))
    genes_sorted = sort_feature(genes,system_uncertainty,threshold)
    #which(genes_sorted%in%gene1)
    gene_feature = c(gene_feature,genes_sorted)
  }
  
  print("Gene Rank Calculated!")
  selected_gene =  gsub("_","-",feature_gene_shrink(table(mem)[which(table(mem)>group_min)],feature_gene_num,gene_feature))
  if(Denoise){
    return(list(denoise_dat = data_iter,HIG = selected_gene,membership = mem, sorted_gene = gene_feature))
  }else{
    return(list(HIG = selected_gene,membership = mem, sorted_gene = gene_feature))
  }
}
##Evaluation Functions
confusion.table <- function(mat){
  r.total <- rowSums(mat)
  c.total <- colSums(mat)
  N <- sum(r.total)
  sq.mat <- mat^2
  sq.sum <- sum(rowSums(sq.mat))
  r.sq.sum <- sum(r.total^2)
  c.sq.sum <- sum(c.total^2)
  a <- 0.5*(sq.sum-N)
  b <- 0.5*(r.sq.sum-sq.sum)
  c <- 0.5*(c.sq.sum-sq.sum)
  d <- 0.5*(sq.sum+N^2-r.sq.sum-c.sq.sum)
  c.mat <- matrix(c(a,b,c,d),nrow = 2,byrow = T)
  return(c.mat)
}
cluster.metric <- function(true,pred,method,beta=NULL){
  true <- as.factor(true)
  pred <- as.factor(pred)
  mat <- table(true,pred)
  c.mat <- confusion.table(mat)
  a <- c.mat[1,1]
  b <- c.mat[1,2]
  c <- c.mat[2,1]
  d <- c.mat[2,2]
  s <- a + b + c + d
  
  if (method == "Purity"){
    N <- length(true)
    l.pred <- length(mat[1,])
    int <- numeric()
    for (i in 1:l.pred) {
      int[i] <- max(mat[,i])
    }
    score <- sum(int)/N
  }
  else if (method == "F"){
    precision <- a/(a+c)
    recall <- a/(a+b)
    score <- (beta^2+1)*precision*recall/(beta^2*precision+recall)
  }
  else if (method == "Rand"){
    score <- (a+d)/(a+b+c+d)
  }
  else if (method == "Jaccard"){
    score <- a/(a+b+c)
  }
  else if (method == "Fowlkes-Mallows"){
    score <- a/sqrt((a+b)*(a+c))
  }
  else if (method == "Hubert-Arabie"){
    numerator <- s*(a+d)-((a+b)*(a+c)+(c+d)*(b+d))
    denominator <- s^2 - ((a+b)*(a+c)+(c+d)*(b+d))
    score <- numerator/denominator
  }
  return(score)
}
score_function<-function(data){
  resol = seq(0.1,1,0.1)
  score_matrix=matrix(0,0,4)
  for(i in 1:length(resol)){
    data<-FindClusters(data,resolution=resol[i],verbose = FALSE)
    score_matrix = rbind(score_matrix,c(cluster.metric(data$cell,data$seurat_clusters,method = "Purity"),
                                        cluster.metric(data$cell,data$seurat_clusters,method = "F", beta = 1),
                                        cluster.metric(data$cell,data$seurat_clusters,method = "Rand"),
                                        cluster.metric(data$cell,data$seurat_clusters,method = "Hubert-Arabie")))
  }
  return(score_matrix)
}
score_function_pure<-function(true_label,predicted_label){
  score  = c(cluster.metric(true_label,predicted_label,method = "Purity"),
             cluster.metric(true_label,predicted_label,method = "F", beta = 1),
             cluster.metric(true_label,predicted_label,method = "Rand"),
             cluster.metric(true_label,predicted_label,method = "Hubert-Arabie"))
  return(score)
}
counter = function(i, c = 100){
  if(i %% c ==0){
    print(i)
  }
}
GetFeatureGene<- function(counts,feature_gene_num = 2000,min.cell = 3){
  ##M3Drop
  counts = counts[which(rowSums(counts!=0)>=min.cell),]
  if(!is.integer(counts)){
    print("Transformation")
    counts2 = NBumiConvertToInteger(counts)
  }else{
    counts2 = counts
  }
  DANB_fit <- NBumiFitModel(counts2)
  NBgene <- NBumiFeatureSelectionCombinedDrop(DANB_fit,ntop = feature_gene_num)$Gene
  norm <- M3DropConvertData(counts2, is.counts=TRUE) 
  MDgene = M3DropFeatureSelection(norm)
  MDgene = rownames(MDgene)[1:feature_gene_num]
  rm(norm)
  gc()
  scran_sce <- SingleCellExperiment(list(counts=counts))
  scran_sce <- computeSumFactors(scran_sce)
  scran_sce <- logNormCounts(scran_sce)
  dec <- modelGeneVar(scran_sce)
  hvg = getTopHVGs(dec,n=feature_gene_num)
  rm(scran_sce)
  gc()
  data <- CreateSeuratObject(counts)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data  <- FindVariableFeatures(data , selection.method = "vst", nfeatures = 2000)
  seurat_gene = VariableFeatures(object = data)
  return(list(M3Drop = MDgene, scran = hvg, seurat = seurat_gene,NBDrop = NBgene))
}
DenoiseEvaluation<-function(counts,counts_denoise,feature_gene,meta,n_dim = 10,res = 0.3){
  score = matrix(0,0,4)
  colnames(score) = c("Purity","F","Rand","ARI")
  for(i in 1:length(feature_gene)){
    suppressMessages({
      data <- CreateSeuratObject(counts)
    })
    data$celltype = meta
    data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
    #data  <- FindVariableFeatures(data , selection.method = "vst", nfeatures = 2000)
    data <- ScaleData(data, features = gsub("_","-",feature_gene[[i]]),verbose = FALSE)
    data <- RunPCA(data, features =gsub("_","-",feature_gene[[i]]),verbose = FALSE)
    data <- FindNeighbors(data, dims= 1:n_dim,verbose = FALSE)
    #score_function(data_seurat)
    data <- RunUMAP(data,dims = 1:n_dim)
    data$cell = meta
    data<-FindClusters(data,resolution=res,verbose = FALSE)
    score_inter1 = score_function_pure(data$cell,data$seurat_clusters)
    #score_inter1 = score_inter1[5,]
    suppressMessages({
      data2 <- CreateSeuratObject(counts_denoise)
    })
    data2$celltype = meta
    #data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    if(names(feature_gene)[i] == "seurat"){
      data2  <- FindVariableFeatures(data2 , selection.method = "vst", nfeatures = 2000,verbose = FALSE)
      data2 <- ScaleData(data2, features = VariableFeatures(object = data2),verbose = FALSE)
      data2 <- RunPCA(data2, features = VariableFeatures(object = data2),verbose = FALSE)
    }
    else{
      data2 <- ScaleData(data2, features = gsub("_","-",feature_gene[[i]]),verbose = FALSE)
      data2 <- RunPCA(data2, features = gsub("_","-",feature_gene[[i]]),verbose = FALSE)
    }
    data2$cell = meta
    data2 <- FindNeighbors(data2, dims= 1:n_dim,verbose = FALSE)
    data2<-FindClusters(data2,resolution=res,verbose = FALSE)
    data2 <- RunUMAP(data2,dims = 1:n_dim)
    score_inter2 = score_function_pure(data2$cell,data2$seurat_clusters)
    #score_function(data2_seurat)
    # data2$cell = meta
    # score_inter2 = score_function_pure()
    # score_inter2 = score_inter2[5,]
    score_inter = rbind(score_inter1,score_inter2)
    colnames(score_inter) = c("Purity","F","Rand","ARI")
    rownames(score_inter) = c(paste(names(feature_gene)[i],c("orig","denoise"),sep="_"))
    score = rbind(score,score_inter)
    print(names(feature_gene)[i])
  }
  return(score)
}
DenoiseEvaluation_COMSE<-function(counts,counts_denoise,feature_gene,meta,n_dim = 10,res=0.3){
  data = CreateSeuratObject(counts)
  data$celltype = meta
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  #data  <- FindVariableFeatures(data , selection.method = "vst", nfeatures = 2000)
  data <- ScaleData(data, features = gsub("_","-",feature_gene[[1]]),verbose = FALSE)
  data <- RunPCA(data, features =gsub("_","-",feature_gene[[1]]),verbose = FALSE)
  data <- FindNeighbors(data, dims= 1:n_dim,verbose = FALSE)
  #score_function(data_seurat)
  data$cell = meta
  data<-FindClusters(data,resolution=res,verbose = FALSE)
  data <- RunUMAP(data,dims = 1:n_dim)
  score_inter1 = score_function_pure(data$cell,data$seurat_clusters)
  #score_inter1 = score_function(data)
  #score_inter1 = score_inter1[5,]
  suppressMessages({
    data2 <- CreateSeuratObject(counts_denoise)
  })
  data2$celltype = meta
  data2 <- ScaleData(data2, features = gsub("_","-",feature_gene[[2]]),verbose = FALSE)
  data2 <- RunPCA(data2, features = gsub("_","-",feature_gene[[2]]),verbose = FALSE)
  data2 <- FindNeighbors(data2, dims= 1:n_dim,verbose = FALSE)
  #score_function(data2_seurat)
  data2$cell = meta
  data2<-FindClusters(data2,resolution=res,verbose = FALSE)
  data2 <- RunUMAP(data2,dims = 1:n_dim)
  score_inter2 = score_function_pure(data2$cell,data2$seurat_clusters)
  #score_inter2 = score_function(data2)
  #score_inter2 = score_inter2[5,]
  score_inter = rbind(score_inter1,score_inter2)
  colnames(score_inter) = c("Purity","F","Rand","ARI")
  rownames(score_inter) = c(paste("COMSE",c("orig","denoise"),sep="_"))
  return(score_inter)
}
InternalClusterValidation = function(data, cluster, method = "CDbw", 
                                     distance = "euclidean", p = 1, 
                                     neighbSize = 10, minkowski_p = 1, 
                                     threads = 1){
  map = function(x, mat){
    y = numeric()
    n = length(x)
    for (i in 1:n) {
      if (x[i] %in% mat[,1]){
        j = min(which(x[i]==mat[,1]))
        y[i] = mat[j,2]
      }
      else {
        y[i] = NA
      }
    }
    return(y)
  }
  
  cl = levels(as.factor(cluster))
  K = length(cl)
  cluster = map(cluster, data.frame(cl, 1:K))
  d = ncol(data)
  N = nrow(data)
  
  center = matrix(NA, nrow = K, ncol = d)
  for (i in 1:K) {
    center[i,] = colMeans(data[which(cluster == i),])
  }
  C = colMeans(data)
  
  dist = ClusterR::distance_matrix(data, method = distance, upper = T, diagonal = T, 
                                   minkowski_p = minkowski_p, threads = threads)
  dist_c = ClusterR::distance_matrix(center, method = distance, upper = T, diagonal = T, 
                                     minkowski_p = minkowski_p, threads = threads)
  
  data_cen = data
  for (j in 1:N) {
    data_cen[j, ] = center[cluster[j], ]
  }
  
  dist_cen = ClusterR::distance_matrix(data_cen, method = distance, upper = T, diagonal = T, 
                                       minkowski_p = minkowski_p, threads = threads)
  
  sigma = function(data){
    s = colSums((data - colMeans(data))^2)
    return(s)
  }
  centerdist = function(data, distance){
    data = rbind(data, colMeans(data))
    dist = ClusterR::distance_matrix(data, method = distance, upper = T, diagonal = T, 
                                     minkowski_p = minkowski_p, threads = threads)
    return(dist[nrow(dist), ncol(dist) - 1])
  }
  
  # icv for compactness & separatness
  
  if (method == "RMSSTD" || method == "all"){
    Wk = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i), ]
      Wk = Wk + sum(sigma(subdata))
    }
    rmsstd = Wk / sqrt(N*(N-K))
    rmsstd = c(rmsstd, rmsstd, NA, NA)
  }
  else{
    rmsstd = c(NA, NA, NA, NA)
  }
  if (method == "RS" || method == "all"){
    Wk = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i), ]
      Wk = Wk + sum(sigma(subdata))
    }
    W = sum(sigma(data))
    rs = 1 - Wk / W
    rs = c(rs, NA, rs, NA) 
  }
  else{
    rs = c(NA, NA, NA, NA)
  }
  if (method == "MHG" || method == "all"){
    mhg = sum(dist * dist_cen) * (2/(N*(N-1)))
    mhg = c(mhg, NA, mhg, NA) 
  }
  else{
    mhg = c(NA, NA, NA, NA)
  }
  if (method == "CH" || method == "all"){
    #icv = fpc::calinhara(data, cluster)
    Wk = 0
    for (i in 1:K) {
      subdata = rbind(data[which(cluster == i), ], center[i, ])
      dist_sub = ClusterR::distance_matrix(subdata, method = distance, upper = T, diagonal = T, 
                                           minkowski_p = minkowski_p, threads = threads)
      Wk = Wk + sum((dist_sub[nrow(dist_sub), ])^2)
    }
    Wk = Wk / (N-K)
    data_cen_c = rbind(data_cen, C)
    dist_cen_c = ClusterR::distance_matrix(data_cen_c, method = distance, upper = T, diagonal = T, 
                                           minkowski_p = minkowski_p, threads = threads)
    Bk = sum((dist_cen_c[nrow(dist_cen_c), ])^2) / (K-1)
    ch = Bk / Wk
    ch = c(ch, Wk, Bk, NA) 
  }
  else{
    ch = c(NA, NA, NA, NA)
  }
  if (method == "I" || method == "all"){
    WD_v = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_v = c(WD_v, WD_vi)
    }
    TD_v = centerdist(data, distance = distance)
    I_c = sum(TD_v)/sum(WD_v)
    I_s = max(dist_c)
    I = (I_c * I_s / K)^p
    I = c(I, I_c, I_s, NA) 
  }
  else{
    I = c(NA, NA, NA, NA)
  }
  if (method == "D" || method == "all"){
    D = clValid::dunn(distance = dist, clusters = cluster)
    D = c(D, NA, NA, NA) 
  }
  else{
    D = c(NA, NA, NA, NA)
  }
  if (method == "S" || method == "all"){
    S = cluster::silhouette(cluster, dist = dist)
    S_sum = summary(S)
    s = c(S_sum$avg.width, NA, NA, NA) 
  }
  else{
    s = c(NA, NA, NA, NA)
  }
  if (method == "DB" || method == "all"){
    WD_s = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_s[i] = sum(WD_vi) / length(WD_vi)
    }
    WD_s = outer(WD_s, WD_s, "+")
    WD_q = WD_s / dist_c
    WD_q_max = numeric()
    for (i in 1:K) {
      WD_q_max[i] = max(WD_q[-i,i])
    }
    db = sum(WD_q_max) / K
    db = c(db, NA, NA, NA)
  }
  else{
    db = c(NA, NA, NA, NA)
  }
  if (method == "XB" || method == "all"){
    WD_v = numeric()
    for (i in 1:K) {
      WD_vi = centerdist(data[which(cluster == i),], distance = distance)
      WD_v = c(WD_v, WD_vi)
    }
    xb_c = sum(WD_v^2)
    xb_s = numeric()
    for (i in 1:K) {
      xb_s[i] = min(dist_c[-i,i])^2
    }
    xb_s = min(xb_s)
    xb = xb_c/(N*xb_s)
    xb = c(xb, xb_c, xb_s, NA)
  }
  else{
    xb = c(NA, NA, NA, NA)
  }
  if (method == "SD" || method == "all"){
    
    dist_c_max = max(dist_c)
    dist_c_min = numeric()
    for (i in 1:K) {
      dist_c_min[i] = min(dist_c[-i,i])
    }
    dist_c_min = min(dist_c_min)
    alpha = dist_c_max / dist_c_min
    Dis = alpha * sum(1/(colSums(dist_c)))
    
    dist_max = max(dist)
    dist_min = numeric()
    for (i in 1:K) {
      dist_min[i] = min(dist[-i,i])
    }
    dist_min = min(dist_min)
    beta = dist_max / dist_min
    Dis_max = beta * sum(1/(colSums(dist)))
    
    sigma_i = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      sigma_i = sigma_i + sqrt(sum((sigma(subdata))^2))/(nrow(subdata)-1)
    }
    sigma_data = sqrt(sum((sigma(data))^2))/(nrow(data)-1)
    Scat = 1/K * sigma_i / sigma_data 
    
    sd = Dis_max * Scat + Dis
    sd = c(sd, Scat, Dis, NA) 
  }
  else{
    sd = c(NA, NA, NA, NA)
  }
  if (method == "SDbw" || method == "all"){
    
    sigma_i = 0
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      sigma_i = sigma_i + sqrt(sum((sigma(subdata))^2))/(nrow(subdata)-1)
    }
    sigma_data = sqrt(sum((sigma(data))^2))/(nrow(data)-1)
    Scat = 1/K * sigma_i / sigma_data 
    
    stdev = 1/K * sqrt(sigma_i)
    sdbw_density = function(x, y){
      d = sqrt(sum((x-y)^2))
      return(ifelse(d > stdev, 0, 1))
    }
    density_ci = numeric()
    for (i in 1:K) {
      subdata = data[which(cluster == i),]
      ni = nrow(subdata)
      density_ci[i] = 0
      for (j in 1:ni) {
        density_ci[i] = density_ci[i] + sdbw_density(subdata[j,], center[i,])
      }
    }
    density_cij = matrix(NA, nrow = K, ncol = K)
    for (i in 1:(K-1)) {
      subdata_i = data[which(cluster == i),]
      for (j in (i+1):K) {
        subdata_j = data[which(cluster == j),]
        subdata_ij = rbind(subdata_i, subdata_j)
        nij = nrow(subdata_ij)
        density_cij[i,j] = 0
        cij = (center[i,] + center[j,]) / 2
        for (l in 1:nij) {
          density_cij[i,j] = density_cij[i,j] + sdbw_density(subdata_ij[l,], cij)
        }
      }
    }
    
    Dens_bw = 0
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        Dens_bw = Dens_bw + density_cij[i,j]/max(density_ci[i], density_ci[j])
      }
    }
    Dens_bw = 2/(K*(K-1)) * Dens_bw
    
    sdbw = Scat + Dens_bw
    
    sdbw = c(sdbw, Scat, Dens_bw, NA)
  }
  else{
    sdbw = c(NA, NA, NA, NA)
  }
  if (method == "CDbw" || method == "all"){
    cdbw = fpc::cdbw(data, cluster)
    cdbw = c(cdbw$cdbw, cdbw$compactness, cdbw$sep, NA)
  }
  else{
    cdbw = c(NA, NA, NA, NA)
  }
  
  # icv for connectivity
  
  if (method == "connect" || method == "all"){
    cnnt = clValid::connectivity(distance = dist, clusters = cluster, neighbSize = neighbSize)
    cnnt = c(cnnt, NA, NA, cnnt)
  }
  else{
    cnnt = c(NA, NA, NA, NA)
  }
  
  index = data.frame(rmsstd,rs,mhg,ch,I,D,s,db,xb,sd,sdbw,cdbw,cnnt)
  colnames(index) = c("RMSSTD", "RS", "MHG", "CH", "I", "D", "S",
                      "DB", "XB", "SD", "S_dbw", "C_dbw", "Connectedness")
  rownames(index) = c("Index", "Compactness", "Separatness", "Connectedness")
  return(index)
}
# Function
Partition <- function(data,k_neighbors = 30, min.cells = 3,ncomponents = 5){
  #print("Patition begin!")
  #print(dim(data))
  pca = prcomp(data, scale=TRUE)
  #print("PCA Completed")
  pca.var <- pca$sdev^2 
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  #saveRDS(pca$x[,1:100],"data_pca.rds")
  data_pca = pca$x[,1:ncomponents]
  graph = BuildSimKNN(data_pca,k = k_neighbors)
  graph[graph!=0] = 1
  G = graph.adjacency(graph,mode = "undirected")
  #print("Graph Construction")
  fc <-multilevel.community(G,weights = E(G)$weight)
  mem = fc$membership
  names(mem) = rownames(data_pca)
  return(mem)
}
COMSE<-function(counts,group_min = 30 , k_neighbors = 30, min.cells = 3, feature_gene_num = 2000, t= 200, q= 0.95,ssize = 100, sslices = 100,ncomponents = 5){
  mem1 = mem_used()
  time1 = Sys.time()
  if(min.cells>0){
    counts = counts[which(rowSums(counts!=0)>=min.cells),] 
  }
  print(dim(counts))
  data = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
  mem = Partition(data,ncomponents = ncomponents)
  names(mem) = rownames(counts)
  group = table(mem)[which(table(mem)>=group_min)]
  print(group)
  cor_matrix_list = list()
  for(i in 1:length(group)){
    data_sub = data[which(mem==names(group)[i]),]
    cor_matrix_list[[i]] = cor(t(data_sub))
  }
  print("Similarity Matrix Calculated!")
  data = log1p(sweep(counts,2,colSums(counts),FUN = "/")*1e4)
  gene_feature = character()
  #weight = numeric()
  print("Gene Rank start calculation!")
  for(i in 1:length(group)){
    print(i)
    data_sub = data[which(mem==names(group)[i]),]
    score = RMLaplacianScore(data_sub,t=t,sample_size = ssize,slices = sslices, k =k_neighbors)
    names(score) = rownames(data)[which(mem==names(group)[i])]
    genes = names(sort(score))
    system_uncertainty = cor_matrix_list[[i]]
    system_total = system_uncertainty[upper.tri(system_uncertainty)]
    threshold = max(0.5,quantile(system_total, q))
    genes_sorted = sort_feature(genes,system_uncertainty,threshold)
    #which(genes_sorted%in%gene1)
    gene_feature = c(gene_feature,genes_sorted)
  }
  print("Gene Rank Calculated!")
  selected_gene =  gsub("_","-",feature_gene_shrink(table(mem)[which(table(mem)>group_min)],feature_gene_num,gene_feature))
  result = list(HIG = selected_gene, membership = mem, sorted_gene = gene_feature)
  mem2 = mem_used()
  time2 = Sys.time()
  print("Time Consuming:")
  print(time2-time1)
  print("Memory Consuming:")
  print(mem2-mem1)
  return(result)
}
Partition_v2 <- function(data_pca,k_neighbors = 30,ncomponents = 5){
  print("Patition begin!")
  graph = BuildSimKNN(data_pca,k = k_neighbors)
  graph[graph!=0] = 1
  G = graph.adjacency(graph,mode = "undirected")
  print("Graph Construction")
  fc <-multilevel.community(G,weights = E(G)$weight)
  mem = fc$membership
  names(mem) = rownames(data_pca)
  return(mem)
}
COMSE_step2<-function(counts,data_pca,group_min = 30 , k_neighbors = 30, min.cells = 3, feature_gene_num = 2000, t= 200, q= 0.95,ssize = 100, sslices = 100,ncomponents = 5){
  counts = counts[which(rowSums(counts!=0)>=min.cells),]
  #print(dim(counts))
  data = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
  mem = Partition_v2(data_pca,ncomponents = ncomponents)
  names(mem) = rownames(counts)
  group = table(mem)[which(table(mem)>=group_min)]
  #print(group)
  cor_matrix_list = list()
  for(i in 1:length(group)){
    data_sub = data[which(mem==names(group)[i]),]
    cor_matrix_list[[i]] = cor(t(data_sub))
  }
  #print("Similarity Matrix Calculated!")
  data = log1p(sweep(counts,2,colSums(counts),FUN = "/")*1e4)
  gene_feature = character()
  #weight = numeric()
  #print("Gene Rank start calculation!")
  for(i in 1:length(group)){
    print(i)
    data_sub = data[which(mem==names(group)[i]),]
    score = RMLaplacianScore(data_sub,t=t,sample_size = ssize,slices = sslices, k = k_neighbors)
    names(score) = rownames(data)[which(mem==names(group)[i])]
    genes = names(sort(score))
    system_uncertainty = cor_matrix_list[[i]]
    system_total = system_uncertainty[upper.tri(system_uncertainty)]
    threshold = max(0.5,quantile(system_total, q))
    genes_sorted = sort_feature(genes,system_uncertainty,threshold)
    #which(genes_sorted%in%gene1)
    gene_feature = c(gene_feature,genes_sorted)
  }
  #print("Gene Rank Calculated!")
  selected_gene =  gsub("_","-",feature_gene_shrink(table(mem)[which(table(mem)>group_min)],feature_gene_num,gene_feature))
  result = list(HIG = selected_gene, membership = mem, sorted_gene = gene_feature)
  return(result)
}

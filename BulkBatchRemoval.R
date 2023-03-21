Regression <- function(counts, regress_data = regresult, methods = c("lr","zero-inflated","tweedie"),regularization = FALSE,genes=NA){
  if(is.na(genes)){
    genes = rownames(counts)
  }
  logdata = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
  res = data.frame(matrix(0,dim(counts)[1],dim(counts)[2]))
  if(methods == "lr"){
    # fit ln model 
    # 1. please use logdata instead of data for model fitting
    # 2. please chaneg the fit formula according to the model
    # 3. counter print the index for cycles
    for (i in 1:dim(logdata)[1]) {
      if(rownames(logdata)[i]%in%genes){
        #inter = my_subset(logdata[i,])
        data_inter = cbind(y=as.numeric(logdata[i,]),regress_data)
        data_inter = as.data.frame(data_inter)
        #data_inter = as.data.frame(data_inter[inter$index,])
        fit = lm(y ~ 0+.,data=data_inter) # X[,1] + X[,2] + X[,3] + X[,4]
        res[i,] = -predict(fit,as.data.frame(regress_data)) + logdata[i,]
      }
      else{
        res[i,] = logdata[i,]
      }
      # report cycle counts
      counter(i, c = 10000)
    }
  }else if(methods == "zero-inflated"){
    for(i in 1:dim(counts)[1]){
      if(rownames(counts)[i]%in%genes){
        #thres = quantile(counts[i,],0.95)
        #index = which(counts[i,]<thres)
        #zero_prop = length(which(as.numeric(counts[i,index])==0))/length(index)
        zero_prop = length(which(as.numeric(counts[i,])==0))/dim(counts)[2]
        data_inter = cbind(y=as.numeric(counts[i,]),regress_data)
        #data_inter = as.data.frame(data_inter[index,])
        data_inter = as.data.frame(data_inter)
        if(zero_prop > 0.05){
          fit <- zeroinfl(y ~ 0+.|.,data=data_inter,dist="poisson")
        }else{
          fit<-glm(y ~ 0+.,data=data_inter,family=poisson())
        }
        # report cycle counts
        res[i,] = -predict(fit,as.data.frame(regress_data)) + counts[i,]
      }else{
        res[i,] = counts[i,]
      }
      counter(i, c = 1000)
    }
  }else if(methods == "tweedie"){
    for(i in 1:dim(counts)[1]){
      if(rownames(counts)[i]%in%genes){
        #thres = quantile(counts[i,],0.95)
        #index = which(counts[i,]<thres)
        #zero_prop = length(which(as.numeric(counts[i,index])==0))/length(index)
        zero_prop = length(which(as.numeric(counts[i,])==0))/dim(counts)[2]
        data_inter = cbind(y=as.numeric(counts[i,]),regress_data)
        #data_inter = as.data.frame(data_inter[index,])
        data_inter = as.data.frame(data_inter)
        if(zero_prop > 0.05){
          fit=glm(y~0+.,family=statmod::tweedie(var.power=(1+zero_prop),link.power = 0.5),data=data_inter)
        }else{
          fit<-glm(y ~0+ .,data=data_inter,family=poisson())
        }
        # report cycle counts
        res[i,] = -predict(fit,as.data.frame(regress_data)) + counts[i,]
      }else{
        res[i,] = counts[i,]
      }
      counter(i, c = 10000)
    }
  }
  rownames(res) = rownames(counts)
  colnames(res) = colnames(counts)
  return(res)
}
Regression_Elastic <- function(counts, regress_data = regresult,
                               methods = c("lr","tweedie"),alpha_given=0,lambda_given = 0.1,regularization = FALSE,genes=NA){
  if(is.na(genes)){
    genes = rownames(counts)
  }
  logdata = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
  print(class(logdata))
  print(class(regress_data))
  res = data.frame(matrix(0,dim(counts)[1],dim(counts)[2]))
  if(methods == "lr"){
    # fit ln model 
    # 1. please use logdata instead of data for model fitting
    # 2. please chaneg the fit formula according to the model
    # 3. counter print the index for cycles
    for (i in 1:dim(logdata)[1]) {
      if(rownames(logdata)[i]%in%genes){
        fit = glmnet(x=as.matrix(regress_data),y = t(logdata[i,]) ,lambda = lambda_given,alpha = alpha_given) # X[,1] + X[,2] + X[,3] + X[,4]
        res[i,] = -predict(fit,as.matrix(regress_data)) + logdata[i,]
      }
      else{
        res[i,] = logdata[i,]
      }
      # report cycle counts
      counter(i, c = 10000)
    }
  }else if(methods == "tweedie"){
    for(i in 1:dim(counts)[1]){
      if(rownames(counts)[i]%in%genes){
        #thres = quantile(counts[i,],0.95)
        #index = which(counts[i,]<thres)
        #zero_prop = length(which(as.numeric(counts[i,index])==0))/length(index)
        zero_prop = length(which(as.numeric(counts[i,])==0))/dim(counts)[2]
        if(zero_prop > 0.05){
          #fit=glmnet(x=regress_data[index,],y = counts[i,index],family=statmod::tweedie(var.power=(1+zero_prop),link.power=0),data=data_inter,alpha = alpha_given,lambda = lambda_given)
          fit=glmnet(x=regress_data,y = t(counts[i,]),family=statmod::tweedie(var.power=(1+zero_prop),link.power=0.5),alpha = alpha_given,lambda = lambda_given)
        }else{
          #fit<-glmnet(x=regress_data[index,],y = counts[i,index],data=data_inter,family=poisson(),alpha=alpha_given,lambda = lambda_given)
          fit<-glmnet(x=regress_data,y = t(counts[i,]),family=poisson(),alpha=alpha_given,lambda = lambda_given)
        }
        # report cycle counts
        res[i,] = -predict(fit,as.matrix(regress_data)) + counts[i,]
      }else{
        res[i,] = counts[i,]
      }
      counter(i, c = 10000)
    }
  }
  rownames(res) = rownames(counts)
  colnames(res) = colnames(counts)
  return(res)
}

Regression_Eval<-function(counts,regress_data,reg_genes = NA, meta,p.adjust = 0.01){
  #var1 = Regression(counts,regress_data,methods = "zero-inflated")
  var1 = Regression(counts,regress_data,methods = "lr")
  var2 = Regression(counts,regress_data,methods = "tweedie")
  var3 = Regression(counts,regress_data,methods = "zero-inflated")
  var4 = Regression_Elastic(counts, regress_data,methods = "lr")
  var5 = Regression_Elastic(counts, regress_data,methods = "tweedie")
  if(!is.na(reg_genes)){
    #var6 = Regression(counts,regress_data,methods = "zero-inflated",genes=reg_genes)
    #print(2)
    var6 = Regression(counts,regress_data,methods = "lr",genes=reg_genes)
    var7 = Regression(counts,regress_data,methods = "tweedie",genes=reg_genes)
    var8 = Regression(counts,regress_data,methods = "zero-inflated",genes=reg_genes)
    var9 = Regression_Elastic(counts, regress_data,methods = "lr",genes=reg_genes)
    var10 = Regression_Elastic(counts, regress_data,methods = "tweedie",genes=reg_genes)
  }
  gene_list = list()
  for(i in 1:length(grep("var",ls()))){
    variable_name = paste("var",i,sep="")
    gene_list[[i]] = mylimma(get(variable_name), meta,p.adj = p.adjust) 
  }
  gene_list[[11]] = mylimma(counts, meta,p.adj = p.adjust)
  #inter = list(var1,var2,var3,var4,counts)
  #names(inter) = c("lr","tweedie","lrr","tweedier","lrv","tweediev","lrrv","tweedierv","orig")
  #names(inter) = c("lr","tweedie","lrr","tweedier","orig")
  #saveRDS(inter,paste("~/hnsc_sub/subsample_res",index,".RDS",sep = ""))
  #return(inter)
  names(gene_list) = c("lr","tweedie","zinb","lrr","tweedier","lrv","tweediev","zinbv","lrrv","tweedierv","orig")
  return(gene_list)
}
Regression_Eval2<-function(counts,regress_data,reg_genes = NA, meta,p.adjust = 0.01){
  #var1 = Regression(counts,regress_data,methods = "zero-inflated")
  var1 = Regression(counts,regress_data,methods = "lr")
  var2 = Regression_Elastic(counts, regress_data,methods = "lr")
  var3 = Regression_Elastic(counts, regress_data,methods = "lr",alpha_given=0.5)
  var4 = Regression_Elastic(counts, regress_data,methods = "lr",alpha_given=1)
  if(!is.na(reg_genes)){
    #var6 = Regression(counts,regress_data,methods = "zero-inflated",genes=reg_genes)
    #print(2)
    var5 = Regression(counts,regress_data,methods = "lr",genes=reg_genes)
    var6 = Regression_Elastic(counts, regress_data,methods = "lr",genes=reg_genes)
    var7 = Regression_Elastic(counts, regress_data,methods = "lr",alpha_given=0.5,genes=reg_genes)
    var8 = Regression_Elastic(counts, regress_data,methods = "lr",alpha_given=1,genes=reg_genes)
  }
  gene_list = list()
  for(i in 1:length(grep("var",ls()))){
    variable_name = paste("var",i,sep="")
    gene_list[[i]] = mylimma(get(variable_name), meta,p.adj = p.adjust) 
  }
  gene_list[[9]] = mylimma(counts, meta,p.adj = p.adjust)
  #inter = list(var1,var2,var3,var4,counts)
  #names(inter) = c("lr","tweedie","lrr","tweedier","lrv","tweediev","lrrv","tweedierv","orig")
  #names(inter) = c("lr","tweedie","lrr","tweedier","orig")
  #saveRDS(inter,paste("~/hnsc_sub/subsample_res",index,".RDS",sep = ""))
  #return(inter)
  names(gene_list) = c("lr","lr_ridge","lr_elastic","lr_lasso","lrv","lrv_ridge","lrv_elastic","lrv_lasso","orig")
  return(gene_list)
}
Paraller_eval<-function(counts_inter,meta,i){
  ## Filtering
  keep.exprs<-filterByExpr(counts_inter,group = meta)
  counts_inter = counts_inter[keep.exprs,]
  counts_inter = counts_inter[which(rowSums(counts_inter!=0)>=3),]
  data = log1p(sweep(counts_inter,1,rowSums(counts_inter),FUN = "/")*1e4)
  ##saveRDS(data,"~/hnsc_sub/test1_1.rds")
  #marker <- c("CD19","IL7R", "CCR7", "S100A4", "CD8A", "FCGR3A", "CD14", "CD27", "GNLY", "NKG7", "MS4A1", "LYZ")
  #result_inter = IdentImmuno(data,marker)
  #mem = result_inter[[1]]
  #idx = result_inter[[2]]
  #regresult = RegressionData(counts_inter,mem,idx)
  regresult = BulkBatchIdent(counts_inter, meta)
  #regresult2 = RegressionData_NMF(counts_inter,mem,idx,rank = 3)
  reg_num = round(dim(counts_inter)[1]/2)
  data_seurat = CreateSeuratObject(counts_inter, min.cells = 3)
  data_seurat = NormalizeData(data_seurat)
  reg_genes = VariableFeatures(FindVariableFeatures(data_seurat, nfeatures = reg_num))
  #reg_genes  = NA
  #result  = Regression_Eval(counts_inter, regresult, reg_genes, meta, index = i)
  #Regression_Eval(counts_inter, regresult, reg_genes, meta,i)
  gene_list_auc = Regression_Eval(counts_inter, as.matrix(regresult), reg_genes, meta)
  #gene_list_nmf = Regression_Eval(counts_inter, regresult2, reg_genes, hnsc$meta[sample_list[[i]]])
  saveRDS(gene_list_auc,paste("~/hnsc_exp1/DEG_List_auc_sub",i,".RDS",sep = ""))
  #saveRDS(gene_list_nmf,paste("~/hnsc_sub/DEG_List_nmf_BRCA_sub",i,".RDS",sep = ""))
  print(i)
  return(gene_list_auc)
}
Paraller_eval2<-function(counts_inter,meta,i){
  ## Filtering
  keep.exprs<-filterByExpr(counts_inter,group = meta)
  counts_inter = counts_inter[keep.exprs,]
  counts_inter = counts_inter[which(rowSums(counts_inter!=0)>=3),]
  data = log1p(sweep(counts_inter,1,rowSums(counts_inter),FUN = "/")*1e4)
  ##saveRDS(data,"~/hnsc_sub/test1_1.rds")
  #marker <- c("CD19","IL7R", "CCR7", "S100A4", "CD8A", "FCGR3A", "CD14", "CD27", "GNLY", "NKG7", "MS4A1", "LYZ")
  #result_inter = IdentImmuno(data,marker)
  #mem = result_inter[[1]]
  #idx = result_inter[[2]]
  #regresult = RegressionData(counts_inter,mem,idx)
  regresult = BulkBatchIdent(counts_inter, meta)
  #regresult2 = RegressionData_NMF(counts_inter,mem,idx,rank = 3)
  reg_num = round(dim(counts_inter)[1]/2)
  data_seurat = CreateSeuratObject(counts_inter, min.cells = 3)
  data_seurat = NormalizeData(data_seurat)
  reg_genes = VariableFeatures(FindVariableFeatures(data_seurat, nfeatures = reg_num))
  #reg_genes  = NA
  #result  = Regression_Eval(counts_inter, regresult, reg_genes, meta, index = i)
  #Regression_Eval(counts_inter, regresult, reg_genes, meta,i)
  gene_list_auc = Regression_Eval2(counts_inter, as.matrix(regresult), reg_genes, meta)
  #gene_list_nmf = Regression_Eval(counts_inter, regresult2, reg_genes, hnsc$meta[sample_list[[i]]])
  saveRDS(gene_list_auc,paste("~/hnsc_exp2/DEG_List_auc_sub",i,".RDS",sep = ""))
  #saveRDS(gene_list_nmf,paste("~/hnsc_sub/DEG_List_nmf_BRCA_sub",i,".RDS",sep = ""))
  print(i)
  return(gene_list_auc)
}
BulkBatchIdent <- function(counts,meta,k_neighbors = 30,group_min = 30, n_components=5){
  keep.exprs<-filterByExpr(counts,group = meta)
  counts = counts[keep.exprs,]
  data = counts[which(rowSums(counts!=0)>=3),]
  data = log1p(sweep(data,1,rowSums(data),FUN = "/")*1e4)
  pca = prcomp(data, scale=TRUE)
  print("PCA Completed!")
  #pca.var <- pca$sdev^2 
  #pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
  #plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
  data_pca = pca$x[,1:n_components]
  graph = BuildSimKNN(data_pca,k = k_neighbors)
  graph[graph!=0] = 1
  G = graph.adjacency(graph,mode = "undirected")
  print("Graph Construction Completed!")
  fc <-multilevel.community(G,weights = E(G)$weight)
  print("Graph Partition Completed!")
  mem = fc$membership
  names(mem) = rownames(data_pca)
  group = table(mem)[which(table(mem)>=group_min)]
  auc_matrix = data.frame(matrix(0,dim(counts)[2],length(group)))
  cells_rankings <- AUCell_buildRankings(as.matrix(counts))
  for(i in 1:length(group)){
    genes <-rownames(counts)[which(mem==names(group)[i])]
    geneSets <- GSEABase::GeneSet(genes, setName=paste("community",names(group)[i],sep = "")) # alternative
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
    auc_matrix[,i] = as.numeric(getAUC(cells_AUC)) 
  }
  colnames(auc_matrix) = paste("auc",1:dim(auc_matrix)[2],sep="")
  
  prob = numeric()
  for(i in 1:length(group)){
    aov <- aov(auc_matrix[,i]~meta)  
    inter = summary(aov)
    prob[i] = inter[[1]]$`Pr(>F)`[1]
  }
  idx = which(prob > 0.05)
  print(group)
  print(idx)
  if(!length(idx)){
    idx = which.max(inter[[1]]$`Pr(>F)`)
  }
  return(auc_matrix[,idx])
}
hnsc <- readRDS("~/HNSC/hnsc.RDS")
counts = hnsc$count
meta = hnsc$meta
sample_list = Subsample(size = 200,data=counts,meta=hnsc$meta)

cl <- parallel::makeCluster(20)
doParallel::registerDoParallel(cl)
test= foreach(i = 1:length(sample_list),.combine =list,.packages = c("limma","pscl","MASS","AUCell","edgeR","dplyr","igraph","Seurat","glmnet"),.export = GetGlobalFunction(),.errorhandling = "pass")%dopar% Paraller_eval(counts[,sample_list[[i]]],meta[sample_list[[i]]],i)
stopCluster(cl)

cl <- parallel::makeCluster(20)
doParallel::registerDoParallel(cl)
test= foreach(i = 1:length(sample_list),.combine =list,.packages = c("limma","pscl","MASS","AUCell","edgeR","dplyr","igraph","Seurat","glmnet"),.export = GetGlobalFunction(),.errorhandling = "pass")%dopar% Paraller_eval2(counts[,sample_list[[i]]],meta[sample_list[[i]]],i)
stopCluster(cl)

path = "/home/luoqinhuan/hnsc_exp2"
file = dir(path)
gene_methods = list()
GetGene<-function(genelist){
  for(i in 1:length(genelist)){
    genelist[[i]] = rownames(genelist[[i]])
  }
  return(genelist)
}

for(i in 1:length(file)){
  file_inter = readRDS(paste(path,"/",file[i],sep = ""))
  gene_methods = c(gene_methods,GetGene(file_inter)) 
}

methods = names(table(names(gene_methods)))
jaccard_similarity <- function(string1,string2){
  return(length(intersect(string1,string2))/max(length(union(string1,string2)),1))
}

cor = numeric()
methods_name = character()
for(i in 1:length(methods)){
  gene_methods_sub = gene_methods[which(names(gene_methods)==methods[i])]
  cor_inter = numeric()
  print(length(gene_methods_sub))
  for(j in 1:(length(gene_methods_sub)-1)){
    for(k in (j+1):length(gene_methods_sub)){
      cor_inter = c(cor_inter,jaccard_similarity(gene_methods_sub[[j]],gene_methods_sub[[k]]))
    }
  }
  cor = c(cor,cor_inter)
  methods_name = c(methods_name, rep(methods[i],length(cor_inter)))
}

del_ind = which(is.na(cor))
cor = cor[-del_ind]
methods_name = methods_name[-del_ind]

data_plot = data.frame(Jaccard_Similarity = cor, methods = methods_name)


ggplot(data_plot,aes(x = methods,y=Jaccard_Similarity,fill =methods))+
  scale_fill_npg()+
  geom_boxplot()+
  theme_bw()


counts_diff = unlist(lapply(gene_methods,length))
data_plot = data.frame(counts = as.numeric(counts_diff), methods = names(counts_diff))

ggplot(data_plot,aes(x = methods,y = counts,fill=methods)) + 
  scale_fill_npg()+
  geom_boxplot()
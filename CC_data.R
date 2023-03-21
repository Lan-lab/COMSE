GO_Analysis<-function(gene){
  gene_id<-mapIds(x=org.Mm.eg.db,
                  keys = gene,
                  keytype = "SYMBOL",
                  column = "ENTREZID")
  gene_id<-na.omit(gene_id)
  
  erich.go.BP<-enrichGO(gene = gene_id,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        readable = T)
  return(erich.go.BP)
}
library(RColorBrewer)
library(gplots)
library(VennDiagram)
Gene_Relevance_Calculation<-function(gene,sets=signature_cycle){
  return(length(grep(gene,sets)))
}
counts <- readRDS("~/public_data_true_label/CC_data.rds")
meta = c(rep("G1",96),rep("G2M",96),rep("S",96))
counts = counts[which(rowSums(counts!=0)>=3),]
data = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
k_neighbors = 30
group_min = 30

scran_sce <- SingleCellExperiment(list(counts=counts))
scran_sce <- computeSumFactors(scran_sce)
scran_sce <- logNormCounts(scran_sce)
dec <- modelGeneVar(scran_sce)
hvg = getTopHVGs(dec,n=2000)

data = log1p(sweep(counts,1,rowSums(counts),FUN = "/")*1e4)
pca = prcomp(data)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
data_pca = pca$x[,1:5]
graph = BuildSimKNN(data_pca,k = k_neighbors)
graph[graph!=0] = 1
G = graph.adjacency(graph,mode = "undirected")
fc <-multilevel.community(G,weights = E(G)$weight)
mem = fc$membership
group = table(mem)[which(table(mem)>=group_min)]
saveRDS(mem,"~/CellCycle/mem.RDS")
rm(graph)
cor_matrix_list = list()
for(i in 1:length(group)){
  data_sub = data[which(mem==names(group)[i]),]
  cor_matrix_list[[i]] = cor(t(data_sub))
}

data = log1p(sweep(counts,2,colSums(counts),FUN = "/")*1e4)
t = 200
q = 0.95
gene_feature = character()
weight = numeric()
for(i in 1:length(group)){
  print(i)
  data_sub = data[which(mem==names(group)[i]),]
  
  score = RMLaplacianScore(data_sub,t=t,sample_size = 30,slices = 100, k =10)
  names(score) = rownames(data)[which(mem==names(group)[i])]
  genes = names(sort(score))
  system_uncertainty = cor_matrix_list[[i]]
  system_total = system_uncertainty[upper.tri(system_uncertainty)]
  threshold = max(0.5,quantile(system_total, q))
  print(threshold)
  genes_sorted = sort_feature(genes,system_uncertainty,threshold)
  #which(genes_sorted%in%gene1)
  gene_feature = c(gene_feature,genes_sorted)
  weight[i] = quantile(system_total, q)
}

saveRDS(gene_feature,"~/CellCycle/feature_gene.RDS")
saveRDS(weight,"~/CellCycle/weight.RDS")

gene_feature = readRDS("~/CellCycle/feature_gene.RDS")
test = feature_gene_shrink(table(mem)[which(table(mem)>group_min)],2000,gene_feature)
feature_gene_num = 2000

DANB_fit <- NBumiFitModel(counts)
NBgene <- NBumiFeatureSelectionCombinedDrop(DANB_fit,ntop = feature_gene_num)$Gene

data_seurat= CreateSeuratObject(counts)
data_seurat$cell = meta
data_seurat = NormalizeData(data_seurat)
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = dim(data_seurat)[1])
feature_gene_seurat <- VariableFeatures(data_seurat)[1:2000]
all.genes <- rownames(data_seurat)
data_seurat <- ScaleData(data_seurat, features = all.genes)

data_seurat <- RunPCA(data_seurat, features =feature_gene_seurat)
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
score <- score_function(data_seurat)
index = which.max(rowMeans(score))
data_seurat <- FindClusters(data_seurat, resolution = index * 0.1)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
score_seurat = score[index,]

dataplot_seurat = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
dataplot_seurat$cell = data_seurat$cell
dataplot_seurat$predicted = data_seurat$seurat_clusters


data_seurat <- RunPCA(data_seurat, features =gsub("_","-",NBgene))
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
score <- score_function(data_seurat)
index = which.max(rowMeans(score))
data_seurat <- FindClusters(data_seurat, resolution = index * 0.1)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
score_NB = score[index,]

dataplot_NB  = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
dataplot_NB$cell = data_seurat$cell
dataplot_NB$predicted = data_seurat$seurat_clusters



data_seurat <- RunPCA(data_seurat, features = gsub("_","-",hvg))
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
score <- score_function(data_seurat)
index = which.max(rowMeans(score))
data_seurat <- FindClusters(data_seurat, resolution = index * 0.1)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
score_scran = score[index,]

dataplot_scran  = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
dataplot_scran$cell = data_seurat$cell
dataplot_scran$predicted = data_seurat$seurat_clusters



data_seurat <- RunPCA(data_seurat, features = gsub("_","-",test))
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
score <- score_function(data_seurat)
score <- score_function(data_seurat)
index = which.max(rowMeans(score))
data_seurat <- FindClusters(data_seurat, resolution = index * 0.1)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
score_COMSE = score[index,]

dataplot_COMSE  = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
dataplot_COMSE$cell = data_seurat$cell
dataplot_COMSE$predicted = data_seurat$seurat_clusters


col5<-pal_npg( "nrc")(10)[c(1,3,6,9)]

p1_1=ggplot(dataplot_seurat, aes(UMAP_1, UMAP_2, color = cell)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


p1_2=ggplot(dataplot_seurat, aes(UMAP_1, UMAP_2, color = predicted)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="Predicted Label",title = "Seurat") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 
p1 = p1_2+p1_1

p2_1=ggplot(dataplot_NB, aes(UMAP_1, UMAP_2, color = cell)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


p2_2=ggplot(dataplot_NB, aes(UMAP_1, UMAP_2, color = predicted)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="Predicted Label",title = "NBDrop") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p2 = p2_2+p2_1


p3_1=ggplot(dataplot_scran, aes(UMAP_1, UMAP_2, color = cell)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


p3_2=ggplot(dataplot_scran, aes(UMAP_1, UMAP_2, color = predicted)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="Predicted Label",title = "Scran") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p3 = p3_2+p3_1


p4_1=ggplot(dataplot_COMSE, aes(UMAP_1, UMAP_2, color = cell)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


p4_2=ggplot(dataplot_COMSE, aes(UMAP_1, UMAP_2, color = predicted)) +
  geom_point(size =1.5, shape = 19, alpha = 0.8) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values=col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="Predicted Label",title = "COMSE") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p4 = p4_2+p4_1

ggsave("~/CellCycle/result1.jpg",p1/p2,width = 8,height = 6)
ggsave("~/CellCycle/result2.jpg",p3/p4,width = 8,height = 6)


score = rbind(score_COMSE,score_NB,score_scran,score_seurat)
rownames(score) = c("COMSE","NBDrop","Scran","Seurat")
colnames(score) = c("Purity","F score","RI","ARI") 

score_melt =  reshape2::melt(score)
colnames(score_melt)[1:2] = c("Method","Validations")

p=ggplot(score_melt, aes(x=Validations, y=value, fill = Method)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85) +
  scale_fill_npg()+
  labs(x = "", y = "Validations")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 16),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 20) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 28),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 20),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


ggsave("~/CellCycle/validations.jpg",p,width = 14,height = 12)


scran_dis = rep(0,length(group))
scran_dis[as.numeric(names(table(mem[which(rownames(data)%in%hvg)])))] = as.numeric(table(mem[which(rownames(data)%in%hvg)]))
scran_dis = scran_dis/as.numeric(group)

nb_dis = rep(0,length(group))
nb_dis[as.numeric(names(table(mem[which(rownames(data)%in%NBgene)])))] = as.numeric(table(mem[which(rownames(data)%in%NBgene)]))
nb_dis = nb_dis/as.numeric(group)

seurat_dis = rep(0,length(group))
seurat_dis[as.numeric(names(table(mem[which(rownames(data)%in%feature_gene_seurat)])))] = as.numeric(table(mem[which(rownames(data)%in%feature_gene_seurat)]))
seurat_dis = seurat_dis/as.numeric(group)

comse_dis = as.numeric(table(mem[which(rownames(data)%in%test)]))/as.numeric(group)

dis = rbind(comse_dis,nb_dis,scran_dis,seurat_dis)
rownames(dis) = c("COMSE","NBDrop","Scran","Seurat")
colnames(dis) = paste("subgraph",1:28,sep="")


signature_cycle <- readRDS("~/CellCycle/signature_cycle.RDS")

genes = rownames(data)
score_gene = sapply(genes,Gene_Relevance_Calculation)


score_mean = aggregate(score_gene,by=list(mem),mean)
data_plot_score  =  data.frame(score = score_mean$x,subgraph = as.factor(names(group)))
data_plot_score$subgraph = factor(data_plot_score$subgraph, levels=1:28)
col5<-colorRampPalette((pal_npg( "nrc")(10)))(28)
p1=ggplot(data_plot_score,aes(x = subgraph,y=score,fill=subgraph))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = col5)+
  labs(x = "Subgraph", y = "Cell Cycle Gene Relevance Score")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 14) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 10),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


rownames(dis) = c("COMSE","NBDrop","Scran","Seurat")
colnames(dis) = as.factor(1:28)
col5<-colorRampPalette((pal_npg( "nrc")(10)))(11)
dis_part = dis[,data_plot_score$subgraph[which(data_plot_score$score>mean(data_plot_score$score))]]
dis_melt = reshape2::melt(dis_part)
colnames(dis_melt)[1:2]=c("Method","subgraph")
p2=ggplot(dis_melt,aes(x = as.factor(subgraph),y=value,fill=Method))+
  geom_bar(stat = "identity", position = "dodge", width = .5,colour = "black", alpha = .85) +
  scale_fill_manual(values = col5)+
  labs(x = "Subgraph", y = "Proportion of gene in CellCycle-related Subgraph")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        #axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 14) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 10),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 


ggsave("~/CellCycle/cellcycle.jpg",p1+p2,dpi=300,width = 20,height = 10)


HVG_detection <- function(counts,feature_gene_num=2000){
  scran_sce <- SingleCellExperiment(list(counts=counts))
  scran_sce <- computeSumFactors(scran_sce)
  scran_sce <- logNormCounts(scran_sce)
  dec <- modelGeneVar(scran_sce)
  hvg = getTopHVGs(dec,n=feature_gene_num)
  return(hvg)
}
NBG_detection <- function(counts,feature_gene_num=2000){
  DANB_fit <- NBumiFitModel(counts)
  NBgene <- NBumiFeatureSelectionCombinedDrop(DANB_fit,ntop = feature_gene_num)$Gene
}
Seurat_detection <- function(counts,feature_gene_num=2000){
  data_seurat= CreateSeuratObject(counts)
  data_seurat = NormalizeData(data_seurat)
  data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = dim(data_seurat)[1])
  feature_gene_seurat <- VariableFeatures(data_seurat)[1:feature_gene_num]
  return(feature_gene_seurat)
}
## HVG compairson
counts_G1 = counts[,1:96]
counts_G2M = counts[,97:192]
counts_S = counts[,193:288]

hvg1 = HVG_detection(counts_G1[which(rowSums(counts_G1)>=3),])
hvg2 = HVG_detection(counts_G2M[which(rowSums(counts_G2M)>=3),])
hvg3 = HVG_detection(counts_S[which(rowSums(counts_S)>=3),])



p = venn.diagram(
  x = list(hvg,hvg1,hvg2,hvg3),
  category.names = c("Total","G1_only","G2M_only","S_only"),
  # filename = 'venn.png',
  filename = NULL,
  output=TRUE,
  fill = pal_npg( "nrc")(4),
  col = pal_npg( "nrc")(4),
  fontface = "bold",
  #cat.col = c("#2E8B57","orange","#483D8B",""),
  cat.fontface = "bold",
  main = "Scran",
  main.fontface = "bold",
  main.cex = 2
)
png("~/CellCycle/Venn_scran.png")
grid.draw(p)
dev.off()




NBgene1 = NBG_detection(counts_G1[which(rowSums(counts_G1)>=3),])
NBgene2 = NBG_detection(counts_G2M[which(rowSums(counts_G2M)>=3),])
NBgene3 = NBG_detection(counts_S[which(rowSums(counts_S)>=3),])



p = venn.diagram(
  x = list(NBgene,NBgene1,NBgene2,NBgene3),
  category.names = c("Total","G1_only","G2M_only","S_only"),
  # filename = 'venn.png',
  filename = NULL,
  output=TRUE,
  fill = pal_npg( "nrc")(4),
  col = pal_npg( "nrc")(4),
  fontface = "bold",
  #cat.col = c("#2E8B57","orange","#483D8B",""),
  cat.fontface = "bold",
  main = "NBDrop",
  main.fontface = "bold",
  main.cex = 2
)
png("~/CellCycle/Venn_NBDrop.png")
grid.draw(p)
dev.off()


Seuratgene= Seurat_detection(counts[which(rowSums(counts)>=3),])
Seuratgene1 = Seurat_detection(counts_G1[which(rowSums(counts_G1)>=3),])
Seuratgene2 = Seurat_detection(counts_G2M[which(rowSums(counts_G2M)>=3),])
Seuratgene3 = Seurat_detection(counts_S[which(rowSums(counts_S)>=3),])



p = venn.diagram(
  x = list(feature_gene_seurat,Seuratgene1,Seuratgene2,Seuratgene3),
  category.names = c("Total","G1_only","G2M_only","S_only"),
  # filename = 'venn.png',
  filename = NULL,
  output=TRUE,
  fill = pal_npg( "nrc")(4),
  col = pal_npg( "nrc")(4),
  fontface = "bold",
  #cat.col = c("#2E8B57","orange","#483D8B",""),
  cat.fontface = "bold",
  main = "Seurat",
  main.cex = 2,
  main.fontface = "bold"
)
png("~/CellCycle/Venn_Seurat.png")
grid.draw(p)
dev.off()

counts_G1 = counts_G1[which(rowSums(counts_G1)>=3),]
data = log1p(sweep(counts_G1,1,rowSums(counts_G1),FUN = "/")*1e4)
pca = prcomp(data)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
data_pca = pca$x[,1:2]
graph = BuildSimKNN(data_pca,k = k_neighbors)
graph[graph!=0] = 1
G = graph.adjacency(graph,mode = "undirected")
fc <-multilevel.community(G,weights = E(G)$weight)
mem = fc$membership
group = table(mem)[which(table(mem)>=group_min)]

cor_matrix_list = list()
for(i in 1:length(group)){
  data_sub = data[which(mem==names(group)[i]),]
  cor_matrix_list[[i]] = cor(t(data_sub))
}

data = log1p(sweep(counts_G1,2,colSums(counts_G1),FUN = "/")*1e4)
t = 200
q = 0.95
gene_feature = character()
weight = numeric()
for(i in 1:length(group)){
  print(i)
  data_sub = data[which(mem==names(group)[i]),]
  
  score = RMLaplacianScore(data_sub,t=t,sample_size = 30,slices = 100, k =10)
  names(score) = rownames(data)[which(mem==names(group)[i])]
  genes = names(sort(score))
  system_uncertainty = cor_matrix_list[[i]]
  system_total = system_uncertainty[upper.tri(system_uncertainty)]
  threshold = max(0.5,quantile(system_total, q))
  print(threshold)
  genes_sorted = sort_feature(genes,system_uncertainty,threshold)
  #which(genes_sorted%in%gene1)
  gene_feature = c(gene_feature,genes_sorted)
  weight[i] = quantile(system_total, q)
}
saveRDS(mem,"~/CellCycle/mem_G1.RDS")
saveRDS(gene_feature,"~/CellCycle/feature_gene_G1.RDS")

test1 = feature_gene_shrink(table(mem)[which(table(mem)>group_min)],2000,gene_feature)
saveRDS(test1,"~/CellCycle/feature_gene_G1_2000.RDS")



counts_G2M = counts_G2M[which(rowSums(counts_G2M)>=3),]
data = log1p(sweep(counts_G2M,1,rowSums(counts_G2M),FUN = "/")*1e4)
pca = prcomp(data)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
data_pca = pca$x[,1:2]
graph = BuildSimKNN(data_pca,k = k_neighbors)
graph[graph!=0] = 1
G = graph.adjacency(graph,mode = "undirected")
fc <-multilevel.community(G,weights = E(G)$weight)
mem = fc$membership
group = table(mem)[which(table(mem)>=group_min)]

cor_matrix_list = list()
for(i in 1:length(group)){
  data_sub = data[which(mem==names(group)[i]),]
  cor_matrix_list[[i]] = cor(t(data_sub))
}

data = log1p(sweep(counts_G2M,2,colSums(counts_G2M),FUN = "/")*1e4)
t = 200
q = 0.95
gene_feature = character()
weight = numeric()
for(i in 1:length(group)){
  print(i)
  data_sub = data[which(mem==names(group)[i]),]
  
  score = RMLaplacianScore(data_sub,t=t,sample_size = 30,slices = 100, k =10)
  names(score) = rownames(data)[which(mem==names(group)[i])]
  genes = names(sort(score))
  system_uncertainty = cor_matrix_list[[i]]
  system_total = system_uncertainty[upper.tri(system_uncertainty)]
  threshold = max(0.5,quantile(system_total, q))
  print(threshold)
  genes_sorted = sort_feature(genes,system_uncertainty,threshold)
  #which(genes_sorted%in%gene1)
  gene_feature = c(gene_feature,genes_sorted)
  weight[i] = quantile(system_total, q)
}
saveRDS(mem,"~/CellCycle/mem_G2M.RDS")
saveRDS(gene_feature,"~/CellCycle/feature_gene_G2M.RDS")

test2 = feature_gene_shrink(table(mem)[which(table(mem)>group_min)],2000,gene_feature)
saveRDS(test2,"~/CellCycle/feature_gene_G2M_2000.RDS")




counts_S= counts_S[which(rowSums(counts_S)>=3),]
data = log1p(sweep(counts_S,1,rowSums(counts_S),FUN = "/")*1e4)
pca = prcomp(data)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
data_pca = pca$x[,1:2]
graph = BuildSimKNN(data_pca,k = k_neighbors)
graph[graph!=0] = 1
G = graph.adjacency(graph,mode = "undirected")
fc <-multilevel.community(G,weights = E(G)$weight)
mem = fc$membership
group = table(mem)[which(table(mem)>=group_min)]

cor_matrix_list = list()
for(i in 1:length(group)){
  data_sub = data[which(mem==names(group)[i]),]
  cor_matrix_list[[i]] = cor(t(data_sub))
}

data = log1p(sweep(counts_S,2,colSums(counts_S),FUN = "/")*1e4)
t = 200
q = 0.95
gene_feature = character()
weight = numeric()
for(i in 1:length(group)){
  print(i)
  data_sub = data[which(mem==names(group)[i]),]
  
  score = RMLaplacianScore(data_sub,t=t,sample_size = 30,slices = 100, k = 10)
  names(score) = rownames(data)[which(mem==names(group)[i])]
  genes = names(sort(score))
  system_uncertainty = cor_matrix_list[[i]]
  system_total = system_uncertainty[upper.tri(system_uncertainty)]
  threshold = max(0.5,quantile(system_total, q))
  print(threshold)
  genes_sorted = sort_feature(genes,system_uncertainty,threshold)
  #which(genes_sorted%in%gene1)
  gene_feature = c(gene_feature,genes_sorted)
  weight[i] = quantile(system_total, q)
}
saveRDS(mem,"~/CellCycle/mem_S.RDS")
saveRDS(gene_feature,"~/CellCycle/feature_gene_S.RDS")

test3 = feature_gene_shrink(table(mem)[which(table(mem)>group_min)],2000,gene_feature)
saveRDS(test3,"~/CellCycle/feature_gene_S_2000.RDS")



p = venn.diagram(
  x = list(test,test1,test2,test3),
  category.names = c("Total","G1_only","G2M_only","S_only"),
  # filename = 'venn.png',
  filename = NULL,
  output=TRUE,
  fill = pal_npg( "nrc")(4),
  col = pal_npg( "nrc")(4),
  fontface = "bold",
  #cat.col = c("#2E8B57","orange","#483D8B",""),
  cat.fontface = "bold",
  main = "COMSE",
  main.fontface = "bold",
  main.cex = 2
)
png("~/CellCycle/Venn_COMSE.png")
grid.draw(p)
dev.off()

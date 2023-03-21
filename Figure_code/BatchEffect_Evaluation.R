library(SeuratDisk)
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)


Convert('/home/luoqinhuan/HumanPancreas/human_pancreas_norm_complexBatch.h5ad', "h5seurat",
        overwrite = TRUE,assay = "RNA")
seurat_obj <- LoadH5Seurat("/home/luoqinhuan/HumanPancreas/human_pancreas_norm_complexBatch.h5seurat",assays ="RNA")
cell_type_retained = names(table(seurat_obj$celltype))[which(as.numeric(table(seurat_obj$celltype))>10)]
data_seurat = seurat_obj[,which(seurat_obj$celltype%in%cell_type_retained)]
data_seurat = subset(data_seurat,subset=tech%in%c("smarter","smartseq2"))
data_seurat = data_seurat[which(rowSums(data_seurat@assays$RNA@counts!=0)>=3),]
data_seurat = NormalizeData(data_seurat)
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = dim(data_seurat)[1])
feature_gene_seurat <- VariableFeatures(data_seurat)[1:2000]
all.genes <- rownames(data_seurat)
data_seurat <- ScaleData(data_seurat, features = all.genes)
#test = gsub("_","-",test)
data_seurat <- RunPCA(data_seurat, features = feature_gene_seurat)
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
data_seurat <- FindClusters(data_seurat, resolution = 0.5)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
data_seurat$cell = data_seurat$celltype
#score_function(data_seurat)
DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "seurat_clusters")+
  DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "celltype")+
  DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "tech")

cell_embedding = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
data_plot = cbind(cell_embedding,tech = data_seurat$tech, CellType = data_seurat$celltype,Cluster = data_seurat$seurat_clusters)
col5<-colorRampPalette((pal_npg( "nrc")(10)))(13)
p1 = ggplot(data_plot, aes(UMAP_1, UMAP_2, color = CellType)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  labs(title = "CellType") +
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values = col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(family = "Times New Roman", face = "bold", size = 16),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title= element_blank(),
        #legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p2 = ggplot(data_plot, aes(UMAP_1, UMAP_2, color =Cluster)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values = col5) +
  labs(x = "UMAP_1", y = "UMAP_2", title ="SeuratClusters") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(family = "Times New Roman", face = "bold", size = 16),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title= element_blank(),
        #legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank())

p3 = ggplot(data_plot, aes(UMAP_1, UMAP_2, color =tech)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values = col5) +
  labs(x = "UMAP_1", y = "UMAP_2", title ="Tech") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(family = "Times New Roman", face = "bold", size = 16),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title= element_blank(),
        #legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank())
p4 = ElbowPlot(data_seurat,n=30)
ggsave("~/COMSE_FigS3.jpg",(p1+p2)/(p3+p4),height = 8,width = 14)

mem = readRDS("/home/luoqinhuan/HumanPancreas/mem.RDS")
names(mem) = rownames(data_seurat)
markers = FindMarkers(data_seurat,group.by = "tech",ident.1 = "smarter",ident.2 = "smartseq2",only.pos = TRUE)
gene_num = c(50,100,200,400,600,800,1000,1200,1500)
diff_gene_prop = data.frame(matrix(0,length(table(mem)),length(gene_num)))
colnames(diff_gene_prop) = gene_num
rownames(diff_gene_prop) = paste("subgraph",names(table(mem)),sep="")

for(i in 1:length(gene_num)){
  count_inter = table(mem[which(names(mem)%in%rownames(markers)[1:gene_num[i]])])
  diff_gene_prop[paste("subgraph",names(count_inter),sep = ""),i] = as.numeric(count_inter)
}
diff_gene_prop = diff_gene_prop[rowSums(diff_gene_prop)>0,]
value = unlist(flatten(as.data.frame(t(diff_gene_prop))))
subgraph = rep(rownames(diff_gene_prop),each = length(gene_num))
number = rep(gene_num,length(rownames(diff_gene_prop)))
data_plot2 = data.frame(value,subgraph,number)
p1= ggplot(data_plot2,aes(x=number,y= value,shape = subgraph,color=subgraph))+
  geom_point(size=3)+
  geom_line(size=1)+
  scale_color_npg()+
  labs(x = "Differential Gene Number", y = "counts") +
  theme(axis.line = element_blank(),
        #axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(family = "Times New Roman", face = "bold", size = 16),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title= element_blank(),
        #legend.title = element_text(family = "Times New Roman", face = "bold", size = 14),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank())

data_plot3 = data.frame(subgraph = paste("subgraph",names(table(mem)),sep=""),count = as.numeric(table(mem)))
col5<-colorRampPalette((pal_npg( "nrc")(10)))(21)
p2 = ggplot(data_plot3,aes(x=subgraph,y=count,fill=subgraph))+
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85)+
  scale_fill_manual(values = col5)+
  labs(x = "Subgraph", y = "Counts",title = "Genes in Different Subgraph")+
  theme(axis.line = element_blank(),
        #axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        axis.text.x =  element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 12),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 10),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 
ggsave("~/COMSE_FigS4.jpg",p2/p1,width = 12,height = 12)
counts = data_seurat@assays$RNA@counts
counts = counts[which(rowSums(counts!=0)>=3),]
counts_norm = t(t(counts)/colSums(counts))
data_plot = data.frame(Total_counts = colSums(counts),tech = data_seurat$tech)
ggplot(data_plot,aes(x = tech,y = Total_counts))+
  geom_boxplot(aes(color=tech))+
  scale_color_npg()


data = log1p(sweep(counts_norm,1,rowSums(counts_norm),FUN = "/")*1e4)
#cor_matrix = cor(t(data))
pca = prcomp(data, scale=TRUE)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) 
plot(pca.var.per[1:30], main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")  
data_pca = pca$x[,1:3]
#data_pca = pca$x[,1:3]
graph = BuildSimKNN(data_pca,k = k_neighbors)
graph[graph!=0] = 1
G = graph.adjacency(graph,mode = "undirected")
fc <-multilevel.community(G,weights = E(G)$weight)
mem = fc$membership
group = table(mem)[which(table(mem)>=group_min)]
saveRDS(mem,"/home/luoqinhuan/HumanPancreas/mem.RDS")
#rm(graph)
cor_matrix_list = list()
for(i in 1:length(group)){
  data_sub = data[which(mem==names(group)[i]),]
  cor_matrix_list[[i]] = cor(t(data_sub))
}
saveRDS(cor_matrix_list,"/home/luoqinhuan/HumanPancreas/cor_matrix_list.RDS")

data = log1p(sweep(counts,2,colSums(counts),FUN = "/")*1e4)
t = 200
q = 0.95
gene_feature = character()
weight = numeric()
for(i in 1:length(group)){
  print(i)
  data_sub = data[which(mem==names(group)[i]),]
  score = RMLaplacianScore(data_sub,t=t,sample_size = 500,slices = 100)
  names(score) = rownames(data)[which(mem==names(group)[i])]
  genes = names(sort(score))
  system_uncertainty = cor_matrix_list[[i]]
  system_total = system_uncertainty[upper.tri(system_uncertainty)]
  threshold = max(0.5,quantile(system_total, q))
  weight[i] = quantile(system_total, q)
  print(threshold)
  genes_sorted = sort_feature(genes,system_uncertainty,threshold)
  #which(genes_sorted%in%gene1)
  gene_feature = c(gene_feature,genes_sorted)
}
saveRDS(gene_feature,"/home/luoqinhuan/HumanPancreas/gene_feature.RDS")

gene_feature = readRDS("/home/luoqinhuan/HumanPancreas/gene_feature.RDS")
mem =  readRDS("/home/luoqinhuan/HumanPancreas/mem.RDS")

names(mem) = rownames(data_seurat)
weight_art = rep(1,length(group))
weight_art[c(16,21)]=0
test = feature_gene_shrink(table(mem)[which(table(mem)>group_min)],2000,gene_feature,weight_art)
data_seurat$auc1 = data_auc1
data_seurat$auc2 = data_auc2
data_seurat$auc3 = data_auc3
data_seurat <- ScaleData(data_seurat,vars.to.regress = c("auc1","auc2","auc3"),model.use = "negbinom",features = intersect(test,feature_gene_seurat))
saveRDS(data_seurat,"~/HumanPancreas/data_seurat_auc1_auc2_auc3.rds")
data_seurat <- RunPCA(data_seurat, features = test)
data_seurat <- FindNeighbors(data_seurat, dims= 1:10)
#score_function(data_seurat)
data_seurat <- FindClusters(data_seurat, resolution = 0.5)
data_seurat <- RunUMAP(data_seurat, dims = 1:10)
data_seurat$cell = data_seurat$celltype

score_function(data_seurat)
score_COMSE = score_function(data_seurat)
score_COMSE = score_COMSE[which.max(rowSums(score_COMSE)),]

DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "seurat_clusters")+
  DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "celltype")+ DimPlot(data_seurat, reduction = "umap",label = TRUE,group.by = "tech")


marker = FindMarkers(data_seurat,ident.1 = "smarter",ident.2 = "smartseq2",only.pos = TRUE,group.by = "tech")

cells_rankings <- AUCell_buildRankings(as.matrix(counts))
genes <-rownames(data_seurat)[which(mem==1)]
geneSets <- GSEABase::GeneSet(genes, setName="member1") # alternative
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
par(mfrow=c(1,1))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
data_auc3 <- as.numeric(getAUC(cells_AUC)) 

genes <-rownames(data_seurat)[which(mem==16)]
geneSets <- GSEABase::GeneSet(genes, setName="member16") # alternative
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
par(mfrow=c(1,1))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
data_auc1 <- as.numeric(getAUC(cells_AUC)) 

genes <-rownames(data_seurat)[which(mem==21)]
geneSets <- GSEABase::GeneSet(genes, setName="member21") # alternative
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
par(mfrow=c(1,1))
set.seed(123)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
data_auc2 <- as.numeric(getAUC(cells_AUC)) 

data_plot = data.frame(mem16=data_auc1,mem21=data_auc2,mem1=data_auc3,seq = data_seurat$tech)

mu <- ddply(data_plot[,c(1,4)], "seq", summarise, grp.mean=mean(mem16))
p1_1=ggplot(data_plot, aes(x=mem16, fill=seq,color=seq)) +
  geom_histogram(binwidth=0.0005, aes(y=..density..), alpha=0.5, position="identity")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=seq), linetype="dashed") +
  scale_fill_npg() + 
  scale_color_npg() +
  theme_classic()+
  labs(x="Subgraph16")+
  theme(legend.position="top")

mu <- ddply(data_plot[,c(2,4)], "seq", summarise, grp.mean=mean(mem21))
p1_2=ggplot(data_plot, aes(x=mem21, fill=seq,color=seq)) +
  geom_histogram(binwidth=0.01, aes(y=..density..), alpha=0.5, position="identity")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=seq), linetype="dashed") +
  scale_fill_npg() + 
  scale_color_npg() +
  theme_classic()+
  labs(x="Subgraph21")+
  theme(legend.position="top")

mu <- ddply(data_plot[,c(3,4)], "seq", summarise, grp.mean=mean(mem1))
p1_3=ggplot(data_plot, aes(x=mem1, fill=seq,color=seq)) +
  geom_histogram(binwidth=0.002, aes(y=..density..), alpha=0.5, position="identity")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=seq), linetype="dashed") +
  scale_fill_npg() + 
  scale_color_npg() +
  theme_classic()+
  labs(x="Subgraph1")+
  theme(legend.position="top")


ggsave("~/HumanPancreas/Batch.jpg",p1_1/p1_2/p1_3+plot_layout(guides = "collect"),width=16,height = 12)
library(MASS)

fit_lda = lda(seq~.,data_plot)
x = predict(fit_lda,data_plot[,1:3])
lda_score = x$x

## Harmony
data_seurat_harmony = RunHarmony(data_seurat,group.by.vars = "tech")
data_seurat_harmony <- RunUMAP(data_seurat_harmony, reduction = "harmony", dims = 1:30)
data_seurat_harmony <- FindNeighbors(data_seurat_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters()
data_seurat_harmony$cell = data_seurat$celltype
score_function(data_seurat_harmony)
score_harmony = score_function(data_seurat_harmony)
score_harmony = score_harmony[which.max(rowSums(score_harmony)),]
#group_by_cluster
plot1 = DimPlot(data_seurat_harmony, reduction = "umap", label=T) + ggtitle("Tech (Harmony)")
#group_by_sample
plot1 = DimPlot(data_seurat_harmony, reduction = "umap", group.by='tech',cols = col5)+  ggtitle("Tech (Harmony)")
plot2 = DimPlot(data_seurat_harmony, reduction = "umap", group.by='celltype',cols = col5) + ggtitle("Cell Type (Harmony)")
#combinate
ggsave("~/HumanPancreas/harmony.jpg",plot1+plot2,dpi = 300,height = 8,width = 14)

pancreas.list <- SplitObject(data_seurat, split.by = "tech")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst",  nfeatures = 2000, verbose = FALSE)
}
reference.list <- pancreas.list[c("smarter", "smartseq2")]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims= 1:10)
pancreas.integrated$cell = pancreas.integrated$celltype
score_CCA = score_function(pancreas.integrated)
score_CCA = score_CCA[which.max(rowSums(score_CCA)),]
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech",cols = col5)+ ggtitle("Tech (Seurat-CCA)")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE,cols = col5) + ggtitle("Cell Type (Seurat-CCA)")
p1+p2

ggsave("~/HumanPancreas/seurat_cca.jpg",p1+p2,dpi = 300,height = 8,width = 14)

dataplot_graph = as.data.frame(data_seurat@reductions$umap@cell.embeddings)
dataplot_graph$CellType = data_seurat$celltype
dataplot_graph$Tech = data_seurat$tech
col5<-colorRampPalette((pal_npg( "nrc")(10)))(13)
p2_1=ggplot(dataplot_graph, aes(UMAP_1, UMAP_2, color = CellType)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  #scale_colour_manual(my36colors) + 
  scale_color_manual(values = col5) +
  labs(x = "UMAP_1", y = "UMAP_2", colour ="CellType") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p2_2=ggplot(dataplot_graph, aes(UMAP_1, UMAP_2, color = Tech)) +
  geom_point(size =1, shape = 19, alpha = 0.5) + 
  #scale_colour_manual(my36colors) + 
  scale_color_npg()+
  labs(x = "UMAP_1", y = "UMAP_2", colour ="Tech") +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.text = element_text(family = "Times New Roman",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

ggsave("~/HumanPancreas/Regressed_out_clustering.jpg",p2_1+p2_2,height = 6,width = 14)


score = rbind(score_COMSE,score_CCA,score_harmony)
rownames(score) = c("COMSE","Seurat-CCA","Harmony")
colnames(score) = c("Purity","F score","RI","ARI")

score_melt =  reshape2::melt(score)
colnames(score_melt)[1:2] = c("Method","Validations")

p_cluster=ggplot(score_melt, aes(x=Validations, y=value, fill = Method)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85) +
  scale_fill_npg()+
  labs(x = "", y = "Validations")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

library(scPOP)

# Batch Effect Evaluation
## ASW
ASW_comse = silhouette_width(data_seurat[,which(data_seurat$celltype%in%c("alpha","beta","delta","gamma"))]@reductions$umap@cell.embeddings, data_seurat[,which(data_seurat$celltype%in%c("alpha","beta","delta","gamma"))]@meta.data, 'tech')
ASW_harmony = silhouette_width(data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("alpha","beta","delta","gamma"))]@reductions$umap@cell.embeddings,data_seurat_harmony[,which(data_seurat$celltype%in%c("alpha","beta","delta","gamma"))]@meta.data, 'tech')
ASW_seurat = silhouette_width(pancreas.integrated[,which(pancreas.integrated$celltype%in%c("alpha","beta","delta","gamma"))]@reductions$umap@cell.embeddings, pancreas.integrated[,which(data_seurat$celltype%in%c("alpha","beta","delta","gamma"))]@meta.data, 'tech')

asw_score = data.frame(ASW = c(ASW_comse,ASW_harmony,ASW_seurat),Method=c("COMSE","Harmony","Seurat-CCA"))
asw_score$Method = factor(asw_score$Method,levels=c(c("COMSE","Seurat-CCA","Harmony")))
p_asw = ggplot(asw_score, aes(x=Method, y=ASW, fill = Method)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5,
           colour = "black", alpha = .85) +
  scale_fill_npg()+
  labs(x = "Method", y = "ASW-Batch")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 12),
        #legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

## LISI
lisi_comse_alpha = lisi(data_seurat[,which(pancreas.integrated$celltype%in%c("alpha"))]@reductions$umap@cell.embeddings, data_seurat[,which(pancreas.integrated$celltype%in%c("alpha"))]@meta.data, 'tech')$tech
lisi_harmony_alpha = lisi(data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("alpha"))]@reductions$umap@cell.embeddings, data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("alpha"))]@meta.data, 'tech')$tech
lisi_seurat_alpha = lisi(pancreas.integrated[,which(pancreas.integrated$celltype%in%c("alpha"))]@reductions$umap@cell.embeddings, pancreas.integrated[,which(data_seurat$celltype%in%c("alpha"))]@meta.data, 'tech')$tech


lisi_comse_beta = lisi(data_seurat[,which(pancreas.integrated$celltype%in%c("beta"))]@reductions$umap@cell.embeddings, data_seurat[,which(pancreas.integrated$celltype%in%c("beta"))]@meta.data, 'tech')$tech
lisi_harmony_beta = lisi(data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("beta"))]@reductions$umap@cell.embeddings, data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("beta"))]@meta.data, 'tech')$tech
lisi_seurat_beta = lisi(pancreas.integrated[,which(pancreas.integrated$celltype%in%c("beta"))]@reductions$umap@cell.embeddings, pancreas.integrated[,which(data_seurat$celltype%in%c("beta"))]@meta.data, 'tech')$tech

lisi_comse_delta = lisi(data_seurat[,which(pancreas.integrated$celltype%in%c("delta"))]@reductions$umap@cell.embeddings, data_seurat[,which(pancreas.integrated$celltype%in%c("delta"))]@meta.data, 'tech')$tech
lisi_harmony_delta = lisi(data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("delta"))]@reductions$umap@cell.embeddings, data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("delta"))]@meta.data, 'tech')$tech
lisi_seurat_delta = lisi(pancreas.integrated[,which(pancreas.integrated$celltype%in%c("delta"))]@reductions$umap@cell.embeddings, pancreas.integrated[,which(data_seurat$celltype%in%c("delta"))]@meta.data, 'tech')$tech


lisi_comse_gamma = lisi(data_seurat[,which(pancreas.integrated$celltype%in%c("gamma"))]@reductions$umap@cell.embeddings, data_seurat[,which(pancreas.integrated$celltype%in%c("gamma"))]@meta.data, 'tech')$tech
lisi_harmony_gamma = lisi(data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("gamma"))]@reductions$umap@cell.embeddings, data_seurat_harmony[,which(pancreas.integrated$celltype%in%c("gamma"))]@meta.data, 'tech')$tech
lisi_seurat_gamma = lisi(pancreas.integrated[,which(pancreas.integrated$celltype%in%c("gamma"))]@reductions$umap@cell.embeddings, pancreas.integrated[,which(data_seurat$celltype%in%c("gamma"))]@meta.data, 'tech')$tech


lisi_score = data.frame(score=c(lisi_comse_alpha,lisi_harmony_alpha,lisi_seurat_alpha,
                                lisi_comse_beta,lisi_harmony_beta,lisi_seurat_beta,
                                lisi_comse_delta,lisi_harmony_delta,lisi_seurat_delta,
                                lisi_comse_gamma,lisi_harmony_gamma,lisi_seurat_gamma)
                        ,Method = c(rep(c("COMSE","Seurat-CCA","Harmony"),c(length(lisi_comse_alpha),length(lisi_comse_alpha),length(lisi_comse_alpha))),
                                    rep(c("COMSE","Seurat-CCA","Harmony"),c(length(lisi_comse_beta),length(lisi_comse_beta),length(lisi_comse_beta))),
                                    rep(c("COMSE","Seurat-CCA","Harmony"),c(length(lisi_comse_delta),length(lisi_comse_delta),length(lisi_comse_delta))),
                                    rep(c("COMSE","Seurat-CCA","Harmony"),c(length(lisi_comse_gamma),length(lisi_comse_gamma),length(lisi_comse_gamma)))
                                    ),
                        celltype = c(rep("alpha",3*length(which(data_seurat$celltype=="alpha"))),
                                     rep("beta",3*length(which(data_seurat$celltype=="beta"))),
                                     rep("delta",3*length(which(data_seurat$celltype=="delta"))),
                                     rep("gamma",3*length(which(data_seurat$celltype=="gamma"))))
                        )
p_lisi=ggplot(lisi_score,aes(y=score,x=Method,fill=Method))+
  geom_boxplot()+
  scale_fill_npg()+
  facet_grid(.~celltype)+
  labs(x = "Method", y = "LISI")+
  theme(axis.line = element_blank(),
        axis.text = element_text(family = "Times New Roman", face = "bold", size = 12),
        axis.ticks = element_blank(),
        #axis.title = element_blank(),
        title =element_text(family = "Times New Roman", face = "bold", size = 16) ,
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.text = element_text(family = "Times New Roman",face = "bold",size = 12),
        legend.position = "none",
        #panel.background = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = .5),
        panel.grid = element_blank()) 

p_evaluate = (p_cluster+p_asw)+plot_layout(guides = "collect")
ggsave("~/HumanPancreas/Evaluations.jpg",p_evaluate,dpi=300,height=8,width = 14)

library(Seurat)
library(ggplot2)

load("qsc_wt_ezh1ko.RData")

DimPlot(object = qsc_wt_ezh1ko, reduction = 'umap',label=T)+labs(title = "qsc_wt_ezh1ko")


cell_id0<-WhichCells(qsc_wt_ezh1ko, idents = "0")
cell_id1<-WhichCells(qsc_wt_ezh1ko, idents = "1")
cell_id2<-WhichCells(qsc_wt_ezh1ko, idents = "2")
cell_id10<-WhichCells(qsc_wt_ezh1ko, idents = "10")
num_id0<-c(length(grep("-1",cell_id0)),length(grep("-2",cell_id0)))
num_id1<-c(length(grep("-1",cell_id1)),length(grep("-2",cell_id1)))
num_id10<-c(length(grep("-1",cell_id10)),length(grep("-2",cell_id10)))
num_id2<-c(length(grep("-1",cell_id2)),length(grep("-2",cell_id2)))

qsc_cell_id<-rbind(
  data.frame(cell_id=cell_id0,cluster=rep(0,time=length(cell_id0))),
  data.frame(cell_id=cell_id1,cluster=rep(1,time=length(cell_id1))),
  data.frame(cell_id=cell_id2,cluster=rep(2,time=length(cell_id2)))
)

qsc_wt_ezh1ko.filtered <- subset(qsc_wt_ezh1ko, idents = c("0", "1", "2"))
which.max(unlist(FetchData(qsc_wt_ezh1ko.filtered, "UMAP_2")))
#save(qsc_wt_ezh1ko,file="qsc_wt_ezh1ko.RData")
#qsc_wt_ezh1ko.filtered.mtx <- as.matrix(GetAssayData(qsc_wt_ezh1ko.filtered,slot = "data"))

DimPlot(object = qsc_wt_ezh1ko.filtered, reduction = 'umap',label=F)+labs(title = "qsc_wt_ezh1ko")
new.cluster.ids <-  c("Ezh1KO","WT","WT+Ezh1KO")
names(new.cluster.ids) <- levels(qsc_wt_ezh1ko.filtered)
qsc_wt_ezh1ko.filtered<- RenameIdents(qsc_wt_ezh1ko.filtered, new.cluster.ids)
DimPlot(object = qsc_wt_ezh1ko.filtered, reduction = 'umap',label=F)+labs(title = "qsc_wt_ezh1ko")
#figure1 G
#compare clusters
library(dplyr)
library(cowplot)
qsc_wt_ezh1ko.only<- subset(qsc_wt_ezh1ko.filtered, idents = c("WT", "Ezh1KO"))
qsc_wt_ezh1ko.filtered.markers <- FindAllMarkers(qsc_wt_ezh1ko.only, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
#write.csv(qsc_wt_ezh1ko.filtered.markers,file="qsc_wt_ezh1ko_filtered_markers.csv")

top40 <- qsc_wt_ezh1ko.filtered.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)

heat_genes<- c('Rpl29', 'Rpl38','Rps25','Rps28','Hspa5', 'Hsp90b1', 'Dnajc3', 'Spry1',
               'Gas1', 'Pax7', 'Hey1', 'Hes1','Myod1','Notch3', 'Rbpj','Heyl','Sdc4','Meg3','Myod1',
               'Txnip', 'Ccnd3', 'Fgfr4', 'Calcr','Hist1h2bc','Bmyc','Pde10a','Rp9','Smc6', 'Col3a1', 'Idh2', 'Usp2')

DoHeatmap(qsc_wt_ezh1ko.only, features = top40$gene) +theme(text = element_text(size = 8))

#figure1 H
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Rpl29', 'Rpl38','Rps25','Rps28'),min.cutoff = "q8")
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Hspa5', 'Hsp90b1', 'Dnajc3', 'Spry1'),min.cutoff = "q8")
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Gas1', 'Pax7', 'Hey1', 'Hes1'))
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Myod1'),min.cutoff = "q8")
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Notch3', 'Rbpj','Heyl'))
FeaturePlot(qsc_wt_ezh1ko.filtered, features = c('Sdc4','Meg3','Myod1'),min.cutoff = "q8")

#monocle
library(monocle3)
library(SeuratWrappers)
qsc_wt_ezh1ko.filtered_cds<-as.cell_data_set(qsc_wt_ezh1ko.filtered)
#qsc_wt_ezh1ko.filtered_cds<- cluster_cells(qsc_wt_ezh1ko.filtered_cds)
#qsc_wt_ezh1ko.filtered_cds <- learn_graph(qsc_wt_ezh1ko.filtered_cds)
#qsc_wt_ezh1ko.filtered_cds <- learn_graph(qsc_wt_ezh1ko.filtered_cds,close_loop = F,use_partition = T,learn_graph_control =list(minimal_branch_len=9))
qsc_wt_ezh1ko.filtered_cds<- cluster_cells(qsc_wt_ezh1ko.filtered_cds,k = 20,resolution = 0.00002)
q#sc_wt_ezh1ko.filtered_cds <- learn_graph(qsc_wt_ezh1ko.filtered_cds,close_loop = F,use_partition = T)
qsc_wt_ezh1ko.filtered_cds <- learn_graph(qsc_wt_ezh1ko.filtered_cds,close_loop = F,use_partition = T,learn_graph_control =list(minimal_branch_len=5))
plot_cells(qsc_wt_ezh1ko.filtered_cds, label_groups_by_cluster = T, label_leaves = F, label_branch_points = T,,graph_label_size = 3)

wt_qsc.min.umap <- which.min(unlist(FetchData(qsc_wt_ezh1ko.filtered, "UMAP_2")))
wt_qsc.min.umap <- colnames(qsc_wt_ezh1ko.filtered)[wt_qsc.min.umap]
qsc_wt_ezh1ko.filtered_cds <- order_cells(qsc_wt_ezh1ko.filtered_cds, root_cells = wt_qsc.min.umap)
plot_cells(qsc_wt_ezh1ko.filtered_cds, color_cells_by = "pseudotime", label_cell_groups =T, label_leaves = F, 
           label_branch_points = T,show_trajectory_graph = T,graph_label_size = 3,label_groups_by_cluster = T)

save(qsc_wt_ezh1ko.filtered_cds,file="qsc_wt_ezh1ko.filtered_cds.RData")

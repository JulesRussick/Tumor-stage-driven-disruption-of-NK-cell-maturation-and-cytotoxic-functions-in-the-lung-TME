library(dplyr)
library(data.table)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(ggplot2)
library(BiocManager)
library(MAST)
library(EnhancedVolcano)
library(harmony)
library(Nebulosa)
library(scater)
library(scran)
library(ggalt)
library(DropletUtils)
library(GEOquery)
library(Biobase)
library(Matrix)
library(stringr)
library(dittoSeq)
library(readxl)
library(tidyverse)


#To download saved files
combined_seurat_object <- readRDS("pathway/2023-12-12_combined_seurat_object.rds")
NK_seurat_object <- readRDS("pathway/2023-12-23_Subset_NK_seurat_object.rds")

#-----------------------------------------------------------------------------------------------------------------
#-------------------------------------------------- Import data -----------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# Downloading data
combined_seurat_object <- load("pathway/lung_ldm.rd")
combined_seurat_object <- CreateSeuratObject(lung_ldm$dataset$umitab)
# Downsampling
combined_seurat_object <- CreateSeuratObject(lung_ldm$dataset$ds[[1]])
set.seed(123)
combined_seurat_object <- NormalizeData(combined_seurat_object)
combined_seurat_object <- FindVariableFeatures(combined_seurat_object)
combined_seurat_object <- SketchData(combined_seurat_object,
                                     ncells = 50000,
                                     method = "LeverageScore",
                                     sketched.assay = "sketch")

# Switch to analyze the sketched dataset (in-memory)
DefaultAssay(combined_seurat_object) <- "sketch"
# Clustering on a sketched dataset
combined_seurat_object <- FindVariableFeatures(combined_seurat_object)
combined_seurat_object <- ScaleData(combined_seurat_object)
combined_seurat_object <- RunPCA(combined_seurat_object)
combined_seurat_object <- FindNeighbors(combined_seurat_object, dims = 1:40)
combined_seurat_object <- FindClusters(combined_seurat_object, resolution = 2)

#Extend results to the full datasets
combined_seurat_object <- ProjectData(
  object = combined_seurat_object,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:40,
  refdata = list(cluster_full = "seurat_clusters")
)

DefaultAssay(combined_seurat_object) <- "RNA"
DimPlot(combined_seurat_object, reduction = "umap",label = T,label.size = 5,label.box = F, label.color = "black",pt.size = 1) 
#+ NoLegend()

# Add Tum_stage to metadata according to sampleID
combined_seurat_object@meta.data <- combined_seurat_object@meta.data %>%   
  mutate(Tum_stage = case_when(
    sampleID %in% c('38','40','43','45','47','57','58','60','62','95','97','297','298','313','317','318','324','328','332','334','342','415','851') ~ 'Stage_I',  
    sampleID %in% c('42','48','51','52','54','99','103','107','111','311','326','340','398','414','448','452','456','460','464','468','472') ~ 'Stage_II',             
    sampleID %in% c('315','320','338') ~ 'Stage_III',                  
    sampleID %in% '36' ~ 'Unknown',
    sampleID %in% c('37','39','41','44','46','49','50','53','55','56','59','61','63','96','98','100','104','108','112','310','312','314','316','319','323','325','327','331','333','335','337','339','341','399','413','449','453','457','461','465','469','473') ~ 'Stage_NT')
  )

saveRDS(combined_seurat_object,paste0("/Users/jules.russick/Desktop/Analyse_papier_NK/Docs/Human data/231128 Merad/",Sys.Date(),"_combined_seurat_object.rds"))


# Identifying NK cell cluster
FeaturePlot(
  object = combined_seurat_object,
  reduction = "umap",
  features = c("gene_names"),
  ncol = 3,
  cols = c("green","purple"),
  pt.size = 0.1,
  raster = F
)

# ----------------------------------------- subsetting on selected clusters --------------------------------------------

# Subset the Seurat object to select specific clusters
NK_seurat_object <- subset(subsubset_seurat_object, ident = c(selected_clusters))
DimPlot(NK_seurat_object, reduction = "umap",label = T,label.size = 8,label.box = T, label.color = "black",pt.size = 1) 

### harmony
DefaultAssay(NK_seurat_object) <- 'RNA'
NK_seurat_object <- NK_seurat_object %>%
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(vars.to.regress = c('nCount_RNA','percent.mt','S.Score','G2M.Score')) %>% 
  RunPCA(npcs = 40) %>% 
  RunHarmony('sampleID', reduction = 'pca', reduction.save = "harmony", plot_convergence = T) %>% 
  RunUMAP(reduction = "harmony", dims = 1:10, reduction.name = 'umap') %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 1)

DimPlot(NK_seurat_object, reduction = "umap",label = T, label.size = 5,label.box = F, label.color = "black", pt.size = 1)
#+ NoLegend()


#-----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------- Visualizations --------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# UMAPs
FeaturePlot(NK_seurat_object, features = "gene_name",cols = c("grey20","green"),pt.size = 2,order = T, raster = F,reduction = "umap")
DimPlot(immune.combined, reduction = "umap",label = T,label.size = 10,label.box = T, label.color = "black",pt.size = 1)
plot_density(immune.combined, "gene_name",size = 2, reduction = "umap")

# Violins Plots
VlnPlot(immune.combined, features = "gene_name",group.by = "seurat_clusters",sort = F) + NoLegend()
VlnPlot(immune.combined, features = c("gene_names"),
        group.by = "orig.ident",
        log = F, 
        sort = "increasing"
)


# BubblePlot
genes_of_interest <- c("gene_names")
dittoDotPlot(
  NK_seurat_object,
  group.by = "seurat_clusters",
  vars = genes_of_interest,
  #assay = assay_name,
  size = 14,
  min.color = "blue",
  max.color = "red"
)

# OR
DotPlot(NK_seurat_object, features = genes_of_interest)

#-----------------------------------------------------------------------------------------------------------------
#--------------------------------------------------- Calcul DEG --------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#1-Non-Tumoral, 2-Tumoral, 3-Grosse-Tumeur
diff_markers.allclust <- FindMarkers(subset_refquery_final,
                                     test.use = "MAST",
                                     ident.1 = c(gp_1_clusters),
                                     ident.2 = c(gp_2_clusters),
                                     logfc.threshold = 0.1,
                                     min.pct = 0.1
)

diff_markers.allclust

#-----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ Pseudotime -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggpubr)
library(reshape2)

cds <- readRDS("pathway/2023-12-26_monocle_all.rds")

#-------

cds <- as.cell_data_set(NK_seurat_object)
## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)
## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(immune.combined[["RNA"]])

cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
#str(cds)
plot_cells(cds, color_cells_by = "partition")

# each partition correspond to a different trajectory
cds <- learn_graph(cds,close_loop = F)
plot_cells(cds, color_cells_by = "cluster",label_branch_points = T,label_leaves = T,label_roots = T)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds,color_cells_by = "pseudotime", label_branch_points = T,label_leaves = T,label_roots = T,show_trajectory_graph = T,trajectory_graph_segment_size = 3,trajectory_graph_color = "black",cell_size = 1)

# #Subset cells by branch
cds_sub <- choose_graph_segments(cds)

#Different expression analysis - Which genes change as a function of pseudotime
NK_test <- graph_test(cds,neighbor_graph = "principal_graph")
pr_deg_ids <- row.names(subset(NK_test, q_value < 0.05))

#-----------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ Visuatizations -------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

#UMAP 3D Pseudotime
cds_3D <- reduce_dimension(cds, max_components = 3)
cds_3D <- cluster_cells(cds_3D)
cds_3D <- learn_graph(cds_3D,close_loop = F)
cds_3D <- order_cells(cds_3D)
cds_3D_plot_obj <- plot_cells_3d(cds_3D, color_cells_by = "pseudotime")
cds_3D_plot_obj


#Genes on UMAP
plot_cells(cds, genes=c("gene_names"),
           show_trajectory_graph=T,
           label_cell_groups=F,
           label_leaves=F)


# Gene expression accross pseudotime - Escalière's code
library("KernSmooth")
library(ggpubr)

pData(cds)$pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime
exprData <- t(as.matrix(GetAssayData(immune.combined[c("gene_names"),])))
pseudotime_expr.df <- data.frame(group = pData(cds)$orig.ident, pseudotime = pData(cds)$pseudotime, as.data.frame(exprData))
pseudotime_expr.df$cell <- rownames(pseudotime_expr.df)

# All groups
pseudotime_expr.melt.df <- melt(pseudotime_expr.df, id.vars = c('group','cell', 'pseudotime'), variable.name = 'gene', value.name = 'expression')
pseudotime_expr.melt.df <- subset(pseudotime_expr.melt.df, ! is.infinite(pseudotime))
max_pseudotime <- max(pseudotime_expr.melt.df[,"pseudotime"])
min_pseudotime <- min(pseudotime_expr.melt.df[,"pseudotime"])

kernel_smoother_bandwith <- max_pseudotime/5
kernel_smoother.df = data.frame()
genes <- unique(pseudotime_expr.melt.df$gene)

for (agene in genes){
  genes_kernel_smoother_list <- locpoly(x=subset(pseudotime_expr.melt.df, gene == agene)$pseudotime, y=subset(pseudotime_expr.melt.df, gene == agene)$expression, bandwidth = kernel_smoother_bandwith, gridsize=1000)
  genes_kernel_smoother.df <- data.frame(pseudotime=genes_kernel_smoother_list[[1]], expression=genes_kernel_smoother_list[[2]])
  genes_kernel_smoother.df$gene <- agene
  kernel_smoother.df <- rbind(kernel_smoother.df, genes_kernel_smoother.df)
}

min_max_df <- kernel_smoother.df %>%
  group_by(gene) %>%
  summarise(min_expr = min(expression),
            max_expr = max(expression))
kernel_smoother.df <- kernel_smoother.df %>%
  left_join(min_max_df, by = "gene") %>%
  group_by(gene) %>% 
  mutate(scaled_expression = (expression - min_expr) / (max_expr - min_expr))
ggscatter(kernel_smoother.df, 
          x = 'pseudotime', 
          y = 'scaled_expression', 
          color = 'gene')

# Per group
pseudotime_expr.melt.df <- melt(pseudotime_expr.df, id.vars = c('group','cell', 'pseudotime'), variable.name = 'gene', value.name = 'expression')
calcKernelSmoother <- function(df, aGroup){
  group.df <- subset(df, ! is.infinite(pseudotime) & group == aGroup)
  max_pseudotime <- max(group.df[,"pseudotime"])
  min_pseudotime <- min(group.df[,"pseudotime"])
  kernel_smoother_bandwith <- max_pseudotime/5
  kernel_smoother.df = data.frame()
  genes <- unique(group.df$gene)
  for (agene in genes){
    genes_kernel_smoother_list <- locpoly(x=subset(group.df, gene == agene)$pseudotime, y=subset(group.df, gene == agene)$expression, bandwidth = kernel_smoother_bandwith, gridsize=1000)
    genes_kernel_smoother.df <- data.frame(pseudotime=genes_kernel_smoother_list[[1]], expression=genes_kernel_smoother_list[[2]])
    genes_kernel_smoother.df$gene <- agene
    kernel_smoother.df <- rbind(kernel_smoother.df, genes_kernel_smoother.df)
  }
  min_max_df <- kernel_smoother.df %>%
    group_by(gene) %>%
    summarise(min_expr = min(expression),
              max_expr = max(expression))
  kernel_smoother.df <- kernel_smoother.df %>%
    left_join(min_max_df, by = "gene") %>%
    group_by(gene) %>% 
    mutate(scaled_expression = (expression - min_expr) / (max_expr - min_expr))
  p <- ggscatter(kernel_smoother.df, 
                 x = 'pseudotime', 
                 y = 'scaled_expression', 
                 color = 'gene',
                 title = aGroup)
  
  return(p)
}


p1 <- calcKernelSmoother(pseudotime_expr.melt.df, '1-Non-Tumoral')
p1 <- p1 + labs(x = "Pseudotime", y = "Relative expression")
#p1 <- p1 + NoLegend()
p2 <- calcKernelSmoother(pseudotime_expr.melt.df, '2-Tumoral')
p2 <- p2 + labs(x = "Pseudotime", y = "Relative expression")
p2 <- p2 + NoLegend()
p3 <- calcKernelSmoother(pseudotime_expr.melt.df, '3-Grosse-Tumeur')
p3 <- p3 + labs(x = "Pseudotime", y = "Relative expression")
p3 <- p3 + NoLegend()

p <- p1 + p2 + p3
p

#-----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------- MiloDE --------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

library("miloDE")
library("miloR")
library("Seurat")
library("uwot")
library("ggplot2")
library("SingleCellExperiment")
library("scater")
library("scran")
library("dplyr")
library("patchwork")


# Rename clusters
clus_id <- c("Cytotoxic", "Cytokine", "Cytotoxic", "Cytokine", "CXCR5+", "Cytokine", "Proliferative", "Cytokine", "Cytotoxic")
names(clus_id) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined,clus_id)
DimPlot(immune.combined,
        reduction = "umap",
        label = T,
        label.size = 7,
        label.box = T,
        repel = T,
        pt.size = 0.3) + NoLegend()

# Add to metadata
immune.combined$clus_id = immune.combined$seurat_clusters
immune.combined$clus_id = ifelse(immune.combined$clus_id%in%c("0","2","11"),"Cytotoxic",ifelse(immune.combined$clus_id%in%c("1","3","6","8"),"Cytokine",ifelse(immune.combined$clus_id%in%c("4"),"CXCR5+","Proliferative")))
sce_mouseEmbryo <- as.SingleCellExperiment(immune.combined, assay = "RNA")

# Add reducedDim
umap_embeddings <- Embeddings(immune.combined, "umap")
reducedDim(sce_mouseEmbryo) <- umap_embeddings
reducedDimNames(sce_mouseEmbryo) -> umap_embeddings
table(colData(sce_mouseEmbryo)[,c('orig.ident','clus_id')])
sce_mouseEmbryo = assign_neighbourhoods(sce_mouseEmbryo, k = 30, order = 2,
                                        filtering = TRUE,
                                        reducedDim_name = "UMAP"
)

# See how neighbourhood size evolves acoording to k value
stat_k = estimate_neighbourhood_sizes(sce_mouseEmbryo, k_grid = seq(10,40,5) ,
                                      order = 2, prop = 0.1 , filtering = TRUE,
                                      reducedDim_name = "UMAP" , plot_stat = TRUE)

de_stat = de_test_neighbourhoods(sce_mouseEmbryo , sample_id = "orig.ident",
                                 design = ~clus_id, covariates = c("clus_id"))

# Annotate neighboorhoods
de_stat = miloR::annotateNhoods(sce_mouseEmbryo , de_stat , coldata_col = "clus_id")

#-----------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------- MiloR --------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

library("miloDE")
library("miloR")
library("Seurat")
library("uwot")
library("ggplot2")
library("SingleCellExperiment")
library("scater")
library("scran")
library("dplyr")
library("patchwork")
library("viridis")
library("viridisLite")
library("ggpubr")
library("Hmisc")
library("reshape2")
library("knitr")

# Renommer clusters
clus_id <- c("Cytotoxic", "Cytokine", "Cytotoxic", "Cytokine", "CXCR5+", "Cytokine", "Proliferative", "Cytokine", "Cytotoxic")
names(clus_id) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined,clus_id)
DimPlot(immune.combined,
        reduction = "umap",
        label = T,
        label.size = 7,
        label.box = T,
        repel = T,
        pt.size = 1.5) + NoLegend()

#Add to metadata 
immune.combined$clus_id = immune.combined$seurat_clusters
immune.combined$clus_id = ifelse(immune.combined$clus_id%in%c("0","2","11"),"Cytotoxic",ifelse(immune.combined$clus_id%in%c("1","3","6","8"),"Cytokine",ifelse(immune.combined$clus_id%in%c("4"),"CXCR5+","Proliferative")))
p_z<-DimPlot(immune.combined,group.by = "clus_id",
             reduction = "umap",
             label = T,
             label.size = 7,
             label.box = T,
             repel = T,
             pt.size = 0.3) + NoLegend()
p_z
sce_NK <- as.SingleCellExperiment(immune.combined, assay = "RNA")

# Add reducedDim
umap_embeddings <- Embeddings(immune.combined, "umap")
reducedDim(sce_NK) <- umap_embeddings
reducedDimNames(sce_NK) -> umap_embeddings

print(sce_NK)

# Check which dimensions we have in the object
reducedDimNames(sce_NK)
#[1] "PCA"       "UMAP"

# Rstimate_neighbourhood_sizes allows to gauge how neighbourhood size distribution
# changes as a function of (order,k). It might be useful to run it first in order to 
# determine optimal range that will return desired neighbourhood sizes.
stat_k = estimate_neighbourhood_sizes(sce_NK, 
                                      k_grid = seq(10,50,5),
                                      order = 2, 
                                      prop = 0.1,
                                      filtering = TRUE,
                                      reducedDim_name = "UMAP" , plot_stat = TRUE)
kable(stat_k , caption = "Neighbourhood size distribution ~ k")
sce_NK_milo = assign_neighbourhoods(sce_NK, k = 50, order = 2, 
                                    filtering = TRUE, 
                                    reducedDim_name = "UMAP"
)

nhoods_sce = nhoods(sce_NK_milo)
nhood_stat_ct = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
nhood_stat_ct = miloR::annotateNhoods(sce_NK_milo , nhood_stat_ct , coldata_col = "clus_id")
head(nhood_stat_ct)
p = plot_milo_by_single_metric(sce_NK_milo, nhood_stat_ct, colour_by = "clus_id" , 
                               layout = "UMAP" , size_range = c(5,5) , edge_width = c(0.2,0.5),node_stroke = 0.5) 
p = p + scale_fill_manual(values = c("#00B9E3","#93AA00","#F8766D","#DB72FB")) 
p

#DE testing
#Calculate AUC per neighbourhood
stat_auc = suppressWarnings(calc_AUC_per_neighbourhood(sce_NK_milo , sample_id = "orig.ident" , condition_id = "clus_id", min_n_cells_per_sample = 3,
                                                       n_threads = 12))
p2 = plot_milo_by_single_metric(sce_NK_milo, stat_auc, colour_by = "auc" , 
                                layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "AUC")
p2
de_stat = de_test_neighbourhoods(sce_NK_milo , sample_id = "orig.ident", 
                                 design = ~clus_id, covariates = c("clus_id"),
                                 plot_summary_stat=T,
                                 layout = "UMAP")
de_stat$clus_id<-nhood_stat_ct$clus_id[match(de_stat$Nhood,nhood_stat_ct$Nhood)]
stat_de_magnitude = rank_neighbourhoods_by_DE_magnitude(de_stat)
p2 = plot_milo_by_single_metric(sce_NK_milo, stat_de_magnitude, colour_by = "n_specific_DE_genes" , 
                                layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# specific\nDE genes" , option = "inferno")
p2

# Visualizing individual genes
genes<-c("gene_names")
plots = lapply(genes , function(gene){
  p = plot_DE_single_gene(sce_NK_milo, de_stat , gene = genes , layout = "UMAP" , set_na_to_0 = TRUE) + 
    ggtitle(gene)
  return(p)
})
p_genes = ggarrange(plotlist = plots,p,p_z , ncol = 2, nrow = 3)
p_only_genes = ggarrange(plotlist = plots,p, ncol = 3, nrow = 3)
umap_annot = p_z
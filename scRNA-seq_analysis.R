rm(list=ls())

library(Seurat)
library(Matrix)
library(dplyr) # Dataframe manipulation
library(monocle)
library(DDRTree)
library(data.table)
library(pheatmap)
library(corrplot)
library(dplyr)
library(ggsci)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(patchwork)
library(magrittr)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(dorothea)
library(progeny)
#############################
#######sample-object construct########
###### use P1_Dx as an example######
#load data from cellranger pipeline
P1_Dx_data <- Read10X("./P1_Dx/outs/filtered_feature_bc_matrix/")
P1_Dx_meta <- read.csv("./P1_Dx/P1_Dx_doublet_annot.csv",header = T,row.names = 1)
#####Creat P1_Dx project######
P1_Dx <- CreateSeuratObject(counts = P1_Dx_data, project = "P1_Dx", 
                                    meta.data = P1_Dx_meta)
P1_Dx <- CreateSeuratObject(counts = P1_Dx_data[[ "Gene Expression" ]], project = "P1_Dx", min.cells = 0, min.features = 0,meta.data = P1_Dx_meta)
rownames(P1_Dx_data[[ "Antibody Capture" ]]) <- sub("_TotalC", "", rownames(P1_Dx_data[[ "Antibody Capture" ]]))
P1_Dx[[ "Antibody" ]] <- CreateAssayObject(counts = P1_Dx_data[[ "Antibody Capture" ]][, colnames(P1_Dx)] )
##### singleR prediction #####
ALL_sce_i <-  as.SingleCellExperiment(P1_Dx)
ref_MonacoImmune <- celldex::MonacoImmuneData()
ref.data <- ref_MonacoImmune
new.data <- ALL_sce_i
predictions_label.main <- SingleR(test=new.data, assay.type.test=1,   ref=ref.data, labels=ref.data$label.main)
table(predictions_label.main$pruned.labels)
predictions_label.fine <- SingleR(test=new.data, assay.type.test=1,   ref=ref.data, labels=ref.data$label.fine)
table(predictions_label.fine$pruned.labels)
P1_Dx$prediction_singleR_lable.main <- predictions_label.main$pruned.labels[match(colnames(P1_Dx),row.names(predictions_label.main))]
P1_Dx$prediction_singleR_lable.fine <- predictions_label.fine$pruned.labels[match(colnames(P1_Dx),row.names(predictions_label.fine))]
###### remove potential doubles #####
P1_Dx <- subset(x = P1_Dx,subset= (doublet_detected == 0))
##### remove low quality cells #####
P1_Dx[["percent.mt"]] <- PercentageFeatureSet(P1_Dx, pattern = "^MT-")
P1_Dx[["percent.rp"]] <- PercentageFeatureSet(P1_Dx, pattern = "^(MRP|RP)[LSFN]")
P1_Dx[["percent.hb"]] <- PercentageFeatureSet(P1_Dx, pattern = "^(HBA1|HBA2|HBB|HBG1|HBG2|HBQ1|HBD|HBM|HBE1|HBZ)$")
P1_Dx <- subset(x = P1_Dx,subset= (nCount_RNA > 800 ) & (nFeature_RNA > 400) &  (percent.mt < 10 ) )

#############################
#######Figure1######
###### pseudotime analysis with monocle2#######
load("./pseudotime.Rdata")
CD8_T_Pathway_cds <- setOrderingFilter(CD8_T_Pathway_cds, ordering_genes = ordering_genes)
CD8_T_Pathway_cds <- reduceDimension(CD8_T_Pathway_cds, max_components = 2, method = 'DDRTree')
CD8_T_Pathway_cds <- orderCells(CD8_T_Pathway_cds)
plot_cell_trajectory(CD8_T_Pathway_cds,  color_by = "cell_types") + scale_color_manual(values=c(
    '#5F559BFF','#E7BA52EF','#E377C2FF','#8C564BFF','#BB0021FF'))
plot_cell_trajectory(CD8_T_Pathway_cds,  color_by = "group") + scale_color_manual(values=c(
    '#636363FF','#D41159'))
#####DEGs identification##########
DefaultAssay(seuratObj_all) <- "RNA"
seuratObj_all  %<>% NormalizeData(.) %<>% ScaleData()
DefaultAssay(seuratObj_all) <- 'RNA'
Idents(seuratObj_all) <- "celltypes_use"
AllMarkers <- FindAllMarkers(seuratObj_all, only.pos =FALSE,logfc.threshold=0 ,min.cells.feature=0,min.pct = 0,return.thresh = 0.05)

#############################
#######Figure2-Sfigure4######
##### TCR analysis of bulk RNA-seq #####
# run run-trust4 of software TRUST first
load("./TCR.Rdata")
library(RColorBrewer)
mycolors <- colorRampPalette(col_vector)(nrow(df))
fish = setCol(fish, mycolors)
suppressWarnings(
fishPlot(fish,shape="spline",title.btm=title,
         cex.title=2, vlines=c(0,150),
         cex.vlab=1,
         vlab=c(vlab))
)

#############################
#######Figure3######
##### GO analysis #####
xx <- compareCluster(gcSample,
                     fun="enrichGO",#One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .
                     #organism="hsa",
                     OrgDb = org.Hs.eg.db,# org.Hs.eg.db,org.Mm.eg.db,
                     ont = "ALL" ,
                     pool = TRUE ,
                     keyType = "SYMBOL",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.2,
                     minGSSize = 10,
                     maxGSSize = 500,
                     readable = FALSE
                     #universe
                    )
##### KEGG analysis #####
y <- gseKEGG(geneList     = geneList_sort_entrez,
               organism     = 'hsa',
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               verbose      = FALSE) %>% setReadable( OrgDb = org.Hs.eg.db, keyType="ENTREZID")
##### cell cycle analysis #####
P1_Dx <- CellCycleScoring(
  object = P1_Dx,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
##### gene signatures scoring of bulk RNA-seq #####
tf_activities_counts <-
    dorothea::run_viper(counts, regulons,
    options =  list(minsize = 5, eset.filter = FALSE,
    cores = 1, verbose = FALSE, method = c("scale")))

#############################
#######Figure5######
##### Calculate the activity of pathways in bulk RNA-seq #####
PathwayActivity_counts <- progeny(counts, scale=T, z_scores =F,
    organism="Human", top = 100)

#############################
#######Sfigure5##############
##### cell interaction analysis #####
load("./cell_interaction.Rdata")
# Define the contrasts and covariates of interest for the DE analysis
contrasts_oi = c("'Dx-HD','HD-Dx'")
contrast_tbl = tibble(contrast = c('Dx-HD','HD-Dx'), group = c("Dx",'HD'))
# Define the weights of the prioritization of both expression, differential expression and NicheNet activity information
# We will set our preference for this dataset as follows - and recommend the user to use the same weights by default:
prioritizing_weights_DE = c("de_ligand" = 1,
                         "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                         "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                         "abund_receiver" = 0)
prioritizing_weights = c(prioritizing_weights_DE,
                         prioritizing_weights_activity,
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency,
                         prioritizing_weights_relative_abundance)
abundance_expression_info = get_abundance_expression_info(sce = sce,
                                                          sample_id = sample_id,
                                                          group_id = group_id,
                                                          celltype_id = celltype_id,
                                                          min_cells = min_cells,
                                                          senders_oi = senders_oi,
                                                          receivers_oi = receivers_oi,
                                                          lr_network = lr_network,
                                                          batches = batches)
DE_info = get_DE_info(sce = sce,
                      sample_id = sample_id,
                      group_id = group_id,
                      celltype_id = celltype_id,
                      batches = batches,
                      covariates = covariates,
                      contrasts_oi = contrasts_oi,
                      min_cells = min_cells)
DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
# Run the NicheNet ligand activity analysis
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose,
  n.cores = n.cores
)))
# Make necessary grouping data frame
# Crucial note: grouping_tbl: group should be the same as in the contrast_tbl, and as in the expression info tables! Rename accordingly
# if this would not be the case. If you followed the guidelines of this tutorial closely, there should be no problem.
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}
# Run the prioritization
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))
lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
                                                    abundance_expression_info,
                                                    celltype_de,
                                                    grouping_tbl,
                                                    prioritization_tables,
                                                    ligand_target_matrix,
                                                    logFC_threshold = logFC_threshold,
                                                    p_val_threshold = p_val_threshold,
                                                    p_val_adj = p_val_adj)
prioritized_tbl_oi_Dx_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 20, groups_oi = group_oi,senders_oi = senders_oi, receivers_oi = receivers_oi,)
plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_Dx_50)
plot_oi


		
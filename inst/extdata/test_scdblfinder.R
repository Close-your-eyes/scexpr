library(scexpr)
library(scDblFinder)
library(ggplot2)
so <- readRDS("/Volumes/CMS_SSD_2TB/R_scRNAseq/2020_10XGenomics_PBMCs/data/SO_processed/full_objects/SO_5k_pbmc_v3_SCT_none_1_500_10_220617-154349.rds")

so@meta.data$prim_cell_atlas_labels
feature_plot2(so, "prim_cell_atlas_labels")

so@assays$RNA@counts <- scexpr::reverse_lognorm(so)
cdsdbl <- scDblFinder::addDoublets(scexpr::get_layer(so, layer = "counts"),
                                   dbr = ncol(scexpr::get_layer(so, layer = "counts"))/1e4*2,
                                   clusters = so@meta.data$prim_cell_atlas_labels)
cdsdbl@colData@listData$type <- ifelse(cdsdbl@colData@listData$type == "singlet", "single_transcriptome_as_is", "doublet_simul")
table(cdsdbl@colData@listData$type)

ncol(scexpr::get_layer(so, layer = "counts"))*0.008*4
sodbl <- Seurat::CreateSeuratObject(assays(cdsdbl)[["counts"]])



for (i in names(cdsdbl@colData@listData)) {
  sodbl <- Seurat::AddMetaData(object = sodbl,
                               metadata = cdsdbl@colData@listData[[i]],
                               col.name = i)
}

sonew <- scexpr::SO_prep02(SO_unprocessed = list(seurat = sodbl),
                           interactive_varfeat_selection = F,
                           interactive_pc_selection = F,
                           nhvf = 600,
                           npcs = 15)

scf <- scDblFinder::scDblFinder(sce = get_layer(obj = sonew,
                                                assay = "RNA",
                                                layer = "counts"),
                                nfeatures = Seurat::VariableFeatures(sonew),
                                dims = ncol(sonew@reductions$pca@cell.embeddings))
sonew@meta.data$dbl_score <- scf@colData$scDblFinder.score
sonew@meta.data$dbl_class <- scf@colData$scDblFinder.class
sonew@meta.data$dbl_class <- ifelse(sonew@meta.data$dbl_class == "singlet", "singlet_infer", "doublet_calc")
sonew <- Seurat::AddMetaData(sonew, so@meta.data[,setdiff(names(so@meta.data), names(sonew@meta.data))])

feature_plot2(sonew, features = c("type", "cluster", "dbl_score", "dbl_class"))
feature_plot_stat(sonew, meta_col = "dbl_class", features = "dbl_score")
feature_plot_stat(sonew, meta_col = "type", features = "dbl_score")
feature_plot_stat(sonew, meta_col = "cluster", features = "dbl_score", split.by = "type") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  facet_wrap(vars(split_feature))

tab <- table(sonew@meta.data$dbl_class, sonew@meta.data$type)
tab

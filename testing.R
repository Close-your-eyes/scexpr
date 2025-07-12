library(scexpr)
#source("/Users/vonskopnik/Documents/R_packages/scexpr/feature_plot2.R")
#source("/Users/vonskopnik/Documents/R_packages/scexpr/helpers.R")
# source("~/Documents/R_packages/scexpr/R/get_data.R")
# source("~/Documents/R_packages/scexpr/R/feature_plot_data.R")
# source("~/Documents/R_packages/scexpr/R/feature_plot_gene.R")
# source("~/Documents/R_packages/scexpr/R/feature_plot_meta.R")
SO1 <- readRDS("/Users/vonskopnik/Documents/scRNAseq/2019_SLE_LN/data/SO_processed/SO_urine_hg38_CR7_wo_Introns_non_contam_rep1_SCT_harmony_1_800_12_230509-104557.rds")
SO1 <- Seurat::UpdateSeuratObject(SO1)

SO2 <- readRDS("/Users/vonskopnik/Documents/scRNAseq/2019_SLE_LN/data/SO_processed/SO_urine_blood_hg38_CR7_non_contam_rep1_wo_Introns_strict_QC_patient_integration_CD8_eff_SCT_harmony_1_300_12_231016-163715.rds")
SO2 <- Seurat::UpdateSeuratObject(SO2)


feature_plot2(SO = list(SO1, SO2), features = c("CX3CR1", "CD8A"))

data <- get_data(SO = list(SO1, SO2), feature = c("CX3CR1", "CD8A", "clusters"))
attributes(data[[1]])


feature_plot2(SO = list(SO1, SO2), features = c("clusters"), get_data_args = list(qmax = 0.9, qmin = 0.05, label_feature = "clusters"))

feature_plot2(SO = list(SO1), features = c("clusters"), get_data_args = list(qmax = 0.9, qmin = 0.05, label_feature = "clusters_short"),
              contour_feature = "clusters_short", theme = colrr::theme_material(white = T))

feature_plot2(SO = list(SO1), features = c("CX3CR1", "CD8A"), get_data_args = list(qmax = 0.9, qmin = 0.05, label_feature = "clusters_short"),
              contour_feature = "clusters_short", theme = colrr::theme_material(white = F),
              col_legend_args = list(barwidth = 1,
                                     barheight = 10,
                                     override.aes = list(size = 4),
                                     title.theme = ggtext::element_markdown(),
                                     title = "..auto..",
                                     order = 1))

pp <- feature_plot2(SO = SO1, features = c("clusters"), contour_feature = "clusters")
feature_plot2(SO = SO1, features = c("CX3CR1"), contour_feature = "clusters", contour_rm_outlier = T)
pp <- feature_plot2(SO = SO1, features = c("CX3CR1"), contour_feature = "clusters",
                    contour_fun = geomtextpath::geom_labelpath, contour_path_label = "pct",
                    contour_rm_outlier = T)

feature_plot2(SO1, c("CD8A", "CD8B"),
              col_legend_args = list(barwidth = 1,
                                     barheight = 8,
                                     override.aes = list(size = 4),
                                     title.theme = ggtext::element_markdown(),
                                     title = "..auto..",
                                     order = 1,
                                     position = "inside"),
              theme_args = list(axis.ticks = ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(),
                                axis.title.x = ggplot2::element_blank(),
                                axis.title.y = ggplot2::element_blank(),
                                panel.grid = ggplot2::element_blank(),
                                legend.background = ggplot2::element_blank(),
                                legend.key.size = ggplot2::unit(0.3, "cm"),
                                legend.key = ggplot2::element_blank(),
                                legend.position.inside = c(0.2,0.2)),
              col_steps = NULL)

feature_plot2(SO1, c("CD8A", "CD8B"),
              col_legend_args = list(barwidth = 1,
                                     barheight = 8,
                                     override.aes = list(size = 4),
                                     title.theme = ggtext::element_markdown(),
                                     title = "..auto..",
                                     order = 1,
                                     position = "inside"),
              theme_args = list(axis.ticks = ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(),
                                axis.title.x = ggplot2::element_blank(),
                                axis.title.y = ggplot2::element_blank(),
                                panel.grid = ggplot2::element_blank(),
                                legend.background = ggplot2::element_blank(),
                                legend.key.size = ggplot2::unit(0.3, "cm"),
                                legend.key = ggplot2::element_blank(),
                                legend.position.inside = c(0.2,0.2)))

pp + ggplot2::guides(x = ggh4x::guide_axis_truncated(trunc_lower = ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1], trunc_upper = ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1] + abs(min(c(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])) - max(c(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range[2], ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$x.range[1])))/4 ),
                     y = ggh4x::guide_axis_truncated(trunc_lower = ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1], trunc_upper = ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1] + abs(min(c(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])) - max(c(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range[2], ggplot2::ggplot_build(p_blood)$layout$panel_params[[1]]$y.range[1])))/4 ))


names(formals(feature_plot_data)) %in% names(formals(feature_plot2))
names(formals(feature_plot2))[which(!names(formals(feature_plot2)) %in% names(formals(feature_plot_data)))]

#check feature cutoff direction



# BiocManager::install("splatter")
library(splatter)
ncol(SO1)
params <- splatter::splatEstimate(as.matrix(SO1@assays$RNA$counts[SO1@assays[["SCT"]]@var.features, sample.int(ncol(SO1), 2000)]))
params2 <- splatter::kersplatEstimate(as.matrix(SO1@assays$RNA$counts[SO1@assays[["SCT"]]@var.features, sample.int(ncol(SO1), 2000)]))
sim <- splatter::splatSimulate(params = params, method = "groups", group.prob = rep(0.2,5))
sim2 <- splatter::kersplatSimulate(params = params2)
SOy <- prep_SO(SO_unprocessed = list(xxx = Seurat::CreateSeuratObject(counts = sim2@assays@data@listData[["counts"]])),
               reductions = "umap", normalization = "LogNormalize", verbose = T, npcs = 10, nhvf = 400)
feature_plot2(SOy, c("Gene2", "RNA_snn_res.1") ,reduction = "umap", theme = colrr::theme_material())
names(SOy@meta.data)

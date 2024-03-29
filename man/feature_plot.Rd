% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_plot.R
\name{feature_plot}
\alias{feature_plot}
\title{Plot features of single cell transcriptomes on a dimension reduction map}
\usage{
feature_plot(
  SO,
  features,
  assay = c("RNA", "SCT"),
  dims = c(1, 2),
  cells = NULL,
  downsample = 1,
  pt.size = 1,
  pt.size.expr.factor = 1,
  order = T,
  order.abs = T,
  shuffle = T,
  order.rev = F,
  min.q.cutoff = 0,
  max.q.cutoff = 1,
  reduction = "umap",
  split.by = NULL,
  split.by.scales = "fixed",
  shape.by = NULL,
  combine = T,
  ncol.combine = NULL,
  nrow.combine = NULL,
  nrow.inner = NULL,
  ncol.inner = NULL,
  feature.aliases = NULL,
  binary = F,
  title = NULL,
  title.font.size = 14,
  cutoff.feature = NULL,
  cutoff.expression = 0,
  exclusion.feature = NULL,
  legend.position = "right",
  legend.title.text.size = 14,
  legend.text.size = 10,
  legend.barheight = 3,
  legend.barwidth = 0.5,
  legend.nrow = NULL,
  legend.ncol = NULL,
  legend.key.size = 0.3,
  legend.shape.size = 3,
  legend.col.size = 3,
  hide.shape.legend = F,
  font.family = "sans",
  col.pal.c = "spectral",
  col.pal.d = "custom",
  col.excluded.cells = "grey95",
  col.non.expresser = "grey85",
  col.expresser = "tomato2",
  col.pal.rev = F,
  theme = ggplot2::theme_bw(),
  plot.axis.labels = F,
  plot.panel.grid = F,
  plot.freq.title = NULL,
  plot.freq = NULL,
  plot.legend.title = F,
  plot.title = T,
  order.discrete = T,
  freq.position = c(-Inf, Inf),
  freq.font.size = 4,
  plot.cutoff.feature.inset = F,
  inset.position = c(0.75, 0.1),
  inset.size = c(0.15, 0.15),
  strip.font.size = 14,
  strip.selection = NA,
  plot.labels = NULL,
  label.size = 12,
  na.rm = F,
  inf.rm = F,
  bury_NA = T,
  trajectory.slot = NULL,
  trajectory.color = "grey30",
  trajectory.size = 0.75,
  trajectory.linetype = "solid",
  contour_feature = NULL,
  col.pal.contour = "custom",
  contour_args = list(contour_var = "ndensity", breaks = 0.3, linewidth = 1),
  plot.expr.freq.by.contour.group = F,
  use_ggnewscale_for_contour_colors = T,
  expand_limits = list(),
  ...
)
}
\arguments{
\item{SO}{one or more Seurat object(s)}

\item{features}{vector of features to plot (genes or column names in meta data)}

\item{assay}{which assay to get expression data from}

\item{dims}{which dimensions of the selected dimension reduction to plot}

\item{cells}{a vector of cell names to include (not selected ones are plotted with color col.excluded.cells)}

\item{downsample}{downsample the number of cells (intended to speed up test plottings)}

\item{pt.size}{dot size per cells}

\item{pt.size.expr.factor}{factor of increased dot size for expressing cells}

\item{order}{for meta.col: remains T if var is continuous but becomes F if var is integer (~probably discrete)}

\item{order.abs}{do use absolute values for ordering (any value away from zero (+/) is treated equally)}

\item{shuffle}{do shuffle if order if FALSE; allows to define a definite order for plotting if set to to F}

\item{order.rev}{reverse the ordering to have lowest values on top (or zeros if order.abs = T)}

\item{min.q.cutoff}{decimal number (> 1, < 0) of lower quantile limit where to cut the color scale, intended to squish extremes}

\item{max.q.cutoff}{decimal number (> 1, < 0) of upper quantile limit where to cut the color scale, intended to squish extremes}

\item{reduction}{which reduction to use for plotting}

\item{split.by}{column in meta data to use to split plots}

\item{shape.by}{column in meta data to shape dots (cells) b<}

\item{combine}{combine multiple features to one plot (TRUE) or return a list with one entry per feature (FALSE)}

\item{ncol.combine}{number of columns in combined graphic (feature combined)}

\item{nrow.combine}{number of rows in combined graphic (feature combined)}

\item{nrow.inner}{number of plots per row within one feature (originating from multiple SO and/or split.by)}

\item{ncol.inner}{number of plots per column within one feature (originating from multiple SO and/or split.by)}

\item{feature.aliases}{vector aliases for features; e.g. c("FAIM3" = "FCMR", "MS4A1" = "CD20")}

\item{title}{force title for the plot instead of feature name which is the default}

\item{title.font.size}{font size of the title}

\item{cutoff.feature}{select a feature for cells to consider in plotting}

\item{cutoff.expression}{select the cutoff feature expression levels for filtering}

\item{legend.position}{values between 0 and -1}

\item{trajectory.slot}{name of Misc slot that contains information on how to plot a trajectory}

\item{expand_limits}{}
}
\description{
Plot features of single cell transcriptomes on a dimension reduction map
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_plot.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Title}
\usage{
volcano_plot(
  SO,
  assay = "RNA",
  volcano.data = NULL,
  negative.group.cells,
  positive.group.cells,
  negative.group.name = "negative.group",
  positive.group.name = "positive.group",
  x.axis.symmetric = T,
  x.axis.extension = 0,
  y.axis.pseudo.log = F,
  pseudo.log.sigma = 1,
  inf.fc.shift = 2,
  pt.size = 1,
  pt.alpha = 0.8,
  font.size = 14,
  pval.tick = 0.01,
  min.pct = 0.1,
  max.iter = 10000,
  max.overlaps = 50,
  label.neg.pos.sep = T,
  label.col = "black",
  label.face = "bold",
  font.family = "Courier",
  label.size = 4,
  labels.topn = 30,
  label.features = NULL,
  topn.metric = "p.value",
  nudge.x = 0,
  nudge.y = 0,
  p.plot = "adj.p.val",
  p.adjust = "bonferroni",
  p.cut = NA,
  p.signif = 0.001,
  fc.cut = NA,
  features.exclude = NULL,
  save.path.interactive = NULL,
  ...
)
}
\arguments{
\item{...}{}
}
\value{

}
\description{
Title
}
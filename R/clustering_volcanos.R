#' Title
#'
#' @param SO
#' @param cluster.meta.cols
#' @param reduction
#' @param clustree.prefix
#' @param save.path
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; e.g. mc.cores may be passed when parallel::mclapply is chosen above
#' @param assay
#'
#' @return
#' @export
#'
#' @examples
clustering_volcanos <- function(SO,
                                assay = c("RNA", "SCT"),
                                cluster.meta.cols,
                                reduction = "tSNE",
                                clustree.prefix = NULL,
                                save.path = NULL,
                                lapply_fun = lapply,
                                ...) {
  # meta.cols - keeping additional columns

  if (!requireNamespace("clustree", quietly = T)) {
    utils::install.packages("clustree")
  }

  if (missing(cluster.meta.cols)) {
    stop("Please provide cluster.meta.cols from SO@meta.data.")
  }

  SO <- .check.SO(SO, length = 1)
  reduction <- .check.reduction(SO, reduction = reduction)
  cluster.meta.cols <- .check.features(SO, features = cluster.meta.cols, rownames = F)
  assay <- match.arg(assay, c("RNA", "SCT"))
  lapply_fun <- match.fun(lapply_fun)

  if (is.null(save.path)) {
    print(paste0("save.path not provided, set to getwd(): ", getwd()))
    save.path <- getwd()
  }
  dir.create(save.path, showWarnings = FALSE, recursive = TRUE)

  if (is.null(clustree.prefix) && length(cluster.meta.cols) > 1) {
    splits <- strsplit(cluster.meta.cols, "")
    max_char <- min(nchar(cluster.meta.cols))
    equal <- T
    i <- 0
    while (equal && i < max_char) {
      i <- i + 1
      equal <- nlevels(as.factor(sapply(splits, "[", i))) == 1 && all(suppressWarnings(is.na(as.numeric(sapply(splits, "[", i)))))
    }
    clustree.prefix <- substr(cluster.meta.cols[1], 1, i-1)
    print(paste0("Tried to assume clustree.prefix: '", clustree.prefix, "'. Provide manually if wrong."))
  } else if (is.null(clustree.prefix) && length(cluster.meta.cols) == 1) {
    stop("clustree.prefix cannot be guessed with one meta.col only.")
  }

  # remove not relevant cluster.meta.cols for clustree plot
  SO@meta.data <- SO@meta.data[,c(intersect(which(grepl(clustree.prefix, names(SO@meta.data))), which(grepl(paste(cluster.meta.cols, collapse = "|"), names(SO@meta.data))))),drop=F]
  if (length(cluster.meta.cols) > 1) {
    clust <- clustree::clustree(SO,
                                prefix = clustree.prefix,
                                exprs = "data",
                                show_axis = T,
                                node_colour = "grey80")
  } else {
    clust <- NULL
  }

  data_list <- lapply(cluster.meta.cols, function(y) {
    comb <- expand.grid(unique(SO@meta.data[,y]),unique(SO@meta.data[,y]))
    comb <- dplyr::arrange(comb, Var1, Var2)
    comb <- comb[which(!duplicated(apply(comb, 1, function(x) paste(sort(c(x[1], x[2])), collapse = "")))),]
    #comb <- comb[which(comb$Var1 != comb$Var2),]
    comb$Var1 <- as.character(comb$Var1)
    comb$Var2 <- as.character(comb$Var2)

    vd <- lapply_fun(split(comb, 1:nrow(comb)), function(x) {
      if (x[1,1] == x[1,2]) {
        pgc <- rownames(SO@meta.data[,y,drop=F][which(SO@meta.data[,y] != x[1,2]),,drop=F])
        pgn <- "all_other"
      } else {
        pgc <- rownames(SO@meta.data[,y,drop=F][which(SO@meta.data[,y] == x[1,2]),,drop=F])
        pgn <- x[1,2]
      }
      vd <- volcano_plot(SO,
                         assay = assay,
                         negative.group.cells = rownames(SO@meta.data[,y,drop=F][which(SO@meta.data[,y] == x[1,1]),,drop=F]),
                         positive.group.cells = pgc,
                         negative.group.name = x[1,1],
                         positive.group.name = pgn,
                         interactive.only = T)
      return(vd)
    }, ...)



    all_feats <- unique(unlist(sapply(sapply(vd, "[", "vd"), rownames)))
    #cell_idents <- stats::setNames(rownames(SO@meta.data), nm = as.character(SO@meta.data[,y]))
    #vd <- sapply(vd, function(x) x[-which(names(x) %in% c("non.aggr.data", "negative.group.cells", "positive.group.cells"))], simplify = F)

    comb[which(comb$Var1 == comb$Var2),"Var2"] <- "all_other"
    vd <- c(vd, list(scexpr::feature_plot(SO, features = y, reduction = reduction, plot.title = F)), list(all_feats = all_feats))
    names(vd) <- c(apply(comb, 1, function(x) paste0(x[1], " vs ", paste(x[2]))), "reduction_plot", "all_feats")
    return(vd)
  })

  SO <- Seurat::DietSeurat(SO, assays = assay, counts = F, features = unique(unlist(sapply(data_list, "[", "all_feats"))))
  SO@meta.data <- SO@meta.data[,cluster.meta.cols,drop=F]
  data_list <- sapply(data_list, function(x) x[-which(names(x) %in% c("all_feats"))], simplify = F)

  data_list <- c(data_list, list(clust), list(SO))
  names(data_list) <- c(cluster.meta.cols, "clustree_plot", "Seurat_object")
  file.copy(from = system.file("extdata", "app.R", package = "scexpr"), to = file.path(save.path, "app.R"))
  saveRDS(data_list, file = file.path(save.path, "data.rds"))

  return(data_list)
}

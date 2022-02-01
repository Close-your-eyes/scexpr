#' Title
#'
#' @param SO
#' @param meta.cols
#' @param reduction
#' @param clustree.prefix
#' @param save.path
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param ... additional argument to the lapply function; e.g. mc.cores may be passed when parallel::mclapply is chosen above
#'
#' @return
#' @export
#'
#' @examples
clustering_volcanos <- function(SO,
                                meta.cols,
                                reduction = "tSNE",
                                clustree.prefix = NULL,
                                save.path = NULL,
                                lapply_fun = lapply,
                                ...) {

  library(clustree) # library is required!! https://github.com/lazappi/clustree/issues/14
  if (missing(meta.cols)) {
    stop("Please provide meta.cols from SO@meta.data.")
  }

  SO <- .check.SO(SO, max.length = 1)
  reduction <- .check.reduction(SO, reduction = reduction)
  meta.cols  <- .check.features(SO, features = meta.cols, rownames = F)
  lapply_fun <- match.fun(lapply_fun)

  if (is.null(save.path)) {
    print(paste0("save.path not provided, set to getwd(): ", getwd()))
    save.path <- getwd()
  }
  dir.create(save.path, showWarnings = FALSE, recursive = TRUE)

  if (is.null(clustree.prefix) && length(meta.cols) > 1) {
    splits <- strsplit(meta.cols, "")
    max_char <- min(nchar(meta.cols))
    equal <- T
    i <- 0
    while (equal && i < max_char) {
      i <- i + 1
      equal <- nlevels(as.factor(sapply(splits, "[", i))) == 1 && all(suppressWarnings(is.na(as.numeric(sapply(splits, "[", i)))))
    }
    clustree.prefix <- substr(meta.cols[1], 1, i-1)
    print(paste0("Tried to assume clustree.prefix: '", clustree.prefix, "'. Provide manually if wrong."))
  } else if (is.null(clustree.prefix) && length(meta.cols) == 1) {
    stop("clustree.prefix cannot be guessed with one meta.col only.")
  }

  # remove not relevant meta.cols for clustree plot
  SO@meta.data <- SO@meta.data[,c(intersect(which(grepl(clustree.prefix, names(SO@meta.data))), which(grepl(paste(meta.cols, collapse = "|"), names(SO@meta.data))))),drop=F]
  if (length(meta.cols) > 1) {
    clust <- clustree::clustree(SO,
                                prefix = clustree.prefix,
                                exprs = "data",
                                show_axis = T,
                                node_colour = "grey80")
  } else {
    clust <- NULL
  }


  data_list <- lapply(meta.cols, function(y) {
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

      volcano_plot(SO,
                   negative.group.cells = rownames(SO@meta.data[,y,drop=F][which(SO@meta.data[,y] == x[1,1]),,drop=F]),
                   positive.group.cells = pgc,
                   negative.group.name = x[1,1],
                   positive.group.name = pgn,
                   interactive.only = T)
    }, ...)

    all_feats <- unique(unlist(sapply(sapply(vd, "[", "non.aggr.data"), rownames)))
    #all_cells <- unique(unlist(sapply(sapply(vd, "[", "non.aggr.data"), colnames)))
    vd <- sapply(vd, function(x) x[-which(names(x) == "non.aggr.data")], simplify = F)

    comb[which(comb$Var1 == comb$Var2),"Var2"] <- "all_other"
    vd <- c(vd, list(scexpr::feature_plot(SO, features = y, reduction = reduction, plot.title = F)), list(all_feats = all_feats))
    names(vd) <- c(apply(comb, 1, function(x) paste0(x[1], " vs ", paste(x[2]))), "reduction_plot", "all_feats")
    return(vd)
  })


  non.aggr.data <- Seurat::GetAssayData(SO, slot = "data", assay = "RNA")[unique(unlist(sapply(data_list, "[", "all_feats"))), ] #unique(unlist(sapply(data_list, "[", "all_cells")))
  data_list <- sapply(data_list, function(x) x[-which(names(x) %in% c("all_feats", "all_cells"))], simplify = F)

  data_list <- c(data_list, list(clust), list(non.aggr.data))
  names(data_list) <- c(meta.cols, "clustree_plot", "non.aggr.data")
  file.copy(from = system.file("extdata/clustering_volcanos", "app.R", package = "scexpr"), to = file.path(save.path, "app.R"))
  saveRDS(data_list, file = file.path(save.path, "data.rds"))
  return(data_list)
}

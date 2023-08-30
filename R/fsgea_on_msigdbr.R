#' Convenient wrapper around fgsea function
#'
#' @param gene.ranks this becomes stats in fgsea_fun
#' @param gene.sets this becomes pathways in fgsea_fun
#' @param min.padj
#' @param use.msigdbr
#' @param msigdbr_args
#' @param fgsea_fun
#' @param fgsea_args
#' @param return.gene.sets set to FALSE in order to save memory when the same gene sets are used multiple times
#' @param return.subset.gene.sets return the subset of each gene set that is actually found in gene.ranks
#' @param ... arguments to .plotEnrichment
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' # get gene sets from msigdbr, split to a named list
#' gene.sets <- scexpr:::.get.split.msigdbr()
#' # provide gene.sets manually and leave use.msigdbr = F in case fgsea_on_msigdbr is called many times with the same gene.sets
#' }
fgsea_on_msigdbr <- function(gene.ranks = NULL,
                             gene.sets = NULL,
                             return.subset.gene.sets = T,
                             return.gene.sets = T,
                             min.padj = 0, # which plots to generate; by default: none
                             use.msigdbr = F,
                             msigdbr_args = list(species = "Homo sapiens", category = NULL, subcategory = NULL),
                             fgsea_fun = fgsea::fgseaMultilevel,
                             fgsea_args = list(stats = gene.ranks, pathways = gene.sets),
                             ...) {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    BiocManager::install("fgsea")
  }
  if (is.null(gene.ranks)) {
    stop("gene.ranks has to be provided.")
  }
  if (!is.numeric(gene.ranks)) {
    stop("gene.ranks has to be a numeric metric.")
  }
  if (is.null(names(gene.ranks))) {
    stop("gene.ranks has to have names which should be genes.")
  }
  if (!is.null(gene.sets)) {
    if (!is.list(gene.sets) || is.null(names(gene.sets))) {
      stop("gene.sets has to be a named list of gene sets.")
    }
  }
  fgsea_fun <- match.fun(fgsea_fun)

  if (use.msigdbr) {
    if (!requireNamespace("msigdbr", quietly = T)) {
      utils::install.packages("msigdbr")
    }
    gene.sets.list <- c(gene.sets, .get.split.msigdbr(msigdbr_args = msigdbr_args))
    gene.sets <- gene.sets.list[["sets"]]
  }
  if (length(gene.sets) == 0) {
    stop("No gene.sets provided and use.msigdbr=F. Either manually provide gene.sets or set use.msigdbr=T to select sets from their.")
  }
  names(gene.sets) <- make.unique(names(gene.sets))

  if (!"pathways" %in% names(fgsea_args)) {
    fgsea_args <- c(list(pathways = gene.sets), fgsea_args)
  }
  if (!"stats" %in% names(fgsea_args)) {
    fgsea_args <- c(list(stats = gene.ranks), fgsea_args)
  }
  results <- as.data.frame(Gmisc::fastDoCall(fgsea_fun, args = fgsea_args))
  if (use.msigdbr) {
    results$pathway_cat <- gene.sets.list[["cats"]][results$pathway]
  }

  results$leadingEdge_sorted_chr <- purrr::map_chr(purrr::map(results$leadingEdge, sort), paste, collapse = ",")
  results$leadingEdge_size <- lengths(results$leadingEdge)
  results$leadingEdge_size_rel <- results$leadingEdge_size/results$size
  results$leadingEdge_rank <- purrr::map_int(results$pathway, function(x) ifelse(results[["ES"]][[which(results[["pathway"]] == x)]] > 0,
                                                                                 max(match(results[["leadingEdge"]][[which(results[["pathway"]] == x)]], names(rev(gene.ranks)))),
                                                                                 min(match(results[["leadingEdge"]][[which(results[["pathway"]] == x)]], names(rev(gene.ranks))))))
  gsea_plots <- lapply(.self_naming(results[which(results$padj <= min.padj),"pathway",drop=T]), function(x) {
    gsea_data <- prep_gsea(name = x,
                           gene.set = gene.sets[[x]],
                           gene.ranks = gene.ranks,
                           gseaParam = 1,
                           pval = signif(results[which(results$pathway == x), "padj"], 2),
                           NES = signif(results[which(results$pathway == x), "NES"], 2))
    # leadingEdge_rank differs by 1 between methods ?!
    # leadingEdge_rank = results[which(results$pathway == x), "leadingEdge_rank"]

    return(list(plot = plot_gsea(data = gsea_data), data = gsea_data))
  })


  return(list(data = results,
              plots = gsea_plots,
              gene.sets = if (return.gene.sets) {gene.sets} else {NULL},
              gene.sets.subset = if(return.subset.gene.sets) {lapply(gene.sets, function(x) intersect(x, names(gene.ranks)))} else {NULL},
              gene.ranks = gene.ranks))
}




.get.split.msigdbr <- function(msigdbr_args = list(species = "Homo sapiens", category = NULL, subcategory = NULL)) {
  sets <- Gmisc::fastDoCall(msigdbr::msigdbr, args = msigdbr_args)[,c("gs_name", "gene_symbol", "gs_cat")]
  sets <- unique(sets) # some gene symbols exist in duplicates, but ENSEMBLE differs; make gene symbols unique
  sets_cat <- unique(sets[,which(names(sets) %in% c("gs_cat", "gs_name"))])
  return(list(sets = split(sets$gene_symbol, sets$gs_name), cats = stats::setNames(sets_cat$gs_cat, sets_cat$gs_name)))
}



.self_naming <- function(x) {
  return(stats::setNames(x, x))
}

prep_gsea <- function(name,
                      gene.set,
                      gene.ranks,
                      gseaParam = 1,
                      pval,
                      NES) {
  data <- .plotEnrichmentData(pathway = gene.set,
                              stats = gene.ranks,
                              gseaParam = gseaParam)

  colorbar_df <- .prep_gsea_colorbar(x = as.data.frame(data[["stats"]]))
  rank_df <- data.frame(gene = names(rev(gene.ranks)[data[["ticks"]][["rank"]]]), rank = data[["ticks"]][["rank"]])

  leadingEdge_rank <- as.numeric(data[["curve"]][which(data[["curve"]]$ES == c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))]), "rank"])
  if (abs(data[["posES"]]) > abs(data[["negES"]])) {
    rank_df$leadingEdge <- rank_df$rank <= leadingEdge_rank
  } else if (abs(data[["posES"]]) < abs(data[["negES"]])) {
    rank_df$leadingEdge <- rank_df$rank >= leadingEdge_rank
  } else {
    stop("posES equal to negES. What now?")
  }
  leadingEdge_size = sum(rank_df$leadingEdge)



  data <- c(data, list(name = name,
                       rank_df = rank_df,
                       colorbar_df = colorbar_df,
                       leadingEdge_rank = leadingEdge_rank,
                       leadingEdge_size = leadingEdge_size,
                       pval = pval,
                       ES = c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))],
                       NES = NES))
  return(data)
}

.prep_gsea_colorbar <- function(x) {
  browser()
  zscore_cuts <- NULL
  x$stat_scale <- scale(x$stat)
  if (is.null(zscore_cuts)) {
    zscore_cuts <- seq(floor(min(x$stat_scale)), ceiling(max(x$stat_scale)), 1)


    unique(c(seq(min(zscore_cuts), max(zscore_cuts), 4), max(zscore_cuts)))

    x <- unique(c(min(zscore_cuts), seq(min(zscore_cuts), 0, 2), 0))
    y <- unique(c(min(zscore_cuts), seq(min(x), 0, 3), 0))

    z <- unique(c(seq(0, max(zscore_cuts), 2), max(zscore_cuts)))
    z1 <- unique(c(0, seq(0, max(z), 3), max(z)))

    split_func(zscore_cuts[which(zscore_cuts <= 0)],4)
  }
  data <-
    x %>%
    dplyr::mutate(zscore = cut(stat_scale, breaks = zscore_cuts)) %>% # as.numeric(as.factor())
    dplyr::group_by(zscore) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank))
  data$zscore <- factor(data$zscore, levels = data %>% dplyr::arrange(min_rank) %>% dplyr::pull(zscore))
  data <- dplyr::arrange(data, zscore)
  data$fill_col <- ""

  # manually take care, that positive stats become reddish color and negative bluish
  pos_col <- RColorBrewer::brewer.pal(10, "RdBu")[1:5]
  neg_col <- RColorBrewer::brewer.pal(10, "RdBu")[6:10]
  data[which(grepl("^\\(-", data$zscore)),"fill_col"] <- grDevices::colorRampPalette(neg_col, interpolate = "linear")(length(which(grepl("^\\(-", data$zscore))))
  data[which(grepl("^\\([[:digit:]]", data$zscore)),"fill_col"] <- grDevices::colorRampPalette(pos_col, interpolate = "linear")(length(which(grepl("^\\([[:digit:]]", data$zscore))))


  return(list(data = data, breaks = zscore_cuts))
}

.plotEnrichmentData <- function(pathway,
                                stats,
                                gseaParam=1) {
  # copied from fgsea pkg

  if (any(!is.finite(stats))){
    stop("Not all stats values are finite numbers")
  }

  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)

  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  pathway <- unique(pathway)

  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.table::data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
  ticks <- data.table::data.table(rank=pathway, stat=statsAdj[pathway])
  stats <- data.table::data.table(rank=seq_along(stats), stat=statsAdj)

  res <- list(
    curve=toPlot,
    ticks=ticks,
    stats=stats,
    posES=max(tops),
    negES=min(bottoms),
    spreadES=max(tops)-min(bottoms),
    maxAbsStat=max(abs(statsAdj)))
}

split_func <- function(x, by) {
  r <- diff(range(x))
  out <- seq(0, r - by - 1, by = by)
  c(round(min(x) + c(0, out - 0.51 + (max(x) - max(out)) / 2), 0), max(x))
}

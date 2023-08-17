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
                             min.padj = 0.001, # which plots to generate
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
    fgsea_args <- c(list(stats = gene.sets), fgsea_args)
  }
  results <- as.data.frame(Gmisc::fastDoCall(fgsea_fun, args = fgsea_args))
  if (use.msigdbr) {
    results$pathway_cat <- gene.sets.list[["cats"]][results$pathway]
  }

  results$leadingEdge_sorted_chr <- purrr::map_chr(purrr::map(results$leadingEdge, sort), paste, collapse = ",")
  results$leadingEdge_size <- lengths(results$leadingEdge)
  results$leadingEdge_size_rel <- results$leadingEdge_size/results$size
  results$leadingEdge_rank <- purrr::map_int(results$pathway, function(x) max(match(results[["leadingEdge"]][[which(results[["pathway"]] == x)]], names(rev(gene.ranks)))))

  gsea_plots <- lapply(.self_naming(results[which(results$padj <= min.padj),"pathway",drop=T]), function(x) {
    .plotEnrichment(pathway = gene.sets[[x]],
                    pathway_name = x,
                    p_val = signif(results[which(results$pathway == x), "padj"], 2),
                    ES = signif(results[which(results$pathway == x), "ES"], 2),
                    NES = signif(results[which(results$pathway == x), "NES"], 2),
                    stats = gene.ranks,
                    color_ES_line = "black",
                    color_min_max_ES_line = "black",
                    leadingEdge_rank = results[which(results$pathway == x), "leadingEdge_rank"],
                    ...)
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



.plotEnrichment <- function(pathway,
                            pathway_name = NULL,
                            p_val = "",
                            ES = "",
                            NES = "",
                            stats,
                            gseaParam=1,
                            ticksSize=0.2,
                            color_ES_line = "black",
                            color_min_max_ES_line = "black",
                            color_leadingEdge = "hotpink2",
                            leadingEdge_rank = NULL,
                            leadingEdge_size = T,
                            label_genes = NULL, # leave NULL be default because it may take some time in case many gene sets are tested
                            theme = ggplot2::theme_bw()) {

  pd <- .plotEnrichmentData(
    pathway = pathway,
    stats = stats,
    gseaParam = gseaParam)


  stats_df <-
    dplyr::bind_rows(list(as.data.frame(pd[["stats"]]) %>%
                            dplyr::filter(stat > 0) %>%
                            dplyr::mutate(group = as.numeric(as.factor(cut(rank,4)))),
                          as.data.frame(pd[["stats"]]) %>%
                            dplyr::filter(stat < 0) %>%
                            dplyr::mutate(group = -as.numeric(as.factor(cut(rank,4)))),
                          as.data.frame(pd[["stats"]]) %>%
                            dplyr::filter(stat == 0) %>%
                            dplyr::mutate(group = 0))) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank))
  stats_df$group <- factor(stats_df$group, levels = stats_df %>% dplyr::arrange(min_rank) %>% dplyr::pull(group))
  stats_df <- dplyr::arrange(stats_df, group)
  stats_df$fill_col <- RColorBrewer::brewer.pal(nrow(stats_df), "RdBu")

  p2 <- ggplot2::ggplot(as.data.frame(pd[["stats"]]), ggplot2::aes(x = rank, y = stat)) +
    ggplot2::geom_col() +
    theme +
    ggplot2::labs(x="gene rank", y="ranking metric")

  p <- with(pd,
            p <- ggplot2::ggplot(data=curve) +
              ggplot2::geom_line(ggplot2::aes(x=rank, y=ES), color = color_ES_line) +
              ggplot2::geom_segment(data=ticks,
                                    mapping=ggplot2::aes(x=rank, y=-spreadES/10,
                                                         xend=rank, yend=spreadES/10),
                                    linewidth=ticksSize) +
              ggplot2::geom_hline(yintercept=posES, colour=color_min_max_ES_line, linetype="dashed") +
              ggplot2::geom_hline(yintercept=negES, colour=color_min_max_ES_line, linetype="dashed") +
              ggplot2::geom_hline(yintercept=0, colour="black") +
              theme +
              ggplot2::labs(x="gene rank", y="enrichment score (ES)")
  )
  p <- p + ggplot2::geom_rect(data = stats_df, aes(xmin = min_rank, xmax = max_rank, fill = I(fill_col)), ymin = -pd[["spreadES"]]/16, ymax = pd[["spreadES"]]/16, alpha = 0.8)

  le_ranks <- pd[["ticks"]][["rank"]][which(pd[["ticks"]][["rank"]] <= leadingEdge_rank)]
  le_gene_df <- data.frame(gene = names(rev(stats)[le_ranks]), x = le_ranks)

  all_ranks <- pd[["ticks"]][["rank"]]
  all_gene_df <- data.frame(gene = names(rev(stats)[all_ranks]), x = all_ranks)
  if (!is.null(label_genes)) {
    p <- p + ggrepel::geom_text_repel(data = if (label_genes == "leadingEdge") le_gene_df else all_gene_df[which(all_gene_df$gene %in% label_genes)], ggplot2::aes(label = gene, x = x, y = 0),
                                      max.overlaps = length(le_ranks), max.time = 5, nudge_y = pd[["spreadES"]]/2.5, nudge_x = leadingEdge_rank*6,
                                      segment.size = 0.2, segment.color = "grey80")
  }

  if (!is.null(leadingEdge_rank)) {
    p <- p + ggplot2::geom_segment(x = leadingEdge_rank,
                                   xend = leadingEdge_rank,
                                   y = 0,
                                   yend = pd[["posES"]],
                                   color = color_leadingEdge,
                                   linetype="dashed")
    if (leadingEdge_size) {
      p <- p + ggplot2::annotate("text",
                                 label = paste0("n = ", which(pd[["ticks"]]$rank == leadingEdge_rank)),
                                 color = color_leadingEdge,
                                 y = pd[["spreadES"]]/12,
                                 x = leadingEdge_rank*1.1,
                                 hjust = 0)
    }
  }

  p <- p + ggplot2::labs(title = paste0(pathway_name, " (n = ", nrow(pd[["ticks"]]), ")"), subtitle = paste0("p = ", p_val, "\nES = ", ES, "\nNES = ", NES))

  return(list(plot = p, metric_plot = p2, data = c(pd, list(pathway_name = pathway_name, ES = ES, NES = NES, p_val = p_val,
                                                            leadingEdge_rank = leadingEdge_rank,
                                                            leadingEdge_size = which(pd[["ticks"]]$rank == leadingEdge_rank),
                                                            leadingEdge_genes = le_gene_df,
                                                            geneSet_genes = all_gene_df,
                                                            color_stats_df = stats_df))))
}

.plotEnrichmentData <- function(pathway, stats,
                                gseaParam=1) {

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


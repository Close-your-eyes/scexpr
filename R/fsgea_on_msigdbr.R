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
    plot_gsea(data = gsea_data)
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
  data <- fgsea::plotEnrichmentData(pathway = gene.set,
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

plot_gsea <- function(data,
                      ticksSize=0.2,
                      color_ES_line = "black",
                      color_min_max_ES_line = "black",
                      color_leadingEdge = "hotpink2",
                      label_genes = NULL, # leave NULL be default because it may take some time in case many gene sets are tested
                      theme = ggplot2::theme_bw(),
                      plot_leadingEdge_rank = T,
                      plot_leadingEdge_size = T) {

  p <- with(data,
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
  p <- p + ggplot2::geom_rect(data = data[["colorbar_df"]], aes(xmin = min_rank,
                                                   xmax = max_rank,
                                                   fill = I(fill_col)),
                              ymin = -data[["spreadES"]]/16,
                              ymax = data[["spreadES"]]/16, alpha = 0.8)


  if (!is.null(label_genes)) {
    p <- p + ggrepel::geom_text_repel(data = if (label_genes == "leadingEdge") le_gene_df else all_gene_df[which(all_gene_df$gene %in% label_genes)], ggplot2::aes(label = gene, x = x, y = 0),
                                      max.overlaps = length(le_ranks),
                                      max.time = 5,
                                      nudge_y = data[["spreadES"]]/2.5,
                                      nudge_x = data[["leadingEdge_rank"]]*6,
                                      segment.size = 0.2,
                                      segment.color = "grey80")
  }

  if (plot_leadingEdge_rank) {
    p <- p + ggplot2::geom_segment(x = data[["leadingEdge_rank"]],
                                   xend = data[["leadingEdge_rank"]],
                                   y = 0,
                                   yend = c(data[["posES"]], data[["negES"]])[which.max(abs(c(data[["posES"]], data[["negES"]])))],
                                   color = color_leadingEdge,
                                   linetype="dashed")
    if (plot_leadingEdge_size) {
      p <- p + ggplot2::annotate("text",
                                 label = paste0("n = ", data[["leadingEdge_size"]]), #  which(data[["ticks"]]$rank == data[["leadingEdge_rank"]])
                                 color = color_leadingEdge,
                                 y = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["spreadES"]]/8, -data[["spreadES"]]/8),
                                 x = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), data[["leadingEdge_rank"]]*1.1, data[["leadingEdge_rank"]]*0.98),
                                 hjust = ifelse(abs(data[["posES"]]) > abs(data[["negES"]]), 0, 1))
    }
  }
  p <- p + ggplot2::labs(title = paste0(data[["name"]], " (n = ", nrow(data[["ticks"]]), ")"), subtitle = paste0("p = ", data[["pval"]], "\nES = ", signif(data[["ES"]], 2), "\nNES = ", data[["NES"]]))

  p_metric <- ggplot2::ggplot(as.data.frame(data[["stats"]]),
                              ggplot2::aes(x = rank, y = stat)) +
    ggplot2::geom_col() +
    theme +
    ggplot2::labs(x="gene rank", y="ranking metric")
  return(list(plot = p, metric_plot = p_metric))
}

.prep_gsea_colorbar <- function(x,
                                cut_level_pos_stat = 4,
                                cut_level_neg_stat = 4) {
  data <-
    dplyr::bind_rows(list(x %>%
                            dplyr::filter(stat > 0) %>%
                            dplyr::mutate(group = as.numeric(as.factor(cut(rank,cut_level_pos_stat)))),
                          x %>%
                            dplyr::filter(stat < 0) %>%
                            dplyr::mutate(group = -as.numeric(as.factor(cut(rank,cut_level_neg_stat)))),
                          x %>%
                            dplyr::filter(stat == 0) %>%
                            dplyr::mutate(group = 0))) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(min_rank = min(rank), max_rank = max(rank))
  data$group <- factor(data$group, levels = data %>% dplyr::arrange(min_rank) %>% dplyr::pull(group))
  data <- dplyr::arrange(data, group)

  ## add interpolation in case more levels than colors in brewer.pal
  data$fill_col <- RColorBrewer::brewer.pal(nrow(data), "RdBu")
  return(data)
}

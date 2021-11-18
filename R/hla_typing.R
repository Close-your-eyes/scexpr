hla_typing <- function(hla_ref,
                       reads,
                       diff_allele_hits_single_results = 5,
                       lapply_fun = lapply) {

  lapply_fun <- match.fun(lapply_fun)

  single_res <- Biostrings::vcountPDict(subject = Biostrings::DNAStringSet(hla_ref$seq_Exon2_3), pdict = Biostrings::PDict(reads$seq))
  colnames(single_res) <- hla_ref$allele # ncol(result) - hla_ref
  rownames(single_res) <- reads$readName # nrow(result) - reads
  n_noHit <- sum(Matrix::rowSums(single_res) == 0)
  n_Hit <- sum(Matrix::rowSums(single_res) > 0)

  top_single_res <- single_res[Matrix::rowSums(single_res) > 0, Matrix::colSums(single_res) >= max(Matrix::colSums(single_res))/diff_allele_hits_single_results]
  top_single_res_df <-
    data.frame(n_Hit = Matrix::colSums(top_single_res)) %>%
    tibble::rownames_to_column("allele") %>%
    dplyr::left_join(hla_ref, by = "allele") %>%
    dplyr::mutate(rank = dplyr::row_number(-n_Hit))

  col.combs <- utils::combn(ncol(top_single_res), 2, simplify = F)
  print(paste0("Calculating top pairwise matches. Combinations: ", length(col.combs), "."))

#lapply_fun
  pairwise_results <- lapply_fun(col.combs, function(x) {
    c(sum(top_single_res[,x[1]] + top_single_res[,x[2]] == 1),
      sum(top_single_res[,x[1]] + top_single_res[,x[2]] == 2))
  }, ...)

  pairwise_results_df <-
    data.frame(unique_explained_reads = sapply(pairwise_results, "[", 1), double_explained_reads = sapply(pairwise_results, "[", 2), allele.1 = colnames(top_single_res)[sapply(col.combs, "[", 1)], allele.2 = colnames(top_single_res)[sapply(col.combs, "[", 2)]) %>%
    dplyr::mutate(total_explained_reads = unique_explained_reads + double_explained_reads)  %>%
    dplyr::mutate(total_explained_reads_rank = dplyr::dense_rank(-total_explained_reads), unique_explained_reads_rank = dplyr::dense_rank(-unique_explained_reads)) %>%
    dplyr::mutate(rank.sum = total_explained_reads_rank + unique_explained_reads_rank)

  top_pairwise_results_df <-
    dplyr::bind_rows(pairwise_results_df %>% dplyr::top_n(-200, unique_explained_reads_rank), pairwise_results_df %>% dplyr::top_n(-200, total_explained_reads_rank)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(combined.rank = dplyr::dense_rank(base::interaction(-total_explained_reads, -unique_explained_reads, lex.order = TRUE)))

  top_pairwise_results_df$allele.1.2 <- apply(top_pairwise_results_df[,c("allele.1", "allele.2")], 1, function(x) paste(sort(c(x[1], x[2])), collapse = "_"))

  top_pairwise_results_df <-
    top_pairwise_results_df %>%
    dplyr::group_by(allele.1.2) %>%
    dplyr::filter(combined.rank == min(combined.rank)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(allele.1.2, .keep_all = T) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(allele, p_group, g_group), by = c("allele.1" = "allele")) %>% dplyr::rename("p_group.1" = p_group, "g_group.1" = g_group) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(allele, p_group, g_group), by = c("allele.2" = "allele")) %>% dplyr::rename("p_group.2" = p_group, "g_group.2" = g_group) %>%
    dplyr::left_join(top_single_res_df %>% dplyr::select(allele, n_Hit), by = c("allele.1" = "allele")) %>% dplyr::rename("explained_reads_allele.1" = n_Hit) %>%
    dplyr::left_join(top_single_res_df %>% dplyr::select(allele, n_Hit), by = c("allele.2" = "allele")) %>% dplyr::rename("explained_reads_allele.2" = n_Hit) %>%
    dplyr::mutate(n_Hit = n_Hit) %>%
    dplyr::mutate(n_noHit = n_noHit) %>%
    dplyr::mutate(explained.reads.diff = abs(explained_reads_allele.1 - explained_reads_allele.2)) %>%
    dplyr::mutate(percent_matching_reads_explained = total_explained_reads/n_Hit) %>%
    dplyr::mutate(percent_all_reads_explained = total_explained_reads/(n_Hit + n_noHit)) %>%
    dplyr::mutate(p_group.1.2 = paste0(p_group.1, "_", p_group.2)) %>%
    dplyr::mutate(g_group.1.2 = paste0(g_group.1, "_", g_group.2)) %>%
    dplyr::mutate(allele_group.1 = stringr::str_extract(allele.1, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group.2 = stringr::str_extract(allele.2, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group.1.2 = paste0(allele_group.1, "_", allele_group.2))


  ## plotting
  top_pairwise_results_df.plot <-
    top_pairwise_results_df %>%
    dplyr::group_by(p_group.1.2) %>%
    dplyr::top_n(-1, combined.rank) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(row.number.rank = dplyr::row_number(combined.rank)) %>%
    dplyr::top_n(-50, row.number.rank) %>%
    dplyr::arrange(row.number.rank) %>%
    dplyr::mutate(plot.color = as.factor(ifelse(row.number.rank %% 2 != 0, 1, 2)))

  ## factor 8
  single.plot <- ggplot2::ggplot(top_single_res_df, ggplot2::aes(x = reorder_within(allele, n_Hit, allele_group), y = n_Hit)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("allele") +
    ggplot2::ylab("n matching read hits") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::facet_wrap(ggplot2::vars(allele_group), scales = "free_x")

  rank.plot.p1 <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = stats::reorder(p_group.1, row.number.rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p_group.1") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pairwise_results_df), 10)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none")

  rank.plot.p2 <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = stats::reorder(p_group.2, row.number.rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p_group.2") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pairwise_results_df), 10)) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none")

  rank.read.plot <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = total_explained_reads, fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("total explained reads") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pairwise_results_df), 10))

  n_Hit <- as.numeric(levels(as.factor(top_pairwise_results_df.plot$n_Hit)))
  if ((max(top_pairwise_results_df.plot$total_explained_reads) - min(top_pairwise_results_df.plot$total_explained_reads)) > 10) {
    rank.read.plot <- rank.read.plot + ggplot2::scale_y_continuous(breaks = seq(max(top_pairwise_results_df.plot$total_explained_reads), min(top_pairwise_results_df.plot$total_explained_reads), by = -floor_any(max(top_pairwise_results_df.plot$total_explained_reads) - min(top_pairwise_results_df.plot$total_explained_reads), 10)), sec.axis = ggplot2::sec_axis(~ . / n_Hit * 100, name = "% of matching reads that\nmatched min. one ref. allele"))
  } else if ((max(top_pairwise_results_df.plot$total_explained_reads) - min(top_pairwise_results_df.plot$total_explained_reads)) > 5) {
    rank.read.plot <- rank.read.plot + ggplot2::scale_y_continuous(breaks = seq(max(top_pairwise_results_df.plot$total_explained_reads), min(top_pairwise_results_df.plot$total_explained_reads), by = -5), sec.axis = ggplot2::sec_axis(~ . / n_Hit * 100, name = "% of matching reads that\nmatched min. one ref. allele"))
  } else {
    rank.read.plot <- rank.read.plot + ggplot2::scale_y_continuous(breaks = seq(max(top_pairwise_results_df.plot$total_explained_reads), min(top_pairwise_results_df.plot$total_explained_reads), by = -1), sec.axis = ggplot2::sec_axis(~ . / n_Hit * 100, name = "% of matching reads that\nmatched min. one ref. allele"))
  }

  height.1 <- nlevels(as.factor(top_pairwise_results_df.plot$p_group.1))
  height.2 <- nlevels(as.factor(top_pairwise_results_df.plot$p_group.2))
  if (height.1 / height.2 < 0.1) {
    height.1 <- 10
    height.2 <- 100 - height.1
  }
  if (height.2 / height.1 < 0.1) {
    height.2 <- 10
    height.1<- 100 - height.2
  }
  sum.1.2 <- height.1 + height.2
  height.3 <- sum.1.2*0.15
  height.1 = (sum.1.2 - height.3)*height.1/sum.1.2
  height.2 = (sum.1.2 - height.3)*height.2/sum.1.2
  total <- height.1 + height.2 + height.3
  height.1 <- height.1/total
  height.2 <- height.2/total
  height.3 <- height.3/total

  pairwise.plot <- cowplot::plot_grid(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, align = "v", rel_heights = c(height.1,height.2,height.3))


  #results <- plot_hla_typing_results(single_res, top_single_res_df, pairwise_results, pairwise_results_df, top_pairwise_results_df)
## check what to return
  #return(results)
}


integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(base::pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

ceiling_any = function(x, accuracy, f = ceiling){f(x/ accuracy) * accuracy}

floor_any = function(x, accuracy, f = floor){f(x/ accuracy) * accuracy}


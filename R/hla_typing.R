#' Align reads from (sc)RNAseq to known HLA reference alleles, count the number of matches and infer the HLA type
#'
#' This function tries to infer the HLA type by plotting the number of matching reads. Inspection of output-plots is seen as an alternative to statistical analysis.
#' Please provide data (read and hla_ref) for only one gene at a time, e.g. only HLA-A or HLA-B or HLA-C. Run the function sequentially for every single gene.
#' Not every allele-combination is easily inferred. Results may vary with respect to persuasiveness. It is assumed that inferring the p-group is sufficient. Also deeper typing becomes uncertain.
#' The results of pairwise matches are hence plotted by p-groups.
#'
#' Reads are checked for perfect matches (hits) in all provided hla reference alleles (or a sub-sequence of them, e.g. only exon 2 and 3 which are subject to highest variation).
#' For every hla allele the number of hits are counted - the allele can 'explain' a number of reads. This is done for every allele on its own.
#' Then pairwise combinations of alleles are checked for the number of reads that they explain redundantly (double_explained_reads) or uniquely (uniquely_explained_reads) and in total (total_explained_reads).
#' The allele-combinations with the highest ranks (combination of total_explained_reads and uniquely_explained_reads) likely reflect the cells' HLA type.
#'
#'
#' @param hla_ref a data frame preferentially prepared with scexpr::hla_df_from_xml
#' @param reads a data frame preferentially prepared with scexpr::reads_from_bam
#' @param allele_diff maximum allowed difference (factor) of read abundances between two alleles;
#' e.g. if highest abundant HLA-A allele 1 has 20 matched reads and diff_allele_hits_single_results = 5, then other alleles need to have
#' at least 4 (=20/5) matched reads to be considered in the subsequent pairwise matching. In other words, what is the maximum expected/allowed
#' expression difference of alleles HLA-A from father and HLA-A from mother. Intended to speed up the subsequent pairwise matching.
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param top_n_pairwise_results top ranks of unique_explained_reads_rank and total_explained_reads_rank of pairwise matches used for plotting
#' @param hla_seq_colName column name of hla sequences in hla_ref
#' @param read_seq_colName column name read sequences in reads
#' @param hla_allele_colName column name of allele names in hla_ref
#' @param read_name_colName column name of read names in reads
#' @param p_group_colName column name of p_group in hla_ref
#' @param g_group_colName column name of g_group in hla_ref
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#'
#' @importFrom magrittr %>%
#'
#' @return list of data frames and ggplot2 objects which, upon visual inspection, may allow to infer hla type
#' @export
#'
#' @examples
hla_typing <- function(hla_ref,
                       reads,
                       allele_diff = 5,
                       top_n_pairwise_results = 50,
                       hla_seq_colName = "seq_Exon2_3",
                       read_seq_colName = "seq",
                       hla_allele_colName = "allele",
                       read_name_colName = "readName",
                       p_group_colName = "p_group",
                       g_group_colName = "g_group",
                       lapply_fun = lapply,
                       ...) {

  # check for unique read names
  # check for columns in data frames

  if (!requireNamespace("BiocManager", quietly = T)) {
    utils::install.packages("BiocManager")
  }
  if (!requireNamespace("Biostrings", quietly = T)) {
    BiocManager::install("Biostrings")
  }

  lapply_fun <- match.fun(lapply_fun)

  print("Calculating single matches.")
  single_res <- Biostrings::vcountPDict(subject = Biostrings::DNAStringSet(hla_ref[,hla_seq_colName,drop=T]), pdict = Biostrings::PDict(reads[,read_seq_colName,drop=T]))
  colnames(single_res) <- hla_ref[,hla_allele_colName,drop=T]
  rownames(single_res) <- reads[,read_name_colName,drop=T]
  n_noHit <- sum(Matrix::rowSums(single_res) == 0)
  n_Hit <- sum(Matrix::rowSums(single_res) > 0)

  top_single_res <- single_res[Matrix::rowSums(single_res) > 0, Matrix::colSums(single_res) >= max(Matrix::colSums(single_res))/allele_diff]
  top_single_res_df <-
    data.frame(n_Hit = Matrix::colSums(top_single_res)) %>%
    tibble::rownames_to_column(hla_allele_colName) %>%
    dplyr::left_join(hla_ref, by = hla_allele_colName) %>%
    dplyr::mutate(rank = dplyr::row_number(-n_Hit))

  col.combs <- utils::combn(ncol(top_single_res), 2, simplify = F)
  print(paste0("Calculating pairwise matches. Combinations: ", length(col.combs), "."))
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
    dplyr::bind_rows(pairwise_results_df %>% dplyr::top_n(-top_n_pairwise_results*2, unique_explained_reads_rank), pairwise_results_df %>% dplyr::top_n(-top_n_pairwise_results*2, total_explained_reads_rank)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(combined.rank = dplyr::dense_rank(base::interaction(-total_explained_reads, -unique_explained_reads, lex.order = TRUE)))

  top_pairwise_results_df$allele.1.2 <- apply(top_pairwise_results_df[,c("allele.1", "allele.2")], 1, function(x) paste(sort(c(x[1], x[2])), collapse = "_"))

  top_pairwise_results_df <-
    top_pairwise_results_df %>%
    dplyr::group_by(allele.1.2) %>%
    dplyr::filter(combined.rank == min(combined.rank)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(allele.1.2, .keep_all = T) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(!!rlang::sym(hla_allele_colName), !!rlang::sym(p_group_colName), !!rlang::sym(g_group_colName)), by = c("allele.1" = hla_allele_colName)) %>% dplyr::rename("p_group.1" = p_group_colName, "g_group.1" = g_group_colName) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(!!rlang::sym(hla_allele_colName), !!rlang::sym(p_group_colName), !!rlang::sym(g_group_colName)), by = c("allele.2" = hla_allele_colName)) %>% dplyr::rename("p_group.2" = p_group_colName, "g_group.2" = g_group_colName) %>%
    dplyr::left_join(top_single_res_df %>% dplyr::select(!!rlang::sym(hla_allele_colName), n_Hit), by = c("allele.1" = hla_allele_colName)) %>% dplyr::rename("explained_reads_allele.1" = n_Hit) %>%
    dplyr::left_join(top_single_res_df %>% dplyr::select(!!rlang::sym(hla_allele_colName), n_Hit), by = c("allele.2" = hla_allele_colName)) %>% dplyr::rename("explained_reads_allele.2" = n_Hit) %>%
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
    dplyr::top_n(-top_n_pairwise_results, row.number.rank) %>%
    dplyr::arrange(row.number.rank) %>%
    dplyr::mutate(plot.color = as.factor(ifelse(row.number.rank %% 2 != 0, 1, 2)))


  single.plot <- ggplot2::ggplot(top_single_res_df, ggplot2::aes(x = reorder_within(!!rlang::sym(hla_allele_colName), n_Hit, allele_group), y = n_Hit)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("allele") +
    ggplot2::ylab("n explained reads") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::facet_wrap(ggplot2::vars(allele_group), scales = "free_x")

  rank.plot.p1 <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = stats::reorder(p_group.1, row.number.rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 1") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pairwise_results_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.plot.p2 <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = stats::reorder(p_group.2, row.number.rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 2") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pairwise_results_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.read.plot <- ggplot2::ggplot(top_pairwise_results_df.plot, ggplot2::aes(x = as.factor(row.number.rank), y = total_explained_reads, fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("total explained reads") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
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

  pairwise.plot <- cowplot::plot_grid(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, align = "v", rel_heights = c(height.1,height.2,height.3)) # check how to replace with patchwork

  return(list(top_single_res_df = top_single_res_df,
              top_single_res_matrix = top_single_res,
              top_pairwise_res_df = top_pairwise_results_df,
              top_pairwise_res_for_plotting_df = top_pairwise_results_df.plot,
              top_single_res_plot = single.plot,
              top_pairwise_res_plot = pairwise.plot,
              rank_plot_p1 = rank.plot.p1,
              rank_plot_p2 = rank.plot.p2,
              rank_read_plot = rank.read.plot))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

ceiling_any = function(x, accuracy, f = ceiling) {
  f(x / accuracy) * accuracy
}

floor_any = function(x, accuracy, f = floor) {
  f(x / accuracy) * accuracy
}


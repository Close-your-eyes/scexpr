#' Align reads from (sc)RNAseq to known HLA reference alleles, count the number of matches and infer the HLA type
#'
#' This function tries to infer the HLA type by plotting the number of matching reads. Inspection of output-plots is seen as an alternative to statistical analysis.
#' Please provide data (read and hla_ref) for only one gene at a time, e.g. only HLA-A or HLA-B or HLA-C. Run the function sequentially for every single gene.
#' Not every allele-combination is easily inferred. Results may vary with respect to persuasiveness. It is assumed that inferring the p-group is sufficient. Also deeper typing becomes uncertain.
#' The results of pairwise matches are hence plotted by p-groups.
#'
#' Reads are checked for perfect matches (hits) in all provided hla reference alleles (or a sub-sequence of them, e.g. only exon 2 and 3 which are subject to highest variation).
#' For every hla allele the number of hits are counted - the allele can 'explain' a number of reads. This is done for every allele on its own.
#' Then pairwise combinations of alleles are checked for the number of reads that they explain redundantly (dbl_expl_reads) or uniquely (uniquely_explained_reads) and in total (tot_expl_reads).
#' The allele-combinations with the highest ranks (combination of tot_expl_reads and uniquely_explained_reads) likely reflect the cells' HLA type.
#'
#'
#' @param hla_ref a data frame preferentially prepared with scexpr::hla_df_from_xml
#' @param reads a data frame preferentially prepared with scexpr::reads_from_bam
#' @param allele_diff maximum allowed difference (factor) of read abundances between two alleles;
#' e.g. if highest abundant HLA-A allele 1 has 20 matched reads and diff_allele_hits_single_results = 5, then other alleles need to have
#' at least 4 (=20/5) matched reads to be considered in the subsequent pairwise matching. In other words, what is the maximum expected/allowed
#' expression difference of alleles HLA-A from father and HLA-A from mother. Intended to speed up the subsequent pairwise matching.
#' @param lapply_fun function name without quotes; lapply, pbapply::pblapply or parallel::mclapply are suggested
#' @param top_n_pairwise_results top ranks of uni_expl_reads_rank and tot_expl_reads_rank of pairwise matches used for plotting
#' @param hla_seq_colName column name of hla sequences in hla_ref
#' @param read_seq_colName column name read sequences in reads
#' @param hla_allele_colName column name of allele names in hla_ref
#' @param read_name_colName column name of read names in reads
#' @param p_group_colName column name of p_group in hla_ref
#' @param g_group_colName column name of g_group in hla_ref
#' @param ... additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen
#' @param maxmis
#'
#' @importFrom magrittr %>%
#' @import Matrix
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
                       hla_seq_colName2 = "seq",
                       read_seq_colName = "seq",
                       hla_allele_colName = "allele",
                       read_name_colName = "readName",
                       p_group_colName = "p_group",
                       g_group_colName = "g_group",
                       lapply_fun = lapply,
                       maxmis = 0,
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
  arg_list <- list(...)

  if (missing(hla_ref)) {
    stop("hla_ref missing.")
  }
  if (missing(reads)) {
    stop("reads missing.")
  }
  if (allele_diff < 1 || allele_diff == 1) {
    stop("allele_diff has to be larger than 1.")
  }
  if (top_n_pairwise_results < 1) {
    stop("top_n_pairwise_results has to be at least 1.")
  }

  if (!read_seq_colName %in% names(reads)) {
    stop(read_seq_colName, " not found in names of reads.")
  }
  if (!read_name_colName %in% names(reads)) {
    stop(read_name_colName, " not found in names of reads.")
  }
  if (!hla_seq_colName %in% names(hla_ref)) {
    stop(hla_seq_colName, " not found in names of hla_ref")
  }
  if (!p_group_colName %in% names(hla_ref)) {
    stop(p_group_colName, " not found in names of hla_ref")
  }
  if (!g_group_colName %in% names(hla_ref)) {
    stop(g_group_colName, " not found in names of hla_ref")
  }
  if (!hla_allele_colName %in% names(hla_ref)) {
    stop(hla_allele_colName, " not found in names of hla_ref")
  }

  if (anyDuplicated(reads[,read_name_colName,drop=T])) {
    message("Duplicated ", read_name_colName, " found. Will make them unique.")
    reads[,read_name_colName] <- make.unique(reads[,read_name_colName,drop=T])
  }

  if (length(invalid <- which(grepl("[^ACTGU]", reads[,read_seq_colName,drop=T]))) > 0) {
    message(length(invalid), " sequences with non-DNA or non-RNA characters detected. Those are excluded.")
    reads <- reads[-invalid,]
    if (nrow(reads) == 0) {
      stop("No reads left after filtering for ones with valid DNA/RNA characters only. Please fix the sequences.")
    }
  }
  # remove columns with all NA, e.g. for HLA-B there is no exon8 --> remove the column
  hla_ref <- hla_ref[,colSums(is.na(hla_ref)) < nrow(hla_ref)]

  complete_alleles <-
    hla_ref %>%
    dplyr::select(dplyr::matches("Exon[[:digit:]]$"), dplyr::contains("UTR"), allele) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::pull(!!rlang::sym(hla_allele_colName))

  hla_ref <-
    hla_ref %>%
    dplyr::mutate(complete_seq = allele %in% complete_alleles)

  hla_ref2 <-
    hla_ref %>%
    dplyr::group_by(allele_coding) %>%
    dplyr::mutate(allele_coding_Exon2_3_max_length = seq_Exon2_3_length == max(seq_Exon2_3_length))

  hla_ref3 <-
    hla_ref2 %>%
    dplyr::filter(allele_coding_Exon2_3_max_length) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(g_group, seq_Exon2_3_length, .keep_all = TRUE)


  first_round_results <- run_read_matching_and_report_results(hla_ref = hla_ref,
                                                              reads = reads,
                                                              allele_diff = allele_diff,
                                                              top_n_pairwise_results = top_n_pairwise_results,
                                                              hla_seq_colName = hla_seq_colName,
                                                              read_seq_colName = read_seq_colName,
                                                              hla_allele_colName = hla_allele_colName,
                                                              read_name_colName = read_name_colName,
                                                              p_group_colName = p_group_colName,
                                                              g_group_colName = g_group_colName,
                                                              lapply_fun = lapply_fun,
                                                              maxmis = maxmis,
                                                              arg_list = arg_list,
                                                              ...)

  allele_group12_medians <-
    first_round_results[["pair_res1_df"]] %>%
    dplyr::group_by(allele_group12) %>%
    dplyr::summarise(median_frac_match_reads_expl = median(frac_match_reads_expl)) %>%
    dplyr::arrange(dplyr::desc(median_frac_match_reads_expl))

  top_alleles <-
    first_round_results[["pair_res1_df"]] %>%
    dplyr::filter(frac_match_reads_expl >= allele_group12_medians[1,2,drop=T]) %>%
    dplyr::select(allele1, allele2)

  hla_ref_top <-
    hla_ref %>%
    dplyr::filter(!!rlang::sym(hla_allele_colName) %in% unique(c(top_alleles$allele1, top_alleles$allele2)))

  hla_ref_top_complete <-
    hla_ref %>%
    dplyr::filter(!!rlang::sym(p_group_colName) %in% unique(hla_ref_top[[p_group_colName]])) %>%
    dplyr::select(dplyr::matches("Exon[[:digit:]]$"), dplyr::contains("UTR"), dplyr::all_of(c(hla_allele_colName, p_group_colName))) %>%
    dplyr::filter(stats::complete.cases(.))

  second_round_results <- NULL
  if (!is.null(hla_seq_colName2)) {
    second_round_results <- run_read_matching_and_report_results(hla_ref = hla_ref_top_complete,
                                                                 reads = reads,
                                                                 allele_diff = allele_diff,
                                                                 top_n_pairwise_results = top_n_pairwise_results,
                                                                 hla_seq_colName = hla_seq_colName2,
                                                                 read_seq_colName = read_seq_colName,
                                                                 hla_allele_colName = hla_allele_colName,
                                                                 read_name_colName = read_name_colName,
                                                                 p_group_colName = p_group_colName,
                                                                 g_group_colName = g_group_colName,
                                                                 lapply_fun = lapply_fun,
                                                                 maxmis = maxmis,
                                                                 arg_list = arg_list,
                                                                 ...)
  }


  return(list(hla_seq_colName = first_round_results, hla_seq_colName2 = second_round_results))

}


run_read_matching_and_report_results <- function(hla_ref,
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
                                                 maxmis = 0,
                                                 arg_list = arg_list,
                                                 ...) {

  library(Matrix) # required for sparseMatrix below; saves memory
  print("Calculating single matches.")
  if ("strand" %in% names(reads)) {
    single_res <- lapply(split(reads, as.character(reads$strand)),
                         single_matching,
                         lapply_fun = lapply_fun,
                         hla_ref = hla_ref,
                         maxmis = maxmis,
                         hla_seq_colName = hla_seq_colName,
                         hla_allele_colName = hla_allele_colName,
                         read_name_colName = read_name_colName,
                         read_seq_colName = read_seq_colName,
                         ...)
    for (i in names(single_res)) {
      message("Strand (", i, "):")
      reads_w_no_match_sum <- sum(Matrix::rowSums(single_res[[i]]) == 0)
      reads_w_min_one_match_sum <- sum(Matrix::rowSums(single_res[[i]]) > 0)
      message(reads_w_min_one_match_sum, " reads with at least one match/hit, (", round(reads_w_min_one_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
      message(reads_w_no_match_sum, " reads with no match/hit, (", round(reads_w_no_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
    }
    single_res <- do.call(rbind, single_res)
  } else {
    single_res <- lapply(reads,
                         single_matching,
                         lapply_fun = lapply_fun,
                         hla_ref = hla_ref,
                         maxmis = maxmis,
                         hla_seq_colName = hla_seq_colName,
                         hla_allele_colName = hla_allele_colName,
                         ...)
  }

  reads_w_min_one_match <- Matrix::rowSums(single_res) > 0
  expl_reads_per_allele <- Matrix::colSums(single_res)

  reads_w_no_match_sum <- sum(Matrix::rowSums(single_res) == 0)
  reads_w_min_one_match_sum <- sum(reads_w_min_one_match)

  message("Total:")
  if (reads_w_min_one_match_sum == 0) {
    stop("No matches/hits determined.")
  }
  message(reads_w_min_one_match_sum, " reads with at least one match/hit, (", round(reads_w_min_one_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")
  message(reads_w_no_match_sum, " reads with no match/hit, (", round(reads_w_no_match_sum/(reads_w_min_one_match_sum+reads_w_no_match_sum)*100), " %)")

  # as.matrix here as this will speed up iteration over col.combs below!
  top_single_res <- as.matrix(single_res[reads_w_min_one_match, expl_reads_per_allele >= max(expl_reads_per_allele)/allele_diff])
  top_sin_res_df <-
    data.frame(expl_reads = Matrix::colSums(top_single_res)) %>%
    tibble::rownames_to_column(hla_allele_colName) %>%
    dplyr::left_join(hla_ref, by = hla_allele_colName) %>%
    dplyr::mutate(rank = dplyr::row_number(-expl_reads))


  # could also be made mapply or purrr::map2 with utils::combn(ncol(top_single_res), 2, simplify = T)
  # but would require another argument to define with mapply fun to use
  #col.combs <- utils::combn(1:ncol(top_single_res), 2, simplify = F)
  #print(paste0("Calculating pairwise matches. Combinations: ", length(col.combs), "."))

  col.combs <- t(utils::combn(1:ncol(top_single_res), 2, simplify = T))
  print(paste0("Calculating pairwise matches. Combinations: ", nrow(col.combs), "."))

  # Call the C++ function to calculate row sums
  # Source the C++ code
  '
  sourceCpp("/Users/vonskopnik/Documents/R_packages/scexpr/src/calculateRowSumsInCpp2.cpp")
  system.time(results <- countOccurrencesInCpp(top_single_res, col.combs)) # col.combs[1:2,,drop=F]
  col.combs2 <- split_mat(col.combs, n_chunks = 4, byrow = T)
  system.time(
    pairwise_results <- lapply_fun(col.combs2, function(x) {
      countOccurrencesInCpp(top_single_res, x)
    }, ...)
  )'

  # doing this in R was too slow. other packages did not have the functionality
  # so, written in c++
  if (identical(lapply_fun, parallel::mclapply) && "mc.cores" %in% names(arg_list)) {
    # split_mat fun from scexpr package
    # Split the matrix into chunks for multithreading
    pairwise_results <- lapply_fun(scexpr::split_mat(col.combs, n_chunks = arg_list[["mc.cores"]], byrow = T), function(x) {
      scexpr:::countOccurrencesInCpp(top_single_res, x)
    }, ...)
    pairwise_results <- do.call(rbind, pairwise_results)
  } else {
    pairwise_results <- scexpr:::countOccurrencesInCpp(top_single_res, col.combs)
  }

  # other attempts with collapse package
  # however, it does not support subsetting in c++
  '  #top_single_res[1:3,1:5]
  pairwise_results <- matrix(0, nrow = length(col.combs), ncol = 2)
  # Iterate over each combination of columns
  for (i in seq_along(col.combs)) {
    # Compute row sums for the current combination using fsum()
    temp <- collapse::fsum(t(top_single_res[, col.combs[[i]]]))
    # Count occurrences of 1 and 2 in row sums
    pairwise_results[i, 1] <- sum(temp == 1)
    pairwise_results[i, 2] <- sum(temp == 2)
  }
'

  '
  tt <- collapse::fsubset(top_single_res, subset = list(c(1:2), c(4:5)))
  out <- collapse::fsum(x = t(top_single_res), vars = c(1,2))
  collapse::collap(t(top_single_res), by = c(1,2), FUN = collapse::fsum)
  m <- qM(mtcars)
  nrow(m)
  tt <- fsum(m, g = list(c(rep(1,16),rep(2,16)),
                         c(rep(3,10),rep(4,22)),
                         c(rep(5,5),rep(6,27))))

  fsum(m, g)
  pairwise_results <- purrr::map(col.combs, function(x) {
    temp <- Matrix::rowSums(top_single_res[,x])
    return(c(sum(temp == 1), sum(temp == 2)))
  }, .progress = T)'


  ## old version, working:
  'col.combs <- utils::combn(1:ncol(top_single_res), 2, simplify = F)
  pairwise_results <- lapply_fun(col.combs, function(x) {
    temp <- matrixStats::rowSums2(top_single_res, cols = x)
    #temp <- collapse::fsum(t(top_single_res[,x]))
    #temp <- Matrix::rowSums(top_single_res[,x])
    return(c(sum(temp == 1), sum(temp == 2)))
  }, ...)'

  pair_res_df <-
    data.frame(uni_expl_reads = pairwise_results[,1], # sapply(pairwise_results, "[", 1),
               dbl_expl_reads = pairwise_results[,2], # sapply(pairwise_results, "[", 2),
               allele1 = colnames(top_single_res)[col.combs[,1]],
               allele2 = colnames(top_single_res)[col.combs[,2]]) %>%
    dplyr::mutate(tot_expl_reads = uni_expl_reads + dbl_expl_reads, .after = dbl_expl_reads)  %>%
    dplyr::mutate(uni_expl_reads_rank = dplyr::dense_rank(-uni_expl_reads),
                  tot_expl_reads_rank = dplyr::dense_rank(-tot_expl_reads)) %>%
    dplyr::mutate(rank_sum = (tot_expl_reads_rank + uni_expl_reads_rank)/2)  %>%
    dplyr::mutate(allele_comb = scexpr:::orderAndConcatenateStrings(as.matrix(.[c("allele1", "allele2")])), .after = "allele2") %>%
    #dplyr::mutate(rank_int = dplyr::dense_rank(base::interaction(-tot_expl_reads, -uni_expl_reads, lex.order = TRUE))) %>%
    dplyr::group_by(allele_comb) %>%
    #dplyr::filter(rank_int == min(rank_int)) %>%
    dplyr::filter(rank_sum == min(rank_sum)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(allele_comb, .keep_all = T) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele1" = hla_allele_colName)) %>%
    dplyr::rename("p_group1" = p_group_colName, "g_group1" = g_group_colName) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele2" = hla_allele_colName)) %>%
    dplyr::rename("p_group2" = p_group_colName, "g_group2" = g_group_colName) %>%
    dplyr::left_join(top_sin_res_df %>% dplyr::select(dplyr::all_of(hla_allele_colName), expl_reads), by = c("allele1" = hla_allele_colName)) %>%
    dplyr::rename("allele1_expl_reads" = expl_reads) %>%
    dplyr::left_join(top_sin_res_df %>% dplyr::select(dplyr::all_of(hla_allele_colName), expl_reads), by = c("allele2" = hla_allele_colName)) %>%
    dplyr::rename("allele2_expl_reads" = expl_reads) %>%
    dplyr::mutate(expl_reads_overall = !!reads_w_min_one_match_sum) %>%
    dplyr::mutate(non_expl_reads_overall = !!reads_w_no_match_sum) %>%
    dplyr::mutate(allele12_expl_read_diff = abs(allele1_expl_reads - allele2_expl_reads)) %>%
    dplyr::mutate(frac_match_reads_expl = tot_expl_reads/expl_reads_overall) %>%
    dplyr::mutate(frac_all_reads_expl = tot_expl_reads/(expl_reads_overall + non_expl_reads_overall)) %>%
    dplyr::mutate(p_group12 = paste0(p_group1, "_", p_group2)) %>%
    dplyr::mutate(g_group12 = paste0(g_group1, "_", g_group2)) %>%
    dplyr::mutate(allele_group1 = stringr::str_extract(allele1, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group2 = stringr::str_extract(allele2, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group12 = paste0(allele_group1, "_", allele_group2))


  #pair_res_df <- dplyr::mutate(pair_res_df, allele_comb = scexpr:::orderAndConcatenateStrings(as.matrix(pair_res_df[,c("allele1", "allele2")])), .after = "allele2")
  #top_pair_res_df$allele_comb <- apply(top_pair_res_df[,c("allele1", "allele2")], 1, function(x) paste(sort(c(x[1], x[2])), collapse = "_"))
  #anyDuplicated(pair_res_df$allele_comb)

  '  top_pair_res_df <-
    dplyr::bind_rows(pair_res_df %>%
                       dplyr::slice_min(n = top_n_pairwise_results*2, order_by = uni_expl_reads_rank, with_ties = F),
                     pair_res_df %>%
                       dplyr::slice_min(n = top_n_pairwise_results*2, order_by = tot_expl_reads_rank, with_ties = F)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(rank_int = dplyr::dense_rank(base::interaction(-tot_expl_reads, -uni_expl_reads, lex.order = TRUE))) %>%
    dplyr::group_by(allele_comb) %>%
    dplyr::filter(rank_int == min(rank_int)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(allele_comb, .keep_all = T) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele1" = hla_allele_colName)) %>% dplyr::rename("p_group1" = p_group_colName, "g_group1" = g_group_colName) %>%
    dplyr::left_join(hla_ref %>% dplyr::select(dplyr::all_of(c(hla_allele_colName, p_group_colName, g_group_colName))), by = c("allele2" = hla_allele_colName)) %>% dplyr::rename("p_group2" = p_group_colName, "g_group2" = g_group_colName) %>%
    dplyr::left_join(top_sin_res_df %>% dplyr::select(dplyr::all_of(hla_allele_colName), n_Hit), by = c("allele1" = hla_allele_colName)) %>% dplyr::rename("allele1_expl_reads" = n_Hit) %>%
    dplyr::left_join(top_sin_res_df %>% dplyr::select(dplyr::all_of(hla_allele_colName), n_Hit), by = c("allele2" = hla_allele_colName)) %>% dplyr::rename("allele2_expl_reads" = n_Hit) %>%
    dplyr::mutate(n_Hit = n_Hit) %>%
    dplyr::mutate(n_noHit = n_noHit) %>%
    dplyr::mutate(allele12_expl_read_diff = abs(allele1_expl_reads - allele2_expl_reads)) %>%
    dplyr::mutate(frac_match_reads_expl = tot_expl_reads/n_Hit) %>%
    dplyr::mutate(frac_all_reads_expl = tot_expl_reads/(n_Hit + n_noHit)) %>%
    dplyr::mutate(p_group12 = paste0(p_group1, "_", p_group2)) %>%
    dplyr::mutate(g_group12 = paste0(g_group1, "_", g_group2)) %>%
    dplyr::mutate(allele_group1 = stringr::str_extract(allele1, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group2 = stringr::str_extract(allele2, "[:alpha:]\\*[:digit:]{2}")) %>%
    dplyr::mutate(allele_group12 = paste0(allele_group1, "_", allele_group2))

'
  ## plotting
  top_pair_res_plot1 <-
    #top_pair_res_df
    pair_res_df %>%
    dplyr::group_by(p_group12) %>%
    dplyr::slice_min(order_by = rank_sum, n = 1) %>%
    dplyr::ungroup()

  allele_group12_medians <-
    top_pair_res_plot1 %>%
    dplyr::group_by(allele_group12) %>%
    dplyr::summarise(median_frac_match_reads_expl = median(frac_match_reads_expl)) %>%
    dplyr::arrange(dplyr::desc(median_frac_match_reads_expl))

  overview_plot <- ggplot2::ggplot(top_pair_res_plot1, ggplot2::aes(x = reorder_within(allele_group2, frac_match_reads_expl, allele_group1), y = frac_match_reads_expl)) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    scale_x_reordered() +
    ggplot2::labs(x = "allele_group2") +
    ggplot2::geom_hline(yintercept = allele_group12_medians[1,2,drop=T], color = "tomato2") +
    ggplot2::geom_hline(yintercept = max(top_pair_res_plot1$frac_match_reads_expl), color = "forestgreen") +
    ggplot2::facet_wrap(ggplot2::vars(allele_group1), nrow = 1, scales = "free_x")

  top_pair_res_plot2 <-
    top_pair_res_plot1 %>%
    dplyr::mutate(row_num_rank = dplyr::row_number(rank_sum)) %>%
    dplyr::slice_min(order_by = row_num_rank, n = top_n_pairwise_results) %>%
    dplyr::arrange(row_num_rank) %>%
    dplyr::mutate(plot.color = as.factor(ifelse(row_num_rank %% 2 != 0, 1, 2)))


  sin_plot <- ggplot2::ggplot(top_sin_res_df, ggplot2::aes(x = reorder_within(!!rlang::sym(hla_allele_colName), expl_reads, allele_group), y = expl_reads)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::xlab("allele") +
    ggplot2::ylab("n explained reads") +
    ggplot2::theme_bw() +
    scale_x_reordered() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::facet_wrap(ggplot2::vars(allele_group), scales = "free_x")

  rank.plot.p1 <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(row_num_rank), y = stats::reorder(p_group1, row_num_rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 1") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.plot.p2 <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(row_num_rank), y = stats::reorder(p_group2, row_num_rank), fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::ylab("p group 2") +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier"))

  rank.read.plot <- ggplot2::ggplot(top_pair_res_plot2, ggplot2::aes(x = as.factor(row_num_rank), y = tot_expl_reads, fill = plot.color)) +
    ggplot2::geom_point(size = 2, shape = 21) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("total\nexplained\nreads") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), text = ggplot2::element_text(family = "Courier")) +
    ggplot2::scale_x_discrete(breaks = seq(0, nrow(pair_res_df), 10))

  expl_reads_overall <- as.numeric(levels(as.factor(top_pair_res_plot2$expl_reads_overall)))
  if ((max(top_pair_res_plot2$tot_expl_reads) - min(top_pair_res_plot2$tot_expl_reads)) > 10) {
    rank.read.plot <-
      rank.read.plot +
      ggplot2::scale_y_continuous(breaks = seq(max(top_pair_res_plot2$tot_expl_reads), min(top_pair_res_plot2$tot_expl_reads),
                                               by = -floor_any((max(top_pair_res_plot2$tot_expl_reads) - min(top_pair_res_plot2$tot_expl_reads))/3, 10)),
                                  sec.axis = ggplot2::sec_axis(~ . / expl_reads_overall * 100, name = "% of matching\nreads that\nmatched min. one\nref. allele"))
  } else if ((max(top_pair_res_plot2$tot_expl_reads) - min(top_pair_res_plot2$tot_expl_reads)) > 5) {
    rank.read.plot <-
      rank.read.plot +
      ggplot2::scale_y_continuous(breaks = seq(max(top_pair_res_plot2$tot_expl_reads), min(top_pair_res_plot2$tot_expl_reads), by = -5),
                                  sec.axis = ggplot2::sec_axis(~ . / expl_reads_overall * 100, name = "% of matching\nreads that\nmatched min. one\nref. allele"))
  } else {
    rank.read.plot <-
      rank.read.plot +
      ggplot2::scale_y_continuous(breaks = seq(max(top_pair_res_plot2$tot_expl_reads), min(top_pair_res_plot2$tot_expl_reads), by = -1),
                                  sec.axis = ggplot2::sec_axis(~ . / expl_reads_overall * 100, name = "% of matching\nreads that\nmatched min. one\nref. allele"))
  }

  height.1 <- nlevels(as.factor(top_pair_res_plot2$p_group1))
  height.2 <- nlevels(as.factor(top_pair_res_plot2$p_group2))
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

  #pair_plot <- cowplot::plot_grid(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, align = "v", rel_heights = c(height.1,height.2,height.3)) # check how to replace with patchwork
  pair_plot <- patchwork::wrap_plots(rank.plot.p1, rank.plot.p2, rank.read.plot, ncol = 1, heights = c(height.1,height.2,height.3))

  return(list(top_sin_res_df = top_sin_res_df,
              top_sin_res_mat = top_single_res,
              pair_res_df = pair_res_df,
              pair_res1_df = top_pair_res_plot1,
              pair_res2_df = top_pair_res_plot2,
              plot_sin_res = sin_plot,
              plot_pair_res1 = overview_plot,
              plot_pair_res2 = pair_plot,
              plot_rank_1 = rank.plot.p1,
              plot_rank_2 = rank.plot.p2,
              plot_rank_3 = rank.read.plot))
}



single_matching <- function(reads,
                            lapply_fun,
                            hla_ref,
                            maxmis,
                            hla_seq_colName,
                            hla_allele_colName,
                            read_name_colName,
                            read_seq_colName,
                            ...) {

  reads <- setNames(reads[,read_seq_colName,drop=T], reads[,read_name_colName,drop=T])
  single_res <- do.call(rbind, lapply_fun(split(reads,ceiling(seq_along(reads)/1e3)), function(rrr) {

    'hit_matrix <- methods::as(Biostrings::vcountPDict(subject = Biostrings::DNAStringSet(hla_ref[,hla_seq_colName,drop=T]),
                                                      pdict = Biostrings::PDict(reads[rows,read_seq_colName,drop=T], max.mismatch = 0),
                                                      max.mismatch = 0),
                              "sparseMatrix")'
    # vwhichPDict allows for max.mismatch, but vcountPDict does not
    hit_matrix <- Biostrings::vwhichPDict(subject = Biostrings::DNAStringSet(hla_ref[,hla_seq_colName,drop=T]),
                                          pdict = Biostrings::PDict(rrr, max.mismatch = maxmis),
                                          max.mismatch = maxmis)
    hit_matrix <- lapply(hit_matrix, function(z) replace(rep(0, length(rrr)), z, 1))
    hit_matrix <- methods::as(do.call(cbind, hit_matrix), "sparseMatrix")
    colnames(hit_matrix) <- hla_ref[,hla_allele_colName,drop=T]
    rownames(hit_matrix) <- names(rrr)
    return(hit_matrix)
  }, ...))
  return(single_res)
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

ceiling_any = function(x, accuracy, f = ceiling) {
  f(x / accuracy) * accuracy
}

floor_any = function(x, accuracy, f = floor) {
  f(x / accuracy) * accuracy
}






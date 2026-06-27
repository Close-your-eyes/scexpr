#' Summarize DEG from multiple comparisons for volcano plot
#'
#' Makes a data frame to use with volcano02_plot.
#'
#' @param volc01_df_list list of data frames from scexpr::volcano01_calc
#' @param neg_name Name of the negative group. If `NULL`, the name is
#'   inferred from the first data frame in `volc01_df_list`.
#' @param pos_name see neg_name
#' @param logfc_avg_fun function for averaging log fold changes
#'
#' @returns data frame
#' @export
#'
#' @examples
volcano03_avg <- function(volc01_df_list,
                          neg_name,
                          pos_name,
                          logfc_avg_fun = mean) {

  id <- "volc_ind"
  logfc <- purrr::map_dfr(volc01_df_list, ~dplyr::select(.x, feature, avg_log2FC), .id = id) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(id), values_from = avg_log2FC) |>
    tidyr::drop_na()
  pval <- purrr::map_dfr(volc01_df_list, ~dplyr::select(.x, feature, p_val_adj), .id = id) |>
    tidyr::pivot_wider(names_from = !!rlang::sym(id), values_from = p_val_adj) |>
    tidyr::drop_na() |>
    dplyr::mutate(dplyr::across(-1, ~scales::squish(.x, range = c(1e-300, 1))))
  avglogfc <- apply(logfc[,-1], 1, logfc_avg_fun)
  signcon <- apply(logfc[,-1], 1, function(x) mean(sign(x)))
  avgpval <- metapod::combineParallelPValues(as.list(pval[,-1]), method = "simes")
  avgdf <- data.frame(feature = logfc$feature,
                      avg_log2FC = avglogfc,
                      p_val_adj = avgpval$p.value,
                      signcon = signcon)

  #dplyr::filter(abs(signcon) == 1)
  try(expr = {
    negpos_name <- gsub("^pct\\.", "", grep("pct", names(volc01_df_list[[1]]), value = T)[1:2])
    if (missing(neg_name)) {
      neg_name <- negpos_name[2]
    }
    if (missing(pos_name)) {
      pos_name <- negpos_name[1]
    }
    attr(avgdf, "neg_name") <- neg_name
    attr(avgdf, "pos_name") <- pos_name
  }, silent = T)

  return(avgdf)
}

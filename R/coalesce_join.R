#' Coalescing join function which will update NA values in the left-hand data frame
#'
#' This is a combination of join functions and the coalesce function from dplyr.
#' It is a convenient way to solve the generic task of updating a data frame (replacing NAs)
#' with another one that holds additional information. In its current form the information (non-NA)
#' in the left-hand data frame (x) will be prioritized over that information in the right-hand
#' data frame (y).
#'
#' originally from: https://alistaire.rbind.io/blog/coalescing-joins/
#'
#' @param x left-hand data frame
#' @param y right-hand data frame
#' @param by which column(s) to join by
#' @param suffix
#' @param join
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
coalesce_join <- function(x, y,
                          by = NULL, suffix = c(".x", ".y"),
                          join = dplyr::left_join, ...) {

  # copied originally from: https://alistaire.rbind.io/blog/coalescing-joins/
  # ideas:
  # https://stackoverflow.com/questions/33954292/merge-two-data-frame-and-replace-the-na-value-in-r#33954334
  # https://community.rstudio.com/t/merging-2-dataframes-and-replacing-na-values/32123/2
  # https://github.com/WinVector/rqdatatable

  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))

  to_coalesce <- names(joined)[!names(joined) %in% cols]
  if (length(to_coalesce) > 0) {
    suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
    # remove suffixes and deduplicate
    to_coalesce <- unique(substr(
      to_coalesce,
      1,
      nchar(to_coalesce) - nchar(suffix_used)
    ))

    coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
      joined[[paste0(.x, suffix[1])]],
      joined[[paste0(.x, suffix[2])]]
    ))
    names(coalesced) <- to_coalesce

    dplyr::bind_cols(joined, coalesced)[cols]
  } else {
    joined
  }

}

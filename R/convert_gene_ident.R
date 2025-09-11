#' Convert gene identifiers
#'
#' Uses biomaRt. See biomaRt::listFilters() and biomaRt::listAttributes()
#'
#' @param x input identifiers
#' @param input input ident, filters in biomaRt
#' @param output output ident, attributes in biomaRt
#'
#' @returns list()
#' @export
#'
#' @examples
#' convert_gene_ident(x = c("CD8A", "CD4"))
convert_gene_ident <- function(x,
                               input = c("hgnc_symbol",
                                         "external_synonym",
                                         "external_gene_name"),
                               output = c("ensembl_gene_id",
                                          "hgnc_symbol",
                                          "external_gene_name",
                                          "external_synonym",
                                          "gene_biotype",
                                          "description")) {
  # make sure all input are in output
  output <- unique(c(output, input))

  mart <- biomaRt::useEnsembl("genes", dataset="hsapiens_gene_ensembl")
  res <- data.frame()
  i <- 1

  while(length(x) > 0 && i <= length(input)) {
    message(input[i])
    res <- dplyr::bind_rows(res,
                            biomaRt::getBM(
                              attributes = output,
                              filters = input[i],
                              values = x,
                              mart = mart) |>
                              dplyr::mutate(dplyr::across(dplyr::everything(), as.character))) |>
      dplyr::mutate(input = !!rlang::sym(input[i]))
    x <- setdiff(x, res[[input[i]]])
    i <- i + 1
  }
  res$GRCh <- "38"

  if (length(x) > 0) {
    message("GRCh37")
    mart <- biomaRt::useEnsembl("genes", dataset="hsapiens_gene_ensembl", GRCh=37)
    i <- 1

    while(length(x) > 0 && i <= length(input)) {
      message(input[i])
      res <- dplyr::bind_rows(res,
                              biomaRt::getBM(
                                attributes = output,
                                filters    = input[i],
                                values     = x,
                                mart       = mart) |>
                                dplyr::mutate(dplyr::across(dplyr::everything(), as.character))) |>
        dplyr::mutate(input = !!rlang::sym(input[i]))
      x <- setdiff(x, res[[input[i]]])
      i <- i + 1
    }
    res$GRCh[which(is.na(res$GRCh))] <- "37"
  }

  res$description <- sapply(strsplit(res$description, " \\[Source"), "[", 1)

  res2 <- res |>
    dplyr::group_by(input, hgnc_symbol, description, GRCh) |>
    dplyr::summarise(ensembl_gene_id = paste(unique(ensembl_gene_id), collapse = ","),
                     external_gene_name = paste(unique(external_gene_name), collapse = ","),
                     external_synonym = paste(unique(external_synonym), collapse = ","),
                     gene_biotype = paste(unique(gene_biotype), collapse = ","),
                     .groups = "drop")

  return(list(long = res, short = res2, unmapped = x))
}

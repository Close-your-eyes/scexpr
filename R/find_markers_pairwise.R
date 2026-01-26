#' Find marker genes in pairwise comparisons
#'
#' @param obj Seurat object
#' @param group column in meta data
#' @param mc.cores multiple threads?
#' @param ... arguments to Seurat::FindMarkers
#' @param split column in obj to split by before marker calculation
#' @param redundant_return with utils::combn unique combinations of ident.1 and
#' ident.2 are calculation and returned as list of lists; by default this is
#' made redundant by adding the flip order of ident.1 and ident.2 to create
#' returned list of lists; with
#'
#' @returns
#' @export
#'
#' @examples
#'\dontrun{
#' markers_pairwise <- scexpr::find_markers_pairwise(so, group = "label_main", mc.cores = 6)
#'
#' # turn this into a result like from SingleR::getClassicMarkers but based on more statistic
#' markers_pairwise_seurat <- markers_pairwise |>
#'   dplyr::slice_min(order_by = p_val, n = 100, by = c(ident.1, ident.2), with_ties = F) |>
#'   dplyr::arrange(p_val)
#' markers_pairwise_seurat <- split(markers_pairwise_seurat, markers_pairwise_seurat$ident.1)
#' markers_pairwise_seurat <- purrr::map(markers_pairwise_seurat, ~split(.x$gene, .x$ident.2))
#' markers_pairwise_seurat_full <- scexpr:::make_list_of_lists_redundant(markers_pairwise_seurat)
#'}
find_markers_pairwise <- function(obj,
                                  group,
                                  #assay = "RNA",
                                  #layer = "data",
                                  split = NULL,
                                  mc.cores = 1,
                                  redundant_return = T,
                                  ...) {

  lvls <- as.character(unique(obj@meta.data[[group]]))
  all_pairs <- utils::combn(lvls, 2, simplify = F)
  Seurat::Idents(obj) <- obj@meta.data[[group]]
  if (!is.null(split)) {
    obj <- Seurat::SplitObject(obj, split.by = split)
  } else {
    obj <- list(seurat = obj)
  }

  out <- purrr::map_dfr(obj, function(obj) {
    out <- parallel::mclapply(all_pairs, function(z) {
      tryCatch(expr = {
        out <- Seurat::FindMarkers(obj, ident.1 = z[1], ident.2 = z[2], ...) |>
          dplyr::mutate(ident.1 = z[1], ident.2 = z[2]) |>
          tibble::rownames_to_column("gene")
      }, error = function(err) {
        out <- NULL
      })
      return(out)
    }, mc.cores = mc.cores)
    out <- dplyr::bind_rows(out)
    return(out)
  }, .id = "split")

  return(out)

  # fun <- if (assay %in% c("RNA", "SCT") && layer == "data") expm1 else identity
  # out <- parallel::mclapply(all_pairs[1:2], function(z) {
  #   out <- presto::wilcoxauc(X = fun(get_layer(obj, assay = assay, layer = layer)),
  #                            y = obj@meta.data[[group]],
  #                            groups_use = z) |>
  #     dplyr::select(-statistic, -pval)
  #   dplyr::left_join(out |> dplyr::filter(group == z[1]) |>
  #                      dplyr::rename("in" = group, "avgExpr_in" = avgExpr),
  #                    out |> dplyr::filter(group == z[2]) |>
  #                      dplyr::rename("out" = group, "avgExpr_out" = avgExpr) |>
  #                      dplyr::select(feature, out, avgExpr_out),
  #                    by = "feature")
  # }, mc.cores = mc.cores)
  #
}

make_list_of_lists_redundant <- function(x) {
  # assumes list depth of 2
  all_names <- unique(c(names(x), unlist(purrr::map(x, names))))
  sep <- brathering::find_sep(all_names)
  flat <- purrr::list_flatten(x, name_spec = paste0("{outer}", sep, "{inner}"))
  flat2 <- flat
  rev_names <- purrr::map_chr(purrr::map(strsplit(names(flat2), sep), rev), paste, collapse = sep)
  names(flat2) <- rev_names
  flat_join <- c(flat, flat2)

  return(brathering::list_unflatten(flat_join, sep = sep))
}


#
# out <- presto::wilcoxauc(X = get_layer(so, assay = "RNA", layer = "data"),
#                          y = so@meta.data$disease_main,
#                          groups_use = c("AKI", "SLE"))
#
# # pairwise
# lvls <- as.character(unique(so@meta.data$disease_main))
# all_pairs <- combn(lvls, 2, simplify = F)
#
# # one vs rest
# lvls <- as.character(unique(so@meta.data$disease_main))
# all_onevsrest <- purrr::map_dfc(lvls, function(x) {
#   not_lvl <- brathering::random_varname(lvls, prefix = paste0("not_", x), n_digits = 0)
#   stats::setNames(data.frame(ifelse(so@meta.data$disease_main == x, x, not_lvl)), x)
# })
#
#
#
# sosub <- subset(so, subset = final_annotation == "B_lymphocytes")
# Seurat::Idents(sosub) <- sosub$disease_main
#
# out <- find_markers_pairwise(sosub, group = "disease_main")
# zz0 <- Seurat::FindMarkers(sosub, ident.1 = "SLE", ident.2 = "AKI", test.use = "MAST")
#
# all <- Seurat::FindAllMarkers(sosub, group.by = "disease_main")
#
# zz <- avg_log2fc_seuratlike(tt$AKI, tt$SLE)
# zz <- stack(zz)
#
#
# find_markers <- function(mat, group, groups_use = NULL) {
#   presto::wilcoxauc(X = mat,
#                     y = group,
#                     groups_use = groups_use)
# }

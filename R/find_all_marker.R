#' Calculate markers with presto
#'
#'  Column logFC is difference in group means, not logFC. But avg_log2FC
#'  in Seurat style can be added.
#'
#' @param obj seurat object
#' @param meta_col column name in meta.data
#' @param assay assay in obj
#' @param calc_avglog2fc calculate avglog2fc like Seurat, prestos logFC is diff
#' of means (wrong naming)
#' @param eps pseudo count to avoid div by zero
#' @param features features selection to increase speed further
#' @param layer which layer to pull expression matrix from
#' @param mc.cores multicores for calc_avglog2fc
#'
#' @returns data frame with DEG
#' @export
#'
#' @examples
find_all_marker <- function(obj,
                            meta_col = NULL,
                            assay = "RNA",
                            layer = "data",
                            calc_avglog2fc = T,
                            eps = 1e-6,
                            features = NULL,
                            na_rm = F,
                            mc.cores = 1) {

  obj <- scexpr:::check.SO(
    obj,
    assay = assay)

  ## get shared features with expr>0
  presto_feat <- purrr::map(obj, ~Matrix::rowSums(get_layer(
    obj = .x,
    layer = layer,
    assay = assay
  )))

  common_feat <- purrr::reduce(purrr::map(presto_feat, names), intersect)
  presto_feat <- purrr::map(presto_feat, ~.x[common_feat])
  sums <- Matrix::colSums(do.call(rbind, presto_feat))
  presto_feat <- names(which(sums>0))

  if (is.null(meta_col)) {
    message("using Idents.")
    obj <- purrr::map(obj, ~SeuratObject::AddMetaData(.x, Seurat::Idents(.x), col.name = "meta_col"))
    meta_col <- "meta_col"
  }

  # presto_feat <- unique(unlist(purrr::map(obj, ~names(which(Matrix::rowSums(get_layer(
  #   obj = .x,
  #   layer = layer,
  #   assay = assay
  # )) > 0)))))


  if (!is.null(features)) {
    features <- scexpr:::check.features(SO = obj, features = unique(features), meta.data = F)
    if (any(!features %in% presto_feat)) {
      message("No expressers found for: ", paste(features[which(!features %in% presto_feat)], collapse = ","), ". Will not be plotted.")
    }
    presto_feat <- intersect(presto_feat, features) # change that?
  }

  X = Gmisc::fastDoCall(cbind, # get_layer orders features, so cbind works
                        purrr::map(obj, ~get_layer( # expm1
                          obj = .x,
                          layer = layer,
                          assay = assay,
                          features = presto_feat)))

  y <- unlist(purrr::map2(obj, meta_col, ~.x@meta.data[,.y,drop=T]), use.names = F)

  if (anyNA(y)) {
    if (na_rm) {
      X <- X[, which(!is.na(y))]
      y <- y[which(!is.na(y))]
    } else {
      y[which(is.na(y))] <- "NA"
    }
  }

  wil_auc_raw <- presto::wilcoxauc(X = X, y = y)

  if (calc_avglog2fc) {
    groups <- as.character(unique(y))
    # logfc calc

    # logfcs <- purrr::map_dfr(purrr::set_names(groups), function(z) {
    #   y <- as.character(y)
    #   z <- as.character(z)
    #   newname <- paste0("not_", z) # not used before?
    #   y[which(y != z)] <- newname
    #   X2 <- brathering::split_mat(X, y, byrow = F)
    #   X2 <- purrr::map(X2, Matrix::rowMeans)
    #   avglog2fc <- log2(X2[[z]])-log2(X2[[newname]])
    #
    #   avglog2fc <- stack(avglog2fc)
    #   names(avglog2fc) <- c("avg_log2FC", "feature")
    #   avglog2fc$feature <- as.character(avglog2fc$feature)
    #
    #   return(avglog2fc)
    # }, .id = "group")


    logfcs <- parallel::mclapply(groups, function(g) {
      idx_g <- y == g
      idx_bg <- y != g

      mu_g <- Matrix::rowMeans(X[, idx_g, drop = FALSE])
      mu_bg <- Matrix::rowMeans(X[, idx_bg, drop = FALSE])

      avg_log2FC <- log2((mu_g + eps) / (mu_bg + eps))

      df <- tibble::tibble(
        group = g,
        feature = rownames(X),
        avg_log2FC = unname(avg_log2FC)
      )

      return(df)
    }, mc.cores = mc.cores)

    wil_auc_raw <- dplyr::left_join(wil_auc_raw,
                                    dplyr::bind_rows(logfcs),
                                    by = dplyr::join_by(group, feature))
  }

  return(wil_auc_raw)
}

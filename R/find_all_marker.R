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
#' @param na_rm Logical. If `TRUE`, exclude cells with missing group labels.
#' Otherwise, keep them in a separate group labeled `"NA"` or a suffixed
#' variant if that label is already in use.
#'
#' @returns data frame with DEG
#' @export
#'
#' @examples
#' \dontrun{
#' marker <- scexpr::find_all_marker(so, meta_col = "orig.ident")
#' }
find_all_marker <- function(obj,
                            meta_col = NULL,
                            assay = "RNA",
                            layer = "data",
                            calc_avglog2fc = T,
                            eps = 1e-6,
                            features = NULL,
                            na_rm = F,
                            mc.cores = 1) {

  scexpr:::.ensure_package("presto")

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

  if (!length(presto_feat)) {
    stop("No expressed features remain for marker calculation.", call. = FALSE)
  }

  layer_matrices <- purrr::map(obj, ~get_layer(
    obj = .x,
    layer = layer,
    assay = assay,
    features = presto_feat
  ))
  group_labels <- purrr::map2(obj, meta_col, function(so, column) {
    metadata <- so@meta.data
    stats::setNames(metadata[, column, drop = TRUE], rownames(metadata))
  })
  layer_labels <- purrr::map2(group_labels, layer_matrices, function(labels, mat) {
    cells <- colnames(mat)
    if (is.null(cells) || any(!cells %in% names(labels))) {
      stop("Could not match assay layer cells to metadata.", call. = FALSE)
    }
    labels[cells]
  })

  # Unlisting a mixture of factors and strings would use the factor codes.
  factor_labels <- vapply(layer_labels, is.factor, logical(1))
  if (any(factor_labels) && !all(factor_labels)) {
    layer_labels[factor_labels] <- lapply(layer_labels[factor_labels], as.character)
  }
  y <- unlist(layer_labels, use.names = FALSE)
  if (is.factor(y) && anyNA(levels(y))) {
    y <- droplevels(y, exclude = NA)
  }
  X <- Gmisc::fastDoCall(cbind, layer_matrices)
  rm(layer_matrices)

  if (anyNA(y)) {
    labeled <- !is.na(y)
    if (na_rm) {
      X <- X[, labeled, drop = FALSE]
      y <- y[labeled]
    } else {
      # Reserve labels from metadata even when their cells are outside the layer.
      used_labels <- unique(unlist(lapply(group_labels, as.character), use.names = FALSE))
      missing_group <- utils::tail(make.unique(c(used_labels[!is.na(used_labels)], "NA")), 1L)
      if (is.factor(y)) {
        levels(y) <- union(levels(y), missing_group)
      }
      y[!labeled] <- missing_group
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

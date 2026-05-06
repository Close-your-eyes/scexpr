#' Find marker genes in expressing cells; dichotomous
#'
#' @param obj seurat object or gene x cell matrix (cells as columns)
#' @param group meta column in object or group vector of length ncol(obj);
#' must have two levels
#' @param get_layer_args args to get_layer
#' @param mc.cores parallel computing: mc.cores
#'
#' @returns
#' @export
#'
#' @examples
find_marker_in_expr_dich <- function(obj,
                                     group,
                                     get_layer_args = list(),
                                     mc.cores = 4) {

  # dich  =  dichotomous

  if (length(group) == 1 && !methods::is(obj, "Seurat")) {
    stop("when obj is a matrix, group has to be a vector.")
  } else if (length(group) == 1) {
    group <- obj@meta.data[[group]]
  }
  group <- as.character(group)

  if (methods::is(obj, "Seurat")) {
    obj <- Gmisc::fastDoCall(get_layer, args = c(list(obj = obj), get_layer_args))
  }

  if (ncol(obj) != length(group)) {
    stop("column number of obj not equal to length of group.")
  }

  if (length(unique(group)) != 2) {
    stop("group must have two levels.")
  }

  obj <- brathering::split_mat(x = obj, f = group, byrow = F)

  idx <- stats::setNames(seq_len(nrow(obj[[1]])), rownames(obj[[1]]))
  out <- parallel::mclapply(idx, function(i) {

    v1 <- obj[[1]][i, ]
    v2 <- obj[[2]][i, ]

    z1 <- v1[v1 > 0]
    z2 <- v2[v2 > 0]

    if (length(z1) == 0L || length(z2) == 0L) return(NULL)
    m1 <- mean(z1)
    m2 <- mean(z2)

    df <- data.frame(
      pval = wilcox.test(z1, z2, exact = FALSE)$p.value,
      mean_diff = m1 - m2,
      log2fc = log2(m1) - log2(m2),
      pct.1 = length(z1)/length(v1),
      pct.2 = length(z2)/length(v2)
    )
    return(df)
  }, mc.cores = mc.cores)

  out <- purrr::discard(out, is.null)
  out <- dplyr::bind_rows(out, .id = "feature")
  out$padj <- stats::p.adjust(out$pval, method = "BH")
  out$group.1 <- names(obj)[1]
  out$group.2 <- names(obj)[2]

  return(out)
}


#' Find marker genes in expressing cells
#'
#' @param obj seurat object or gene x cell matrix (cells as columns)
#' @param group meta column in object or group vector of length ncol(obj)
#' @param get_layer_args args to get_layer
#' @param mc.cores parallel computing: mc.cores
#'
#' @returns
#' @export
#'
#' @examples
find_marker_in_expr_all <- function(obj,
                                    group,
                                    get_layer_args = list(),
                                    mc.cores = 4) {


  if (length(group) == 1 && !methods::is(obj, "Seurat")) {
    stop("when obj is a matrix, group has to be a vector.")
  } else if (length(group) == 1) {
    group <- obj@meta.data[[group]]
  }
  group <- as.character(group)

  if (methods::is(obj, "Seurat")) {
    obj <- Gmisc::fastDoCall(get_layer, args = c(list(obj = obj), get_layer_args))
  }

  if (ncol(obj) != length(group)) {
    stop("column number of obj not equal to length of group.")
  }

  out <- purrr::map_dfr(unique(group), function(x) {
    group_dich <- ifelse(group == x, x, paste0("not_", x))
    y <- find_marker_in_expr_dich(
      obj,
      group = group_dich,
      get_layer_args = get_layer_args,
      mc.cores = mc.cores
    )
    if (y$group.1[1] != x) {
      y <- dplyr::rename(y, "group.1" = group.2, "group.2" = group.1)
    }
    return(y)
  })

  return(out)
}


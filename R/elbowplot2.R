#' Plot variance along PCs
#'
#' @param obj single seurat or list thereof
#' @param assay assay
#' @param npcs how many PC to print
#' @param theme ggplot theme
#'
#' @returns list of data and plots
#' @export
#'
#' @examples
elbowplot2 <- function(obj,
                       assay = "RNA",
                       npcs = min(30, length(Seurat::VariableFeatures(obj[[1]]))),
                       theme = ggplot2::theme_grey()) {

  if (!is.list(obj)) {
    obj <- list("1" = obj)
  } else {
    if (is.null(names(obj))) {
      names(obj) <- as.character(seq_along(obj))
    }
  }

  out <- purrr::map(names(obj), function(x) {

    pca <- tryCatch(expr = {
      stats::prcomp(scexpr::get_layer(obj[[x]],
                                      assay = assay,
                                      layer = "scale.data",
                                      transpose = T),
                    center=F,
                    scale.=F)
    },
    error = function(err) {
      stats::prcomp(scexpr::get_layer(obj[[x]],
                                      assay = assay,
                                      layer = "data",
                                      features = Seurat::VariableFeatures(obj[[x]]),
                                      transpose = T),
                    center=T,
                    scale.=T)
    })

    stdev <- pca$sdev
    var <- stdev^2
    df <- data.frame(pc = seq_along(stdev),
                     stdev = stdev,
                     var = var,
                     var_pct = var/sum(var)*100,
                     var_cum = cumsum(var),
                     var_cum_pct = cumsum(var)/sum(var)*100) |>
      dplyr::filter(pc <= npcs)

    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = pc, y = var)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_sqrt(sec.axis = ggplot2::sec_axis(~ . / sum(var)*100, name = "var_pct")) +
      ggplot2::scale_x_continuous() +
      theme +
      ggplot2::labs(subtitle = x)

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = pc, y = var_cum)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_sqrt(sec.axis = ggplot2::sec_axis(~ . / sum(var)*100, name = "var_cum_pct")) +
      ggplot2::scale_x_continuous() +
      theme

    p12 <- patchwork::wrap_plots(p1,p2,ncol = 1, axes = "collect_x")

    return(list(data = df, plot = p12))
  })

  plist <- sapply(out, "[", 2)
  df <- dplyr::bind_rows(sapply(out, "[", 1), .id = "obj")

  return(list(data = df, plot = patchwork::wrap_plots(plist, nrow = 1, axis_titles = "collect_y")))

  # fit <- nls.multstart::nls_multstart(
  #   var ~ a * exp(-b * pc) + c,
  #   data = df,
  #   iter = 200,
  #   start_lower = list(a = 0, b = 1e-4, c = 0),
  #   start_upper = list(a = max(df$var), b = 1, c = median(tail(df$var, 10))),
  #   lower = c(a = 0, b = 0, c = 0)
  # )
  # predict(fit, 1:100)
  #
  # # derivative
  # a <- coef(fit)["a"]
  # b <- coef(fit)["b"]
  # c <- coef(fit)["c"]
  # df$var_fit <- predict(fit)
  # df$var_deriv <- -a * b * exp(-b * df$pc)

  # model can be improved
  #   geom_function(fun = function(pc) eval(parse(text = get_nls_fun(fit))), color = "blue", linewidth = 1)

  # if (length(obj) > 1) {
  #   mapping <- ggplot2::aes(color = obj)
  # } else {
  #   mapping <- NULL
  # }
  #
  # p1 <- ggplot2::ggplot(df, ggplot2::aes(x = pc, y = var)) +
  #   ggplot2::geom_point(mapping = mapping) +
  #   ggplot2::scale_y_sqrt(sec.axis = ggplot2::sec_axis(~ . / sum(var)*100, name = "var_pct")) +
  #   ggplot2::scale_x_continuous(limits = c(1,npcs)) +
  #   theme +
  #   ggplot2::facet_wrap(vars(obj), scales = "free_y")
  #
  # p2 <- ggplot2::ggplot(df, ggplot2::aes(x = pc, y = var_cum)) +
  #   ggplot2::geom_point(mapping = mapping) +
  #   ggplot2::scale_y_sqrt(sec.axis = ggplot2::sec_axis(~ . / sum(var)*100, name = "var_cum_pct")) +
  #   ggplot2::scale_x_continuous(limits = c(1,npcs)) +
  #   theme
  #ggplot2::scale_x_log10(breaks = c(1,5,10,12,15,20,25,30,40,50,100), limits = c(1,npcs))


  # return(patchwork::wrap_plots(p1,p2,ncol = 1, axes = "collect_x"))

}

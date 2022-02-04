#' Title
#'
#' @param d
#' @param feat
#' @param pt.size
#' @param geom1
#' @param geom2
#' @param font.size
#' @param filter.non.expr
#' @param plot.expr.freq
#' @param label.size
#' @param ngc
#' @param pgc
#' @param ngn
#' @param pgn
#'
#' @return
#' @export
#'
#' @examples
.expr_jitter <- function (d,
                          feat,
                          pt.size,
                          geom1 = c("jitter", "point", "dotplot"),
                          geom2 = c("none", "violin", "boxplot"),
                          font.size = 14,
                          filter.non.expr = F,
                          plot.expr.freq = F,
                          label.size = 4,
                          ngc,
                          pgc,
                          ngn,
                          pgn) {

  #geom <- match.arg(geom, c("jitter", "point"))
  #my_geom <- switch(geom, jitter = geom_jitter, point = geom_point)
  #my_geom(...)

  geom1 <- match.arg(geom1, c("jitter", "point", "dotplot"))
  geom2 <- match.arg(geom2, c("violin", "boxplot", "none"))
  ccc <- stats::setNames(list(ngc, pgc), nm = c(ngn, pgn))

  dd <- do.call(rbind,
                lapply(names(ccc),
                       function(x) data.frame(stats::setNames(reshape2::melt(as.matrix(d[feat,ccc[[x]],drop=F]),value.name = "expr")[,-2], nm = c("Feature", "expr")), group = x)))

  if (plot.expr.freq) {
    stat <-
      dd %>%
      dplyr::group_by(Feature) %>%
      dplyr::mutate(max.feat.expr = max(expr)) %>%
      dplyr::group_by(Feature, group, max.feat.expr) %>%
      dplyr::summarise(pct.expr = sum(expr > 0)/dplyr::n(), .groups = "drop")
    '    stat1 <-
      dd %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarise(max.feat.expr = max(expr), .groups = "drop")

    stat2 <- do.call(rbind,lapply(as.character(unique(dd$Feature)), function(x) {
      do.call(rbind,lapply(unique(dd$group), function(y) {
        print(x)
        print(y)
        print(sum(dd[intersect(which(dd$Feature == x), which(dd$group == y)),"expr"] > 0))
        print(length(intersect(which(dd$Feature == x), which(dd$group == y))))
        data.frame(Feature = x,
                   group = y,
                   pct.expr = sum(dd[intersect(which(dd$Feature == x), which(dd$group == y)),"expr"] > 0)/length(intersect(which(dd$Feature == x), which(dd$group == y)))
        )
      }))
    }))

    stat <- dplyr::left_join(stat2, stat1, by = "Feature")'

  }

  if (filter.non.expr) {
    dd <- dd[which(dd$expr > 0),]
  }

  if (geom2 != "split violin") {
    plot <- ggplot2::ggplot(dd, ggplot2::aes(x = group, y = expr))
    if (geom1 == "jitter") {
      plot <- plot + ggplot2::geom_jitter(width = 0.2, size = pt.size, color = "tomato2")
    }

    if (geom1 == "point") {
      plot <- plot + ggplot2:: geom_point(size = pt.size, color = "tomato2")
    }

    if (geom1 == "dotplot") {
      plot <- plot + ggplot2::geom_dotplot(binaxis = "y", binwidth = (max(dd$expr) - min(dd$expr))/40, stackdir = "center", fill = "tomato2", stackratio = 0.7)
    }

    if (geom2 == "violin") {
      plot <- plot + ggplot2::geom_violin(alpha = 0.5)
    }
    if (geom2 == "boxplot") {
      plot <- plot + ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA)
    }

    if (plot.expr.freq) {
      plot <- plot + ggplot2::geom_text(data = stat, ggplot2::aes(label = round(pct.expr, 2), y = max.feat.expr + 0.4), size = label.size)
      #expand_limits
    }

    plot <-
      plot +
      ggplot2::theme_bw(base_size = font.size) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), legend.position = "none", panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), text = ggplot2::element_text(family = "Courier"))

  } else {
    plot <- ggplot2::ggplot(dd, ggplot2::aes(x = 0, y = expr, fill = group)) +
      .geom_split_violin() +
      ggplot2::theme_bw(base_size = font.size) +
      ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "right", legend.title = ggplot2::element_blank(),  panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "white"), text = ggplot2::element_text(family = "Courier"))
  }

  plot <- plot + ggplot2::facet_wrap(ggplot2::vars(Feature), scales = "free_y")

  return(plot)
}


.GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", #GeomViolin
                                     draw_group = function(self, data, ..., draw_quantiles = NULL) {
                                       data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                                       grp <- data[1, "group"]
                                       newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                                       newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                                       newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                                       if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                         stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                                   1))
                                         quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                         aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                         aesthetics$alpha <- rep(1, nrow(quantiles))
                                         both <- cbind(quantiles, aesthetics)
                                         quantile_grob <- GeomPath$draw_panel(both, ...)
                                         ggplot2:::ggname(".geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                                       }
                                       else {
                                         ggplot2:::ggname(".geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                                       }
                                     })
.geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                               draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                               show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = .GeomSplitViolin,
                 position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

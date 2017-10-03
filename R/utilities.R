#' pal
#'
#' @param  col A vector of hexcodes
#'
#' @return This function displays a colour pallete, you can change the colour of the border but its probably not neccessary.
#' Written by Achim Zeileis (in the colorspace package apparently) but I found it here:
#' http://www.r-bloggers.com/the-paul-tol-21-color-salute/
#' @examples pal(g_colours)
#' @export
pal <- function(col, border = "light gray", ...){

  n <- length(col)

  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)

  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)

}

#' pickRandomRows
#'
#' @param data A dataframe
#' @param numberOfRows Number of rows to sample
#'
#' This function randomly selects rows from a data frame. I got it from here/; http://www.markhneedham.com/blog/2014/11/26/r-dplyr-select-random-rows-from-a-data-frame/
#' @export
pickRandomRows = function(data, numberOfRows = 10) {

  require(dplyr)
  require(magrittr)

  df %>% slice(runif(numberOfRows,0, length(data[,1])))
}

#' read.nano
#'
#' @param Path The path to the Nanophotometer output file
#' @param n_samples The number of Nanophotometer measurements (i.e. if you did one twice count it twice)
#'
#' This reads the excel file produced by the Nanophotometer and outputs a markdown table.
#' @examples Path = "~/Documents/Nanophotometer/21_7_16_.xlsx"
#' read_nano(Path, 7)
#' @export
read.nano <- function(Path){

  require(readxl)
  require(dplyr)
  require(magrittr)
  require(knitr)

  dat <- read_excel(Path, skip = 20)
  dat <- dplyr::select(dat, 2,3,11,12)
  dat <- magrittr::set_colnames(dat, c("Sample", "Concentration", "A260/A280", "A260/A230"))

  kable(dat, format = "markdown", row.names = FALSE,  padding = 2, align = 'c')
}

#' moe
#'
#' @param data A vector of numeric type
#' @param alpha The significance level (default = 0.05)
#'
#' This calculates a margin of error (half a confidence interval) with the default being 95%.
#'
#' @export
moe <- function(data, alpha = 0.05){
  moe <- qt(1 - (alpha / 2), sum(!is.na(data)) - 1) * (sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data))))
  return(moe)
}

#' px_to_mm
#'
#' @export
px_to_mm <- function(n){
  n*0.264583
}


ggarrange <- function(..., plotlist = NULL, ncol = NULL, nrow = NULL,
                      labels = NULL,
                      align = c("none", "h", "v", "hv"),
                      widths = 1, heights = 1,
                      legend = NULL, common.legend = FALSE,
                      label_size = 14, label_fontfamily = NULL,
                      label_fontface = "bold", label_colour = NULL,
                      label_x = 0, label_y = 1,
                      hjust = -0.5, vjust = 1.5)
{
  align <- match.arg(align)
  plots <- c(list(...), plotlist)
  nb.plots <- length(plots)
  nb.plots.per.page <- .nbplots_per_page(ncol, nrow)

  if(is.null(legend) & common.legend)
    legend <- "top"
  legend <- .check_legend(legend)
  if(!is.null(legend))
    plots <- purrr::map(plots, function(x) x + theme(legend.position = legend))

  # Split plots over multiple pages
  if(nb.plots > nb.plots.per.page){
    plots <- split(plots, ceiling(seq_along(plots)/nb.plots.per.page))
  }

  # One unique page
  else plots <- list(plots)

  res <- purrr::map(plots, .plot_grid,
                    ncol = ncol, nrow = nrow, labels = labels, align = align,
                    rel_widths = widths, rel_heights = heights,
                    legend = legend, common.legend = common.legend,
                    label_size = label_size, label_fontfamily = label_fontfamily,
                    label_fontface = label_fontface, label_colour = label_colour,
                    label_x = label_x, label_y = label_y,
                    hjust = hjust, vjust = vjust)

  if(length(res) == 1) res <- res[[1]]

  class(res) <- c(class(res), "ggarrange")
  res
}




# Compute number of plots per page
.nbplots_per_page <- function(ncol = NULL, nrow = NULL){

  if(!is.null(ncol) & !is.null(nrow))
    ncol * nrow
  else if(!is.null(ncol))
    ncol
  else if(!is.null(nrow))
    nrow
  else Inf
}


.plot_grid <- function(plotlist, legend = "top", common.legend = FALSE,  ... ){


  if(common.legend){
    # Legend infos
    leg <- get_legend(plotlist[[1]])
    lheight <- sum(leg$height)
    lwidth <- sum(leg$width)
    plotlist <- purrr::map(plotlist, function(x) x + theme(legend.position = "none"))
  }

  res <- cowplot::plot_grid(plotlist = plotlist, ...)
  if(!common.legend) return(res)

  arrangeGrob <- gridExtra::arrangeGrob
  unit.c <- grid::unit.c
  .unit <- grid::unit(1, "npc")

  res <- switch(legend,
                top = arrangeGrob(leg, res, ncol = 1,
                                  heights = unit.c(lheight, .unit - lheight)),
                bottom = arrangeGrob(res, leg, ncol = 1,
                                     heights = unit.c(unit(1, "npc") - lheight, lheight)),
                left = arrangeGrob(leg, res, ncol = 2,
                                   widths = unit.c(lwidth, .unit - lwidth)),
                right = arrangeGrob(res, leg, ncol = 2,
                                    widths = unit.c(.unit - lwidth, lwidth))
  )

  p <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(res))
  p

}

.check_legend <- function(legend){

  allowed.values <- c("top", "bottom", "left", "right", "none")

  if(is.null(legend) | is.numeric(legend))
    return(legend)
  else if(is.logical(legend)){
    if(legend) legend <- "top"
    else legend <- "none"
  }
  else if(is.character(legend)){
    legend <- legend[1]
    if(!legend %in% allowed.values)
      stop("Argument legend should be one of ", .collapse(allowed.values, sep = ", "))
  }
  return (legend)
}




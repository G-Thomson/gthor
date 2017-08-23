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




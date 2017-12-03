#' ggtimeline
#'
#'
#' @param x is a data.frame with three columns
#' @param begin is the date the timeline begins in the form "dd/mm/yyyy"
#' @param end is the date the timeline ends in the form "dd/mm/yyyy"
#' @param group Boolean denoting whether activities should be grouped in terms of colour
#'
#' @return  This function draws a timeline given a data.frame with three or four columns.
#' These columns should correspond to activity, start and end. The start and end columns
#' should be dates in the form of "dd/mm/yyyy". You also need to specific when
#' you want the time table to begin and end. Only have a "group" column if you have group = TRUE.
#' @examples timeline(timeline_dates, "1/1/2015", "1/3/2018")
#' @export
ggtimeline <- function(x, begin, end, group = F){

  if(group == F){
    if(ncol(x) != 3){stop("\n\nToo many columns!! You can only have three.\n\n")}
    colnames(x) <- c("activity", "start", "halt")
  }
  else{
    if(!"group" %in% colnames(x))stop("\n\nIf you want to group tasks you need a column named 'group'.\n\n")
    colnames(x) <- c("activity", "start", "halt", "group")
  }

  require(scales)
  require(ggplot2)

  x$start = as.Date(x$start, "%d/%m/%Y")
  x$halt = as.Date(x$halt, "%d/%m/%Y")

  if(group == F){
    chart <- ggplot(x, aes(start, activity, colour = activity))
  }
  else{
    chart <- ggplot(x, aes(start, activity, colour = group))
  }

chart +
  geom_errorbarh(aes(xmin = start, xmax = halt, height = 0.4), size = 3) +
      scale_x_date(labels = date_format("%b %Y"),
                   limits = c(as.Date(begin, "%d/%m/%Y"), as.Date(end, "%d/%m/%Y")),
                   breaks = seq(as.Date(begin, "%d/%m/%Y"), as.Date(end, "%d/%m/%Y"), by = '2 month')) +
      theme(legend.position = "none",
            axis.text = element_text(size = 16),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()) +
      scale_color_manual(values = g_colours)

}

#' ggplotRegression
#'
#' @param fit is a linear model from the lm function
#'
#' @return This plots a linear model with ggplot2. It is a function written by Susan
#' Johnston (https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/).
#' @export
ggplotRegression <- function (fit) {

  require(ggplot2)

  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))

}

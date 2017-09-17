#' my_theme
#'
#' @return This is the theme I use for my plots. Simple and based on theme_wsj().
#' Use it with theme_set or as a parameter in a ggplot call.
#' @examples theme_set(my_theme)
#' @export
my_theme <- function(base_size = 12,
                     base_family = "sans") {
  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
      theme(line = element_line(linetype = 1, colour = "black"),
            rect = element_rect(fill = "transparent", linetype = 0, colour = NA),
            text = element_text(colour = "black"),
            title = element_text(size = rel(3)),
            axis.title = element_text(vjust = 4),
            axis.title.x = element_blank(),
            axis.title.y = element_text(vjust = 4),
            axis.text = element_text(face = "bold", size = rel(2)),
            axis.text.x = element_text(colour = NULL),
            axis.text.y = element_text(colour = NULL),
            axis.ticks = element_line(colour = NULL),
            axis.ticks.x = element_line(colour = NULL),
            axis.line = element_line(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            panel.grid = element_line(colour = NULL, linetype = 3),
            panel.grid.major = element_line(colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA),
            plot.title = element_text(hjust = 0, face = "bold"),
            plot.margin = unit(c(1, 1, 1, 1), "lines"),
            strip.background = element_rect(fill = "transparent", colour = NA),
            strip.text = element_text(hjust=0, face = "bold", size = rel(2))))
}

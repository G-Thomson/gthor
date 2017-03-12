#' flr.summary
#'
#' @param data This is a dataframe of flowering time fata with the following columns; Sowing.Number,
#' Conditions, Genotype, Plant.Number, Date.Sown, Date.Flowered, Number.of.Days.to.First.Floral.Bud,
#' Position.of.1st.Floral.Bud, Node.Number.of.Main.Axis.at.1st.Floral.Bud and Comments.
#' @param LoI These are the Lines of Interest to be examined. Enter these as a vector of sowing
#'  numbers as characters.
#'
#' @return This function summarises a flowering time data frame producing a list of the
#' summarised data and relevant raw data
#' @examples data = Mid_2016
#' LoI = c("P912", "P913") # Lines of interest
#' d <- flr.summary(data, LoI)
#' @export
flr.summary <- function(data, LoI){

  require(dplyr)

  colnames(data) <- make.names(colnames(data))

  # Rename data
  data <- rename(data, SN = Sowing_Number,
                 Conditions = Conditions,
                 Genotype = Genotype,
                 Plant_num = Plant_Number,
                 Sowing_date = Date_Sown,
                 Flowering_date = Date_Flowered,
                 Days_to_flower = Number_of_Days_to_First_Floral_Bud,
                 First_bud_loc = Position_of_1st_Floral_Bud,
                 Node_num = Node_Number_of_Main_Axis_at_1st_Floral_Bud,
                 Mutant = Mutant)

  # Relevent summary data
  rel.sum <- filter(data, SN %in% LoI) %>%
    select(SN, Days_to_flower, Node_num) %>%
    group_by(SN) %>%
    summarise_each(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE), len = length)) %>%
    mutate(Day_se = Days_to_flower_sd/sqrt(Days_to_flower_len),
           Node_se = Node_num_sd/sqrt(Node_num_len))

  Gen.sum <- distinct(data, SN, .keep_all = TRUE) %>%
    select(SN, Genotype, Conditions)

  rel.sum <- inner_join(rel.sum, Gen.sum)

  # Relevent raw data
  rel.raw <- filter(data, SN %in% LoI) %>%
    select(SN, Genotype, Conditions, Plant_num, Days_to_flower, Node_num, Mutant) %>%
    filter(!is.na(Days_to_flower))

  res = list(Summary = rel.sum, Raw = rel.raw)
  return(res)
}

#' flr.plot
#'
#' @param data This is a dataframe from flr.summary
#'
#' @return This function plots flowering time data. If you a just comparing one line and its wild-type control it
#'  will colour and plot these side by side. If there are multiple conditions then it will colour these and facet by genotype.
#' @examples flr.plot(data)
#' @export
flr.plot <- function(data){

  require(tidyverse)
  require(cowplot)
  require(ggthemes)

  # Number of unique conditions
  cond_cnt <- data[["Summary"]] %>%
    select_("Conditions") %>%
    unique(.) %>%
    length(.)

  # Plot y axis limits
  flr_ax_lim = data[["Raw"]]  %>%
    select_("Days_to_flower")  %>%
    max() %>%
    plyr::round_any(10, ceiling)
  node_ax_lim = data[["Raw"]]  %>%
    select_("Node_num")  %>%
    max() %>%
    plyr::round_any(5, ceiling)

  if(cond_cnt == 1){
    x = "Genotype"
    Colour = "Genotype"
  }
  else{
    x = "Conditions"
    Colour = "Conditions"
    # facet
  }

  Flr <-  ggplot(data[["Summary"]], aes_string(x = x, y = "Days_to_flower_mean", fill = Colour)) +
            geom_point(data = data[["Raw"]], aes_string(x = x, y = "Days_to_flower")) +
            geom_bar(data = data[["Summary"]], stat = "identity") +
            geom_errorbar(aes(ymin=Days_to_flower_mean-Day_se,
                              ymax=Days_to_flower_mean+Day_se), width=.2)  +
            theme_wsj() +
            theme(legend.position="none",
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  strip.background = element_rect(fill = "transparent",colour = NA),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(),
                  axis.title.x = element_text(vjust = -2),
                  axis.title.y = element_text(vjust = 4)) +
            scale_fill_manual(values = g_colours) +
            labs(y = "Mean days to flower") +
            scale_y_continuous(expand=c(0,0), limits = c(0, flr_ax_lim)) +
            ggtitle("Days to Flower")

  Node <- ggplot(data[["Summary"]], aes_string(x = x, y = "Node_num_mean", fill = Colour)) +
    geom_point(data = data[["Raw"]], aes_string(x = x, y = "Node_num")) +
    geom_bar(data = data[["Summary"]], stat = "identity") +
    geom_errorbar(aes(ymin=Node_num_mean - Node_se,
                      ymax=Node_num_mean + Node_se), width=.2)  +
    theme_wsj() +
    theme(legend.position="none",
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          strip.background = element_rect(fill = "transparent",colour = NA),
          axis.text = element_text(size = 14),
          axis.title = element_text(),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 4)) +
    scale_fill_manual(values = g_colours) +
    labs(y = "Mean nodes at time of flowering") +
    scale_y_continuous(expand=c(0,0), limits = c(0, node_ax_lim)) +
    ggtitle("Nodes at Flowering")

  if(cond_cnt > 1){
    Flr <- Flr + facet_grid(.~Genotype, scales="free")
    Node <- Node + facet_grid(.~Genotype, scales="free")
  }
  else{}

  plot_grid(Flr, Node)

}

#' flr.hist.graph
#'
#' @param data This is a a flr.summary object used to draw the graphs.
#' @param meas The measurement of flowering time. This is by default "fl" for days to flowering but could also be "nd" for nodes at the time of flowering.
#' @param new.titles By default this function labels each facet with the entry in the Genotype column of the flr.summary object. This can be changed by inputting a named vector of names which map the genotypes to the new names.
#'
#' @return This is a function for drawing a histogram of flowering time measurments from a screen. It takes a flr.summary object as input and outputs a graph for each line highlighting homozygotes. At this stage the Mutant column needs to be either "WT" or "Homo"
#'
#' @examples new_labs <- c(`NF5076` = "NF5076 - LD+V", `R108` = "R108 - LD+V")
#' @examples flr.hist.graph(res, new.titles = new_labs)
flr.hist.graph <- function(data, meas = "fl", new.titles = NULL){

  require(tidyverse)
  require(scales)
  require(ggthemes)

  # This takes the raw data from the flr.summary object and adds a colour column
  dat <- data$Raw %>% mutate(colr = ifelse(Genotype == "R108", "#00A08A",
                                           ifelse(Mutant == "Homo", "#F2AD00", "#FF0000")))

  # Creates the basic plots based on measurment of flowering
  if(meas == "nd"){
    plot <- ggplot(data=dat, aes(Node_num, fill = Mutant)) +
      geom_histogram(data=subset(dat, Genotype!="R108"),aes(y=..count..), binwidth = 1.5) +
      geom_histogram(data=subset(dat, Genotype=="R108"),aes(x = Node_num, y=..count..), binwidth = 1.5)
  }
  else if(meas == "fl"){
    plot <- ggplot(data=dat, aes(Days_to_flower, fill = Mutant)) +
      geom_histogram(data=subset(dat, Genotype!="R108"),aes(y=..count.., fill = colr), binwidth = 1.5) +
      geom_histogram(data=subset(dat, Genotype=="R108"),aes(x = Days_to_flower, y=..count.., fill = colr), binwidth = 1.5)
  }
  else{ stop("ERROR: Please specify a correct measurment of flowering. See docs for details.") }

  # Generic plot aesthetics
  plot <- plot +
    scale_fill_manual(values = g_colours[c(2,3,1)]) +
    theme_wsj() +
    theme(legend.position="none",
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.title = element_text(size = 20),
          axis.title.y = element_text(angle = 90),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          # panel.margin = unit(3, "lines"),
          strip.background = element_blank(),
          strip.text = element_text(hjust=0, size = 16)) +
    scale_y_continuous(expand=c(0,0), breaks=pretty_breaks())

  # Flowering or node specific graph elements - axis titles, facet titles
  if(meas == "nd"){
    plot <- plot + labs(x = "Nodes at Flowering", y = "Number of plants", size = 40) +
      theme(strip.text = element_text(hjust=0, size = 16)) +
      facet_wrap(~Genotype, nrow=length(unique(dat$Genotype)))
  }
  else if(meas == "fl"){
    plot <- plot + labs(x = "Days to Flowering", y = "Number of plants", size = 40) +
      theme(strip.text = element_text(hjust=0, size = 16))
  }

  # Faceting with or without new titles
  if(!is.null(new.titles)){

    if(mode(new.titles) != "character" | is.null(names(new.titles)) | length(new.titles) != length(unique(dat$Genotype))){
      stop("ERROR: new.titles must be a named character vector. See docs for details.")
    }

    else{
      plot <- plot + facet_wrap(~Genotype, nrow=length(unique(dat$Genotype)), labeller = as_labeller(new.titles)) #scales = "free_x"
    }

  }

  else{
    plot <- plot + facet_wrap(~Genotype, nrow=length(unique(dat$Genotype)))
  }

  # Calls plot
  plot
}

#' flr.hist
#'
#' @param data This is a a flr.summary object used to draw the graphs
#' @param new.titles By default this function labels each facet with the entry in the Genotype column of the flr.summary object. This can be changed by inputting a named vector of names which map the genotypes to the new names.
#'
#' @return This is a wrapper function for flr.hist.graph in that it calls it twice. Once for days to flower and once for nodes at the time of flowering. It then arranges them.
#' @examples new_labs <- c(`NF5076` = "NF5076 - LD+V", `R108` = "R108 - LD+V")
#' @examples flr.hist(res, new.titles = new_labs)
#' @export
flr.hist <- function(data, new.titles = NULL){

  require(cowplot)

  # Creates each graph
  flr <- flr.hist.graph(data, new.titles = new.titles)
  nd <- flr.hist.graph(data, meas = "nd", new.titles = new.titles)

  # Calls and plots them
  plot_grid(flr, nd, labels = "AUTO", align = 'h', rel_widths = c(1, 0.4))
}

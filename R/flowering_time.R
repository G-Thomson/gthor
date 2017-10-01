#' flr.test
#'
#' @export
flr.test <- function(data, LoI, group = conditions){

  require(dplyr)
  require(magrittr)
  require(broom)

  # Make group a quosure
  group <- enquo(group)

  # This just checks that the grouping argument is valid
  if(!quo_name(group) %in% c("sowing_number", "conditions", "mutant")) stop("You can only group by sowing_number or conditions")

  # Rename data
  data <- rename(data,
                 sowing_number = .data$Sowing_Number,
                 conditions = .data$Conditions,
                 genotype = .data$Genotype,
                 individual = .data$Individual,
                 plant_num = .data$Plant_Number,
                 sowing_date = .data$Date_Sown,
                 flowering_date = .data$Date_Flowered,
                 days_to_flower = .data$Number_of_Days_to_First_Floral_Bud,
                 first_bud_loc = .data$Position_of_1st_Floral_Bud,
                 node_num = .data$Node_Number_of_Main_Axis_at_1st_Floral_Bud,
                 mutant = .data$Mutant)

  # Ensure there are no NA entries
  data <- filter(data,
                 !is.na(flowering_date))

  # Ensure flowering time data is numeric
  data <- mutate(data,
                 days_to_flower = as.numeric(days_to_flower),
                 node_num = as.numeric(node_num))

  # Relevent summary data
  rel.sum <- filter(data, sowing_number %in% LoI) %>%
    select(sowing_number, genotype, conditions, plant_num, days_to_flower, node_num, mutant)

  # When grouping by condition I would like to preserve the sowing number in the
  # nested data so I can then extract it as a new column. This is not necessary
  # when grouping by sowing number as each has a single condition
  if(quo_name(group) == "conditions"){
    rel.sum <- rel.sum %>%
      select(-mutant) %>%
      group_by(!!group) %>%
      nest(sowing_number, plant_num, genotype, days_to_flower, node_num)  %>%
      mutate(sowing_numbers = map(data, ~ unique(.x$sowing_number)),
             data = map(data, ~ select(.x, -sowing_number)),

             data = map(data, ~ gather(.x, measurment, value, -genotype, -plant_num)),
             data = map(data, ~ unite(.x, geno_meas, genotype, measurment)),
             data = map(data, ~ group_by(.x, geno_meas)),
             data = map(data, ~ mutate(.x, counter = rank(as.numeric(plant_num), ties.method= "first"))),
             data = map(data, ~ select(.x, -plant_num)),
             data = map(data, ~ spread(.x, geno_meas, value))
      )

  }
  else if(quo_name(group) == "mutant"){
    rel.sum <- rel.sum %>%
      group_by(!!group, conditions) %>%
      nest(sowing_number, genotype, days_to_flower, node_num) %>%
      mutate(sowing_numbers = map(data, ~ unique(.x$sowing_number)))
  }

  # Data checks
  day_check <- rel.sum$data %>%
    map_lgl(~ ncol(select(.x, contains("days_to_flower"))) == 2)
  node_check <- rel.sum$data %>%
    map_lgl(~ ncol(select(.x, contains("node_num"))) == 2)

  if(length(day_check) != sum(day_check) | length(node_check) != sum(node_check)) stop("Check your data. There needs to be TWO sets of observations for each test")

  # Here I conduct t.tests between genotypes
  rel.sum <- rel.sum %>%
    mutate(test_days = map(data, ~ tidy(t.test(pull(.x[,2]), pull(.x[,4])))),
           test_days = map(test_days, ~ select(.x, estimate, statistic, p.value, conf.low, conf.high)),

           test_nodes = map(data, ~ tidy(t.test(pull(.x[,3]), pull(.x[,5])))),
           test_nodes = map(test_nodes, ~ select(.x, estimate, statistic, p.value, conf.low, conf.high)),

           days_p.value = map_dbl(test_days, ~.x$p.value),
           days_sig_diff = days_p.value < 0.05,
           nodes_p.value = map_dbl(test_nodes, ~.x$p.value),
           nodes_dig_diff = nodes_p.value < 0.05)

  return(rel.sum)
}


#' flr.summary
#'
#' @param data This is a dataframe of flowering time fata with the following columns; Sowing.Number,
#'  Conditions, Genotype (of parent seed), Individual, Plant.Number, Date.Sown, Date.Flowered,
#'  Number.of.Days.to.First.Floral.Bud, Position.of.1st.Floral.Bud,
#'  Node.Number.of.Main.Axis.at.1st.Floral.Bud, Comments and Mutant.
#'
#' @param LoI These are the Lines of Interest to be examined. Enter these as a vector of sowing
#'  numbers as characters.
#' @param group Data can be grouped either by by sowing number (sowing_number), growth conditions
#' (conditions) or mutantion genotype (mutant). Only these options are possible. When data is grouped by conditions the sowing
#' numbers column is a list in case multiple sowing numbers are included.
#'
#' @return This function summarises a flowering time data frame producing a dataframe of the
#' summarised data grouped either by sowing number or growth conditions. Means and 95% confidence
#' intervals are calculated. A list column containing the raw data and mutant status if
#' applicable is also included.
#'
#' @examples data = Mid_2016
#' # Lines of interest
#' LoI <- c("P880", "P881", "P883", "Q022", "Q023", "P935", "P938", "P936", "P939")
#' flr.summary(flowering_data , LoI)
#' @export
flr.summary <- function(data, LoI, group = sowing_number, WT = "R108"){

  require(dplyr)
  require(forcats)
  require(magrittr)

  # Make group a quosure
  group <- enquo(group)

  # This just checks that the grouping argument is valid
  if(!quo_name(group) %in% c("sowing_number", "conditions", "mutant")) stop("You can only group by sowing_number or conditions")

  # Rename data
  data <- rename(data,
                 sowing_number = .data$Sowing_Number,
                 conditions = .data$Conditions,
                 genotype = .data$Genotype,
                 individual = .data$Individual,
                 plant_num = .data$Plant_Number,
                 sowing_date = .data$Date_Sown,
                 flowering_date = .data$Date_Flowered,
                 days_to_flower = .data$Number_of_Days_to_First_Floral_Bud,
                 first_bud_loc = .data$Position_of_1st_Floral_Bud,
                 node_num = .data$Node_Number_of_Main_Axis_at_1st_Floral_Bud,
                 mutant = .data$Mutant)

  # Ensure there are no NA entries
  data <- filter(data,
                 !is.na(flowering_date))

  # Ensure flowering time data is numeric
  data <- mutate(data,
                 days_to_flower = as.numeric(days_to_flower),
                 node_num = as.numeric(node_num))

  # Relevent summary data
  rel.sum <- filter(data, sowing_number %in% LoI) %>%
    select(sowing_number, genotype, conditions, days_to_flower, node_num, mutant)

  # When grouping by condition I would like to preserve the sowing number in the
  # nested data so I can then extract it as a new column. This is not necessary
  # when grouping by sowing number as each has a single condition
  if(quo_name(group) == "sowing_number"){
    rel.sum <- rel.sum %>%
      group_by(!!group, conditions, genotype) %>%
      nest(days_to_flower, node_num, mutant)
  }
  else if(quo_name(group) == "conditions"){
    rel.sum <- rel.sum %>%
      group_by(!!group, genotype) %>%
      nest(sowing_number, days_to_flower, node_num, mutant) %>%
      mutate(sowing_numbers = map(data, ~ unique(.x$sowing_number)))
  }
  else if(quo_name(group) == "mutant"){
    rel.sum <- rel.sum %>%
      group_by(!!group, conditions, genotype) %>%
      nest(sowing_number, days_to_flower, node_num) %>%
      mutate(sowing_numbers = map(data, ~ unique(.x$sowing_number)))
  }

  # Here I calculate confidence intervals
  rel.sum <- rel.sum %>%
    mutate(days_to_flower_mean = map_dbl(data, ~ mean(.x$days_to_flower, na.rm = T)),
           node_num_mean = map_dbl(data, ~ mean(.x$node_num, na.rm = T)),
           days_to_flower_moe = map_dbl(data, ~ moe(.x$days_to_flower)),
           node_num_moe = map_dbl(data, ~ moe(.x$node_num)),

           days_to_flower_lower_ci = days_to_flower_mean - days_to_flower_moe,
           days_to_flower_upper_ci = days_to_flower_mean + days_to_flower_moe,
           node_num_lower_ci = node_num_mean - node_num_moe,
           node_num_upper_ci = node_num_mean + node_num_moe) %>%
    select(-days_to_flower_moe, -node_num_moe)

  # This needs to be generalised
  if(WT == "R108"){
    rel.sum <- rel.sum %>%
      mutate(genotype = factor(genotype),
             genotype = fct_relevel(genotype, "R108"))
  }

  return(rel.sum)
}

#' flr.summary.old
#'
#' @param data This is a dataframe of flowering time fata with the following columns; Sowing.Number,
#'  Conditions, Genotype, Plant.Number, Date.Sown, Date.Flowered,
#'  Number.of.Days.to.First.Floral.Bud, Position.of.1st.Floral.Bud,
#'  Node.Number.of.Main.Axis.at.1st.Floral.Bud and Comments.
#'
#' @param LoI These are the Lines of Interest to be examined. Enter these as a vector of sowing
#'  numbers as characters.
#'
#' @return This function summarises a flowering time data frame producing a list of the
#' summarised data and relevant raw data
#' @examples data = Mid_2016
#' LoI = c("P912", "P913") # Lines of interest
#' d <- flr.summary(data, LoI)
#' @export
flr.summary.old <- function(data, LoI){

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
            geom_bar(data = data[["Summary"]], stat = "identity") +
            geom_point(data = data[["Raw"]], aes_string(x = x, y = "Days_to_flower")) +
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
    geom_bar(data = data[["Summary"]], stat = "identity") +
    geom_point(data = data[["Raw"]], aes_string(x = x, y = "Node_num")) +
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

  # Subset list
  dat <- data$Raw

  # Creates the basic plots based on measurment of flowering
  if(meas == "nd"){
    plot <- ggplot(data=dat, aes(Node_num, fill = Mutant)) +
      geom_histogram(data=subset(dat, Genotype!="R108"),aes(y=..count..), binwidth = 1) +
      geom_histogram(data=subset(dat, Genotype=="R108"),aes(x = Node_num, y=..count..), binwidth = 1)
  }
  else if(meas == "fl"){
    plot <- ggplot(data=dat, aes(Days_to_flower, fill = Mutant)) +
      geom_histogram(data=subset(dat, Genotype!="R108"),aes(y=..count..), binwidth = 1) +
      geom_histogram(data=subset(dat, Genotype=="R108"),aes(x = Days_to_flower, y=..count..), binwidth = 1)
  }
  else{ stop("ERROR: Please specify a correct measurment of flowering. See docs for details.") }

  # Generic plot aesthetics
  plot <- plot +
    scale_fill_manual(values = g_colours[c(3,2)]) +
    ggthemes::theme_wsj() +
    theme(legend.position="none",
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.title = element_text(size = 14),
          axis.title.y = element_text(angle = 90),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          strip.background = element_blank(),
          strip.text = element_text(hjust=0, size = 16)) +
    scale_y_continuous(expand=c(0,0), breaks=scales::pretty_breaks()) +
    scale_x_continuous(breaks=scales::pretty_breaks())

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

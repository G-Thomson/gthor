#' read.sds
#'
#' @param path This is a character string to point you to the file
#' @param col These are the columns you want to keep, by default it selects
#' "Well", "Sample",  "Detector", "Reporter" and "Ct". If you want to change
#' it it takes a character vector
#'
#' @return This function imports SDS 2.4 files from a RT-qPCR machine into
#' R as a dataframe. It does this by cutting off the header, selecting the
#' relevant columns and ordering the rows.
#' @examples read.sds(Path)
#' @export
read.sds <- function(path, col = c("Well", "Sample",  "Detector", "Reporter",  "Ct")){

  require(dplyr)

  # Read raw table ignoring the header
  df <- read.delim(Path, skip=10)
  colnames(df)[2:3] <- c("Sample", "Detector")

  # Select the relevant rows and columns
  df <- select(df, Well, Sample, Detector, Reporter, Ct) %>%
    filter(Well %in% as.character(1:384)) %>%
    mutate(Well = as.numeric(as.character(Well)))

  # Add in omitted rows up to your last used cell - I do this so it is easy to add
  # in well coordinates
  miss <- which(!seq_len(max(df$Well, na.rm = TRUE)) %in% df$Well)
  empt_rows <- cbind(miss, data.frame(matrix(ncol = ncol(df)-1, nrow = length(miss))))
  colnames(empt_rows) <- colnames(df)
  df <- bind_rows(df, empt_rows)

  # Sort by Well column
  df <- arrange(df, Well)

  # Modify the Well column to plate co-ordinates
  Pos <- sapply(LETTERS[1:16], paste0, "%d") %>%
    lapply(sprintf,  seq_len(24)) %>%
    unlist()
  df["Well"] <- Pos[1:nrow(df)]

  # Remove empty rows
  df <- df[complete.cases(df), ]

  # Make Cts numeric - there will be a warning about "Undefined" being coerced to
  # NA. I supress this.
  df["Ct"] <- suppressWarnings(as.numeric(as.character(df$Ct)))

  return(df)
}

#' sds.summary
#'
#' @param data This is a dataframe created by the read.sds function
#' @param Experimental_conditions If you have multiple experimental conditions that you want to compare this
#' is a character vector of those conditions. By default it is NULL
#' @param Problems This is a character vector of the wells (specified by coordinate e.g. "A1")
#' which have been identified as problematic by looking at the Ct curves etc. They get
#' removed from further analysis
#'
#' @return This function takes the SDS data and does the munging to get it in a useful
#' form for your experiment, especially if multiple experiments are done on one plate.
#' Essentially it averages out technical replicates for a given experiment and groups
#' biological replicates of the same measurements for downstram analysis.
#' It returns a list with each entry being a dataframe for a gene/condition combination.
#' @examples Experimental_conditions <- c("SD", "LD")
#' Problems <- c("A14", "B1", "B12", "I4", "I5", "I6", "J1")
#' sds.summary(test, Experimental_conditions = Experimental_conditions, Problems = Problems)
#' @export
sds.summary <- function(data, Experimental_conditions  = NULL, Problems = NULL){

  require(dplyr)
  require(stringr)

  # This is where the input data is
  df <- data
  # .. and this is where the output data will end up
  res <- list()

  # Subset data by particular experiment, if required
  if(!is.null(Experimental_conditions)){
    cond <- Experimental_conditions %>%
      paste(collapse = '|')
    df <- filter(df, grepl(cond, Sample))
  }

  # Eliminate problems, if required
  if(!is.null(Problems)){
    df <- filter(df, !Well %in% Problems)
  }

  # Put the genes in this experiment in a vector
  genes <- unique(df$Detector)

  # For each gene...
  for(i in 1:length(genes)){
    # Average out technical replicates
    sub <- filter(df, Detector == genes[i]) %>%
      select(Sample, Ct) %>%
      group_by(Sample)  %>%
      summarise_each(funs(mean(., na.rm = TRUE)))

    # If there are different conditions subset by them
    if(!is.null(Experimental_conditions)){
      sub <- mutate(sub, Condition = str_match(Sample, cond)[1])

      # Extract each condition-gene table in a list
      for(j in 1:length(Experimental_conditions)){
        name <- paste0(genes[i], "_", Experimental_conditions[j])
        sub.cond <- filter(sub, grepl(Experimental_conditions[j], Sample)) %>%
          select(Sample, Ct)
        res[[name]] <- sub.cond
      }
    }

    # Otherwise place each gene table in a list
    else{
      print(as.character(genes[i]))
      res[[as.character(genes[i])]] <- sub
    }
  }

  # Return the list with epty entries eliminated
  res <- res[which(lapply(res, nrow) != 0)]
  return(res)
}

#' sds.test
#'
#' @param data This is a list from the sds.summary function
#' @param Reference_gene The name of your reference gene that you are measuring relative expression against
#' as a character
#' @param Measured_genes A vector of the genes you are measuing as characters
#' @param Experimental_conditions If you have multiple experimental conditions that you want to compare this
#' is a character vector of those conditions. By default it is NULL
#'
#' @return This is the function that produces the results. Essentially it uses the ddCT
#' method to produce a list of dataframes, one for each gene. Each dataframe contains the
#' mean, sd and se for each sample. There is also a "Gene and Condition" tag to make plotting easy.
#' @examples Measured_genes <- c("GOI_1", "GOI_2")
#' sds.test(test2, "PDF2", Measured_genes, Experimental_conditions = Experimental_conditions)
#' @export
sds.test <- function(data, Reference_gene, Measured_genes, Experimental_conditions = NULL, Replicate_id = "_"){

  require(dplyr)
  require(tidyr)
  require(stringr)

  # For the output data
  res <- list()

  # What is in the data
  data.names <- names(data)

  # Analyse one gene at a time
  for(g in 1:length(Measured_genes)){

    # This is where all the gene's data will go
    gene.data <- data.frame()

    # If there are different experimental conditions
    if(!is.null(Experimental_conditions)){

      # Assess one condition at a time and make a table with the corresponding reference gene data
      for(e in 1:length(Experimental_conditions)){
        rel.ref <-  intersect(grep(Reference_gene, data.names, value=TRUE), grep(Experimental_conditions[e], data.names, value=TRUE))
        rel.data <- intersect(grep(Measured_genes[g], data.names, value=TRUE), grep(Experimental_conditions[e], data.names, value=TRUE))
        gene.data <- bind_rows(gene.data, full_join(data[[rel.ref]], data[[rel.data]], by = "Sample"))
      }
    }

    # Otherwise just pick out the relevant gene data
    else{
      rel.ref <-  grep(Reference_gene, data.names, value=TRUE)
      rel.data <- grep(Measured_genes[g], data.names, value=TRUE)
      gene.data <- bind_rows(gene.data, full_join(data[[rel.ref]], data[[rel.data]], by = "Sample"))
   # grep("P951", data.names, value=TRUE)
    }

    # Once a genes data for all conditions (if required) is collected do the analysis, this allows relative expression to be observed
    # over multiple conditions
    gene.res <- mutate(gene.data, dCt = Ct.y - Ct.x) %>%
      select(Sample, dCt) %>%
      mutate(ddCt = dCt-min(dCt, na.rm = TRUE)) %>%
      mutate(comp_exp = 2^-ddCt) %>%
      separate(Sample, c("Sample", "Bio_Rep"), Replicate_id) %>%
      select(Sample, comp_exp) %>%
      group_by(Sample) %>%
      summarise_each(funs(mean = mean(., na.rm = TRUE), sd = sd, len = length)) %>%
      mutate(se = sd/sqrt(len)) %>%
      select(Sample, mean, sd, se)

    # If there are different experimental conditions
    if(!is.null(Experimental_conditions)){

      # Add in a extra descriptive column which makes plotting the data easy
      cond <- Experimental_conditions %>%
        paste(collapse = '|')
      cond <- paste0(Measured_genes[g], "_in_", str_extract(gene.res$Sample, cond))
      gene.res["Gene_and_Condition"] <- cond
    }

    # Save the gene's data to a final output list
    res[[Measured_genes[g]]] <- gene.res
  }

  # Return the list
  return(res)
}

#' read.eds
#'
#' @param path This is a character string to point you to the file
#' @param in_file This is the type of input file. By default it is "Results" for Ct results
#' but you can also set it to "Melt" for melting curve data.
#'
#' @return This function imports csv files from a eds file produced by the Quant 5
#' RT-qPCR machine into R as a dataframe. It does this by cutting off the header, selecting the
#' relevant columns and ordering the rows.
#' @examples read.eds(Path)
#' @examples read.eds(Path, in_file == "Melt")
#' @export
read.eds <- function(Path, in_file = "Results"){

  require(dplyr)

  # Read raw table ignoring the header
  df <- read.delim(Path, skip=17,sep=",")

  # Select the relevant rows and columns
  if(in_file == "Results"){
    df <- select(df, Well, Sample, Target, Ct)
  }
  else if(in_file == "Melt"){
    df <- select(df, Well, Sample, Target, Temperature, Derivative)
  }

  # Make Cts numeric - there will be a warning about "Undefined" being coerced to
  # NA. I supress this.
  if(in_file == "Results"){
    df["Ct"] <- suppressWarnings(as.numeric(as.character(df$Ct)))
  }
  return(df)
}

#' eds.summary
#'
#' @param data This is a dataframe created by the read.sds function
#' @param Experimental_conditions If you have multiple experimental conditions that you want to compare this
#' is a character vector of those conditions. By default it is NULL
#' @param Problems This is a character vector of the wells (specified by coordinate e.g. "A1")
#' which have been identified as problematic by looking at the Ct curves etc. They get
#' removed from further analysis
#'
#' @return This function takes the eds data and does the munging to get it in a useful
#' form for your experiment, especially if multiple experiments are done on one plate.
#' Essentially it averages out technical replicates for a given experiment and groups
#' biological replicates of the same measurements for downstram analysis.
#' It returns a list with each entry being a dataframe for a gene/condition combination.
#' @examples Experimental_conditions <- c("SD", "LD")
#' Problems <- c("A14", "B1", "B12", "I4", "I5", "I6", "J1")
#' eds.summary(test, Experimental_conditions = Experimental_conditions, Problems = Problems)
#' @export
eds.summary <- function(data, Experimental_conditions  = NULL, Problems = NULL){

  require(tidyverse)

  # This is where the input data is
  df <- data
  # .. and this is where the output data will end up
  res <- list()

  # Subset data by particular experiment, if required
  if(!is.null(Experimental_conditions)){
    cond <- Experimental_conditions %>%
      paste(collapse = '|')
    df <- filter(df, grepl(cond, Sample))
  }

  # Eliminate problems, if required
  if(!is.null(Problems)){
    df <- filter(df, !Well %in% Problems)
  }

  # Put the genes in this experiment in a vector
  genes <- unique(df$Target)

  # For each gene...
  for(i in 1:length(genes)){
    # Average out technical replicates
    sub <- filter(df, Target == genes[i]) %>%
      select(Sample, Ct) %>%
      group_by(Sample)  %>%
      summarise_each(funs(mean(., na.rm = TRUE)))

    # If there are different conditions subset by them
    if(!is.null(Experimental_conditions)){
      sub <- mutate(sub, Condition = str_match(Sample, cond)[1])

      # Extract each condition-gene table in a list
      for(j in 1:length(Experimental_conditions)){
        name <- paste0(genes[i], "_", Experimental_conditions[j])
        sub.cond <- filter(sub, grepl(Experimental_conditions[j], Sample)) %>%
          select(Sample, Ct)
        res[[name]] <- sub.cond
      }
    }

    # Otherwise place each gene table in a list
    else{
      print(as.character(genes[i]))
      res[[as.character(genes[i])]] <- sub
    }
  }

  # Return the list with epty entries eliminated
  res <- res[which(lapply(res, nrow) != 0)]
  return(res)
}

#' eds.test
#'
#' @param data This is a list from the eds.summary function
#' @param Reference_gene The name of your reference gene that you are measuring relative expression against
#' as a character
#' @param Measured_genes A vector of the genes you are measuing as characters
#' @param Experimental_conditions If you have multiple experimental conditions that you want to compare this
#' is a character vector of those conditions. By default it is NULL
#'
#' @return This is the function that produces the results. Essentially it uses the ddCT
#' method to produce a list of dataframes, one for each gene. Each dataframe contains the
#' mean, sd and se for each sample. There is also a "Gene and Condition" tag to make plotting easy.
#' @examples Measured_genes <- c("GOI_1", "GOI_2")
#' eds.test(test2, "PDF2", Measured_genes, Experimental_conditions = Experimental_conditions)
#' @export
eds.test <- function(data, Reference_gene, Measured_genes, Experimental_conditions = NULL, Replicate_id = "_"){

  require(tidyverse)

  # For the output data
  res <- list()

  # What is in the data
  data.names <- names(data)

  # Analyse one gene at a time
  for(g in 1:length(Measured_genes)){

    # This is where all the gene's data will go
    gene.data <- data.frame()

    # If there are different experimental conditions
    if(!is.null(Experimental_conditions)){

      # Assess one condition at a time and make a table with the corresponding reference gene data
      for(e in 1:length(Experimental_conditions)){
        rel.ref <-  intersect(grep(Reference_gene, data.names, value=TRUE), grep(Experimental_conditions[e], data.names, value=TRUE))
        rel.data <- intersect(grep(Measured_genes[g], data.names, value=TRUE), grep(Experimental_conditions[e], data.names, value=TRUE))
        gene.data <- bind_rows(gene.data, inner_join(data[[rel.ref]], data[[rel.data]], by = "Sample"))
      }
    }

    # Otherwise just pick out the relevant gene data
    else{
      rel.ref <-  grep(Reference_gene, data.names, value=TRUE)
      rel.data <- grep(Measured_genes[g], data.names, value=TRUE)
      gene.data <- bind_rows(gene.data, inner_join(data[[rel.ref]], data[[rel.data]], by = "Sample"))
    }

#     return(Measured_genes[g])}
# }
# eds.test(data_res_New_FTa1, "MtPDF2", "New_FTa1", Replicate_id = "_B")

    # Once a genes data for all conditions (if required) is collected do the analysis, this allows relative expression to be observed
    # over multiple conditions
    gene.res <- mutate(gene.data, dCt = Ct.y - Ct.x) %>%
      select(Sample, dCt) %>%
      mutate(ddCt = dCt-min(dCt, na.rm = TRUE)) %>%
      mutate(comp_exp = 2^-ddCt) %>%
      separate(Sample, c("Sample", "Bio_Rep"), Replicate_id) %>%
      select(Sample, comp_exp) %>%
      group_by(Sample) %>%
      summarise_each(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE), len = length)) %>%
      mutate(se = sd/sqrt(len)) %>%
      select(Sample, mean, sd, se)

    # If there are different experimental conditions
    if(!is.null(Experimental_conditions)){

      # Add in a extra descriptive column which makes plotting the data easy
      cond <- Experimental_conditions %>%
        paste(collapse = '|')
      cond <- paste0(Measured_genes[g], "_in_", str_extract(gene.res$Sample, cond))
      gene.res["Gene_and_Condition"] <- cond
    }

    # Save the gene's data to a final output list
    res[[Measured_genes[g]]] <- gene.res
  }

  # Return the list
  return(res)
}

#' rt.summary
#'
#' @export
rt.summary <- function(data){

  require(dplyr)

  data %>%
    select(date, sample, genotype, condition, bio_replicate, target, Ct) %>%
    group_by(date, target, sample, genotype, condition, bio_replicate) %>%
    summarise(mean_ct = mean(Ct), sd_div_mean = sd(Ct)/mean(Ct)) %>%
    mutate(sd_check = sd_div_mean < 0.015)
}

#' rt.test
#'
#' @export
rt.test <- function(data, reference = MtPDF2){

  require(dplyr)
  require(tidyr)
  require(stringr)
  require(purrr)

  # Make reference a quosure
  reference <- enquo(reference)

  # Check data format
  if(sum(colnames(data) == c("date",
                             "target",
                             "sample",
                             "genotype",
                             "condition",
                             "bio_replicate",
                             "mean_ct",
                             "sd_div_mean",
                             "sd_check")) != 9){

    stop("ERROR: Input must be a data frame from rt.summary().")}

  # Check the reference gene is there
  if(!quo_name(reference) %in% (data %>% pull(target))){

    stop("ERROR: Your reference gene isn't in the dataset.")
  }

  # Perform the ddCt method
  data %>%
    ungroup() %>%
    select(-sd_div_mean, -sd_check) %>%
    spread(target, mean_ct) %>%
    mutate(sample = as.character(sample)) %>% # for future sanity
    mutate_at(vars(-one_of(c("date","sample", "genotype", "condition", "bio_replicate", "MtPDF2"))),
              funs(dCt = . - MtPDF2)) %>%
    mutate_at(vars(contains("dCt")),
              funs(ddCt = . - min(.,  na.rm = TRUE))) %>%
    mutate_at(vars(contains("ddCt")),
              funs(comp_exp = 2 ^ - .)) %>%
    filter(!is.na(ddCt)) %>% # problem with multiple genes at once.
    nest(-sample, -genotype, -condition) %>%



    mutate(mean = map(data, ~ summarise_at(.x, vars(contains("comp_exp")),
                                           funs(mean_expr = mean(., na.rm = T) ))), # T NOT TRUE for some reason
           se   = map(data, ~ summarise_at(.x, vars(contains("comp_exp")),
                                           funs(se = sd(., na.rm = T)/sqrt(sum(!is.na(.))) ))), # T NOT TRUE for some reason

           # Shorten names
           data = map(data, ~ set_colnames(.x, str_replace(names(.), "_dCt_", "_"))),
           data = map(data, ~ set_colnames(.x, str_replace(names(.), "_ddCt_", "_"))),
           mean = map(mean, ~ set_colnames(.x, str_replace(names(.), "_dCt_ddCt_comp_exp_", "_"))),
           se = map(se, ~ set_colnames(.x, str_replace(names(.), "_dCt_ddCt_comp_exp_", "_")))
    ) %>%
    unnest(mean, se) %>%
    select(1:4, order(colnames(.)))
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qrt_pcr.R
\name{sds.summary}
\alias{sds.summary}
\title{sds.summary}
\usage{
sds.summary(data, Experimental_conditions = NULL, Problems = NULL)
}
\arguments{
\item{data}{This is a dataframe created by the read.sds function}

\item{Experimental_conditions}{If you have multiple experimental conditions that you want to compare this
is a character vector of those conditions. By default it is NULL}

\item{Problems}{This is a character vector of the wells (specified by coordinate e.g. "A1")
which have been identified as problematic by looking at the Ct curves etc. They get
removed from further analysis}
}
\value{
This function takes the SDS data and does the munging to get it in a useful
form for your experiment, especially if multiple experiments are done on one plate.
Essentially it averages out technical replicates for a given experiment and groups
biological replicates of the same measurements for downstram analysis.
It returns a list with each entry being a dataframe for a gene/condition combination.
}
\description{
sds.summary
}
\examples{
Experimental_conditions <- c("SD", "LD")
Problems <- c("A14", "B1", "B12", "I4", "I5", "I6", "J1")
sds.summary(test, Experimental_conditions = Experimental_conditions, Problems = Problems)
}

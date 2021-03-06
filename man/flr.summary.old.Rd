% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flowering_time.R
\name{flr.summary.old}
\alias{flr.summary.old}
\title{flr.summary.old}
\usage{
flr.summary.old(data, LoI)
}
\arguments{
\item{data}{This is a dataframe of flowering time fata with the following columns; Sowing.Number,
Conditions, Genotype, Plant.Number, Date.Sown, Date.Flowered,
Number.of.Days.to.First.Floral.Bud, Position.of.1st.Floral.Bud,
Node.Number.of.Main.Axis.at.1st.Floral.Bud and Comments.}

\item{LoI}{These are the Lines of Interest to be examined. Enter these as a vector of sowing
numbers as characters.}
}
\value{
This function summarises a flowering time data frame producing a list of the
summarised data and relevant raw data
}
\description{
flr.summary.old
}
\examples{
data = Mid_2016
LoI = c("P912", "P913") # Lines of interest
d <- flr.summary(data, LoI)
}

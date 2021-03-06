% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nof1.summary.R
\name{raw_table}
\alias{raw_table}
\title{Summary data table for nof1}
\usage{
raw_table(nof1)
}
\arguments{
\item{nof1}{nof1 object created using nof1.data}
}
\description{
Summary data table for nof1
}
\examples{
Y <- laughter$Y
Treat <- laughter$Treat
nof1 <- nof1.data(Y, Treat, ncat = 11, baseline = "Usual Routine", response = "ordinal")
raw_table(nof1)
}

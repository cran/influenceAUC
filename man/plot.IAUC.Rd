% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporter.R
\name{plot.IAUC}
\alias{plot.IAUC}
\title{Visualize IAUC result}
\usage{
\method{plot}{IAUC}(x, ...)
}
\arguments{
\item{x}{An IAUC class object for `plot` method}

\item{...}{Not used directly}
}
\description{
Visualize IAUC output sequentially
}
\examples{
library(ROCR)
data("ROCR.simple")
Ioutput <- IAUC(ROCR.simple$predictions, ROCR.simple$labels)
plot(Ioutput)
}
\seealso{
\link{IAUC}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporter.R
\name{plot.LAUC}
\alias{plot.LAUC}
\title{Visualize LAUC results}
\usage{
\method{plot}{LAUC}(x, ...)
}
\arguments{
\item{x}{An LAUC class object for `plot` method}

\item{...}{Not used directly}
}
\description{
Visualize LAUC output sequentially
}
\examples{
library(ROCR)
data("ROCR.simple")
Loutput <- LAUC(ROCR.simple$predictions, ROCR.simple$labels)
plot(Loutput)
}
\seealso{
\link{LAUC}
}

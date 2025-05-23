% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stars.R
\name{flr2stars}
\alias{flr2stars}
\title{flr2stars()}
\usage{
flr2stars(object, uncertainty = NULL, quantiles = c(0.05, 0.95))
}
\arguments{
\item{object}{of class FLStockR (MLE)}

\item{uncertainty}{of class FLStock with iters}

\item{quantities}{default is 90CIs as c(0.05,0.95)}
}
\value{
STARS list with $timeseris and $refpts
}
\description{
flr2stars()
}

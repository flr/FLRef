% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ices.R
\name{ss2ices}
\alias{ss2ices}
\title{ss2ices()}
\usage{
ss2ices(mvln, quantiles = c(0.05, 0.95))
}
\arguments{
\item{mvln}{output of ssmvln()}

\item{quantiles}{default is 95CIs as c(0.025,0.975)}

\item{output}{choice c("iters","mle")[1]}

\item{Fmsy}{if specified the ratio F/Fmsy is calculated (required for ensembles)}

\item{Btrigger}{if specified the ratio SSB/Btrigger is calculated (required for ensembles)}
}
\value{
ICES list with $timeseris and $refpts
}
\description{
ss2ices()
}

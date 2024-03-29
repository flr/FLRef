% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotDyns.R
\name{plotbioage}
\alias{plotbioage}
\title{plotbioage()
Plots stock N_a, W_a, M_a and Mat_a by year}
\usage{
plotbioage(stk, ncol = 2)
}
\arguments{
\item{stk}{stock object class FLStock}

\item{ncol}{number of columns in multiplot}
}
\value{
ggplot
}
\description{
plotbioage()
Plots stock N_a, W_a, M_a and Mat_a by year
}
\examples{
data(ple4)
plotbioage(ple4)
}

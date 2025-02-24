% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesisTestingTools.R
\name{twoSidedZ}
\alias{twoSidedZ}
\title{Two Sided Z function}
\usage{
twoSidedZ(z)
}
\arguments{
\item{z}{vector of quantiles.}
}
\description{
This function converts a test statistic x into the area under the z-density for values >|x|
}
\details{
Refer to the \code{pnorm} function for more details
}
\examples{
twoSidedZ(2)

}

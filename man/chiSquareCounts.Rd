% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesisTestingTools.R
\name{chiSquareCounts}
\alias{chiSquareCounts}
\title{Chi Square Test}
\usage{
chiSquareCounts(data, var1, var2)
}
\arguments{
\item{tib}{data supplied as a tibble}

\item{x, y}{variable as arguments}
}
\value{
It returns the test statistics, the p-value and the degree of freedom
}
\description{
Calculates a chi-square test for count data supplied as a tibble
}
\details{
You may supply variable names as strings,
or with syntax similar to janitor::tabyl using the curly-curly operator in R ({{}}).
}

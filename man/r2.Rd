% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesisTestingTools.R
\name{r2}
\alias{r2}
\title{R Squared Function}
\usage{
r2(y_pred, y_true)
}
\arguments{
\item{y_pred}{array-like of shape (n_samples,) Estimated target values}

\item{y_true}{array-like of shape (n_samples,) Ground truth (correct) target values.}
}
\value{
outputs the r-square value
}
\description{
Calculates the r-squared value
}
\examples{
y_true = c(4, -1.5, 5, 9)
y_pred = c(3.5, 1.0, 5, 10)
r2(y_pred, y_true)

}

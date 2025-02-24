% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classificationMetricsTools.R
\name{specificity}
\alias{specificity}
\title{Specificity Function}
\usage{
specificity(y_pred, y_true)
}
\arguments{
\item{y_pred}{array-like of shape (n_samples,) Estimated targets as returned by a classifier.}

\item{y_true}{array-like of shape (n_samples,) Ground truth (correct) target values.}
}
\value{
Outputs the specificity
}
\description{
Calculates the specificity  given the true and predicted values
}
\details{
Specificity is defined as True Negatives / All Negatives
}
\examples{
Simulated y_true and y_predicted values
y_true = c(1,1,1,0,0,1,0,1,1,0)
y_pred = c(1,1,1,1,0,0,1,0,1,1)
specificity(y_pred, y_true)

}

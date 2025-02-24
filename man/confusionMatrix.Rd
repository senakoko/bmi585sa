% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classificationMetricsTools.R
\name{confusionMatrix}
\alias{confusionMatrix}
\title{Confusion Matrix Function}
\usage{
confusionMatrix(y_pred, y_true)
}
\arguments{
\item{y_pred}{array-like of shape (n_samples,) Estimated targets as returned by a classifier.}

\item{y_true}{array-like of shape (n_samples,) Ground truth (correct) target values.}
}
\value{
Outputs the confusion matrix
}
\description{
Compute confusion matrix to evaluate the accuracy of a classification.
It is used to calculate other metrics like precision. sensitivity and f1_score.
}
\examples{
Simulated y_true and y_predicted values
y_true = c(1,1,1,0,0,1,0,1,1,0)
y_pred = c(1,1,1,1,0,0,1,0,1,1)
confusionMatrix(y_pred,y_true)

}

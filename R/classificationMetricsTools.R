#' Confusion Matrix Function
#' 
#' Compute confusion matrix to evaluate the accuracy of a classification.
#' It is used to calculate other metrics like precision. sensitivity and f1_score.
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the confusion matrix
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' confusionMatrix(y_pred,y_true)
#' 
confusionMatrix = function(y_pred, y_true){
  
  # Convert y_true and y_pred to integer 
  y_true = as.integer(y_true)
  y_pred = as.integer(y_pred) 
  # find unique values
  uniq = unique(y_true)
  uniq = uniq[order(uniq)]
  n = length(uniq)
  c_matrix = matrix(0,nrow = n, ncol = n)
  # loop to generate confusion matrix
  for (i in 1:length(y_true)) {
    if (y_true[i] == y_pred[i]) {
      v = which(y_true[i] == uniq)
      c_matrix[v,v] = c_matrix[v,v] + 1
    }
    else{
      v = which(y_true[i] == uniq)
      u = which(y_pred[i] == uniq)
      c_matrix[v,u] = c_matrix[v,u] + 1
    }
  }
  return(c_matrix)
}

#' Sensitivity Function
#' 
#' Compute the sensitivity given the true and predicted values
#' 
#' @details 
#' Sensitivity is defined as 
#' True Positives / All Positives = True Positives / ( True Positives + False Negatives). 
#' Sensitivity is also referred to as recall.
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the sensitivity
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' sensitivity(y_pred,y_true)
#' 
sensitivity = function(y_pred,y_true){
  
  cmatrix = confusionMatrix(y_pred, y_true)
  sens_rate = cmatrix[2,2]/(cmatrix[2,2] + cmatrix[2,1] )
  return(sens_rate)
}

#' Specificity Function
#' 
#' Calculates the specificity  given the true and predicted values
#' 
#' @details
#' Specificity is defined as True Negatives / All Negatives
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the specificity
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' specificity(y_pred, y_true)
#' 
specificity = function(y_pred, y_true){
  
  cmatrix = confusionMatrix(y_pred, y_true)
  spec_rate = cmatrix[1,1]/(cmatrix[1,1] + cmatrix[1,2])
  return(spec_rate)
}

#' Accuracy Function
#' 
#' Calculates the accuracy  given the true and predicted values
#' 
#' @details 
#' Accuracy is defined as (TP +  TN) / (All Positives + All Negatives).
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the Acccuracy
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' accuracy(y_pred, y_true)
#' 
accuracy = function(y_pred, y_true){
  
  cmatrix = confusionMatrix(y_pred, y_true)
  accu_rate = (cmatrix[2,2] + cmatrix[1,1])/((cmatrix[2,1]+cmatrix[2,2]) + (cmatrix[1,1]+ cmatrix[1,2]))
  return(accu_rate)
}

#' Positive Predictive Value or Precision Function
#' 
#' Calculates the Positive Predictive Value or precision given the true and predicted values
#' 
#' @details 
#' Positive Predictive Value, or PPV, is defined as True Positives / ( True Positives + False Positives). 
#' PPV is also referred to as precision.
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the precision
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' ppv(y_pred, y_true)
#' 
ppv = function(y_pred, y_true){
  
  cmatrix = confusionMatrix(y_pred, y_true)
  ppv_rate = cmatrix[2,2] / (cmatrix[2,2] + cmatrix[1,2])
  return(ppv_rate)
}

#' F1 Score Function
#' 
#' Calculates the F1-Score given the true and predicted values
#' 
#' @details 
#' F1 score is defined as 2 x (precision x recall) / (precision + recall)
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated targets as returned by a classifier.
#'
#' @return 
#' Outputs the F1-score
#' 
#' @export
#' @examples 
#' Simulated y_true and y_predicted values
#' y_true = c(1,1,1,0,0,1,0,1,1,0)
#' y_pred = c(1,1,1,1,0,0,1,0,1,1)
#' f1(y_pred, y_true)
#' 
f1 = function(y_pred, y_true){
  
  precision = ppv(y_pred, y_true)
  recall = sensitivity(y_pred,y_true)
  f1_score = 2 * ((precision * recall) / (precision + recall))
  return(f1_score)
}




#' Unscale Matrix-like Objects
#' 
#' This function takes a numerical object from scale and reverses the centering/scaling. 
#' It assumes the input is the output of the scale function in R.
#' 
#' @param X  vector or matrix
#' 
#' @return 
#' An unscaled vector or matrix
#' 
#' @details 
#' Refer to the scale function in R for more details about it. 
#' 
#' @export
#' @examples  
#' x = c(2,4,3,5)
#' X = matrix(x, nrow = 2)
#' X_scaled = scale(X)
#' unscale(X_scaled)
#' 
unscale = function(X){
  center = attr(X,"scaled:center")
  scale=attr(X,"scaled:scale")
  X_unscale = sweep(X,2,scale,FUN =  "*")
  X_unscale = sweep(X_unscale, 2, center, FUN =   "+")
  return(X_unscale)
}


#' Principal Component Approximation
#' 
#' This function returns an approximation to the data x based on the number of Principal 
#' Components (PCs). 
#' 
#' @param X    tibble or matrix
#' @param npc  the number of pca components
#' 
#' @details 
#' The closer the number of principal components to the number of features in the original data, 
#' the better approximated the input value.
#' 
#' @return 
#' An approximated tibble of the original data
#' 
#' @export
#' @examples 
#' npc = 6
#' pcApprox(mtcars[,1:npc],npc)
#' 
pcApprox = function(X, npc){
  X_scaled = scale(X)
  X2D = princomp(X_scaled, scores = TRUE)
  scores = X2D$scores[,1:npc]
  loading = X2D$loadings[,1:npc]
  X_orig = as.matrix(scores) %*% as.matrix(loading)
  center = attr(X_scaled,"scaled:center")
  scale=attr(X_scaled,"scaled:scale")
  X_unscale = sweep(X_orig,2,scale,FUN =  "*")
  X_unscale = sweep(X_orig, 2, center, FUN =   "+")
  return(X_unscale)
}

#' Lollipop Plot
#' 
#' This function creates a lollipop plot of the principal 
#' component loadings of the (potentially unscaled/uncentered) data x
#' 
#' @param X   tibble or matrix
#' 
#' @details 
#' The function scales the data before performing the principal component analysis.
#' 
#' @export
#' @returns 
#' A lollipop plot of the principal component loadings
#'
#' @examples 
#' pcLollipop(mtcars[,1:6])
#'
pcLollipop = function(X){
  X_scaled = scale(X)
  X2D = princomp(X_scaled, scores = TRUE)
  scores = X2D$scores
  loading = X2D$loadings[,1:dim(X)[2]]
  loading = as.matrix(loading)
  feat_names = rownames(loading)
  loading_new = as_tibble(cbind(feat_names, round(loading,1)))
  loading_long = loading_new |> pivot_longer(-feat_names, names_to = "features", 
                                             values_to = "comp_values")
  loading_long$comp_values = as.numeric(loading_long$comp_values)
  ll_plot = loading_long |> ggplot(aes(x=feat_names, y=comp_values)) +
    geom_point() + geom_segment(aes(x=feat_names, xend=feat_names, 
                                    yend=comp_values,y=0, color=features)) +
    facet_wrap(~features, scales = "fixed", as.table = FALSE)
  print(ll_plot)
}
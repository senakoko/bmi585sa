#' BoxMuller Function
#'
#' The Box-Muller transformation takes two samples from the uniform distribution 
#' on the interval [0, 1] and maps them to two standard, normally distributed samples:
#' 𝑋1=−2log(𝑈1)‾‾‾‾‾‾‾‾‾‾√cos(2𝜋𝑈2)
#' 𝑋2=−2log(𝑈1)‾‾‾‾‾‾‾‾‾‾√sin(2𝜋𝑈2)
#' 
#' @param n the number of samples to create normal distribution
#' 
#' @return 
#' a normal distribution
#' 
#' @export
#' @examples 
#' n = 500
#' BM_data = boxMuller(n)
#' find means of boxMuller
#' round(mean(BM_data[,1]),2)
#' find variance of boxMuller  
#' round(var(BM_data[,1]),2) 
#' 
boxMuller = function(n=500){
  bm_mat = matrix(nrow = n, ncol = 2)
  for (i in 1:n) {
    u1 =  runif(1)
    u2 =  runif(1)
    exp1 =  cos(2 * pi * u2)
    exp2 =  sin(2 * pi * u2)
    nom1 =  sqrt(-2 * log(u1))
    nom2 =  sqrt(-2 * log(u1))
    X1 =  nom1 * exp1
    X2 =  nom2 * exp2
    bm_mat[i,1] = X1
    bm_mat[i,2] = X2
  }
  return(bm_mat)
}

#' Two Sided T function
#'
#' This function converts a test statistic x into the area under the t-density for values >|x|
#' 
#' @details 
#' Refer to the `pt` function for more details 
#' 
#' @param t  vector of quantiles.
#' 
#' @param df degrees of freedom (> 0, maybe non-integer). df = Inf is allowed.
#' 
#' @export
#' @examples 
#' twoSidedT(t=5, df=2)
#' 
twoSidedT = function(t,df){
  two_value = 2 * pt(q = t, df, lower.tail = FALSE)
  two_value = round(two_value*100,2)
  return(two_value)
}

#' Two Sided Z function
#'
#' This function converts a test statistic x into the area under the z-density for values >|x|
#' 
#' @details 
#' Refer to the `pnorm` function for more details 
#' 
#' @param z  vector of quantiles.
#' 
#' @export
#' @examples 
#' twoSidedZ(2)
#' 
twoSidedZ = function(z){
  two_value = 2 * pnorm(q = z, lower.tail=FALSE)
  two_value = round(two_value*100,2)
  return(two_value)
}

#' Simulate data 
#' 
#' This function is used to generate data for the effect size function. 
#' It creates a simulated data with treatment and control groups
#' 
#' @param N  the number of samples 
#' 
#' @return 
#' a simulated with treatment and control groups
#' 
#' @export
#' @examples 
#' simulated_cancer = simulateTrial(50)
#' 
simulateTrial = function(N=50){
  b_not = 4
  b_treat = -2
  treat = rep(1,N)
  control = rep(0,N)
  treat_cont = c(treat,control)
  epi = rnorm(length(treat_cont),0,1)
  tm = b_not + b_treat*treat_cont + epi
  result = tibble(treat_grp = treat_cont, tm_size = tm)
  return(result)
}
#' Effect Size Function
#' 
#' A function to calculate Cohen’s d effect size.
#' 
#' @details 
#' This is the difference in means between two groups 
#' divided by some measure of the standard deviation.
#' 
#' @param t - vector of data points
#' 
#' @param g - groups for which each value in parameter t belongs
#' 
#' @return 
#' outputs the Cohen’s d effect size
#' 
#' @export
#' @examples 
#' simulated_cancer = simulateTrial(50)
#' effectSize(simulated_cancer, simulated_cancer[,1]) 
#' 
effectSize = function(t,g){
  g_v = unique(g)
  tm_mean = rep(NA,dim(g_v)[1])
  tm_sd = rep(NA,dim(g_v)[1])
  tm_n = rep(NA,dim(g_v)[1])
  for (i in 1:dim(g_v)[1]) {
    g_value = g_v[[i,1]]
    tm_size = t$tm_size[t$treat_grp == g_value]
    tm_mean[i] = mean(tm_size, na.rm=T)
    tm_sd[i] = sd(tm_size, na.rm = T)
    tm_n[i] = length(tm_size)
  }
  tm_diff = abs(diff(tm_mean))
  sd_pooled = sqrt(((tm_n[1]-1)*tm_sd[1]**2)+ ((tm_n[2]-1)*tm_sd[2]**2))
  cohen_d = tm_diff / sd_pooled
  return(cohen_d)
}


#' Welch T-test function
#' 
#' Perform a welch t-test
#' 
#' @details 
#' The output from this function should be closer to the value of the t.test function in R.
#' 
#' @param x  vector of values
#' 
#' @param y  vector of values. (Should be the same length as x)
#' 
#' @return 
#' It returns the test statistics, the p-value and the degree of freedom
#' 
#' @export
#' @examples 
#' x = rnorm(10,0,1)
#' y = rnorm(10,1,1)
#' welchT(x,y)
#' 
welchT = function(x,y){
  mean_x = mean(x, na.rm=T)
  sd_x = sd(x, na.rm = T)
  mean_y = mean(y, na.rm=T)
  sd_y = sd(y, na.rm = T)
  nx = length(x)
  ny = length(y)
  
  mean_diff = mean_x - mean_y
  SD = sqrt((sd_x/sqrt(length(x)))^2 + (sd_y/sqrt(length(y)))^2)
  t_value = mean_diff / SD
  deg_n = (sd_x^2/length(x) + sd_y^2/length(y))^2
  deg_d = sd_x^4/(length(x)^2*length(x)-1) + sd_y^4/(length(y)^2*length(y)-1)
  deg = deg_n/deg_d
  
  p_value = 2*pt(t_value,df=deg, lower.tail = FALSE)
  
  result = data.frame(t_value = t_value,p_value=p_value, df = deg)
  return(result)
}

#' Minimum Sample Size Function
#' 
#' This function finds the minimum (per group) sample size when given the effect size
#'
#' @param d  the effect size
#' 
#' @return 
#' the minimum sample size
#' 
#' @export
#' @examples 
#' minimumN(0.6)
#' 
minimumN = function(d=0.8){
  result = power.t.test(power = 0.8, delta = d)
  return(ceiling(result$n))
}


#' Chi Square Test
#'
#' calculates a chi-square test for count data supplied as a contingency table
#' 
#' @param tib data supplied as a contingency table
#' 
#' @export
#' @return 
#' It returns the test statistics, the p-value and the degree of freedom
#' 
chiSquareCont = function(tib){
  dim_d = dim(a)
  expected = matrix(nrow = dim_d[1], ncol = dim_d[2])
  for (r in 1:dim_d[1]) {
    for (c in 1:dim_d[2]) {
      expected[r,c] = (sum(a[r,])*sum(a[,c]))/sum(a)
    }
  }
  chi_values = matrix(nrow = dim_d[1], ncol = dim_d[2])
  for (r in 1:dim_d[1]) {
    for (c in 1:dim_d[2]) {
      chi_values[r,c] = (a[r,c] - expected[r,c])^2/expected[r,c]
    }
  }
  df = (dim_d[1]-1)*(dim_d[2]-1)
  t_stats = sum(chi_values)
  p_value = 1-pchisq(t_stats, df=df)
  return(c(t_stats,df,p_value))
}

#' Chi Square Test
#'
#' Calculates a chi-square test for count data supplied as a tibble
#' 
#' @param tib data supplied as a tibble
#' 
#' @param x,y variable as arguments
#' 
#' @details 
#' You may supply variable names as strings, 
#' or with syntax similar to janitor::tabyl using the curly-curly operator in R ({{}}).
#' 
#' @export
#' @return 
#' It returns the test statistics, the p-value and the degree of freedom
#' 
chiSquareCounts = function(data, var1, var2){
a = janitor::tabyl(data, {{var1}}, {{var2}})
dim_d = dim(a)
r_n = dim_d[1]
c_n = dim_d[2] 
b = a[1:r_n,2:c_n]
results = chiSquareCont(b)
return(results)
}

#' Post Hoc Power Analysis
#'
#' This function simulates post-hoc power
#' 
#' @param effect_size the effect size
#' 
#' @param n1,n2 the number of samples per group
#' 
#' @export
#' @return 
#' outputs the post hoc power
#' 
#' @examples 
#' N = 1000
#' data_power = rep(NA,N)
#' for (i in 1:N) {data_power[i] = postHocPower(i)}
#' 
#' or
#' res = replicate(100,postHocPower(effect_size=1, n1=30, n2=30))
#' mean(res)
#' 
postHocPower <- function(effect_size=1, n1=30, n2=30) {
  x1 <- rnorm(n1, 0, 1)
  x2 <- rnorm(n2, effect_size, 1)
  
  t <- t.test(x1, x2, var.equal=TRUE)  # run t-test on generated data
  stat <- t$statistic
  p <- t$p.value
  
  # return(c(t=stat, p=p, sig=(p < .05)))
  return(sig=(p < .05))
  # return a named vector with the results we want to keep
}

#' Bonferroni Holm Function
#' 
#' Performs an adjusted bonferroni correction
#' 
#' @param p_values a vector of p_value
#' 
#' @param alpha alpha value to use for adjustment
#' 
#' @param method method for adjustment. Default is holm
#' 
#' @export
#' 
#' @return 
#' a logical vector indicating which p-values are statistically significant.
#' 
#' @examples 
#' Lets say we have performed a series of statistical tests 
#' and obtained the following p-values:
#' p.vals = seq(0.0025,0.025,0.0025)
#' 
#' We could give them names, to identify the results of each test.
#' names(p.vals) = paste("Test",LETTERS[1:10])
#'
#' bhAdjust(p.vals)
#' 
bhAdjust = function(p_values, alpha=0.5, method= "holm"){
  new_p_values = p.adjust(p_values, method = method)
  logi_values = new_p_values < alpha
  return(logi_values)
}


#' FDR Function
#' 
#' Performs an adjusted FDR correction
#' 
#' @param p_values a vector of p_value
#' 
#' @param alpha alpha value to use for adjustment
#' 
#' @param method method for adjustment. Default is FDR
#' 
#' @return 
#' a logical vector indicating which p-values are statistically significant.
#' 
#' @export
#' @examples 
#' Lets say we have performed a series of statistical tests 
#' and obtained the following p-values:
#' p.vals = seq(0.0025,0.025,0.0025)
#' 
#' We could give them names, to identify the results of each test.
#' names(p.vals) = paste("Test",LETTERS[1:10])
#'
#' fdrAdjust(p_values = p.vals, alpha = 0.5, method = "fdr")
#' 
fdrAdjust = function(p_values, alpha=0.5, method= "fdr"){
  new_p_values = p.adjust(sort(p_values) , method = method)
  new_p_values = sort(new_p_values)
  logi_values = new_p_values < alpha
  return(logi_values)
}


#' R Squared Function
#'
#' Calculates the r-squared value 
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated target values
#' 
#' @return 
#' outputs the r-square value
#' 
#' 
#' @export
#' @examples 
#' y_true = c(4, -1.5, 5, 9)
#' y_pred = c(3.5, 1.0, 5, 10)
#' r2(y_pred, y_true)
#' 
r2 = function(y_pred, y_true){
  rss = sum((y_true - y_pred)**2)
  tss = sum((y_true - mean(y_true))**2)
  results = 1 - (rss/tss)
  return(results)
}

#' Adjusted R Squared Function
#'
#' Calculates the adjusted r-squared value 
#' 
#' @param y_true array-like of shape (n_samples,) Ground truth (correct) target values.
#' 
#' @param y_pred array-like of shape (n_samples,) Estimated target values
#' 
#' @return 
#' outputs the r-square value
#' 
#' 
#' @export
#' @examples 
#' y_true = c(4, -1.5, 5, 9)
#' y_pred = c(3.5, 1.0, 5, 10)
#' adjR2(y_pred, y_true)
#' 
adjR2 = function(y_true, y_pred){
  rss = sum((y_true - y_pred)**2)
  tss = sum((y_true - mean(y_true))**2)
  results = 1 - ((rss/(length(y_true)-1-1))/(tss/(length(y_true)-1)))
  return(results)
}


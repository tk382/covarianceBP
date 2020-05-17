get_grad = function(x, y1, y2){
  a = mean(y1^2)
  b = mean(y2^2)
  d = mean(y1*y2)
  xmat = cbind(rep(1, length(y1)), x)
  first = (a*b*d - d^3) * colSums(xmat)
  second = b * d * colSums(xmat * (y1^2))
  third = a * d * colSums(xmat * (y2^2))
  fourth = (a*b+d^2)*colSums(xmat * y1 * y2)
  out = (first - second - third + fourth)
  return(out)
}

get_fisher = function(x, y1, y2){
  a = mean(y1^2)
  b = mean(y2^2)
  d = mean(y1*y2)
  xmat = cbind(rep(1, length(y1)), x)
  first = (a*b+d^2)*(a*b-d^2)*(t(xmat) %*% xmat)
  second = 4*a*b*d*d*
    (colSums(xmat) %*% t(colSums(xmat)))/nrow(x)
  return(solve(first-second))
}


get_score = function(x, y1, y2){
  a = mean(y1^2)
  b = mean(y2^2)
  d = mean(y1*y2)
  const = 1/(a*b-d^2)
  return(const* 
           t(get_grad(x, y1, y2))  %*%
           get_fisher(x,y1,y2) %*%
           get_grad(x,y1,y2))
}

get_degree = function(x, y, Y){
  d = 0
  for (k in 1:ncol(Y)){
    d = d + get_score(x, y, Y[,k])
  }
  return(d)
}

get_eta = function(rho12, rho23, rho13, rho11, rho22, rho33){
  # num = (rho11*rho23 + rho12*rho13)
  # denom = sqrt(rho11*rho22-rho12^2)*sqrt(rho11*rho33-rho13^2)
  # return(num/denom)
  a = rho11
  b = rho22
  c = rho33
  d = rho12
  e = rho13
  f = rho23
  num = 2*a^2*b*c*d*e + 3*a*c*d^3*e + 3*a*b*d*e^3-2*a^2*b^2*c*e
  num = num - a^2*b*e^2*f - a*b^2*e^3 + a^3*b*c*f
  num = num - 2*a*b*c*d^2*e - b*d^2*e^3 - a^2*c*d^2*f+2*d^3*e^3
  num = num - a^2*b*c*d^2 - a*b*d^2*e^2 + 2*a^2*d*e*f^2
  denom = sqrt((a*b+d^2)*(a*b-d^2)^2*(a*b+e^2)*(a*b-e^2)^2)
  return(num/denom)
}

#' Estimate H
#'
#' @param Sigma True or estimated variance of all the variables of interest.
#' The first row / column must be the main variable of interest
#' @return Esitmated covariance matrix of the part of score statistics
#' @export
get_H = function(Sigma){
  K = nrow(Sigma)
  est_H = matrix(NA, K-1, K-1)
  for (i in 2:(K-1)){
    for (j in (i+1):K){
      est_H[i-1, j-1] = get_eta(Sigma[1,i], Sigma[i,j], Sigma[j,1],
                                Sigma[1,1], Sigma[i,i], Sigma[j,j])
      est_H[j-1, i-1] = est_H[i-1, j-1]
    }
  }
  diag(est_H) = 1
  return(est_H)
}


#' Estimate p-value from the degree statistic
#'
#' @param y vector: main variable of interest
#' @param Y Matrix: other variables of interest to measure the correlation with y
#' @param d double: computed degree statistic
#' @param numsim integer: number of simulations to draw the null distribution
#' @return p value of the degree statistic
#' @export
#' @examples
#' n = 30
#' k = 3
#' nullX = as.matrix(rnorm(30), ncol=1)
#' nullY = MASS::mvrnorm(n, rep(0,3), diag(3))
#' d = get_degree(as.matrix(nullX, ncol=1), nullY[,1], nullY[,2:3])
#' p = get_p_from_degree(nullY[,1], nullY[,2:3], d)
#' print(paste("The degree statistic is", d, "and the p-value is", p))
get_p_from_degree = function(y, Y, d, numsim = 5000){
  bigy = cbind(y, Y)
  K = ncol(bigy)
  est_Sigma = stats::cor(bigy)
  est_H = get_H(est_Sigma)
  lambda = eigen(est_H)$values
  U      = eigen(est_H)$vectors
  null_d = matrix(NA, numsim, K-1)
  for (k in 1:(K-1)){
    null_d[,k] = rgamma(numsim, 1/2, 1/(2*lambda[k]))
  }
  p = sum(rowSums(null_d) > as.numeric(d))/numsim
  return(p)
}



#' Computes the small sample correction
#'
#' @param score score statistic computed from get_score
#' @param coef coefficients for the cubic function
#' @return small-sample adjusted score statistic
#' @export
post_score = function(score, coef) {
  roots = polyroot(c(-score, coef))
  return(Re(roots)[abs(Im(roots)) < 1e-06][1])
}


#' Helper function for get_est_H
#' Returns each element of matrix H
#'
#' @param rho12 correlation between variable 1 and variable 2
#' @param rho23 correlation between variable 2 and variable 3
#' @param rho13 correlation between variable 1 and variable 3
#' @return element (j,k) of H from elements (i,j), (j,k), (i,k) of Sigma
#' @export
get_cor_r = function(rho12, rho23, rho13) {
  num = (rho23 + 2 * rho12 * rho23) *
    (rho12^2 + 1) * (rho13^2 + 1)
  num = num + rho12 * rho13 * (6 + 2 *rho12 + 2 * rho13 + 2 * rho23)
  num = num - rho12 * (rho13^2 + 1) * (3 * rho13 + rho13 + 2 * rho12 *rho23)
  num = num - rho13 * (rho12^2 + 1) * (3 * rho12 + rho12 + 2 * rho13 * rho23)
  denom = (1 - rho12^2) * (1 - rho13^2) * sqrt(1 + rho12^2) * sqrt(1 + rho13^2)
  return(num/denom)
}


#' Computes the correlation matrix of the combined test statistics r
#' Takes input of the estimated or true covariance matrix
#'
#' @param Sigma estimate or true covariance matrix of $K$ variables
#' @return estimated H from the covariance matrix Sigma
#' @export
get_est_H = function(Sigma) {
  K = nrow(Sigma)
  est_H = matrix(0, K - 1, K - 1)
  diag(est_H) = 1
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      eta = get_cor_r(Sigma[1, i],
                      Sigma[i, j], Sigma[j, 1])
      est_H[i - 1, j - 1] = est_H[j -
                                    1, i - 1] = eta
    }
  }
  return(est_H)
}






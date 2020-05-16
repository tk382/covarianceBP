#### rho generating functions ####
h1_fisher = function(X, alpha) {
  tmp = X %*% alpha
  return((exp(tmp) - 1)/(exp(tmp) + 1))
}

h6_quadratic = function(X, alpha) {
  sigma1 = (X %*% alpha - 0.1)^2 - 0.99
  return(sigma1)
}

####

X = matrix(rnorm(n), nrow = n)
rho = matrix(0, n, 4)
alphas = c(.1, .3, .5, .7)
for (i in 1:4){
  rho[,i] = h6_quadratic(x, alphas[i])
}
rho = pmin(rho, 0.9)
par(mfrow = c(2,2))
plot(rho[,1] ~ x)
plot(rho[,2] ~ x)
plot(rho[,3] ~ x)
plot(rho[,4] ~ x)

set.seed(2020)
B = 1000
s1 = s2 = matrix(NA, B, 4)
for (a in 1:4){
  for (b in 1:B){
    Y = matrix(NA, n, 2)
    for (i in 1:n){
      Sigma = matrix(c(1, rho[i,a], rho[i,a], 1), nrow=2)
      Y[i,] = MASS::mvrnorm(n = 1, mu = rep(0,2), Sigma = Sigma)
    }
    y1 = Y[,1]
    y2 = Y[,2]
    s1[b,a] = get_score(x, y1, y2)
    s2[b,a] = get_q(y1, y2, x, correction = FALSE)
  }
  print(a)
}

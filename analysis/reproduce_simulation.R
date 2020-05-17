library(diversitree)
library(mclust)
library(ggplot2)
library(reshape2)
Rcpp::sourceCpp("src/utilities.cpp")

#### rho generating functions ####
h1_fisher = function(X, alpha) {
  X = cbind(rep(1, nrow(X)), X)
  tmp = X %*% alpha
  return((exp(tmp) - 1)/(exp(tmp) + 1))
}

h6_quadratic = function(X, alpha) {
  X = cbind(rep(1, nrow(X)), X)
  sigma1 = (X %*% alpha - 0.1)^2 - 0.99
  return(sigma1)
}


#### Likelihood function ####
likelihood_fisher = function(alpha, opt = list(X, Y)){
  rho = h1_fisher(opt$X, alpha)
  # sigma1 = mean(opt$Y[,1]^2)
  # sigma2 = mean(opt$Y[,2]^2)
  l = 0
  for (i in 1:nrow(opt$X)){
    l = l + dmvnorm(matrix(opt$Y[i,], nrow=1), 
                    rep(0,2), 
                    matrix(c(1, rho[i], rho[i], 1), nrow=2),
                    log=TRUE)
  }
  return(-l)
}
likelihood_null = function(alpha, opt = list(X, Y)){
  newrho = h1_fisher(X, c(alpha, 0))[1]
  sigma1 = mean(Y[,1]^2)
  sigma2 = mean(Y[,2]^2)
  l = sum(dmvnorm(Y, 
                  rep(0,2), 
                  matrix(c(1, newrho, newrho, 1), nrow=2),
                  log=TRUE))
  return(-l)
}

#### run simulation for FISHER ####
## set up
set.seed(20200516)
N = 70
B= 1000
alphalist = c(0, 0.25, 0.5, 0.75, 1)
X = matrix(rnorm(N), ncol = 1)
shuffle = matrix(NA, N, 1000)
for (d in 1:1000){
  shuffle[,d] = sample(1:N)
}

## output
scorep = lap = lrp = matrix(0, B, length(alphalist))
scoret = lat = lrt = matrix(0, B, length(alphalist))

for (a in 1:length(alphalist)){
  alpha = c(0, alphalist[a])
  print(alpha)
  # generate rho
  rho = h1_fisher(X, alpha)
  Y = matrix(NA, N, 2)
  for (b in 1:B){
    # generate Y
    for (c in 1:N){
      Sigma = matrix(rho[c], 2, 2); diag(Sigma) = 1
      Y[c,] = mvrnormArma(1, rep(0,2), Sigma)
    }
    ## LRT
    # print("LR..")
    t = Sys.time()
    mle_fisher = optim(c(0, alpha[1]), 
                       likelihood_fisher,
                       opt = list(X = X, Y = Y),
                       method = "L-BFGS-B",
                       lower=-10, upper=10)
    mle_null = optim(mean(Y[,1]*Y[,2]),
                     likelihood_null,
                     opt = list(X = X, Y = Y),
                     method = "Brent",
                     lower = -5, upper= 5)
    llr_fisher = 2 * (-mle_fisher$value + mle_null$value)
    lrp[b, a] = pchisq(llr_fisher, 1, lower.tail=FALSE)
    lrt[b, a] = Sys.time() - t
    
    ## LM
    # print("LM..")
    t = Sys.time()
    score = get_score_c(x = X, y1 = Y[,1], y2 = Y[,2])
    scorep[b, a] = pchisq(score, 1, lower.tail = FALSE)
    scoret[b, a] = Sys.time() - t
    
    ## LA
    # print("LA..")
    t = Sys.time()
    truela = mean(Y[,1]*Y[,2]*X)
    la_null = rep(0, 1000)
    for (d in 1:1000){
      la_null[d] = mean(Y[,1] * Y[,2] * X[shuffle[,d]])
    }
    lap[b,a] = sum(abs(la_null) > abs(truela)) / 1000
    lat[b,a] = Sys.time() - t
  }
}

## analyze time
lrt = as.data.frame(lrt)
lat = as.data.frame(lat)
scoret = as.data.frame(scoret)
colnames(lrt) = colnames(lat) = colnames(scoret) = paste0("alpha=",alphalist)
lrt$test = "Likelihood Ratio"
lat$test = "Liquid Association"
scoret$test = "Score"
fisher_time = rbind(lrp, lap, scorep)
fisher_time2 = reshape2::melt(fisher_time, id = "test")
colnames(fisher_time2) = c("test", "alpha", "pvalue")
fisher_time2$generate = "HyperbolicTangent"
write.table(fisher_time2, "output/simulation_fishertime.txt",
            col.names=TRUE, row.names=FALSE, quote=TRUE)

##analyze p-value
lrp = as.data.frame(lrp)
lap = as.data.frame(lap)
scorep = as.data.frame(scorep)
colnames(lrp) = colnames(lap) = colnames(scorep) = paste0("alpha=",alphalist)
lrp$test = "Likelihood Ratio"
lap$test = "Liquid Association"
scorep$test = "Score"
fisher_result = rbind(lrp, lap, scorep)
fisher_result2 = reshape2::melt(fisher_result, id = "test")
colnames(fisher_result2) = c("test", "alpha", "pvalue")
fisher_result2$generate = "HyperbolicTangent"
write.table(fisher_result2, "output/simulation_fisherresult.txt",
            col.names=TRUE, row.names=FALSE, quote=TRUE)
ggplot(fisher_result2, aes(x = alpha, y = pvalue, col = test, fill = test)) + 
  geom_violin() + 
  geom_jitter(aes(col = test), size = 0.1, alpha = 0.1) 
  
library(dplyr)
fisher_result2 %>% group_by(test, alpha) %>%
  summarize(power = sum(pvalue < 0.05)/B)


#### run simulation for QUADRATIC ####
## set up
set.seed(20200516)
N = 70
B= 1000
alphalist = c(0.2, 0.3, 0.4, 0.5)
X = matrix(runif(N), ncol = 1)
shuffle = matrix(NA, N, 1000)
for (d in 1:1000){
  shuffle[,d] = sample(1:N)
}
## output
scorep = lap = lrp = matrix(0, B, length(alphalist))
scoret = lat = lrt = matrix(0, B, length(alphalist))

for (a in 1:length(alphalist)){
  alpha = c(0, alphalist[a])
  print(alpha)
  # generate rho
  rho = h6_quadratic(X, alpha)
  Y = matrix(NA, N, 2)
  for (b in 1:B){
    # generate Y
    for (c in 1:N){
      Sigma = matrix(rho[c], 2, 2); diag(Sigma) = 1
      Y[c,] = mvrnormArma(1, rep(0,2), Sigma)
    }
    ## LRT
    # print("LR..")
    t = Sys.time()
    mle_fisher = optim(c(0, alpha[1]), 
                       likelihood_fisher,
                       opt = list(X = X, Y = Y),
                       method = "L-BFGS-B",
                       lower=-10, upper=10)
    mle_null = optim(mean(Y[,1]*Y[,2]),
                     likelihood_null,
                     opt = list(X = X, Y = Y),
                     method = "Brent",
                     lower = -5, upper= 5)
    llr_fisher = 2 * (-mle_fisher$value + mle_null$value)
    lrp[b, a] = pchisq(llr_fisher, 1, lower.tail=FALSE)
    lrt[b, a] = Sys.time() - t
    
    ## LM
    # print("LM..")
    t = Sys.time()
    score = get_score_c(x = X, y1 = Y[,1], y2 = Y[,2])
    scorep[b, a] = pchisq(score, 1, lower.tail = FALSE)
    scoret[b, a] = Sys.time() - t
    
    ## LA
    # print("LA..")
    t = Sys.time()
    truela = mean(Y[,1]*Y[,2]*X)
    la_null = rep(0, 1000)
    for (d in 1:1000){
      la_null[d] = mean(Y[,1] * Y[,2] * X[shuffle[,d]])
    }
    lap[b,a] = sum(abs(la_null) > abs(truela)) / 1000
    lat[b,a] = Sys.time() - t
  }
}


## analyze time
lrt = as.data.frame(lrt)
lat = as.data.frame(lat)
scoret = as.data.frame(scoret)
colnames(lrt) = colnames(lat) = colnames(scoret) = paste0("alpha=",alphalist)
lrt$test = "Likelihood Ratio"
lat$test = "Liquid Association"
scoret$test = "Score"
quadratic_time = rbind(lrt, lat, scoret)
quadratic_time2 = reshape2::melt(quadratic_time, id = "test")
colnames(quadratic_time2) = c("test", "alpha", "time", "generate")
quadratic_time2$generate = "Quadratic"
write.table(quadratic_time2, "output/simulation_quadratictime.txt",
            col.names=TRUE, row.names=FALSE, quote=TRUE)
ggplot(quadratic_time2, aes(x = test, y = time, col = alpha)) + 
  geom_violin() + 
  geom_jitter(aes(col = alpha), size = 0.1, alpha = 0.3)  + 
  scale_y_continuous(trans='log',
                     breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1)) + 
  ggsci::scale_color_startrek() +
  theme_bw() 



##analyze p-value
lrp = as.data.frame(lrp)
lap = as.data.frame(lap)
scorep = as.data.frame(scorep)
colnames(lrp) = colnames(lap) = colnames(scorep) = paste0("alpha=",alphalist)
lrp$test = "Likelihood Ratio"
lap$test = "Liquid Association"
scorep$test = "Score"
quadratic_result = rbind(lrp, lap, scorep)
quadratic_result2 = reshape2::melt(quadratic_result, id = "test")
colnames(quadratic_result2) = c("test", "alpha", "pvalue")
quadratic_result2$generate = "HyperbolicTangent"
write.table(quadratic_result2, "output/simulation_quadraticresult.txt",
            col.names=TRUE, row.names=FALSE, quote=TRUE)
ggplot(quadratic_result2, aes(x = alpha, y = pvalue, col = test, fill = test)) + 
  geom_violin() + 
  geom_jitter(aes(col = test), size = 0.1, alpha = 0.1) 

library(dplyr)
quadratic_result2 %>% group_by(test, alpha) %>%
  summarize(power = sum(pvalue < 0.05)/B)

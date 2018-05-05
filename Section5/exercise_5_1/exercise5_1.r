setwd("C:\\Users\\zc3763\\Dropbox\\Courses\\18SP\\SDS383D\\Section5\\exercise_5_1")
library(mvtnorm)
library(MASS)
profitData = read.csv("restaurants.csv",header = TRUE)

N = 1000
X = as.matrix(profitData$SeatingCapacity)
Y = as.matrix(profitData$Profit)
intercept = rep(1,nrow(X))
X = cbind(intercept,X)
X = cbind(X,profitData$DinnerService)


mu_n = solve(t(X) %*% X + diag(1,ncol(X))) %*% (t(X) %*% Y)
Sigma_nw = solve(t(X) %*% X + diag(1,ncol(X)))
orate = 1 + 0.5 * (t(Y) %*% Y - t(mu_n) %*% (t(X) %*% X + diag(1,ncol(X),ncol(X))) %*% mu_n)

beta = matrix(0,ncol(X),N + 1)
omega = rgamma(N,shape = (1 + 0.5 * (nrow(X) + ncol(X))), rate = orate)

for (index in 1:N){
  
  beta[,index + 1] = t(rmvnorm(1, mean = mu_n, sigma = Sigma_nw / omega[index]))

}
beta = beta[,500:ncol(beta)]
betaSelect = rowMeans(beta)
yBayes = X %*% betaSelect
residBayes = Y - yBayes
mean(residBayes ** 2)
hist(residBayes,breaks = 50)
setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section2\\R")
library(mvtnorm)
library(MASS)
dentalData = read.csv("dental.csv",header = TRUE)

N = 3000
X = as.matrix(dentalData$age)
Y = as.matrix(dentalData$distance)
intercept = rep(1,nrow(X))
X = cbind(intercept,X)
dentalData$intercept = intercept
dentalData$genderFlag = rep(1,nrow(dentalData))
dentalData$genderFlag[dentalData$Sex == 'Female'] = 0
X = cbind(X,dentalData$genderFlag)

mu_n = solve(t(X) %*% X + diag(1,ncol(X))) %*% (t(X) %*% Y)
Sigma_nw = solve(t(X) %*% X + diag(1,ncol(X)))
orate = 1 + 0.5 * (t(Y) %*% Y - t(mu_n) %*% (t(X) %*% X + diag(1,ncol(X),ncol(X))) %*% mu_n)

beta = matrix(0,ncol(X),N + 1)
beta[,1] = 1
omega = matrix(1,ncol(X),N + 1)
Lambda = matrix(1,nrow(X),N + 1)

for (index in 1:N){

	mu_n = solve(t(X) %*% diag(Lambda[,index]) %*% X + diag(1,ncol(X))) %*% (t(X) %*% diag(Lambda[,index]) %*% Y)
	Sigma_nw = solve(t(X) %*% diag(Lambda[,index]) %*% X + diag(1,ncol(X)))
	orate = 1 + 0.5 * (t(Y) %*% diag(Lambda[,index]) %*% Y - t(mu_n) %*% (t(X) %*% diag(Lambda[,index]) %*% X + diag(1,ncol(X),ncol(X))) %*% mu_n)
	omega[index + 1] = rgamma(1,shape = (1 + 0.5 * (nrow(X) + ncol(X))), rate = orate)
	Lambda[,index + 1] = rgamma(nrow(X), shape = 1.5, rate = 0.5 + omega[index] / 2 * (Y - X %*% beta[,index]) ** 2)
  
	beta[,index + 1] = t(rmvnorm(1, mean = mu_n, sigma = Sigma_nw / omega[index]))

}
beta = beta[,501:ncol(beta)]
betaSelect = rowMeans(beta)
yBayes = X %*% betaSelect
residBayes = Y - yBayes
sqrt(mean(residBayes ** 2))
hist(residBayes)

ols_mdl = lm("distance ~ age", data = dentalData)
residOLS = ols_mdl$residuals
hist(residOLS)

ridge_mdl = lm.ridge("distance ~ intercept + age + 0", lambda = 0.2, data = dentalData)
yRidge = X %*% ridge_mdl$coef
residRidge = Y - yRidge
#hist(residRidge)
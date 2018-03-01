setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section2\\R")
library(mvtnorm)
library(MASS)
dentalData = read.csv("dental.csv",header = TRUE)

N = 1000

intercept = rep(1,nrow(dentalData))

dentalData$intercept = intercept
dentalData$genderFlag = rep(1,nrow(dentalData))
dentalData$genderFlag[dentalData$Sex == 'Female'] = 0

maleData = dentalData[dentalData$Sex == 'Male',]
femaleData = dentalData[dentalData$Sex == 'Female',]

X = cbind(maleData$intercept,maleData$age)
Y = maleData$distance

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
hist(residBayes)

ols_mdl = lm("distance ~ age", data = maleData)
residOLS = ols_mdl$residuals
mean(residOLS ** 2)
hist(residOLS)

ridge_mdl = lm.ridge("distance ~ intercept + age + 0", lambda = 0.2, data = dentalData)
yRidge = X %*% ridge_mdl$coef
residRidge = Y - yRidge
hist(residRidge)
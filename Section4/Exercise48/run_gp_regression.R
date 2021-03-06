library(readr)
library(rstan)
setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section4\\Exercise48")

data(faithful)

gp_data<-read_rdump('faithful_data.R')

filename = "gp_regression.stan"

stan_code <-readChar(filename, file.info(filename)$size)

fit <- stan(model_code=stan_code,data=gp_data)

params<-extract(fit)

par(mfrow=c(1, 3))

alpha_breaks=10 * (0:50) / 50 - 5
hist(log(params$alpha), main="", xlab="log(alpha)", col="red", yaxt='n')
quantile(params$alpha, probs = c(0.025,0.0975))

beta_breaks=10 * (0:50) / 50 - 5
hist(log(params$rho), main="", xlab="log(rho)", col="red", yaxt='n')
quantile(params$rho, probs = c(0.025,0.0975))

sigma_breaks=5 * (0:50) / 50
hist(params$sigma, main="", xlab="sigma", col="red", yaxt='n')
quantile(params$sigma, probs = c(0.025,0.0975))

readline(prompt="Press [enter] to continue")

probs = c(0.1,0.5,0.9)
cred <- sapply(1:length(gp_data$x_predict), function(n) quantile(params$y_predict[,n], probs=probs))

plot(1, type="n", main=title,
     xlab="x", ylab="y", xlim=c(0, 120), ylim=c(-10, 10))
polygon(c(gp_data$x_predict, rev(gp_data$x_predict)), c(cred[1,], rev(cred[3,])),
        col = c_light, border = NA)
lines(gp_data$x_predict, cred[2,], col=c_dark, lwd=2)

points(gp_data$x,gp_data$y)
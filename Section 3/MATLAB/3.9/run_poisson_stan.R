library(readr)
library(rstan)
library(Metrics)
setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section3\\MATLAB\\3.9")
tea_discipline_oss <- read_csv("tea_discipline_oss.csv") 
#View(tea_discipline_oss)

uncensored_data = subset(tea_discipline_oss,ACTIONS>0)
tea <-data.frame(grade=uncensored_data$GRADE,y=uncensored_data$ACTIONS)
tea$intercept =1
tea<-as.list(tea)
tea$N<-nrow(uncensored_data)
fileName <- "poisson.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
resStan<-stan(model_code=stan_code,data=tea,chains=3,iter=3000,warmup=1000,thin=10)
traceplot(resStan, pars = c("beta"), inc_warmup = FALSE) #set inc_warmup = TRUE to see burn in
beta1 = mean(resStan@sim[["samples"]][[1]][["beta[1]"]]+resStan@sim[["samples"]][[2]][["beta[1]"]]+resStan@sim[["samples"]][[3]][["beta[1]"]])/3
beta2 = mean(resStan@sim[["samples"]][[1]][["beta[2]"]]+resStan@sim[["samples"]][[2]][["beta[2]"]]+resStan@sim[["samples"]][[3]][["beta[2]"]])/3
yhat = beta2*uncensored_data$GRADE +beta1
rmse(uncensored_data$ACTIONS,yhat)
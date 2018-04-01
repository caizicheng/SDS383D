library(readr)
library(rstan)
library(Metrics)
setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section3\\MATLAB\\3.9")
tea_discipline_oss <- read_csv("tea_discipline_oss.csv") 
#View(tea_discipline_oss)

uncensored_data = subset(tea_discipline_oss,ACTIONS>0)
gender = uncensored_data$SEXX
gender[gender == 'FEMALE'] = 0
gender[gender == 'MALE'] = 1
gender = as.integer(gender)
genderattendance = gender*uncensored_data$SE_ATTEND
gendergrade = gender*uncensored_data$GRADE
tea <-data.frame(grade=uncensored_data$GRADE,se_attend=uncensored_data$SE_ATTEND,gender=gender,y=uncensored_data$ACTIONS,genderattendance=genderattendance,gendergrade=gendergrade)
tea$intercept =1
tea<-as.list(tea)
tea$N<-nrow(uncensored_data)
fileName <- "poisson_2.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
resStan<-stan(model_code=stan_code,data=tea,chains=1,iter=1000,warmup=100,thin=10)
traceplot(resStan, pars = c("beta"), inc_warmup = FALSE) #set inc_warmup = TRUE to see burn in
beta1 = mean(resStan@sim[["samples"]][[1]][["beta[1]"]][50:100])
beta2 = mean(resStan@sim[["samples"]][[1]][["beta[2]"]][50:100])
beta3 = mean(resStan@sim[["samples"]][[1]][["beta[3]"]][50:100])
beta4 = mean(resStan@sim[["samples"]][[1]][["beta[4]"]][50:100])
beta5 = mean(resStan@sim[["samples"]][[1]][["beta[5]"]][50:100])
beta6 = mean(resStan@sim[["samples"]][[1]][["beta[6]"]][50:100])
yhat = beta1+beta2*uncensored_data$GRADE+beta3*gender+beta4*uncensored_data$SE_ATTEND+beta5*genderattendance+beta6*gendergrade
rmse(uncensored_data$ACTIONS,yhat)
library(readr)
library(rstan)
setwd("C:\\Users\\Zicheng Cai\\Dropbox\\Courses\\18SP\\SDS383D\\Section3\\MATLAB\\3.9")
tea_discipline_oss <- read_csv("tea_discipline_oss.csv") 
#View(tea_discipline_oss)

uncensored_data = subset(tea_discipline_oss,ACTIONS>0)
gender = uncensored_data$SEXX
gender[gender == 'FEMALE'] = 0
gender[gender == 'MALE'] = 1
gender = as.integer(gender)
tea <-data.frame(grade=uncensored_data$GRADE,se_attend=uncensored_data$SE_ATTEND,gender=gender,y=uncensored_data$ACTIONS)
tea$intercept =1
tea<-as.list(tea)
tea$N<-nrow(uncensored_data)
fileName <- "poisson_1.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)
resStan<-stan(model_code=stan_code,data=tea,chains=3,iter=3000,warmup=1000,thin=10)
traceplot(resStan, pars = c("beta"), inc_warmup = FALSE) #set inc_warmup = TRUE to see burn in
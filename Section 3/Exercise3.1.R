setwd("C:\\Users\\zc3763\\Dropbox\\Courses\\18SP\\SDS383D\\git\\SDS383D\\Section 3")
rawData = read.csv("pima.csv")
rawData$intercept = rep(1,nrow(rawData))
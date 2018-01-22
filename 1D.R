n = 60
numCycles = 50
for (index in 1:numCycles){
  
  samples = ceiling(runif(n,min = 1,max = 3000))
  # Twice the sample mean
  estimator1[index] = 2 * mean(samples)
  
  # Sample mean plus 3x sample standard deviation
  estimator2[index] = mean(samples) + 3*sd(samples)
  
  # (n+1)/n times the sample maximum
  estimator3[index] = (n+1)/n * max(samples)
  
}
mean(estimator1)
var(estimator1)
t.test(estimator1,mu=3000)
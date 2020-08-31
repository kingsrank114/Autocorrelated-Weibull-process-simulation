#Estimation of Autocorrelation function up to lag d of sample Y with standard deviation sigma
acf_method <- function(y, sigma, d){
  
  acf <- rep(NA, d)
  
  acf[1] <- 1
  
  l <- length(y)
  
  for(i in 2:d){
    
    acf[i] <- cov(y[i:l], y[1:(l-i+1)])/sigma^2 
    
  }
  
  return(acf=acf)
  
}


#Generate first-order autocorrelated Weibull Distribution process

#Input
#n --- Number of sample size desired
#k --- Weibull distribution shape parameter (k > 0)
#lambda --- Weibull distribution scale parameter (lambda > 0)
#rho --- first-order autocorrelation 

#Output
#mean_parameter --- The theoretical mean for Weibull distribution   
#sd_parameter --- The theoretical standard deviation for Weibull distribution 
#mean_simulation --- The mean of simulated sample
#sd_parameter --- The standard deviation of simulated sample
#rho1_est --- Estimated first-order autocorrelation

#Conditional distribution

#CDF function of first-order autocorrelated Weibull distribution
library(zipfR)

F1 <- function(y, k, lambda, rho, y_previous, u){
  
  mean <- lambda * gamma(1 + 1/k)
  
  term <-  lambda  * (gamma(1 + 2/k) - (gamma(1+1/k))^2) 
  
  pweibull(y, k, lambda, lower.tail = TRUE) + (rho * (y_previous - mean) / term) * 
    ((Igamma((k+1)/k, (y/lambda)^k) - Igamma((k+1)/k, 0))  -   gamma(1+1/k) * 
       pweibull(y, k, lambda, lower.tail = TRUE)) - u
  
}

First_order_simulation <- function(n, k, lambda, rho){
  
  Y <- rep(NA, n)
  term <- rep(NA, n-1)
  Y[1] <- rweibull(1, k, lambda)
  
  #Data simulation process
  for(i in 2:n){
    
    u <- runif(1, 0, 1)
    
    Y[i] <- uniroot(F1, c(0, 1000), k, lambda, rho, Y[i-1],u)$root
    
  }
  
  #Checking for the reliability of the simulation!!
  for(i in 2:n){
    
    mean_PARAMETER <- lambda * gamma(1+1/k)
    mean_simulation <- mean(Y)
    
    sd_PARAMETER <- sqrt(lambda^2* (gamma(1+2/k) - (gamma(1+1/k))^2))
    sd_simulation <- sd(Y)
    
    Zt <- (Y[i] - mean_PARAMETER)/sd_PARAMETER
    Zt1 <- (Y[i-1] - mean_PARAMETER)/ sd_PARAMETER
    
    term[i-1] <- Zt * Zt1
  }
  
  list(Y=Y, rho1_est = mean(term), mean_PARAMETER=mean_PARAMETER, mean_simulation=mean_simulation,
       sd_PARAMETER=sd_PARAMETER, sd_simulation=sd_simulation)
  
}

n = 10000 
lambda = 1
k = 0.5
rho = 0.15
sim <- First_order_simulation(n, k, lambda, rho)
Y <- sim$Y
rho1_est <- sim$rho1_est
rho1_est
mean_PARAMETER <- sim$mean_PARAMETER
mean_PARAMETER
mean_simulation <- sim$mean_simulation
mean_simulation
sd_PARAMETER <- sim$sd_PARAMETER
sd_PARAMETER
sd_simulation <- sim$sd_simulation
sd_simulation

#Picture
n = 10000
rho = 0.15 
sim1 <- First_order_simulation(n, 0.5, 1, rho)
Y1 <- sim1$Y
sim2 <- First_order_simulation(n, 1, 1, rho)
Y2 <- sim2$Y
sim3 <- First_order_simulation(n, 1.5, 1, rho)
Y3 <- sim3$Y
sim4 <- First_order_simulation(n, 5, 1, rho)
Y4 <- sim4$Y

par(mfrow=c(2,2))
hist(Y1, main = "lambda = 1, k = 0.5, rho = 0.15", prob=TRUE, ylim = c(0,1))
lines(density(Y1), col = "purple", lwd=2)

hist(Y2, main = "lambda = 1, k = 1, rho = 0.15", prob=TRUE, ylim = c(0,1))
lines(density(Y2), col = "red", lwd=2)

hist(Y3, main = "lambda = 1, k = 1.5, rho = 0.15", prob=TRUE, ylim = c(0,1))
lines(density(Y3), col = "pink", lwd=2)

hist(Y4, main = "lambda = 1, k = 5, rho = 0.15", prob=TRUE, ylim = c(0,2 ))
lines(density(Y4), col = "green", lwd=2)


#Check the joint distribution a density
lambda = 1
k=1
mean_PARAMETER <- lambda * gamma(1+1/k)
sd_PARAMETER <- sqrt(lambda^2 * (gamma(1+2/k) - (gamma(1+1/k))^2))
integrate(function(y){
  
  sapply(y, function(y){
    
    integrate(function(y_previous){dweibull(y, k, lambda) * dweibull(y_previous, k, lambda) *
        (1 + rho * ((y-mean_PARAMETER)/sd_PARAMETER) * ((y_previous - mean_PARAMETER)/ sd_PARAMETER))}, 0, Inf )$value
    
  })
  
}, 0, Inf)

#Integration Cheking.................................................................
f <- function(y){-lambda*gamma(1+1/k) * y^(k-1) * exp(-(y/lambda)^k)}
integrate(f, 0, 100)
y=100
-lambda^(k+1) * gamma(1+1/k)*  pweibull(y, k, lambda, lower.tail = TRUE) / k

f <- function(y){y^(k) * exp(-(y/lambda)^k)}
integrate(f, 0, 1)
y=1
lambda^(k+1) * (Igamma((k+1)/k, (y/lambda)^k) - Igamma((k+1)/k, 0)) / k


mean = 3
var = 1
y_previous = 2.8
f <- function(y){
  mean <- lambda * gamma(1 + 1/k)
  
  var <-  lambda^2 * (gamma(1 + 2/k) - (gamma(1+1/k))^2) 
  
  zt <- (y - mean) / sqrt(var)
  zt1 <- (y_previous - mean) / sqrt(var)
  dweibull(y, k, lambda) * (1 + rho * zt * zt1)
  
}
integrate(f, 0, 1)

y=1
pweibull(y, k, lambda, lower.tail = TRUE) + (rho * (k/lambda) * (y_previous - mean) ) / (lambda^(k-1) * var) *
  (lambda^(k+1) * (-Igamma((k+1)/k, (y/lambda)^k) + Igamma((k+1)/k, 0)) / k -lambda^(k+1) * gamma(1+1/k)*  pweibull(y, k, lambda, lower.tail = TRUE) / k )

#..............................................................................................................................
#.............................................................................................................................................................................................
#Generate Second-order autocorrelated Weibull Distribution

#Input
#n --- Number of sample size desired
#k --- Weibull distribution shape parameter (k > 0)
#lambda --- Weibull distribution scale parameter (lambda > 0)
#rho1 --- first-order autocorrelation 
#rho2 --- second-order autocorrelation
#rho_star --- three-way mixed autocorrelation, i.e., rho_012

#Output
#mean_parameter --- The theoretical mean for Weibull distribution   
#sd_parameter --- The theoretical standard deviation for Weibull distribution 
#mean_simulation --- The mean of simulated sample
#sd_parameter --- The standard deviation of simulated sample
#rho1_est --- Estimated first-order autocorrelation
#rho2_est --- Estimated second-order autocorrelation
#rho_star_est --- Estimated three-way autocorrelation

#CDF function of first-order autocorrelated Weibull distribution
F2 <- function(y, k, lambda, rho1, rho2, rho_star, y_p1, y_p2, u){
  
  mean <- lambda * gamma(1 + 1/k)
  
  sd <-  sqrt(lambda^2  * (gamma(1 + 2/k) - (gamma(1+1/k))^2) )
  
  Zt1 <- (y_p1 - mean)/sd
  Zt2 <- (y_p2 - mean)/sd
  
  pweibull(y, k, lambda, lower.tail = TRUE) + 
    ((rho1 * Zt1 + rho2 * Zt2 + rho_star * Zt1 * Zt2) / (1 + rho1 * Zt1 * Zt2)) * 
    (1 / sqrt(gamma(1 + 2/k) - (gamma(1+1/k))^2)) *
    ((Igamma((k+1)/k, (y/lambda)^k) - Igamma((k+1)/k, 0)) - gamma(1+1/k) * 
       pweibull(y, k, lambda, lower.tail = TRUE)) - u
  
}

second_order_simulation <- function(n, k, lambda, rho1, rho2, rho_star){
  
  Y <- rep(NA, n)
  term1 <- rep(NA, n-2)
  term2 <- rep(NA, n-2)
  term3 <- rep(NA, n-2)
  Y[1] <- rweibull(1, k, lambda)
  u <- runif(1, 0, 1)
  Y[2] <- uniroot(F1, c(0, 1000), k, lambda,rho1,Y[1],u)$root
  
  #Data simulation process
  for(i in 3:n){
    
    u <- runif(1, 0, 1)    
    
    Y[i] <- uniroot(F2, c(0, 1000), k, lambda, rho1, rho2, rho_star, Y[i-1], Y[i-2], u)$root
    
  }
  
  #Checking for the reliability of the simulation!!
  for(i in 3:n){
    
    mean_simulation <- mean(Y)
    sd_simulation <- sd(Y)
    
    mean <- lambda * gamma(1 + 1/k)
    sd <-  sqrt(lambda^2  * (gamma(1 + 2/k) - (gamma(1+1/k))^2) )
    
    Zt <- (Y[i] - mean)/sd
    Zt1 <- (Y[i-1] - mean)/ sd
    Zt2 <- (Y[i-2] - mean)/ sd
    
    term1[i-2] <- Zt * Zt1
    term2[i-2] <- Zt * Zt2
    term3[i-2] <- Zt * Zt1 * Zt2
  }
  
  list(Y=Y, rho1_est = mean(term1), rho2_est = mean(term2), rho_star_est = mean(term3), 
       mean_PARAMETER = mean, mean_simulation=mean_simulation,sd_PARAMETER = sd, 
       sd_simulation=sd_simulation)
  
}


n = 100000
k = 5
lambda= 2
rho1 = 0.3
rho2 = 0.22
rho_star= 0.18

sim2 <- second_order_simulation(n, k, lambda, rho1, rho2, rho_star)
Y <- sim2$Y
rho1_est <- sim2$rho1_est
rho1_est
rho2_est <- sim2$rho2_est
rho2_est
rho_star_est <- sim2$rho_star_est
rho_star_est
mean_PARAMETER <- sim2$mean_PARAMETER
mean_PARAMETER
mean_simulation2 <- sim2$mean_simulation
mean_simulation2
sd_PARAMETER <- sim2$sd_PARAMETER
sd_PARAMETER
sd_simulation2 <- sim2$sd_simulation
sd_simulation2


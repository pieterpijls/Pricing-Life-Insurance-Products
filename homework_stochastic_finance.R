#####################################
### Question 1
#####################################

## Simulate Sample Paths ##

## define model parameters
mu0 <- 0.0015
theta_mu <- 0
k_mu <- -0.08
beta_mu <- 0.008
C = 100000
interest = 0.01

## simulate short rate paths
n <- 1000    # MC simulation trials
T <- 30    # total time
m <- 2000   # subintervals
dt <- T/m  # difference in time each subinterval

mu <- matrix(0,m+1,n)  # matrix to hold short rate paths
mu[1,] <- mu0
set.seed(1)

e1 = matrix(rnorm(m*n,mean=0,sd=1),m,n) #for mortality
e2 = matrix(rnorm(m*n,mean=0,sd=1),m,n) #for short term interest rates
e3 = matrix(rnorm(m*n,mean=0,sd=1),m,n) #for F
e4 = matrix(rnorm(m*n,mean=0,sd=1),m,n) #for G


for(j in 1:n){
  for(i in 2:(m+1)){
    dmu <- k_mu*(theta_mu-mu[i-1,j])*dt + beta_mu*sqrt(dt)*e1[i-1,j]
    mu[i,j] <- abs(mu[i-1,j] + dmu)
  }
} 

### calculate premium
premium_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  premium_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*(1/(1+interest)^T)
}

mean(premium_matrix)

## define model parameters
r0 <- 0.01
theta_r <- 0.02
k_r <- 0.2
beta_r <- 0.012

r <- matrix(0,m+1,n)  # matrix to hold short rate paths
r[1,] <- r0

for(j in 1:n){
  for(i in 2:(m+1)){
    dr <- k_r*(theta_r-r[i-1,j])*dt + beta_r*sqrt(dt)*e2[i-1,j]
    r[i,j] <- r[i-1,j] + dr
  }
} 

### calculate fair value
FV_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*exp(-sum(r[,i]*dt))
}

mean(FV_matrix)


###find equilibrium
rsum=0
for(i in 1:n){
  rsum <- exp(-sum(r[,i]*dt)) + rsum
}
rmean = rsum/n
i = (1/(rmean^(1/T))-1)

#####################################
### Question 2
#####################################

### calculate premium
premium <- C*(1/(1+interest))^T
premium

## define model parameters
p <- 0.85
F0 <- 1
theta_F <- 0
k_F <- -0.04
beta_F <- 0.05
rho <- -0.5

F <- matrix(0,m+1,n)  # matrix to hold short rate paths
F[1,] <- F0

for(j in 1:n){
  for(i in 2:(m+1)){
    dF <- k_F*(theta_F-F[i-1,j])*dt + beta_F*sqrt(dt)*(e2[i-1,j]*rho+sqrt(1-rho^2)*e3[i-1,j])
    F[i,j] <- F[i-1,j] + dF
  }
} 

### calculate fair value
FV_guar_matrix <- matrix(nrow=1,ncol=n)
FV_bonus_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_guar_matrix[1,i] <- C*exp(-sum(r[,i]*dt))*(1+interest)^T
  FV_bonus_matrix[1,i] <- p*premium*exp(-sum(r[,i]*dt))*max(0,((F[2001,i])-(1+interest)^T)) # minus one in this equation??
}

mean(FV_guar_matrix)
mean(FV_bonus_matrix)
FV <- mean(FV_guar_matrix) + mean(FV_bonus_matrix)
FV

###equilibrated contract

i <- (C/FV)^(1/T)-1
#####################################
### Question 3
#####################################
### calculate premium
premium_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  premium_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*(1/(1+interest)^T)
}

premium <- mean(premium_matrix)
premium

### calculate fair value
s <- 0.999441704
g <- 0.999733441
c <- 1.101077536

tpx <- (s^T)*(g^(c^(35+T)-c^35))


FV_guar_matrix <- matrix(nrow=1,ncol=n)
FV_bonus_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_guar_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*exp(-sum(r[,i]*dt)) 
  FV_bonus_matrix[1,i] <- p*C*max(0,(tpx/(exp(-sum(mu[,i]*dt))))-1)*exp(-sum(mu[,i]*dt))*exp(-sum(r[,i]*dt))
}
mean(FV_guar_matrix)
mean(FV_bonus_matrix)
FV <- mean(FV_guar_matrix) + mean(FV_bonus_matrix)
FV

i <- (C/FV)^(1/T)-1
i

#####################################
### Question 4
#####################################

### calculate premium
premium_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  premium_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*(1/(1+interest)^T)
}

premium <- mean(premium_matrix)
premium

### calculate fair value
G0 <- 1
theta_G <- 0
k_G <- -0.08
beta_G <- 0.19
rho_2 <- 0.25

G <- matrix(0,m+1,n)  # matrix to hold short rate paths
G[1,] <- G0

for(j in 1:n){
  for(i in 2:(m+1)){
    e5 <- e2[i-1,j]*rho+sqrt(1-rho^2)*e3[i-1,j]
    dG <- k_G*(theta_G-G[i-1,j])*dt + beta_G*sqrt(dt)*(e5*rho_2+sqrt(1-rho_2^2)*e4[i-1,j])
    G[i,j] <- G[i-1,j] + dG
  }
} 

FV_guar_matrix <- matrix(nrow=1,ncol=n)
FV_bonus_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_guar_matrix[1,i] <- C*exp(-sum(mu[,i]*dt))*exp(-sum(r[,i]*dt)) 
  FV_bonus_matrix[1,i] <- p*C*exp(-sum(r[,i]*dt))*exp(-sum(mu[,i]*dt))*max(0,((G[2001,i])/((1+interest)^T))-1)
}

mean(FV_guar_matrix)
mean(FV_bonus_matrix)
FV <- mean(FV_guar_matrix) + mean(FV_bonus_matrix)
FV

i <- (C/FV)^(1/T)-1
i

#####################################
### Question 5
#####################################

guarantee <- 50000

FV_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_matrix[1,i] <- guarantee*max(1,(G[2001,i]))*exp(-sum(r[,i]*dt))*exp(-sum(mu[,i]*dt))
}

mean(FV_matrix)

#####################################
### Question 6
#####################################

FV_matrix <- matrix(nrow=1,ncol=n)
for(i in 1:n){
  FV_matrix[1,i] <- guarantee*max((F[2001,i]),(G[2001,i]))*exp(-sum(r[,i]*dt))
}

mean(FV_matrix)



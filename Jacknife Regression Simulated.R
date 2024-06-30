rm(list = ls())
n=100
set.seed(2)
# Generate X1 regressor
x=rnorm(n,0,1)
# Specify noncontant variance of error
si2=(x-2)^2/4
max(si2)

# Generate error with oncontant variance
e=0
for(i in 1:n)e[i]=rnorm(1,0,si2[i])

# Generate y with beta_0=2 and beta_1=5
y=2+5*x+e
########################################

Data=data.frame(y,x)
Data

plot(Data)

# Fit OLS regression model
OLS=lm(y~x,data=Data)
summary(OLS)

beta_hat=OLS$coefficients
# Model residuals 
e_hat=OLS$residuals
y_hat=OLS$fitted.values
k=2 # no of parameters
# Estimate of sigma^2
sigma2_hat=sum(e_hat^2)/(n-k)
sigma2_hat
# Design Matrix
X=cbind(rep(1,n),x)
# Variance of beta hat
var_beta=sigma2_hat*solve(t(X)%*%X)
var_beta

# Hat matrix
H=X%*%solve(t(X)%*%X)%*%t(X)
var_e=sigma2_hat*(diag(n)-H)

cbind(diag(var_e), si2, diag(var_e)-si2 )



# Residual v/s Fitted plot
plot(y_hat,e_hat)
###### Jackknife regression
Beta_mat=matrix(0,n,k)

for (i in 1:n) 
{
  model=lm(y~x,data=Data[-i,])
  Beta_mat[i,]=model$coefficients
}
Beta_mat

Beta_bar=colMeans(Beta_mat)



# Difference between beta_i and beta_hat

diff=matrix(0,n,k)
for(i in 1:n) diff[i,]=Beta_mat[i,]-beta_hat

# jackknife variance estimator  
sum1=matrix(0,k,k)
for (i in 1:n) {
  sum1=(t(t(diff[i,])) %*% t(diff[i,]))
}

Var_J=sum1*(n-k)/n
Var_J


# weighted jackknife variance estimator  
sum2=matrix(0,k,k)
for (i in 1:n) {
  sum2=sum2+(1-H[i,i])^2 * (t(t(diff[i,])) %*% t(diff[i,]))
}

Var_H=sum2*n/(n-k)
Var_H

Var_J
var_beta

# Actual variance
Var_beta_act=solve(t(X)%*%X) %*% t(X) %*% diag(e_hat^2) %*% X %*% solve(t(X)%*%X)
Var_beta_act

Var_H

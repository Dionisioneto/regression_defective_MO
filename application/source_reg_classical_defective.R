## ---
## Source functions
## Classical Defective Models
## ---

# 1. Base functions - Gompertz

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}


Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*t)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}


## log-likelihood function for Regression Gompertz

loglik_reg_gpz = function(par,time,delta,xa,xb){
  #parameters
  as = par[1:(dim(xa)[2]+1)] # b0, b1, b2
  bs = par[(dim(xa)[2]+2):(dim(xa)[2]+1+dim(xb)[2]+1)] # d0, d1, d2
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  
  den = as.vector(ft_Gompertz(t=time,alpha=a,beta=b))
  surv = as.vector(St_Gompertz(t=time,alpha=a,beta=b))
  
  log_lik = sum(delta*log(den) + (1-delta)*log(surv))
  return(-1*log_lik)
}


# 2. Base functions - Inverse Gaussian


ft_IG = function(t,alpha,beta){
  den = (1/sqrt(2*pi*beta*t^3))*exp((-(1/(2*beta*t)))*(1-alpha*t)^2)
  return(den)
}

Ft_IG = function(t,alpha,beta){
  arg1 = (alpha*t-1)/sqrt(beta*t)
  arg2 = (-alpha*t-1)/sqrt(beta*t)
  Ft = pnorm(q=arg1,mean=0,sd=1) + (exp(2*alpha/beta)*pnorm(q=arg2,mean=0,sd=1))
  return(Ft)
}

St_IG = function(t,alpha,beta){
  arg1 = (alpha*t-1)/sqrt(beta*t)
  arg2 = (-alpha*t-1)/sqrt(beta*t)
  Ft = pnorm(q=arg1,mean=0,sd=1) + (exp(2*alpha/beta)*pnorm(q=arg2,mean=0,sd=1))
  return(1-Ft)
}


## log-likelihood function for Regression Inverse Gaussian

loglik_regig = function(par,time,delta,xa,xb){
  #parameters
  as = par[1:(dim(xa)[2]+1)] # b0, b1, b2
  bs = par[(dim(xa)[2]+2):(dim(xa)[2]+1+dim(xb)[2]+1)] # d0, d1, d2
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  
  den = as.vector(ft_IG(t=time,alpha=a,beta=b))
  surv = as.vector(St_IG(t=time,alpha=a,beta=b))
  
  log_lik = sum(delta*log(den) + (1-delta)*log(surv))
  return(-1*log_lik)
}






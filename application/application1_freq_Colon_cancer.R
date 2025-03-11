### -----------
### Estudo inicial do modelo frequentista
### de regressão MO Inverse Gaussian e MO Gompertz
### Data: Colon Cancer
### -----------

library(survival)

source('https://raw.githubusercontent.com/Dionisioneto/regression_defective_MO/refs/heads/main/application/source_reg_MO_Gompertz.R')
source('https://raw.githubusercontent.com/Dionisioneto/regression_defective_MO/refs/heads/main/application/source_reg_MO_IG.R')
source('https://raw.githubusercontent.com/Dionisioneto/regression_defective_MO/refs/heads/main/application/source_reg_MO_Gompertz.R')

## Funções de comparação de modelos (Inferência Clássica)

AICm = function(loglik,nparams){-2*loglik+2*nparams}
AICcm = function(loglik,nparams,n){(-2*loglik+2*nparams)+(2*nparams*(nparams+1))/(n-nparams-1)}
BICm = function(loglik,nparams,n){-2*loglik + nparams*log(n)}
HQICm = function(loglik,nparams,n){-2*loglik + (2*nparams*log(log(n)))}
CAICm = function(loglik,nparams,n){-2*loglik+nparams*(log(n)+1)}

colon = survival::colon

condition = (colon$etype==2) 
colon=colon[condition,] ## filtrando para o evento de interesse ser a morte.

head(colon)
dim(colon)

## Organizando o banco de dados

surv_time = colon$time/365 ## transforming from days to weeks
censoring = colon$status

# age = ifelse(colon$age>=70,1,0)
# cov = as.matrix(cbind(age))
covn = colon$node4
cov = as.matrix(cbind(covn))

## Kaplan-Meier suvival curve

kpfit_colon = survfit(Surv(surv_time, censoring) ~ covn)

plot(kpfit_colon , conf.int=F) # Existe a presença de fração de cura


## 1: Modelo Marshall-Olkin Gompertz

start = c(-0.1,-0.1,-0.1,0.1,1) # a0, a1, b0, b1, lambda

maxlike_regmog = optim(par = start,
                       fn = loglik_regmog,
                       gr = NULL,
                       hessian = T,
                       method = "BFGS",
                       time = surv_time,
                       delta = censoring,
                       xa=cov,
                       xb=cov)

maxlike_regmog$par


se_remog = sqrt(diag(solve(maxlike_regmog$hessian)))

maxlike_regmog$par - 1.96*se_remog
maxlike_regmog$par + 1.96*se_remog



AICm(loglik = -maxlike_regmog$value, nparams = length(maxlike_regmog$par))
AICcm(loglik = -maxlike_regmog$value, nparams = length(maxlike_regmog$par), n = length(surv_time))
BICm(loglik = -maxlike_regmog$value, nparams = length(maxlike_regmog$par), n = length(surv_time))
HQICm(loglik = -maxlike_regmog$value, nparams = length(maxlike_regmog$par), n = length(surv_time))
CAICm(loglik = -maxlike_regmog$value, nparams = length(maxlike_regmog$par), n = length(surv_time))


## 2: Modelo Marshall-Olkin Inverse Gaussian

start = c(-0.5,0.1,-0.1,-0.1,2) # a0, a1, b0, b1, lambda


maxlike_regmoig = optim(par = start,
                       fn = loglik_regmoig,
                       gr = NULL,
                       hessian = T,
                       method = "BFGS",
                       time = surv_time,
                       delta = censoring,
                       xa=cov,
                       xb=cov)

maxlike_regmoig$par

se_remoig = sqrt(diag(solve(maxlike_regmoig$hessian)))

maxlike_regmoig$par - 1.96*se_remoig
maxlike_regmoig$par + 1.96*se_remoig

AICm(loglik = -maxlike_regmoig$value, nparams = length(maxlike_regmoig$par))
AICcm(loglik = -maxlike_regmoig$value, nparams = length(maxlike_regmoig$par), n = length(surv_time))
BICm(loglik = -maxlike_regmoig$value, nparams = length(maxlike_regmoig$par), n = length(surv_time))
HQICm(loglik = -maxlike_regmoig$value, nparams = length(maxlike_regmoig$par), n = length(surv_time))
CAICm(loglik = -maxlike_regmoig$value, nparams = length(maxlike_regmoig$par), n = length(surv_time))



## 3: Modelo Defeituosos Clássicos (com estrutura regressora)

## 3.1: Gompertz regression

start = c(-0.5,0.1,-0.1,-0.1) # a0, a1, b0, b1

maxlike_reg_gpz = optim(par = start,
                        fn = loglik_reg_gpz,
                        gr = NULL,
                        hessian = T,
                        method = "BFGS",
                        time = surv_time,
                        delta = censoring,
                        xa=cov,
                        xb=cov)


maxlike_reg_gpz$par

AICm(loglik = -maxlike_reg_gpz$value, nparams = length(maxlike_reg_gpz$par))
AICcm(loglik = -maxlike_reg_gpz$value, nparams = length(maxlike_reg_gpz$par), n = length(surv_time))
BICm(loglik = -maxlike_reg_gpz$value, nparams = length(maxlike_reg_gpz$par), n = length(surv_time))
HQICm(loglik = -maxlike_reg_gpz$value, nparams = length(maxlike_reg_gpz$par), n = length(surv_time))
CAICm(loglik = -maxlike_reg_gpz$value, nparams = length(maxlike_reg_gpz$par), n = length(surv_time))


## 3.2: Gaussian inverse regression

start = c(-0.5,0.1,-0.1,-0.1) # a0, a1, b0, b1

maxlike_regig = optim(par = start,
                        fn = loglik_regig,
                        gr = NULL,
                        hessian = T,
                        method = "BFGS",
                        time = surv_time,
                        delta = censoring,
                        xa=cov,
                        xb=cov)


maxlike_regig$par

AICm(loglik = -maxlike_regig$value, nparams = length(maxlike_regig$par))
AICcm(loglik = -maxlike_regig$value, nparams = length(maxlike_regig$par), n = length(surv_time))
BICm(loglik = -maxlike_regig$value, nparams = length(maxlike_regig$par), n = length(surv_time))
HQICm(loglik = -maxlike_regig$value, nparams = length(maxlike_regig$par), n = length(surv_time))
CAICm(loglik = -maxlike_regig$value, nparams = length(maxlike_regig$par), n = length(surv_time))


## O modelo de fração de cura weibull
library(rstpm2)
library(cuRe)

fit.wei = fit.cure.model(Surv(time/365, status) ~ node4,
                          #formula.surv = list(~ 1, ~ 1),  # Sem covariáveis na sobrevivência
                          data = colon,
                          bhazard = NULL,
                          type = "mixture",
                          dist = "weibull",
                          link = "logit")


summary(fit.wei)

fit.wei$coefs

fit.wei$optim$value
length(fit.wei$optim$par)

AICm(loglik = -fit.wei$optim$value, nparams = length(fit.wei$optim$par))
AICcm(loglik = -fit.wei$optim$value, nparams = length(fit.wei$optim$par), n = length(surv_time))
BICm(loglik = -fit.wei$optim$value, nparams = length(fit.wei$optim$par), n = length(surv_time))
HQICm(loglik = -fit.wei$optim$value, nparams = length(fit.wei$optim$par), n = length(surv_time))
CAICm(loglik = -fit.wei$optim$value, nparams = length(fit.wei$optim$par), n = length(surv_time))


## Na Abordagem Clássica o Modelo Marshall-Olkin Gompertz ganhou!  

# Para os modelos Marshall-Olkin Gompertz e Marshall-Olkin Inverse
# as covariáveis foram significativas.

## ---------
## Gráficos da curva de sobrevivência
## ---------
kpfit_colon = survfit(Surv(surv_time, censoring) ~ 1)

plot(kpfit_colon , conf.int=F)




## ---------
## Gráficos da curva de sobrevivência
## ---------

## Gráficos
# Ajuste do modelo base [Gompertz]

as_gptz = maxlike_reg_gpz$par[1:2]
bs_gptz = maxlike_reg_gpz$par[3:4]

t0=seq(0,10,length=length(aest[covn==0]))
t1=seq(0,10,length=length(aest[covn==1]))

aest_gptz = cbind(1,cov) %*% as_gptz
best_gptz = exp(cbind(1,cov) %*% bs_gptz)

st0_gptz = St_Gompertz(t=t0,alpha=aest_gptz[covn==0],beta=best_gptz[covn==0])
st1_gptz = St_Gompertz(t=t1,alpha=aest_gptz[covn==1],beta=best_gptz[covn==1])


## ajuste do modelo [Marshall-Olkin Gompertz]
kpfit_colon = survfit(Surv(surv_time, censoring) ~ covn)

plot(kpfit_colon, conf.int=F, bty = "n", xlim=range(surv_time)+c(0,1),
     xlab="Time (Years)", ylab = "Survival function",lwd=1.5, lty = c(6, 1))

lines(t0,st0_gptz,type='l', col = "red", lwd=2, lty = 6)
lines(t1,st1_gptz,type='l', col = "red", lwd=2, lty = 1)



# Marshall-Olkin Gompertz 
as = maxlike_regmog$par[1:2]
bs = maxlike_regmog$par[3:4]
ls = maxlike_regmog$par[5]

aest = cbind(1,cov) %*% as
best = exp(cbind(1,cov) %*% bs)

t0=seq(0,10,length=length(aest[covn==0]))
t1=seq(0,10,length=length(aest[covn==1]))

st0 = Stmo_gompertz(t=t0,alpha=aest[covn==0],beta=best[covn==0],lambda=ls)
st1 = Stmo_gompertz(t=t1,alpha=aest[covn==1],beta=best[covn==1],lambda=ls)


lines(t0,st0,type='l', col = "steelblue", lwd=2, lty = 6)
lines(t1,st1,type='l', col = "steelblue", lwd=2, lty = 1)

legend("topright", legend = c("<=4 Lymph nodes",">4 Lymph nodes","Gompertz", "Marshall-Olkin Gompertz"),
       col=c("black","black","red", "steelblue"), lty = c(6,1,1,1), cex=1, box.lty=0)



# Ajuste do modelo base [Inverse Gaussian]

as_ig = maxlike_regig$par[1:2]
bs_ig = maxlike_regig$par[3:4]

aest_ig = cbind(1,cov) %*% as_ig
best_ig = exp(cbind(1,cov) %*% bs_ig)

st0_ig = St_IG(t=t0,alpha=aest_ig[covn==0],beta=best_ig[covn==0])
st1_ig = St_IG(t=t1,alpha=aest_ig[covn==1],beta=best_ig[covn==1])

plot(kpfit_colon, conf.int=F, bty = "n", xlim=range(surv_time)+c(0,1),
     xlab="Time (Years)", ylab = "Survival function",lwd=1.5, lty = c(6, 1))

lines(t0,st0_ig,type='l', col = "red", lwd=2, lty = 6)
lines(t1,st1_ig,type='l', col = "red", lwd=2, lty = 1)


# Ajuste do modelo base [Marshall-Olkin Inverse Gaussian]

as_moig = maxlike_regmoig$par[1:2]
bs_moig = maxlike_regmoig$par[3:4]
ls_moig = maxlike_regmoig$par[5]

aest_moig = cbind(1,cov) %*% as_moig
best_moig = exp(cbind(1,cov) %*% bs_moig)

st0_moig = Stmo_IG(t=t0,alpha=aest_moig[covn==0],beta=best_moig[covn==0],lambda=ls_moig)
st1_moig = Stmo_IG(t=t1,alpha=aest_moig[covn==1],beta=best_moig[covn==1],lambda=ls_moig)

lines(t0,st0_moig,type='l', col = "steelblue", lwd=2, lty = 6)
lines(t1,st1_moig,type='l', col = "steelblue", lwd=2, lty = 1)


legend("topright", 
       legend = c("<=4 Lymph nodes",">4 Lymph nodes","Inverse Gaussian", "Marshall-Olkin Inv. Gauss."),
       col=c("black","black","red", "steelblue"), lty = c(6,1,1,1), cex=1, box.lty=0)





## ---------
## Analise residual do modelo 
## Marshall-Olkin Gompertz
## ---------

## 1. Martingales residuals


as = maxlike_regmog$par[1:2]
bs = maxlike_regmog$par[3:4]
ls = maxlike_regmog$par[5]

aest = cbind(1,cov) %*% as
best = cbind(1,cov) %*% bs

rm_mog = censoring + log(Stmo_gompertz(t=surv_time,alpha = aest,beta=best,lambda=ls))
#r_mog = censoring*(-log(Stmo_gompertz(t=surv_time,alpha = aest,beta=best,lambda=ls))) + (1-censoring) * (1-log(Stmo_gompertz(t=surv_time,alpha = aest,beta=best,lambda=ls)))
rd_mog = sign(rm_mog)*(-2*(rm_mog + censoring*log(censoring-rm_mog)))^(1/2)

## 2. Deviance residuals

rd_mog = sign(rm_mog)*(-2*(rm_mog + censoring*log(censoring-rm_mog)))^(1/2)

par(mfrow=c(1,2))
plot(1:length(rm_mog),rm_mog, col = ifelse(censoring == 1, "black", "red"),
     xlab = "Index", ylab = "Martingale residuals")

plot(1:length(rd_mog),rd_mog, ylim=c(-4,4), col = ifelse(censoring == 1, "black", "red"),
     xlab = "Index", ylab = "Deviance residuals")

abline(h=-3,lwd=2,col='red', lty=6)
abline(h=3,lwd=2,col='red', lty=6)


## ---
## Global diagnosis
## Modelo Marshall-Olkin Gompertz
## ---

## Distancia genreralizada de Cook

theta_hat = maxlike_regmog$par
mfisher_hat = maxlike_regmog$hessian
log_vero = -maxlike_regmog$value

n_amostral = length(surv_time)
gd = rep(0,n_amostral)
ldist = rep(0,n_amostral)
  
start = c(-0.1,-0.1,-0.1,0.1,1) # a0, a1, b0, b1, lambda

for(i in 1:n_amostral){
  
  maxlike_regmogi = optim(par = start,
                         fn = loglik_regmog,
                         gr = NULL,
                         hessian = T,
                         method = "BFGS",
                         time = surv_time[-i],
                         delta = censoring[-i],
                         xa=as.matrix(cov[-i,]),
                         xb=as.matrix(cov[-i,]))
  
  esti = maxlike_regmogi$par
  gd[i] = t(theta_hat-esti) %*%  mfisher_hat %*%  (theta_hat-esti)
  ldist[i] = 2*(log_vero-(-maxlike_regmogi$value))
}

#par(mfrow=c(1,2))

## Distancia de Cook
mod_cook = abs(gd)

plot(mod_cook, type = "h", lwd = 2, ylim = c(0, max(mod_cook) + 15000),
     xlim=c(0,length(mod_cook)+100),
     ylab = "|Generalized Cook Distance|",  bty = "n")

above_threshold = which(mod_cook > 2000)

if (length(above_threshold) > 0) {
  for (i in seq_along(above_threshold)) {
    text(above_threshold[i], mod_cook[above_threshold[i]] + 5000 + (i * 1000),  # Incremento para separar os rótulos
         labels = above_threshold[i], pos = 1, offset = -0.5, col = "red", cex = 0.8)
  }
}


## Distancia de verossimilhancas
mod_ldist = abs(ldist)

plot(mod_ldist, type = "h", lwd = 2, ylab = "|Likelihood Distance|",
     ylim = c(0,max(mod_ldist)+20), xlim=c(0,length(mod_ldist)+100),  bty = "n")

above_th2d = which(mod_ldist > 20)

if (length(above_th2d) > 0) {
  for (i in seq_along(above_th2d)) {
    text(above_th2d[i], mod_ldist[above_th2d[i]] + 5 + (i * 2),  # Incremento de 2 unidades para cada rótulo
         labels = above_th2d[i], pos = 1, col = "red", cex = 0.8)
  }
}


## Cálculos de medidas
## relative change

rc_est = function(estset,est){abs((est-estset)/est)*100}
rcse_est = function(seset,se){abs((se-seset)/se)*100}


sets = list(c(14), c(21), c(116), c(139), c(172), c(228),
             c(235), c(14, 21, 116, 172, 228, 235),
             c(14, 21, 116, 139, 172, 228, 235))

matriz_rc_summary = array(data = NA, 
                           dim = c(2, length(maxlike_regmog$par) ,length(sets)),
                           dimnames = list(c("RC", "RCse"), c("a0", "a1", "b0", "b1", "lambda"), 
                                           paste0("Set", seq_along(sets))))

matriz_rc_ic = array(data = NA, 
                          dim = c(length(maxlike_regmog$par), 2, length(sets)),
                          dimnames = list(c("a0", "a1", "b0", "b1", "lambda"), c("inf","sup"), 
                                          paste0("Set", seq_along(sets))))


start = c(-0.1,-0.1,-0.1,0.1,1) # a0, a1, b0, b1, lambda

for(j in 1:length(sets)){
  maxlike_regmog_set = optim(par = start,
                               fn = loglik_regmog,
                               gr = NULL,
                               hessian = T,
                               method = "BFGS",
                               time = surv_time[-sets[j][[1]]],
                               delta = censoring[-sets[j][[1]]],
                               xa=as.matrix(cov[-sets[j][[1]],]),
                               xb=as.matrix(cov[-sets[j][[1]],]))
  
  est_set = maxlike_regmog_set$par
  se_set = sqrt(diag(solve(maxlike_regmog_set$hessian)))
  
  matriz_rc_summary[1,,j]=round(rc_est(estset=est_set,est=maxlike_regmog$par),4)
  matriz_rc_summary[2,,j]=round(rcse_est(seset=se_set,se=se_remog),4)
  
  matriz_rc_ic[,1,j] = round(est_set - 1.96*se_set,4)
  matriz_rc_ic[,2,j] = round(est_set + 1.96*se_set,4)
  
}


head(matriz_rc_summary)
head(matriz_rc_ic)

# em um CSV
write.csv2(as.data.frame(matriz_rc_summary), "mog_matriz_rc_summary.csv", row.names = FALSE)
write.csv2(as.data.frame(matriz_rc_ic), "mog_matriz_rc_ic.csv", row.names = FALSE)


## Teste da razão de verossimilhanças


## H0: O modelo restrito é preferível do que o geral;
## H1: O modelo geral é preferível do que o modelo restrito.

## MO-Gompertz versus Gompertz
lrv_g = -2*(-maxlike_reg_gpz$value-(-maxlike_regmog$value))

1-pchisq(q=lrv_g,df = 1)

## MO-IG versus inverse Gaussian
lrv_ig = -2*(-maxlike_regig$value-(-maxlike_regmoig$value))

1-pchisq(q=lrv_ig,df = 1)


qchisq(p=1-0.05,df=1)





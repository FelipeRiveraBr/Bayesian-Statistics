##PARCIAL 3 BAYESIANA
Y<-c(1,3,2,12,1,1)
X<-c(33,12,27,90,12,17)


### PUNTO 3


###  MODELO 1  ###
library(coda)

#Datos

x=c(33,14,27,90,12,17)
y=c(1,3,2,12,1,1)

#Hiperparametros

adealpha=1
bdealpha=1
adebeta=10
bdebeta=1

#Valores iniciales

set.seed(6)
nsamples=50000
alpha=rgamma(1,adealpha,bdealpha)
beta=rgamma(1,adebeta,bdebeta)
theta=rgamma(6,alpha,beta)
theta=t(theta)
delta=1.5
estimadores=matrix(ncol=8,nrow=nsamples)
vero=c()
inicial=0

# Cadena

for(n in 1:nsamples){
  for(i in 1:6){
    theta[i]=rgamma(1,alpha+y[i],beta+x[i])
  }
  beta=rgamma(1,6*alpha+adebeta,sum(theta)+bdebeta)
  alphapropuesta=abs(runif(1,alpha-delta,alpha+delta))
  tasa=((beta^(6*(alphapropuesta-alpha)))*(prod(theta)^(alphapropuesta-alpha))*((alphapropuesta/alpha)^(adealpha-1))*exp(-bdealpha*(alphapropuesta-alpha))*(gamma(alpha)^6))/(gamma(alphapropuesta)^6)
  if(runif(1) < tasa){
    alpha=alphapropuesta
    inicial=inicial+1
  }
  estimadores[n,]=c(theta,beta,alpha)
  vero[n]=sum(dpois(y,theta*x,log = T))
}
inicial/nsamples
step=seq(10,50000,by=10)
estimadores=estimadores[step,]
vero=vero[step]

colMeans(estimadores)


par(mfrow=c(3,2))
plot(estimadores[,1],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[1]),
     main=expression(paste(theta[1]," Modelo 1")))
plot(estimadores[,2],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[2]),
     main=expression(paste(theta[2]," Modelo 1")))
plot(estimadores[,3],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[3]),
     main=expression(paste(theta[3]," Modelo 1")))
plot(estimadores[,4],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[4]),
     main=expression(paste(theta[4]," Modelo 1")))
plot(estimadores[,5],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[5]),
     main=expression(paste(theta[5]," Modelo 1")))
plot(estimadores[,6],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[6]),
     main=expression(paste(theta[6]," Modelo 1")))

LLPROM1<-mean(vero)

tamanoefectivo_1=apply(X=estimadores, MARGIN = 2, FUN = effectiveSize)
tamanoefectivo_1
TEPROM1<-mean(tamanoefectivo_1)

desviaciones_1=apply(X=estimadores, MARGIN = 2, FUN = sd)
montecarloerrors_1=desviaciones_1/sqrt(tamanoefectivo_1)
montecarloerrors_1
EMC1<-mean(montecarloerrors_1)

###  MODELO 2  ###


B<-50000
n<-length(Y)

#hiperparametros
ALPHA_alpha <- 1
BETA_alpha <- 1
ALPHA_beta <- 10
BETA_beta <- 1
ALPHA_lambda <- 0
BETA_lambda <- 10



#Valores Iniciales
thetai<- rep(0.5,6)
alpha <-1
beta<-1
lambda  <-1

#logaritmos
logtheta <-log(thetai)
logalpha<-log(alpha)
logbeta<-log(beta)
loglambda<-log(lambda)

#parametros de ajuste
delta_theta <-2*c(1,1,1,1,1,1)
delta_alpha <-1
delta_lambda <-3



logtheta_star<-NULL
logalpha_star<-NULL
loglambda_star<-NULL
THETA <- matrix(NA,ncol =9,nrow = B)
LL<- NULL
accept_theta   <- rep(0,6)
accept_alpha   <- 0
accept_lambda   <- 0



#MCMC
set.seed(1234)

for(b in 1:B) {

#Muestrear Beta
alpha_b <- n*exp(logalpha) + ALPHA_beta
beta_b  <- sum(exp(logtheta)) + BETA_beta
beta    <- rgamma(1, alpha_b, beta_b)
  
#Muestrear theta por algoritmo de metropolis
  for (i in 1:n) {
    logtheta_star[i]<-runif(1,logtheta[i]-delta_theta[i],
                            logtheta[i]+delta_theta[i])
    
    #tasa de aceptacion
    lhr.t<-dnbinom(Y[i],(exp(logtheta_star[i])*X[i])/exp(loglambda) ,
                   1/(1+exp(loglambda)) , log = T) -
      dnbinom(Y[i],(exp(logtheta[i])*X[i])/exp(loglambda) ,1/(1+exp(loglambda)) 
              ,log = T) +
      dgamma(exp(logtheta_star[i]), exp(logalpha), beta, log = T) -
      dgamma(exp(logtheta[i]), exp(logalpha), beta, log = T)+
      logtheta_star[i]-logtheta[i]
    

    if (log(runif(1)) < lhr.t) {
      logtheta[i] <- logtheta_star[i]
      accept_theta[i]   <- accept_theta[i] + 1 
    }
  }
  

#Muestrear alpha por algoritmo de metropolis
  logalpha_star <-runif(1,logalpha-delta_alpha,logalpha+delta_alpha)
  
  # tasa de aceptacion
  lhr.a <- sum(dgamma(exp(logtheta), exp(logalpha_star), beta, log = T)) -
    sum(dgamma(exp(logtheta), exp(logalpha), beta, log = T)) +
    dgamma(exp(logalpha_star), ALPHA_alpha, BETA_alpha, log = T) -
    dgamma(exp(logalpha), ALPHA_alpha, BETA_alpha, log = T)+
    logalpha_star-logalpha
  
  
  if (log(runif(1)) < lhr.a) {
    logalpha <- logalpha_star
    accept_alpha   <- accept_alpha + 1
  }

  
  
#Muestrear lambda por algoritmo de metropolis
  loglambda_star <-runif(1,loglambda-delta_lambda,loglambda+delta_lambda)
  
  # tasa de aceptacion
  lhr.l <- sum(dnbinom(Y,(exp(logtheta)*X)/exp(loglambda_star) ,
                       1/(1+exp(loglambda_star)) , log = T)) -
    sum(dnbinom(Y, (exp(logtheta)*X)/exp(loglambda) ,1/(1+exp(loglambda))  , 
                log = T)) +
    dunif(exp(loglambda_star), ALPHA_lambda, BETA_lambda) -
    dunif(exp(loglambda), ALPHA_lambda, BETA_lambda)+
    loglambda_star-loglambda
 
  
  if (log(runif(1)) < lhr.l) { 
    loglambda <- loglambda_star 
    accept_lambda   <- accept_lambda + 1 
  }
  
  #almacenar los parámetros
  THETA[b,]<- c(logtheta,beta,logalpha,loglambda)

  #LL
  LL[b]<-sum(dnbinom(Y,(exp(logtheta)*X)/exp(loglambda),1/(1+exp(loglambda)), 
                      log = T))
}




#Tasas de aceptación
accept_lambda/B
accept_theta/B
accept_alpha/B

#transformar los parametros
THETA[,c(1,2,3,4,5,6,8,9)]<-exp(THETA[,c(1,2,3,4,5,6,8,9)])

#Muestreo aleatorio sobre las cadenas
set.seed(1)
step<-sample(1:20000,replace = F)
THETA<-THETA[step,]
LL<-LL[step]


#log verosimilitud
par(mfrow=c(2,1))
plot(vero,type = "p", pch = ".", main="Log-verosimilitud Modelo 1", 
     cex.axis = 0.8 , xlab = "Iteración",ylab = "Log-verosimilitud")
plot(LL,type = "p", pch = ".", main="Log-verosimilitud  Modelo 2", 
     cex.axis = 0.8 , xlab = "Iteración",ylab = "Log-verosimilitud")
par(mfrow=c(3,2))
plot(THETA[,1],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[1]),
     main=expression(paste(theta[1],"  Modelo 2")))
plot(THETA[,2],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[2]),
     main=expression(paste(theta[2],"  Modelo 2")))
plot(THETA[,3],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[3]),
     main=expression(paste(theta[3],"  Modelo 2")))
plot(THETA[,4],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[4]),
     main=expression(paste(theta[4],"  Modelo 2")))
plot(THETA[,5],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[5]),
     main=expression(paste(theta[5],"  Modelo 2")))
plot(THETA[,6],type = "p", pch = ".", cex.axis = 0.8 , 
     xlab = "Iteración", ylab = expression(theta[6]),
     main=expression(paste(theta[6],"  Modelo 2")))


LLPROM2<-mean(LL)

tamanoefectivo_2=apply(X=THETA, MARGIN = 2, FUN = effectiveSize)
tamanoefectivo_2
TEPROM2<-mean(tamanoefectivo_2)

desviaciones_2=apply(X=THETA, MARGIN = 2, FUN = sd)
montecarloerrors_2=desviaciones_2/sqrt(tamanoefectivo_2)
montecarloerrors_2
EMC2<-mean(montecarloerrors_2)




## PUNTO 4

##Graficas posteriores
prevobs<-Y/X
prevobs

par(mfrow=c(3,2))

#Theta 1
plot(density(THETA[,1]),ylim=c(0,15),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[1])),xlab=expression(theta[1]),
     ylab=expression(paste("p(",theta[1],"|y)")))
lines(density(estimadores[,1]),col="orange2", lwd=2)
abline(v=prevobs[1],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)

#Theta 2
plot(density(THETA[,2]),ylim=c(0,6),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[2])),xlab=expression(theta[2]),
     ylab=expression(paste("p(",theta[2],"|y)")))
lines(density(estimadores[,2]),col="orange2", lwd=2)
abline(v=prevobs[2],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)

#Theta 3
plot(density(THETA[,3]),ylim=c(0,10),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[3])),xlab=expression(theta[3]),
     ylab=expression(paste("p(",theta[3],"|y)")))
lines(density(estimadores[,3]),col="orange2", lwd=2)
abline(v=prevobs[3],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)

#Theta 4
plot(density(THETA[,4]),ylim=c(0,12),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[4])),xlab=expression(theta[4]),
     ylab=expression(paste("p(",theta[4],"|y)")))
lines(density(estimadores[,4]),col="orange2", lwd=2)
abline(v=prevobs[4],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)

#Theta 5
plot(density(THETA[,5]),ylim=c(0,8),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[5])),xlab=expression(theta[5]),
     ylab=expression(paste("p(",theta[5],"|y)")))
lines(density(estimadores[,5]),col="orange2", lwd=2)
abline(v=prevobs[5],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)

#Theta 6
plot(density(THETA[,6]),ylim=c(0,9),col="blue4", lwd=2,
     main=expression(paste("distribución posterior ",theta[6])),xlab=expression(theta[6]),
     ylab=expression(paste("p(",theta[6],"|y)")))
lines(density(estimadores[,6]),col="orange2", lwd=2)
abline(v=prevobs[6],lty=2,col="red3")
legend("topright",legend = c("Modelo 1","Modelo 2","Prev. Observada"), 
       col=c("blue4", "orange2","red3"), lty=c(1,1,2),lwd = c(2,2,1),
       xpd = TRUE,seg.len = 3,cex=1)



## PUNTO 5

par(mfrow=c(2,3),mar=c(4, 4, 4, 4))
for (i in 1:6) {
  if (i==2) {
  }else{
    plot(estimadores[,i],estimadores[,2],xlim = c(0,0.7),
         main = bquote(theta[2]~vs ~theta[.(i)]),xlab=bquote(~theta[.(i)]),
         ylab=bquote(~theta[.(2)]),cex.lab=1.1)
    abline(0,1,lty=1,lwd=1, col= "orange1")
    pr<-round(mean(estimadores[,2]>estimadores[,i]),1)
    legend("topright",legend=c(expression(paste("Pr(",theta[2],">",theta[i],")",
                                            " ")),pr))
    legend("right",legend=expression(paste(theta[2],"=",theta[i])),
           col="orange2",lty=1)
  }
}

par(mfrow=c(2,3),mar=c(4, 4, 4, 4))
for (i in 1:6) {
  if (i==2) {
  }else{
    plot(THETA[,i],THETA[,2],xlim = c(0,2),
         main = bquote(theta[2]~vs ~theta[.(i)]),xlab=bquote(~theta[.(i)]),
         ylab=bquote(~theta[.(2)]),cex.lab=1.1)
    abline(0,1,lty=1,lwd=1, col= "orange1")
    pr<-round(mean(THETA[,2]>THETA[,i]),1)
    legend("topright",legend=c(expression(paste("Pr(",theta[2],">",theta[i],")",
                                                " ")),pr))
    legend("right",legend=expression(paste(theta[2],"=",theta[i])),
           col="orange2",lty=1)
    }
}

##PUNTO 6

maximo=c()
for(i in 1:dim(estimadores)[1]){
  maximo[i]<-max(estimadores[i,1:6])
}
mean(estimadores[,2]==maximo)

maximo=c()
for(i in 1:dim(THETA)[1]){
  maximo[i]<-max(THETA[i,1:6])
}
mean(THETA[,2]==maximo)


## PUNTO 7

##DIC
#Modelo 1
theta_hat  <- as.vector(colMeans(estimadores[,1:6]))
lpyth_m1   <- sum(dpois(y,theta_hat*x, log = T))
pDIC_m1    <- 2*(lpyth_m1 - mean(vero))
dic_m1     <- -2*lpyth_m1 + 2*pDIC_m1
dic_m1

#Modelo 2
t2_hat <- as.vector(colMeans(THETA[,1:6]))
l2_hat <- mean(THETA[,9])
lpyth_m2 <- sum(dnbinom(Y, t(t2_hat*x/l2_hat),1/(1+l2_hat), log = T))
pDIC_m2 <- 2*(lpyth_m2 - mean(LL))
dic_m2 <- -2*lpyth_m2 + 2*pDIC_m2
dic_m2




#WAIC
# WAIC M1
lppd_m1  <- 0
pWAIC_m1 <- 0
for (i in 1:n) {
  # lppd
  tmp1    <- dpois(x = y[i], lambda = estimadores[,i]*x[i])
  lppd_m1 <- lppd_m1 + log(mean(tmp1))
  # pWAIC
  tmp2 <- dpois(x = y[i], lambda = estimadores[,i]*x[i], log = T)
  pWAIC_m1 <- pWAIC_m1 + 2*(log(mean(tmp1)) - mean(tmp2))
}

waic_m1 <- -2*lppd_m1 + 2*pWAIC_m1
waic_m1

# WAIC M2
lppd_m2  <- 0
pWAIC_m2 <- 0
for (i in 1:n) {
  # lppd
  tmp1    <- dnbinom(Y[i], X[i]*THETA[,i]/THETA[,9],1/(1+THETA[,9]))
  lppd_m2 <- lppd_m2 + log(mean(tmp1))
  # pWAIC
  tmp2 <-  dnbinom(Y[i], X[i]*THETA[,i]/THETA[,9],1/(1+THETA[,9]),log=T)
  pWAIC_m2 <- pWAIC_m2 + 2*(log(mean(tmp1)) - mean(tmp2))
}

waic_m2 <- -2*lppd_m2 + 2*pWAIC_m2
waic_m2



# BIC
# BIC M1
k_m1 <- 8
bic_m1 <- -2*lpyth_m1 + k_m1*log(n)
bic_m1
# BIC M2
k_m2 <- 9
bic_m2 <- -2*lpyth_m2 + k_m2*log(n)
bic_m2


#AIC
#AIC M1
maxv_1<-max(vero)
AIC1<-2*k_m1-2*maxv_1
#AIC M2
maxv_2<-max(LL)
AIC2<-2*k_m2-2*maxv_2



## ppp

# ppp M1

M <-5000
ts_obs<-c(mean(y),sd(y))
TS1<-matrix(NA,ncol = 2,nrow = 5000)
yp <-matrix(NA,ncol = 1,nrow = 6)
set.seed(123)
for (m in 1:M) {
  for (i in 1:6) {
    yp[i]<-rpois(1,estimadores[m,i]*x[i])    
  }
  TS1[m,]<-c(mean(yp),sd(yp))
  
}

pppmean<-mean(TS1[,1]>ts_obs[1])
pppSD<-mean(TS1[,2]>ts_obs[2])

par(mfrow=c(1,2))
den_m1 <- density(TS1[,1], adjust = 2)
plot(den_m1,    
     xlab = expression(paste('t*=',bar(bold(y))[b],sep='')),
     ylab = expression(paste('P(t*|',bold(y),')')),
     main ="ppp Media Modelo 1"
)
abline(v = ts_obs[1], col = "blue4",lty=2 ,lwd = 2)
legend("topright", legend =expression(bold(bar(y)[obs])) , col = "blue4", 
       lwd = 2, cex=0.8, lty=2)
legend("right", paste("ppp = ",pppmean), adj = 0.2,cex = 0.8)


den_sd1 <- density(TS1[,2], adjust = 1.5)
plot(den_sd1,   
     xlab = expression(paste('t*=',bold(sd(y[b])),sep='')),
     ylab = expression(paste('P(t*|',bold(y),')')),
     main ="ppp SD Modelo 1"
)
abline(v = ts_obs[2], col = "orange2", lwd = 2,lty=2)
legend("topright", legend = expression(bold(sd(y[obs]))),col = "orange2", 
       lwd = 2, cex=0.8, lty=2)
legend("right", paste("ppp = ",pppSD), adj = 0.2,cex = 0.8)



# ppp M2
M <-5000
ts_obs<-c(mean(y),sd(y))
TS1<-matrix(NA,ncol = 2,nrow = 5000)
yp <-matrix(NA,ncol = 1,nrow = 6)
set.seed(123)
for (m in 1:M) {
  for (i in 1:6) {
    yp[i]<-rpois(1,THETA[m,i]*x[i])    
  }
  TS1[m,]<-c(mean(yp),sd(yp))
  
}

pppmean<-mean(TS1[,1]>ts_obs[1])
pppSD<-mean(TS1[,2]>ts_obs[2])

par(mfrow=c(1,2))
den_m1 <- density(TS1[,1], adjust = 2)
plot(den_m1,    
     xlab = expression(paste('t*=',bar(bold(y))[b],sep='')),
     ylab = expression(paste('P(t*|',bold(y),')')),
     main ="ppp Media Modelo 2"
)
abline(v = ts_obs[1], col = "blue4",lty=2 ,lwd = 2)
legend("topright", legend =expression(bold(bar(y)[obs])) , col = "blue4", 
       lwd = 2, cex=0.8, lty=2)
legend("right", paste("ppp = ",pppmean), adj = 0.2,cex = 0.8)


den_sd1 <- density(TS1[,2], adjust = 1.5)
plot(den_sd1,   
     xlab = expression(paste('t*=',bold(sd(y[b])),sep='')),
     ylab = expression(paste('P(t*|',bold(y),')')),
     main ="ppp SD Modelo 2"
)
abline(v = ts_obs[2], col = "orange2", lwd = 2,lty=2)
legend("topright", legend = expression(bold(sd(y[obs]))),col = "orange2", 
       lwd = 2, cex=0.8, lty=2)
legend("right", paste("ppp = ",pppSD), adj = 0.2,cex = 0.8)














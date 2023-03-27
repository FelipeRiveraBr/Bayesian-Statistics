#PARCIAL 1
set.seed(123)
#datos
Y<-c(495, 541, 1461, 1555 ,1603, 2201, 2750, 3468, 3516, 4319, 
     6622, 7728, 13159 ,21194)
Y #vector de información
n<-length(Y)
s<- sum(Y) 
s #estadístico suficiente

##PUNTO 4
#Distribución de lambda
#hiperparámetros previa
Mu0<-4500
sig0<-1800
alpha1<-((Mu0/sig0)^2)+2
beta1<-((alpha1-1)*Mu0)
alpha1
beta1
#hiperparámetros posterior
alpha2<-alpha1+n
beta2<-beta1+s
#función de densidad Gamma Inversa
dGI<-function (x, alpha, beta) {
  dens<- alpha* log (beta)-lgamma(alpha)-(alpha+1)*log(x)-beta/x
  return (exp(dens))
}

lambdapre<-seq(0,Mu0+5*sig0,length=100001)
lambdapre<-lambdapre[-1]
plambdapre<-dGI(lambdapre,alpha1,beta1)
prelambda<-cbind(lambdapre,plambdapre)

lambdapost<-seq(0,Mu0+5*sig0,length=100001)
lambdapost<-lambdapost[-1]
plambdapost<-dGI(lambdapost,alpha2,beta2)
postlambda<-cbind(lambdapost,plambdapost)


min(lambdapre)<min(lambdapost)
max(lambdapre)<max(lambdapost)

max(plambdapre)<max(plambdapost)

plot(NA,NA,xlim=c(min(lambdapre),max(lambdapre)),ylim=c(0,max(plambdapost)),
     xlab=expression ( lambda ),
     ylab=expression(paste('p(',lambda,'|y)',' , p(',lambda,')')),
     cex.lab=1.1)

lines(prelambda,col='royalblue4',lwd=2)
lines(postlambda,col='goldenrod1',lwd=2)
legend(x = "topright", inset = 0.05,
     legend = c(expression(paste('p(',lambda,')')),
     expression(paste('p(',lambda,'|y)'))), 
     lty = c(1, 1),bty = "n",
     col = c("royalblue4", "goldenrod1"),
     lwd=c(2,2))

##PUNTO 6
##BAYESIANO
#Simulación
#Simulamos de una gamma
set.seed(123)
invlambdasim<-rgamma(100000,alpha2,beta2)
lambdasim<-1/invlambdasim
Mupostsim<-mean(lambdasim)
Mupostsim
varpostsim<-var(lambdasim)
sigmapostsim<-sqrt(varpostsim)
CVpost<-sigmapostsim/Mupostsim
CVpost
ls<-quantile(lambdasim, 0.975)
li<-quantile(lambdasim, 0.025)
intbay<-c(li,ls)
intbay

##Frecuentista
Ybar<-mean(Y)
Ybar
IF<-n/(mean(Y)^2)  
varfrec<-1/IF  
sdfrec<-sqrt(varfrec)
CVfreq<-sdfrec/Ybar
CVfreq
intfreq<-qnorm(c(.025,.975),Ybar,sdfrec)
intfreq

#Bootstrap
set.seed(123)
out <- NULL
for (i in 1:100000) {
  yy <- sample(x = Y, size = length(Y), replace = T)
  out[i] <- mean(yy)
}
est_boot <- mean(out)
est_boot
CV_boot <- sd(out)/mean(out)
CV_boot
ic_boot  <- quantile(out, probs = c(.025, .975))
ic_boot

##PUNTO7
#Pr(lambda<4000|Y) y Pr(y*<4000|Y)
#Para Pr(lambda<4000|Y) utilizamos la simulación obtenida
mean(lambdasim<4000)
#Para Pr(y*<4000|Y) debemos aproximar la distribución
#predictiva posterior
#simulamos Y de la dist muestral usando los lambda
#obtenidos de la dist posterior como parámetros
set.seed(123)
ystar<- rexp(100000,rate=1/lambdasim)
mean(ystar<4000)
mean(lambdasim==4000)

##PUNTO8
##Factor de Bayes
a<-alpha1
b<-beta1
lambda0<-4000
B10<-exp(a*log(b)+lgamma(a+n)-lgamma(a)-(a+n)*log(b+s)+n*log(lambda0)
         +(s/lambda0))
B10

##PUNTO 9
#Farctor de Bayes 
s1<-s
n1<-n
Y1<-Y
Y2<-c(294, 569, 766, 1576, 1602, 2015, 2166, 3885, 8141, 10285)
s2<-sum(Y2)
n2<-length(Y2)
B10_2<-exp(a*log(b)+lgamma(a+n1)+lgamma(a+n2)+(a+n1+n2)*log(b+s1+s2)
           -lgamma(a)-(a+n1)*log(b+s1)-(a+n2)*log(b+s2)-lgamma(a+n1+n2))
B10_2


##PUNTO 10
#Bondad de ajuste
#Tipo 1
tobs_1<-s/n
t_1<-NULL
set.seed(123)
for (i in 1:100000) {
  t_1[i]<-mean(rexp(100000,rate=1/lambdasim[i]))
}
hist(x=t_1, freq=F, col="gray90", border="gray90", xlab="t_1",
     ylab="p(t_1 | y)", main="Bondad de ajuste Alambre tipo 1",
     ylim = c(0,0.000405))
lines(density(t_1), col = "goldenrod1", lwd = 2)

abline(v = tobs_1, col = "royalblue4", lwd = 2, lty = 1)
abline(v = quantile(x = t_1, probs = c(0.025, 0.975)), lty = 2, 
       lwd = 1, col = "gray")
legend("topright", legend = c("Posterior", "IC 95%", "t obs"), 
       col = c("goldenrod1", "gray", "royalblue4"),
       lty = c(1,2,1), lwd = 1, bty = "n")
#pr(t>tobs|y)
ppp_1<-mean(t_1>tobs_1)



#Tipo 2

alpha2_2<-alpha1+n2
beta2_2<-beta1+s2
set.seed(123)
invlambdasim_2<-rgamma(100000,alpha2_2,beta2_2)
lambdasim_2<-1/invlambdasim_2
tobs_2<-s2/n2
t_2<-NULL
set.seed(123)
for (i in 1:100000) {
  t_2[i]<-mean(rexp(100000,rate=1/lambdasim_2[i]))
}
hist(x=t_2, freq=F, col="gray90", border="gray90", xlab="t_2",
     ylab="p(t_2 | y)", main="Bondad de ajuste Alambre tipo 2",
     ylim = c(0,0.000505))
lines(density(t_2), col = "goldenrod1", lwd = 2)

abline(v = tobs_2, col = "royalblue4", lwd = 2, lty = 1)
abline(v = quantile(x = t_2, probs = c(0.025, 0.975)), lty = 2, 
       lwd = 1, col = "gray")
legend("topright", legend = c("Posterior", "IC 95%", "t obs"), 
       col = c("goldenrod1", "gray", "royalblue4"),
       lty = c(1,2,1), lwd = 1, bty = "n")

#pr(t>tobs|y)
ppp_2<-mean(t_2>tobs_2)
ppp_2

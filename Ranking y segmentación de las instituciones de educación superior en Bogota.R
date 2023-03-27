library(readxl)


###
##MUESTREADORES


#Muestreador de Theta j
sample_thetaj <- function (nj, ybj, Muh, tau2, sig2j, thetaj, xi) 
{
  tau2n <- 1/(1/tau2 + nj/sig2j)
  mn  <- tau2n*(Muh[xi]/tau2 + nj*ybj/sig2j)
  thetaj <- rnorm(n = m, mean = mn, sd = sqrt(tau2n))
  return(thetaj)
}

#Muesteador de sigma2j
sample_sig2j <- function (nj, nu0, sig2, s2j, thetaj, ybj,sig2j) 
{
  alpha <- (nj + nu0)/2
  beta  <- (nu0*sig2 + (nj-1)*s2j + nj*(ybj - thetaj)^2)/2
  sig2j  <- 1/rgamma(n = m, shape = alpha, rate = beta)
  return(sig2j)
}

#Muestreador de Muh
sample_Muh <- function (nh, thetabh, mu0, gamma02, tau2, Muh) 
{
  for (h in 1:length(nh)) {
    if (nh[h] > 0) {
      gamman2 <- 1/(1/gamma02 + nh[h]/tau2)
      mu0n  <- gamman2*(mu0/gamma02 + nh[h]*thetabh[h]/tau2)
    } else {
      gamman2 <- gamma02
      mu0n  <- mu0
    }
    Muh[h] <- rnorm(n = 1, mean = mu0n, sd = sqrt(gamman2))
  }
  return(Muh)
}

#Muestreador de xi
sample_xi <- function (omega, tau2, xi, Muh, thetaj) 
{
  H <- length(omega)
  for (j in 1:length(thetaj)) {
    lp <- NULL
    for (h in 1:H) {
      lp[h] <- log(omega[h]) + dnorm(x = thetaj[j], mean = Muh[h], 
                                     sd = sqrt(tau2), log = T)
    }
    xi[j] <- sample(x = 1:H, size = 1, replace = F, prob = exp(lp - max(lp)))
  }
  return(xi)
}

#Muestreador de tau2
sample_tau2 <- function (eta0, tau02, tau2, xi, Muh, thetaj) 
{
  alpha <- (eta0 + length(thetaj))/2
  beta  <- (eta0*tau02 + sum((thetaj - Muh[xi])^2))/2
  tau2  <- 1/rgamma(n = 1, shape = alpha, rate = beta)
  return(tau2)
}

#Muestreador de sig2
sample_sig2<-function(alpha0,beta0,nu0,sig2j)
{
  alpha<-alpha0 +length(sig2j)*nu0/2
  beta<-beta0 + (nu0*sum(1/sig2j))/2
  sig2<-rgamma(n = 1, shape = alpha, rate = beta)
  return(sig2)  
}

#Muestreador de omega
sample_omega <- function (nh, alphasup0, omega)
{
  omega <- c(gtools::rdirichlet(n = 1, alpha = alphasup0 + nh))
  return(omega)
}



###
#Importar datos
data <- read_excel("C:/Users/Usuario/Downloads/data.xlsx")
y<-data$PUNT_GLOBAL
#estadisticos por institución
stats <- read_excel("C:/Users/Usuario/Downloads/data.xlsx", 
                    sheet = "Hoja2")

nj<-stats$`Cuenta de PUNT_GLOBAL2`
ybj<-stats$`Promedio de PUNT_GLOBAL`
s2j<-stats$`Var de PUNT_GLOBAL`
universidades<-stats$`Etiquetas de fila`
universidades<-gsub("-BOGOTÁ D.C.", "",universidades)


H<-3 #número de clusters
m<-length(nj) #número de instituciones




###
##MUESTREADOR DE GIBBS
MCMC <- function(n_sams,n_burn, nj, ybj, s2j, H) {

#hiperparametros
  
  nu0<-1
  mu0<-150
  gamma02<-30^2
  eta0<-1
  tau02<-30^2
  alpha0<-1
  beta0<-1/(30^2)
  alphasup0<-c(1/3,1/3,1/3)
  
  
# numero de iteraciones total
  B <- n_sams + n_burn

# valores iniciales
  set.seed(1)
  thetaj<-ybj
  cl<-kmeans(thetaj,3)
  sig2j<-s2j
  Muh<-cl$centers
  xi  <- cl$cluster
  omega <- table(xi)/m
  tau2  <- 1
  sig2  <- 100
  
# almacenamiento
  THETA <- matrix(data = NA, nrow = B, ncol = 3*m+2*H+2)
  LL    <- matrix(data = NA, nrow = B, ncol = 1)
  
# cadena
set.seed(1)
for (b in 1:B) {
  # actualizar estadisticos suficientes
  nh  <- as.numeric(table(factor(xi, levels = 1:H)))
  thetabh <- rep(NA, H)
  for (h in 1:H) if (nh[h] > 0) thetabh[h] <- mean(thetaj[xi == h])
  
  # actualizar parámetros
  Muh <- sample_Muh(nh, thetabh, mu0, gamma02, tau2, Muh)
  thetaj<-sample_thetaj(nj, ybj, Muh, tau2, sig2j, thetaj, xi) 
  sig2j<-sample_sig2j(nj, nu0, sig2, s2j, thetaj, ybj,sig2j) 
  tau2<-sample_tau2(eta0, tau02, tau2, xi, Muh, thetaj) 
  sig2<-sample_sig2(alpha0,beta0,nu0,sig2j)
  omega<-sample_omega (nh, alphasup0, omega)
  xi<-sample_xi(omega, tau2, xi, Muh, thetaj)

  # almacenar y log-verosimilitud
    THETA [b,]<- c(thetaj, sig2j,Muh, xi,tau2,sig2,omega)
    LL[b,] <- sum(dnorm(x = y, mean = rep(thetaj, nj), sd = sqrt(rep(sig2j, nj)), log = T))
}
  # fin de la cadena
  # salida
  colnames(THETA) <- c(paste0("theta"," ", universidades), 
                       paste0("sig2j"," ", universidades), paste0("mu", 1:H),
                       paste0("xi"," ", universidades), "tau2", "sig2", 
                       paste0("omega", 1:H))
  colnames(LL) <- c("ll")
  THETA <- as.data.frame(THETA)[(n_burn+1):B,]
  LL    <- as.data.frame(LL)[(n_burn+1):B,]
  return(list(THETA = THETA, LL = LL))
}



###
##Ajuste del modelo
set.seed(1)
Cadena<-MCMC(n_sams = 5000,n_burn = 0, nj, ybj, s2j, H)
plot(Cadena$LL,type = "p", pch = ".", cex.axis = 0.8, main = "", xlab = "Iteración", ylab = "Log-verosimilitud")
Cadena<-MCMC(n_sams = 5000,n_burn = 500, nj, ybj, s2j, H)
plot(Cadena$LL,type = "p", pch = ".", cex.axis = 0.8, main = "", xlab = "Iteración", ylab = "Log-verosimilitud")


THETA <- THETA[,]

###
##Rankimg
mm<-20
#Parámetros
THETA<-Cadena$THETA
ids2 <-universidades
that  <- colMeans(THETA[,1:m])
ic1   <- apply(X = THETA[,1:m], MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
ic2   <- apply(X = THETA[,1:m], MARGIN = 2, FUN = function(x) quantile(x, c(0.005,0.995)))

coef_var <- function(x) {
  sd(x) / mean(x)
}

CV   <- apply(X = THETA[,1:m], MARGIN = 2, FUN = function(x) coef_var(x))

indices <- tail(order(that),m) 
ranking<-tail(indices,mm) 
indices<-rev(indices)
#ranking <- order(that)
ids2 <- ids2[ ranking]
that <- that[ ranking]
CV<-round(CV[ ranking]*100,2)
CV<-paste("CV:",CV,"%")
dim(ic1)
ic1  <- ic1 [,ranking]
ic2  <- ic2 [,ranking]




par(mar=c(5,25,5,4))
plot(NA, NA, xlab = "Puntaje", ylab = "", main = "Ranking Bayesiano",
     xlim = c(100,200), ylim = c(1,mm), cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:mm, labels = ids2, las = 1,cex.axis = 0.8)
abline(v = 150,  col = "gray", lwd = 3)
abline(h = 1:mm, col = "lightgray", lwd = 1)
for (j in 1:mm) {
  segments(x0 = ic2[1,j], y0 = j, x1 = ic2[2,j], y1 = j, col = "orange", lwd = 3)
  segments(x0 = ic1[1,j], y0 = j, x1 = ic1[2,j], y1 = j, col = "royalblue", lwd = 3)
  lines(x = that[j], y = j, type = "p", pch = 16, cex = 0.8, col = "black")
}
#axis(side = 4, at = 1:mm, labels = CV, las = 1,cex.axis = 0.5)
text(that, (1:mm),
     labels = CV,
     cex = 0.9, pos = 3, col = "grey49")

legend(x = "topright",
       inset = c(0, -0.12), 
       legend = c("IC (95%)", "IC (99%)"), 
       lty = c(1, 1),
       col = c("royalblue", "orange"),
       lwd = 2,
       xpd = TRUE,
       seg.len = 3,
       cex=0.8)

 
###
#MATRIZ DE INCIDENCIA
A <- matrix(data = 0, nrow = m, ncol = m)
XI<-THETA[,(2*m+H+1):(3*m+H)]

B <- nrow(XI)
B_grid <- seq(from = 10, to =  B, by = 10)
B <- length(B_grid)
for (b in B_grid) {
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      if (XI[b,i] == XI[b,j]) {
        A[i,j] <- A[i,j] + 1/B
      } 
    }
  }
}
A <- A + t(A)
diag(A) <- 1


A  <- A[indices,indices]

# funcion para graficar las probabilidades por medio de un diagrama de calor
heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis,...)
{ 
  JJ <- dim(mat)[1]
  colorscale <- c("white", rev(heat.colors(100)))
  if(missing(labs))     labs <- 1:JJ
  if(missing(col.axis)) col.axis <- rep("black", JJ)
  if(missing(cex.axis)) cex.axis <- 0.5
  if(missing(tick))     tick <- TRUE
  ## adjacency matrix
  image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "Posición en el ranking", ylab = "Posición en el ranking",
        col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
  for(j in 1:JJ){
    axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick = tick, col.axis = col.axis[j])
    axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick = tick, col.axis = col.axis[j])
  }
  box()
  if(show.grid) grid(nx = JJ, ny = JJ)
}
par(mar=c(5,5,3,1))
heat.plot0(mat = A,main="Matriz de incidencia")


#######BONDAD DE AJUSTE
data<-as.data.frame(data)
B<-dim(THETA)[1]
MIN<-MAX<-IQR<-MEDIA<-MEDIANA<-SD<-matrix(NA,nrow=B,ncol = m)
pppMIN<-pppMAX<-pppIQR<-pppMEDIA<-pppMEDIANA<-pppSD<-NULL
unis<-stats$`Etiquetas de fila`
set.seed(1)
for (j in 1:m) {
for (i in 1:B) {
  ystar<-rnorm(1000,THETA[i,j],sqrt(THETA[i,m+j]))
  MIN[i,j]<-min(ystar)
  MAX[i,j]<-max(ystar)
  IQR[i,j]<-quantile(ystar)[4]-quantile(ystar)[2]
  MEDIA[i,j]<-mean(ystar)
  MEDIANA[i,j]<-median(ystar)
  SD[i,j]<-sd(ystar)
}
  pppMIN[j]<-mean(MIN[,j]<min(data[data[1] ==unis[j],2]))
  pppMAX[j]<-mean(MAX[,j]<max(data[data[1] ==unis[j],2]))
  pppIQR[j]<-mean(IQR[,j]<quantile(data[data[1] ==unis[j],2])[4]-
                    quantile(data[data[1] ==unis[j],2])[2])
  pppMEDIA[j]<-mean(MEDIA[,j]<mean(data[data[1] ==unis[j],2]))
  pppMEDIANA[j]<-mean(MEDIANA[,j]<median(data[data[1] ==unis[j],2]))
  pppSD[j]<-mean(SD[,j]<sd(data[data[1] ==unis[j],2]))
}
Estadistio<-rep(c("Min","Max","IQR","Media","Median","SD"),rep(72,6))
Estadistio<-as.factor(ppp_stat)
ppp<-c(pppMIN,pppMAX,pppIQR,pppMEDIA,pppMEDIANA,pppSD)
ppp<-data.frame(ppp,Estadistio)

gris<-rgb(0.6,0.6,0.7,0.3)
boxplot(ppp~Estadistio,data = ppp)
stripchart(ppp$ppp ~ ppp$Estadistio, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = gris)




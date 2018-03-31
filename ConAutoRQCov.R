
AS.CAViaR.Variance <- function(beta,y,theta,VaR){
  # y <- ins
  # VaR <- SAV.VaR01[1:5298]
  # beta <- SAV.beta01
  residuals <- y - VaR
  N <- length(residuals)
  sortedRes <- sort(abs(residuals))
  k <- ifelse(theta==0.01,40,60)
  bandwidth <- sortedRes[k]
  t=0
  derivative1 = rep(0,N)
  derivative2 = rep(0,N)
  derivative3 = rep(0,N)
  derivative4 = rep(0,N)
  gradient = matrix(0,N,4)
  D <- matrix(0,4,4)
  A = D
  for(i in 2:N){
    derivative1[i] <- 1+beta[2]*derivative1[i-1]
    derivative2[i] <- VaR[i-1]+beta[2]*derivative2[i-1]
    derivative3[i] <- beta[2]*derivative3[i-1]+max(y[i-1],0)
    derivative4[i] <- beta[2]*derivative4[i-1]-min(y[i-1],0)
    gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i],derivative4[i])
    A = A + (gradient[i,])%*%t(gradient[i,])
    if(abs(residuals[i])<=bandwidth){
      t=t+1
      D = D + (gradient[i,])%*%t(gradient[i,])
    }
  }
  tStrErr = t
  A = A/N
  D = D/(2*bandwidth*N)
  VCmatrix = theta*(1-theta)*solve(D)%*%A%*%solve(D)/N
  return(VCmatrix)
}
AS.Var <- AS.CAViaR.Variance(AS.beta01,ins,0.01,AS.VaR01[1:5298])
AS.std <- sqrt(diag(AS.Var))
AS.Pvalue <- pnorm(-abs(AS.beta01/AS.std))

setwd("C:\\Users\\samsung\\Desktop")
dat <- read.csv("SVAR.csv",sep='',header=T)
dat <- na.omit(dat)
st <- dat[,1]
dt <- dat[,2]
lr <- dat[,3]
md <- dat[,4]
iv <- dat[,5]

y <- iv
m <- st
start <- c(0.5,0.5,0.5,0.5)
theta <- 0.01 #0.05
res <- SAV_estim(y,m,theta,start)

SAV_estim <- function(y,m,theta,start){
  ini <- quantile(y,theta)
  n <- length(y)
  SAVobjective <- function(beta){
    q <- rep(NA,n)
    q[1] <- ini
    for(i in 2:n){
      q[i] <- beta[1]+beta[2]*q[i-1]+beta[3]*abs(y[i-1])+beta[4]*abs(m[i-1])
    }
    Hit <- ifelse(y<q,theta-1,theta)
    return(sum(Hit*(y-q))/n)
  }
  ui <- rbind(c(0,1,0,0),c(0,-1,0,0))
  ci <- c(0,-1)
  m <- constrOptim(start, SAVobjective, NULL,ui=ui,ci=ci,method = "Nelder-Mead")
  beta <- m$par
  q <- rep(NA,n)
  q[1] <- ini
  for(i in 2:n){
    q[i] <- beta[1]+beta[2]*q[i-1]+beta[3]*abs(y[i-1])
  }
  prop <- mean(y<q)
  return(list(prop=prop,q=q,beta=beta))
}

SAV_loop <- function(y,m,beta,theta){
  empirical <- quantile(y,theta)
  n <- length(y)
  q <- rep(NA,n)
  q[i] <- empirical
  for(i in 2:n){
    q[i] <- beta[1]+beta[2]*q[i-1]+beta[3]*abs(y[i-1])+beta[4]*abs(m[i-1])
  }
  return(q)
}

VaR <- SAV_loop(y,m,res$beta,theta)


SAV.CAViaR.Variance <- function(beta,y,theta,VaR,m){
  # y <- ins
  # VaR <- SAV.VaR01[1:5298]
  # beta <- SAV.beta01
  residuals <- y - VaR
  N <- length(residuals)
  sortedRes <- sort(abs(residuals))
  k <- ifelse(theta==0.01,40,60)
  bandwidth <- sortedRes[k]
  t=0
  derivative1 = rep(0,N)
  derivative2 = rep(0,N)
  derivative3 = rep(0,N)
  derivative4 = rep(0,N)
  gradient = matrix(0,N,4)
  D <- matrix(0,4,4)
  A = D
  for(i in 2:N){
    derivative1[i] <- 1+beta[2]*derivative1[i-1]
    derivative2[i] <- VaR[i-1]+beta[2]*derivative2[i-1]
    derivative3[i] <- beta[2]*derivative3[i-1]+abs(y[i-1])
    derivative4[i] <- beta[2]*derivative4[i-1]+abs(m[i-1])
    gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i],derivative4[i])
    A = A + (gradient[i,])%*%t(gradient[i,])
    if(abs(residuals[i])<=bandwidth){
      t=t+1
      D = D + (gradient[i,])%*%t(gradient[i,])
    }
  }
  tStrErr = t
  A = A/N
  D = D/(2*bandwidth*N)
  VCmatrix = theta*(1-theta)*solve(D)%*%A%*%solve(D)/N
  std <- diag(VCmatrix)
  Pvalue <- pnorm(-abs(beta/SAV.std))
  return(list(std=std,Pvalue=Pvalue))
}

y <- iv
M <- cbind(st,dt,lr,md)
theta <- c(0.01,0.05)
start <- c(0.5,0.5,0.5,0.5)
res <- simu.SAV(y,M,theta,start)

simu.SAV <- function(y,M,theta,start){
  kk <- ncol(M)
  tt <- length(theta)
  
  psi <- matrix(NA,24,4)
  
  for(i in 1:tt){
    for(j in 1:kk){
      m <- M[,j]
      est <- SAV_estim(y,m,theta[i],start)
      beta <- est$beta
      VaR <- est$q
      res <- SAV.CAViaR.Variance(beta,y,theta[i],VaR,m)
      
      std <- res$std
      Pvalue <- res$Pvalue
      
      
      esti1 <- rbind(beta,std,Pvalue)
      
      colres <- as.vector(esti1)
      
      if(i == 1){psi[1:12,j]<- colres}
      else{psi[13:24,j] <- colres}
      
      
    }
    
  }
  
  return(psi)
  
  
}
write.csv(round(res,6), file = "C:/Users/samsung/Desktop/huang1.csv", row.names = F, quote = F)













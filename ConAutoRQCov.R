SAV.CAViaR.Variance <- function(beta,y,theta,VaR){
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
  gradient = matrix(0,N,3)
  D <- matrix(0,3,3)
  A = D
  for(i in 2:N){
    derivative1[i] <- 1+beta[2]*derivative1[i-1]
    derivative2[i] <- VaR[i-1]+beta[2]*derivative2[i-1]
    derivative3[i] <- beta[2]*derivative3[i-1]+abs(y[i-1])
    gradient[i,] = c(derivative1[i],derivative2[i],derivative3[i])
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
SAV.Var <- SAV.CAViaR.Variance(SAV.beta01,ins,0.01,SAV.VaR01[1:5298])
SAV.std <- sqrt(diag(SAV.Var))
SAV.Pvalue <- pnorm(-abs(SAV.beta01/SAV.std))


write.table(VaR, file = "C:/Users/samsung/Desktop/paper2/VaR.txt", append =FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".",row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
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










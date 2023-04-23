library(resample); library(matlib)

dimred <- function(lambda, threshold=0.95){
  y <- cumsum(lambda)/sum(lambda)
  lambda <- lambda[y <= threshold]
  a <- length(lambda)
  a
}
dimred <- compiler:::cmpfun(dimred)

PCA_train <- function(X_train){
  m <- ncol(X_train); n_train <- dim(X_train)[1]
  b <- t(as.matrix(colMeans(X_train)))
  I <- matrix(1,nrow = n_train, ncol = 1)
  Sig <- matrix(0, m, m)
  vars <- colVars(X_train)
  for (i in 1:m){
    Sig[i,i] <- sqrt(vars[i])
  }
  S <- cor(X_train)
  S[is.na.data.frame(as.data.frame(S))] <- 0
  eig <- eigen(S)
  a <-dimred(eig$values, 0.95)
  D <- matrix(0, m, m)
  for (i in 1:m){
    D[i,i] <- eig$values[i]
  }
  V <- eig$vectors
  P <- V[,1:a]
  Da <- D[1:a,1:a]
  v <-0
  # Calculating T^2 alpha
  Fo <- qf(0.95,a,n_train-a)
  T2a <- ((a*(n_train-1)*(n_train+1))/(n_train*(n_train-a)))*Fo
  # Calculating Q alpha
  Theta1 <- Theta2 <- Theta3 <- 0
  alpha <- 0.05
  if(a == n_train){
    a <- n_train-1
  }
  for (i in (a+1):m){
    Theta1 <- Theta1+D[i,i]^1
    Theta2 <- Theta2+D[i,i]^2
    Theta3 <- Theta3+D[i,i]^3
  }
  h0 <- 1-(2*Theta1*Theta3)/(3*Theta2^2)
  Ca <- qnorm(1-alpha,0,1)
  Qa <- Theta1*(((h0*Ca*(2*Theta2)^0.5)/(Theta1))+1+((Theta2*h0*(h0-1))/(Theta1^2)))^(1/h0)
  Sig0 <- D^0.5;
  Sig_a <- Sig0[1:a ,1:a]
  out = list(T2a, Qa, a, P, Sig, b, Sig_a, S, Da, D)
  names(out) = c("T2a", "Qa", "a", "P", "Sig", "b", "Sig_a", "S", "Da", "D")
  out
}
PCA_train <- compiler:::cmpfun(PCA_train)

PCA_test <- function(X_train, X_test){
  out <- PCA_train(X_train)
  T2a <- out$T2a; Qa <- out$Qa; Sig <- out$Sig; b <- out$b; Sig_a <- out$Sig_a; 
  P <- out$P; a <- out$a
  n_test <- nrow(X_test); m = ncol(X_test)
  I <- matrix(1,nrow = n_test, ncol = 1)
  X <- (as.matrix(X_test - I%*%b))%*%inv(Sig)
  T2 <- rep(0, n_test); T2a <- T2a*rep(1, n_test); Q = T2; Qa <- Qa*rep(1, n_test)
  K=P%*%((inv(Sig_a))^2)%*%(t(P))
  K2=diag(1,m,m)-P%*%(t(P))
  for (i in 1:n_test) {
    T2[i] <- X[i,]%*%K%*%(matrix(data = X[i,], nrow = ncol(X), ncol = 1))
    r <- K2%*%(matrix(data = X[i,], nrow = ncol(X), ncol = 1))
    Q[i] <- (matrix(data = r, nrow = 1, ncol = ncol(X)))%*%r
  }
  output <- list(T2,T2a,Q,Qa,a)
  names(output) <- c("T2","T2a","Q","Qa","a")
  output
}
PCA_test <- compiler:::cmpfun(PCA_test)

closure <- function(x,k){
  out <- k*x/(sum(x))
}
closure <- compiler:::cmpfun(closure)

contribution <- function(XT,train_range,test_range){
  n <- nrow(XT);  m <- ncol(XT)
  X_train <- XT[train_range,]; n_train <- dim(X_train)[1]
  out <- PCA_train(X_train)
  T2a<-out$T2a; Qa<-out$Qa; a<-out$a; P<-out$P; Sig<-out$Sig; b<-out$b; Sig_a<-out$Sig_a
  S <- out$S; Da <- out$Da; D <- out$D; v <-0
  I <- matrix(1,nrow = n, ncol = 1)
  XT <- (as.matrix(XT - I%*%b))%*%inv(Sig) #(as.matrix(X_test - I%*%b))%*%inv(Sig)
  Contribution_T <- matrix(0, m,dim(as.matrix(test_range))[1])
  G <- P%*%((inv(Da))^2)%*%(t(P))
  G_Q <- diag(x = 1,m,m)-P%*%(t(P))
  for (k in test_range){
    r <- 0 
    tk <- (t(P)%*%(as.matrix(XT[k,])))
    v<-v+1
    size_tk_r <- sum((tk^2/as.matrix(diag(Da))) > (T2a[1])^(1/a))
    no <- FALSE
    if (size_tk_r == 0){
      size_tk_r <- r <- 1
      no <- TRUE
    }
    tk_r <- s <- index <- matrix(0, 1, size_tk_r)
    for (i in 1:a){
      if(tk[i]^2/D[i,i] > (T2a[1])^(1/a)){
        r<-r+1
        tk_r[r] <- tk[i]
        s[r] <- D[i,i]
        index[r] <- i
      }
    }
    for (j in 1:m){
      Cont_ij <- 0
      for(i in 1:r){
        if (no == TRUE){
          suma <- 0
        }else{
          suma <- (tk_r[i]/s[i])*P[j, index[i]]*XT[k,j]
        }
        if (suma >=0){
          Cont_ij <- Cont_ij + suma
        }
      }
      Contribution_T[j,v] <- Cont_ij
    }
  }
  rownames(Contribution_T) <- colnames(Data)
  output = Contribution_T
}
contribution <- compiler:::cmpfun(contribution)

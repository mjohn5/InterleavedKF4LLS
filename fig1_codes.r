



####################################################################################################################
## copy/pasting the entire set of codes below will generate a figure like Figure-1 in the paper (John and Wu, 2023)
####################################################################################################################

library(MASS)

M <- 5000 ## number of simulation-iterations

  n <- 100; p <- 6  ## number of rows and columns of 'covariate' matrices in each data-block


bias.1.mat <- bias.2.mat <- bias.3a.mat <- bias.3b.mat <- bias.4.mat <- matrix(, nrow = M, ncol = p)
mse.1.mat  <- mse.2.mat  <- mse.3a.mat  <- mse.3b.mat  <- mse.4.mat  <- matrix(, nrow = M, ncol = p)



#  Sigma <- matrix(rep(0.2, p*p), p, p); diag(Sigma) <- 1
  Sigma <- matrix(rep(0.0, p*p), p, p); diag(Sigma) <- 1

  alpha <- 0.25 ## try 0.5, 0.4, 0.25, 0.05

  mu.1 <- runif(1, -10, 10); mu.2 <- runif(1, -10, 10); mu.12 <- alpha*mu.1 + (1-alpha)*mu.2



  r.b <- as.matrix(runif(6, -5, 5), nrow = p, ncol = 1)
  r.f <- as.matrix(runif(6, -5, 5), nrow = p, ncol = 1)


for(m in 1:M){ if(m%%1000 == 0){print(m)} 



  U1  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))
  U2  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))
  U3  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))
  U4  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))
  U5  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))
  U6  <- as.matrix(mvrnorm(n, mu = rep(mu.1, p), Sigma = Sigma))

  V1  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))
  V2  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))
  V3  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))
  V4  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))
  V5  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))
  V6  <- as.matrix(mvrnorm(n, mu = rep(mu.2, p), Sigma = Sigma))


Z1 <- alpha*U1 + (1-alpha)*V1; Z2 <- alpha*U2 + (1-alpha)*V2; Z3 <- alpha*U3 + (1-alpha)*V3; Z4 <- alpha*U4 + (1-alpha)*V4; Z5 <- alpha*U5 + (1-alpha)*V5; Z6 <- alpha*U6 + (1-alpha)*V6;



  e.b <- as.matrix(rnorm(6*n, 0, 1), nrow = 6*n, ncol = 1);
  e.b1 <- e.b[1:n,]; e.b2 <- e.b[(n+1):(2*n),]; e.b3 <- e.b[((2*n)+1):(3*n),];  e.b4 <- e.b[((3*n)+1):(4*n),]; e.b5 <- e.b[((4*n)+1):(5*n),]; e.b6 <- e.b[((5*n)+1):(6*n),] 
  e.f <- as.matrix(rnorm(6*n, 0, 1), nrow = 6*n, ncol = 1)
  e.f1 <- e.f[1:n,]; e.f2 <- e.f[(n+1):(2*n),]; e.f3 <- e.f[((2*n)+1):(3*n),];  e.f4 <- e.f[((3*n)+1):(4*n),]; e.f5 <- e.f[((4*n)+1):(5*n),]; e.f6 <- e.f[((5*n)+1):(6*n),] 
  e.p <- as.matrix(rnorm(6*n, 0, 1), nrow = 6*n, ncol = 1)
  e.p1 <- e.p[1:n,]; e.p2 <- e.p[(n+1):(2*n),]; e.p3 <- e.p[((2*n)+1):(3*n),];  e.p4 <- e.p[((3*n)+1):(4*n),]; e.p5 <- e.p[((4*n)+1):(5*n),]; e.p6 <- e.p[((5*n)+1):(6*n),] 

  b1 <- U1%*%r.b + e.b1; b2 <- U2%*%r.b + e.b2; b3 <- U3%*%r.b + e.b3;  b4 <- U4%*%r.b + e.b4;  b5  <- U5%*%r.b + e.b5;  b6  <- U6%*%r.b + e.b6
  f1 <- V1%*%r.f + e.f1; f2 <- V2%*%r.f + e.f2; f3 <- V3%*%r.f + e.f3;  f4 <- V4%*%r.f + e.f4;  f5 <-  V5%*%r.f + e.f5;  f6 <-  V6%*%r.f + e.f6

  r.p <- alpha*r.b + (1-alpha)*r.f

  p1 <- Z1%*%r.p + e.p1; p2 <- Z2%*%r.p + e.p2; p3 <- Z3%*%r.p + e.p3;  p4 <- Z4%*%r.p + e.p4;  p5  <- Z5%*%r.p + e.p5;  p6  <- Z6%*%r.p + e.p6


##############
### KF4LLS
########

### (1) training with U's and b's only (i.e. only with 'bird' data blocks)

psi.0 <- as.matrix(rnorm(p), nrow = p, ncol = 1)

H1 <-      t(U1)%*%U1; H1.inv <- solve(H1); psi.1 <- psi.0 + H1.inv%*%t(U1)%*%(b1 - U1%*%psi.0)
H2 <- H1 + t(U2)%*%U2; H2.inv <- solve(H2); psi.2 <- psi.1 + H2.inv%*%t(U2)%*%(b2 - U2%*%psi.1)
H3 <- H2 + t(U3)%*%U3; H3.inv <- solve(H3); psi.3 <- psi.2 + H3.inv%*%t(U3)%*%(b3 - U3%*%psi.2)
H4 <- H3 + t(U4)%*%U4; H4.inv <- solve(H4); psi.4 <- psi.3 + H4.inv%*%t(U4)%*%(b4 - U4%*%psi.3)
H5 <- H4 + t(U5)%*%U5; H5.inv <- solve(H5); psi.5 <- psi.4 + H5.inv%*%t(U5)%*%(b5 - U5%*%psi.4)
H6 <- H5 + t(U6)%*%U6; H6.inv <- solve(H6); psi.6 <- psi.5 + H6.inv%*%t(U6)%*%(b6 - U6%*%psi.5)
  bias.1.mat[m,] <- c(mean(abs(psi.1 - r.p)), mean(abs(psi.2 - r.p)), mean(abs(psi.3 - r.p)), mean(abs(psi.4 - r.p)), mean(abs(psi.5 - r.p)), mean(abs(psi.6 - r.p)))
  res.sq.1 <-    (p1 - Z1%*%psi.1)^2
  res.sq.2 <- c( (p1 - Z1%*%psi.2)^2, (p2 - Z2%*%psi.2)^2 )
  res.sq.3 <- c( (p1 - Z1%*%psi.3)^2, (p2 - Z2%*%psi.3)^2,  (p3 - Z3%*%psi.3)^2);
  res.sq.4 <- c( (p1 - Z1%*%psi.4)^2, (p2 - Z2%*%psi.4)^2,  (p3 - Z3%*%psi.4)^2,  (p4 - Z4%*%psi.4)^2);
  res.sq.5 <- c( (p1 - Z1%*%psi.5)^2, (p2 - Z2%*%psi.5)^2,  (p3 - Z3%*%psi.5)^2,  (p4 - Z4%*%psi.5)^2, (p5 - Z5%*%psi.5)^2);
  res.sq.6 <- c( (p1 - Z1%*%psi.6)^2, (p2 - Z2%*%psi.6)^2,  (p3 - Z3%*%psi.6)^2,  (p4 - Z4%*%psi.6)^2, (p5 - Z5%*%psi.6)^2, (p6 - Z6%*%psi.6)^2);
  mse.1.mat[m,] <- c(mean(res.sq.1), mean(res.sq.2), mean(res.sq.3), mean(res.sq.4), mean(res.sq.5), mean(res.sq.6))
  if(m == 1){ beta.1.mat <- cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }
  if(m > 1){ beta.1.mat <- beta.1.mat + cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }


### (2) training with V's and f's only (i.e. only with 'fish' data blocks)

H1 <-      t(V1)%*%V1; H1.inv <- solve(H1); psi.1 <- psi.0 + H1.inv%*%t(V1)%*%(f1 - V1%*%psi.0)
H2 <- H1 + t(V2)%*%V2; H2.inv <- solve(H2); psi.2 <- psi.1 + H2.inv%*%t(V2)%*%(f2 - V2%*%psi.1)
H3 <- H2 + t(V3)%*%V3; H3.inv <- solve(H3); psi.3 <- psi.2 + H3.inv%*%t(V3)%*%(f3 - V3%*%psi.2)
H4 <- H3 + t(V4)%*%V4; H4.inv <- solve(H4); psi.4 <- psi.3 + H4.inv%*%t(V4)%*%(f4 - V4%*%psi.3)
H5 <- H4 + t(V5)%*%V5; H5.inv <- solve(H5); psi.5 <- psi.4 + H5.inv%*%t(V5)%*%(f5 - V5%*%psi.4)
H6 <- H5 + t(V6)%*%V6; H6.inv <- solve(H6); psi.6 <- psi.5 + H6.inv%*%t(V6)%*%(f6 - V6%*%psi.5)
  bias.2.mat[m,] <- c(mean(abs(psi.1 - r.p)), mean(abs(psi.2 - r.p)), mean(abs(psi.3 - r.p)), mean(abs(psi.4 - r.p)), mean(abs(psi.5 - r.p)), mean(abs(psi.6 - r.p)))
  res.sq.1 <-    (p1 - Z1%*%psi.1)^2
  res.sq.2 <- c( (p1 - Z1%*%psi.2)^2, (p2 - Z2%*%psi.2)^2 )
  res.sq.3 <- c( (p1 - Z1%*%psi.3)^2, (p2 - Z2%*%psi.3)^2,  (p3 - Z3%*%psi.3)^2);
  res.sq.4 <- c( (p1 - Z1%*%psi.4)^2, (p2 - Z2%*%psi.4)^2,  (p3 - Z3%*%psi.4)^2,  (p4 - Z4%*%psi.4)^2);
  res.sq.5 <- c( (p1 - Z1%*%psi.5)^2, (p2 - Z2%*%psi.5)^2,  (p3 - Z3%*%psi.5)^2,  (p4 - Z4%*%psi.5)^2, (p5 - Z5%*%psi.5)^2);
  res.sq.6 <- c( (p1 - Z1%*%psi.6)^2, (p2 - Z2%*%psi.6)^2,  (p3 - Z3%*%psi.6)^2,  (p4 - Z4%*%psi.6)^2, (p5 - Z5%*%psi.6)^2, (p6 - Z6%*%psi.6)^2);
  mse.2.mat[m,] <- c(mean(res.sq.1), mean(res.sq.2), mean(res.sq.3), mean(res.sq.4), mean(res.sq.5), mean(res.sq.6))
  if(m == 1){ beta.2.mat <- cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }
  if(m > 1){ beta.2.mat <- beta.2.mat + cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }


### (4) training with Z's and p's only (i.e. only with 'penguin' data blocks)

H1 <-      t(Z1)%*%Z1; H1.inv <- solve(H1); psi.1 <- psi.0 + H1.inv%*%t(Z1)%*%(p1 - Z1%*%psi.0)
H2 <- H1 + t(Z2)%*%Z2; H2.inv <- solve(H2); psi.2 <- psi.1 + H2.inv%*%t(Z2)%*%(p2 - Z2%*%psi.1)
H3 <- H2 + t(Z3)%*%Z3; H3.inv <- solve(H3); psi.3 <- psi.2 + H3.inv%*%t(Z3)%*%(p3 - Z3%*%psi.2)
H4 <- H3 + t(Z4)%*%Z4; H4.inv <- solve(H4); psi.4 <- psi.3 + H4.inv%*%t(Z4)%*%(p4 - Z4%*%psi.3)
H5 <- H4 + t(Z5)%*%Z5; H5.inv <- solve(H5); psi.5 <- psi.4 + H5.inv%*%t(Z5)%*%(p5 - Z5%*%psi.4)
H6 <- H5 + t(Z6)%*%Z6; H6.inv <- solve(H6); psi.6 <- psi.5 + H6.inv%*%t(Z6)%*%(p6 - Z6%*%psi.5)
  bias.4.mat[m,] <- c(mean(abs(psi.1 - r.p)), mean(abs(psi.2 - r.p)), mean(abs(psi.3 - r.p)), mean(abs(psi.4 - r.p)), mean(abs(psi.5 - r.p)), mean(abs(psi.6 - r.p)))
  res.sq.1 <-    (p1 - Z1%*%psi.1)^2
  res.sq.2 <- c( (p1 - Z1%*%psi.2)^2, (p2 - Z2%*%psi.2)^2 )
  res.sq.3 <- c( (p1 - Z1%*%psi.3)^2, (p2 - Z2%*%psi.3)^2,  (p3 - Z3%*%psi.3)^2);
  res.sq.4 <- c( (p1 - Z1%*%psi.4)^2, (p2 - Z2%*%psi.4)^2,  (p3 - Z3%*%psi.4)^2,  (p4 - Z4%*%psi.4)^2);
  res.sq.5 <- c( (p1 - Z1%*%psi.5)^2, (p2 - Z2%*%psi.5)^2,  (p3 - Z3%*%psi.5)^2,  (p4 - Z4%*%psi.5)^2, (p5 - Z5%*%psi.5)^2);
  res.sq.6 <- c( (p1 - Z1%*%psi.6)^2, (p2 - Z2%*%psi.6)^2,  (p3 - Z3%*%psi.6)^2,  (p4 - Z4%*%psi.6)^2, (p5 - Z5%*%psi.6)^2, (p6 - Z6%*%psi.6)^2);
  mse.4.mat[m,] <- c(mean(res.sq.1), mean(res.sq.2), mean(res.sq.3), mean(res.sq.4), mean(res.sq.5), mean(res.sq.6))
  if(m == 1){ beta.4.mat <- cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }
  if(m > 1){ beta.4.mat <- beta.4.mat + cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }


## centering all columns in all data blocks and target vectors, needed for the remaining set of codes 
U1 <- U1 - apply(U1,2,mean); U2 <- U2 - apply(U2,2,mean); U3 <- U3 - apply(U3,2,mean); U4 <- U4 - apply(U4,2,mean); U5 <- U5 - apply(U5,2,mean); U6 <- U6 - apply(U6,2,mean);  
V1 <- V1 - apply(V1,2,mean); V2 <- V2 - apply(V2,2,mean); V3 <- V3 - apply(V3,2,mean); V4 <- V4 - apply(V4,2,mean); V5 <- V5 - apply(V5,2,mean); V6 <- V6 - apply(V6,2,mean);  
b1 <- b1 - mean(b1); b2 <- b2 - mean(b2); b3 <- b3 - mean(b3);  b4 <- b4 - mean(b4);    b5 <- b5 - mean(b5);    b6 <- b6 - mean(b6)
f1 <- f1 - mean(f1); f2 <- f2 - mean(f2); f3 <- f3 - mean(f3);  f4 <- f4 - mean(f4);    f5 <- f5 - mean(f5);    f6 <- f6 - mean(f6)


## in the paper, after multiplying by sqrt(alpha) and sqrt(1-alpha), respectively, U's and V's were labelled with an alpha suffix attached.
## But, in the codes we don't have that suffix.
U1 <- sqrt(alpha)*U1;   U2 <- sqrt(alpha)*U2;   U3 <- sqrt(alpha)*U3;   U4 <- sqrt(alpha)*U4;   U5 <- sqrt(alpha)*U5;   U6 <- sqrt(alpha)*U6; 
V1 <- sqrt(1-alpha)*V1; V2 <- sqrt(1-alpha)*V2; V3 <- sqrt(1-alpha)*V3; V4 <- sqrt(1-alpha)*V4; V5 <- sqrt(1-alpha)*V5; V6 <- sqrt(1-alpha)*V6; 

b1 <- sqrt(alpha)*b1;   b2 <- sqrt(alpha)*b2;   b3 <- sqrt(alpha)*b3;   b4 <- sqrt(alpha)*b4;      b5 <- sqrt(alpha)*b5;     b6 <- sqrt(alpha)*b6; 
f1 <- sqrt(1-alpha)*f1; f2 <- sqrt(1-alpha)*f2; f3 <- sqrt(1-alpha)*f3; f4 <- sqrt(1-alpha)*f4;    f5 <- sqrt(1-alpha)*f5;   f6 <- sqrt(1-alpha)*f6; 

#
### Interleaved KF4LLS
#

### (3a) training with 'bird' blocks and 'fish' blocks interleaved, starting with a 'bird' block.


H1 <-      t(U1)%*%U1; H1.inv <- solve(H1); psi.1 <- psi.0 + H1.inv%*%t(U1)%*%(b1 - U1%*%psi.0)
H2 <- H1 + t(V1)%*%V1; H2.inv <- solve(H2); psi.2 <- psi.1 + H2.inv%*%t(V1)%*%(f1 - V1%*%psi.1)
H3 <- H2 + t(U2)%*%U2; H3.inv <- solve(H3); psi.3 <- psi.2 + H3.inv%*%t(U2)%*%(b2 - U2%*%psi.2)
H4 <- H3 + t(V2)%*%V2; H4.inv <- solve(H4); psi.4 <- psi.3 + H4.inv%*%t(V2)%*%(f2 - V2%*%psi.3)
H5 <- H4 + t(U3)%*%U3; H5.inv <- solve(H5); psi.5 <- psi.4 + H5.inv%*%t(U3)%*%(b3 - U3%*%psi.4)
H6 <- H5 + t(V3)%*%V3; H6.inv <- solve(H6); psi.6 <- psi.5 + H6.inv%*%t(V3)%*%(f3 - V3%*%psi.5)
  bias.3a.mat[m,] <- c(mean(abs(psi.1 - r.p)), mean(abs(psi.2 - r.p)), mean(abs(psi.3 - r.p)), mean(abs(psi.4 - r.p)), mean(abs(psi.5 - r.p)), mean(abs(psi.6 - r.p)))
  res.sq.1 <-    (p1 - Z1%*%psi.1)^2
  res.sq.2 <- c( (p1 - Z1%*%psi.2)^2, (p2 - Z2%*%psi.2)^2 )
  res.sq.3 <- c( (p1 - Z1%*%psi.3)^2, (p2 - Z2%*%psi.3)^2,  (p3 - Z3%*%psi.3)^2);
  res.sq.4 <- c( (p1 - Z1%*%psi.4)^2, (p2 - Z2%*%psi.4)^2,  (p3 - Z3%*%psi.4)^2,  (p4 - Z4%*%psi.4)^2);
  res.sq.5 <- c( (p1 - Z1%*%psi.5)^2, (p2 - Z2%*%psi.5)^2,  (p3 - Z3%*%psi.5)^2,  (p4 - Z4%*%psi.5)^2, (p5 - Z5%*%psi.5)^2);
  res.sq.6 <- c( (p1 - Z1%*%psi.6)^2, (p2 - Z2%*%psi.6)^2,  (p3 - Z3%*%psi.6)^2,  (p4 - Z4%*%psi.6)^2, (p5 - Z5%*%psi.6)^2, (p6 - Z6%*%psi.6)^2);
  mse.3a.mat[m,] <- c(mean(res.sq.1), mean(res.sq.2), mean(res.sq.3), mean(res.sq.4), mean(res.sq.5), mean(res.sq.6))
  if(m == 1){ beta.3a.mat <- cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }
  if(m > 1){ beta.3a.mat <- beta.3a.mat + cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }


### (3b) training with 'bird' blocks and 'fish' blocks interleaved, starting with a 'fish' block.


H1 <-      t(V1)%*%V1; H1.inv <- solve(H1); psi.1 <- psi.0 + H1.inv%*%t(V1)%*%(f1 - V1%*%psi.0)
H2 <- H1 + t(U1)%*%U1; H2.inv <- solve(H2); psi.2 <- psi.1 + H2.inv%*%t(U1)%*%(b1 - U1%*%psi.1)
H3 <- H2 + t(V2)%*%V2; H3.inv <- solve(H3); psi.3 <- psi.2 + H3.inv%*%t(V2)%*%(f2 - V2%*%psi.2)
H4 <- H3 + t(U2)%*%U2; H4.inv <- solve(H4); psi.4 <- psi.3 + H4.inv%*%t(U2)%*%(b2 - U2%*%psi.3)
H5 <- H4 + t(V3)%*%V3; H5.inv <- solve(H5); psi.5 <- psi.4 + H5.inv%*%t(V3)%*%(f3 - V3%*%psi.4)
H6 <- H5 + t(U3)%*%U3; H6.inv <- solve(H6); psi.6 <- psi.5 + H6.inv%*%t(U3)%*%(b3 - U3%*%psi.5)
  bias.3b.mat[m,] <- c(mean(abs(psi.1 - r.p)), mean(abs(psi.2 - r.p)), mean(abs(psi.3 - r.p)), mean(abs(psi.4 - r.p)), mean(abs(psi.5 - r.p)), mean(abs(psi.6 - r.p)))
  res.sq.1 <-    (p1 - Z1%*%psi.1)^2
  res.sq.2 <- c( (p1 - Z1%*%psi.2)^2, (p2 - Z2%*%psi.2)^2 )
  res.sq.3 <- c( (p1 - Z1%*%psi.3)^2, (p2 - Z2%*%psi.3)^2,  (p3 - Z3%*%psi.3)^2);
  res.sq.4 <- c( (p1 - Z1%*%psi.4)^2, (p2 - Z2%*%psi.4)^2,  (p3 - Z3%*%psi.4)^2,  (p4 - Z4%*%psi.4)^2);
  res.sq.5 <- c( (p1 - Z1%*%psi.5)^2, (p2 - Z2%*%psi.5)^2,  (p3 - Z3%*%psi.5)^2,  (p4 - Z4%*%psi.5)^2, (p5 - Z5%*%psi.5)^2);
  res.sq.6 <- c( (p1 - Z1%*%psi.6)^2, (p2 - Z2%*%psi.6)^2,  (p3 - Z3%*%psi.6)^2,  (p4 - Z4%*%psi.6)^2, (p5 - Z5%*%psi.6)^2, (p6 - Z6%*%psi.6)^2);
  mse.3b.mat[m,] <- c(mean(res.sq.1), mean(res.sq.2), mean(res.sq.3), mean(res.sq.4), mean(res.sq.5), mean(res.sq.6))
  if(m == 1){ beta.3b.mat <- cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }
  if(m > 1){ beta.3b.mat <- beta.3b.mat + cbind(r.p, r.b, r.f, psi.1, psi.2, psi.3, psi.4, psi.5, psi.6) }


   }



mean.bias.1.vec <- apply(apply(as.matrix(beta.1.mat[,4:9]/M) - matrix(rep(r.p, length(r.p)), nrow = length(r.p), ncol = length(r.p)), 2, abs), 2, mean)
bmin <- min(mean.bias.1.vec);          bmax <- max(mean.bias.1.vec)

mean.bias.2.vec <- apply(apply(as.matrix(beta.2.mat[,4:9]/M) - matrix(rep(r.p, length(r.p)), nrow = length(r.p), ncol = length(r.p)), 2, abs), 2, mean)
bmin <- min(c(bmin, mean.bias.2.vec)); bmax <- max(c(bmax, mean.bias.2.vec))

mean.bias.3a.vec <- apply(apply(as.matrix(beta.3a.mat[,4:9]/M) - matrix(rep(r.p, length(r.p)), nrow = length(r.p), ncol = length(r.p)), 2, abs), 2, mean)
bmin <- min(c(bmin, mean.bias.3a.vec)); bmax <- max(c(bmax, mean.bias.3a.vec))

mean.bias.3b.vec <- apply(apply(as.matrix(beta.3b.mat[,4:9]/M) - matrix(rep(r.p, length(r.p)), nrow = length(r.p), ncol = length(r.p)), 2, abs), 2, mean)
bmin <- min(c(bmin, mean.bias.3b.vec)); bmax <- max(c(bmax, mean.bias.3b.vec))

mean.bias.4.vec <- apply(apply(as.matrix(beta.4.mat[,4:9]/M) - matrix(rep(r.p, length(r.p)), nrow = length(r.p), ncol = length(r.p)), 2, abs), 2, mean)
bmin <- min(c(bmin, mean.bias.4.vec)); bmax <- max(c(bmax, mean.bias.4.vec))


   mean.mse.1.vec  <- apply(mse.1.mat,  2, mean); mmin <- min(mean.mse.1.vec);           mmax <- max(mean.mse.1.vec)
   mean.mse.2.vec  <- apply(mse.2.mat,  2, mean); mmin <- min(c(mmin, mean.mse.2.vec));  mmax <- max(c(mmax, mean.mse.2.vec))
   mean.mse.3a.vec <- apply(mse.3a.mat, 2, mean); mmin <- min(c(mmin, mean.mse.3a.vec)); mmax <- max(c(mmax, mean.mse.3a.vec))
   mean.mse.3b.vec <- apply(mse.3b.mat, 2, mean); mmin <- min(c(mmin, mean.mse.3b.vec)); mmax <- max(c(mmax, mean.mse.3b.vec))
   mean.mse.4.vec  <- apply(mse.4.mat,  2, mean); mmin <- min(c(mmin, mean.mse.4.vec));  mmax <- max(c(mmax, mean.mse.4.vec))



print(beta.1.mat/M)
print(beta.2.mat/M)
print(beta.3a.mat/M)
print(beta.3b.mat/M)
print(beta.4.mat/M)


#################################################################################################################
######## Figure - 1 ( we didn't set the seed, so the results will be (very) slightly different than in the paper 

par(mfrow = c(1,2))
#   plot(c(1:p), mean.bias.4.vec, type = "b", col = "green", xlab = "iterations", ylab = "bias over iterations", ylim =c(bmin, bmax), lwd = 2, pch = 16, cex.lab = 1.75, cex.axis = 1.75)
   plot(c(1:p), mean.bias.1.vec, type = "b", col = "red", xlab = "iterations", ylab = "", ylim =c(bmin, bmax), lwd = 2, pch = 16, cex.lab = 1.75, cex.axis = 1.75)
   mtext(expression(paste( plain("Bias"))),side=2,line=3.5, padj=1,at=1.5,cex=1.75)
   points(c(1:p), mean.bias.2.vec, col = "orange", pch = 16)
   lines(c(1:p),  mean.bias.2.vec, col = "orange", lwd = 2)
   points(c(1:p), mean.bias.3a.vec, col = "blue")
   lines(c(1:p),  mean.bias.3a.vec, col = "blue", lwd = 2, lty = 2)
   points(c(1:p), mean.bias.3b.vec, col = "blue")
   lines(c(1:p),  mean.bias.3b.vec, col = "blue", lwd = 2, lty = 2)
   points(c(1:p), mean.bias.4.vec, col = "green")
   lines(c(1:p),  mean.bias.4.vec, col = "green", lwd = 2)


   plot(c(1:p), mean.mse.1.vec, type = "b", col = "red", xlab = "iterations", ylab = "", ylim =c(0, mmax), lwd = 2, pch = 16, cex.lab = 1.75, cex.axis = 1.75)
   mtext(expression(paste( plain("MSE"))),side=2,line=3.5, padj=1, cex=1.75)
   points(c(1:p), mean.mse.2.vec, col = "orange", pch = 16)
   lines(c(1:p),  mean.mse.2.vec, col = "orange", lwd = 2)
   points(c(1:p), mean.mse.3a.vec, col = "blue")
   lines(c(1:p),  mean.mse.3a.vec, col = "blue", lwd = 2, lty = 2)
   points(c(1:p), mean.mse.3b.vec, col = "blue")
   lines(c(1:p),  mean.mse.3b.vec, col = "blue", lwd = 2, lty = 2)
   points(c(1:p), mean.mse.4.vec, col = "green")
   lines(c(1:p),  mean.mse.4.vec, col = "green", lwd = 2)
 


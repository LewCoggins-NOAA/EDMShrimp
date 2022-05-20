######################################################################################################
# Cheng-Han Tsai, Stephan Munch 10/Mar/2022 hGP-EDM  (refs Munch et al. 2017, Rogers and Munch 2020)
# 
# Warning : The example code and required packages are still improving and may be re-packaged in future
#           Contact me for any bugs !
#
# Issues to do: 1. Comparing the performance of R built-in optim or nlminb vs. rprop algorithm
#               2. Implementing Laplace Approximation method
#               3. Comparing the MAP estimates vs. MCMC and importance sampling estimates
#               4. Re-packaging and probably merging with the R package maintained by Tanya Rogers on github 
#                  or potentially by Steve's lab hubs
# 
#
#####################################################################################################
#####################################################################################################
# R-TMB hierarchical sepGP-EDM code example (version 1.2)
# Cheng-Han Tsai 20/Feb/2022
# 1. R-TMB package is required for automatic differentiation (AD)
# 2. The prior density uses TMB build-in R style functions of beta-distribution "dbeta()"
# 3. The hyper-parameter settings for priors are referenced to Munch et al. 2017, Rogers&Munch 2020
####################################################################################################
#setwd("~/Desktop/")
hgpEDMtmb_v1.cpp <-"
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // rr must match the dimension of X
  // upper bound of exp(logtau), exp(logg) < ysigma e.g. the var(y); exp(logrho) < 1.01
  // the mean of half-normal prior of W is set to 1 and thus the highest prob density is 0
  //
  DATA_MATRIX(X);
  DATA_VECTOR(y);
  DATA_INTEGER(nsite);
  DATA_INTEGER(nrblock);
  DATA_MATRIX(rr);
  DATA_SCALAR(ysigma2);
  PARAMETER_VECTOR(logW);
  PARAMETER(logtau);
  PARAMETER(logg);
  PARAMETER(logrho);
  //
  vector<Type> W = exp(logW);
  Type tau = exp(logtau);
  Type g = exp(logg);
  Type rho = exp(logrho);
  //
  vector<Type> v1(X.cols());
  vector<Type> v2(X.cols());
  vector<Type> r(X.cols());
  matrix<Type> K(X.rows(),X.rows());
  //
  int idx = 0;
  for(int i=0; i<K.rows(); i++) {
    for(int j=0; j<K.rows(); j++) {
      idx = (i-i%nrblock)/nrblock;
      r = rr.row(idx);
      v1 = X.row(i)-X.row(j);
      v1 = v1.abs();
      v2 = (W*v1/r)*(W*v1/r);
      K(i,j) = (1-rho)*tau*exp(-sum(v2));
    }
  }
  matrix<Type> G(K);
  vector<Type> gvec(K.rows());
  G.fill(0.0);
  gvec.fill(g);
  G.diagonal() = gvec;
  K = K+G;
  //
  matrix<Type> I1(nsite,nsite);
  matrix<Type> I2(nrblock,nrblock);
  vector<Type> ii(nsite);
  I1.fill(0.0);
  I2.fill(1.0);
  ii.fill(1.0);
  I1.diagonal() = ii;
  matrix<Type> B = kronecker(I1,I2);
  matrix<Type> Kf = K.array()*B.array();
  //
  vector<Type> unconstrained_params(nsite*(nsite-1)/2);
  matrix<Type> R(nsite,nsite);
  matrix<Type> S(nsite,nsite);
  matrix<Type> RS(nsite,nsite);
  unconstrained_params.fill(1.0);
  R = density::UNSTRUCTURED_CORR(unconstrained_params).cov();
  REPORT(R);
  S.fill(0.0);
  vector<Type> sdvec(nsite);
  sdvec.fill(sqrt(rho*tau));
  S.diagonal() = sdvec;
  RS = S*R*S;
  matrix<Type> Ks = kronecker(RS,I2);
  //
  matrix<Type> Kc = Kf + Ks;
  REPORT(Kc);
  Type nll = density::MVNORM_t<Type>(Kc)(y);
  vector<Type> v3 = W*W;
  Type v4 = (-0.5)*sum(v3)/(3.14159/2);
  Type v5 = dbeta(tau/ysigma2,Type(2.0),Type(2.0),true);
  Type v6 = dbeta(g/ysigma2,Type(2.0),Type(2.0),true);
  Type v7 = dbeta(rho/1.01,Type(2.0),Type(2.0),true);
  nll -= v4;
  nll -= v5;
  nll -= v6;
  nll -= v7;
  return nll;
}"
writeLines(hgpEDMtmb_v1.cpp, con="hgpEDMtmb_v1.cpp")
library(TMB)
compile("hgpEDMtmb_v1.cpp")
dyn.load(dynlib("hgpEDMtmb_v1"))
#########################################################################################
  
###############################################################################################
# R functions required for distance calculation and for hierarchical GP prediction
###############################################################################################
# Functions to compute distance between coordinate vectors with weighted length-scale parameter
dist.fn <- function(X,Y,rr,nrblock,w=1){
  D <- matrix(NA,nrow(X),nrow(Y))
  for(i in 0:(nrow(X)-1)){
    for(j in 0:(nrow(Y)-1)){
      idx <- (j-j%%nrblock)/nrblock+1
      r <- rr[idx,]
      v1 <- abs(X[i+1,]-Y[j+1,])
      v2 <- w*(v1/r)*w*(v1/r)
      D[i+1,j+1] <- sum(v2)
    }
  }
  return(D)
}
###########################################################################################################
# Function to compute MAP-based predictions
# Parameters: est (MAP estimates), 
#             X (time delay coordinates), y (vector of observations), 
#             XX (predicted time delay coordinates), KC (covariance matrix estimates)
#             R (dynamic correlation estimates)
#             rr (pre-processing scaling factors for time delay coordinates)
#             nrblock.pred (predicted data size for each spatial block)
#             nrblock (data size for each block)
#             
predstackGPsep <- function(est,X,y,XX,Kc=NULL,R,nsite,rr,nrblock.pred,nrblock) {
  phi.mle <- exp(est[1:ncol(X)])
  tau.mle <- exp(est[ncol(X)+1])
  g.mle <- exp(est[ncol(X)+2])
  rho.mle <- exp(est[ncol(X)+3]) # beware of parameter order 
  if(is.null(Kc)){ 
    D <- dist.fn(X=X,Y=X,rr=rr,nrblock=nrblock,w=phi.mle)
    K <- (1-rho.mle)*tau.mle*exp(-D) + diag(g.mle,nrow(D))
    B <-  kronecker(matrix(1,nsite,nsite),matrix(1,nrblock,nrblock))
    Kf <- K*B
    S <- diag(sqrt(rho.mle*tau.mle),nsite)
    RS <- S%*%R%*%S
    Ks <- kronecker(RS,matrix(1,nrblock,nrblock))
    Kc <- Kf+Ks
  }
  Ki <- solve(Kc)
  DXtheta <- dist.fn(X=XX,Y=X,rr=rr,nrblock=nrblock,w=phi.mle)
  KX <- (1-rho.mle)*tau.mle*exp(-DXtheta)
  KXf <- KX*kronecker(diag(1,nsite),matrix(1,nrblock.pred,nrblock))
  S <- diag(sqrt(rho.mle*tau.mle),nsite)
  RS <- S%*%R%*%S
  KXs <- kronecker(RS,matrix(1,nrblock.pred,nrblock))
  KXc <- KXf+KXs
  mup <- KXc%*%Ki%*%y
  #
  if(nrblock.pred==1){
    DXXtheta <- dist.fn(X=XX,Y=XX,rr=rr,nrblock=nrblock,w=phi.mle)
    KXX <- (1-rho.mle) * tau.mle * exp(-DXXtheta)
    Sigma <- (KXX+RS) - KXc%*%Ki%*%t(KXc)
    return(list(mup=mup, Sigma=Sigma))
  }
  else{
    DXXtheta <- dist.fn(XX,XX,rr,nrblock,w=phi.mle)
    KXX <- (1-rho.mle) * tau.mle * exp(-DXXtheta)
    KXXf <- KXX*kronecker(diag(1,nsite),matrix(1,nrblock.pred,nrblock))
    KXXc <- KXXf+KXs
    Sigma <- KXXc - KXc%*%Ki%*%t(KXc)
    return(list(mup=mup, Sigma=Sigma))
  }
}
############################################################################################




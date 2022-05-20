
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
}

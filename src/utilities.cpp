#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace std;
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; 
//           compile-command: "gcc -s -Wall -O3 -I/usr/share/R/include 
//                             -o rmath_qbeta rmath_qbeta.c -lRmath -lm" -*-


// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_grad_c(const arma::mat x,
                     const arma::vec y1,
                     const arma::vec y2){
  double a = mean(y1%y1);
  double b = mean(y2%y2);
  double d = mean(y1%y2);
  arma::mat intercept = arma::ones<arma::mat>(x.n_rows, 1);
  arma::mat x2 = join_rows(intercept, x);
  int P = x.n_cols;
  arma::vec first(P+1);
  arma::vec second(P+1);
  arma::vec third(P+1);
  arma::vec fourth(P+1);
  for (int p = 0; p < P+1; ++p){
    first(p) = sum(x2.col(p)) * (a*b*d - d*d*d);
    second(p) = sum(x2.col(p) % y1 % y1) * b * d;
    third(p) = sum(x2.col(p) % y2 % y2) * a * d;
    fourth(p) = sum(x2.col(p) % y1 % y2) * (a*b+d*d);
  }
  arma::vec out = first - second - third + fourth;
  return out;
}


// [[Rcpp::export]]
arma::mat get_fisher_c(const arma::mat x,
                       const arma::vec y1,
                       const arma::vec y2){
  double a = mean(y1%y1);
  double b = mean(y2%y2);
  double d = mean(y1%y2);
  int N = y1.size();
  arma::mat intercept = arma::ones<arma::mat>(x.n_rows, 1);
  arma::mat x2 = join_rows(intercept, x);
  arma::mat first = (a*b+d*d)*(a*b-d*d)*(x2.t() * x2);
  arma::rowvec colsumx = sum(x2, 0);
  arma::mat second = 4 * a * b * d * d * colsumx.t() * colsumx/N;
  return (first-second).i();
}

// [[Rcpp::export]]
double get_score_c(const arma::mat x,
                      const arma::vec y1,
                      const arma::vec y2){
  double a = mean(y1%y1);
  double b = mean(y2%y2);
  double d = mean(y1%y2);
  double mult = 1/(a*b-d*d);
  arma::vec grad = get_grad_c(x, y1, y2);
  arma::mat fisher = get_fisher_c(x, y1, y2);
  arma::mat val = mult * grad.t() * fisher * grad;
  return val(0,0);
}

// [[Rcpp::export]]
double get_degree_c(const arma::mat x,
                  const arma::vec y,
                  const arma::mat Y){
  int K = Y.n_cols;
  double d = 0;
  for (int k = 0; k < K; ++k){
    d = d + get_score_c(x, y, Y.col(k));
  }
  return d;
}

// [[Rcpp::export]]
double get_eta_c(const double rho12,
                    const double rho23,
                    const double rho13,
                    const double rho11,
                    const double rho22,
                    const double rho33){
  double a = rho11;
  double b = rho22;
  double c = rho33;
  double d = rho12;
  double e = rho13;
  double f = rho23;
  double num = 2*a*a*b*c*d*e + 3*a*c*d*d*d*e + 3*a*b*d*e*e*e;
  num = num - 2*a*a*b*b*c*e - a*a*b*e*e*f - a*b*b*e*e*e;
  num = num + a*a*a*b*c*f - 2*a*b*c*d*d*e - b*d*d*e*e*e;
  num = num - a*a*c*d*d*f + 2*d*d*d*e*e*e - a*a*b*c*d*d;
  num = num - a*b*d*d*e*e + 2*a*a*d*e*f*f;
  double denom = sqrt((a*b+d*d)*(a*b-d*d)*(a*b-d*d));
  denom = denom * sqrt((a*b+e*e)*(a*b-e*e)*(a*b-e*e));
  return num/denom;
}

// [[Rcpp::export]]
arma::mat get_H_c(const arma::mat Sigma){
  int K = Sigma.n_rows;
  arma::mat est_H = arma::zeros<arma::mat>(K-1, K-1);
  double eta;
  for (int i=1; i < K-1; ++i){
    for (int j=(i+1); j < K; ++j){
      eta = get_eta_c(Sigma(0,i), Sigma(i,j), Sigma(j,0),
                      Sigma(0,0), Sigma(i,i), Sigma(j,j));
      est_H(i-1, j-1) = eta;
      est_H(j-1, i-1) = eta;
    }
  }
  est_H.diag().ones();
  return est_H;
}

// [[Rcpp::export]]
arma::mat shuffle_row(arma::mat A) {
  for (int i=0; i < A.n_rows; ++i) {
    A.row(i) = shuffle( A.row(i), 1 );
  }
  return A;
}

// [[Rcpp::export]]
arma::cube shuffle_x_c(arma::mat x, 
                      int B){
  arma::cube X(x.n_rows, x.n_cols, B);
  for (int b=0; b < B; ++b){
    X.slice(b) = shuffle_row(x);
  }
  return X;
}

// [[Rcpp::export]]
arma::mat bootstrap_c(const arma::mat x,
                      const int B,
                      const arma::mat Y){
  arma::cube Xb = shuffle_x_c(x, B);
  const int K = Y.n_cols;
  arma::mat out(B, K);
  for (int b = 0; b < B; ++b){
    for (int k = 1; k < K; ++k){
      out(b, k) = get_score_c(Xb.slice(b), Y.col(0), Y.col(k));
    }
  }
  return out;
}

// [[Rcpp::export]]
double get_p_from_degree_c(const arma::vec y,
                           const arma::mat Y,
                           const double d,
                           const int numsim = 1000){
  arma::mat y2(y);
  y2.reshape(Y.n_rows, 1);
  arma::mat bigY = join_rows(y2, Y);
  int K = bigY.n_cols;
  arma::mat est_Sigma = cor(bigY);
  arma::mat est_H = get_H_c(est_Sigma);
  arma::vec lambda;
  arma::mat eigvec;
  eig_sym(lambda, eigvec, est_Sigma);
  arma::mat null_d(numsim, K-1);
  for (int k = 0; k < K-1; ++k){
    arma::vec tmp(numsim);
    for (int i = 0; i < numsim; ++i){
      double param1 = 0.5;
      double param2 = 2*lambda(k);
      tmp(i) = R::rgamma(param1, param2);
    }
    null_d.col(k) = tmp;
  }
  arma::vec rowsums = sum(null_d, 1);
  double p = 0;
  for (int i = 0; i < numsim; ++i){
    if (rowsums(i) > d){
      p += 1;
    }
  }
  p = p/numsim;
  return p;
}
  
  
// [[Rcpp::export]]
double generate_rgamma(int numsim){
  double out;
  for (int i = 0; i < numsim; ++i){
    out = R::rgamma(0.5, 3);
    cout << out << " ";
  }
  return out;
}


double get_A1_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = 24 * (k-1) * (p-1);
  double second = -24 * n * sum(TT.diag() % RR.diag());
  double third = 0;
  for (int i=0; i < n; ++i){
    for (int j=0; j < n; ++j){
      double tmp = RR(i,i) * RR(j,j);
      double tmp2 = RR(i,j) * RR(i,j) * 2;
      third += TT(i,j) * (tmp + tmp2);
    }
  }
  third *= 6 * n;
  return first + second + third;
}



double get_A2_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = -24 * (p-1) * (p-1);
  double second = 36 * n * accu(TT.diag() % TT.diag());
  double third = -48 * accu(TT%TT);
  double fourth = 0;
  for (int i=0; i < (n-1); ++i){
    for (int j=(i+1); j < n; ++j){
      fourth += TT(i,j) * (TT(j,j) * RR(i,i) + TT(i,i) * RR(j,j));
    }
  }
  for (int i=0; i < n; ++i){
    fourth += TT(i,i) * TT(i,i) * RR(i,i);
  }
  fourth *= (-24) * n;
  return (first + second + third + fourth);
}


double get_A3_c(int k,
                int p,
                arma::mat TT,
                arma::mat RR,
                int n){
  double first = 0;
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < n; ++j){
      first += TT(i,i) * TT(j,j) * TT(i,j);
    }
  }
  first *= 24 * n;
  double second = 16 * n * accu(TT%TT%TT);
  return first + second;
}


Rcpp::List get_TT_RR_c(arma::mat A){
  int n = A.n_rows;
  int P = A.n_cols;
  arma::mat X(n,P+1);
  arma::vec onevec = arma::ones<arma::vec>(n);
  X = join_rows(onevec, A);
  arma::rowvec barA = mean(A, 0);
  arma::mat H = arma::zeros<arma::mat>(P, P);
  for (int i=0; i < n; ++i){
    H += (A.row(i)-barA).t() * (A.row(i)-barA);
  }
  arma::mat Hinv = H.i();
  arma::mat TT(n,n);
  for (int i = 0; i < (n-1); ++i){
    for (int j = (i+1); j < n; ++j){
      arma::mat tmp = (A(i) - barA) * Hinv * (A(j) - barA).t();
      TT(i,j) =tmp(0,0);
      TT(j,i)= TT(i,j);
    }
  }
  for (int i=0; i < n; ++i){
    arma::mat tmp = (A.row(i)- barA) * Hinv * (A.row(i) - barA).t();
    TT(i,i) = tmp(0,0);
  }
  arma::mat RR = X * (X.t() * X).i() * X.t();
  return Rcpp::List::create(Rcpp::Named("TT") = TT,
                            Rcpp::Named("RR") = RR);
}



arma::vec cubic_coeff(arma::mat x){
  int p = x.n_cols + 1;
  int k = p;
  arma::mat A = orth(x) * sqrt(x.n_rows);
  Rcpp::List TTRR = get_TT_RR_c(A);
  arma::mat TT = as<arma::mat>(TTRR["TT"]);
  arma::mat RR = as<arma::mat>(TTRR["RR"]);
  int n = A.n_rows;
  double A1 = get_A1_c(k,p,TT,RR,n);
  double A2 = get_A2_c(k,p,TT,RR,n);
  double A3 = get_A3_c(k,p,TT,RR,n);
  arma::vec coef = arma::zeros<arma::vec>(3);
  coef(0) = (A3-A2+A1) / (12*n*(p-1))  + 1;
  coef(1) = (A2-2*A3)/(12*n*(p-1)*(p+1));
  coef(2) = A3/(12*n*(p-1)*(p+1)*(p+3));
  return coef;
}


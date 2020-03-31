#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

//const double log2pi = std::log(2.0 * M_PI);

arma::vec get_grad_c(const double sigw,
                      const double sigv,
                      const arma::vec uw,
                      const arma::vec uv,
                      const arma::mat x){

  arma::vec tmp1 = uw%uw / (sigw * sigw);
  arma::vec tmp2 = uv%uv / (sigv * sigv);
  arma::vec tmp3 = arma::zeros<arma::vec>(uw.size());
  tmp3.fill(1/sigw);
  arma::vec tmp4 = arma::zeros<arma::vec>(uw.size());
  tmp4.fill(1/sigv);
  arma::vec tmp = tmp3 - tmp1 - tmp4 + tmp2;
  int P = x.n_cols;
  arma::vec grad = arma::zeros<arma::vec>(P);
  for (int p = 0; p < P; ++p){
    grad(p) = sum(x.col(p) % tmp);
  }
  return grad;
}


// [[Rcpp::export]]
double const_c(const double sigw,
        const double sigv){
  double tmp = 2 * (1/(sigw*sigw) + 1/(sigv*sigv));
  return 1/tmp;
}

// //' get_score
// //' The order of y1 and y2 does not matter.
// //'
// //' @param y1 observation of the first variable
// //' @param y2 observation of the second variable
// //' @param x covariate vector
// //'
// //' @export
// //' @example
// //' x = rnorm(100)
// //' y1 = rnorm(100)
// //' y2 = rnorm(100)
// //' q = get_score(x, y1, y2)
// // [[Rcpp::export]]
// double get_score(const arma::mat x,
//             const arma::vec y1,
//             const arma::vec y2){
//   arma::mat B = orth(x) * sqrt(x.n_rows);
//   arma::vec w = y1 + y2;
//   arma::vec v = y1 - y2;
//   int n = w.size();
//   double sigw = sum(w%w)/n;
//   double sigv = sum(v%v)/n;
//   arma::vec grad = get_grad_c(sigw, sigv, w, v, B);
//   double fisherinf = const_c(sigw, sigv);
//   double q = sum(grad%grad) / n * fisherinf;
//   return q;
// }

// //' get_score_wv_c
// //'
// //' This function takes the sum of two variables, w, the difference of them, v, and one dimensional covariate x
// //' and returns the score test statistic
// //'
// //' @param x Covariate vector
// //' @param w sum of variable1 and variable2
// //' @param v difference of variable1 and variable2
// //' @export
// // [[Rcpp::export]]
// double get_score_wv_c(const arma::mat x,
//                      const arma::vec w,
//                      const arma::vec v) {
//   int n = w.size();
//   double sigw = sum(w%w)/n;
//   double sigv = sum(v%v)/n;
//   arma::vec grad = get_grad_c(sigw, sigv, w, v, x);
//   double fisherinf = const_c(sigw, sigv);
//   double q = sum(grad%grad) / n * fisherinf;
//   return q;
// }

//' mvrnormArma
//'
//' This function returns n samples of multivariate normal distribution with
//' mean mu and variance Sigma.
//'
//' @param n Sample size
//' @param mu Mean vector
//' @param Sigma Variance-covariance matrix
//' @export
// [[Rcpp::export]]
arma::mat mvrnormArma(const int n,
                      const arma::vec mu,
                      const arma::mat Sigma) {
  //returns random multivariate normal vectors with mean mu and covariance Sigma
  //input : integer n for the number of vectors you'd like to draw
  //      : vector mu for the mean
  //      : matrix Sigma for the covariance - needs to be psd
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}


// //' get_degree_c
// //'
// //' This function returns a sum statistic d for vector y tested with all other variables in matrix Y against covariate x
// //'
// //' @param x Covariate vector
// //' @param y Variable of interest
// //' @param Y All other variables to test with y
// //' @export
// // [[Rcpp::export]]
// double get_degree(const arma::mat x,
//                     const arma::vec y,
//                     const arma::mat Y){
//   //degree statistic for i'th gene based on matrix Y
//   arma::mat B = orth(x) * sqrt(x.n_rows);
//   int K = Y.n_cols;
//   double d = 0;
//   for (int i=0; i < K; ++i){
//     d = d + get_score(B, y, Y.col(i));
//   }
//   return d;
// }


//
// //' @export
// // [[Rcpp::export]]
// double get_eta(const double rho12,
//                  const double rho23,
//                  const double rho13){
//   double num = (rho23 + 2 * rho12 * rho23)* (rho12*rho12+1) * (rho13*rho13 + 1);
//   num = num + rho12 * rho13* (6 + 2 * rho12 + 2 * rho13 + 2 * rho23);
//   num = num - rho12 * (rho13*rho13 + 1) * (3*rho13 + rho13 + 2 * rho12*rho23);
//   num = num - rho13 * (rho12*rho12 + 1) * (3*rho12 + rho12 + 2 * rho13 * rho23);
//   double denom = (1-rho12*rho12)*(1-rho13*rho13) * sqrt(1+rho12*rho12) * sqrt(1+rho13*rho13);
//   return num/denom;
// }


// //' @export
// // [[Rcpp::export]]
// arma::mat get_H(const arma::mat Sigma){
//   int K = Sigma.n_rows;
//   arma::mat est_H = arma::zeros<arma::mat>(K-1, K-1);
//   double eta;
//   for (int i=1; i < K-1; ++i){
//     for (int j=(i+1); j < K; ++j){
//       eta = get_eta(Sigma(0,i), Sigma(i,j), Sigma(j,0));
//       est_H(i-1, j-1) = eta;
//       est_H(j-1, i-1) = eta;
//     }
//   }
//   est_H.diag().ones();
//   return est_H;
// }
//


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




arma::mat shuffle_x_c(arma::vec x, int B){
  arma::mat X(x.size(), B);
  for (int b=0; b < B; ++b){
    X.col(b) = shuffle(x);
  }
  return X;
}


arma::mat store_W_c(const arma::vec y,
                    const arma::mat smallY){
  int n = smallY.n_rows;
  int K = smallY.n_cols;
  arma::mat W(n, K);
  for (int k = 0; k < K; ++k){
    W.col(k) = y + smallY.col(k);
  }
  return W;
}

arma::mat store_V_c(const arma::vec y,
                    const arma::mat smallY){
  int n = smallY.n_rows;
  int K = smallY.n_cols;
  arma::mat V(n, K);
  for (int k = 0; k < K; ++k){
    V.col(k) = y - smallY.col(k);
  }
  return V;
}


// //' bootstrap_c
// //'
// //' This function performs permutation test for given covariate x B times
// //' Also takes as input W instead of Y matrix
// //'
// //' @param x Covariate vector
// //' @param B number of permutations
// //' @param W sum of the main variable and each target
// //' @param V difference of the main variable and each target
// //' @export
// // [[Rcpp::export]]
// arma::mat bootstrap_c(const arma::vec x,
//                       const int B,
//                       const arma::mat W,
//                       const arma::mat V){
//   arma::mat Xb = shuffle_x_c(x, B);
//   const int K = W.n_cols+1;
//   arma::mat out(B, K-1);
//   for (int b = 0; b < B; ++b){
//     for (int k = 0; k < (K-1); ++k){
//       out(b, k) = get_score_wv_c(Xb.col(b), W.col(k), V.col(k));
//     }
//   }
//   return out;
// }
//
//
//
//







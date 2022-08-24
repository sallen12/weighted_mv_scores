// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;     

////////////////////////////////////////////////////////////////////////////////
// kernel functions

// Euclidean distance
// [[Rcpp::export]]
double euclnormC(arma::colvec x){
  double out = sqrt(sum(square(x)));
  return(out);
}

// ims kernel
// [[Rcpp::export]]
double imskernelC(arma::colvec x){
  double out = -pow(1 + pow(euclnormC(x), 2.0), -0.50);
  return(out);
}

// variogram score kernel
// [[Rcpp::export]]
double vskernelC(arma::colvec x1, arma::colvec x2, arma::mat w_vs, double p){
  
  double out = 0;
  double d = x1.size();
  for (int i = 1; i < (d+1); i++) {
    for (int j = i; j < (d+1); j++) {
      double vx1 = pow(abs(x1[i-1] - x1[j-1]), p);
      double vx2 = pow(abs(x2[i-1] - x2[j-1]), p);
      out += 2*w_vs(i-1,j-1)*pow(vx1 - vx2, 2.0);
    }
  }
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// multivariate scores

// inverse multiquadric score
// [[Rcpp::export]]
double imsC(arma::colvec y, arma::mat dat, NumericVector w){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*imskernelC(dat.col(i-1) - y);
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = 1; j < (m+1); j++) {
      s2 += w[i-1]*w[j-1]*imskernelC(dat.col(i-1) - dat.col(j-1));
    }
  }
  
  double s3 = imskernelC(y - y);
  
  double out = s1 - s2 / 2 - s3 / 2;
  return (out);
  
}

// variogram score (with weights)
// [[Rcpp::export]]
double vsC_w(arma::colvec y, arma::mat dat, arma::mat w_vs, arma::colvec w, double p){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*vskernelC(dat.col(i-1), y, w_vs, p);
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = 1; j < (m+1); j++) {
      s2 += w[i-1]*w[j-1]*vskernelC(dat.col(i-1), dat.col(j-1), w_vs, p);
    }
  }
  double out = s1 - s2 / 2;
  return (out);
}

////////////////////////////////////////////////////////////////////////////////
// vertically re-scaled multivariate scores

// vertically re-scaled energy score
// [[Rcpp::export]]
double vresC(arma::colvec y, arma::mat dat, double w_y, NumericVector w_dat, arma::colvec x0, NumericVector w){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*euclnormC(dat.col(i-1) - y)*w_y*w_dat[i-1];
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = i; j < (m+1); j++) {
      s2 += 2*w[i-1]*w[j-1]*euclnormC(dat.col(i-1) - dat.col(j-1))*w_dat[i-1]*w_dat[j-1];
    }
  }
  
  double s3 = 0;
  for (int i = 1; i < (m+1); i++) {
    s3 += w[i-1]*euclnormC(dat.col(i-1) - x0)*w_dat[i-1];
  }
  s3 += -euclnormC(y - x0)*w_y;
  
  double s4 = 0;
  for (int i = 1; i < (m+1); i++) {
    s4 += w[i-1]*w_dat[i-1];
  }
  s4 += -w_y;
  
  double out = s1 - s2 / 2 + (s3 * s4);
  return (out);
  
}

// vertically re-scaled variogram score (with weights)
// [[Rcpp::export]]
double vrvsC(arma::colvec y, arma::mat dat, double w_y, NumericVector w_dat, arma::colvec x0, arma::mat w_vs, arma::colvec w, double p){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*vskernelC(dat.col(i-1), y, w_vs, p)*w_y*w_dat[i-1];
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = 1; j < (m+1); j++) {
      s2 += w[i-1]*w[j-1]*vskernelC(dat.col(i-1), dat.col(j-1), w_vs, p)*w_dat[i-1]*w_dat[j-1];
    }
  }
  
  double s3 = 0;
  for (int i = 1; i < (m+1); i++) {
    s3 += w[i-1]*vskernelC(dat.col(i-1), x0, w_vs, p)*w_dat[i-1];
  }
  s3 += -vskernelC(y, x0, w_vs, p)*w_y;
  
  double s4 = 0;
  for (int i = 1; i < (m+1); i++) {
    s4 += w[i-1]*w_dat[i-1];
  }
  s4 += -w_y;
  
  double out = s1 - s2 / 2 + (s3 * s4);
  return (out);
}

// vertically re-scaled inverse multiquadric score
// [[Rcpp::export]]
double vrimsC(arma::colvec y, arma::mat dat, double w_y, NumericVector w_dat, NumericVector w){
  
  double s1 = 0;
  double m = dat.n_cols;
  for (int i = 1; i < (m+1); i++) {
    s1 += w[i-1]*imskernelC(dat.col(i-1) - y)*w_y*w_dat[i-1];
  }
  
  double s2 = 0;
  for (int i = 1; i < (m+1); i++) {
    for (int j = 1; j < (m+1); j++) {
      s2 += w[i-1]*w[j-1]*imskernelC(dat.col(i-1) - dat.col(j-1))*w_dat[i-1]*w_dat[j-1];
    }
  }
  
  double s3 = imskernelC(y - y)*w_y*w_y;
  
  double out = s1 - s2 / 2 - s3 / 2;
  return (out);
  
}


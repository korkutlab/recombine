#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;


// products of sample pairs for each feature
// [[Rcpp::export]]
IntegerMatrix multfun(IntegerMatrix x_) {
  // Declarations
  int n = x_.nrow();
  int p = x_.ncol();
  int n2 = n*(n - 1)/2;

  // fill na with 0
  IntegerMatrix x(n, p);
  for (int i = 0; i < n*p; i++) {
    if (IntegerVector::is_na(x_[i])) x[i] = 0;
    else x[i] = x_[i];
  }

  // create containers for R output
  IntegerMatrix d(n2, p);

  // calculation
  int ii = 0;
  for (int i = 0; i < n-1; i++) {
    for (int ip = i+1; ip < n; ip++) {
      for (int j = 0; j < p; j++)
        d[ii + j*n2] = x[i + j*n]*x[ip + j*n]; // R matrix is column-major internally
      ii += 1;
    }
  }

  return d;
}


// absolute distance of sample pairs for each feature
// [[Rcpp::export]]
NumericMatrix distfun(NumericMatrix x_) {
  // Declarations
  int n = x_.nrow();
  int p = x_.ncol();
  int n2 = n*(n - 1)/2;

  // fill na with 0
  NumericMatrix x(n, p);
  for (int i = 0; i < n*p; i++) {
    if (NumericVector::is_na(x_[i])) x[i] = 0;
    else x[i] = x_[i];
  }

  // create containers for R output
  NumericMatrix d(n2, p);

  // calculation
  int ii = 0;
  for (int i = 0; i < n-1; i++) {
    for (int ip = i+1; ip < n; ip++) {
      for (int j = 0; j < p; j++)
        d[ii + j*n2] = fabs(x[i + j*n] - x[ip + j*n]); // R matrix is column-major internally
      ii += 1;
    }
  }

  return d;
}


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "lasso.h"


// Lasso -------

double soft(double x, double d) {
  double value = fabs(x) - d > 0.0 ? fabs(x) - d : 0.0; // max
  return copysign(1.0, x)*value;
}

double sum(double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += x[i];
  }
  return s;
}

double max(double *x, int n) {
  double m = -DBL_MAX;
  for (int i = 0; i < n; i++) {
    m = m > x[i] ? m : x[i];
  }
  return m;
}

void pmax(double *x, int n, double cutoff) {
  for (int i = 0; i < n; i++) {
    x[i] = x[i] > cutoff ? x[i] : cutoff;
  }
}

double max_fabs(double *x, int n) {
  double m = 0;
  for (int i = 0; i < n; i++) {
    m = m > fabs(x[i]) ? m : fabs(x[i]);
  }
  return m;
}

double l2n(double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += pow(x[i], 2);
  }
  return sqrt(s);
}

double l1n(double *x, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += fabs(x[i]);
  }
  return s;
}

double diff_l1n(double *x1, double *x2, int n) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += fabs(x2[i] - x1[i]);
  }
  return s;
}

void normalize_by_l2n(double *x, int n) {
  double x_l2n = l2n(x, n);
  if (x_l2n > 1e-10) {
    for (int i = 0; i < n; i++) {
      x[i] /= x_l2n;
    }
  } else {
    for (int i = 0; i < n; i++) {
      x[i] = 0.0;
    }
  }
}

double l1n_after_normalized_by_l2n(double *x, int n) {
  double x_l2n = l2n(x, n);
  double s = 0.0;
  if (x_l2n > 1e-10) {
    for (int i = 0; i < n; i++) {
      s += fabs(x[i]/x_l2n);
    }
  }
  return s;
}

double BinarySearch(double *argu, int n, double sumabs) {
  if(l2n(argu, n) == 0 || l1n_after_normalized_by_l2n(argu, n) <= sumabs) return 0;
  double lam1 = 0;
  double lam2 = max_fabs(argu, n) - 1e-5;
  int iter = 1;
  double *su = new double[n];
  while(iter <= 15 && (lam2 - lam1) > 1e-4) {
    for (int i = 0; i < n; i++) {
      su[i] = soft(argu[i], (lam1 + lam2)/2);
    }
    if(l1n_after_normalized_by_l2n(su, n) < sumabs) {
      lam2 = (lam1 + lam2)/2;
    } else {
      lam1 = (lam1 + lam2)/2;
    }
    iter += 1;
  }
  delete [] su;
  return (lam1 + lam2)/2;
}

// double rand_gen() {
//   // return a uniformly distributed random value
//   return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
// }
// 
// double normalRandom() {
//   // return a normally distributed random value
//   double v1=rand_gen();
//   double v2=rand_gen();
//   return cos(2*3.14*v2)*sqrt(-2.*log(v1));
// }

void matrix_multiply(double *A, double *B, int m, int n, int p, double *C) {
  int i, j, k;
  for (i = 0; i < m; i++) {
    for (j = 0; j < p; j++) {
      C[i + j*m] = 0; // R matrix is column-major internally
      for (k = 0; k < n; k++)
        C[i + j*m] += A[i + k*m] * B[k + j*n];
    }
  }
}

double cross_prod(double *x1, double *x2, int n) {
  double val = 0;
  for (int i = 0; i < n; i++) val += x1[i]*x2[i];
  return(val);
}


// Spike-and-slab Lasso ------
// This code has been adapted from the SSLASSO package (Rockova, et al, 2018)

double expectation_approx(double *beta, double a, double b, int p, int l) {
  int sum = 0;
  for (int i = 0; i < p; i++) {
    if(beta[l*p + i] != 0) sum++;
  }

  return (sum+a)/(a+b+p);
}

double pstar(double x, double theta, double lambda1, double lambda0) {
  if (lambda1 == lambda0) return 1;

  double value;
  value = (1-theta)/theta*lambda0/lambda1*exp(-fabs(x)*(lambda0-lambda1));
  value += 1;
  value = 1/value;

  return value;
}

double lambdastar(double x, double theta, double lambda1, double lambda0){
  if (lambda1 == lambda0) return lambda1;

  double aux;
  aux = pstar(x,theta,lambda1,lambda0);

  return aux*lambda1 + (1-aux)*lambda0;
}

double g(double x, double theta, double lambda1, double lambda0){
  double value = lambdastar(x,theta,lambda1,lambda0);

  return pow((value-lambda1),2) + 2*log(pstar(x,theta,lambda1,lambda0));
}

double threshold(double theta, double lambda1, double lambda0){
  if (lambda0 == lambda1) return lambda1;

  if ( g(0,theta,lambda1,lambda0) > 0) {
    return sqrt(2*log(1/pstar(0,theta,lambda1,lambda0))) + lambda1;
  } else {
    return lambdastar(0,theta,lambda1,lambda0);
  }
}

// Cross product of u with jth column of D
double crossprod_Dj_u(double *D, double *u, int n2, int j) {
  int nn = n2*j;

  double val = 0;
  for (int i = 0; i < n2; i++) val += D[nn+i]*u[i]; // R matrix is column-major internally

  return val;
}

double SSL(double dtu, double w, double lambda0, double lambda1, double theta, double delta) {
  double s = 0;
  double lambda;

  if (dtu > 0) s = 1;
  else if (dtu < 0) s = -1;

  if (fabs(dtu) <= delta) return 0;
  else {
    lambda = lambdastar(w, theta, lambda1, lambda0);

    double temp;
    temp = fabs(dtu) - lambda;

    if (temp > 0) return temp*s;
    else return 0;
  }
}

int checkConvergence(double *beta, double *beta_old, double eps, int p) {
  int converged = 1;

  for (int j=0; j<p; j++) {
    if (fabs((beta[j]-beta_old[j])/beta_old[j]) > eps) {
      converged = 0;
      break;
    }
  }

  return converged;
}


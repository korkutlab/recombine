#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <Rcpp.h>
#include "lasso.h"
#include "flsa.h"

using namespace Rcpp;


// Lasso as a constrained optimization problem
// [[Rcpp::export]]
List SHC_lasso_getuw(NumericMatrix ds_,
                     NumericVector wbounds_,
                     int max_iter,
                     bool init_random,
                     bool silent,
                     bool warm_start) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int L = wbounds_.length();
  //int n = int(sqrt(2 * n2)) + 1;
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *wbounds = new double[L];
  for (int i = 0; i < L; i++)
    wbounds[i] = wbounds_[i];
  
  double wbound;
  
  // create containers for R output
  double *u = new double[L*n2];
  double *w = new double[L*p];
  double *crit = new double[L];
  int *iter = new int[L];
  
  // init
  double *pu = NULL;
  
  if (init_random) {
    NumericVector rand_numbers = Rcpp::runif(L*p);
    for (int i = 0; i < L*p; i++) {
      w[i] = rand_numbers[i]*2.0/sqrt(p);
    }
  } else {
    for (int i = 0; i < L*p; i++) {
      w[i] = 1.0/sqrt(p);
    }
  }
  
  double *pw = NULL;
  
  double *w_old = new double[p];
  
  double *dw = new double[n2];
  
  for (int i = 0; i < L; i++) {
    iter[i] = 0;
  }
  
  double *dtu = new double[p];
  
  double lam;
  
  // Regularization Path in a Reverse Way (starting from less penalty)
  for (int l = L-1; l > -1; l--) {
    Rcpp::checkUserInterrupt();
    
    if(!silent) Rcpp::Rcout << L-l;
    R_FlushConsole();
    
    wbound = wbounds[l];
    // printf("wbound = %f\n", wbound);
    
    // set pointers
    pu = u + l*n2;
    pw = w + l*p;
    
    // warm start
    if (warm_start) {
      if (l != L-1) {
        for (int i = 0; i < p; i++) {
          w[l*p + i] = w[(l+1)*p + i];
        }
      }
    }
    
    // if w are all zeros due to warm start, no iteration wil be performed
    
    // main loop
    while(l1n(pw, p) > 1e-10 && iter[l] <= max_iter) {
      iter[l]++;
      
      for (int i = 0; i < p; i++) {
        w_old[i] = pw[i];
      }
      
      // update u
      matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
      normalize_by_l2n(pu, n2);
      
      // update w
      matrix_multiply(pu, ds, 1, n2, p, dtu); // D'u = (u'D)'
      pmax(dtu, p, 0.0);
      lam = BinarySearch(dtu, p, wbound);
      for (int i = 0; i < p; i++) {
        pw[i] = soft(dtu[i], lam);
      }
      normalize_by_l2n(pw, p);
      
      // check convergence
      if (diff_l1n(w_old, pw, p)/l1n(w_old, p) < 1e-4) {
        break;
      }
    }
    
    // update u
    matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
    normalize_by_l2n(pu, n2);
    
    // crit = u'Dw
    matrix_multiply(ds, pw, n2, p, 1, dw);
    crit[l] = cross_prod(pu, dw, n2);
  }
  if(!silent) Rcpp::Rcout << std::endl;
  
  // prepare results
  NumericVector u_(L*n2), w_(L*p), crit_(L);
  IntegerVector iter_(L);
  
  for (int i = 0; i < L*n2; i++)
    u_[i] = u[i];
  for (int i = 0; i < L*p; i++)
    w_[i] = w[i];
  for (int i = 0; i < L; i++)
    crit_[i] = crit[i];
  for (int i = 0; i < L; i++)
    iter_[i] = iter[i];
  
  List res;
  res["u"] = u_;
  res["w"] = w_;
  res["crit"] = crit_;
  res["iter"] = iter_;
  
  // cleanup
  delete [] ds;
  delete [] wbounds;
  
  delete [] u;
  delete [] w;
  delete [] crit;
  delete [] iter;
  
  delete [] w_old;
  delete [] dtu;
  delete [] dw;
  
  return res;
}


// Lasso as a penalty in the form of Lagrangian
// [[Rcpp::export]]
List SHC_lasso_lagrange_getuw(NumericMatrix ds_,
                              NumericVector lambdas_,
                              int max_iter,
                              bool init_random,
                              bool silent,
                              bool warm_start) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int L = lambdas_.length();
  //int n = int(sqrt(2 * n2)) + 1;
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *lambdas = new double[L];
  for (int i = 0; i < L; i++)
    lambdas[i] = lambdas_[i];
  
  double lambda;
  
  // create containers for R output
  double *u = new double[L*n2];
  double *w = new double[L*p];
  double *crit = new double[L];
  int *iter = new int[L];
  
  // init
  double *pu = NULL;
  
  if (init_random) {
    NumericVector rand_numbers = Rcpp::runif(L*p);
    for (int i = 0; i < L*p; i++) {
      w[i] = rand_numbers[i]*2.0/sqrt(p);
    }
  } else {
    for (int i = 0; i < L*p; i++) {
      w[i] = 1.0/sqrt(p);
    }
  }
  
  double *pw = NULL;
  
  double *w_old = new double[p];
  
  double *dw = new double[n2];
  
  for (int i = 0; i < L; i++) {
    iter[i] = 0;
  }
  
  double *dtu = new double[p];
  
  // Regularization Path
  for (int l = 0; l < L; l++) {
    Rcpp::checkUserInterrupt();
    
    if(!silent) Rcpp::Rcout << l+1;
    R_FlushConsole();
    
    lambda = lambdas[l];
    // printf("lambda = %f\n", lambda);
    
    // set pointers
    pu = u + l*n2;
    pw = w + l*p;
    
    // warm start
    if (warm_start) {
      if (l != 0) {
        for (int i = 0; i < p; i++) {
          w[l*p + i] = w[(l-1)*p + i];
        }
      }
    }
    
    // if w are all zeros due to warm start, no iteration wil be performed
    
    // main loop
    while (l1n(pw, p) > 1e-10 && iter[l] < max_iter) {
      iter[l]++;
      
      for (int i = 0; i < p; i++) {
        w_old[i] = pw[i];
      }
      
      // update u
      matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
      normalize_by_l2n(pu, n2);
      
      // update w
      matrix_multiply(pu, ds, 1, n2, p, dtu); // D'u = (u'D)'
      pmax(dtu, p, 0.0);
      for (int i = 0; i < p; i++) {
        pw[i] = soft(dtu[i], lambda);
      }
      normalize_by_l2n(pw, p);
      
      // check convergence
      if (diff_l1n(w_old, pw, p)/l1n(w_old, p) < 1e-4) {
        break;
      }
    }
    
    // update u
    matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
    normalize_by_l2n(pu, n2);
    
    // crit = u'Dw
    matrix_multiply(ds, pw, n2, p, 1, dw);
    crit[l] = cross_prod(pu, dw, n2);
  }
  if(!silent) Rcpp::Rcout << std::endl;
  
  // prepare results
  NumericVector u_(L*n2), w_(L*p), crit_(L);
  IntegerVector iter_(L);
  
  for (int i = 0; i < L*n2; i++)
    u_[i] = u[i];
  for (int i = 0; i < L*p; i++)
    w_[i] = w[i];
  for (int i = 0; i < L; i++)
    crit_[i] = crit[i];
  for (int i = 0; i < L; i++)
    iter_[i] = iter[i];
  
  List res;
  res["u"] = u_;
  res["w"] = w_;
  res["crit"] = crit_;
  res["iter"] = iter_;
  
  // cleanup
  delete [] ds;
  delete [] lambdas;
  
  delete [] u;
  delete [] w;
  delete [] crit;
  delete [] iter;
  
  delete [] w_old;
  delete [] dtu;
  delete [] dw;
  
  return res;
}


// Fused Lasso ------
// [[Rcpp::export]]
List SHC_FL_getuw(NumericMatrix ds_,
                  NumericVector lambda1s_,
                  NumericVector lambda2s_,
                  int max_iter,
                  bool init_random,
                  bool silent) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int pm1 = p - 1;
  int Len1 = lambda1s_.length();
  int Len2 = lambda2s_.length();
  int L = Len1*Len2;
  //int n = int(sqrt(2 * n2)) + 1;
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *lambda1s = new double[Len1];
  for (int i = 0; i < Len1; i++)
    lambda1s[i] = lambda1s_[i];
  
  double *lambda2s = new double[Len2];
  for (int i = 0; i < Len2; i++)
    lambda2s[i] = lambda2s_[i];
  
  double lambda1, lambda2;
  
  // create containers for R output
  double *u = new double[L*n2];
  double *w = new double[L*p];
  double *crit = new double[L];
  int *iter = new int[L];
  int *flsa_converged = new int[L];
  
  // init
  double *pu = NULL;
  
  if (init_random) {
    NumericVector rand_numbers = Rcpp::runif(L*p);
    for (int i = 0; i < L*p; i++) {
      w[i] = rand_numbers[i]*2.0/sqrt(p);
    }
  } else {
    for (int i = 0; i < L*p; i++) {
      w[i] = 1.0/sqrt(p);
    }
  }
  
  double *pw = NULL;
  
  double *w_old = new double[p];
  
  // flsa
  double *z = new double[L*pm1];
  for (int i = 0; i < L*pm1; i++) {
    z[i] = 0.0;
  }
  double *pz = NULL;
  
  double *z0 = new double[pm1];
  for (int i = 0; i < pm1; i++) {
    z0[i] = 0.0;
  }
  
  int n_infor = 4;
  double *infor = new double[L*n_infor];
  for (int i = 0; i < L*n_infor; i++) {
    infor[i] = 0.0;
  }
  double *p_infor = NULL;
  
  // parameters for sfa in flsa
  int maxStep = 1000;
  double tol = 1e-8;
  int tau = 1;
  
  double *dw = new double[n2];
  
  for (int i = 0; i < L; i++) {
    iter[i] = 0;
  }
  
  for (int i = 0; i < L; i++) {
    flsa_converged[i] = 0;
  }
  
  double *dtu = new double[p];
  
  
  // Regularization Path
  for (int l1 = 0; l1 < Len1; l1++) {
    for (int l2 = 0; l2 < Len2; l2++) {
      int l = l1*Len2 + l2;
      Rcpp::checkUserInterrupt();
      
      if(!silent) Rcpp::Rcout << l+1;
      R_FlushConsole();
      
      lambda1 = lambda1s[l1];
      lambda2 = lambda2s[l2];
      // printf("lambda1 = %f\n", lambda1);
      // printf("lambda2 = %f\n", lambda2);
      
      // set pointers
      pu = u + l*n2;
      pw = w + l*p;
      pz = z + l*pm1;
      p_infor = infor + l*n_infor;
      
      /*
      // warm start
      if (l != 0) {
      if (l % Len2 != 0) { // same lambda1 as previous run
      for (int i = 0; i < p; i++) {
      w[l*p + i] = w[(l-1)*p + i];
      }
      } else { // start over to the first value in lambda1s
      for (int i = 0; i < p; i++) {
      w[l*p + i] = w[(l-Len2)*p + i];
      }
      }
      }
      
      // if w are all zeros due to warm start, no iteration wil be performed
      
      The warm start makes flsa be numerically unstable and the program crushes undeterministically.
      Thus it is not used in this case.
      */
      
      // main loop
      while (l1n(pw, p) > 1e-10 && iter[l] < max_iter) {
        iter[l]++;
        
        for (int i = 0; i < p; i++) {
          w_old[i] = pw[i];
        }
        
        // update u
        matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
        normalize_by_l2n(pu, n2);
        
        // update D'u
        matrix_multiply(pu, ds, 1, n2, p, dtu); // D'u = (u'D)'
        
        // z0 is uded for warm start in flsa
        for (int i = 0; i < pm1; i++) {
          z0[i] = pz[i];
        }
        
        // printf("l1n w dtu z0: %f %f %f\n", l1n(pw, p), l1n(dtu, p), l1n(z0, pm1));
        // R_FlushConsole();
        
        // update w and z
        flsa(pw, pz, p_infor, dtu, z0, lambda1, lambda2, p, maxStep, tol, tau);
        
        // update w due to squared l2-norm constraint
        normalize_by_l2n(pw, p);
        
        // check convergence
        if (diff_l1n(w_old, pw, p)/l1n(w_old, p) < 1e-4) {
          break;
        }
      }
      
      // update u
      matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
      normalize_by_l2n(pu, n2);
      
      // crit = u'Dw
      matrix_multiply(ds, pw, n2, p, 1, dw);
      crit[l] = cross_prod(pu, dw, n2);
      
      // convergence status of the last flsa
      flsa_converged[l] = (p_infor[1] != maxStep);
    }
  }
  if(!silent) Rcpp::Rcout << std::endl;
  
  
  // ////////////////////////
  // printf("\n*Debugging here*\n");
  // SEXP tmp_;
  // PROTECT(tmp_ = allocVector(REALSXP, 1));
  // REAL(tmp_)[0] = 1;
  // UNPROTECT(1);
  // return tmp_;
  // ////////////////////////
  
  
  // prepare results
  NumericVector u_(L*n2), w_(L*p), crit_(L);
  IntegerVector iter_(L), flsa_converged_(L);
  
  for (int i = 0; i < L*n2; i++)
    u_[i] = u[i];
  for (int i = 0; i < L*p; i++)
    w_[i] = w[i];
  for (int i = 0; i < L; i++)
    crit_[i] = crit[i];
  for (int i = 0; i < L; i++)
    iter_[i] = iter[i];
  for (int i = 0; i < L; i++)
    flsa_converged_[i] = flsa_converged[i];
  
  List res;
  res["u"] = u_;
  res["w"] = w_;
  res["crit"] = crit_;
  res["iter"] = iter_;
  res["flsa_converged"] = flsa_converged_;
  
  // cleanup
  delete [] ds;
  delete [] lambda1s;
  delete [] lambda2s;
  
  delete [] u;
  delete [] w;
  delete [] crit;
  delete [] iter;
  delete [] flsa_converged;
  
  delete [] z;
  delete [] z0;
  delete [] infor;
  
  delete [] w_old;
  delete [] dtu;
  delete [] dw;
  
  return res;
}


// Spike-and-slab Lasso ------
// [[Rcpp::export]]
List SHC_SSL_getuw(NumericMatrix ds_,
                   CharacterVector penalty_,
                   double lambda1,
                   NumericVector lambda0s_,
                   double theta,
                   double aa,
                   double bb,
                   double eps,
                   int max_iter,
                   bool init_random,
                   bool silent,
                   bool warm_start) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int L = lambda0s_.length();
  //int n = int(sqrt(2 * n2)) + 1;
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *lambda0s = new double[L];
  for (int i = 0; i < L; i++)
    lambda0s[i] = lambda0s_[i];
  
  const char *penalty = penalty_[0];
  
  double lambda0;
  
  // create containers for R output
  double *u = new double[L*n2];
  double *w = new double[L*p];
  double *crit = new double[L];
  int *iter = new int[L];
  double *thetas = new double[L];
  double *deltas = new double[L];
  
  // init
  double *pu = NULL;
  
  if (init_random) {
    NumericVector rand_numbers = Rcpp::runif(L*p);
    for (int i = 0; i < L*p; i++) {
      w[i] = rand_numbers[i]*2.0/sqrt(p);
    }
  } else {
    for (int i = 0; i < L*p; i++) {
      w[i] = 1.0/sqrt(p);
    }
  }
  
  double *pw = NULL;
  
  double *w_old = new double[p];
  
  double *dw = new double[n2];
  
  for (int i = 0; i < L; i++) {
    iter[i] = 0;
  }
  
  double delta = 0;
  
  // dtu before soft and normalization
  double *dtu = new double[p];
  for (int j=0; j<p; j++) {
    dtu[j] = 0;
  }
  
  int converged = 0;
  //int counter = 0;
  //int violations = 0;
  
  // Regularization Path
  for (int l = 0; l < L; l++) {
    Rcpp::checkUserInterrupt();
    
    if(!silent) Rcpp::Rcout << l+1;
    R_FlushConsole();
    
    lambda0 = lambda0s[l];
    // printf("lambda1 = %f\n", lambda1);
    // printf("lambda0 = %f\n", lambda0);
    
    // set pointers
    pu = u + l*n2;
    pw = w + l*p;
    
    // warm start
    if (warm_start) {
      if (l != 0) {
        for (int i = 0; i < p; i++) {
          w[l*p + i] = w[(l-1)*p + i];
        }
      }
    }
    
    // calculate theta and delta
    if (strcmp(penalty, "adaptive")==0) {
      theta = expectation_approx(w, aa, bb, p, l);
    }
    delta = threshold(theta, lambda1, lambda0);
    
    // if w are all zeros due to warm start, no iteration wil be performed
    
    // main loop
    while (l1n(pw, p) > 1e-10 && iter[l] < max_iter) {
      iter[l]++;
      
      for (int j = 0; j < p; j++) {
        w_old[j] = pw[j];
      }
      
      // printf("theta = %f\n", theta);
      // printf("delta = %f\n", delta);
      
      // update u
      matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
      normalize_by_l2n(pu, n2);
      
      // update w
      // we can't use coordinate-wise update b/c w need to be updated simutaneously
      matrix_multiply(pu, ds, 1, n2, p, dtu); // D'u = (u'D)'
      for (int j=0; j<p; j++) {
        pw[j] = SSL(dtu[j], w_old[j], lambda0, lambda1, theta, delta);
      }
      normalize_by_l2n(pw, p);
      
      // update theta and delta
      if(strcmp(penalty, "adaptive")==0) {
        theta = expectation_approx(w, aa, bb, p, l);
        delta = threshold(theta, lambda1, lambda0);
      }
      
      // check for convergence
      converged = checkConvergence(pw, w_old, eps, p);
      
      if (converged) {
        break;
      }
    }
    
    // update u
    matrix_multiply(ds, pw, n2, p, 1, pu); // u = Dw
    normalize_by_l2n(pu, n2);
    
    // crit = u'Dw
    matrix_multiply(ds, pw, n2, p, 1, dw);
    crit[l] = cross_prod(pu, dw, n2);
    
    // out theta and delta
    thetas[l] = theta;
    deltas[l] = delta;
  } // end of regularization loop
  if(!silent) Rcpp::Rcout << std::endl;
  
  // prepare results
  NumericVector u_(L*n2), w_(L*p), crit_(L), thetas_(L), deltas_(L);
  IntegerVector iter_(L);
  
  for (int i = 0; i < L*n2; i++)
    u_[i] = u[i];
  for (int i = 0; i < L*p; i++)
    w_[i] = w[i];
  for (int i = 0; i < L; i++)
    crit_[i] = crit[i];
  for (int i = 0; i < L; i++)
    iter_[i] = iter[i];
  for (int i = 0; i < L; i++)
    thetas_[i] = thetas[i];
  for (int i = 0; i < L; i++)
    deltas_[i] = deltas[i];
  
  List res;
  res["u"] = u_;
  res["w"] = w_;
  res["crit"] = crit_;
  res["iter"] = iter_;
  res["thetas"] = thetas_;
  res["deltas"] = deltas_;
  
  // cleanup
  delete [] ds;
  delete [] lambda0s;
  
  delete [] u;
  delete [] w;
  delete [] crit;
  delete [] iter;
  delete [] thetas;
  delete [] deltas;
  
  delete [] w_old;
  delete [] dtu;
  delete [] dw;
  
  return res;
}


// Get crit given D and w
// [[Rcpp::export]]
List SHC_get_crit(NumericMatrix ds_,
                  NumericVector w_) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int L = w_.length() / p;
  //int n = int(sqrt(2 * n2)) + 1;
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *w = new double[L*p];
  for (int i = 0; i < L*p; i++)
    w[i] = w_[i];
  
  // create containers for R output
  double *crit = new double[L];
  
  double *u = new double[L*n2];
  double *dw = new double[n2];
  
  double *pw = NULL;
  
  // Regularization Path
  for (int l = 0; l < L; l++) {
    Rcpp::checkUserInterrupt();
    
    pw = w + l*p;
    
    // update u
    matrix_multiply(ds, pw, n2, p, 1, u); // u = Dw
    normalize_by_l2n(u, n2);
    
    // crit = u'Dw
    matrix_multiply(ds, pw, n2, p, 1, dw);
    crit[l] = cross_prod(u, dw, n2);
  }
  
  // prepare results
  NumericVector crit_(L);
  
  for (int i = 0; i < L; i++)
    crit_[i] = crit[i];
  
  List res;
  res["crit"] = crit_;
  
  // cleanup
  delete [] ds;
  delete [] w;
  
  delete [] crit;
  
  delete [] u;
  delete [] dw;
  
  return res;
}


// Get u given D and w
// [[Rcpp::export]]
List SHC_get_u(NumericMatrix ds_,
               NumericVector w_) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  
  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];
  
  double *w = new double[p];
  for (int i = 0; i < p; i++)
    w[i] = w_[i];
  
  // create containers for R output
  double *u = new double[n2];
  
  // update u
  matrix_multiply(ds, w, n2, p, 1, u); // u = Dw
  normalize_by_l2n(u, n2);
  
  // prepare results
  NumericVector u_(n2);
  
  for (int i = 0; i < n2; i++)
    u_[i] = u[i];
  
  List res;
  res["u"] = u_;
  
  // cleanup
  delete [] ds;
  delete [] w;
  
  delete [] u;
  
  return res;
}


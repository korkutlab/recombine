#include <math.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;


// convert ds matrix (n2*p) and var weights to distance matrix (n*n, symmetric)
void ds2dist(double *ds,
             double *w,
             double *dist,
             int n,
             int p) {
  int n2 = n*(n - 1)/2;

  // calculation
  int ii = 0;
  double d;
  for (int i = 0; i < n-1; i++) {
    dist[i*n + i] = 0.0; // diagnal is 0
    for (int ip = i+1; ip < n; ip++) {
      d = 0.0;
      for (int j = 0; j < p; j++)
        d += ds[ii + j*n2]*w[j]; // R matrix is column-major internally
      dist[i*n + ip] = d;
      dist[ip*n + i] = d; // dist is symmetric
      ii += 1;
    }
  }
}


// convert ds matrix (n2*p) to distance matrix (n*n, symmetric)
void ds2dist_unweighted(double *ds,
                        double *dist,
                        int n,
                        int p) {
  int n2 = n*(n - 1)/2;

  // calculation
  int ii = 0;
  double d;
  for (int i = 0; i < n-1; i++) {
    dist[i*n + i] = 0.0; // diagnal is 0
    for (int ip = i+1; ip < n; ip++) {
      d = 0.0;
      for (int j = 0; j < p; j++)
        d += ds[ii + j*n2]; // R matrix is column-major internally
      dist[i*n + ip] = d;
      dist[ip*n + i] = d; // dist is symmetric
      ii += 1;
    }
  }
}


// get kNN distances and ids from distance matrix
void kNN(double *dist,
         int n,
         int k,
         double *dist_knn,
         int *id_knn) {
  for (int ii = 0; ii < n; ii++) {
    // build container
    // A priority queue is a container adaptor that provides constant time lookup of the largest (by default) element.
    // A user-provided Compare can be supplied to change the ordering, e.g. using std::greater<T> would cause the smallest element to appear as the top().
    std::priority_queue<std::pair<double, int>,
                        std::vector<std::pair<double, int>>,
                        std::greater<std::pair<double, int>>> q;
    for (int i = 0; i < n; ++i) {
      if (i != ii) q.push(std::pair<double, int>(dist[ii*n + i], i));
    }

    // get knn distance and index
    for (int i = 0; i < k; ++i) {
      dist_knn[ii*k + i] = q.top().first;
      id_knn[ii*k + i] = q.top().second;
      q.pop();
    }
  }
}


// get local outlier factors from distance matrix
// lof depends on actual values of distance matrix
// square.distance, absolute.value, and eucledian lead to different lof for same data
void get_lof(double *dist,
             int n,
             int k,
             double *lof) {
  double *dist_knn = new double[n*k];
  int *id_knn = new int[n*k];
  double *lrd = new double[n];

  // get k nearest neighbors + distances
  kNN(dist, n, k, dist_knn, id_knn);

  // std::cout << " " << std::endl;
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < k; ++j)
  //     std::cout << dist_knn[i*k + j] << " ";
  //   std::cout << std::endl;
  // }
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < k; ++j)
  //     std::cout << id_knn[i*k + j] << " ";
  //   std::cout << std::endl;
  // }

  double d;

  // calculate local reachability density
  for (int i = 0; i < n; ++i) {
    d = 0.0;
    // walk through kNN
    for (int j = 0; j < k; ++j) {
      int ip = id_knn[i*k + j];
      d += std::max(dist_knn[i*k + j], dist_knn[ip*k + (k-1)]); // reachability distance of object i w.r.t. object ip
    }
    lrd[i] = k / d;
  }

  // calculate lof
  for (int i = 0; i < n; ++i) {
    d = 0.0;
    // walk through kNN
    for (int j = 0; j < k; ++j) {
      int ip = id_knn[i*k + j];
      d += lrd[ip];
    }
    lof[i] = d / k / lrd[i];
  }

  // for (int i = 0; i < n; ++i) {
  //   std::cout << lrd[i] << " ";
  // }
  // std::cout << std::endl;
  // for (int i = 0; i < n; ++i) {
  //     std::cout << lof[i] << " ";
  // }
  // stop("test");

  delete [] dist_knn;
  delete [] id_knn;
  delete [] lrd;
}


// get median of 1D array
double median(double *arr,
              int n) {
  if(n == 0) {
    return 0.0;
  }

  // std::nth_element and std::max_element has post-condition:
  // After calling, the elements in vector may be reordered and the resulting order is implementation defined.
  std::vector<double> v;
  for (int i = 0; i < n; ++i)
    v.push_back(arr[i]);

  auto n_med = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n_med, v.end());
  auto med = v[n_med];
  if(!(v.size() & 1)) { //If the set size is even
    auto max_it = std::max_element(v.begin(), v.begin()+n);
    med = (*max_it + med) / 2.0;
  }

  return med;
}


// scale 1D array to have zero mean and unit variance
void scale(double *arr,
           int n) {
  if(n < 2) {
    return;
  }

  double mean = 0.0;
  for (int i = 0; i < n; ++i) {
    mean += arr[i];
  }
  mean /= n;

  double sd = 0.0;
  for (int i = 0; i < n; ++i) {
    arr[i] -= mean;
    sd += arr[i]*arr[i];
  }
  sd = sqrt(sd / (n-1));

  for (int i = 0; i < n; ++i) {
    arr[i] /= sd;
  }
}


// get local outlier probability from distance matrix
// loop depends on actual values of distance matrix
// square.distance, absolute.value, and eucledian lead to different loop for same data
void get_loop(double *dist,
              int n,
              int k,
              double lambda,
              double *loop) {
  double *dist_knn = new double[n*k];
  int *id_knn = new int[n*k];
  double *pdist = new double[n];
  double *plof = new double[n];

  // get k nearest neighbors + distances
  kNN(dist, n, k, dist_knn, id_knn);

  // std::cout << " " << std::endl;
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < k; ++j)
  //     std::cout << dist_knn[i*k + j] << " ";
  //   std::cout << std::endl;
  // }
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < k; ++j)
  //     std::cout << id_knn[i*k + j] << " ";
  //   std::cout << std::endl;
  // }

  // calculate probabilistic set distance of kNN with significance lambda
  for (int i = 0; i < n; ++i) {
    double d = 0.0;
    // walk through kNN
    for (int j = 0; j < k; ++j) {
      d += dist_knn[i*k + j]*dist_knn[i*k + j];
    }

    pdist[i] = lambda*sqrt(d/k);
  }

  // calculate probabilistic local outlier factor
  for (int i = 0; i < n; ++i) {
    // calculate expected pdist of kNN
    double epdist = 0.0;
    // walk through kNN
    for (int j = 0; j < k; ++j) {
      int ip = id_knn[i*k + j];
      epdist += pdist[ip];
    }
    epdist /= k;

    plof[i] = pdist[i]/epdist - 1;
  }

  // calculate nplof
  double nplof = 0.0;
  for (int i = 0; i < n; ++i) {
    nplof += plof[i]*plof[i];
  }
  nplof = lambda*sqrt(nplof/n);

  // calculate loop
  for (int i = 0; i < n; ++i) {
    double d = erf(plof[i]/(sqrt(2.0)*nplof));
    if (d > 0)
      loop[i] = d;
    else
      loop[i] = 0;
  }

  // for (int i = 0; i < n; ++i) {
  //     std::cout << loop[i] << " ";
  // }
  // stop("test");

  delete [] dist_knn;
  delete [] id_knn;
  delete [] pdist;
  delete [] plof;
}


// get outlier indicators from distance matrix
// distance matrix is obtained iternally from ds matrix (n2*p) and var weights
void get_outliers(double *dist,
                  int n,
                  int loop_k,
                  double loop_lambda,
                  double loop_threshold,
                  bool outlier_on,
                  double *loop,
                  int *is_outlier) {
  // reset is_outlier
  for (int i = 0; i < n; ++i)
    is_outlier[i] = 0;

  if (!outlier_on)
    return;

  // get loop
  get_loop(dist, n, loop_k, loop_lambda, loop);

  // assign outlier based on loop_threshold
  for (int i = 0; i < n; ++i)
    if (loop[i] > loop_threshold) is_outlier[i] = 1;
}


// An R wrapper of get_outliers() given distance matrix
// [[Rcpp::export]]
List get_outliers_from_dist(NumericMatrix dist_,
                            int loop_k,
                            double loop_lambda,
                            double loop_threshold,
                            bool outlier_on) {
  // Declarations
  int n = dist_.nrow();

  // copy input vectors
  double *dist = new double[n*n];
  for (int i = 0; i < n*n; i++) {
    dist[i] = dist_[i];
  }

  // create containers for R output
  int *is_outlier = new int[n];
  double *loop = new double[n];

  // get outlier indicators
  get_outliers(dist, n,
               loop_k,
               loop_lambda,
               loop_threshold,
               outlier_on,
               loop,
               is_outlier);

  // prepare results
  IntegerVector is_outlier_(n);
  NumericVector loop_(n);

  for (int i = 0; i < n; i++)
    is_outlier_[i] = is_outlier[i];
  for (int i = 0; i < n; i++)
    loop_[i] = loop[i];

  List res;
  res["is_outlier"] = is_outlier_;
  res["loop"] = loop_;

  // cleanup
  delete [] dist;
  delete [] is_outlier;
  delete [] loop;

  return res;
}


// An R wrapper of get_outliers() given ds matrix
// [[Rcpp::export]]
List get_outliers_from_ds(NumericMatrix ds_,
                          int loop_k,
                          double loop_lambda,
                          double loop_threshold,
                          bool outlier_on) {
  // Declarations
  int n2 = ds_.nrow();
  int p = ds_.ncol();
  int n = int(sqrt(2 * n2)) + 1;

  // copy input vectors
  double *ds = new double[n2*p];
  for (int i = 0; i < n2*p; i++)
    ds[i] = ds_[i];

  // get distance matrix (symmetric)
  // ds_all is unmasked ds
  double *dist = new double[n*n];
  ds2dist_unweighted(ds, dist, n, p);

  // create containers for R output
  int *is_outlier = new int[n];
  double *loop = new double[n];

  // get outlier indicators
  get_outliers(dist, n,
               loop_k,
               loop_lambda,
               loop_threshold,
               outlier_on,
               loop,
               is_outlier);

  // prepare results
  IntegerVector is_outlier_(n);
  NumericVector loop_(n);

  for (int i = 0; i < n; i++)
    is_outlier_[i] = is_outlier[i];
  for (int i = 0; i < n; i++)
    loop_[i] = loop[i];

  List res;
  res["is_outlier"] = is_outlier_;
  res["loop"] = loop_;

  // cleanup
  delete [] ds;
  delete [] dist;

  delete [] is_outlier;
  delete [] loop;

  return res;
}


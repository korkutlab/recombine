getuw_ssl <- function(ds,
                      penalty = c("adaptive", "separable"),
                      lambda0s = NULL,
                      nlambda0 = 100,
                      lambda1 = NULL,
                      theta = 0.5,
                      a = 1,
                      b,
                      eps = 0.001,
                      max.iter = 500,
                      iseed = NULL,
                      silent = FALSE,
                      warm_start = TRUE) {
  penalty <- match.arg(penalty)
  
  # dim
  n2 <- nrow(ds)
  p <- ncol(ds)
  n <- ceiling(sqrt(2 * n2)) # n2 = n(n-1)/2
  
  if (is.null(lambda0s)) {
    lambda0s <- seq(1, n, length = nlambda0)
    lambda1 <- lambda0s[1]
  } else {
    nlambda0 <- length(lambda0s)
    if (is.null(lambda1)) {
      lambda1 <- lambda0s[1]
    }
  }
  
  # Lambda0s should be an increasing sequence
  monotone <- sum((lambda0s[-1] - lambda0s[-nlambda0]) > 0)
  if (monotone != nlambda0 - 1){
    stop("lambda0s must be a monotone increasing sequence")
  }
  if (lambda1 > min(lambda0s) ) {
    stop("lambda1 must be smaller than lambda0s")
  }
  
  if(missing(b)) {
    b <- p
  }
  
  if (is.null(iseed)) {
    init_random <- FALSE
    iseed <- 42L
  } else {
    init_random <- TRUE
  }
  
  # call
  set.seed(iseed)
  res <- SHC_SSL_getuw(ds,
                       penalty,
                       as.double(lambda1),
                       as.numeric(lambda0s),
                       as.double(theta),
                       as.double(a),
                       as.double(b),
                       eps,
                       as.integer(max.iter),
                       init_random,
                       silent,
                       warm_start)
  
  u <- res[["u"]]
  w <- res[["w"]]
  crit <- res[["crit"]]
  iter <- res[["iter"]]
  
  # # check convergence
  # if (any(iter == max.iter)) {
  #   message("Algorithm did not converge for the following regularization values")
  #   print(lambda0s[iter == max.iter])
  # }
  #
  # if (iter[nlambda0] == max.iter) {
  #   message("Algorithm did not converge at the last regularization value")
  # }
  
  # postprocess u
  u_list <- list()
  for (i in c(1:nlambda0)) {
    u_list[[i]] <- post_u(u[((i-1)*n2+1):(i*n2)])
  }
  
  # postprocess w
  w_list <- list()
  for (i in c(1:nlambda0)) {
    w_list[[i]] <- w[((i-1)*p+1):(i*p)]
  }
  
  w_l0norm <- sapply(w_list, l0norm)
  nonzero_w_indices_list <- lapply(w_list, function(x) which(x > 0))
  
  return(list(lassotype = 'SSL',
              lassotype_full = 'Spike-and-slab LASSO',
              u = u_list, w = w_list,
              crit = crit, iter = iter,
              w_l0norm = w_l0norm,
              nonzero_w_indices = nonzero_w_indices_list,
              hyps = lambda0s))
}


#' orSHC_SSL
#'
#' @description
#' Outlier-robust sparse hierarchical clustering with spike-and-slab lasso
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param penalty type of spike-and-slab lasso penalty. This should be either "adaptive" (with a random theta, where theta ~ B(a,b)) or "separable" (with a fixed theta). The default is "adaptive".
#' @param lambda0 spike variance parameter. 
#' @param lambda1 slab variance parameter. Note lambda1 should be no greater than lambda0. 
#' @param theta prior mixing proportion. For "separable" penalty, this value is fixed. For "adaptive" penalty, this value is used as a starting value. The default is 0.5.
#' @param a parameter of the beta prior B(a,b) for the adaptive penalty. The default is 1.
#' @param b parameter of the beta prior B(a,b) for the adaptive penalty. The default is \code{ncol(x)}.
#' @param eps convergence criterion. Converged when all component-wise differences in w are less than eps. The default is 0.001.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param loop_k size of k-nearest neighbors defining local contexts in local outlier probability (LoOP) calculation. The default is 20.
#' @param loop_lambda scaling parameter in local outlier probability (LoOP) calculation. The default is 3.
#' @param loop_threshold threshold of LoOP above which LoOP values imply outliers. The default is 0.5.
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#'
#' @details
#' This function, orSHC_SSL, which robustifies SHC_SSL, is robust to the existence of outliers. 
#' It trims outliers during feature selection and perform final clustering using all samples. 
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{u}{the optimal weighted dissimilarity matrix.}
#' \item{w}{the optimal weights.}
#' \item{crit}{the optimal objective.}
#' \item{iter}{number of iteration steps.}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{nonzero_w_indices}{indices of nonzero w's.}
#' \item{hyp}{hyperparameter used (i.e., lambda0).}
#' \item{is_outlier}{outlier indicator of each sample.}
#' \item{loop}{LoOP value of each sample.}
#' \item{hc}{an object of class hclust by running hierarchical clustering on u.}
#' 
#' @examples
#' source(system.file("scripts", "gen_sim_data.R", package = "recombine"))
#' d <- gen_sim_data(out_pct = 0.1, iseed = 1)
#' result <- orSHC_SSL(d$x,
#'                     lambda0 = 200,
#'                     lambda1 = 0.001)
#' result
#' 
#' @export
#'
orSHC_SSL <- function(
  x,
  penalty = c("adaptive", "separable"),
  lambda0,
  lambda1,
  theta = 0.5,
  a = 1,
  b,
  eps = 0.001,
  standardize.arrays = FALSE,
  dissimilarity = c("squared.distance",
                    "absolute.value"),
  method = c("average",
             "complete",
             "single",
             "centroid",
             "ward.D",
             "ward.D2",
             "mcquitty",
             "median"),
  loop_k = 20,
  loop_lambda = 3,
  loop_threshold = 0.5,
  max.iter = 500,
  iseed = NULL,
  silent = FALSE)
{
  lassotype <- "ssl"
  dissimilarity <- match.arg(dissimilarity)
  method <- match.arg(method)
  if (is.null(iseed)) iseed <- as.integer(Sys.time())
  
  # determine outlier flags
  if (nrow(x) < loop_k + 1) stop("loop_k needs to be smaller than x")
  if (!(loop_lambda > 0)) stop("loop_lambda needs to be positive")
  
  outlier_on <- !is.null(loop_threshold)
  
  if (outlier_on) {
    if (!(loop_threshold > 0 &
          loop_threshold < 1))
      stop("loop_threshold needs to be positive")
  }
  
  if (standardize.arrays) {
    x <- sweep(x,2,apply(x,2,mean,na.rm=TRUE),"-")
    x <- sweep(x,2,apply(x,2,sd,na.rm=TRUE),"/")
  }
  
  if(missing(b)) {
    b <- ncol(x)
  }
  
  ds <- x2ds_nona(x, dissimilarity=dissimilarity)
  
  # step 1: remove outliers
  if (outlier_on) {
    res <- get_outliers_from_ds(ds,
                                as.integer(loop_k),
                                as.double(loop_lambda),
                                as.double(loop_threshold),
                                outlier_on)
    
    is_outlier <- res[["is_outlier"]]
    is_outlier <- as.logical(is_outlier)
    
    loop <- res[["loop"]]
    
    # remove outliers
    ds_all  <- ds
    ds <- x2ds_nona(x[!is_outlier, ], dissimilarity=dissimilarity)
  } else {
    ds_all <- NULL
    is_outlier <- NULL
    loop = NULL
  }
  
  # step 2: run sparse clustering
  out <- getuw_ssl(ds = ds,
                   penalty = penalty,
                   lambda0s = lambda0,
                   lambda1 = lambda1,
                   theta = theta,
                   a = a,
                   b = b,
                   eps = eps,
                   max.iter = max.iter,
                   iseed = iseed,
                   silent = silent)
  
  # step 3: get u for all data using outlier-free w
  # wrap out
  class_name <- "orSHC_SSL"
  result <- extract_unperm_out(index = 1L,
                               unperm_outs = out,
                               ds = ds_all,
                               is_outlier,
                               loop,
                               dissimilarity = dissimilarity,
                               method = method,
                               labels = rownames(x),
                               class_name = class_name)
  
  return(result)
}


#' SHC_SSL
#'
#' @description
#' Sparse hierarchical clustering with spike-and-slab lasso
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param penalty type of spike-and-slab lasso penalty. This should be either "adaptive" (with a random theta, where theta ~ B(a,b)) or "separable" (with a fixed theta). The default is "adaptive".
#' @param lambda0 spike variance parameter. 
#' @param lambda1 slab variance parameter. Note lambda1 should be no greater than lambda0. 
#' @param theta prior mixing proportion. For "separable" penalty, this value is fixed. For "adaptive" penalty, this value is used as a starting value. The default is 0.5.
#' @param a parameter of the beta prior B(a,b) for the adaptive penalty. The default is 1.
#' @param b parameter of the beta prior B(a,b) for the adaptive penalty. The default is \code{ncol(x)}.
#' @param eps convergence criterion. Converged when all component-wise differences in w are less than eps. The default is 0.001.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#'
#' @details
#' SHC_SSL employs a spike-and-slab lasso penalty to approach the theoretically ideal feature selection coupled with clustering. 
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{u}{the optimal weighted dissimilarity matrix.}
#' \item{w}{the optimal weights.}
#' \item{crit}{the optimal objective.}
#' \item{iter}{number of iteration steps.}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{nonzero_w_indices}{indices of nonzero w's.}
#' \item{hyp}{hyperparameter used (i.e., lambda0).}
#' \item{hc}{an object of class hclust by running hierarchical clustering on u.}
#' 
#' @examples
#' d <- example_sim_data
#' result <- SHC_SSL(d$x,
#'                   lambda0 = 200,
#'                   lambda1 = 0.001)
#' result
#' 
#' @export
#'
SHC_SSL <- function(
  x,
  penalty = c("adaptive", "separable"),
  lambda0,
  lambda1,
  theta = 0.5,
  a = 1,
  b,
  eps = 0.001,
  standardize.arrays = FALSE,
  dissimilarity = c("squared.distance",
                    "absolute.value"),
  method = c("average",
             "complete",
             "single",
             "centroid",
             "ward.D",
             "ward.D2",
             "mcquitty",
             "median"),
  max.iter = 500,
  iseed = NULL,
  silent = FALSE)
{
  if(missing(b)) {
    b <- ncol(x)
  }
  
  result <- orSHC_SSL(
    x = x,
    penalty = penalty,
    lambda0 = lambda0,
    lambda1 = lambda1,
    theta = theta,
    a = a,
    b = b,
    eps = eps,
    standardize.arrays = standardize.arrays,
    dissimilarity = dissimilarity,
    method = method,
    loop_threshold = NULL,
    max.iter = max.iter,
    iseed = iseed,
    silent = silent)
  
  result$is_outlier <- NULL
  result$loop <- NULL
  result$class_name <- "SHC_SSL"
  
  return(result)
}


#' orSHC_SSL_gapstat
#'
#' @description
#' This function computes gap statistic profile for a series of hyperparameters (lambda0's) by running orSHC_SSL on both original data and permuted data.
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param penalty type of spike-and-slab lasso penalty. This should be either "adaptive" (with a random theta, where theta ~ B(a,b)) or "separable" (with a fixed theta). The default is "adaptive".
#' @param lambda0s a vector of spike variance parameters. If Null, lambda0s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 100)}. The default is NULL.
#' @param lambda1 slab variance parameter. Note lambda1 should be no greater than lambda0. If Null, lambda1 is assigned as the smallest value of lambda0s. The default is NULL.
#' @param theta prior mixing proportion. For "separable" penalty, this value is fixed. For "adaptive" penalty, this value is used as a starting value. The default is 0.5.
#' @param a parameter of the beta prior B(a,b) for the adaptive penalty. The default is 1.
#' @param b parameter of the beta prior B(a,b) for the adaptive penalty. The default is \code{ncol(x)}.
#' @param eps convergence criterion. Converged when all component-wise differences in w are less than eps. The default is 0.001.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param loop_k size of k-nearest neighbors defining local contexts in local outlier probability (LoOP) calculation. The default is 20.
#' @param loop_lambda scaling parameter in local outlier probability (LoOP) calculation. The default is 3.
#' @param loop_threshold threshold of LoOP above which LoOP values imply outliers. The default is 0.5.
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#' @param warm_start “warm start” strategy to accelerate convergence for a series of hyperparameters. Since solutions for similar hyperparameters are close, w is initialized as a solution for a nearby hyperparameter so that its solution can be found relatively quickly. The default is TRUE.
#' @param return_all_results should orSHC_SSL results for all hyperparameters be returned? If False, only the result for the best hyperparameter is returned. The default is FALSE.
#' @param nperms number of permutations. The default is 10.
#' @param sel_rule rule for choosing the best hyperparameter. This should be either "max" (maximum gap statistic) or "1se.rule" (having a low complexity but within 1 standard error of the maximum gap statistic). The default is "max".
#'
#' @details
#' For a hyperparameter, the gap statistic measures the strength of the clustering based on real data with respect to the one based on randomly permuted data that are supposed to have no cluster. 
#' This function runs orSHC_SSL for every hyperparameter, calculates gap statistic, and get the best hyperparamter that corresponds to the maximum gap statistic.
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{tots}{the strength of the clustering using real data (i.e., maximized objective function values).}
#' \item{permtots}{the strength of the clustering using permuted data (i.e., maximized objective function values).}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{gaps_mean}{mean values of gap statistic.}
#' \item{gaps_se}{standard errors of gap statistic.}
#' \item{hyps}{hyperparameters (i.e., lambda0s).}
#' \item{best_hyp}{best hyperparameter (i.e., best lambda0).}
#' \item{best_index}{index of best hyperparameter (i.e., best lambda0).}
#' \item{result}{orSHC_SSL result. If return_all_results is TRUE, results of all hyperparameters are stored as a list, in which each element is a result of orSHC_SSL; otherwise, it only contains the result corresponding to the best hyperparameter.}
#' 
#' @examples
#' source(system.file("scripts", "gen_sim_data.R", package = "recombine"))
#' d <- gen_sim_data(out_pct = 0.1, iseed = 1)
#' out <- orSHC_SSL_gapstat(d$x,
#'                          lambda0s = seq(0.01, 0.5*max(nrow(d$x), ncol(d$x)), length = 20),
#'                          lambda1 = 0.001,
#'                          nperms = 3)
#' plot(out)
#' 
#' @export
#'
orSHC_SSL_gapstat <- function(
  x,
  penalty = c("adaptive", "separable"),
  lambda0s = NULL,
  lambda1 = NULL,
  theta = 0.5,
  a = 1,
  b,
  eps = 0.001,
  standardize.arrays = FALSE,
  dissimilarity = c("squared.distance",
                    "absolute.value"),
  method = c("average",
             "complete",
             "single",
             "centroid",
             "ward.D",
             "ward.D2",
             "mcquitty",
             "median"),
  loop_k = 20,
  loop_lambda = 3,
  loop_threshold = 0.5,
  max.iter = 500,
  iseed = NULL,
  silent = FALSE,
  warm_start = TRUE,
  return_all_results = FALSE,
  nperms = 10,
  sel_rule = c("max", "1se.rule"))
{
  lassotype <- "ssl"
  dissimilarity <- match.arg(dissimilarity)
  if (is.null(iseed)) iseed <- as.integer(Sys.time())
  sel_rule <- match.arg(sel_rule)
  
  # determine outlier flags
  if (nrow(x) < loop_k + 1) stop("loop_k needs to be smaller than x")
  if (!(loop_lambda > 0)) stop("loop_lambda needs to be positive")
  
  outlier_on <- !is.null(loop_threshold)
  
  if (outlier_on) {
    if (!(loop_threshold > 0 &
          loop_threshold < 1))
      stop("loop_threshold needs to be positive")
  }
  
  if (standardize.arrays) {
    x <- sweep(x,2,apply(x,2,mean,na.rm=TRUE),"-")
    x <- sweep(x,2,apply(x,2,sd,na.rm=TRUE),"/")
  }
  
  if(is.null(lambda0s)) lambda0s <- seq(1, max(nrow(x), ncol(x)), length = 100)
  if(length(lambda0s) < 2) stop("Lambda0s should be a vector with at least 2 elements.")
  n_regul <- length(lambda0s)
  
  if(missing(b)) {
    b <- ncol(x)
  }
  
  ds <- x2ds_nona(x, dissimilarity=dissimilarity)
  
  if (outlier_on) {
    res <- get_outliers_from_ds(ds,
                                as.integer(loop_k),
                                as.double(loop_lambda),
                                as.double(loop_threshold),
                                outlier_on)
    
    is_outlier <- res[["is_outlier"]]
    is_outlier <- as.logical(is_outlier)
    
    loop <- res[["loop"]]
    
    # remove outliers
    ds_all <- ds
    ds <- x2ds_nona(x[!is_outlier, ], dissimilarity=dissimilarity)
  } else {
    ds_all <- NULL
    is_outlier <- NULL
    loop = NULL
  }
  
  # run sparse clustering
  cat("Running on unpermuted data", fill=TRUE)
  out <- getuw_ssl(ds = ds,
                   penalty = penalty,
                   lambda0s = lambda0s,
                   lambda1 = lambda1,
                   theta = theta,
                   a = a,
                   b = b,
                   eps = eps,
                   max.iter = max.iter,
                   iseed = iseed,
                   silent = silent,
                   warm_start = warm_start)
  
  lassotype_full <- out$lassotype_full
  w_l0norm <- out$w_l0norm
  tots <- out$crit
  
  permtots <- matrix(NA, nrow=nperms, ncol=n_regul)
  
  cat("Running on permuted data", fill=TRUE)
  permds <- ds
  for(k in 1:nperms) {
    cat("Permutation ", k, " of ", nperms, fill=TRUE)
    # Oooohhhh.. It turns out that rather than permuting the columns of x and then computing a dist matrix, we can simply permute
    #  the columns of the (n choose 2)xp dist matrix.
    for(j in 1:ncol(permds)) permds[,j] <- sample(permds[,j])
    perm.out <- getuw_ssl(ds = permds,
                          penalty = penalty,
                          lambda0s = lambda0s,
                          lambda1 = lambda1,
                          theta = theta,
                          a = a,
                          b = b,
                          eps = eps,
                          max.iter = max.iter,
                          iseed = iseed,
                          silent = silent,
                          warm_start = warm_start)
    permtots[k,] <- perm.out$crit
  }
  
  gaps_mean <- log(tots + 1) - colMeans(log(permtots + 1))
  gaps_se <- apply(log(permtots + 1), 2, sd) / sqrt(nperms)
  
  if (sel_rule == "max") {
    best_index <- get_best_index_at_max(gaps_mean)
  } else {
    best_index <- get_best_index_1se_rule_right(gaps_mean, gaps_se)
  }
  
  # get all results or best result using best hyp
  class_name <- "orSHC_SSL"
  if (return_all_results) {
    result <- lapply(c(1:n_regul),
                     extract_unperm_out,
                     unperm_outs = out,
                     ds = ds_all,
                     is_outlier,
                     loop,
                     dissimilarity = dissimilarity,
                     method = method,
                     labels = rownames(x),
                     class_name = class_name)
  } else {
    result <- extract_unperm_out(index = best_index,
                                 unperm_outs = out,
                                 ds = ds_all,
                                 is_outlier,
                                 loop,
                                 dissimilarity = dissimilarity,
                                 method = method,
                                 labels = rownames(x),
                                 class_name = class_name)
  }
  
  out <- list(lassotype=lassotype,
              lassotype_full=lassotype_full,
              tots=tots,
              permtots=permtots,
              w_l0norm=w_l0norm,
              gaps_mean=gaps_mean,
              gaps_se=gaps_se,
              hyps=lambda0s,
              best_hyp=lambda0s[best_index],
              best_index=best_index,
              result=result)
  
  class(out) <- "orSHC_SSL_gapstat"
  return(out)
}


#' SHC_SSL_gapstat
#'
#' @description
#' This function computes gap statistic profile for a series of hyperparameters (lambda0's) by running SHC_SSL on both original data and permuted data.
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param penalty type of spike-and-slab lasso penalty. This should be either "adaptive" (with a random theta, where theta ~ B(a,b)) or "separable" (with a fixed theta). The default is "adaptive".
#' @param lambda0s a vector of spike variance parameters. If Null, lambda0s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 100)}. The default is NULL.
#' @param lambda1 slab variance parameter. Note lambda1 should be no greater than lambda0. If Null, lambda1 is assigned as the smallest value of lambda0s. The default is NULL.
#' @param theta prior mixing proportion. For "separable" penalty, this value is fixed. For "adaptive" penalty, this value is used as a starting value. The default is 0.5.
#' @param a parameter of the beta prior B(a,b) for the adaptive penalty. The default is 1.
#' @param b parameter of the beta prior B(a,b) for the adaptive penalty. The default is \code{ncol(x)}.
#' @param eps convergence criterion. Converged when all component-wise differences in w are less than eps. The default is 0.001.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#' @param warm_start “warm start” strategy to accelerate convergence for a series of hyperparameters. Since solutions for similar hyperparameters are close, w is initialized as a solution for a nearby hyperparameter so that its solution can be found relatively quickly. The default is TRUE.
#' @param return_all_results should SHC_SSL results for all hyperparameters be returned? If False, only the result for the best hyperparameter is returned. The default is FALSE.
#' @param nperms number of permutations. The default is 10.
#' @param sel_rule rule for choosing the best hyperparameter. This should be either "max" (maximum gap statistic) or "1se.rule" (having a low complexity but within 1 standard error of the maximum gap statistic). The default is "max".
#'
#' @details
#' For a hyperparameter, the gap statistic measures the strength of the clustering based on real data with respect to the one based on randomly permuted data that are supposed to have no cluster. 
#' This function runs SHC_SSL for every hyperparameter, calculates gap statistic, and get the best hyperparamter that corresponds to the maximum gap statistic.
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{tots}{the strength of the clustering using real data (i.e., maximized objective function values).}
#' \item{permtots}{the strength of the clustering using permuted data (i.e., maximized objective function values).}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{gaps_mean}{mean values of gap statistic.}
#' \item{gaps_se}{standard errors of gap statistic.}
#' \item{hyps}{hyperparameters (i.e., lambda0s).}
#' \item{best_hyp}{best hyperparameter (i.e., best lambda0).}
#' \item{best_index}{index of best hyperparameter (i.e., best lambda0).}
#' \item{result}{SHC_SSL result. If return_all_results is TRUE, results of all hyperparameters are stored as a list, in which each element is a result of SHC_SSL; otherwise, it only contains the result corresponding to the best hyperparameter.}
#' 
#' @examples
#' d <- example_sim_data
#' out <- SHC_SSL_gapstat(d$x,
#'                        lambda0s = seq(0.01, 0.5*max(nrow(d$x), ncol(d$x)), length = 20),
#'                        lambda1 = 0.001,
#'                        nperms = 3)
#' plot(out)
#' 
#' @export
#'
SHC_SSL_gapstat <- function(
  x,
  penalty = c("adaptive", "separable"),
  lambda0s = NULL,
  lambda1 = NULL,
  theta = 0.5,
  a = 1,
  b,
  eps = 0.001,
  standardize.arrays = FALSE,
  dissimilarity = c("squared.distance",
                    "absolute.value"),
  method = c("average",
             "complete",
             "single",
             "centroid",
             "ward.D",
             "ward.D2",
             "mcquitty",
             "median"),
  max.iter = 500,
  iseed = NULL,
  silent = FALSE,
  warm_start = TRUE,
  return_all_results = FALSE,
  nperms = 10,
  sel_rule = c("max", "1se.rule"))
{
  if(missing(b)) {
    b <- ncol(x)
  }
  
  out <- orSHC_SSL_gapstat(
    x = x,
    penalty = penalty,
    lambda0s = lambda0s,
    lambda1 = lambda1,
    theta = theta,
    a = a,
    b = b,
    eps = eps,
    standardize.arrays = standardize.arrays,
    dissimilarity = dissimilarity,
    method = method,
    loop_threshold = NULL,
    max.iter = max.iter,
    iseed = iseed,
    silent = silent,
    warm_start = warm_start,
    return_all_results = return_all_results,
    nperms = nperms,
    sel_rule = sel_rule)
  
  if (return_all_results) {
    out$result <- lapply(out$result, 
                         function(result) {
                           result$is_outlier <- NULL
                           result$loop <- NULL
                           result$class_name <- "SHC_SSL"
                           return(result)
                         }
    )
  } else {
    out$result$is_outlier <- NULL
    out$result$loop <- NULL
    out$result$class_name <- "SHC_SSL"
  }
  
  class(out) <- "SHC_SSL_gapstat"
  return(out)
}


#' @export
print.SHC_SSL <- function(x, ...) {
  print.SHC_base(x, ...)
}


#' @export
print.SHC_SSL_gapstat <- function(x, ...) {
  print.SHC_base_gapstat(x, ...)
}


#' @export
plot.SHC_SSL_gapstat <- function(x, ...) {
  plot.SHC_base_gapstat(x, ...)
}


#' @export
print.orSHC_SSL <- function(x, ...) {
  print.SHC_base(x, ...)
}


#' @export
print.orSHC_SSL_gapstat <- function(x, ...) {
  print.SHC_base_gapstat(x, ...)
}


#' @export
plot.orSHC_SSL_gapstat <- function(x, ...) {
  plot.SHC_base_gapstat(x, ...)
}


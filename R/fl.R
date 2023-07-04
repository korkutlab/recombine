getuw_fl <- function(ds,
                     lambda1s = NULL,
                     lambda2s = NULL,
                     nlambda1 = 10,
                     nlambda2 = 10,
                     max.iter = 500,
                     iseed = NULL,
                     silent = FALSE) {
  # dim
  n2 <- nrow(ds)
  p <- ncol(ds)
  n <- ceiling(sqrt(2 * n2)) # n2 = n(n-1)/2
  
  # lambda1s
  if (is.null(lambda1s)) {
    lambda1s <- seq(1, n, length = nlambda1)
  } else {
    nlambda1 <- length(lambda1s)
    
    if (nlambda1 > 1) {
      # Lambda1s should be an increasing sequence
      monotone <- sum((lambda1s[-1] - lambda1s[-nlambda1]) > 0)
      if (monotone != nlambda1 - 1){
        stop("lambda1s must be a monotone increasing sequence")
      }
    }
  }
  
  # lambda2s
  if (is.null(lambda2s)) {
    lambda2s <- seq(1, n, length = nlambda2)
  } else {
    nlambda2 <- length(lambda2s)
    
    # Lambda2s should be an increasing sequence
    monotone <- sum((lambda2s[-1] - lambda2s[-nlambda2]) > 0)
    if (monotone != nlambda2 - 1){
      stop("lambda2s must be a monotone increasing sequence")
    }
  }
  
  L1 <- length(lambda1s)
  L2 <- length(lambda2s)
  L <- L1*L2
  
  lambda1s_path <- rep(lambda1s, each = L2)
  lambda2s_path <- rep(lambda2s, times = L1)
  
  if (is.null(iseed)) {
    init_random <- FALSE
    iseed <- 42L
  } else {
    init_random <- TRUE
  }
  
  # call
  set.seed(iseed)
  res <- SHC_FL_getuw(ds,
                      as.numeric(lambda1s),
                      as.numeric(lambda2s),
                      as.integer(max.iter),
                      init_random,
                      silent)
  
  u <- res[["u"]]
  w <- res[["w"]]
  crit <- res[["crit"]]
  iter <- res[["iter"]]
  flsa_converged <- res[["flsa_converged"]]
  
  # # check convergence
  # if (any(iter == max.iter)) {
  #   message("Algorithm did not converge for the following regularization values")
  #   print('Lambda1')
  #   print(lambda1s_path[iter == max.iter])
  #   print('Lambda2')
  #   print(lambda2s_path[iter == max.iter])
  # }
  #
  # if (iter[L] == max.iter) {
  #   message("Algorithm did not converge at the last regularization value")
  # }
  
  # postprocess u
  u_list <- list()
  for (i in c(1:L)) {
    u_list[[i]] <- post_u(u[((i-1)*n2+1):(i*n2)])
  }
  
  # postprocess w
  w_list <- list()
  for (i in c(1:L)) {
    w_list[[i]] <- w[((i-1)*p+1):(i*p)]
  }
  
  w_l0norm <- sapply(w_list, l0norm)
  nonzero_w_indices_list <- lapply(w_list, function(x) which(x > 0))
  
  # sucessive differences of w
  w_diff_list <- list()
  for (i in c(1:L)) {
    w_diff_list[[i]] <- w_list[[i]][2:p] - w_list[[i]][1:(p-1)]
  }
  
  w_diff_l0norm <- sapply(w_diff_list, l0norm)
  
  return(list(lassotype = 'FL',
              lassotype_full = 'Fused LASSO',
              u = u_list, w = w_list,
              crit = crit, iter = iter,
              flsa_converged = flsa_converged,
              w_l0norm = w_l0norm,
              w_diff_l0norm = w_diff_l0norm,
              nonzero_w_indices = nonzero_w_indices_list,
              hyp1s = lambda1s_path,
              hyp2s = lambda2s_path))
}


#' orSHC_FL
#'
#' @description
#' Outlier-robust sparse hierarchical clustering with fused lasso
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param lambda1 sparsity parameter. 
#' @param lambda2 successive difference sparsity parameter. 
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
#' This function, orSHC_FL, which robustifies SHC_FL, is robust to the existence of outliers. 
#' It trims outliers during feature selection and perform final clustering using all samples. 
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{u}{the optimal weighted dissimilarity matrix.}
#' \item{w}{the optimal weights.}
#' \item{crit}{the optimal objective.}
#' \item{iter}{number of iteration steps.}
#' \item{flsa_converged}{indicator of convergence of the FLSA (Fused Lasso Signal Approximator) solver.}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{w_diff_l0norm}{L0 norm of successive differences of w's.}
#' \item{nonzero_w_indices}{indices of nonzero w's.}
#' \item{hyp1}{hyperparameter lambda1 used.}
#' \item{hyp2}{hyperparameter lambda2 used.}
#' \item{is_outlier}{outlier indicator of each sample.}
#' \item{loop}{LoOP value of each sample.}
#' \item{hc}{an object of class hclust by running hierarchical clustering on u.}
#' 
#' @examples
#' source(system.file("scripts", "gen_sim_data.R", package = "recombine"))
#' d <- gen_sim_data(out_pct = 0.1, iseed = 1)
#' 
#' # reorder features based on generative process
#' order <- c()
#' for (i in 1:n_clusters) {
#'   order <- c(order, seq(from = i, to = p_inf, by = n_clusters))
#' }
#' order <- c(order, (p_inf+1):p)
#' 
#' mt <- d$x[, order]
#' 
#' # run
#' result <- orSHC_FL(d$x,
#'                    lambda1 = 200,
#'                    lambda2 = 3000)
#' result
#' 
#' @export
#'
orSHC_FL <- function(
  x,
  lambda1,
  lambda2,
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
  max.iter = 100,
  iseed = NULL,
  silent = FALSE)
{
  lassotype <- "fl"
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
  out <- getuw_fl(ds = ds,
                  lambda1s = lambda1,
                  lambda2s = lambda2,
                  max.iter = max.iter,
                  iseed = iseed,
                  silent = silent)
  
  # step 3: get u for all data using outlier-free w
  # wrap out
  class_name <- "orSHC_FL"
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


#' SHC_FL
#'
#' @description
#' Sparse hierarchical clustering with fused lasso
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param lambda1 sparsity parameter. 
#' @param lambda2 successive difference sparsity parameter. 
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#'
#' @details
#' SHC_FL employs a sparse fused lasso penalty to perform feature selection coupled with clustering. 
#' It encourages sparsity of features, and encourages neighboring features to be similar and some to be identical.
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{u}{the optimal weighted dissimilarity matrix.}
#' \item{w}{the optimal weights.}
#' \item{crit}{the optimal objective.}
#' \item{iter}{number of iteration steps.}
#' \item{flsa_converged}{indicator of convergence of the FLSA (Fused Lasso Signal Approximator) solver.}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{w_diff_l0norm}{L0 norm of successive differences of w's.}
#' \item{nonzero_w_indices}{indices of nonzero w's.}
#' \item{hyp1}{hyperparameter lambda1 used.}
#' \item{hyp2}{hyperparameter lambda2 used.}
#' \item{hc}{an object of class hclust by running hierarchical clustering on u.}
#' 
#' @examples
#' d <- example_sim_data
#' 
#' # reorder features based on generative process
#' order <- c()
#' for (i in 1:n_clusters) {
#'   order <- c(order, seq(from = i, to = p_inf, by = n_clusters))
#' }
#' order <- c(order, (p_inf+1):p)
#' 
#' mt <- d$x[, order]
#' 
#' # run
#' result <- SHC_FL(d$x,
#'                  lambda1 = 200,
#'                  lambda2 = 3000)
#' result
#' 
#' @export
#'
SHC_FL <- function(
  x,
  lambda1,
  lambda2,
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
  max.iter = 100,
  iseed = NULL,
  silent = FALSE)
{
  result <- orSHC_FL(
    x = x,
    lambda1 = lambda1,
    lambda2 = lambda2,
    standardize.arrays = standardize.arrays,
    dissimilarity = dissimilarity,
    method = method,
    loop_threshold = NULL,
    max.iter = max.iter,
    iseed = iseed,
    silent = silent)
  
  result$is_outlier <- NULL
  result$loop <- NULL
  result$class_name <- "SHC_FL"
  
  return(result)
}


#' orSHC_FL_gapstat
#'
#' @description
#' This function computes gap statistic profile for a series of hyperparameters (lambda1's and lambda2's) by running orSHC_FL on both original data and permuted data.
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param lambda1s sparsity parameters. If Null, lambda1s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 10)}. The default is NULL.
#' @param lambda2s successive difference sparsity parameters. If Null, lambda2s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 10)}. The default is NULL.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param loop_k size of k-nearest neighbors defining local contexts in local outlier probability (LoOP) calculation. The default is 20.
#' @param loop_lambda scaling parameter in local outlier probability (LoOP) calculation. The default is 3.
#' @param loop_threshold threshold of LoOP above which LoOP values imply outliers. The default is 0.5.
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#' @param return_all_results should orSHC_FL results for all hyperparameters be returned? If False, only the result for the best hyperparameter is returned. The default is FALSE.
#' @param nperms number of permutations. The default is 10.
#' @param sel_rule rule for choosing the best hyperparameter. This should be either "max" (maximum gap statistic) or "1se.rule" (having a low complexity but within 1 standard error of the maximum gap statistic). The default is "max".
#'
#' @details
#' For a hyperparameter, the gap statistic measures the strength of the clustering based on real data with respect to the one based on randomly permuted data that are supposed to have no cluster. 
#' This function runs orSHC_FL for every hyperparameter (i.e., a pair of lambda1 and lambda2), calculates gap statistic, and get the best hyperparamter that corresponds to the maximum gap statistic.
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{tots}{the strength of the clustering using real data (i.e., maximized objective function values).}
#' \item{permtots}{the strength of the clustering using permuted data (i.e., maximized objective function values).}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{gaps_mean}{mean values of gap statistic.}
#' \item{gaps_se}{standard errors of gap statistic.}
#' \item{hyp1s}{hyperparameters lambda1's.}
#' \item{hyp2s}{hyperparameters lambda2's.}
#' \item{best_hyp1}{best hyperparameter lambda1.}
#' \item{best_hyp2}{best hyperparameter lambda2.}
#' \item{best_index}{index of best hyperparameter (i.e., the optimal pair of lambda1 and lambda2). The pairs are stored in a 1-dimensional vector, with a lambda1-major order.}
#' \item{result}{orSHC_FL result. If return_all_results is TRUE, results of all hyperparameters are stored as a list, in which each element is a result of orSHC_FL; otherwise, it only contains the result corresponding to the best hyperparameter.}
#' 
#' @examples
#' source(system.file("scripts", "gen_sim_data.R", package = "recombine"))
#' d <- gen_sim_data(out_pct = 0.1, iseed = 1)
#' 
#' # reorder features based on generative process
#' order <- c()
#' for (i in 1:n_clusters) {
#'   order <- c(order, seq(from = i, to = p_inf, by = n_clusters))
#' }
#' order <- c(order, (p_inf+1):p)
#' 
#' mt <- d$x[, order]
#' 
#' # to save time, run a 5-by-5 grid search; 
#' # in practice, a more refined grid search may be performed.
#' out <- orSHC_FL_gapstat(d$x,
#'                         lambda1s = seq(0, 0.5*max(nrow(mt), ncol(mt)), length = 5),
#'                         lambda2s = seq(1, 10*max(nrow(mt), ncol(mt)), length = 5),
#'                         nperms = 3)
#' out
#' 
#' @export
#'
orSHC_FL_gapstat <- function(
  x,
  lambda1s = NULL,
  lambda2s = NULL,
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
  max.iter = 100,
  iseed = NULL,
  silent = FALSE,
  return_all_results = FALSE,
  nperms = 10,
  sel_rule = c("max", "1se.rule"))
{
  lassotype <- "fl"
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
  
  if(is.null(lambda1s)) lambda1s <- seq(1, max(nrow(x), ncol(x)), length = 10)
  if(is.null(lambda2s)) lambda2s <- seq(1, max(nrow(x), ncol(x)), length = 10)
  if(length(lambda2s) < 2 & length(lambda2s) < 2) stop("Lambda1s or Lambda2s should be a vector with at least 2 elements.")
  n_regul <- length(lambda1s)*length(lambda2s)
  
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
  out <- getuw_fl(ds = ds,
                  lambda1s = lambda1s,
                  lambda2s = lambda2s,
                  max.iter = max.iter,
                  iseed = iseed,
                  silent = silent)
  
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
    perm.out <- getuw_fl(ds = permds,
                         lambda1s = lambda1s,
                         lambda2s = lambda2s,
                         max.iter = max.iter,
                         iseed = iseed,
                         silent = silent)
    permtots[k,] <- perm.out$crit
  }
  
  gaps_mean <- log(tots + 1) - colMeans(log(permtots + 1))
  gaps_se <- apply(log(permtots + 1), 2, sd) / sqrt(nperms)
  
  if (sel_rule == "max") {
    sel <- get_2d_best_index_at_max(gaps_mean,
                                    lambda1s,
                                    lambda2s)
  } else {
    sel <- get_2d_best_index_1se_rule_topright(gaps_mean,
                                               gaps_se,
                                               lambda1s,
                                               lambda2s)
  }
  
  best_index <- sel$best_index
  best_hyp1 <- lambda1s[sel$best_index1]
  best_hyp2 <- lambda2s[sel$best_index2]
  
  # get all results or best result using best hyp
  class_name <- "orSHC_FL"
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
              hyp1s=lambda1s,
              hyp2s=lambda2s,
              best_hyp1=best_hyp1,
              best_hyp2=best_hyp2,
              best_index=best_index,
              result=result)
  
  class(out) <- "orSHC_FL_gapstat"
  return(out)
}


#' SHC_FL_gapstat
#'
#' @description
#' This function computes gap statistic profile for a series of hyperparameters (lambda1's and lambda2's) by running SHC_FL on both original data and permuted data.
#'
#' @param x data matrix. Rows are samples, while columns are features.
#' @param lambda1s sparsity parameters. If Null, lambda1s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 10)}. The default is NULL.
#' @param lambda2s successive difference sparsity parameters. If Null, lambda2s is assigned as \code{seq(1, max(nrow(x), ncol(x)), length = 10)}. The default is NULL.
#' @param standardize.arrays should the data matrix be standardized? The default is FALSE.
#' @param dissimilarity type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".
#' @param method the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".
#' @param max.iter maximum number of iterations. The default is 500.
#' @param iseed seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.
#' @param silent should progress messages be suppressed? The default is FALSE.
#' @param return_all_results should SHC_FL results for all hyperparameters be returned? If False, only the result for the best hyperparameter is returned. The default is FALSE.
#' @param nperms number of permutations. The default is 10.
#' @param sel_rule rule for choosing the best hyperparameter. This should be either "max" (maximum gap statistic) or "1se.rule" (having a low complexity but within 1 standard error of the maximum gap statistic). The default is "max".
#'
#' @details
#' For a hyperparameter, the gap statistic measures the strength of the clustering based on real data with respect to the one based on randomly permuted data that are supposed to have no cluster. 
#' This function runs SHC_FL for every hyperparameter (i.e., a pair of lambda1 and lambda2), calculates gap statistic, and get the best hyperparamter that corresponds to the maximum gap statistic.
#'
#' @return
#' \item{lassotype}{abbreviated name of lasso type.}
#' \item{lassotype_full}{full name of lasso type.}
#' \item{tots}{the strength of the clustering using real data (i.e., maximized objective function values).}
#' \item{permtots}{the strength of the clustering using permuted data (i.e., maximized objective function values).}
#' \item{w_l0norm}{L0 norm of w.}
#' \item{gaps_mean}{mean values of gap statistic.}
#' \item{gaps_se}{standard errors of gap statistic.}
#' \item{hyp1s}{hyperparameters lambda1's.}
#' \item{hyp2s}{hyperparameters lambda2's.}
#' \item{best_hyp1}{best hyperparameter lambda1.}
#' \item{best_hyp2}{best hyperparameter lambda2.}
#' \item{best_index}{index of best hyperparameter (i.e., the optimal pair of lambda1 and lambda2). The pairs are stored in a 1-dimensional vector, with a lambda1-major order.}
#' \item{result}{SHC_FL result. If return_all_results is TRUE, results of all hyperparameters are stored as a list, in which each element is a result of SHC_FL; otherwise, it only contains the result corresponding to the best hyperparameter.}
#' 
#' @examples
#' d <- example_sim_data
#' 
#' # reorder features based on generative process
#' order <- c()
#' for (i in 1:n_clusters) {
#'   order <- c(order, seq(from = i, to = p_inf, by = n_clusters))
#' }
#' order <- c(order, (p_inf+1):p)
#' 
#' mt <- d$x[, order]
#' 
#' # to save time, run a 5-by-5 grid search; 
#' # in practice, a more refined grid search may be performed.
#' out <- SHC_FL_gapstat(d$x,
#'                       lambda1s = seq(0, 0.5*max(nrow(mt), ncol(mt)), length = 5),
#'                       lambda2s = seq(1, 10*max(nrow(mt), ncol(mt)), length = 5),
#'                       nperms = 3)
#' out
#' 
#' @export
#'
SHC_FL_gapstat <- function(
  x,
  lambda1s = NULL,
  lambda2s = NULL,
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
  max.iter = 100,
  iseed = NULL,
  silent = FALSE,
  return_all_results = FALSE,
  nperms = 10,
  sel_rule = c("max", "1se.rule"))
{
  out <- orSHC_FL_gapstat(
    x = x,
    lambda1s = lambda1s,
    lambda2s = lambda2s,
    standardize.arrays = standardize.arrays,
    dissimilarity = dissimilarity,
    method = method,
    loop_threshold = NULL,
    max.iter = max.iter,
    iseed = iseed,
    silent = silent,
    return_all_results = return_all_results,
    nperms = nperms,
    sel_rule = sel_rule)
  
  if (return_all_results) {
    out$result <- lapply(out$result, 
                         function(result) {
                           result$is_outlier <- NULL
                           result$loop <- NULL
                           result$class_name <- "SHC_FL"
                           return(result)
                         }
    )
  } else {
    out$result$is_outlier <- NULL
    out$result$loop <- NULL
    out$result$class_name <- "SHC_FL"
  }
  
  class(out) <- "SHC_FL_gapstat"
  return(out)
}


#' @export
print.SHC_FL <- function(x, ...) {
  print.SHC_hyp2d(x, ...)
}


#' @export
print.SHC_FL_gapstat <- function(x, ...) {
  mat <- round(cbind(x$hyp1s, x$hyp2s, x$w_l0norm, x$gaps_mean, x$gaps_se), 5)
  dimnames(mat) <- list(1:length(x$gaps_mean), c("Hyperparameter1", "Hyperparameter2", "# Non-Zero W's", "Gap Statistic", "Standard Error"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", x$best_hyp1, " ", x$best_hyp2, fill=TRUE)
}


#' @export
print.orSHC_FL <- function(x, ...) {
  print.SHC_hyp2d(x, ...)
}


#' @export
print.orSHC_FL_gapstat <- function(x, ...) {
  mat <- round(cbind(x$hyp1s, x$hyp2s, x$w_l0norm, x$gaps_mean, x$gaps_se), 5)
  dimnames(mat) <- list(1:length(x$gaps_mean), c("Hyperparameter1", "Hyperparameter2", "# Non-Zero W's", "Gap Statistic", "Standard Error"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", x$best_hyp1, " ", x$best_hyp2, fill=TRUE)
}


l0norm <- function(x) {
  return(sum(x > 1e-10))
}

l1norm <- function(x) {
  return(sum(abs(x)))
}

l2norm <- function(x) {
  return(sqrt(sum(x*x)))
}

l2norm_squared <- function(x) {
  return(sum(x*x))
}


# this function accepts missing values
#
# in daisy from cluster:
# If the metric is "euclidean",
# and n_g is the number of columns in which neither row i and j have NAs,
# then the dissimilarity d(i,j) returned is
# sqrt(p/n_g) (p=ncol(x)) times the Euclidean distance between
# the two vectors of length n_g shortened to exclude NAs.
# The rule is similar for the "manhattan" metric,
# except that the coefficient is p/n_g.
#
x2dist <- function(x,
                   dissimilarity = c("squared.distance",
                                     "absolute.value",
                                     "euclidean")) {
  dissimilarity <- match.arg(dissimilarity)
  
  # calculate ds assuming dissimilarity = "absolute.value"
  xnona <- x
  xnona[is.na(x)] <- 0
  ds <- matrix(distfun(xnona), ncol=ncol(x))
  
  # adjust ds if dissimilarity is squared.distance or euclidean
  if (dissimilarity %in% c("squared.distance", "euclidean")) ds <- ds^2
  
  # adjust ds if NA exists
  if (sum(is.na(x)) > 0) {
    # get indicator matrix
    xbinary <- matrix(1, nrow = nrow(x), ncol = ncol(x))
    xbinary[is.na(x)] <- 0
    mult <- matrix(multfun(xbinary), ncol=ncol(x))
    
    # increase distances of non-missing features to account for artificial zero distances of missing features
    ds <- sweep(ds,
                1,
                ncol(ds)/apply(mult != 0, 1, sum),
                "*")
    
    # assign zero to distances of missing features
    ds[mult == 0] <- 0
  }
  
  # get dists by aggregate ds
  dists_flat <- rowSums(ds)
  
  # adjust dists if dissimilarity is euclidean
  if (dissimilarity == "euclidean") dists_flat <- sqrt(dists_flat)
  
  # transform shape
  n <- nrow(x)
  dists <- matrix(0,
                  nrow = n,
                  ncol = n)
  dists[lower.tri(dists)] <- dists_flat
  dists <- as.dist(dists)
  
  return(dists)
}


# this function is used for sparse clustering, which doesn't accept missing values
x2ds_nona <- function(x,
                      dissimilarity = c("squared.distance",
                                        "absolute.value")) {
  dissimilarity <- match.arg(dissimilarity)
  
  if(sum(is.na(x)) > 0) stop("x can't have missing values")
  ds <- matrix(distfun(x), ncol=ncol(x))
  
  if (dissimilarity == "squared.distance") ds <- ds^2
  
  return(ds)
}


# convert ds to distance matrix
ds2dist <- function(ds) {
  n2 <- nrow(ds)
  n <- as.integer(sqrt(2 * nrow(ds))) + 1
  
  # calculate distance
  dists_flat <- rowSums(ds)
  
  dists <- matrix(0,
                  nrow = n,
                  ncol = n)
  dists[lower.tri(dists)] <- dists_flat
  dists <- as.dist(dists)
  
  return(dists)
}


# convert ds and w to distance matrix
ds_w_2_dist <- function(ds, w) {
  n2 <- nrow(ds)
  p <- ncol(ds)
  n <- as.integer(sqrt(2 * nrow(ds))) + 1
  
  # calculate weighted distance
  ds_w <- ds * rep(w, each = n2)
  dists_flat <- rowSums(ds_w)
  
  dists <- matrix(0,
                  nrow = n,
                  ncol = n)
  dists[lower.tri(dists)] <- dists_flat
  dists <- as.dist(dists)
  
  return(dists)
}


post_u <- function(u) {
  # dim
  n2 <- length(u)
  n <- ceiling(sqrt(2 * n2)) # n2 = n(n-1)/2
  
  # format u as a matrix
  u2 <- matrix(0,
               nrow = n,
               ncol = n)
  u2[lower.tri(u2)] <- u
  u2 <- as.matrix(as.dist(u2)) / sqrt(2)
  
  return(u2)
}


u2hclust <- function(u,
                     labels = NULL,
                     method = c("average",
                                "complete",
                                "single",
                                "centroid",
                                "ward.D",
                                "ward.D2",
                                "mcquitty",
                                "median")) {
  method <- match.arg(method)
  hc <- hclust(as.dist(u),
               method = method)
  
  # the original hc has numbers as lables
  # change the labels to sampleid
  if (!is.null(labels)) {
    hc$labels <- labels
  }
  
  return(hc)
}


#' @importFrom stringr str_replace
#' @importFrom magrittr %>%
extract_unperm_out <- function(index,
                               unperm_outs,
                               ds = NULL,
                               is_outlier,
                               loop,
                               dissimilarity,
                               method,
                               labels,
                               class_name) {
  # extract elements from lists
  out <- list()
  vars_identical <- c("lassotype",
                      "lassotype_full")
  for (name in names(unperm_outs)) {
    if (name %in% vars_identical) {
      out[[name]] <- unperm_outs[[name]]
    } else {
      out[[name]] <- unperm_outs[[name]][[index]]
    }
  }
  class(out) <- class(unperm_outs)
  
  out$is_outlier <- is_outlier
  out$loop <- loop
  
  # recalculate u using w for given data
  if (!is.null(ds)) {
    res <- SHC_get_u(ds, out$w)
    
    # postprocess u
    out$u <- post_u(res[["u"]])
  }
  
  names(out) <- names(out) %>%
    str_replace("^hyps$", "hyp") %>%
    str_replace("^hyp1s$", "hyp1") %>%
    str_replace("^hyp2s$", "hyp2")
  
  
  hc <- u2hclust(out$u,
                 labels = labels,
                 method = method)
  hc$dist.method <- dissimilarity
  out[["hc"]] <- hc
  
  class(out) <- class_name
  return(out)
}


get_best_index_at_max <- function(gaps) {
  best_index <- which.max(gaps)
  
  # # the optimal gap needs to be positive, meaning more structure than null
  # if (gaps[best_index] > 0) {
  #   return(best_index)
  # }
  # else {
  #   return(NA)
  # }
  # for synthetic data, random shuffling may increase structures,
  # thus larger crit than the original synthetic data
  
  return(best_index)
}


get_best_index_1se_rule_left <- function(gaps, segaps) {
  gm <- max(gaps)
  sm <- segaps[which.max(gaps)]
  g <- gm - sm
  best_index <- which(gaps > g)[1]
  
  # # the optimal gap needs to be positive, meaning more structure than null
  # if (gaps[best_index] > 0) {
  #   return(best_index)
  # }
  # else {
  #   return(NA)
  # }
  # for synthetic data, random shuffling may increase structures,
  # thus larger crit than the original synthetic data
  
  return(best_index)
}


get_best_index_1se_rule_right <- function(gaps, segaps) {
  gm <- max(gaps)
  sm <- segaps[which.max(gaps)]
  g <- gm - sm
  best_gaps <- which(gaps > g)
  best_index <- best_gaps[length(best_gaps)]
  
  # # the optimal gap needs to be positive, meaning more structure than null
  # if (gaps[best_index] > 0) {
  #   return(best_index)
  # }
  # else {
  #   return(NA)
  # }
  # for synthetic data, random shuffling may increase structures,
  # thus larger crit than the original synthetic data
  
  return(best_index)
}


get_2d_best_index_at_max <- function(gaps, lambda1s, lambda2s) {
  # calculation is done in C
  # lambda1s are rows, lambda2s are columns
  l1 <- length(lambda1s)
  l2 <- length(lambda2s)
  
  # get id of max in flattened
  idmax <- which.max(gaps)
  
  # get position in 2d
  id1 <- (idmax - 1) %/% l2 + 1
  id2 <- (idmax - 1) %% l2 + 1
  
  # # the optimal gap needs to be positive, meaning more structure than null
  # if (gaps[idmax] > 0) {
  #   return(list(best_index = idmax,
  #               best_index1 = id1,
  #               best_index2 = id2))
  # }
  # else {
  #   return(list(best_index = NA,
  #               best_index1 = NA,
  #               best_index2 = NA))
  # }
  # for synthetic data, random shuffling may increase structures,
  # thus larger crit than the original synthetic data
  
  return(list(best_index = idmax,
              best_index1 = id1,
              best_index2 = id2))
}


get_2d_best_index_1se_rule_topright <- function(gaps, segaps, lambda1s, lambda2s) {
  # calculation is done in C
  # lambda1s are rows, lambda2s are columns
  l1 <- length(lambda1s)
  l2 <- length(lambda2s)
  
  # get id of max in flattened
  idmax <- which.max(gaps)
  g <- gaps[idmax] - segaps[idmax]
  
  # get position in 2d
  id1 <- (idmax - 1) %/% l2 + 1
  id2 <- (idmax - 1) %% l2 + 1
  
  # mask bottom-left and uninterested regions and get interested position
  mt_gaps <- matrix(gaps, nrow = l1, ncol = l2, byrow = TRUE) # Row major as in C style
  mt_gaps[1:(id1-1), 1:(id2-1)] <- -Inf
  mt_gaps[mt_gaps > g] <- -Inf
  
  best_index <- which.max(mt_gaps)
  
  # get position in 2d
  id1 <- (best_index - 1) %/% l2 + 1
  id2 <- (best_index - 1) %% l2 + 1
  
  # # the optimal gap needs to be positive, meaning more structure than null
  # if (gaps[best_index] > 0) {
  #   return(list(best_index = best_index,
  #               best_index1 = id1,
  #               best_index2 = id2))
  # }
  # else {
  #   return(list(best_index = NA,
  #               best_index1 = NA,
  #               best_index2 = NA))
  # }
  # for synthetic data, random shuffling may increase structures,
  # thus larger crit than the original synthetic data
  
  return(list(best_index = best_index,
              best_index1 = id1,
              best_index2 = id2))
}


se <- function(x, na.rm = FALSE) {
  if (na.rm) {
    n <- length(x[!is.na(x)])
  } else {
    n <- length(x)
  }
  sqrt(var(x, na.rm = na.rm)/n)
}


print.SHC_base <- function(x, ...) {
  cat("Lasso Type: ", x$lassotype_full, fill=TRUE)
  cat("Hyperparameters: ", x$hyp, fill=TRUE)
  cat(fill=TRUE)
}


print.SHC_hyp2d <- function(x, ...) {
  cat("Lasso Type: ", x$lassotype_full, fill=TRUE)
  cat("Hyperparameters: ", x$hyp1, "\t", x$hyp2, fill=TRUE)
  cat(fill=TRUE)
}


print.SHC_base_gapstat <- function(x, ...) {
  mat <- round(cbind(x$hyps, x$w_l0norm, x$gaps_mean, x$gaps_se), 4)
  dimnames(mat) <- list(1:length(x$gaps_mean), c("Hyperparameter", "# Non-Zero W's", "Gap Statistic", "Standard Error"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", x$best_hyp, fill=TRUE)
}


plot.SHC_base_gapstat <- function(x, ...) {
  plot(x$w_l0norm[x$w_l0norm > 0],
       x$gaps_mean[x$w_l0norm > 0],
       log="x",
       main="Gap Statistics",
       xlab="# Non-zero Wj's",
       ylab="")
  lines(x$w_l0norm[x$w_l0norm > 0],
        x$gaps_mean[x$w_l0norm > 0])
}


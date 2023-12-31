% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fl.R
\name{orSHC_FL}
\alias{orSHC_FL}
\title{orSHC_FL}
\usage{
orSHC_FL(
  x,
  lambda1,
  lambda2,
  standardize.arrays = FALSE,
  dissimilarity = c("squared.distance", "absolute.value"),
  method = c("average", "complete", "single", "centroid", "ward.D", "ward.D2",
    "mcquitty", "median"),
  loop_k = 20,
  loop_lambda = 3,
  loop_threshold = 0.5,
  max.iter = 100,
  iseed = NULL,
  silent = FALSE
)
}
\arguments{
\item{x}{data matrix. Rows are samples, while columns are features.}

\item{lambda1}{sparsity parameter.}

\item{lambda2}{successive difference sparsity parameter.}

\item{standardize.arrays}{should the data matrix be standardized? The default is FALSE.}

\item{dissimilarity}{type of dissimilarity metric. This should be either "squared.distance" or "absolute.value". The default is "squared.distance".}

\item{method}{the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC). The default is "average".}

\item{loop_k}{size of k-nearest neighbors defining local contexts in local outlier probability (LoOP) calculation. The default is 20.}

\item{loop_lambda}{scaling parameter in local outlier probability (LoOP) calculation. The default is 3.}

\item{loop_threshold}{threshold of LoOP above which LoOP values imply outliers. The default is 0.5.}

\item{max.iter}{maximum number of iterations. The default is 500.}

\item{iseed}{seed of random numbers used for initializing w. If Null, all elements of w are initialized as a fixed value of 1.0/sqrt(p). Otherwise, the initial w elements are uniformly distributed in \code{[0, 2.0/sqrt(p)]}. The default is NULL.}

\item{silent}{should progress messages be suppressed? The default is FALSE.}
}
\value{
\item{lassotype}{abbreviated name of lasso type.}
\item{lassotype_full}{full name of lasso type.}
\item{u}{the optimal weighted dissimilarity matrix.}
\item{w}{the optimal weights.}
\item{crit}{the optimal objective.}
\item{iter}{number of iteration steps.}
\item{flsa_converged}{indicator of convergence of the FLSA (Fused Lasso Signal Approximator) solver.}
\item{w_l0norm}{L0 norm of w.}
\item{w_diff_l0norm}{L0 norm of successive differences of w's.}
\item{nonzero_w_indices}{indices of nonzero w's.}
\item{hyp1}{hyperparameter lambda1 used.}
\item{hyp2}{hyperparameter lambda2 used.}
\item{is_outlier}{outlier indicator of each sample.}
\item{loop}{LoOP value of each sample.}
\item{hc}{an object of class hclust by running hierarchical clustering on u.}
}
\description{
Outlier-robust sparse hierarchical clustering with fused lasso
}
\details{
This function, orSHC_FL, which robustifies SHC_FL, is robust to the existence of outliers.
It trims outliers during feature selection and perform final clustering using all samples.
}
\examples{
source(system.file("scripts", "gen_sim_data.R", package = "recombine"))
d <- gen_sim_data(out_pct = 0.1, iseed = 1)

# reorder features based on generative process
order <- c()
for (i in 1:n_clusters) {
  order <- c(order, seq(from = i, to = p_inf, by = n_clusters))
}
order <- c(order, (p_inf+1):p)

mt <- d$x[, order]

# run
result <- orSHC_FL(d$x,
                   lambda1 = 200,
                   lambda2 = 3000)
result

}

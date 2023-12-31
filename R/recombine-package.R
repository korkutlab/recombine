#' recombine package
#'
#' RECOMBINE: Recurrent Composite Markers for Biological Identities with Neighborhood Enrichment
#'
#' This package includes five novel algorithms for feature selection coupled with clustering. Briefly, we advance SHC (sparse hierarchical clustering) and obtain two algorithms: SHC-SSL (sparse hierarchical clustering with spike-and-slab lasso), which replaces lasso with spike-and-slab lasso to debias the effect of lasso without a priori information, and SHC-FL (sparse hierarchical clustering with fused lasso), which replaces lasso with fused lasso to incorporate prior information of feature order. Further, rSHC, rSHC-SSL, and rSHC-FL, are robust versions of SHC, SHC-SSL, and SHC-FL, respectively, in which selected features are robust to the existence of outliers. 
#'
#' @name recombine-package
#' @docType package
#' @references
#' Li X. and Korkut A. (2023) Recurrent composite markers for cell subpopulations with high granularity.
#' @keywords package
#' @useDynLib recombine
#' @importFrom graphics lines plot
#' @importFrom stats as.dist hclust rnorm runif sd var
#' @importFrom Rcpp evalCpp
NULL


#' An example dataset of simulation data
#'
#' This dataset has 3 clusters, each of which has 40 data points with 50 informative features and 950 noise features
#'
#' This dataset can be generated by running: 
#' \code{source(system.file("scripts", "gen_sim_data.R", package = "recombine"))} 
#' followed by
#' \code{example_sim_data <- gen_sim_data(size_clusters = c(40, 40, 40),
#'                                        p_inf = 50,
#'                                        p_noise = 950,
#'                                        iseed = 1)}
#' 
#' @name example_sim_data
#' @docType data
#' @format The format is a list containing the following elements:
#' - x: data matrix of the simulation dataset.
#' - y: integer vector corresponding to the cluster membership.
#'
#' @keywords datasets
#' @examples
#' library(recombine)
#' example_sim_data
#'
NULL


#' scRNA data and RCMs of mouse intestinal organoids
#'
#' The raw single cell expression dataset was obtained from Grün, D., et al. (2015). Single-cell messenger RNA sequencing reveals rare intestinal cell types. Nature, 525(7568), 251–255.
#' RECOMBINE was applied to this dataset to extract recurrent composite markers (RCMs) for both major and rare cell subpopulations.
#' 
#' @name intestine_data
#' @docType data
#' @format The format is a list containing the following elements:
#' - mt_expr: single cell expression data matrix. mt_expr has been normalized and log-transformed. 
#' - w_nonzero: w values of discriminant markers with nonzero w's.
#' - df_cell_subpop: data frame of cell subpopulations.
#' - mt_nn: indicator matrix of nearest neighbors. mt_nn is obtained from the nearest neighbors calculated by Seurat based on discriminant markers.
#' - df_rcm_cell: data frame of RCMs per cell.
#' - df_rcm_subpop: data frame of RCMs per cell subpopulation.
#'
#' @keywords datasets
#' @examples
#' library(recombine)
#' intestine_data
#'
NULL


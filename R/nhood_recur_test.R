#' get_mt_nn
#'
#' @description
#' Get indicator matrix of nearest neighbors
#'
#' @param mt_dist distance matrix of cells. mt_dist needs to be a symmetric matrix or an object of stats::dist. 
#' @param k number of nearest neighbors. The default is 20.
#'
#' @details
#' Given a distance matrix of cells, get_mt_nn extracts K nearest neighbors of each cell to build an indicator matrix where each row represents a cell's neighborhood status, with 1 indicating a member of nearest neighbors.
#'
#' @return
#' \item{indicator matrix of nearest neighbors}{a binary matrix of ncells by ncells, with element \code{mt_nn[i,j]=1} indicating cell j is one of the nearest neighbors of cell i.}
#' 
#' @examples
#' x <- matrix(rnorm(1000), nrow = 100)
#' mt_nn <- get_mt_nn(dist(x))
#' mt_nn[1:10, 1:10]
#' 
#' @export
#'
get_mt_nn <- function(
  mt_dist, 
  k = 20
) 
{
  if (class(mt_dist) == "dist") {
    mt_dist <- as.matrix(mt_dist)
  } else if (!isSymmetric.matrix(mt_dist)) {
    stop("mt_dist needs to be a symmetric matrix or an object of stats::dist.")
  }
  
  mt_nn <- matrix(0, nrow = nrow(mt_dist), ncol = ncol(mt_dist))
  for (i in 1:nrow(mt_dist)) {
    # find k nearest neighbors
    # matches <- setdiff(order(mt_dist[i,],decreasing = F)[1:(k+1)],i) # exclude self
    matches <- order(mt_dist[i,],decreasing = F)[1:k] # include self
    mt_nn[i, matches] <- 1
  }
  dimnames(mt_nn) <- dimnames(mt_dist)
  return(mt_nn)
}


init_mt <- function(cells, markers) {
  mt <- matrix(NA, nrow = length(cells), ncol = length(markers))
  rownames(mt) <- cells
  colnames(mt) <- markers
  return(mt)
}


#' nhood_recur_test
#'
#' @description
#' Neighborhood recurrence test
#'
#' @param mt_expr data matrix. Rows are cells and columns are features.
#' @param mt_nn indicator matrix of nearest neighbors. mt_nn is a binary matrix of ncells by ncells, with \code{mt_nn[i,j]=1} indicating cell j is one of the nearest neighbors of cell i.
#' @param verbose logical variable of progress messages to be printed. The default is FALSE.
#'
#' @details
#' Neighborhood recurrence test determines if markers are recurrently up/down-regulated in the nearest neighborhood of each cell.
#'
#' @return
#' \item{data frame of RCMs per cell}{including columns: cell, marker, nhood_zscore, nhood_recur_pval, nhood_recur_fdr.}
#' 
#' @examples
#' df_rcm_cell <- nhood_recur_test(intestine_data$mt_expr[, names(intestine_data$w_nonzero)], intestine_data$mt_nn)
#' df_rcm_cell
#' 
#' @export
#'
nhood_recur_test <- function(
  mt_expr,
  mt_nn,
  verbose = FALSE
) 
{
  stopifnot("mt_expr need to have both valid rownames (cells) and colnames (markers)" = 
              (!is.null(rownames(mt_expr))) & (!is.null(colnames(mt_expr))))
  
  # each row of mt_nn corresponds to a cell with k nearest neighbors
  # rowSums(mt_nn) %>% summary()
  # colSums(mt_nn) %>% summary()
  
  # check consistency of k
  k_cells <- sum(mt_nn[1, ] > 0)
  stopifnot("mt_nn need to have same # of nn for each row" = 
              all(rowSums(mt_nn > 0) == k_cells))
  
  # check if cell names are consistent
  stopifnot("mt_nn need to have same dimensions of row and column" = 
              all.equal(rownames(mt_nn), colnames(mt_nn)))
  stopifnot("mt_expr and mt_nn need to have same dimension of rows" = 
              all.equal(rownames(mt_expr), rownames(mt_nn)))
  
  # get z-score and p-value
  cells <- rownames(mt_expr)
  markers <- colnames(mt_expr)
  
  mt_z_score <- init_mt(cells, markers)
  mt_p_value <- init_mt(cells, markers)
  
  for (marker in markers) {
    if (verbose) {
      message("processing ", marker)
    }
    
    expr <- mt_expr[, marker]
    
    # null model (CLT): get mean and sd of neighborhood mean expression
    mean_score <- mean(expr)
    sd_score <- sd(expr)/sqrt(k_cells)
    
    # calculate z-score and p-value for each cell
    for (cell in cells) {
      nn_cells <- colnames(mt_nn)[mt_nn[cell, ] > 0]
      nn_mean <- expr[nn_cells] %>% mean()
      z_score <- (nn_mean - mean_score)/sd_score
      p_value <- pnorm(abs(z_score), lower.tail = FALSE)
      mt_z_score[cell, marker] <- z_score
      mt_p_value[cell, marker] <- p_value
    }
  }
  
  df_zscore <- dplyr::tibble(cell = rownames(mt_z_score)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mt_z_score)) %>%
    tidyr::gather("marker", "nhood_zscore", -cell)
  
  df_pval <- dplyr::tibble(cell = rownames(mt_p_value)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mt_p_value))  %>%
    tidyr::gather("marker", "nhood_recur_pval", -cell)
  
  df <- df_zscore %>%
    dplyr::left_join(df_pval, by = c("cell", "marker")) %>%
    dplyr::mutate(nhood_recur_fdr = p.adjust(nhood_recur_pval, method = "fdr"))
  
  df <- df %>%
    dplyr::arrange(cell, nhood_recur_fdr)
  
  return(df) # RCMs per cell
}


#' get_rcm_subpop
#'
#' @description
#' Recurrent composite markers for cell subpopulations
#'
#' @param df_rcm_cell data frame of RCMs per cell, including columns: cell, marker, nhood_zscore, nhood_recur_pval, nhood_recur_fdr.
#' @param df_cell_subpop data frame of cell subpopulations, including columns: cell, name of subpopulations (e.g., cluster).
#' @param subpop_name name of subpopulations. The default is cluster.
#' @param nhood_recur_fdr_threshod threshold for the statistical significance of neighborhood recurrence FDRs. The default is 0.05.
#'
#' @details
#' Given the recurrent composite markers (RCMs) per cell and the membership of cell subpopulations, get_rcm_subpop calculates the average neighborhood Z score and fraction of significant cells for each marker for each cell subpopulation.
#'
#' @return
#' \item{data frame of RCMs per cell subpopulation}{including columns: subpop, marker, fract_signif_cells, avg_nhood_zscore, PR_AUC, ROC_AUC.}
#' 
#' @examples
#' df_rcm_subpop <- get_rcm_subpop(intestine_data$df_rcm_cell, intestine_data$df_cell_subpop)
#' df_rcm_subpop
#' 
#' @export
#'
get_rcm_subpop <- function(
  df_rcm_cell, 
  df_cell_subpop, 
  subpop_name = "cluster",
  nhood_recur_fdr_threshod = 0.05
) 
{
  stopifnot("colnames(df_cell_subpop) need to have cell and name of subpopulations (e.g., cluster)" =
              sum(!(c("cell", subpop_name) %in% colnames(df_cell_subpop))) == 0)
  stopifnot("colnames(df_rcm) need to have cell, marker, nhood_zscore and nhood_recur_fdr" = 
              sum(!(c("cell", "marker", "nhood_zscore", "nhood_recur_fdr") %in% colnames(df_rcm_cell))) == 0)
  
  stopifnot("df_rcm_cell has duplicated items of cell-marker pairs" =
              (nrow(dplyr::distinct(df_rcm_cell, cell, marker))) == (nrow(df_rcm_cell)))
  stopifnot("df_cell_subpop has duplicated cells" =
              (nrow(dplyr::distinct(df_cell_subpop, cell))) == (nrow(df_cell_subpop)))
  
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  
  # rename subpop for programming convenience
  colnames(df_cell_subpop)[colnames(df_cell_subpop) == subpop_name] <- "subpop"
  
  # merge rcm with subpop info and filter cells without subpop info
  df_data <- df_rcm_cell %>%
    dplyr::left_join(df_cell_subpop, by = "cell") %>%
    dplyr::filter(!is.na(subpop))
  
  # fraction of significant cells
  df_n <- df_cell_subpop %>% 
    dplyr::distinct(cell, subpop) %>%
    dplyr::group_by(subpop) %>% 
    dplyr::summarise(n_tot = dplyr::n()) %>% 
    dplyr::ungroup()
  
  df_fract <- df_data %>% 
    dplyr::filter(nhood_recur_fdr < nhood_recur_fdr_threshod) %>% 
    dplyr::group_by(subpop, marker) %>% 
    dplyr::summarise(n = dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(df_n, by = "subpop") %>% 
    dplyr::mutate(fract_signif_cells = n/n_tot) %>% 
    dplyr::select(subpop, marker, fract_signif_cells) %>% 
    dplyr::arrange(subpop, dplyr::desc(fract_signif_cells)) 
  
  # nhood_zscore mean for each marker
  df_zscore <- df_data %>%
    dplyr::group_by(subpop, marker) %>%
    dplyr::summarise(avg_nhood_zscore = mean(nhood_zscore)) %>%
    dplyr::ungroup()
  
  # AUC
  df_auc <- dplyr::tibble()
  markers <- df_data$marker %>% unique()
  for (m in markers) {
    df <- df_data %>%
      dplyr::filter(marker == m)
    z_scores <- df$nhood_zscore
    
    subpops <- df$subpop %>% unique()
    pr_aucs <- c()
    roc_aucs <- c()
    for (s in subpops) {
      class_labels <- ifelse(df$subpop == s, 1, -1)
      mm <- precrec::evalmod(scores = z_scores, labels = class_labels)
      aucs <- precrec::auc(mm) %>% 
        dplyr::select(curvetypes, aucs) %>% 
        tibble::deframe()
      pr_aucs <- c(pr_aucs, aucs["PRC"])
      roc_aucs <- c(roc_aucs, aucs["ROC"])
    }
    
    df_auc <- df_auc %>%
      dplyr::bind_rows(dplyr::tibble(subpop = subpops,
                                     marker = m,
                                     PR_AUC = pr_aucs,
                                     ROC_AUC = roc_aucs))
  }
  
  # merge
  df <- df_fract %>% 
    dplyr::left_join(df_zscore, by = c("subpop", "marker")) %>%
    dplyr::left_join(df_auc, by = c("subpop", "marker")) %>%
    dplyr::arrange(subpop, 
                   dplyr::desc(fract_signif_cells*sign(avg_nhood_zscore)),
                   dplyr::desc(avg_nhood_zscore)) 
  
  # rename subpop_name back
  colnames(df)[colnames(df) == "subpop"] <- subpop_name
  
  return(df)
}


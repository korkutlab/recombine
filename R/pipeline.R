# prep data and pseudo-cells if needed
data_prep <- function(sobj,
                      subpop_name = NULL,
                      max_features = 10000,
                      pseudocell_embedding_nhvgs = 2000,
                      pseudocell_embedding_npcs = 50,
                      pseudocell_embedding_batchname = NULL,
                      pseudocell_embedding_max.iter.harmony = 100,
                      pseudocell_microcluster_res = 50,
                      pseudocell_microcluster_res.iter.rate = 0.8) {
  # check if LogNorm data is provided
  stopifnot(all.equal(dim(sobj@assays$RNA$data),
                      dim(sobj)))
  
  # identify max features by hvg as potential features
  if (nrow(sobj) > max_features) {
    message("Selecting ", max_features, " HVGs as the gene pool for RECOMBINE marker selection")
    sobj <- Seurat::FindVariableFeatures(sobj, selection.method = "vst", nfeatures = max_features)
    features <- Seurat::VariableFeatures(sobj)
    sobj <- sobj[features, ]
  }
  
  # max (pseudo-)cells
  max_nc <- sqrt(2**31/nrow(sobj)) %>% floor()
  
  # check if max_features is too high
  if (max_nc < 100) {
    stop("Too high max_features that leads to <100 (pseudo-)cells. Please reduce max_features below 200000.")
  }
  
  # expressions
  mt <- sobj@assays$RNA$data
  mt <- mt %>%
    as.matrix() %>%
    t()
  tb_data <- dplyr::tibble(cell = rownames(mt)) %>%
    dplyr::bind_cols(dplyr::as_tibble(mt))
  
  if (is.null(subpop_name)) {
    tb_anno <- NULL
  } else {
    tb_anno <- dplyr::tibble(cell = colnames(sobj),
                             subpop = sobj@meta.data[[subpop_name]])
    colnames(tb_anno)[2] <- subpop_name
  }
  
  # generate pseudo-cell data --------
  if (ncol(sobj) > max_nc) {
    message("Clustering with high resolution to find microclusters (pseudo-cells)")
    
    # identify top 2000 hvg
    sobj <- Seurat::FindVariableFeatures(sobj, 
                                         selection.method = "vst", 
                                         nfeatures = pseudocell_embedding_nhvgs)
    
    # Scaling the data
    sobj <- Seurat::ScaleData(sobj)
    
    # Perform linear dimensional reduction
    sobj <- Seurat::RunPCA(sobj, 
                           features = Seurat::VariableFeatures(object = sobj), 
                           npcs = pseudocell_embedding_npcs,
                           verbose = FALSE)
    
    # harmonoy
    if (!is.null(pseudocell_embedding_batchname)) {
      stopifnot(pseudocell_embedding_batchname %in% colnames(sobj@meta.data))
      sobj <- sobj %>%
        harmony::RunHarmony(pseudocell_embedding_batchname,
                            dims.use = 1:pseudocell_embedding_npcs,
                            max.iter.harmony = pseudocell_embedding_max.iter.harmony)
      sobj <- Seurat::FindNeighbors(sobj,
                                    reduction = "harmony", 
                                    dims = 1:pseudocell_embedding_npcs,
                                    verbose = FALSE)
    } else {
      sobj <- Seurat::FindNeighbors(sobj,
                                    dims = 1:pseudocell_embedding_npcs,
                                    verbose = FALSE)
    }
    
    # iterative clustering with adjusted res
    nc <- ncol(sobj)
    res <- pseudocell_microcluster_res
    while(nc > max_nc) {
      message("Running clustering with resolution ", res)
      sobj <- Seurat::FindClusters(sobj,
                                   resolution = res)
      
      # cluster IDs of cells
      cl <- Seurat::Idents(sobj) %>%
        as.character() %>%
        paste0("c", .)
      names(cl) <- names(Seurat::Idents(sobj))
      
      nc <- cl %>% unique() %>% length()
      res <- res*pseudocell_microcluster_res.iter.rate
    }
    
    # get expr for pseudo-cells using mean of contained cells
    mt <- as.matrix(tb_data[, -1])
    rownames(mt) <- tb_data$cell
    
    mt <- scale(mt, scale = FALSE)
    mt <- mt/sd(as.double(mt))
    
    tb_all <- dplyr::tibble(gene = colnames(mt))
    clusters <- unique(cl)
    for (cname in clusters) {
      cells <- names(cl)[cl == cname]
      tb <- dplyr::tibble(gene = colnames(mt)) %>%
        dplyr::bind_cols(dplyr::tibble(expr = colMeans(mt[cells, , drop = FALSE])))
      colnames(tb)[2] <- cname
      
      tb_all <- tb_all %>%
        dplyr::left_join(tb, by = "gene")
    }
    
    # transpose
    mt <- as.matrix(tb_all[, -1])
    rownames(mt) <- tb_all$gene
    
    mt <- t(mt)
    
    tb <- dplyr::tibble(cell = rownames(mt)) %>%
      dplyr::bind_cols(dplyr::as_tibble(mt))
    
    # sort pseudo-cells according to clusters
    tb <- tb %>%
      dplyr::mutate(cluster = stringr::str_replace(cell, "^c", "")) %>%
      # mutate(cluster = ifelse(grepl("singleton", cluster),
      #                         NA,
      #                         cluster)) %>%
      dplyr::mutate(cluster = as.integer(cluster)) %>%
      dplyr::arrange(cluster) %>%
      dplyr::select(-cluster)
    
    tb_pseudocell <- tb %>%
      dplyr::rename(pseudo_cell = cell)
    
    tb_pseudocell_member <- tb_pseudocell %>%
      dplyr::select(pseudo_cell) %>%
      dplyr::left_join(dplyr::tibble(cell = names(cl),
                                     pseudo_cell = cl), by = "pseudo_cell")
  } else {
    tb_pseudocell <- NULL
    tb_pseudocell_member <- NULL
  }
  
  # return preprocessed data
  d <- list(data = tb_data,
            anno = tb_anno,
            data_pseudocell = tb_pseudocell,
            anno_pseudocell_member = tb_pseudocell_member)
  return(d)
}


# embed and knn
embed_knn <- function(sobj, 
                      features,
                      do.scale = FALSE,
                      npcs = 20,
                      k.param = 20) {
  stopifnot(length(setdiff(features, rownames(sobj))) == 0)
  
  # subset features
  sobj <- sobj[features, ]
  
  # Scaling the data
  sobj <- Seurat::ScaleData(sobj, 
                            features = rownames(sobj), 
                            do.scale = do.scale)
  
  # Perform linear dimensional reduction
  sobj <- Seurat::RunPCA(sobj, 
                         features = rownames(sobj), 
                         npcs = npcs,
                         verbose = FALSE)
  
  # knn
  sobj <- Seurat::FindNeighbors(sobj,
                                dims = 1:npcs,
                                k.param = k.param)
  
  return(sobj)
}


#' recombine_pipeline
#'
#' @description
#' Streamlined RECOMBNE pipeline
#'
#' @param sobj Seurat object that includes LogNormalized data in layer sobj@assays$RNA$data.
#' @param subpop_name name of subpopulations that can be accessed in sobj@meta.data. The default is NULL (no subpopulation markers will be in the output).
#' @param fixed.hyperparameter fixed-parameter version (fRECOMBINE).
#' @param n_constraint constrained number of discriminant markers.
#' @param lambda0.fixed lambda0 value for fRECOMBINE.
#' @param lambda0s a sereis of lambda0 values for hyperparameter selection in RECOMBINE.
#' @param lambda1 lambda1 value.
#' @param nperms number of permutation for gap statistic calculation in RECOMBINE.
#' @param max_features max number of features before RECOMBINE analysis. If the total feature number exceeds this, HVGs will be selected to match max_features before RECOMBINE analysis.
#' @param pseudocell_embedding_nhvgs number of HVGs used in embedding during pseudo-cell construction.
#' @param pseudocell_embedding_npcs number of PCs used in embedding during pseudo-cell construction.
#' @param pseudocell_embedding_batchname name of batches. If not NULL, it needs to be a string and can be accessed in sobj@meta.data.
#' @param pseudocell_embedding_max.iter.harmony max number of iteration for running Harmony to correct batch effect indicated by pseudocell_embedding_batchname.
#' @param pseudocell_microcluster_res initial resolution of clustering during pseudo-cell construction.
#' @param pseudocell_microcluster_res.iter.rate iterative damping rate of clustering resolution during pseudo-cell construction.
#' @param allcell_discrim.marker_do.scale scale during contruction of KNN graph using discriminant markers.
#' @param allcell_discrim.marker_npcs number of PCs during contruction of KNN graph using discriminant markers.
#' @param allcell_discrim.marker_neighborhood_k.param K of the KNN graph using discriminant markers.
#' @param nhood_recur_fdr_threshod threshold for the statistical significance of neighborhood recurrence FDRs. The default is 0.05.
#'
#' @details
#' Given a Seurat object with log-normalized data and optionally subpopulation annotation, this pipeline output discriminant markers of the dataset and recurrent composite markers for subpopulations and individual cells.
#'
#' @return
#' \item{list of data frames}{including features with nonzero weights, recurrent composite markers of subpopulations, recurrent composite markers of individual cells, and gap statistic output from a series of lambda0s.}
#' 
#' @examples
#' recombine.out <- recombine_pipeline(sobj, subpop_name = "celltype")
#' 
#' @export
#'
recombine_pipeline <- function(
    sobj,
    subpop_name = NULL,
    fixed.hyperparameter = TRUE,
    n_constraint = 50,
    lambda0.fixed = 1,
    lambda0s = seq(0.001, 4000, length = 100),
    lambda1 = 0.0001,
    nperms = 10,
    max_features = 10000,
    pseudocell_embedding_nhvgs = 2000,
    pseudocell_embedding_npcs = 50,
    pseudocell_embedding_batchname = NULL,
    pseudocell_embedding_max.iter.harmony = 100,
    pseudocell_microcluster_res = 50,
    pseudocell_microcluster_res.iter.rate = 0.8,
    allcell_discrim.marker_do.scale = FALSE,
    allcell_discrim.marker_npcs = 20,
    allcell_discrim.marker_neighborhood_k.param = 20,
    nhood_recur_fdr_threshod = 0.05) 
{
  message("Running RECOMBINE pipeline")
  
  # data prep
  message("Preparing expression matrix and annotation")
  data_list <- data_prep(sobj,
                         subpop_name = subpop_name,
                         max_features = max_features,
                         pseudocell_embedding_nhvgs = pseudocell_embedding_nhvgs,
                         pseudocell_embedding_npcs = pseudocell_embedding_npcs,
                         pseudocell_embedding_batchname = pseudocell_embedding_batchname,
                         pseudocell_embedding_max.iter.harmony = pseudocell_embedding_max.iter.harmony,
                         pseudocell_microcluster_res = pseudocell_microcluster_res,
                         pseudocell_microcluster_res.iter.rate = pseudocell_microcluster_res.iter.rate)
  
  # expr matrix
  mt_expr <- as.matrix(data_list$data[, -1])
  rownames(mt_expr) <- data_list$data$cell
  
  # matrix used for SHC_SSL
  if (is.null(data_list$data_pseudocell)) {
    mt <- mt_expr
  } else {
    mt <- as.matrix(data_list$data_pseudocell[, -1])
    rownames(mt) <- data_list$data_pseudocell$pseudo_cell
  }
  
  # center for each feature
  mt <- scale(mt, scale = FALSE)
  
  # scale so that the overall variance is unit 
  mt <- mt/sd(as.double(mt))
  
  # run SHC_SSL
  message("Running SHC-SSL")
  if (fixed.hyperparameter) {
    res <- SHC_SSL(x = mt,
                   lambda0 = lambda0.fixed,
                   lambda1 = lambda1)
    w <- res$w
    names(w) <- colnames(mt)
    
    # w_nonzero
    w <- w[w > 0]
    tb <- dplyr::tibble(marker = names(w),
                        w = w) %>%
      dplyr::arrange(dplyr::desc(w))
    if (nrow(tb) > n_constraint) {
      tb <- tb[1:n_constraint, ]
    }
    tb_w <- tb
    gapstat_out <- NULL
  } else {
    out <- SelectCF_gapstat(x = mt,
                            nperms = nperms,
                            lambda0s = lambda0s,
                            lambda1 = lambda1)
    w <- out$result$w
    names(w) <- colnames(mt)
    
    # w_nonzero
    w <- w[w > 0]
    tb <- dplyr::tibble(marker = names(w),
                        w = w) %>%
      dplyr::arrange(dplyr::desc(w))
    tb_w <- tb
    gapstat_out <- out
  }
  
  # embed and knn
  message("Embedding cells using discriminant markers and contructing the KNN graph")
  markers <- tb_w$marker
  sobj <- embed_knn(sobj, 
                    markers,
                    do.scale = allcell_discrim.marker_do.scale,
                    npcs = allcell_discrim.marker_npcs,
                    k.param = allcell_discrim.marker_neighborhood_k.param) 
  
  # K nearest neighbors graph from Seurat
  mt_nn <- sobj@graphs$RNA_nn %>% as.matrix()
  
  k_cells <- sum(mt_nn[1, ] > 0)
  stopifnot(all(rowSums(mt_nn > 0) == k_cells))
  
  # check if cell names are consistent
  stopifnot(all.equal(rownames(mt_nn), colnames(mt_nn)))
  stopifnot(all.equal(rownames(mt_expr), rownames(mt_nn)))
  
  # Extract RCMs at the cell level using neighborhood recurrence test
  message("Extracting RCMs at the cell level using neighborhood recurrence test")
  tb_rcm_cell <- nhood_recur_test(mt_expr[, markers], mt_nn, verbose = TRUE)
  
  # Extract RCMs at the subpop level
  message("Extracting RCMs at the subpopulation level")
  if (!is.null(subpop_name)) {
    tb_rcm_subpop <- get_rcm_subpop(tb_rcm_cell, 
                                    data_list$anno,
                                    subpop_name = subpop_name,
                                    nhood_recur_fdr_threshod = nhood_recur_fdr_threshod)
  } else {
    tb_rcm_subpop <- NULL
  }
  
  message("Done.")
  
  # output
  recombine.out <- list(df_w = tb_w,
                        df_rcm_subpop = tb_rcm_subpop,
                        df_rcm_cell = tb_rcm_cell,
                        gapstat_out = gapstat_out)
  return(recombine.out)
}


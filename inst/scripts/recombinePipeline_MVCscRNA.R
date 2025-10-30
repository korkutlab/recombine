rm(list = ls())

library(tidyverse)
library(Seurat)
library(recombine)


# Here we use a scRNA-seq data of the mouse visual cortex from Allen Brain Atlas, which can be downloaded from: 
# http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985

# The following files are used after uncompressing the downloaded file:
# mouse_VISp_2018-06-14_exon-matrix.csv
# mouse_VISp_2018-06-14_genes-rows.csv
# mouse_VISp_2018-06-14_samples-columns.csv


# get expression matrix and cell annotation ------
# expression
tb <- read_csv("mouse_VISp_2018-06-14_exon-matrix.csv")
colnames(tb)[1] <- "gene_entrez_id"
tb_data <- tb

# gene names
tb <- read_csv("mouse_VISp_2018-06-14_genes-rows.csv")
stopifnot(all.equal(tb_data$gene_entrez_id, tb$gene_entrez_id))
tb_data$gene_symbol <- tb$gene_symbol

tb_data <- tb_data %>%
  select(-gene_entrez_id) %>%
  select(gene_symbol, everything())

# transpose
mt <- as.matrix(tb_data[, -1])
rownames(mt) <- tb_data$gene_symbol
mt <- t(mt)
tb_data <- tibble(cell = rownames(mt)) %>%
  bind_cols(as_tibble(mt))

# anno
tb <- read_csv("mouse_VISp_2018-06-14_samples-columns.csv")
stopifnot(all.equal(tb_data$cell, tb$sample_name))
tb <- tb %>%
  rename(cell = sample_name) %>%
  select(cell, class, subclass, cluster)
tb_anno <- tb

# filter low quality cells
low_quality = c('No Class', 'Low Quality')
tb_anno <- tb_anno %>%
  filter(!(class %in% low_quality) &
           !(subclass %in% low_quality) &
           !(cluster %in% low_quality))

tb_data <- tb_data %>%
  filter(cell %in% tb_anno$cell)


# generate log-normalized data --------
# Initialize the Seurat object
mt <- as.matrix(tb_data[, -1])
rownames(mt) <- tb_data$cell
mt <- t(mt)

sobj <- CreateSeuratObject(counts = mt,
                           project = "mvc",
                           min.features = 0,
                           min.cells = 0)

# add cell annotation
stopifnot(all.equal(rownames(sobj@meta.data), tb_anno$cell))
sobj@meta.data <- sobj@meta.data %>%
  cbind(tb_anno[, -1] %>% as.data.frame())

# Normalizing the data
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)


# run RECOMBINE pipeline --------
recombine.out <- recombine_pipeline(sobj, subpop_name = "subclass")

# discriminant markers and weights
recombine.out$df_w

# recurrent composite markers of cell subpopulations
recombine.out$df_rcm_subpop %>%
  filter(fract_signif_cells > 0.5 & avg_nhood_zscore > 2 & PR_AUC > 0.5) %>%
  group_by(subclass) %>%
  slice_max(avg_nhood_zscore, n = 3) %>%
  ungroup()

recombine.out %>%
  saveRDS("recombine.out.rds")


sessionInfo()

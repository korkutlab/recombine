
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RECOMBINE

## Overview

RECOMBINE (REcurrent COmposite Markers for Biological Identities with
Neighborhood Enrichment) is a computational framework for unbiased
selection of discriminant markers that hierarchically distinguish cells
and extraction of recurrent composite markers (RCMs) that characterize
cell types and states with high granularity.

RECOMBINE consists of two steps. The first step employs sparse
hierarchical clustering with spike-and-slab lasso (SHC-SSL) to select
features that discriminate hierarchical cell subpopulations from high
dimensional data such as scRNA-seq. The second step employs neighborhood
recurrence test to extract RCMs at the cell and subpopulation levels.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("korkutlab/recombine")
```

## Example: Running the Streamlined RECOMBINE Pipeline

Here we use a scRNA-seq data of the mouse visual cortex from Allen Brain
Atlas, which can be downloaded from:
<http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985>

The following data will be used after uncompressing the downloaded file:
(1) mouse_VISp_2018-06-14_exon-matrix.csv (2)
mouse_VISp_2018-06-14_genes-rows.csv (3)
mouse_VISp_2018-06-14_samples-columns.csv

``` r
library(tidyverse)
library(Seurat)
library(recombine)
```

First, we extract expression matrix and cell annotation from the above
files, and filter out low-quality cells.

``` r
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
```

Next, we create a Seurat object that contains LogNormalized data and
meta data, which will be the data input of RECOMBINE analysis.

``` r
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
sobj <- NormalizeData(sobj, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
```

Now, the data input has been prepared and we are ready to run the
RECOMBINE pipeline. For simplicity and saving computation time, here we
run fix-hyperparameter RECOMBINE to select top 50 discriminant markers.
To select the full set of discriminant markers, please adjust the
parameters in the recombine_pipeline function or see more examples in
Further Tutorials: Step-by-Step RECOMBINE Workflows.

``` r
# run RECOMBINE pipeline
recombine.out <- recombine_pipeline(sobj, 
                                    subpop_name = "subclass",
                                    n_constraint = 50)

# discriminant markers and weights
recombine.out$df_w
#> # A tibble: 50 × 2
#>    marker     w
#>    <chr>  <dbl>
#> 1  Vip    0.404
#> 2  Gad1   0.264
#> 3  Sst    0.196
#> 4  Synpr  0.193
#> 5  Npy    0.189
#> 6  Cck    0.178
#> 7  Pcp4   0.146
#> 8  Nrgn   0.142
#> 9  Rab3b  0.127
#> 10 Gad2   0.124
#> # ℹ 40 more rows

# recurrent composite markers of cell subpopulations
recombine.out$df_rcm_subpop %>%
  filter(fract_signif_cells > 0.5 & avg_nhood_zscore > 2 & PR_AUC > 0.5) %>%
  group_by(subclass) %>%
  slice_max(avg_nhood_zscore, n = 3) %>%
  ungroup()
#> # A tibble: 28 × 6
#>    subclass marker        fract_signif_cells avg_nhood_zscore PR_AUC ROC_AUC
#>    <chr>    <chr>                      <dbl>            <dbl>  <dbl>   <dbl>
#> 1  Astro    Apoe                       1                22.2   0.704   0.996
#> 2  Astro    Aldoc                      1                19.3   1       1
#> 3  Astro    Slc1a2                     1                13.2   1       1
#> 4  L2/3 IT  Arpp19                     0.999             5.47  0.629   0.950
#> 5  L2/3 IT  Olfm1                      0.999             4.84  0.678   0.951
#> 6  L2/3 IT  Chn1                       1                 4.72  0.532   0.957
#> 7  L4       Rgs4                       1                 7.43  0.779   0.985
#> 8  L4       Arpp21                     1                 6.43  0.552   0.960
#> 9  L4       1110008P14Rik              1                 5.67  0.639   0.943
#> 10 L6 CT    3110035E14Rik              1                 8.66  0.641   0.984
#> # ℹ 18 more rows
```

## Further Tutorials: Step-by-Step RECOMBINE Workflows

Tutorial 1: RECOMBINE identifies recurrent composite markers for common
and rare cell subpopulations of mouse intestinal organoids
(<https://github.com/korkutlab/recombine/blob/main/inst/tutorials/intestine_scRNA/intestine_scRNA.md>)

Tutorial 2: RECOBMINE selects unbiased marker panel from scRNA-seq data
for targeted spatial transcriptomics
(<https://github.com/korkutlab/recombine/blob/main/inst/tutorials/mvc_scRNA_spatial/mvc_scRNA_spatial.md>)

Tutorial 3: RECOMBINE SHC-SSL performs accurate feature selection
coupled with hierarchical clustering
(<https://github.com/korkutlab/recombine/blob/main/inst/tutorials/sim/sim.md>)

<!-- ## Citation -->
<!-- Li X. and Korkut A. (2023) Recurrent composite markers of cell types and states.  -->

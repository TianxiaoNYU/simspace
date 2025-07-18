#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cells_meta_path <- args[1]
cells_omics_path <- args[2]
group_col <- args[3]
spatial_x <- args[4]
spatial_y <- args[5]
seed <- as.numeric(args[6])
if (length(args) != 6) {
  stop("Usage: Rscript scdesign.R <cells_meta_path> <cells_omics_path> <group_col> <spatial_x> <spatial_y> <seed>")
}

# Load required packages
suppressPackageStartupMessages({
  # options(pkgType = "binary")
  # if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  #   if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #     message("Installing 'BiocManager' package...")
  #     install.packages("BiocManager", repos = "https://cloud.r-project.org")
  #     BiocManager::install(version = "3.19")
  #   }
  #   if (!requireNamespace("remotes", quietly = TRUE)) {
  #     install.packages("remotes", repos = "https://cloud.r-project.org")
  #   }
  #   message("Installing 'SparseArray' package from Bioconductor...")
  #   BiocManager::install("SparseArray")
  #   message("Installing 'SingleCellExperiment' package from Bioconductor...")
  #   BiocManager::install("SingleCellExperiment")
  # }
  # if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  #   stop("Package 'SingleCellExperiment' could not be installed or loaded.")
  # }
  library(SingleCellExperiment)

  # if (!requireNamespace("scDesign3", quietly = TRUE)) {
  #   message("Installing 'scDesign3' package from GitHub...")
  #   remotes::install_github("SONGDONGYUAN1994/scDesign3")
  # }
  # if (!requireNamespace("scDesign3", quietly = TRUE)) {
  #   stop("Package 'scDesign3' could not be installed or loaded.")
  # }
  library(scDesign3)
})

# Read the input data
cat('Reading input data...\n')
cells_meta <- read.csv(cells_meta_path, header = TRUE, row.names = 1)
cells_omics <- read.csv(cells_omics_path, header = TRUE, row.names = 1)
if (nrow(cells_meta) != ncol(cells_omics)) {
  if (nrow(cells_meta) == nrow(cells_omics)) {
    cells_omics <- as.data.frame(t(cells_omics))
  } else {
  stop("Number of cells in cells_meta and cells_omics do not match.")
  }
}

if (group_col %in% colnames(cells_meta)) {
  cells_meta$celltype <- cells_meta[[group_col]]
} else {
  stop(paste("Column", group_col, "not found in cells_meta"))
}

if (spatial_x %in% colnames(cells_meta)) {
  cells_meta$row <- cells_meta[[spatial_x]]
} else {
  stop(paste("Column", spatial_x, "not found in cells_meta"))
}
if (spatial_y %in% colnames(cells_meta)) {
  cells_meta$col <- cells_meta[[spatial_y]]
} else {
  stop(paste("Column", spatial_y, "not found in cells_meta"))
}

sce <- SingleCellExperiment(list(counts = cells_omics), colData = cells_meta)
counts(sce) <- as.matrix(counts(sce))
logcounts(sce) <- log1p(counts(sce))

new_meta_raw <- read.csv('./tmp/new_meta.csv', header = TRUE)
new_meta <- new_meta_raw[,c('row', 'col', 'fitted_celltype')]
colnames(new_meta)[3] <- 'celltype'
new_meta$corr_group <- new_meta$celltype

set.seed(seed)

cat('Running scDesign3 simulation...\n')
{
  example_data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "celltype",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = c('row', 'col'), 
    corr_by = "celltype"
  )
  example_marginal <- fit_marginal(
    data = example_data,
    predictor = "gene",
    mu_formula = "celltype",
    sigma_formula = "celltype",
    family_use = "nb",
    n_cores = 6,
    usebam = FALSE,
    parallelization = "pbmcmapply"
  )
  example_copula <- fit_copula(
    sce = sce,
    assay_use = "counts",
    marginal_list = example_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 6,
    input_data = example_data$dat
  )

  example_para <- extract_para(
    sce = sce,
    marginal_list = example_marginal,
    n_cores = 6,
    family_use = "nb",
    new_covariate = new_meta, 
    data = example_data$dat
  )
}

example_newcount <- simu_new(
    sce = sce,
    mean_mat = example_para$mean_mat,
    sigma_mat = example_para$sigma_mat,
    zero_mat = example_para$zero_mat,
    quantile_mat = NULL,
    copula_list = example_copula$copula_list,
    n_cores = 6,
    family_use = "nb",
    input_data = example_data$dat,
    new_covariate = new_meta,
    important_feature = example_copula$important_feature,
    filtered_gene = example_data$filtered_gene
  )

all_sim_counts <- example_newcount
all_sim_meta <- new_meta

all_sim_counts <- as.data.frame(all_sim_counts)
all_sim_meta <- as.data.frame(all_sim_meta)
all_sim_meta$corr_group <- NULL
all_sim_meta$state <- new_meta$state

# Reindex the dataframes
rownames(all_sim_meta) <- paste0("cell_", seq_len(nrow(all_sim_meta)))
colnames(all_sim_counts) <- rownames(all_sim_meta)

cat('Saving simulated counts and metadata to CSV files...\n')
write.csv(all_sim_counts, './tmp/simulated_data.csv', row.names = TRUE, quote = FALSE)
write.csv(all_sim_meta, './tmp/simulated_meta.csv', row.names = TRUE, quote = FALSE)

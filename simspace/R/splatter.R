#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cells_meta_path <- args[1]
n_gene <- as.numeric(args[2])
if (length(args) != 2) {
  stop("Usage: Rscript splatter.R <cells_meta_path> <n_gene>")
}

# Load required packages
suppressPackageStartupMessages({
  if (!requireNamespace("splatter", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("splatter")
  }
  library(splatter)
})

# Read the input data
cells_meta <- read.csv(cells_meta_path)
n_cells <- nrow(cells_meta)
# Calculate the proportion of each cell state
cell_state_counts <- table(cells_meta$state)
group_probs <- as.numeric(cell_state_counts / sum(cell_state_counts))
group_probs <- group_probs / sum(group_probs)

set.seed(1)
params <- newSplatParams()
params <- setParam(params, "nGenes", n_gene)  # Number of genes
params <- setParam(params, "batchCells", n_cells)  # Total number of cells
params <- setParam(params, "group.prob", group_probs)  # Equal probability for 8 cell types
params <- setParam(params, "de.prob", 0.6)  # % of genes are differentially expressed

# Simulate data
sim_data <- splatSimulate(params, method = "groups")

# Check if counts(sim_data) is a dgCMatrix
if (inherits(counts(sim_data), "dgCMatrix")) {
  counts_matrix <- as.matrix(counts(sim_data))
} else {
  counts_matrix <- counts(sim_data)
}

write.csv(as.data.frame(counts_matrix), './tmp/simulated_data.csv', row.names = TRUE, quote = F)
write.csv(as.data.frame(colData(sim_data)), './tmp/simulated_meta.csv', row.names = TRUE, quote = F)
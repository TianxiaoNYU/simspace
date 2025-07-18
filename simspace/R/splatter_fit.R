#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cells_meta_path <- args[1]
cells_omics_path <- args[2]
group_col <- args[3]
n_cell <- as.numeric(args[4])

# Load required packages
suppressPackageStartupMessages({
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("SingleCellExperiment")
  }
  library(SingleCellExperiment)
  
  if (!requireNamespace("splatter", quietly = TRUE)) {
    BiocManager::install("splatter")
  }
  library(splatter)
})

# Read the input data
cat('Reading input data...\n')
cells_meta <- read.csv(cells_meta_path, header = TRUE, row.names = 1)
cells_omics <- read.csv(cells_omics_path, header = TRUE, row.names = 1)
if (nrow(cells_meta) != ncol(cells_omics)) {
  stop("Number of cells in cells_meta and cells_omics do not match.")
}

if (group_col %in% colnames(cells_meta)) {
  cells_meta$state <- cells_meta[[group_col]]
} else {
  stop(paste("Column", group_col, "not found in cells_meta"))
}

n_groups <- length(unique(cells_meta$state))

# Group the metadata and omics data by state
grouped_data <- split(seq_len(nrow(cells_meta)), cells_meta$state)

# Initialize a list to store splatter fits
splatter_fits <- list()

cat('Fitting splatter to each group...\n')
# Fit splatter to each group
for (state in names(grouped_data)) {
  group_indices <- grouped_data[[state]]
  group_omics <- cells_omics[, group_indices, drop = FALSE]
  group_omics <- as.matrix(group_omics)
  
  # Check if the group_omics matrix is empty
  if (nrow(group_omics) == 0) {
    cat(paste("Group", state, "has no data. Skipping...\n"))
    next
  }
  # Check if the group_omics matrix has only one row
  if (nrow(group_omics) == 1) {
    cat(paste("Group", state, "has only one cell. Skipping...\n"))
    next
  }

  # Estimate parameters using splatter
  params <- splatEstimate(group_omics)
  params <- setParam(params, "batchCells", n_cell) # Update the number of cells
  
  sim <- splatSimulate(params)

  if (inherits(counts(sim), "dgCMatrix")) {
    counts_matrix <- as.matrix(counts(sim))
  } else {
    counts_matrix <- counts(sim)
  }

  sim_counts <- as.data.frame(counts_matrix)
  sim_meta <- as.data.frame(colData(sim))
  splatter_fits[[state]] <- list(
    params = params,
    sim_counts = sim_counts,
    sim_meta = sim_meta
  )
}

# Combine sim_counts from all groups into one big dataframe
cat('Combining simulated counts into one dataframe...\n')
combined_sim_counts <- do.call(cbind, lapply(splatter_fits, function(fit) fit$sim_counts))
# Add group information as a column
# group_labels <- unlist(lapply(names(splatter_fits), function(state) {
#   rep(state, ncol(splatter_fits[[state]]$sim_counts))
# }))
# combined_sim_counts <- cbind(Group = group_labels, combined_sim_counts)
# Print a message indicating completion
write.csv(combined_sim_counts, './tmp/simulated_data.csv', row.names = TRUE, quote = FALSE)
cat('Simulated counts combined successfully.\n')



# Combine sim_meta from all groups into one big dataframe
cat('Combining simulated metadata into one dataframe...\n')
combined_sim_meta <- do.call(rbind, lapply(splatter_fits, function(fit) fit$sim_meta))
combined_sim_meta$Group <- unlist(lapply(names(splatter_fits), function(state) {
  rep(state, nrow(splatter_fits[[state]]$sim_meta))
}))
combined_sim_meta <- combined_sim_meta[order(combined_sim_meta$Group), ]
# Print a message indicating completion
cat('Simulated metadata combined successfully.\n')
# Save the combined simulated counts and metadata to CSV files
write.csv(combined_sim_meta, './tmp/simulated_meta.csv', row.names = TRUE, quote = FALSE)
cat('Saving simulated counts and metadata to CSV files...\n')
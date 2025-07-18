#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
cells_meta_path <- args[1]
cells_omics_path <- args[2]
group_col <- args[3]
spatial_x <- args[4]
spatial_y <- args[5]
n_rep <- as.numeric(args[6])
seed <- as.numeric(args[7])
if (length(args) != 7) {
  stop("Usage: Rscript SRTsim.R <cells_meta_path> <cells_omics_path> <group_col> <spatial_x> <spatial_y> <n_rep> <seed>")  
}

if (n_rep <= 0) {
  stop("n_rep must be a positive integer.")
}

suppressPackageStartupMessages({
  if (!requireNamespace("SRTsim", quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github("xzhoulab/SRTsim")
  }
  library(SRTsim)
})

# Read the input data
cat('Reading input data...\n')
cells_meta <- read.csv(cells_meta_path, header = TRUE, row.names = 1)
cells_omics <- read.csv(cells_omics_path, header = TRUE, row.names = 1)
if (nrow(cells_meta) != ncol(cells_omics)) {
  stop("Number of cells in cells_meta and cells_omics do not match.")
}
rownames(cells_meta) <- paste0("Cell_", rownames(cells_meta))
colnames(cells_omics) <- rownames(cells_meta)
cells_omics <- as.matrix(cells_omics)

example_loc <- cells_meta[,c(spatial_x, spatial_y, group_col)]
colnames(example_loc) <- c("x", "y", "label")

simSRT <- createSRT(
  count_in=cells_omics,
  loc_in =example_loc)

set.seed(seed)
## Estimate model parameters for data generation
simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")

for(i in 1:n_rep){
  simSRT1 <- srtsim_count(simSRT1, verbose = T, numCores = 4)
  example_newcount <- as.data.frame(as.matrix(simSRT1@simCounts))
  example_newmeta <- simSRT1@simcolData

  cell_index_range <- seq(1, nrow(example_newmeta)) + 
    (i - 1) * nrow(example_newmeta)
  rownames(example_newmeta) <- paste0("Cell_", cell_index_range)
  colnames(example_newcount) <- rownames(example_newmeta)

  if (i == 1) {
    all_sim_counts <- example_newcount
    all_sim_meta <- example_newmeta
  } else {
    all_sim_counts <- cbind(all_sim_counts, example_newcount)
    all_sim_meta <- rbind(all_sim_meta, example_newmeta)
  }
}

all_sim_counts <- as.data.frame(all_sim_counts)
all_sim_meta <- as.data.frame(all_sim_meta)

# Reindex the dataframes
rownames(all_sim_meta) <- paste0("cell_", seq_len(nrow(all_sim_meta)))
colnames(all_sim_counts) <- rownames(all_sim_meta)

cat('Saving simulated counts and metadata to CSV files...\n')
write.csv(all_sim_counts, './tmp/simulated_data.csv', row.names = TRUE, quote = FALSE)
write.csv(all_sim_meta, './tmp/simulated_meta.csv', row.names = TRUE, quote = FALSE)
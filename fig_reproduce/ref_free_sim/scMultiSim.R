library(scMultiSim)
save_path <- "./"

# Define ligand-receptor interaction parameters, which is used 
# in the original paper to simulate cell-cell interactions
lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(5.2, 5.9)
)

# Load predefined GRN parameters
data(GRN_params_100)

# Run simulations with different random seeds.
# All other parameters are kept the same as in the original paper
for(i in 0:8){
  set.seed(i)
  options_ <- list(
    GRN = GRN_params_100,
    speed.up = TRUE,
    num.genes = 120,
    num.cells = 300,
    num.cifs = 20,
    cif.sigma = 0.2,
    tree = Phyla3(),
    intrinsic.noise = 0.5,
    cci = list(
      params = lig_params,
      max.neighbors = 4,
      grid.size = 30,
      cell.type.interaction = "random",
      step.size = 0.5
    )
  )
  results <- sim_true_counts(options_)
  
  # Access the simulated spatial data
  coordinates <- results$cci_locs
  spatial_data <- results$cell_meta
  spatial_data$x <- coordinates[,1]
  spatial_data$y <- coordinates[,2]
  
  # Save the spatial data to a CSV file
  write.csv(spatial_data, 
            paste0(save_path, "spatial_sim_", i, ".csv"),
            quote = F)
}














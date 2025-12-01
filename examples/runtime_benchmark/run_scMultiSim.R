args <- commandArgs(trailingOnly = TRUE)
nrow <- as.numeric(args[1])
ncol <- as.numeric(args[2])
seed <- as.numeric(args[3])

library(scMultiSim)
save_path <- "/Users/zhaotianxiao/Library/CloudStorage/Dropbox/FenyoLab/Project/Spatialsim/SimSpace/fig_data/scMultiSim/"

lig_params <- data.frame(
  target    = c(101, 102),
  regulator = c(103, 104),
  effect    = c(5.2, 5.9)
)

data(GRN_params_100)
set.seed(seed)

options_ <- list(
  GRN = GRN_params_100,
  speed.up = TRUE,
  num.genes = 220,
  num.cells = nrow * ncol,
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














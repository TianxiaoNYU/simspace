# 📘 SimSpace Tutorials

Welcome to the SimSpace tutorials! This guide walks you through two key modes of spatial omics simulation using SimSpace:

- **Reference-based simulation**: Generate synthetic datasets that mimic the spatial organization of a real reference dataset.
- **Reference-free simulation**: Create synthetic spatial data using randomized parameters and spatial priors without needing real data.

---

## 🧭 Getting Started

Before running the tutorials, make sure you have a proper python environment installed SimSpace:

```bash
pip install simspace  # or pip install -e . if installing from source
```

(Optional) To simulate omics profiles using R packages like scDesign3, splatter, or SRTsim, please follow the R environment setup guide in README.

## Reference-Free Simulation

In this mode, spatial structure and cell types are generated from scratch based on a Markov Random Field (MRF) model and user-defined parameters. Ideal for testing robustness and benchmarking methods under controlled conditions.

### 📂 Inputs
- Number of spatial niches
- Number of cell types
- Shape of spatial grid
- Neighborhood structure
- (Optional) MRF model parameters if one want design their own spatial patterns

### 🛠️ Steps
1. Initialize the parameters and SimSpace object
2. Run MRF-based simulation
3. Visualize the spatial simulation results
4. Simulate omics profiles (if needed)

### 📌 Example

👉 See [examples/reference_free.ipynb](examples/reference_free.ipynb) for code.

## Reference-Based Simulation

In this setting, you provide a real spatial omics dataset as the reference. SimSpace will generate synthetic datasets that preserve the spatial layout and cell type composition of the input.

### 📂 Inputs
- Inputs as in reference-free simulation
- Reference data
  - A cell metadata with spatial coordinates and cell types
  - (Optional) Omics profiles to use as templates

### 🛠️ Steps
1. Load the reference dataset
2. Fit SimSpace spatial parameters from the reference
3. Run SimSpace in reference-based mode
4. Simulate omics profiles (if needed)

### 📌 Example

👉 See [examples/reference_based.ipynb](examples/reference_based.ipynb) and [examples/spatial_fitting.ipynb](examples/spatial_fitting.ipynb) for code.


## 🙋 Need Help?

For issues or feature requests, open an issue on GitHub.
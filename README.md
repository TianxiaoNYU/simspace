# SimSpace

**SimSpace** is a Python framework for simulating spatial omics data with realistic cellular distributions and tissue organization. Designed for benchmarking spatial analysis methods, SimSpace enables generation of synthetic datasets that preserve spatial autocorrelation, cell-cell interactions, and spatial proximities using a Markov Random Field (MRF) model.

## 📦 Installation

To install the latest version of SimSpace, we recommend using conda to setup the environment:

```bash
git clone https://github.com/TianxiaoNYU/simspace.git
```

- Create a conda environment for simspace
```bash
cd simspace
conda env create -f environment.yml
conda activate simspace
```

- Install simspace from PyPi
```bash
pip install simspace
```

### 🧬 Optional: Setting Up the R Environment for Omics Simulation

Besides built-in functions, SimSpace also supports external omics profile simulation via R-based tools including **scDesign3**, **SRTsim**, and **splatter**. You can either install these packages manually or use the [`renv`](https://rstudio.github.io/renv/) package to recreate the exact R environment used by SimSpace.

#### Steps:

1. Ensure that **R (version 4.4 or compatible)** is installed on your system.
2. Navigate to the R project folder and restore the environment:

```bash
cd simspace/R
Rscript -e 'install.packages("renv"); renv::restore()'
```
This will install all required R dependencies in a reproducible, isolated environment.

## 📘 Tutorials

To get started with SimSpace, we provide detailed tutorials covering both reference-based and reference-free simulation modes.

- **Step-by-step tutorials** can be found in [`tutorials.md`](./tutorials.md)
- **Executable notebook examples** are located in the [`examples/`](./examples/) directory

These resources walk through how to configure and run simulations as well as visualize outputs.

## 🚀 Quick Start

Here’s a basic example to simulate a 2D tissue with 3 spatial niches and 8 cell types:

```python
from simspace import util, spatial

# Define simulation parameters
params = util.generate_random_parameters(
    n_group=3,
    n_state=8,
    seed=42)

# Run simulation
sim = util.sim_from_params(
    params,
    shape=(50, 50),    # shape of the simulation grid
    custom_neighbor=spatial.generate_offsets(3, 'manhattan'),
    seed=42
)

# Visualize
sim.plot()
```

## 🙋‍♀️ About

Developed by Tianxiao Zhao at NYU Grossman School of Medicine. Should you have any questions, please contact Tianxiao Zhao at Tianxiao.Zhao@nyulangone.org

## 🔗 References
If you use SimSpace in your work, please cite the work on BioRxiv:



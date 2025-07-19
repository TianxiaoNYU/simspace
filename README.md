# SimSpace

**SimSpace** is a Python framework for simulating spatial omics data with realistic cellular distributions and tissue organization. Designed for benchmarking spatial analysis methods, SimSpace enables generation of synthetic datasets that preserve spatial autocorrelation, cell-cell interactions, and reference-based spatial layouts using a Markov Random Field (MRF) model.

## 📦 Installation

To install the latest version of SimSpace, we recommend using conda to setup the environment:

```bash
git clone https://github.com/TianxiaoNYU/simspace.git
cd simspace
# Create a conda environment for simspace
conda env create -f environment.yml
conda activate simspace
# Install simspace from PyPi
pip install simspace
```

### 🧬 Setting up the R environment

To reproduce the R environment required for SimSpace omics simulation:

1. Make sure you have R (version 4.4 or compatible) installed.
2. From the project folder:

```bash
cd simspace/R
Rscript -e 'install.packages("renv"); renv::restore()'
```

## 🚀 Quick Start

Here’s a basic example to simulate a 2D tissue with 4 cell types:

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
    shape=(50, 50),
    num_iteration=4, 
    n_iter=6, 
    custom_neighbor=spatial.generate_offsets(3, 'manhattan'),
    seed=42
)

# Visualize
sim.plot()
```

## 🙋‍♀️ About

Developed by Tianxiao Zhao at NYU Grossman School of Medicine. Should you have any questions, please contact Tianxiao Zhao at Tianxiao.Zhao@nyulangone.org

## 🔗 References
If you use SimSpace in your work, please cite:



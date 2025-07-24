# SimSpace

[![Docs](https://readthedocs.org/projects/simspace/badge/?version=latest)](https://simspace.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/simspace.svg)](https://pypi.org/project/simspace/)

**SimSpace** is a Python framework for simulating spatial omics data with realistic cellular distributions and tissue organization. Designed for benchmarking spatial analysis methods, SimSpace enables generation of synthetic datasets that preserve spatial autocorrelation, cell-cell interactions, and spatial proximities using a Markov Random Field (MRF) model.

![SimSpace Workflow](images/overview.png)

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

SimSpace supports external omics profile simulation via R-based tools, including **scDesign3**, **SRTsim**, and **splatter**. These tools are optional but recommended if you want to simulate gene expression profiles in addition to spatial patterns.

To enable this functionality, please install the required R packages manually in your system R environment:

Steps:
1.	Ensure that R (version 4.4 or compatible) is installed on your system. You can download it from CRAN.
2.	Open an R session and install the required packages:
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
devtools::install_github("xzhoulab/SRTsim")
```
```R
if (!require("devtools", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("splatter"))
```

Once installed, SimSpace will automatically use these tools when relevant R-based simulations are requested.

## 📖 Documentation

Full documentation for SimSpace is available at:

➡️ [simspace.readthedocs.io](https://simspace.readthedocs.io/en/latest/)

The documentation includes:

- Installation instructions (Python + optional R setup)
- Tutorials (reference-free & reference-based simulation)
- API reference for all modules and functions

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

# Check and save the simulated spatial data
sim.meta.head()
# sim.meta.to_csv('simspace.csv')
```

## 🙋‍♀️ About

Developed by Tianxiao Zhao at NYU Grossman School of Medicine. Should you have any questions, please contact Tianxiao Zhao at Tianxiao.Zhao@nyulangone.org

## 🔗 References
If you use SimSpace in your work, please cite the work on BioRxiv:
[Zhao T, Zhang K, Hollenberg M, Zhou W, Fenyo D. SimSpace: a comprehensive in-silico spatial omics data simulation framework. bioRxiv. 2025:2025.07.18.665587.](https://www.biorxiv.org/content/10.1101/2025.07.18.665587v1)



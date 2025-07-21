# 👋 Welcome to SimSpace Documentation

**SimSpace** is a flexible and extensible simulation framework for generating spatial omics data with customizable spatial structures, cell types, and gene expression profiles. It supports both reference-free and reference-based simulation modes, allowing users to create diverse tissue architectures ranging from well-defined niches to spatially intermixed environments across resolutions. Moreover, SimSpace also supports three-dimensional simulations, capturing the full complexity of tissue structure.

This documentation covers:

- ✅ Installation instructions (Python & R environments)
- 🧪 Tutorials for reference-free and reference-based simulation
- 🛠 API Reference for all SimSpace modules and functions

---

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
- Ensure that R (version 4.4 or compatible) is installed on your system. You can download it from CRAN.
- Open an R session and install the required packages:

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

---

## 🔗 Useful Links

- 📦 [GitHub Repository](https://github.com/TianxiaoNYU/SimSpace)
- 🐍 [PyPI Package](https://pypi.org/project/simspace/)

---

Feel free to explore the sidebar for detailed module-level documentation and usage examples. If you encounter any issues or have suggestions, you're welcome to open an issue on GitHub.
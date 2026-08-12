# SimSpace tutorials

SimSpace supports two primary workflows:

- **Reference-free simulation** generates spatial organization and molecular
  profiles from user-defined or randomly drawn parameters.
- **Reference-based simulation** calibrates or reuses spatial parameters from
  cell-resolved reference metadata, generates new coordinates, and optionally
  fits molecular profiles from a reference matrix.

## Start with the command line

The main [README](README.md) contains complete `simspace reference-free` and
`simspace reference-based` commands, a parameter guide, output definitions,
and background-execution examples. This is the shortest route to a complete
simulation and does not require a notebook.

```bash
simspace reference-free --help
simspace reference-based --help
```

## Step-by-step notebooks

The notebooks remain available for users who want to inspect and modify each
Python step:

- [Reference-free simulation](tutorials/reference_free.ipynb)
- [Reference-based Xenium simulation](tutorials/reference_based_Xenium.ipynb)
- [Reference-based CODEX simulation](tutorials/reference_based_CODEX.ipynb)
- [Spatial parameter fitting](tutorials/spatial_fitting.ipynb)
- [Reference-free 3D simulation](tutorials/3D_simulation.ipynb)
- [Domain design](tutorials/domain_design.ipynb)

Install SimSpace before running a notebook:

```bash
pip install simspace
```

The scDesign3, SRTsim, and Splatter molecular-profile adapters additionally
require their corresponding packages to be installed in the normal local R
environment, as described in the main README. No separate R virtual
environment is needed.

For questions or feature requests, open an issue in the SimSpace GitHub
repository.

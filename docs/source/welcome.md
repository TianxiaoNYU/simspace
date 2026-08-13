# Welcome to SimSpace

SimSpace generates spatial omics data with controllable tissue niches,
cell-type organization, and molecular profiles. Version 0.4.1 supports both a
Python API and a headless command-line interface.

## Install

```bash
pip install simspace
```

## Run without a notebook

```bash
simspace reference-free \
  --output-dir results/reference_free \
  --shape 100 100 \
  --n-niches 3 \
  --n-cell-types 8 \
  --n-genes 100 \
  --seed 42
```

```bash
simspace reference-based \
  --reference-metadata data/reference_metadata.csv \
  --cell-type-column Cluster \
  --x-column x_centroid \
  --y-column y_centroid \
  --params data/fitted_params.json \
  --profile-model native \
  --output-dir results/reference_based \
  --seed 42
```

Both commands write `metadata.csv`, `molecular_profiles.csv.gz`,
`spatial_scatter.png`, and `parameters.json`. They are noninteractive and can
be launched with `nohup`, a shell script, or a cluster scheduler. Run
`simspace reference-free --help` or `simspace reference-based --help` for the
complete parameter list.

The tutorials in the documentation sidebar explain the corresponding Python
steps. They are optional and are not required for command-line execution.
Optional scDesign3, SRTsim, and Splatter adapters use packages installed in the
normal local R environment; no separate R virtual environment is required.

## Links

- [GitHub repository](https://github.com/TianxiaoNYU/simspace)
- [PyPI package](https://pypi.org/project/simspace/)
- [Full CLI and output guide](https://github.com/TianxiaoNYU/simspace#command-line-quick-start)

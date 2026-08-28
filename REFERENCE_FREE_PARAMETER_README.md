# Designing reference-free SimSpace parameter JSON files

This guide explains how to design a readable, reproducible spatial parameter
file for a reference-free SimSpace simulation. It also explains which parts of
the cortex-inspired design in manuscript Figure 3 are controlled by the normal
SimSpace `parameters.json` file and which require custom geometry in Python.

## The important distinction

The JSON accepted by `simspace reference-free --params` controls the two MRFs
and cell retention. It contains:

- the numbers of niches and cell types;
- affinities among niches;
- one cell-type affinity design per niche;
- one retention probability per cell type; and
- the overall MRF smoothness parameter.

Grid shape, neighborhood distance, Gibbs sweeps, coordinate jitter, random
seed, and molecular-profile settings are CLI or Python arguments. Cell-type
names, colors, tissue masks, and manually drawn niche boundaries are also not
part of the standard parameter JSON.

Figure 3 uses reference-free SimSpace, but its four curved bands are created by
a custom geometry function and assigned to `simulation.niche` before the
cell-type MRF is run. Therefore, the archived Figure 3
`cortex_parameters.json` is a complete provenance manifest, not a file that can
be supplied directly to `simspace reference-free --params`.

## Start with a design, not with numbers

A good spatial simulation usually has a short biological or visual design
statement. For example:

> Make three contiguous tissue domains. Cell types 1 and 2 form a superficial
> module, types 2 and 3 form a transition module, and types 3 and 4 form a deep
> module. Retain all four types at similar rates.

Turn that statement into parameters in this order:

1. **Tissue-scale organization:** decide how many niches are needed and which
   niches should form separate domains.
2. **Local cell organization:** for each niche, decide which cell types may mix
   and which should avoid one another.
3. **Density:** choose cell-type retention probabilities that give an adequate
   number of displayed cells.
4. **Texture:** choose the MRF scale, neighborhood, and number of finite Gibbs
   sweeps.
5. **Presentation:** add modest coordinate jitter and use fixed, accessible
   colors when plotting.

This hierarchy is the main lesson of Figure 3: geometry establishes the large
visual structure, niche-conditioned cell-type design establishes biological
organization, and the MRF adds local texture.

## Standard JSON schema

Every top-level field must have a `value` wrapper because this is the format
read by `simspace.util.load_params()` and written by
`simspace.util.save_params()`.

| Field | Meaning | Required size |
|---|---|---|
| `n_group` | Number of niches, `M` | one integer |
| `n_state` | Number of cell types, `K` | one integer |
| `niche_theta` | Off-diagonal niche affinities | `M(M-1)/2` values |
| `theta_list` | Off-diagonal cell-type affinities, one vector per niche | `M` vectors of `K(K-1)/2` values |
| `density_replicates` | Cell-type retention probabilities | `K` values in `[0, 1]` |
| `phi_replicates` | Overall MRF spatial-dependence scale | one positive number |

SimSpace reconstructs a symmetric matrix from each affinity vector and sets
its diagonal to `1`. The vector order is the upper triangle, row by row. For
four cell types, the order is:

```text
(CellType_1, CellType_2)
(CellType_1, CellType_3)
(CellType_1, CellType_4)
(CellType_2, CellType_3)
(CellType_2, CellType_4)
(CellType_3, CellType_4)
```

The niche vector uses the same convention. For three niches, its order is
`(Niche_0, Niche_1)`, `(Niche_0, Niche_2)`, and
`(Niche_1, Niche_2)`. Niche IDs are zero-based in the output, whereas the CLI
names cell states `CellType_1`, `CellType_2`, and so on.

Larger affinity values reward local contact between a pair; smaller or
negative values discourage it. The difference among values matters more than
any single value. Start with moderate contrasts, inspect several independent
seeds, and increase the contrast only when the intended organization is not
visible.

## Complete runnable example

Save the following as `designed_parameters.json`. It specifies three niches
and four cell types. Each niche favors a different overlapping cell-type
module, which creates more coherent structure than filling the affinity
vectors with unrelated random values.

```json
{
  "n_group": {
    "value": 3
  },
  "n_state": {
    "value": 4
  },
  "niche_theta": {
    "value": [-0.45, -0.45, -0.45]
  },
  "theta_list": {
    "value": [
      [0.55, -0.45, -0.55, -0.25, -0.45, 0.20],
      [0.15, -0.25, -0.45, 0.45, -0.15, 0.25],
      [0.05, -0.40, -0.50, -0.25, -0.35, 0.55]
    ]
  },
  "density_replicates": {
    "value": [0.34, 0.34, 0.34, 0.34]
  },
  "phi_replicates": {
    "value": 2.8
  }
}
```

Run it with the remaining design choices supplied explicitly:

```bash
simspace reference-free \
  --params designed_parameters.json \
  --output-dir results/designed_reference_free \
  --shape 104 118 \
  --neighbor-distance 2 \
  --cell-sweeps 2 \
  --niche-sweeps 6 \
  --jitter 0.22 \
  --n-genes 100 \
  --bg-ratio 0.7 \
  --lr-ratio 0 \
  --seed 23
```

This produces an irregular MRF-based tissue rather than the curved cortical
bands in Figure 3. That difference is intentional: the standard JSON describes
affinities, not a tissue boundary or a prescribed niche map.

## Build the JSON from full matrices

For more than a few cell types, do not manually count flattened values. Write
the biologically meaningful symmetric matrices, then let NumPy and SimSpace
serialize them in the correct order:

```python
import numpy as np
from simspace import util


def upper_triangle(matrix):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Affinity matrices must be square")
    if not np.allclose(matrix, matrix.T):
        raise ValueError("Affinity matrices must be symmetric")
    return matrix[np.triu_indices(matrix.shape[0], 1)]


niche_matrix = np.array([
    [1.00, -0.45, -0.45],
    [-0.45, 1.00, -0.45],
    [-0.45, -0.45, 1.00],
])

cell_matrices = [
    np.array([
        [1.00, 0.55, -0.45, -0.55],
        [0.55, 1.00, -0.25, -0.45],
        [-0.45, -0.25, 1.00, 0.20],
        [-0.55, -0.45, 0.20, 1.00],
    ]),
    np.array([
        [1.00, 0.15, -0.25, -0.45],
        [0.15, 1.00, 0.45, -0.15],
        [-0.25, 0.45, 1.00, 0.25],
        [-0.45, -0.15, 0.25, 1.00],
    ]),
    np.array([
        [1.00, 0.05, -0.40, -0.50],
        [0.05, 1.00, -0.25, -0.35],
        [-0.40, -0.25, 1.00, 0.55],
        [-0.50, -0.35, 0.55, 1.00],
    ]),
]

params = {
    "n_group": 3,
    "n_state": 4,
    "niche_theta": upper_triangle(niche_matrix),
    "theta_list": [upper_triangle(matrix) for matrix in cell_matrices],
    "density_replicates": np.repeat(0.34, 4),
    "phi_replicates": 2.8,
}

util.save_params(params, "designed_parameters.json")
```

The diagonal entries shown in the full matrices are reminders; SimSpace fixes
them to `1` when reconstructing the matrices.

## How each choice changes the result

### Niche affinities

- Keep off-diagonal `niche_theta` values well below the fixed diagonal when
  you want recognizable, contiguous domains.
- Bring selected off-diagonal values closer to the diagonal when two niches
  should intermix or share a broad interface.
- Avoid extreme values at first. They can produce nearly isolated islands and
  make seed-to-seed behavior harder to assess.

### Cell-type affinities

- Use positive within-module and negative between-module affinities to express
  a simple local organization.
- Change `theta_list[m]` only when niche `m` should have a distinct local
  pattern.
- Affinities control pairwise proximity, not a directional cell-cell
  interaction and not a causal signaling mechanism.
- The standard JSON has no niche-by-cell-type abundance table. If exact
  niche-specific mixtures are important, initialize the cell grid by those
  mixtures in Python as Figure 3 does.

### Density

`density_replicates[k]` is the probability of retaining a candidate location
whose assigned cell type is `k`. With uniform retention `d`, the expected cell
count is approximately `ROWS × COLS × d` before considering a custom tissue
mask. Unequal retention also changes the displayed cell-type composition, so
use it deliberately rather than as a cosmetic control.

### Smoothness, neighborhoods, and sweeps

- Increasing `phi_replicates` generally strengthens the spatial organization
  encoded by the affinities.
- Increasing `--neighbor-distance` broadens the local cell-type context.
- `--cell-sweeps` and `--niche-sweeps` are finite generative controls, not
  convergence settings. More sweeps allow more reorganization from the
  initialized labels.
- `--jitter` changes coordinate appearance after labels are generated; it does
  not create biological spatial dependence.

These controls interact. Tune them with fixed seeds first, then repeat the
selected design across new seeds to ensure that its qualitative behavior is
not a single-realization accident.

## Reproducing the Figure 3 design strategy

The cortex-inspired Figure 3 anchor separates four design layers:

1. **Geometry:** a `104 × 118` lattice is restricted to a curved sector. An
   inner-radius parameter controls cortical thickness.
2. **Niches:** normalized cortical depth is cut into four curved bands, with
   small smooth perturbations at the interfaces.
3. **Cell types:** nine declared cortical excitatory subclass labels are
   initialized from a different categorical mixture in each band, followed by
   one niche-conditioned cell-type Gibbs sweep.
4. **Molecular profiles:** named markers receive declared
   cell-type-conditional expression programs. Their spatial patterns follow
   the simulated cell labels; no direct coordinate-conditioned expression
   effect is used for the displayed genes.

The selected anchor uses neighbor distance `2`, `phi = 2.6`, uniform retention
`0.34`, coordinate jitter `0.22`, and fixed geometry, initialization, and
simulation seeds. Figure 3B then varies only cortical thickness and
layer-specificity while holding all other settings fixed. This one- or
two-parameter-at-a-time sweep is a useful template for designing any new
reference-free simulation family.

The exact implementation and archived parameters are in the manuscript
reproducibility repository:

- [`cortex_reference_free.py`](https://github.com/TianxiaoNYU/simspace-reproducibility/blob/main/main_figures/Fig3/cortex_reference_free.py)
  builds the geometry, niche mixtures, cell-type affinities, and anchor;
- [`cortex_parameters.json`](https://github.com/TianxiaoNYU/simspace-reproducibility/blob/main/main_figures/Fig3/output/cortex_parameters.json)
  records the complete Figure 3 provenance manifest; and
- [`Fig3.py`](https://github.com/TianxiaoNYU/simspace-reproducibility/blob/main/main_figures/Fig3/Fig3.py)
  generates the controlled parameter sweep and molecular profiles.

Use the script as a pattern when the desired tissue needs a mask, curved or
layered boundaries, named cell types, explicit niche-specific mixtures, or a
custom marker program. Use the smaller standard JSON when an irregular
MRF-generated domain layout is sufficient.

## A practical tuning workflow

1. Generate one scaffold file with the CLI or the matrix-building code above.
2. Fix the random seed and vary one design family at a time: geometry, niche
   affinity, cell-type affinity, MRF scale, density, or jitter.
3. Compare a small grid of candidates with identical plotting limits, point
   sizes, colors, and seeds.
4. Select parameters using declared checks relevant to the intended use, such
   as cell count, cell-type proportions, domain sizes, boundary continuity, or
   median position along a designed axis.
5. Run the selected parameters across several new seeds and report the full
   range rather than retaining only the most attractive realization.
6. Save the selected `parameters.json`, all non-JSON CLI arguments, the seed,
   the SimSpace version, and any custom geometry or expression code.

## Validation checklist

Before running a large simulation, confirm that:

- all six required fields are present and every field has a `value` wrapper;
- `theta_list` contains exactly `n_group` vectors;
- `niche_theta` has `n_group × (n_group - 1) / 2` values;
- every cell-type affinity vector has
  `n_state × (n_state - 1) / 2` values;
- `density_replicates` has `n_state` values, all between `0` and `1`;
- JSON contains no comments, `NaN`, or trailing commas;
- the shape, neighborhood, sweeps, jitter, molecular settings, and seed are
  recorded beside the JSON; and
- conclusions are based on multiple seeds, not appearance alone.


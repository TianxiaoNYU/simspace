import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing
import subprocess
import os
from scipy.spatial import KDTree


def buildOmicsMean(
        omics_meta: pd.DataFrame,
        meta: pd.DataFrame,
        ) -> pd.DataFrame:
    """Build the legacy phenotype-conditioned Poisson mean matrix.

    This helper is deterministic and does not consume random numbers.  It is
    intentionally separate from :func:`simOmics` so the legacy simulation path
    and its random-number ordering remain unchanged.
    """
    if 'state' not in meta.columns:
        raise ValueError('The metadata does not contain the "state" column. Please check the column names.')
    if not isinstance(omics_meta, pd.DataFrame) or not isinstance(meta, pd.DataFrame):
        raise TypeError('The input data must be a pandas DataFrame.')

    gene_columns = [f'Gene_{gene_id}' for gene_id in omics_meta['GeneID']]
    means = np.empty((len(meta), len(gene_columns)), dtype=float)
    states = meta['state'].to_numpy()

    for state in pd.unique(states):
        type_column = f'Type_{state}'
        if type_column not in omics_meta.columns:
            raise ValueError(f'Omics metadata does not contain the required column "{type_column}".')
        means[states == state, :] = omics_meta[type_column].to_numpy(dtype=float)

    return pd.DataFrame(means, columns=gene_columns, index=meta.index)


def _scale_spatial_coordinates(coordinates: np.ndarray, scaling: str) -> np.ndarray:
    """Scale two-dimensional coordinates for an optional spatial basis."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError('Spatial coordinates must be an n_cells by 2 array.')
    if not np.isfinite(coordinates).all():
        raise ValueError('Spatial coordinates must contain only finite values.')

    if scaling == 'none':
        return coordinates.copy()
    if scaling != 'unit_range':
        raise ValueError('coordinate_scaling must be either "unit_range" or "none".')

    lower = coordinates.min(axis=0)
    span = coordinates.max(axis=0) - lower
    span[span == 0] = 1.0
    return (coordinates - lower) / span


def _build_spatial_basis(meta: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, np.ndarray]:
    """Construct a reproducible, named spatial basis from a configuration."""
    coordinate_columns = tuple(config.get('coordinate_columns', ('row', 'col')))
    if len(coordinate_columns) != 2:
        raise ValueError('coordinate_columns must contain exactly two column names.')
    missing = [column for column in coordinate_columns if column not in meta.columns]
    if missing:
        raise ValueError(f'Cell metadata is missing spatial coordinate columns: {missing}.')

    scaling = config.get('coordinate_scaling', 'unit_range')
    coordinates = _scale_spatial_coordinates(
        meta.loc[:, list(coordinate_columns)].to_numpy(dtype=float),
        scaling,
    )
    basis_name = config.get('basis', 'linear')

    if basis_name == 'linear':
        basis = coordinates.copy()
        basis_columns = [f'{coordinate_columns[0]}_linear', f'{coordinate_columns[1]}_linear']
    elif basis_name in {'radial', 'hotspot'}:
        default_center = np.array([0.5, 0.5]) if scaling == 'unit_range' else coordinates.mean(axis=0)
        center = np.asarray(config.get('center', default_center), dtype=float)
        if center.shape != (2,) or not np.isfinite(center).all():
            raise ValueError('center must contain two finite coordinates.')
        distance = np.linalg.norm(coordinates - center[None, :], axis=1)
        if basis_name == 'radial':
            basis = distance[:, None]
            basis_columns = ['radial_distance']
        else:
            bandwidth = float(config.get('bandwidth', 0.2 if scaling == 'unit_range' else 1.0))
            if not np.isfinite(bandwidth) or bandwidth <= 0:
                raise ValueError('bandwidth must be a positive finite number.')
            basis = np.exp(-0.5 * (distance / bandwidth) ** 2)[:, None]
            basis_columns = ['hotspot_basis']
    elif basis_name == 'structure_distance':
        structure_points = np.asarray(config.get('structure_points'), dtype=float)
        if structure_points.ndim != 2 or structure_points.shape[1] != 2 or len(structure_points) == 0:
            raise ValueError('structure_points must be a non-empty n_points by 2 array.')
        if not np.isfinite(structure_points).all():
            raise ValueError('structure_points must contain only finite values.')
        distance = np.linalg.norm(
            coordinates[:, None, :] - structure_points[None, :, :],
            axis=2,
        ).min(axis=1)
        basis = distance[:, None]
        basis_columns = ['structure_distance']
    else:
        raise ValueError(
            'basis must be one of "linear", "radial", "hotspot", or "structure_distance".'
        )

    if bool(config.get('center_basis', True)):
        basis = basis - basis.mean(axis=0, keepdims=True)

    basis_frame = pd.DataFrame(basis, columns=basis_columns, index=meta.index)
    coordinate_frame = pd.DataFrame(
        coordinates,
        columns=[f'{coordinate_columns[0]}_scaled', f'{coordinate_columns[1]}_scaled'],
        index=meta.index,
    )
    return pd.concat([coordinate_frame, basis_frame], axis=1), basis


def _resolve_spatial_genes(omics_meta: pd.DataFrame, genes) -> list[int]:
    """Resolve configured gene IDs or generated names to unique row positions."""
    if genes is None:
        raise ValueError('direct_spatial_effects must declare a non-empty "genes" collection.')
    if isinstance(genes, (str, int, np.integer)):
        genes = [genes]
    genes = list(genes)
    if len(genes) == 0:
        raise ValueError('direct_spatial_effects must declare at least one gene.')

    gene_ids = list(omics_meta['GeneID'])
    gene_names = [f'Gene_{gene_id}' for gene_id in gene_ids]
    resolved = []
    for gene in genes:
        if gene in gene_ids:
            position = gene_ids.index(gene)
        elif isinstance(gene, str) and gene in gene_names:
            position = gene_names.index(gene)
        else:
            raise ValueError(f'Unknown spatial-effect gene: {gene!r}.')
        if position in resolved:
            raise ValueError(f'Spatial-effect gene {gene!r} was specified more than once.')
        resolved.append(position)
    return resolved


def _coefficient_vector(coefficients, gene_name: str, gene_id, n_basis: int) -> np.ndarray:
    """Return one coefficient vector from a shared vector or gene-keyed mapping."""
    coefficient_value = coefficients
    if isinstance(coefficients, dict):
        for key in (gene_name, gene_id, str(gene_id)):
            if key in coefficients:
                coefficient_value = coefficients[key]
                break
        else:
            raise ValueError(f'No coefficient vector was provided for {gene_name}.')

    coefficient_vector = np.asarray(coefficient_value, dtype=float)
    if coefficient_vector.ndim == 0:
        coefficient_vector = coefficient_vector.reshape(1)
    if coefficient_vector.shape != (n_basis,):
        raise ValueError(
            f'Coefficients for {gene_name} must contain {n_basis} value(s); '
            f'received shape {coefficient_vector.shape}.'
        )
    if not np.isfinite(coefficient_vector).all():
        raise ValueError(f'Coefficients for {gene_name} must be finite.')
    return coefficient_vector


def simOmicsWithSpatialEffects(
        omics_meta: pd.DataFrame,
        meta: pd.DataFrame,
        direct_spatial_effects: dict,
        seed: int = 0,
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Simulate counts with opt-in direct spatial effects on Poisson means.

    The existing phenotype-conditioned means are multiplied by
    ``exp(B @ gamma_g)`` only for explicitly selected genes.  A zero
    coefficient vector therefore recovers the legacy mean and, for the same
    seed, the legacy count realization.
    """
    if not isinstance(direct_spatial_effects, dict):
        raise TypeError('direct_spatial_effects must be a dictionary or None.')
    if 'coefficients' not in direct_spatial_effects:
        raise ValueError('direct_spatial_effects must provide "coefficients".')

    latent_mean = buildOmicsMean(omics_meta, meta)
    design_frame, basis = _build_spatial_basis(meta, direct_spatial_effects)
    selected_positions = _resolve_spatial_genes(
        omics_meta,
        direct_spatial_effects.get('genes'),
    )

    truth_rows = []
    for position in selected_positions:
        gene_id = omics_meta['GeneID'].iloc[position]
        gene_name = f'Gene_{gene_id}'
        coefficients = _coefficient_vector(
            direct_spatial_effects['coefficients'],
            gene_name,
            gene_id,
            basis.shape[1],
        )
        linear_predictor = basis @ coefficients
        if not np.isfinite(linear_predictor).all():
            raise ValueError(f'The spatial linear predictor for {gene_name} is not finite.')
        if np.max(np.abs(linear_predictor)) > 50:
            raise ValueError(
                f'The spatial linear predictor for {gene_name} exceeds 50 in absolute value; '
                'reduce the configured coefficient magnitude.'
            )
        multiplier = np.exp(linear_predictor)
        latent_mean.iloc[:, position] *= multiplier

        truth_row = {
            'GeneID': gene_id,
            'Gene': gene_name,
            'spatial_basis': direct_spatial_effects.get('basis', 'linear'),
            'min_log_effect': float(linear_predictor.min()),
            'max_log_effect': float(linear_predictor.max()),
            'mean_multiplier': float(multiplier.mean()),
        }
        for coefficient_index, coefficient in enumerate(coefficients):
            truth_row[f'coefficient_{coefficient_index}'] = float(coefficient)
        truth_rows.append(truth_row)

    np.random.seed(seed)
    omics_data = np.zeros(latent_mean.shape, dtype=float)
    mean_values = latent_mean.to_numpy(dtype=float)
    for cell_index in range(mean_values.shape[0]):
        for gene_index in range(mean_values.shape[1]):
            omics_data[cell_index, gene_index] = np.random.poisson(
                mean_values[cell_index, gene_index]
            )

    simulated_counts = pd.DataFrame(
        omics_data,
        columns=latent_mean.columns,
        index=meta.index,
    )
    gene_truth = pd.DataFrame(truth_rows)
    return simulated_counts, latent_mean, gene_truth, design_frame


def _expand_cell_parameter(value, index: pd.Index, name: str) -> np.ndarray:
    """Expand a scalar or cell-indexed parameter to a finite numeric vector."""
    if np.isscalar(value):
        values = np.full(len(index), float(value), dtype=float)
    elif isinstance(value, pd.Series):
        values = value.reindex(index).to_numpy(dtype=float)
    elif isinstance(value, dict):
        values = pd.Series(value).reindex(index).to_numpy(dtype=float)
    else:
        values = np.asarray(value, dtype=float)
        if values.shape != (len(index),):
            raise ValueError(f'{name} must be a scalar or contain one value per cell.')
    if not np.isfinite(values).all():
        raise ValueError(f'{name} must contain only finite values and cover every cell.')
    return values


def applyCaptureNoise(
        gene_data: pd.DataFrame,
        capture_efficiency=1.0,
        seed: int = 0,
        ) -> tuple[pd.DataFrame, pd.Series]:
    """Apply optional binomial capture thinning to an integer count matrix."""
    efficiencies = _expand_cell_parameter(
        capture_efficiency,
        gene_data.index,
        'capture_efficiency',
    )
    if np.any((efficiencies < 0) | (efficiencies > 1)):
        raise ValueError('capture_efficiency must lie between 0 and 1.')
    efficiency_series = pd.Series(
        efficiencies,
        index=gene_data.index,
        name='capture_efficiency',
    )
    if np.all(efficiencies == 1):
        return gene_data.copy(), efficiency_series

    counts = gene_data.to_numpy(dtype=float)
    if np.any(counts < 0) or not np.allclose(counts, np.rint(counts)):
        raise ValueError('Capture thinning requires a non-negative integer count matrix.')
    rng = np.random.default_rng(seed)
    thinned = rng.binomial(np.rint(counts).astype(np.int64), efficiencies[:, None])
    return pd.DataFrame(thinned, index=gene_data.index, columns=gene_data.columns), efficiency_series


def applyAmbientNoise(
        gene_data: pd.DataFrame,
        ambient_rate=0.0,
        ambient_profile=None,
        seed: int = 0,
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    """Add optional Poisson background contamination to simulated profiles."""
    rates = _expand_cell_parameter(ambient_rate, gene_data.index, 'ambient_rate')
    if np.any(rates < 0):
        raise ValueError('ambient_rate must be non-negative.')

    if ambient_profile is None:
        profile_values = gene_data.sum(axis=0).to_numpy(dtype=float) + 1.0
    elif isinstance(ambient_profile, pd.Series):
        profile_values = ambient_profile.reindex(gene_data.columns).to_numpy(dtype=float)
    elif isinstance(ambient_profile, dict):
        profile_values = pd.Series(ambient_profile).reindex(gene_data.columns).to_numpy(dtype=float)
    else:
        profile_values = np.asarray(ambient_profile, dtype=float)
        if profile_values.shape != (gene_data.shape[1],):
            raise ValueError('ambient_profile must contain one value per gene.')
    if not np.isfinite(profile_values).all() or np.any(profile_values < 0):
        raise ValueError('ambient_profile must contain finite, non-negative values for every gene.')
    if profile_values.sum() <= 0:
        raise ValueError('ambient_profile must have a positive sum.')
    profile_values = profile_values / profile_values.sum()

    rate_series = pd.Series(rates, index=gene_data.index, name='ambient_rate')
    profile_series = pd.Series(
        profile_values,
        index=gene_data.columns,
        name='ambient_profile',
    )
    if np.all(rates == 0):
        additions = pd.DataFrame(
            np.zeros(gene_data.shape, dtype=np.int64),
            index=gene_data.index,
            columns=gene_data.columns,
        )
        return gene_data.copy(), additions, rate_series, profile_series

    rng = np.random.default_rng(seed)
    additions_array = rng.poisson(rates[:, None] * profile_values[None, :])
    additions = pd.DataFrame(
        additions_array,
        index=gene_data.index,
        columns=gene_data.columns,
    )
    observed = gene_data.astype(float) + additions
    return observed, additions, rate_series, profile_series


def applyDropout(
        gene_data: pd.DataFrame,
        dropout: dict | None,
        latent_mean: pd.DataFrame | None = None,
        seed: int = 0,
        ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply optional constant or mean-dependent excess dropout."""
    if dropout is None:
        probabilities = np.zeros(gene_data.shape, dtype=float)
    elif not isinstance(dropout, dict):
        raise TypeError('dropout must be a dictionary or None.')
    else:
        mode = dropout.get('mode', 'constant')
        if mode == 'constant':
            probability = float(dropout.get('probability', 0.0))
            if not np.isfinite(probability) or probability < 0 or probability > 1:
                raise ValueError('Constant dropout probability must lie between 0 and 1.')
            probabilities = np.full(gene_data.shape, probability, dtype=float)
        elif mode == 'mean_dependent':
            if latent_mean is None:
                raise ValueError('mean_dependent dropout requires the latent mean matrix.')
            latent_mean = latent_mean.reindex(index=gene_data.index, columns=gene_data.columns)
            if latent_mean.isna().any().any():
                raise ValueError('The latent mean matrix does not align with the observed count matrix.')
            intercept = float(dropout.get('intercept', 0.0))
            slope = float(dropout.get('slope', -1.0))
            if not np.isfinite(intercept) or not np.isfinite(slope):
                raise ValueError('Dropout intercept and slope must be finite.')
            linear_predictor = intercept + slope * np.log1p(latent_mean.to_numpy(dtype=float))
            linear_predictor = np.clip(linear_predictor, -30, 30)
            probabilities = 1.0 / (1.0 + np.exp(-linear_predictor))
        else:
            raise ValueError('dropout mode must be either "constant" or "mean_dependent".')

    probability_frame = pd.DataFrame(
        probabilities,
        index=gene_data.index,
        columns=gene_data.columns,
    )
    if np.all(probabilities == 0):
        mask = np.zeros(gene_data.shape, dtype=bool)
        return gene_data.copy(), pd.DataFrame(mask, index=gene_data.index, columns=gene_data.columns), probability_frame

    rng = np.random.default_rng(seed)
    mask = rng.random(gene_data.shape) < probabilities
    # Construct the result from an explicitly writable NumPy copy.  With
    # pandas copy-on-write enabled, ``DataFrame.values`` can be read-only even
    # after ``DataFrame.copy()``, so in-place assignment through that view is
    # not portable across supported pandas versions.
    observed_values = gene_data.to_numpy(copy=True)
    observed_values[mask] = 0
    observed = pd.DataFrame(
        observed_values,
        index=gene_data.index,
        columns=gene_data.columns,
    )
    mask_frame = pd.DataFrame(mask, index=gene_data.index, columns=gene_data.columns)
    return observed, mask_frame, probability_frame

def simOmicsMeta(
        meta: pd.DataFrame,
        n_genes: int = 1000,
        bg_ratio: float = 0.5,
        bg_param: typing.Tuple[float, float] = (1, 1),
        marker_param: typing.Tuple[float, float] = (5, 2),
        lr_ratio: float = 0.5,
        seed: int = 0,
        ) -> pd.DataFrame:
    """
    Simulate the metadata of omics data
    Args:
        meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        n_genes (int): The number of genes to simulate.
        bg_ratio (float): The ratio of background genes (non-marker genes) to the total number of genes.
        bg_param (tuple): The parameters for the gamma distribution to simulate the mean expression level of background genes.
        marker_param (tuple): The parameters for the gamma distribution to  simulate the mean expression level of marker genes.
        lr_ratio (float): The ratio of ligand-receptor pairs to the total number of marker genes.
        seed (int): The random seed for reproducibility.
    Returns:
        pd.DataFrame: A DataFrame containing the simulated metadata of the omics data, with columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
            - Type_{cell_type}: The mean expression level of the gene in the corresponding cell type.
    
    Raises:
        ValueError: If the background gene ratio or ligand-receptor pair ratio is not between 0 and 1, or if the metadata does not contain the 'state' column.
    
    Examples:
        >>> meta = pd.DataFrame({'state': ['A', 'B', 'C']})
        >>> omics_meta = simOmicsMeta(meta, n_genes=100, bg_ratio=0.3, lr_ratio=0.2, seed=42)
        >>> print(omics_meta.head())
    """

    np.random.seed(seed)

    if bg_ratio < 0 or bg_ratio > 1:
        raise ValueError('The background gene ratio should be between 0 and 1.')
    if lr_ratio < 0 or lr_ratio > 1:
        raise ValueError('The ligand-receptor pair ratio should be between 0 and 1.')
    if 'state' not in meta.columns:
        raise ValueError('The metadata does not contain the "state" column. Please check the column names.')
    
    type_arr = meta['state'].unique()
    type_arr = sorted(type_arr)
    n_types = len(type_arr)

    # Randomly assign a cell type to each gene
    marker_arr = np.random.choice(type_arr, n_genes)

    # Randomly set a part of genes as background genes (non-marker genes)
    n_bg = int(n_genes * bg_ratio)
    bg_indices = np.random.choice(n_genes, n_bg, replace=False)
    marker_arr[bg_indices] = -1

    # Randomly generates the ligand-receptor pair among genes that are not background genes
    ligand_receptor_arr = np.full(n_genes, -1)
    marker_indices = np.where(marker_arr != -1)[0]
    n_pairs = int(len(marker_indices) * lr_ratio)
    paired_indices = np.random.choice(marker_indices, n_pairs, replace=False)

    for i in paired_indices:
        if ligand_receptor_arr[i] != -1:
            continue
        possible_pairs = np.where((marker_arr != -1) & (ligand_receptor_arr == -1))[0]
        possible_pairs = possible_pairs[possible_pairs != i]
        if len(possible_pairs) > 0:
            tmp = np.random.choice(possible_pairs)
            ligand_receptor_arr[i] = tmp
            ligand_receptor_arr[tmp] = i

    # Create the omics metadata
    omics_meta = pd.DataFrame({
        "GeneID": np.arange(n_genes),
        "Marker": marker_arr,
        "LRindex": ligand_receptor_arr,
        })
    
    # Simulate the mean expression level of each gene
    for celltype in type_arr:
        omics_meta[f'Type_{celltype}'] = np.zeros(n_genes)

    for i in range(n_genes):
        gene_lam = np.random.gamma(bg_param[0], bg_param[1], size=n_types)
        omics_meta.iloc[i, 3:] = gene_lam

        if omics_meta.loc[i, 'Marker'] != -1:
            omics_meta.loc[i, f'Type_{omics_meta.loc[i, "Marker"]}'] = np.random.gamma(marker_param[0], marker_param[1], 1)

    return omics_meta

def simOmics(omics_meta: pd.DataFrame, 
             meta: pd.DataFrame,
             seed: int = 0,
             ) -> pd.DataFrame:
    """
    Simulate the omics data based on the metadata and cell metadata. 

    Args:
        omics_meta (pd.DataFrame): The metadata of the omics data, which should contain the following columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
            - Type_{cell_type}: The mean expression level of the gene in the corresponding cell type.
        meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        seed (int): The random seed for reproducibility.

    Returns:
        pd.DataFrame: A DataFrame containing the simulated omics data, with columns:
            - Gene_{gene_id}: The expression level of the gene in each cell.

    Raises:
        ValueError: If the metadata does not contain the 'state' column.
        TypeError: If the input data is not a DataFrame.

    Examples:
        >>> omics_meta = pd.DataFrame({
        ...     'GeneID': [0, 1, 2],
        ...     'Marker': ['A', 'B', -1],
        ...     'LRindex': [-1, 0, -1],
        ...     'Type_A': [10, 20, 0],
        ...     'Type_B': [0, 30, 0]
        ... })
        >>> meta = pd.DataFrame({'state': ['A', 'B']})
        >>> omics_data = simOmics(omics_meta, meta, seed=42)
        >>> print(omics_data.head())
    """
    if 'state' not in meta.columns:
        raise ValueError('The metadata does not contain the "state" column. Please check the column names.')
    if not isinstance(omics_meta, pd.DataFrame) or not isinstance(meta, pd.DataFrame):
        raise TypeError('The input data must be a pandas DataFrame.')
    
    np.random.seed(seed)
    
    n_genes = omics_meta.shape[0]
    n_cells = len(meta)

    # Simulate the expression level of each gene in each cell
    # The expression level is generated from a poisson distribution
    omics_data = np.zeros((n_cells, n_genes))
    for i in range(n_cells):
        cell_type = meta.loc[i, 'state']
        for j in range(n_genes):
            expr = np.random.poisson(omics_meta.loc[j, f'Type_{cell_type}'])
            omics_data[i, j] = expr

    omics_data = pd.DataFrame(omics_data, columns=[f'Gene_{gene_id}' for gene_id in omics_meta['GeneID']], index=meta.index)

    return omics_data

def simSpatialOmics(gene_data: pd.DataFrame,
                    gene_meta: pd.DataFrame,
                    cell_meta: pd.DataFrame,
                    k_neighors: int = 10,
                    spatial_effect: float = 1.0,
                    se_threshold: float = 1.5,
                    seed: int = 0,
                    ) -> pd.DataFrame:
    """
    Simulate the spatial omics data

    Args:
        gene_data (pd.DataFrame): The gene expression data, with cells as rows and genes as columns.
        gene_meta (pd.DataFrame): The metadata of the genes, which should contain the following columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
        cell_meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        k_neighors (int): The number of nearest neighbors to consider for spatial effects.
        spatial_effect (float): The factor by which to increase or decrease the expression level based on spatial effects.
        se_threshold (float): The threshold for spatial effect application.
        seed (int): The random seed for reproducibility.

    Returns:
        pd.DataFrame: A DataFrame containing the simulated spatial omics data, with cells as rows and genes as columns.

    Raises:
        ValueError: If the spatial effect is not greater than 1, or if the cell metadata does not contain the coordinates or cell types.
    
    Examples:
        >>> gene_data = pd.DataFrame({
        ...     'Gene_0': [10, 20, 30],
        ...     'Gene_1': [5, 15, 25],
        ...     'Gene_2': [0, 10, 20]
        ... }, index=['cell_1', 'cell_2', 'cell_3'])
        >>> gene_meta = pd.DataFrame({
        ...     'GeneID': [0, 1, 2],
        ...     'Marker': ['A', 'B', -1],
        ...     'LRindex': [-1, 0, -1]
        ... })
        >>> cell_meta = pd.DataFrame({
        ...     'state': ['A', 'B', 'A'],
        ...     'row': [0, 1, 0],
        ...     'col': [0, 1, 2]
        ... })
        >>> spatial_omics = simSpatialOmics(gene_data, gene_meta, cell_meta, k_neighors=2, spatial_effect=2, seed=42)
        >>> print(spatial_omics.head())
    """

    assert spatial_effect > 0, 'The spatial effect should be greater than 0.'

    if 'col' not in cell_meta.columns or 'row' not in cell_meta.columns:
        raise ValueError('The cell metadata does not contain the coordinates. Please check the column names.')
    
    if 'state' not in cell_meta.columns:
        raise ValueError('The cell metadata does not contain the cell types. Please check the column names.')
    
    output_data = gene_data.copy()

    np.random.seed(seed)
    n_gene = gene_data.shape[1]
    n_cell = gene_data.shape[0]

    # Extract coordinates
    coordinates = cell_meta[['row', 'col']].values

    # Build a KDTree for efficient neighbor search
    tree = KDTree(coordinates)

    # Find the top 10 nearest neighbors for each cell
    _, neighbors = tree.query(coordinates, k=k_neighors+1)  # k+1 because the closest neighbor is the point itself

    # Create a dictionary to store neighbors for each cell, excluding the cell itself
    cell_neighbors = {k: neighbors[k][1:] for k in range(n_cell)}

    for i in range(n_gene):
        LRindex = gene_meta.loc[i, 'LRindex']
        Marker = gene_meta.loc[i, 'Marker']
        if Marker == -1:
            continue
        if Marker == LRindex:
            continue
        if LRindex != -1:
            mean_expr = gene_data.iloc[:, LRindex].mean()
            for j in range(n_cell):
                if cell_meta.state[j] == gene_meta.loc[i, 'Marker']:
                    neighbor_indices = cell_neighbors[j]
                    neighbor_expr = gene_data.iloc[neighbor_indices, LRindex].values.sum()
                    neighbor_ratio = neighbor_expr / mean_expr / k_neighors
                    if neighbor_ratio > se_threshold:
                        output_data.iloc[j, i] *= spatial_effect * neighbor_ratio
                    else:
                        output_data.iloc[j, i] /= spatial_effect
                
    return output_data

def run_splatter(new_meta, 
                 ngene = 1000, 
                 r_script_path=None):
    """
    Run the splatter simulation to generate synthetic single-cell RNA-seq data.

    Args:
        new_meta (pd.DataFrame): A DataFrame containing simulated spatial metadata for omics simulation, which is derived from the .meta of the simspace object.
        ngene (int): The number of genes to simulate.
        r_script_path (str): The path to the R script that performs the splatter simulation. Default is None, which uses the script of simspace package.
    
    Returns:
        tuple: A tuple containing two DataFrames:
            - splatter_data: The simulated gene expression data.
            - splatter_meta: The metadata of the simulated cells.
    
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    
    Examples:
        >>> splatter_data, splatter_meta = run_splatter(sim.meta, ngene=1000) # sim is simulated simspace object
        >>> print(splatter_data.head())
        >>> print(splatter_meta.head())
    """

    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/splatter.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")
    new_meta.to_csv("./tmp/cells_meta.csv", index=False)
        
    result = subprocess.run(
        ["Rscript", r_script_path, "./tmp/cells_meta.csv", ngene],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("Splatter simulation complete.")
        # Load the simulated data
        splatter_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        splatter_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/cells_meta.csv")
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")

        return splatter_data, splatter_meta

def splatter_fit(count_path, 
                 meta_path, 
                 group_col, 
                 n_cells = 2000,
                 r_script_path=None):
    """
    Fit the splatter model to the given reference data and metadata to simulate new omics data.

    Args:
        count_path (str): The path to the count matrix file.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information.
        n_cells (int): The number of cells to simulate. Should match the number of cells in the spatial simulation results.
        r_script_path (str): The path to the R script that performs the splatter fitting. Default is None, which uses the script of simspace package.
    
    Returns:
        tuple: A tuple containing two DataFrames:
            - splatter_data: The simulated gene expression data.
            - splatter_meta: The metadata of the simulated cells.
    
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    
    Examples:
        >>> splatter_data, splatter_meta = splatter_fit("path/to/count.csv",
        ...                                              "path/to/meta.csv",
        ...                                              "group_column",
        ...                                              n_cells=2000)
        >>> print(splatter_data.head())
        >>> print(splatter_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/splatter_fit.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")

    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    result = subprocess.run(
        ["Rscript", r_script_path, meta_path, count_path, group_col, str(n_cells)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("Splatter fit complete.")
        splatter_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        splatter_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return splatter_data, splatter_meta

def scdesign_fit(count_path, 
                 meta_path, 
                 group_col, 
                 spatial_x,
                 spatial_y,
                 new_meta,
                 seed = 0,
                 r_script_path=None,
                 family="auto"):
    """
    Fit the scdesign model to the given reference data and metadata to simulate new omics data.
    Args:
        count_path (str): Path to a nonnegative reference feature matrix. Integer-valued
            data use negative-binomial marginals and continuous data use log-Gaussian
            marginals when ``family='auto'``.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information.
        spatial_x (str): The column name in the metadata that contains the x-coordinate of the spatial location.
        spatial_y (str): The column name in the metadata that contains the y-coordinate of the spatial location.
        new_meta (pd.DataFrame): A DataFrame containing simulated spatial metadata for omics simulation, which is derived from the .meta of the simspace object.
        seed (int): The random seed for reproducibility.
        r_script_path (str): The path to the R script that performs the scDesign fitting. Default is None, which uses the script of simspace package.
        family (str): Marginal family selection: ``'auto'``, ``'nb'``, or
            ``'gaussian'``. The Gaussian pathway models ``log1p`` intensities and
            returns simulated values on the original scale.
    Returns:
        tuple: A tuple containing two DataFrames:
            - sim_data: The simulated gene expression data.
            - sim_meta: The metadata of the simulated cells.
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    Examples:
        >>> sim_data, sim_meta = scdesign_fit("path/to/count.csv",
        ...                              "path/to/meta.csv",
        ...                              "group_column",
        ...                              "x_coordinate",
        ...                              "y_coordinate",
        ...                              new_meta=sim.meta, # sim is simulated simspace object
        ...                              seed=42)
        >>> print(sim_data.head())
        >>> print(sim_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/scdesign.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    family = str(family).lower()
    if family not in {"auto", "nb", "gaussian"}:
        raise ValueError("family must be one of 'auto', 'nb', or 'gaussian'")
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    new_meta.to_csv("./tmp/new_meta.csv", index=False)

    command = [
        "Rscript", r_script_path, meta_path, count_path, group_col,
        spatial_x, spatial_y, str(seed),
    ]
    # Preserve compatibility with custom six-argument scripts when auto-detection
    # is requested. The bundled script defaults a missing family argument to auto.
    if family != "auto":
        command.append(family)

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stdout)
        print(result.stderr)
        return None, None
    else:
        for line in result.stdout.splitlines():
            if line.startswith("Using scDesign3 family"):
                print(line)
        print("scDesign3 fit complete.")
        sim_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        sim_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        os.remove("./tmp/new_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return sim_data, sim_meta
    

def srtsim_fit(count_path, 
               meta_path, 
               group_col = 'state', 
               spatial_x = 'x',
               spatial_y = 'y',
               n_rep = 1,
               seed = 0,
               r_script_path=None):
    """
    Fit the SRTsim model to the given reference data and metadata to simulate new omics data.
    Args:
        count_path (str): The path to the count matrix file.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information. Default is 'state'.
        spatial_x (str): The column name in the metadata that contains the x-coordinate of the spatial location.
        spatial_y (str): The column name in the metadata that contains the y-coordinate of the spatial location.
        n_rep (int): The number of replicates to simulate. Default is 1. Since SRTsim can only simulate exact same number of cells as the reference, this parameter is used when the number of cells in the reference is less than the number of cells in the spatial simulation results.
        seed (int): The random seed for reproducibility. Default is 0.
        r_script_path (str): The path to the R script that performs the SRTsim fitting. Default is None, which uses the script of simspace package.
    Returns:
        tuple: A tuple containing two DataFrames:
            - sim_data: The simulated gene expression data.
            - sim_meta: The metadata of the simulated cells.
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    Examples:
        >>> sim_data, sim_meta = srtsim_fit("path/to/count.csv",
        ...                              "path/to/meta.csv",
        ...                              group_col='state',
        ...                              spatial_x='x',
        ...                              spatial_y='y',
        ...                              n_rep=1,
        ...                              seed=42)
        >>> print(sim_data.head())
        >>> print(sim_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/srtsim.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    result = subprocess.run(
        ["Rscript", r_script_path, meta_path, count_path, group_col, spatial_x, spatial_y, str(n_rep), str(seed)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("SRTsim fit complete.")
        sim_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        sim_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return sim_data, sim_meta

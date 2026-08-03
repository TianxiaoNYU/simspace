import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from simspace import SimSpace, omics
from simspace.util import sim_from_json
from simspace.spatial import generate_offsets


def _small_omics_inputs():
    meta = pd.DataFrame({
        'state': [0, 0, 1, 1],
        'row': [0.0, 1.0, 0.0, 1.0],
        'col': [0.0, 0.0, 1.0, 1.0],
    })
    gene_meta = pd.DataFrame({
        'GeneID': [0, 1],
        'Marker': [-1, -1],
        'LRindex': [-1, -1],
        'Type_0': [4.0, 7.0],
        'Type_1': [4.0, 7.0],
    })
    return meta, gene_meta


def test_zero_direct_spatial_effect_matches_legacy_counts():
    meta, gene_meta = _small_omics_inputs()
    legacy = omics.simOmics(gene_meta, meta, seed=13)
    extended, latent_mean, gene_truth, design = omics.simOmicsWithSpatialEffects(
        gene_meta,
        meta,
        direct_spatial_effects={
            'genes': ['Gene_0'],
            'basis': 'linear',
            'coefficients': [0.0, 0.0],
        },
        seed=13,
    )

    assert_frame_equal(extended, legacy)
    assert_frame_equal(latent_mean, omics.buildOmicsMean(gene_meta, meta))
    assert gene_truth.loc[0, 'Gene'] == 'Gene_0'
    assert {'row_scaled', 'col_scaled', 'row_linear', 'col_linear'} <= set(design.columns)


def test_direct_spatial_effect_changes_only_selected_gene_mean():
    meta, gene_meta = _small_omics_inputs()
    baseline_mean = omics.buildOmicsMean(gene_meta, meta)
    _, spatial_mean, _, _ = omics.simOmicsWithSpatialEffects(
        gene_meta,
        meta,
        direct_spatial_effects={
            'genes': [0],
            'basis': 'linear',
            'coefficients': [2.0, 0.0],
        },
        seed=9,
    )

    assert spatial_mean.loc[1, 'Gene_0'] > spatial_mean.loc[0, 'Gene_0']
    assert spatial_mean.loc[3, 'Gene_0'] > spatial_mean.loc[2, 'Gene_0']
    assert_frame_equal(spatial_mean[['Gene_1']], baseline_mean[['Gene_1']])


def test_observation_components_are_exact_noops_when_disabled():
    counts = pd.DataFrame(
        [[0.0, 3.0], [4.0, 8.0]],
        columns=['Gene_0', 'Gene_1'],
    )
    captured, _ = omics.applyCaptureNoise(counts, capture_efficiency=1.0, seed=1)
    ambient, additions, _, _ = omics.applyAmbientNoise(
        captured,
        ambient_rate=0.0,
        seed=2,
    )
    observed, mask, probabilities = omics.applyDropout(
        ambient,
        dropout={'mode': 'constant', 'probability': 0.0},
        seed=3,
    )

    assert_frame_equal(captured, counts)
    assert_frame_equal(ambient, counts)
    assert_frame_equal(observed, counts)
    assert additions.to_numpy().sum() == 0
    assert not mask.to_numpy().any()
    assert probabilities.to_numpy().sum() == 0


def test_ambient_noise_adds_nonnegative_counts():
    counts = pd.DataFrame(
        np.zeros((20, 3), dtype=int),
        columns=[f'Gene_{gene_id}' for gene_id in range(3)],
    )
    observed, additions, rates, profile = omics.applyAmbientNoise(
        counts,
        ambient_rate=25.0,
        ambient_profile=[1.0, 2.0, 1.0],
        seed=10,
    )

    assert additions.to_numpy().sum() > 0
    assert (observed.to_numpy() >= counts.to_numpy()).all()
    assert (rates == 25.0).all()
    assert np.isclose(profile.sum(), 1.0)


def test_mean_dependent_dropout_is_lower_for_high_expression():
    counts = pd.DataFrame(
        np.full((2, 2), 5, dtype=int),
        columns=['Gene_0', 'Gene_1'],
    )
    latent_mean = pd.DataFrame(
        [[0.1, 1.0], [10.0, 100.0]],
        columns=counts.columns,
    )
    _, _, probabilities = omics.applyDropout(
        counts,
        dropout={'mode': 'mean_dependent', 'intercept': 1.0, 'slope': -2.0},
        latent_mean=latent_mean,
        seed=12,
    )

    assert probabilities.loc[0, 'Gene_0'] > probabilities.loc[0, 'Gene_1']
    assert probabilities.loc[0, 'Gene_1'] > probabilities.loc[1, 'Gene_0']
    assert probabilities.loc[1, 'Gene_0'] > probabilities.loc[1, 'Gene_1']


def test_capture_and_dropout_change_counts_in_declared_direction():
    counts = pd.DataFrame(
        np.full((30, 4), 10, dtype=int),
        columns=[f'Gene_{gene_id}' for gene_id in range(4)],
    )
    captured, _ = omics.applyCaptureNoise(counts, capture_efficiency=0.4, seed=11)
    dropped, mask, _ = omics.applyDropout(
        captured,
        dropout={'mode': 'constant', 'probability': 1.0},
        seed=12,
    )

    assert captured.to_numpy().sum() < counts.to_numpy().sum()
    assert mask.to_numpy().all()
    assert dropped.to_numpy().sum() == 0


def test_create_omics_default_matches_historical_generator_path():
    meta, _ = _small_omics_inputs()
    sim = SimSpace(shape=(2, 2), num_states=2, random_seed=23)
    sim.meta = meta.copy()
    sim.create_omics(
        n_genes=12,
        bg_ratio=0.5,
        lr_ratio=0.5,
        spatial=True,
        k_neighors=3,
    )

    historical_counts = omics.simOmics(sim.gene_meta, sim.meta, seed=sim.seed)
    historical_spatial_counts = omics.simSpatialOmics(
        historical_counts,
        sim.gene_meta,
        sim.meta,
        k_neighors=3,
        spatial_effect=3,
        se_threshold=1.5,
        seed=1,
    )

    assert_frame_equal(sim.omics, historical_spatial_counts)
    assert not hasattr(sim, 'omics_extension_config')


def test_empty_extension_configs_match_default_output():
    meta, _ = _small_omics_inputs()
    default = SimSpace(shape=(2, 2), num_states=2, random_seed=29)
    default.meta = meta.copy()
    default.create_omics(
        n_genes=12,
        bg_ratio=0.5,
        lr_ratio=0.5,
        spatial=True,
        k_neighors=3,
    )

    disabled = SimSpace(shape=(2, 2), num_states=2, random_seed=29)
    disabled.meta = meta.copy()
    disabled.create_omics(
        n_genes=12,
        bg_ratio=0.5,
        lr_ratio=0.5,
        spatial=True,
        k_neighors=3,
        technical_noise={},
        dropout={},
    )

    assert_frame_equal(disabled.gene_meta, default.gene_meta)
    assert_frame_equal(disabled.omics, default.omics)


def test_optional_pipeline_is_reproducible():
    meta, _ = _small_omics_inputs()
    options = {
        'n_genes': 8,
        'bg_ratio': 0.5,
        'lr_ratio': 0.0,
        'spatial': False,
        'direct_spatial_effects': {
            'genes': ['Gene_0'],
            'basis': 'linear',
            'coefficients': [1.0, -0.5],
        },
        'technical_noise': {
            'capture_efficiency': 0.75,
            'ambient_rate': 2.0,
            'seed': 301,
        },
        'dropout': {
            'mode': 'mean_dependent',
            'intercept': 0.5,
            'slope': -1.0,
            'seed': 302,
        },
        'store_intermediate': True,
    }

    first = SimSpace(shape=(2, 2), num_states=2, random_seed=31)
    first.meta = meta.copy()
    first.create_omics(**options)
    second = SimSpace(shape=(2, 2), num_states=2, random_seed=31)
    second.meta = meta.copy()
    second.create_omics(**options)

    assert_frame_equal(first.omics, second.omics)
    assert_frame_equal(first.omics_latent_mean, second.omics_latent_mean)
    assert_frame_equal(first.omics_dropout_mask, second.omics_dropout_mask)
    assert_frame_equal(first.spatial_gene_truth, second.spatial_gene_truth)


def test_create_omics_stores_optional_truth_layers():
    meta, _ = _small_omics_inputs()
    sim = SimSpace(shape=(2, 2), num_states=2, random_seed=17)
    sim.meta = meta
    sim.create_omics(
        n_genes=6,
        bg_ratio=0.5,
        lr_ratio=0.0,
        spatial=False,
        direct_spatial_effects={
            'genes': ['Gene_0'],
            'basis': 'hotspot',
            'center': [0.5, 0.5],
            'bandwidth': 0.4,
            'coefficients': [1.0],
        },
        technical_noise={
            'capture_efficiency': 0.8,
            'ambient_rate': 1.0,
            'seed': 101,
        },
        dropout={'mode': 'constant', 'probability': 0.1, 'seed': 102},
        store_intermediate=True,
    )

    assert sim.omics.shape == (len(meta), 6)
    assert sim.omics_latent_mean.shape == sim.omics.shape
    assert sim.omics_dropout_mask.shape == sim.omics.shape
    assert sim.spatial_gene_truth.loc[0, 'Gene'] == 'Gene_0'
    assert sim.omics_extension_config['technical_noise']['capture_efficiency'] == 0.8

def test_omics():
    """
    Test function.
    """
    shape = (80, 80)
    # Generate a simulation space
    sim1 = sim_from_json(
        os.path.join(
            os.path.dirname(__file__), 
            '../data', 
            'fitted_params.json'
        ),
        shape=shape, 
        num_iteration=4, 
        n_iter=6, 
        custom_neighbor=generate_offsets(3, 'manhattan'),
        seed=0
    )
    sim1.update_seed(seed=1)

    sim1.fit_scdesign(
        os.path.join(
            os.path.dirname(__file__), 
            '../data', 
            'reference_count.csv'
        ),
        os.path.join(
            os.path.dirname(__file__), 
            '../data', 
            'reference_metadata.csv'
        ),
        'Cluster',      # Cell annotation column in reference meta data to use
        'x_centroid',   # X coordinate column in reference meta data to use
        'y_centroid',   # Y coordinate column in reference meta data to use
        seed=0,
    )

    assert sim1 is not None, "Simulation space should not be None"
    assert sim1.omics['ABCC11'] is not None, "Omics data should not be None"
    assert len(sim1.omics) > 0, "Omics data should not be empty"
    assert 'fitted_celltype' in sim1.meta.columns, "Meta data should contain 'fitted_celltype' column"
    assert round(sim1.meta['row'][0], 2) == 0.35, "First row value should match expected value"
    assert round(sim1.meta['row'][1], 2) == 0.08, "Second row value should match expected value"

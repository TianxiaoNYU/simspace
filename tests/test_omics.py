import os
import warnings
warnings.filterwarnings("ignore")

from simspace.util import sim_from_json
from simspace.spatial import generate_offsets

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
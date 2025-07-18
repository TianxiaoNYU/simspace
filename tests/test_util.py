import pytest
from simspace.spatial import generate_offsets
from simspace import util

def test_save_params(tmp_path):
    params = util.generate_random_parameters(
        n_group=3,
        n_state=8,
        seed=42
    )
    output_file = tmp_path / "params.json"
    util.save_params(params, output_file)
    
    assert output_file.exists()
    
    loaded_params = util.load_params(output_file)
    print(loaded_params['density_replicates'])
    print(params['density_replicates'])
    assert set(params.keys()) == set(loaded_params.keys())

def test_generate_random_parameters():
    params = util.generate_random_parameters(
        n_group=3,
        n_state=5,
        seed=42
    )
    assert isinstance(params, dict)
    assert 'n_group' in params
    assert 'n_state' in params
    assert 'niche_theta' in params
    assert 'theta_list' in params
    assert 'density_replicates' in params
    assert 'phi_replicates' in params

def test_sim_from_params():
    params = util.generate_random_parameters(
        n_group=3,
        n_state=8,
        seed=42
    )
    sim = util.sim_from_params(
        params,
        shape = (100, 100),
        num_iteration=4, 
        n_iter=6, 
        custom_neighbor=generate_offsets(3, 'manhattan'),
        seed=0
    )
    
    assert hasattr(sim, 'meta')
    assert round(sim.meta.row[1], 2) == 0.08
    assert round(sim.meta.row[15], 2) == 1.07
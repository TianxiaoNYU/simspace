import pytest
from simspace import SimSpace

def test_simspace_initialization():
    sim = SimSpace(
        shape=(10, 10), 
        num_states=5, 
        num_iterations=100
    )
    assert sim.shape == (10, 10)
    assert sim.num_states == 5
    assert sim.num_iterations == 100

import warnings
warnings.filterwarnings("ignore")

import os
import json
import numpy as np
import pandas as pd

import simspace

# Function to save parameters to a JSON file
def save_params(params, output_file):
    """
    Save genetic algorithm parameters to a JSON file with annotations.

    Args:
        params (dict): Dictionary containing parameter names and values.
        output_file (str): Path to the output JSON file.
    Raises:
        ValueError: If the parameters are not in the expected format or if required keys are missing.
        FileNotFoundError: If the output file cannot be created.
    Note:        The function converts NumPy arrays to lists before saving to ensure compatibility with JSON format.
        If the 'theta_list' contains NumPy arrays, they are converted to lists.
        If 'theta_list' is not a list of lists or NumPy arrays, a ValueError is raised.
    Returns:
        None
    Example:
        params = generate_random_parameters(
            n_group=3, 
            n_state=5, 
            seed=42
        )
        save_parameters_to_json(params, 'params.json')
    This will save the parameters to 'params.json' with each value converted to a list if it is a NumPy array.
    """

    if isinstance(params['theta_list'][0], np.ndarray):
        params['theta_list'] = [theta.tolist() for theta in params['theta_list']]
    elif not isinstance(params['theta_list'][0], list):
        raise ValueError(f"Theta list should be a list of lists or numpy arrays. Got a list of {type(params['theta_list'][0])}")
    
    param_tolist = {key: {"value": value.tolist() if isinstance(value, np.ndarray) else value} for key, value in params.items()}
    if not isinstance(param_tolist, dict):
        raise ValueError("Parameters must be a dictionary.")
    expected_keys = ['n_group', 'n_state', 'niche_theta', 'theta_list', 'density_replicates', 'phi_replicates']
    for key in expected_keys:
        if key not in param_tolist:
            raise ValueError(f"Missing expected parameter key: {key}")
        
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        with open(output_file, 'w') as json_file:
            json.dump(param_tolist, json_file, indent=4)
    except Exception as e:
        raise FileNotFoundError(f"Could not create or write to file {output_file}: {e}")

# Function to load parameters from a JSON file
def load_params(input_file):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File {input_file} does not exist.")
    try:
        with open(input_file, 'r') as json_file:
            data = json.load(json_file)
        return {key: value['value'] for key, value in data.items()}
    except json.JSONDecodeError as e:
        raise ValueError(f"Error decoding JSON: {e}")
    
# Function to randomly generate the parameters
def generate_random_parameters(
        n_group, 
        n_state, 
        seed=0,
        ):
    """
    Generate random parameters for the simulation.

    Args:
        n_group (int): Number of groups.
        n_state (int): Number of states.
        n_subgroup (int): Number of subgroups.
        n_replicates (int): Number of replicates.
        n_neighbors (int): Number of neighbors.

    Returns:
        dict: Dictionary containing the generated parameters.
    Notes:
        - The function generates random values for niche theta, theta list, density replicates, and
          phi replicates.
        - The theta values are generated uniformly within specified ranges.
        - The function uses a seed for reproducibility.
    """
    np.random.seed(seed)
    theta_list = []
    for _ in range(n_group):
        theta = np.random.uniform(-0.8, 0.8, size=(n_state-1)*n_state//2)
        theta_list.append(theta)
    parameters = {
        'n_group': n_group,
        'n_state': n_state,
        'niche_theta': np.random.uniform(-0.5, 0.5, size=(n_group-1)*n_group//2),
        'theta_list': theta_list,
        'density_replicates': np.random.uniform(0.01, 0.4, size=n_state),
        'phi_replicates': np.random.uniform(4.4, 5),
    }

    return parameters

# Function to simulate from parameters
def sim_from_params(
    parameters: dict, 
    shape: tuple, 
    num_iteration: int, 
    n_iter: int, 
    custom_neighbor: callable, 
    step: float = 0.2,
    seed: int = 0
):
    n_group = parameters['n_group']
    n_state = parameters['n_state']
    
    niche_theta = np.zeros((n_group, n_group))
    niche_theta[np.triu_indices(n_group, 1)] = parameters['niche_theta']
    niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
    np.fill_diagonal(niche_theta, 1)

    theta_list = []
    for i in range(n_group):
        theta_tmp = np.zeros((n_state, n_state))
        theta_tmp[np.triu_indices(n_state, 1)] = parameters['theta_list'][i]
        theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
        np.fill_diagonal(theta_tmp, 1)
        theta_list.append(theta_tmp)
    
    density_replicates = np.array(parameters['density_replicates'])
    density_replicates[density_replicates < 0] = 0
    phi_replicates = parameters['phi_replicates']

    Sim = simspace.SimSpace(
        shape = shape,
        num_states = n_state,
        num_iterations= num_iteration,
        theta=theta_list,
        phi=phi_replicates,
        neighborhood=custom_neighbor, 
        random_seed=seed,
        )
    Sim.initialize()   # Initialize the grid
    Sim.create_niche(num_niches=n_group, n_iter=n_iter, theta_niche=niche_theta)
    Sim.gibbs_sampler()    # Gibbs sampling
    Sim.density_sampler(density_replicates)  # Cell density of each niche
    Sim.perturbation(step = step)     # Perturbation

    return Sim

# Function to simulate from saved json files
def sim_from_json(
    input_file: str, 
    shape: tuple, 
    num_iteration: int, 
    n_iter: int, 
    custom_neighbor: callable, 
    seed: int = 0
):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File {input_file} does not exist.")
    try:
        parameters = load_params(input_file)
    except json.JSONDecodeError as e:
        raise ValueError(f"Error decoding JSON: {e}")
    
    n_group = parameters['n_group']
    diag_length = len(parameters['theta_list'][0])
    n_state = int((1 + np.sqrt(1 + 8 * diag_length)) / 2)
    if n_state * (n_state - 1) != 2 * diag_length:
        raise ValueError("Invalid theta matrix size.")

    niche_theta = np.zeros((n_group, n_group))
    niche_theta[np.triu_indices(n_group, 1)] = parameters['niche_theta']
    niche_theta = niche_theta + niche_theta.T - np.diag(niche_theta.diagonal())
    np.fill_diagonal(niche_theta, 1)

    theta_list = []
    for i in range(n_group):
        theta_tmp = np.zeros((n_state, n_state))
        theta_tmp[np.triu_indices(n_state, 1)] = parameters['theta_list'][i]
        theta_tmp = theta_tmp + theta_tmp.T - np.diag(theta_tmp.diagonal())
        np.fill_diagonal(theta_tmp, 1)
        theta_list.append(theta_tmp)

    density_replicates = np.array(parameters['density_replicates'])
    density_replicates[density_replicates < 0] = 0
    phi_replicates = parameters['phi_replicates']

    Sim = simspace.SimSpace(
        shape = shape,
        num_states = n_state,
        num_iterations= num_iteration,
        theta=theta_list,
        phi=phi_replicates,
        neighborhood=custom_neighbor, 
        random_seed=seed,
        )
    Sim.initialize()   # Initialize the grid
    Sim.create_niche(num_niches=n_group, n_iter=n_iter, theta_niche=niche_theta)
    Sim.gibbs_sampler()    # Gibbs sampling
    Sim.density_sampler(density_replicates)  # Cell density of each niche
    Sim.perturbation(step = 0.2)     # Perturbation

    return Sim

# Function to convolve omics data with a kernel
def convolve(
        simspace: object, 
        kernel: np.array,
        scale: int = 1,
        conv_type = 'average'):
    """
    Convolve the omics data with a kernel by averaging.

    Args:
        simspace (object): SimSpace object containing cell metadata and omics data.
        kernel (tuple): Size of the kernel as a tuple (width, height).
        scale (int): Scaling factor for the omics data. Defaults to 1.
        conv_type (str): Type of convolution to perform ('average' or 'sum'). Defaults to 'average'.
    
    Returns:
        tuple: A tuple containing two DataFrames:
            - spot_meta: Metadata for the spots, including their coordinates and state proportions.
            - spot_omics: Omics data for the spots, either averaged or summed based on conv_type.

    Examples:
        >>> simspace = SimSpace(...)  # Initialize your SimSpace object
        >>> kernel = (5, 5)
        >>> spot_meta, spot_omics = convolve(simspace, kernel, scale=1, conv_type='average')
        >>> print(spot_meta.head())
        >>> print(spot_omics.head())
    """
    cell_meta = simspace.meta
    omics = simspace.omics

    if 'fitted_celltype' in cell_meta.columns:
        group_col = 'fitted_celltype'
    else:
        group_col = 'state_rank'

    cell_meta.index = [f'cell_{i+1}' for i in range(len(cell_meta))]

    spot_meta = []
    spot_omics = []
    state_list = np.sort(cell_meta[group_col].unique())
    max_x = cell_meta['col'].max().astype(int)
    max_y = cell_meta['row'].max().astype(int)
    for i in range(0, max_x, kernel[0]):
        for j in range(0, max_y, kernel[1]):
            region = cell_meta[(cell_meta['col'] >= i) & (cell_meta['col'] < i + kernel[0]) & 
                               (cell_meta['row'] >= j) & (cell_meta['row'] < j + kernel[1])]
            if not region.empty:
                centroid_x = kernel[0] // 2 + i
                centroid_y = kernel[1] // 2 + j
                state_proportions = region[group_col].value_counts(normalize=True, sort=False)
                state_proportions = state_proportions.reindex(state_list, fill_value=0)
                spot_meta.append([centroid_x, centroid_y] + state_proportions.tolist())

                if conv_type == 'average':
                    spot_omic = omics.loc[region.index].sum() * scale
                elif conv_type == 'sum':
                    spot_omic = omics.loc[region.index].sum()
                else:
                    raise ValueError("Invalid conv_type. Choose 'average' or 'sum'.")
                spot_omic = spot_omic.round().astype(int)
                spot_omics.append(spot_omic)
    
    spot_meta = pd.DataFrame(spot_meta, columns=['col', 'row'] + list(state_list))
    spot_omics = pd.DataFrame(spot_omics, columns=omics.columns)

    return spot_meta, spot_omics


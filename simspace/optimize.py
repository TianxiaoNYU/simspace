## Import necessary libraries
import numpy as np
from multiprocessing import Pool
import os
import copy

import simspace as ss

import warnings
warnings.filterwarnings('ignore')


## Define the Fitness Function, which do the simulation first and then calculate the fitness value
#### The fitness value is the distance between the simulation result and the experimental data
def _get_spatial_statistics(simulation, n_state):
    simulation.columns = ['celltype', 'x', 'y', 'celltype_ranked']
    tmp = simulation[simulation['celltype'] == -1].index
    simulation = simulation.drop(tmp)
    simulation = simulation.reset_index(drop=True)
    simulation['celltype'] = simulation['celltype'].astype(int)
    cell_type = simulation['celltype'].unique()
    cell_type = np.sort(cell_type)

    mi_sim = ss.spatial.integrate_morans_I(simulation['celltype'], simulation[['x', 'y']], cell_type)
    if len(mi_sim) < n_state:
        mi_sim = np.pad(mi_sim, (0, n_state - len(mi_sim)), 'constant', constant_values=np.nan)
    cell_counts = simulation['celltype'].value_counts().sort_index()
    mi_sim = [mi_sim[i] for i in np.argsort(-cell_counts)]  
    if simulation[['x', 'y']].shape[0] < 40:
        le_pdf = np.zeros(11)  
    else:  
        le_sim = ss.spatial.calculate_local_entropy(simulation['celltype'], simulation[['x', 'y']], k = 10) # k = 20 in Michelle's data
        le_pdf, _ = np.histogram(le_sim, bins=[x / 4 for x in range(0, 12)], density=True)

    mi_sim = np.array(mi_sim)
    sim_res = np.hstack((le_pdf, mi_sim))
    
    return sim_res

## Define the Crossover Function, which generates new individuals by combining the genes of the selected individuals
def _crossover(parent1, parent2, crossover_rate):
    if parent1['n_group'] != parent2['n_group']:
        warnings.warn("Parents have different number of groups")
        return parent1

    if np.random.rand() < crossover_rate:
        child = {
            "n_group": parent1['n_group'],
            "n_state": parent1['n_state'],
            "niche_theta": np.where(np.random.rand() < 0.5, parent1['niche_theta'], parent2['niche_theta']),
            "theta_list": [np.where(np.random.rand() < 0.5, parent1['theta_list'][i], parent2['theta_list'][i]) for i in range(parent1['n_group'])],
            "density_replicates": np.where(np.random.rand() < 0.5, parent1['density_replicates'], parent2['density_replicates']),
            "phi_replicates": float(np.where(np.random.rand() < 0.5, parent1['phi_replicates'], parent2['phi_replicates']))
        }
        return child
    return parent1

## Define the Mutation Function, which introduces random changes to the genes of the selected individuals
def _mutate(individual, mutation_rate):
    for i in range(len(individual['niche_theta'])):
        if np.random.rand() < mutation_rate:
            individual['niche_theta'][i] += np.random.normal(0, 0.2)
    
    for i in range(len(individual['theta_list'])):
        for j in range(len(individual['theta_list'][i])):
            if np.random.rand() < mutation_rate:
                individual['theta_list'][i][j] += np.random.normal(0, 0.2)

    for i in range(len(individual['density_replicates'])):
        if np.random.rand() < mutation_rate:
            individual['density_replicates'][i] += np.random.normal(0, 0.05)
        if individual['density_replicates'][i] < 0:
            individual['density_replicates'][i] = 0.01

    if np.random.rand() < mutation_rate:
        individual['phi_replicates'] += np.random.normal(0, 0.1)
    return individual

def _initialize_population(population_size, n_group, n_state, group_mean = 3, seed=0):
    """
    Initialize a population of simulation parameters for the evolutionary algorithm.
    Args:
        population_size: Number of individuals in the population.
        n_group: Number of groups in the simulation.
        n_state: Number of states in the simulation.
        group_mean: Mean number of groups (default is 3).
        seed: Random seed for reproducibility (default is 0).
    Returns:
        param_dist_list: List of dictionaries containing the simulation parameters for each individual.
    """

    np.random.seed(seed)

    param_dist_list = []
    if n_group < 0:
        n_group = n_state // group_mean
        if n_group < 1:
            n_group = 1

    for i in range(population_size):
        # subgroup_length = split_group(n_state, group_mean)
        niche_theta = np.random.uniform(-0.5, 0.5, size=(n_group-1)*n_group//2)
        # theta = np.random.uniform(-0.8, 0.8, size=(n_state-1)*n_state//2)
        theta_list = []
        for i in range(n_group):
            theta = np.random.uniform(-0.8, 0.8, size=(n_state-1)*n_state//2)
            theta_list.append(theta)
        density_replicates = np.random.uniform(0.01, 0.4, size=n_state)
        phi_replicates = np.random.uniform(4.4, 5)

        param_dist = {
            "n_group": n_group, 
            "n_state": n_state,
            "niche_theta": niche_theta, 
            "theta_list": theta_list, 
            "density_replicates": density_replicates, 
            "phi_replicates": phi_replicates}

        param_dist_list.append(param_dist)

    return param_dist_list

def _parallel_fitness_evaluation(
        population, 
        target, 
        shape,
        custom_neighbor,
        num_iterations,
        n_iter,
        parallel=True,
        replicate = 1):
    """
    Evaluate the fitness of the entire population in parallel.
    """
    if parallel:
        cpu_count = _get_cpu_count()
    else:
        cpu_count = 1
    with Pool(processes=cpu_count) as pool:
        fitness_scores = pool.starmap(_fitness_function, 
                                      [(ind, shape, custom_neighbor, num_iterations, n_iter, target, replicate) for ind in population])
    return np.array(fitness_scores)

def _get_cpu_count():
    """
    Get the number of CPUs allocated by SLURM or default to all available CPUs.
    """
    return int(os.getenv('SLURM_CPUS_PER_TASK', os.cpu_count()))

def _fitness_function(
        parameters, 
        shape,
        custom_neighbor,
        num_iterations,
        n_iter,
        target_vector, 
        replicate = 1):

    if replicate == 1:
        Sim1 = ss.util.sim_from_params(
            parameters=parameters,
            shape=shape,
            custom_neighbor=custom_neighbor,
            num_iteration=num_iterations,
            n_iter=n_iter,
            seed=0)
        sim_res = _get_spatial_statistics(Sim1.meta, parameters['n_state'])

        if len(sim_res) != len(target_vector):
             fitness = np.linalg.norm(sim_res)
        else:
            fitness = np.linalg.norm(sim_res - target_vector)
    elif replicate > 1:
        assert isinstance(replicate, int), "Replicate should be an integer"
        fitness = []
        for i in range(replicate):
            Sim = ss.util.sim_from_params(
                parameters=parameters,
                shape=shape,
                custom_neighbor=custom_neighbor,
                num_iteration=num_iterations,
                n_iter=n_iter,
                seed=i)
            sim_res = _get_spatial_statistics(Sim.meta, parameters['n_state'])
            if len(sim_res) != len(target_vector):
                fitness.append(np.linalg.norm(sim_res))
            else:
                fitness.append(np.linalg.norm(sim_res - target_vector))
        fitness = np.mean(fitness)
    else:
        raise ValueError("Replicate should be a positive integer")
        
    return fitness

def _tournament_selection(population, fitness_scores, k=3):
    selected = []
    for _ in range(len(population)):
        indices = np.random.choice(len(population), k, replace=False)
        best = indices[np.argmin(fitness_scores[indices])]
        selected.append(population[best])

    return selected


def spatial_fit(
        target,
        population_size: int=50,
        generations: int=20,
        mutation_rate: float=0.2,
        crossover_rate: float=0.6,
        shape: tuple=(50, 50),
        n_group: int=2,
        n_state: int=8,
        custom_neighbor: list=ss.spatial.generate_offsets(3, 'manhattan'),
        num_iterations: int=4,
        n_iter: int=6,
        replicate: int=1,
        seed: int=0,
        parallel: bool=True,
        verbose: bool=True
        ):
    """
    Perform the Evolutionary Algorithm to optimize simulation parameters based on a target vector.
    
    Args:
        target (list): Target vector to optimize against.
        population_size (int): Number of individuals in the population (default is 50).
        generations (int): Number of generations to run the algorithm (default is 20).
        mutation_rate (float): Probability of mutation for each individual (default is 0.2).
        crossover_rate (float): Probability of crossover between individuals (default is 0.6).
        shape (tuple): Shape of the simulation grid (default is (50, 50)).
        n_group (int): Number of groups in the simulation (default is 2).
        n_state (int): Number of states in the simulation (default is 8).
        custom_neighbor (list): Custom neighbor offsets for the simulation (default is None).
        num_iterations (int): Number of iterations for the simulation (default is 4).
        n_iter (int): Number of iterations for the simulation (default is 6).
        replicate (int): Number of replicates for the simulation (default is 1).
        seed (int): Random seed for reproducibility (default is 0).
        parallel (bool): Whether to run the fitness evaluation in parallel (default is True).
        verbose (bool): Whether to print progress information (default is True).
    
    Returns:
        best_solution (dict): The best solution found during the optimization process.
    
    Raises:
        ValueError: If the target vector is not a list or if the population size is not a positive integer.
    """
    if not isinstance(target, list):
        raise ValueError("Target should be a list containing local entropy and Moran's I values.")
    if not isinstance(population_size, int) or population_size <= 0:
        raise ValueError("Population size should be a positive integer.")

    ## Define (or load) the initial population, which is the simulation parameters
    population = _initialize_population(population_size, n_group, n_state, seed=seed)

    ## Call the Evolutionary Algorithm Function to perform the EA
    best_solution = None
    best_fitness = float('inf')

    for generation in range(generations):
        # Evaluate fitness
        fitness_scores = _parallel_fitness_evaluation(
            population=population, 
            target=target, 
            replicate=replicate,
            shape=shape,
            custom_neighbor=custom_neighbor,
            num_iterations=num_iterations,
            n_iter=n_iter,
            parallel=parallel
            )
        
        # Track the best solution
        best_gen_idx = np.argmin(fitness_scores)
        if fitness_scores[best_gen_idx] < best_fitness:
            best_fitness = fitness_scores[best_gen_idx]
            best_solution = copy.deepcopy(population[best_gen_idx])

        if verbose:
            print(f"Generation {generation}: Best Fitness = {best_fitness}")

        # Selection
        selected_population = _tournament_selection(population, fitness_scores)

        # Crossover and mutation
        next_generation = []
        for i in range(0, population_size, 2):
            parent1 = selected_population[i]
            parent2 = selected_population[(i+1) % population_size]
            child1 = _crossover(parent1, parent2, crossover_rate)
            child2 = _crossover(parent2, parent1, crossover_rate)
            next_generation.append(_mutate(child1, mutation_rate))
            next_generation.append(_mutate(child2, mutation_rate))

        population = next_generation

    if verbose:
        print("Optimization complete!")
        print("Best solution:", best_solution)
        print("Best fitness:", best_fitness)

    return best_solution
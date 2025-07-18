import numpy as np

def initialize_configuration(shape, num_states):
    """
    Initialize the configuration z with random states in a 2-D array.
    :param shape: The shape of the 2-D configuration (height, width).
    :param num_states: The number of possible states for each pixel.
    :return: The initial configuration.
    """
    return np.random.randint(num_states, size=shape)

def compute_hk(z, theta, k, i, j, neighbors):
    """
    Compute h_k(z) for a given pixel (i, j).
    :param z: Current configuration.
    :param theta: Parameter matrix.
    :param k: State to compute the sum for.
    :param i: Row index of the current pixel.
    :param j: Column index of the current pixel.
    :param neighbors: List of neighbor indices.
    :return: Computed h_k(z).
    """
    hk = 0
    for ni, nj in neighbors:
        hk += theta[k, z[ni, nj]]
    return hk

def conditional_probability(z, theta, k, i, j, neighbors):
    """
    Compute the conditional probability P(Z_ij = k | Z_-ij = z_-ij).
    :param z: Current configuration.
    :param theta: Parameter matrix.
    :param k: State to compute the probability for.
    :param i: Row index of the current pixel.
    :param j: Column index of the current pixel.
    :param neighbors: List of neighbor indices.
    :return: Conditional probability.
    """
    hk_values = [compute_hk(z, theta, k_prime, i, j, neighbors) for k_prime in range(theta.shape[0])]
    exp_values = np.exp(hk_values)
    return np.exp(hk_values[k]) / np.sum(exp_values)
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import libpysal
from esda.moran import Moran
from esda.moran import Moran_Local
from libpysal.weights import KNN
from sklearn.neighbors import NearestNeighbors

def generate_offsets(distance, 
                     method: str = 'manhattan',
                     linear: bool = False):
    """
    Generate neighbor offsets based on the specified distance and method.

    Args:
        distance (int): Distance parameter.
        method (str): Method to generate offsets ('manhattan' or 'euclidean').
        linear (bool): Whether to include the cell itself in the offsets (default is False).

    Returns:
        list: List of generated neighbor offsets.

    Raises:
        ValueError: If the distance is not an integer or if the method is not recognized.
    
    Examples:
        >>> from simspace.spatial import generate_offsets
        >>> offsets = generate_offsets(1, 'manhattan')
        >>> print(offsets)
        [(-1, 0), (0, -1), (0, 1), (1, 0)]
        >>> offsets = generate_offsets(2, 'euclidean', linear=True)
        >>> print(offsets)
        [(-2, 0), (-1, -1), (-1, 0), (-1, 1), (0, -2), (0, -1), (0, 1), (0, 2), (1, -1), (1, 0), (1, 1), (2, 0), (0, 0)]
    """

    if not isinstance(distance, int):
        raise ValueError("The distance must be an integer.")
    
    if method not in ['manhattan', 'euclidean']:
        raise ValueError("The method should be either 'manhattan' or 'euclidean'.")

    offsets = []
    for i in range(-distance, distance+1):
        for j in range(-distance, distance+1):
            if method == 'manhattan':
                if abs(i) + abs(j) <= distance and (i, j) != (0, 0):
                    offsets.append((i, j))
            elif method == 'euclidean':
                if np.sqrt(i**2 + j**2) <= distance and (i, j) != (0, 0):
                    offsets.append((i, j))
    if linear:
        offsets.append((0, 0))  # Include the cell itself if linear is True
    return offsets

def generate_offsets3D(distance, method, linear=False):
    """
    Generate 3D neighbor offsets based on the specified distance and method.

    Args:
        distance (int): Distance parameter.
        method (str): Method to generate offsets ('manhattan' or 'euclidean').
        linear (bool): Whether to include the cell itself in the offsets (default is False).

    Returns:
        list: List of generated neighbor offsets.

    Raises:
        ValueError: If the distance is not an integer or if the method is not recognized.
    
    Examples:
        >>> offsets = generate_offsets3D(3, 'manhattan')
        >>> print(offsets)
    """

    if not isinstance(distance, int):
        raise ValueError("The distance must be an integer.")
    
    if method not in ['manhattan', 'euclidean']:
        raise ValueError("The method should be either 'manhattan' or 'euclidean'.")

    offsets = []
    for i in range(-distance, distance+1):
        for j in range(-distance, distance+1):
            for k in range(-distance, distance+1):
                if method == 'manhattan':
                    if abs(i) + abs(j) + abs(k)*2 <= distance and (i, j, k) != (0, 0, 0):
                        offsets.append((i, j, k))
                elif method == 'euclidean':
                    if np.sqrt(i**2 + j**2 + 4*k**2) <= distance and (i, j, k) != (0, 0, 0):
                        offsets.append((i, j, k))
    if linear:
        offsets.append((0, 0, 0))
    return offsets

## Calculate Moran's I
def calculate_morans_I(
        data: pd.DataFrame, 
        coordinates: pd.DataFrame,
        k = 5) -> float:
    """
    Calculate Moran's I for a given dataset and spatial weights.

    Args:
        data: pandas DataFrame containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates. Used for libpysal.cg.KDTree()
        k: number of nearest neighbors to consider for spatial weights. Default is 5.

    Returns:
        morans_I: Moran's I value.
    """
    kd = libpysal.cg.KDTree(coordinates)
    weights = KNN(kd, k=k)

    # Calculate Moran's I
    morans_I = Moran(data, weights)

    return morans_I.I

def integrate_morans_I(data: pd.DataFrame, 
                       coordinates: pd.DataFrame, 
                       typelist: list) -> list:
    """
    Calculate Moran's I for a given dataset and spatial weights.

    Args:
        data: pandas DataFrame containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates. Used for libpysal.cg.KDTree()
        typelist: list of types to calculate Moran's I for.
    
    Returns:
        mi_list: List of Moran's I values for each type in typelist.
    
    Raises:
        ValueError: If typelist is empty.
    """
    if not isinstance(typelist, list) or len(typelist) == 0:
        raise ValueError("typelist must be a non-empty list.")
    mi_list = []
    for type in typelist:
        tmp = data == type
        mi = calculate_morans_I(tmp, coordinates)
        mi_list.append(mi)
    return mi_list
    

## Calculate local Moran's I
def calculate_local_morans_I(data: pd.DataFrame, 
                             coordinates: pd.DataFrame, 
                             k: int = 20) -> np.ndarray:
    """
    Calculate local Moran's I for a given dataset and spatial weights.

    Args:
        data: pandas DataFrame or Series containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        k: number of nearest neighbors to consider for spatial weights.

    Returns:
        local_morans_I: Local Moran's I values.
    """
    kd = libpysal.cg.KDTree(coordinates)
    weights = libpysal.weights.KNN(kd, k=k)

    # Calculate local Moran's I
    local_morans_I = Moran_Local(data, weights)

    return local_morans_I.Is

## Plot local Moran's I
def plot_local_morans_I(
        data: pd.DataFrame, 
        coordinates: pd.DataFrame, 
        local_morans_I: np.ndarray, 
        ax=None):

    """
    Plot local Moran's I values on a scatter plot.

    Args:
        data: pandas DataFrame or Series containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        local_morans_I: Local Moran's I values.
        ax: matplotlib axis object to plot on.

    Returns:
        ax: matplotlib axis object.
    """
    if ax is None:
        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    else:
        ax1, ax2 = ax

    # Plot data
    ax1.scatter(coordinates.iloc[:, 0], coordinates.iloc[:, 1], c=data, cmap='viridis', s=5)

    # Plot local Moran's I values
    ax2.scatter(coordinates.iloc[:, 0], coordinates.iloc[:, 1], c=local_morans_I, s=5, alpha=0.75)

    return ax1, ax2

## Calculate the local entropy
def calculate_local_entropy(
        data: pd.DataFrame,
        coordinates: pd.DataFrame,
        k: int = 20) -> np.ndarray:
    """
    Calculate the local entropy for a given dataset and spatial coordinates.

    Args:
        data: pandas DataFrame or Series containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        k: number of nearest neighbors to consider for spatial weights.

    Returns:
        local_entropy: Local entropy values.
    """
    data.reset_index(drop=True, inplace=True)
    nbrs = NearestNeighbors(n_neighbors=k).fit(coordinates)
    _, indices = nbrs.kneighbors(coordinates)

    # Calculate the entropy for each cell
    local_entropy = []
    for i in range(len(data)):
        neighbors = indices[i]
        neighbor_values = data[neighbors]
        _, value_counts = np.unique(neighbor_values, return_counts=True)
        probabilities = value_counts / np.sum(value_counts)
        entropy = -np.sum(probabilities * np.log2(probabilities))
        local_entropy.append(entropy)

    return local_entropy

## histogram of local entropy
def plot_local_entropy(
        local_entropy: np.ndarray, 
        ax=None) -> plt.Axes:
    """
    Plot a histogram of local entropy values.

    Args:
        local_entropy: Local entropy values.
        ax: matplotlib axis object to plot on.

    Returns:
        ax: matplotlib axis object.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    ax.hist(local_entropy, bins=50, alpha=0.75)
    ax.set_xlabel('Local entropy')
    ax.set_ylabel('Frequency')

    return ax

def spatial_stat(
    data: pd.DataFrame,
    coordinates: pd.DataFrame,
    typelist: list) -> np.ndarray:
    """
    Calculate moran's I and local entropy for a given dataset.

    Args:
        data: pandas DataFrame containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        typelist: list of types to calculate Moran's I for. 

    Returns:
        res: numpy array containing moran's I and local entropy values.
    """
    morans_I = integrate_morans_I(data, coordinates, typelist)
    local_entropy = calculate_local_entropy(data, coordinates)

    res = np.array([morans_I, local_entropy])
    return res
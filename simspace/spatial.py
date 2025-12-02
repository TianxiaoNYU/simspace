import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import libpysal
from esda.moran import Moran
from esda.moran import Moran_Local
from esda.geary import Geary
from libpysal.weights import KNN
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import cKDTree
import scipy.spatial.distance as sdist

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
                       typelist) -> list:
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
    if len(typelist) == 0:
        raise ValueError("typelist must be a non-empty list or non-empty array.")
    mi_list = []
    for type in typelist:
        tmp = data == type
        mi = calculate_morans_I(tmp, coordinates)
        mi_list.append(mi)
    return mi_list

## Calculate Geary's C
def calculate_gearys_C(data: pd.DataFrame, 
                       coordinates: pd.DataFrame, 
                       k: int = 20) -> float:
    """
    Calculate Geary's C for a given dataset and spatial weights.
    Args:
        data: pandas DataFrame or Series containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        k: number of nearest neighbors to consider for spatial weights.
    """

    kd = libpysal.cg.KDTree(coordinates)
    weights = KNN(kd, k=k)

    if isinstance(data, np.ndarray) and data.dtype == bool:
        data = data.astype(int)
    elif isinstance(data, pd.Series) and data.dtype == bool:
        data = data.astype(int)
    elif isinstance(data, pd.DataFrame):
        for col in data.columns:
            if data[col].dtype == bool:
                data[col] = data[col].astype(int)

    # Calculate Geary's C
    gearys_C = Geary(data, weights)

    return gearys_C.C

def integrate_gearys_C(data: pd.DataFrame, 
                       coordinates: pd.DataFrame, 
                       typelist) -> list:
    """
    Calculate Geary's C for a given dataset and spatial weights.
    Args:
        data: pandas DataFrame containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates. Used for libpysal.cg.KDTree()
        typelist: list of types to calculate Geary's C for.

    Returns:
        gc_list: List of Geary's C values for each type in typelist.

    Raises:
        ValueError: If typelist is empty.
    """
    if len(typelist) == 0:
        raise ValueError("typelist must be a non-empty list or non-empty array.")
    gc_list = []
    for type in typelist:
        tmp = data == type
        gc = calculate_gearys_C(tmp, coordinates)
        gc_list.append(gc)
    return gc_list

## Calculate variogram
def compute_variogram(
        coords: np.ndarray, 
        labels: np.ndarray, 
        cell_types: list, 
        n_bins: int = 20, 
        max_dist: float = None):
    """
    Compute cross-variograms between all pairs of cell types.

    Args:
        coords: numpy array of shape (n_cells, n_dimensions) containing the spatial coordinates.
        labels: numpy array of shape (n_cells,) containing the cell type labels.
        cell_types: list of unique cell types to consider.
        n_bins: number of distance bins to use for the variogram (default is 20).
        max_dist: maximum distance to consider for the variogram (default is None, which uses half the maximum distance between points).
    
    Returns:
        variograms: dictionary where keys are tuples of cell type pairs and values are tuples of (bin_centers, gamma).
    """
    tree = cKDTree(coords)
    if max_dist is None:
        max_dist = np.linalg.norm(coords.max(0) - coords.min(0)) / 2

    bins = np.linspace(0, max_dist, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    
    variograms = {}
    for i, c1 in enumerate(cell_types):
        I1 = (labels == c1).astype(float)
        for j, c2 in enumerate(cell_types[i:], i):
            I2 = (labels == c2).astype(float)
            diffs = I1[:, None] - I2[None, :]
            dists = tree.sparse_distance_matrix(tree, max_dist=max_dist, output_type='ndarray')
            # compute squared differences only for pairs within max_dist
            gamma = np.zeros(n_bins)
            counts = np.zeros(n_bins)
            for k in range(len(dists)):
                d = dists['dist'][k]
                b = np.searchsorted(bins, d) - 1
                if b >= 0 and b < n_bins:
                    i1, i2 = dists['i'][k], dists['j'][k]
                    val = 0.5 * (I1[i1] - I2[i2])**2
                    gamma[b] += val
                    counts[b] += 1
            gamma /= np.maximum(counts, 1)
            variograms[(c1, c2)] = (bin_centers, gamma)
    return variograms

def integrate_variogram(
        data: pd.DataFrame,
        coordinates: pd.DataFrame, 
        typelist) -> dict:
    """
    Compute the average variogram across all cell type pairs.

    Args:
        data: pandas DataFrame containing the variable of interest.
        coordinates: numpy array or pandas DataFrame containing the spatial coordinates.
        typelist: list of unique cell types to consider.

    Returns:
        variograms: dictionary where keys are tuples of cell type pairs and values are tuples of (bin_centers, gamma).
        
    Raises:
        ValueError: If typelist is empty.
    """
    if len(typelist) == 0:
        raise ValueError("typelist must be a non-empty list or non-empty array.")
    variograms = compute_variogram(coordinates.values, data.values, typelist)
    return variograms

## Calculate pairwise distance between cell types as cell type interaction score
#### Idea was origninally from scFeatures
def calculate_interaction_score(
    data: pd.Series | pd.DataFrame,
    coordinates: pd.DataFrame | np.ndarray,
    typelist: list,
    k: int = 50,
    summary: str = "mean",
    use_knn: bool = True,
) -> pd.DataFrame:
    """
    Compute a type x type matrix of inter-type distances.

    By default uses k-NN distances (fast, memory-safe):
        - For c1 != c2: distances from each c1 cell to its k nearest c2 neighbors.
        - For c1 == c2: distances from each cell to its k nearest *other* cells of the same type
          (self-matches removed).

    If use_knn=False, computes full pairwise distances:
        - For c1 != c2: all |c1| x |c2| distances (uses scipy.spatial.distance.cdist).
        - For c1 == c2: all unique pairs within the type (uses pdist).

    Args:
        data: pd.Series or single-column pd.DataFrame of categorical labels (cell types).
        coordinates: (N, d) array of spatial coordinates.
        typelist: list of types to include in the output matrix (order preserved).
        k: number of neighbors to consider for k-NN (default: 20).
        summary: how to summarize the distances ("mean", "median", "min", "max").
        use_knn: whether to use k-NN distances (True, default) or full pairwise distances (False).  

    Returns:
        M: pd.DataFrame of shape (len(typelist), len(typelist)) with summarized inter-type distances.
    """
    # --- sanitize inputs ---
    if isinstance(data, pd.DataFrame):
        if data.shape[1] != 1:
            raise ValueError("`data` must be a pd.Series or a single-column DataFrame.")
        data = data.iloc[:, 0]
    if not isinstance(data, pd.Series):
        data = pd.Series(np.asarray(data).ravel(), index=getattr(coordinates, "index", None))
    # Align order with coordinates
    if isinstance(coordinates, pd.DataFrame):
        coords = coordinates.values
        labels = data.reindex(coordinates.index).values
    else:
        coords = np.asarray(coordinates)
        labels = np.asarray(data)

    if coords.ndim != 2:
        raise ValueError("`coordinates` must be 2D (N, d).")
    if coords.shape[0] != labels.shape[0]:
        raise ValueError("`coordinates` and `data` must have the same length.")

    # mapping from type -> indices
    type_to_idx = {t: np.where(labels == t)[0] for t in typelist}

    # helper to summarize an array
    def summarize(arr: np.ndarray) -> float:
        arr = np.asarray(arr, dtype=float)
        if arr.size == 0:
            return np.nan
        if summary == "mean":
            return float(np.mean(arr))
        if summary == "median":
            return float(np.median(arr))
        if summary == "min":
            return float(np.min(arr))
        if summary == "max":
            return float(np.max(arr))
        raise ValueError("summary must be one of {'mean','median','min','max'}")

    M = pd.DataFrame(index=typelist, columns=typelist, dtype=float)

    if use_knn:
        # Build one KDTree per c2 (reuse across all c1 queries)
        trees = {}
        for c2 in typelist:
            idx2 = type_to_idx[c2]
            if idx2.size == 0:
                trees[c2] = None
            else:
                trees[c2] = cKDTree(coords[idx2])

        for c1 in typelist:
            idx1 = type_to_idx[c1]
            if idx1.size == 0:
                M.loc[c1, :] = np.nan
                continue

            pts1 = coords[idx1]
            for c2 in typelist:
                idx2 = type_to_idx[c2]
                tree2 = trees[c2]

                if idx2.size == 0 or tree2 is None:
                    M.loc[c1, c2] = np.nan
                    continue

                if c1 == c2:
                    # query same set; remove self-matches
                    # ask for k+1 neighbors then drop the first (distance 0)
                    kk = min(k + 1, max(1, idx2.size))  # guard tiny sets
                    d, _ = tree2.query(pts1, k=kk)
                    d = np.atleast_2d(d)
                    if d.shape[1] == 1:
                        # only self available -> no valid neighbor
                        M.loc[c1, c2] = np.nan
                        continue
                    d = d[:, 1:]  # drop self distances
                else:
                    kk = min(k, max(1, idx2.size))
                    d, _ = tree2.query(pts1, k=kk)
                    d = np.atleast_2d(d)

                # summarize across all queried neighbor distances
                M.loc[c1, c2] = summarize(d.ravel())

    else:
        # Full pairwise (heavier): use scipy distance routines
        for c1 in typelist:
            idx1 = type_to_idx[c1]
            pts1 = coords[idx1]
            if idx1.size == 0:
                M.loc[c1, :] = np.nan
                continue

            for c2 in typelist:
                idx2 = type_to_idx[c2]
                pts2 = coords[idx2]
                if idx2.size == 0:
                    M.loc[c1, c2] = np.nan
                    continue

                if c1 == c2:
                    # all unique within-type pairs
                    if pts1.shape[0] < 2:
                        M.loc[c1, c2] = np.nan
                        continue
                    d = sdist.pdist(pts1)  # condensed vector
                else:
                    # cross-type all pairs (N1 x N2)
                    d = sdist.cdist(pts1, pts2).ravel()

                M.loc[c1, c2] = summarize(d)

    return M


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
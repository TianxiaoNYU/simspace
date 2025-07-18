import numpy as np

def get_neighbors(i, j, shape):
    """
    Get the nearest neighbors for a pixel (i, j) in a 2-D grid.
    :param i: Row index of the current pixel.
    :param j: Column index of the current pixel.
    :param shape: Shape of the 2-D grid.
    :return: List of neighbor indices.
    """
    neighbors = []
    if i > 0: neighbors.append((i-1, j))
    if i < shape[0] - 1: neighbors.append((i+1, j))
    if j > 0: neighbors.append((i, j-1))
    if j < shape[1] - 1: neighbors.append((i, j+1))
    return neighbors

def get_neighbors_8(i, j, shape):
    """
    Get the 8 nearest neighbors for a pixel (i, j) in a 2-D grid.
    :param i: Row index of the current pixel.
    :param j: Column index of the current pixel.
    :param shape: Shape of the 2-D grid.
    :return: List of neighbor indices.
    """
    neighbors = []
    for ni in range(max(0, i-1), min(shape[0], i+2)):
        for nj in range(max(0, j-1), min(shape[1], j+2)):
            if (ni, nj) != (i, j):
                neighbors.append((ni, nj))
    return neighbors

def get_custom_neighbors(i, j, shape, neighborhood):
    """
    Get the custom neighbors for a pixel (i, j) in a 2-D grid based on the specified neighborhood.
    :param i: Row index of the current pixel.
    :param j: Column index of the current pixel.
    :param shape: Shape of the 2-D grid.
    :param neighborhood: List of neighbor offsets relative to the current pixel.
    :return: List of neighbor indices.
    """
    neighbors = []
    for offset in neighborhood:
        ni = i + offset[0]
        nj = j + offset[1]
        if 0 <= ni < shape[0] and 0 <= nj < shape[1]:
            neighbors.append((ni, nj))
    return neighbors

def generate_offsets(distance, method):
    """
    Generate neighbor offsets based on the specified distance and method.
    :param distance: Distance parameter.
    :param method: Method to generate offsets ('manhattan' or 'euclidean').
    :return: List of neighbor offsets.
    """
    offsets = []
    for i in range(-distance, distance+1):
        for j in range(-distance, distance+1):
            if method == 'manhattan':
                if abs(i) + abs(j) <= distance and (i, j) != (0, 0):
                    offsets.append((i, j))
            elif method == 'euclidean':
                if np.sqrt(i**2 + j**2) <= distance and (i, j) != (0, 0):
                    offsets.append((i, j))
    return offsets

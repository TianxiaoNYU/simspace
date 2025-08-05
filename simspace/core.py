import numpy as np
import math
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import colorcet as cc
import typing
import pickle

from . import spatial, omics, niche

import warnings
warnings.filterwarnings('ignore')

class SimSpace:
    """
    SimSpace is a class for simulating spatial omics data in a 2-D or 3-D grid.

    Args:
        shape (tuple): The shape of the grid. Can be a tuple of two integers for 2-D or three integers for 3-D.
        num_states (int): The number of possible states.
        theta (numpy.ndarray): The theta matrix.
        num_iterations (int): The number of iterations for Gibbs sampling.
        phi (int or float): The phi parameter for Gibbs sampling.
        rho (int or float): The rho parameter for Gibbs sampling.
        neighborhood (list, optional): The custom neighborhood offsets. Defaults to an empty list.
        random_seed (int, optional): The random seed. Defaults to 1111.
    """
    def __init__(self, 
                 shape: tuple = (50, 50), 
                 num_states: int = 5, 
                 theta: np.ndarray = [], 
                 num_iterations: int = 5, 
                 phi: float = 5,
                 rho: float = 1,
                 grid: np.ndarray = [], 
                 neighborhood: list = [], 
                 random_seed: int = 1111) -> None:
        """
        Initialize the SimSpace object.
        """
        self.shape = shape
        self.num_states = num_states

        if isinstance(theta, list):
            self.theta = theta
        elif isinstance(theta, (np.ndarray, np.generic)):
            self.theta = [tmp for tmp in theta]
        else:
            raise ValueError("Invalid theta matrix.") 
        
        if not isinstance(num_iterations, int):
            raise ValueError("The number of iterations must be an integer.")
        self.num_iterations = num_iterations

        if isinstance(phi, int) or isinstance(phi, float):
            self.phi = phi
        else:
            raise ValueError("Invalid phi. Should be an integer or float.")
        
        if isinstance(rho, int) or isinstance(rho, float):
            self.rho = rho
        else:
            raise ValueError("Invalid rho. Should be an integer or float.")

        if not isinstance(random_seed, int):
            raise ValueError("The random seed must be an integer.")
        self.seed = random_seed

        if len(neighborhood) == 0:
            self.neighborhood = spatial.generate_offsets(2, 'manhattan') # Default neighborhood
        else:
            self.neighborhood = neighborhood

        if len(grid) == 0:
            self.grid = np.zeros(shape)
        else:
            self.grid = grid.copy()

        self.niche = np.zeros(shape).astype('int')

    ############################ Helper Functions ############################
    def initialize(self):
        """
        Initialize the grid with random states.

        Notes:
            The grid is initialized with random integers representing states.
            The shape of the grid is defined by the `shape` attribute.
            The random seed is set to ensure reproducibility.
        """
        np.random.seed(self.seed)
        init_state = np.random.randint(self.num_states, size=self.shape[0]*self.shape[1])
        self.grid = init_state.reshape(self.shape)

    def initialize3D(self):
        """
        Initialize the 3D grid with random states.

        Raises:
            ValueError: If the shape is not a 3-D tuple.

        Notes:
            The grid is initialized with random integers representing states.
            The shape of the grid is defined by the `shape` attribute.
            The random seed is set to ensure reproducibility.
        """
        if len(self.shape) != 3:
            raise ValueError("The shape must be a 3-D tuple.")

        np.random.seed(self.seed)
        init_state = np.random.randint(self.num_states, size=self.shape[0]*self.shape[1]*self.shape[2])
        self.grid = init_state.reshape(self.shape)

    def _wide_to_long(self):
        """
        Convert wide format to long format.
        """
        self.grid_long = self.grid.reshape(-1, 1)
        self.grid_long = pd.DataFrame(self.grid_long, columns=['state'])
        self.grid_long['row'] = np.repeat(range(self.shape[0]), self.shape[1])
        self.grid_long['col'] = np.tile(range(self.shape[1]), self.shape[0])
        self.grid_long['state'] = self.grid_long['state'].astype('int')

        self.meta = self.grid_long.copy()
        self.meta = self.meta[self.meta['state'] >= 0]
        cell_counts = self.meta['state'].value_counts()
        ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(cell_counts.index, 1)}
        self.meta['state_rank'] = self.meta['state'].map(ranked_cell_types)
        self.meta['state'] = self.meta['state'].astype('category')
        self.ncells = self.meta.shape[0]

    def update_seed(self, 
                    seed: int) -> None:
        """
        Update the random seed.

        Args:
            seed (int): The new random seed.
        
        Raises:
            ValueError: If the seed is not an integer or is negative.
        
        Example:
            >>> sim.update_seed(42)
        """
        if not isinstance(seed, int):
            raise ValueError("The random seed must be an integer.")
        if seed < 0:
            raise ValueError("The random seed must be a non-negative integer.")
        self.seed = seed

    def save(self, 
             path: str, 
             file_name: str = 'simspace.pkl') -> None:
        """
        Save the grid to a file using pickle. 

        Args:
            path (str): The path to save the grid.
            file_name (str): The name of the file to save the grid. Defaults to 'sim_space.pkl'.
        
        Raises:
            ValueError: If the path or file_name is not a string.

        Example:
            >>> sim.save('/path/to/save', 'my_simspace.pkl')
        """
        if not isinstance(path, str):
            raise ValueError("The path must be a string.")
        if not isinstance(file_name, str):
            raise ValueError("The file name must be a string.")

        if not path.endswith('/'):
            path += '/'
        with open(path + file_name, 'wb') as f:
            pickle.dump(self, f)

    ############################ Getter Functions ############################
    def print(self, type: str = 'wide') -> None:
        """
        Print the final grid.

        Args:
            type (str): The type of grid to print. Can be 'long' or 'wide'. Defaults to 'wide'.

        Raises:
            ValueError: If the type is not 'long' or 'wide'.

        Notes:
            This function prints the final grid in its current state.
            The grid is a 2-D or 3-D numpy array representing the simulated spatial data.
        """
        if type not in ['long', 'wide']:
            raise ValueError("Invalid type. Should be 'long' or 'wide'.")
        if type == 'long':
            print(self.grid_long)
        elif type == 'wide':
            print(self.grid)

    def plot_grid(self, 
                  figsize: tuple = (5, 5), 
                  dpi: int = 150) -> None:
        """
        Plot the final grid using seaborn.

        Args:
            figsize (tuple): The size of the figure. Defaults to (5, 5).
            dpi (int): The resolution of the figure. Defaults to 150.
        Raises:
            ValueError: If the grid is not a 2-D numpy array. 
        Notes:
            This function uses seaborn to create a heatmap of the grid.
            The grid should be a 2-D numpy array where each cell represents a state.
            The color palette used is cc.glasbey as default.
        Example:
            >>> sim.plot_grid(figsize=(10, 10), dpi=300)
        """
        if len(self.shape) != 2:
            raise ValueError("The grid must be a 2-D numpy array. Try plot3D() for 3-D grids.")
        cmap = sns.color_palette(cc.glasbey, self.num_states)
        tmp_grid = self.grid.copy()
        if np.min(tmp_grid) < 0:
            if tmp_grid.dtype.kind != 'f':
                tmp_grid = tmp_grid.astype('float')
                tmp_grid[tmp_grid < 0] = np.nan
        plt.figure(figsize=figsize, dpi=dpi)
        sns.heatmap(tmp_grid, cmap=cmap)
        plt.gca().set_aspect('equal')
        plt.show()

    def plot(self, 
             feature: str = 'state', 
             figsize: tuple = (8, 8), 
             dpi: int = 150, 
             size: int = 20, 
             title: str = None, 
             save_path: str = None, 
             legend: bool = True) -> None:
        """
        Plot SimSpace simulation results using seaborn.

        Args:
            feature (str): The feature to plot. Can be 'state', 'celltype', or any other feature in the meta or omics data.
            figsize (tuple): The size of the figure. Defaults to (8, 8).
            dpi (int): The resolution of the figure. Defaults to 150.
            size (int): The size of the points in the scatter plot. Defaults to 20.
            title (str, optional): The title of the plot.
            save_path (str, optional): The path to save the plot. If None, the plot will be shown instead.
            legend (bool): Whether to show the legend. Defaults to True.
        Raises:
            ValueError: If the feature is not found in the meta or omics data.
        Notes:
            This function uses seaborn to create a scatter plot of the specified feature.
            If the feature is 'state' or any features in self.meta.columns, it will plot the metadata with cells' simulated coordinates.
            If the feature is not found in self.meta, it will check if it exists in self.omics data.
        Example:
            >>> sim.plot(feature='state', figsize=(5, 5), dpi=300, size=14)
            >>> # For cell type visualization after omics fitting
            >>> sim.plot(feature='fitted_celltype')
            >>> # For a specific omics feature
            >>> sim.plot(feature='Gene_1')
        """
        if feature in self.meta.columns:
            # cmap = sns.color_palette(cc.glasbey, n_colors=self.meta[feature].nunique())
            cmap = sns.color_palette('tab20', n_colors=self.meta[feature].nunique())

            fig = plt.figure(figsize=figsize, dpi=dpi)
            ax = fig.add_subplot()
            ax.set_aspect('equal')
            ax = sns.scatterplot(data=self.meta, x='col', y='row', hue=feature, s=size, palette=cmap)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=feature)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            plt.tight_layout()
        else:
            if not hasattr(self, 'omics') or self.omics is None:
                raise ValueError(f"Feature '{feature}' not found in meta or omics data. Please check the feature name or if omics data exists.")
            if self.omics is not None and feature in self.omics.columns:
                plt_data = pd.DataFrame({
                    'row': self.meta['row'],
                    'col': self.meta['col'],
                    feature: self.omics[feature].values
                })
                fig = plt.figure(figsize=figsize, dpi=dpi)
                ax = fig.add_subplot()
                ax.set_aspect('equal')
                ax = sns.scatterplot(data=plt_data, x='col', y='row', hue=feature, s=size)
                if legend:
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=feature)
                else:
                    ax.get_legend().remove()
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                plt.tight_layout()
                if title is not None:
                    plt.title(title)
                    
            else:
                raise ValueError(f"Feature '{feature}' not found in meta or omics data. Please check the feature name.")

        if save_path is not None:
            plt.savefig(save_path)
        else:
            plt.show()

    def plot3D(self, 
               axis: str = 'z', 
               pos: int = 0, 
               figsize: tuple = (6, 6), 
               dpi: int = 150, 
               save_path: str = None) -> None:
        """
        Plot the 3D grid along a specified axis at a given position.

        Args:
            axis (str): The axis to plot. Can be 'x', 'y', or 'z'. Defaults to 'z'.
            pos (int): The position along the specified axis to plot. Defaults to 0.
            figsize (tuple): The size of the figure. Defaults to (6, 6).
            dpi (int): The resolution of the figure. Defaults to 150.
            save_path (str, optional): The path to save the plot. If None, the plot will be shown instead.
        Raises:
            ValueError: If the specified axis is invalid.
        Notes:
            This function uses seaborn to create a scatter plot of the specified axis at the given position.
            The plot will show the distribution of cell annotations in the 3D grid at that position.
            The color palette used is cc.glasbey as default.
        Example:
            >>> sim.plot3D(axis='z', pos=5)
        """
        if axis not in ['x', 'y', 'z']:
            raise ValueError("Invalid axis. Should be 'x', 'y', or 'z'.")
        left_axis = [ax for ax in ['x', 'y', 'z'] if ax != axis]
        plot_dat = self.meta[self.meta[axis].round(0) == pos]
        cmap = sns.color_palette(cc.glasbey, n_colors=self.num_states)

        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot()
        # if axis == 'z':
        ax.set_aspect('equal')
        ax = sns.scatterplot(data=plot_dat, x=left_axis[0], y=left_axis[1], hue='state', palette=cmap)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        # ax.invert_yaxis()
        ax.set_title(f"{axis} = {pos}")
        ax.set_xlabel(left_axis[0])
        ax.set_ylabel(left_axis[1])
        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path)
        else:
            plt.show()

    def plot_niche(self, figsize = (5, 5), dpi = 150):
        """
        Plot the niche in the simulation. 

        Args:
            figsize (tuple): The size of the figure. Defaults to (5, 5).
            dpi (int): The resolution of the figure. Defaults to 150.
        Raises:
            ValueError: If the niche is not a 2-D numpy array.
        Notes:
            This function uses seaborn to create a heatmap of the niche.
            The niche should be a 2-D numpy array where each cell represents a niche class.
            The color palette used is the default seaborn palette.
        Example:
            >>> sim.plot_niche(figsize=(5, 5), dpi=150)
        """
        if len(self.shape) != 2:
            raise ValueError("The niche must be a 2-D numpy array. Try plot_niche3D() for 3-D niches.")
        plt.figure(figsize=figsize, dpi=dpi)
        sns.heatmap(self.niche)
        plt.gca().set_aspect('equal')
        plt.gca().invert_yaxis()
        plt.show()
    
    ############################ Spatial Simulation ############################
    def create_niche(self, 
                     num_niches: int = 3, 
                     n_iter: int = 10, 
                     theta_niche: np.ndarray = None,
                     neighborhood: list = None,
                     ) -> None:
        """
        Apply MRF-based niche creation in a 2-D grid.

        Args:
            num_niches (int): The number of niches to create.
            n_iter (int): The number of iterations for the random niche creation.
            theta_niche (np.ndarray): The transition matrix for the niche-level MRF model.
            neighborhood (list): The list of neighboring cells, generated by spatial.generate_offsets. Default is spatial.generate_offsets(5, 'manhattan').
        Raises:
            ValueError: If the theta_niche is not a 2-D numpy array or if the number of niches exceeds the number of theta matrices provided.
        Notes:
            This function uses a Markov Random Field (MRF) approach to create spatial niches. 
            It initializes the niche grid with random integers representing niches and iteratively updates the niches based on their neighbors.
            If no theta_niche is provided, it defaults to a diagonal matrix with ones on the diagonal. This will result in a uniform distribution of niches.
        Examples:
            >>> sim.create_niche(num_niches=3, n_iter=6, theta_niche=theta_niche)
        """
        if neighborhood is None:
            neighborhood = spatial.generate_offsets(5, 'manhattan')

        if not isinstance(theta_niche, (np.ndarray, np.generic)):
            raise ValueError("Invalid theta_niche. Should be a 2-D numpy array.")
        if len(theta_niche.shape) != 2 or theta_niche.shape[0] != theta_niche.shape[1]:
            raise ValueError("Invalid theta_niche. Should be a square matrix.")
        np.random.seed(self.seed)
        self.niche = np.random.randint(num_niches, size=self.shape[0]*self.shape[1]).reshape(self.shape)

        if theta_niche is None:
            theta_niche = np.diag([1]*num_niches).astype(float)

        for _ in range(n_iter):
            for i in np.random.permutation(self.shape[0]):
                for j in np.random.permutation(self.shape[1]):
                    if self.grid[i, j] < 0: continue ## Skip the empty cells
                    neighbors = self.get_custom_neighbors(i, j, 
                                                          neighborhood=neighborhood) 
                    if len(neighbors) == 0: continue ## Skip if there are no neighbors

                    probabilities = self._conditional_probability(self.niche, theta_niche, neighbors) ## Compute the probabilities
                    self.niche[i, j] = np.random.choice(range(num_niches), p=probabilities) ## Update the cell

    def gibbs_sampler(self):
        """
        Perform Gibbs sampling to approximate the field grid in a 2-D grid.
        
        Notes:
            This function uses Gibbs sampling to update the grid cells based on their neighbors.
        Examples:
            >>> sim.gibbs_sampler()
        """
        np.random.seed(self.seed)
        for _ in range(self.num_iterations):
            for i in np.random.permutation(self.shape[0])[0:math.ceil(self.shape[0] * self.rho)]:
                for j in np.random.permutation(self.shape[1])[0:math.ceil(self.shape[1] * self.rho)]:
                    if self.grid[i, j] < 0: continue ## Skip the empty cells
                    neighbors = self.get_custom_neighbors(i, j, neighborhood= self.neighborhood) ## Get the neighbors
                    if len(neighbors) == 0: continue ## Skip if there are no neighbors

                    niche_class = self.niche[i, j]
                    cell_theta = self.theta[niche_class] ## Get the theta matrix for the cell
                    probabilities = self._conditional_probability(self.grid, cell_theta, neighbors) ## Compute the probabilities
                    self.grid[i, j] = np.random.choice(range(self.num_states), p=probabilities) ## Update the cell

    def _conditional_probability(self, 
                                 z: np.ndarray, 
                                 theta: np.ndarray, 
                                 neighbors: list) -> float:
        """
        Compute the conditional probability P(Z_ij = k | Z_-ij = z_-ij) given the neighbors.

        Args:
            z (numpy.ndarray): Current configuration.
            theta (numpy.ndarray): The theta matrix.
            neighbors (list): List of neighbor indices.

        Returns:
            float: The computed conditional probability.
        """
        hk_values = [self._compute_hk(z, k_prime, theta, neighbors) for k_prime in range(theta.shape[0])]
        if np.max(hk_values) > 0:
            hk_values = hk_values / np.max(hk_values) * self.phi
        exp_values = np.exp(hk_values)
        res = exp_values / np.sum(exp_values)
        return res
    
    def _compute_hk(self, 
                    z: np.ndarray, 
                    k: int, 
                    theta: np.ndarray, 
                    neighbors: list) -> float:
        """
        Compute the value of hk based on the given parameters.

        Args:
            z (numpy.ndarray): Current configuration.
            k (int): State to compute the probability for.
            theta (numpy.ndarray): The theta matrix.
            neighbors (list): List of neighbor indices.

        Returns:
            float: The computed value of hk.
        """
        hk = 0
        for ni, nj in neighbors:
            hk += theta[k, z[ni, nj]]
        return hk

    def get_custom_neighbors(self, 
                             i: int, 
                             j: int, 
                             neighborhood: list) -> list:
        """
        Get the custom neighbors for a pixel (i, j) in a 2-D grid based on the specified neighborhood.

        Args:
            i (int): Row index of the current pixel.
            j (int): Column index of the current pixel.
            neighborhood (list): List of neighbor offsets relative to the current pixel.

        Returns:
            list: List of custom neighbor indices.

        Examples:
            >>> sim.get_custom_neighbors(5, 5, neighborhood=spatial.generate_offsets(2, 'manhattan'))
            [(4, 5), (6, 5), (5, 4), (5, 6)]
        """
        neighbors = []
        for offset in neighborhood:
            ni = i + offset[0]
            nj = j + offset[1]
            if 0 <= ni < self.shape[0] and 0 <= nj < self.shape[1]:
                if self.grid[ni, nj] >= 0: ## Skip the empty cells
                    neighbors.append((ni, nj)) ## Add the neighbor
        return neighbors

    def density_sampler(self, 
                        threshold: float | list) -> None:
        """
        Cell density sampler. This sampler will randomly select cells based on the given threshold either uniformly or per state.

        Args:
            threshold (float or list): The threshold for the density sampler. Should be a list which must match either the number of states or the number of niches in the grid. If only one state or niche is available, it can be a float.

        Raises:
            ValueError: If the threshold is not a float or list.
        
        Notes:
            If a list of thresholds is provided, it must match the number of states in the grid.
            If a single float is provided, it will be applied uniformly across all states.
        """

        if isinstance(threshold, float):
            threshold = [threshold]
        elif isinstance(threshold, list):
            pass
        elif isinstance(threshold, np.ndarray):
            threshold = threshold.tolist()
        else:
            raise ValueError("Invalid threshold.")
        
        density_grid = np.random.uniform(0, 1, self.shape)
        self.density_bool = np.zeros_like(density_grid)

        ## If the threshold is provided for each state, then apply the threshold for each state
        if len(threshold) == self.num_states:
            # print("Density sampler applied for each state.")
            for i in range(len(threshold)):
                self.density_bool[self.grid == i] = density_grid[self.grid == i] < threshold[i]
        else:   
            ## If the threshold is not provided for each state, then apply the same threshold for all states in one niche
            # print("Density sampler applied for each niche.")
            if len(threshold) != self.niche.max() + 1:
                raise ValueError("The number of thresholds must match the number of niches in the grid.")
            for i in range(len(threshold)):
                self.density_bool[self.niche == i] = density_grid[self.niche == i] < threshold[i]

        self.grid[self.density_bool == 0] = -1

    def perturbation(self, step: float | int) -> None:
        """
        Perturb the coordinates of the grid. It adds Gaussian noise to the coordinates of the grid, so the cells are randomly displaced instead of being perfectly aligned to the grid.

        Args:
            step (float | int): The standard deviation of the Gaussian noise to be added to the coordinates.
        Raises:
            ValueError: If the step is not a positive number.
        Examples:
            >>> sim.perturbation(step = 0.2)
        """

        if not isinstance(step, (float, int)):
            raise ValueError("Step must be a float or int.")
        if step < 0:
            raise ValueError("Step must be a positive number.")

        np.random.seed(self.seed)
        self._wide_to_long()
        self.meta['row'] = self.meta['row'] + np.random.normal(0, step, self.meta.shape[0])
        self.meta['col'] = self.meta['col'] + np.random.normal(0, step, self.meta.shape[0])
        self.meta.reset_index(drop=True, inplace=True)


    def _add_ellipse(self, niche_idx, center, radius_x, radius_y, angle=0, overlap=False):
        """
        Add a niche to the grid.

        Args:
            center (tuple): The center of the niche.
            radius_x (int): The radius of the niche along the x-axis.
            radius_y (int): The radius of the niche along the y-axis.
            angle (int): The rotation angle of the niche.
            overlap (bool): Whether to allow overlap with existing niches. If True, the niche will be added only if it does not overlap with existing niches.
        """
        ellipse_mask = niche.create_ellipse(self.shape, center, radius_x, radius_y, angle)
        if overlap:
            self.niche = self.niche * (1 - ellipse_mask) + ellipse_mask * niche_idx
        else:
            self.niche += ellipse_mask * niche_idx

    def _add_rectangle(self, center, length, width):
        """
        Add a rectangle to the grid.

        Args:
            center (tuple): The center of the vessel.
            length (int): The length of the vessel.
            width (int): The width of the vessel.
        """
        rectangle_mask = niche.create_rectangle(self.shape, center, length, width)
        self.niche += rectangle_mask * max(self.niche.flatten() + 1)

    def manual_niche(self, 
                     pattern: dict = {}) -> None:
        """
        Manually create niches with given patterns in the grid.

        Args:
            pattern (dict): A dictionary where keys are niche indices and values are tuples containing the shape ('ellipse' or 'rectangle') and parameters for the shape.
                            For 'ellipse': (center_x, center_y, radius_x, radius_y, angle)
                            For 'rectangle': (center_x, center_y, length, width)
        Raises:
            ValueError: If the pattern is not a valid dictionary or if the shapes are not 'ellipse' or 'rectangle'.
        """
        for key in pattern.keys():
            if pattern[key][0] == 'ellipse':
                self._add_ellipse(
                    list(pattern).index(key) + 1, # niche index
                    (pattern[key][1][0], pattern[key][1][1]), # center
                    pattern[key][1][2], pattern[key][1][3], # radius_x, radius_y
                    pattern[key][1][4], # angle
                    ) 
            elif pattern[key][0] == 'rectangle':
                self._add_rectangle(
                    list(pattern).index(key) + 1,
                    (pattern[key][1][0], pattern[key][1][1]), 
                    pattern[key][1][2], 
                    pattern[key][1][3]
                    )
            else:
                raise ValueError("Invalid pattern. Should be either 'ellipse' or 'rectangle'.")

        if len(self.theta) < max(self.niche.flatten()):
            raise ValueError("The number of niches exceeds the number of theta matrices provided. Prohaps there are intercections between niches.")
        self.niche = self.niche.astype('int')
     

    ############################ 3D Spatial Simulation #########################

    def _wide_to_long3D(self):
        """
        Convert wide format to long format for 3-D grids.
        """
        self.grid_long = self.grid.reshape(-1, 1)
        self.grid_long = pd.DataFrame(self.grid_long, columns=['state'])
        self.grid_long['x'] = np.repeat(range(self.shape[0]), self.shape[1] * self.shape[2])
        self.grid_long['y'] = np.tile(np.repeat(range(self.shape[1]), self.shape[2]), self.shape[0])
        self.grid_long['z'] = np.tile(range(self.shape[2]), self.shape[0] * self.shape[1])
        self.grid_long['state'] = self.grid_long['state'].astype('int')

        self.meta = self.grid_long.copy()
        self.meta = self.meta[self.meta['state'] >= 0]
        self.meta['state'] = self.meta['state'].astype('category')

    def create_niche3D(
            self, 
            num_niches: int = 3, 
            n_iter: int = 10, 
            theta_niche: np.ndarray = None,
            neighborhood: np.ndarray = None) -> None:
        """
        Apply MRF-based niche creation in a 3-D grid.

        Args:
            num_niches (int): The number of niches to create.
            n_iter (int): The number of iterations for the random niche creation.
            theta_niche (np.ndarray): The transition matrix for the niche-level MRF model.
            neighborhood (list): The list of neighboring cells, generated by spatial.generate_offsets. Default is spatial.generate_offsets(5, 'manhattan').
        Raises:
            ValueError: If the self.shape is not a 3-D numpy array or if the theta_niche is not a 2-D numpy array.
        Examples:
            >>> sim.create_niche3D(num_niches=3, n_iter=6, theta_niche=theta_niche)
        """
        if neighborhood is None:
            neighborhood = spatial.generate_offsets3D(5, 'manhattan')

        if len(self.shape) != 3:
            raise ValueError("The grid must be a 3-D numpy array. Try create_niche() for 2-D grids.")
        if not isinstance(theta_niche, (np.ndarray, np.generic)):
            raise ValueError("Invalid theta_niche. Should be a 2-D numpy array.")
        if len(theta_niche.shape) != 2 or theta_niche.shape[0] != theta_niche.shape[1]:
            raise ValueError("Invalid theta_niche. Should be a square matrix.")
        
        np.random.seed(self.seed)
        self.niche = np.random.randint(num_niches, size=self.shape[0]*self.shape[1]*self.shape[2]).reshape(self.shape)

        if theta_niche is None:
            theta_niche = np.diag([1]*num_niches).astype(float)

        for _ in range(n_iter):
            for i in np.random.permutation(self.shape[0]):
                for j in np.random.permutation(self.shape[1]):
                    for k in np.random.permutation(self.shape[2]):
                        if self.grid[i, j, k] < 0: continue ## Skip the empty cells
                        neighbors = self.get_custom_neighbors3D(i, j, k,
                                                                neighborhood=neighborhood) 
                        if len(neighbors) == 0: continue ## Skip if there are no neighbors

                        probabilities = self._conditional_probability3D(self.niche, theta_niche, neighbors) ## Compute the probabilities
                        self.niche[i, j, k] = np.random.choice(range(num_niches), p=probabilities) ## Update the cell

    def gibbs_sampler3D(self):
        """
        Perform Gibbs sampling to approximate the field grid in a 3-D grid.

        Notes:
            This function uses Gibbs sampling to update the grid cells based on their neighbors.
        Examples:
            >>> sim = SimSpace(shape=(10, 30, 30), num_states=5, theta=np.random.rand(5, 5), num_iterations=10, seed=42)
            >>> sim.gibbs_sampler3D()
        """
        np.random.seed(self.seed)
        for _ in range(self.num_iterations):
            for i in np.random.permutation(self.shape[0])[0:math.ceil(self.shape[0] * self.rho)]:
                for j in np.random.permutation(self.shape[1])[0:math.ceil(self.shape[1] * self.rho)]:
                    for k in np.random.permutation(self.shape[2])[0:math.ceil(self.shape[2] * self.rho)]:
                        if self.grid[i, j, k] < 0: continue ## Skip the empty cells
                        neighbors = self.get_custom_neighbors3D(i, j, k, neighborhood=self.neighborhood) ## Get the neighbors
                        if len(neighbors) == 0: continue ## Skip if there are no neighbors

                        niche_class = self.niche[i, j, k]
                        cell_theta = self.theta[niche_class] ## Get the theta matrix for the cell
                        probabilities = self._conditional_probability3D(self.grid, cell_theta, neighbors) ## Compute the probabilities
                        self.grid[i, j, k] = np.random.choice(range(self.num_states), p=probabilities) ## Update the cell

    def _compute_hk3D(self, 
                      z: np.ndarray, 
                      k: int, 
                      theta: np.ndarray, 
                      neighbors: list) -> float:
        """
        Compute the value of hk based on the given parameters.

        Args:
            z (numpy.ndarray): Current configuration.
            k (int): State to compute the probability for.
            theta (numpy.ndarray): The theta matrix.
            neighbors (list): List of neighbor indices.

        Returns:
            float: The computed value of hk.
        """
        hk = 0
        for ni, nj, nk in neighbors:
            hk += theta[k, z[ni, nj, nk]]
        return hk

    def _conditional_probability3D(self, 
                                   z: np.ndarray, 
                                   theta: np.ndarray, 
                                   neighbors: list) -> np.ndarray:
        """
        Compute the conditional probability P(Z_ij = k | Z_-ij = z_-ij) given the neighbors.

        Args:
            z (numpy.ndarray): Current configuration.
            neighbors (list): List of neighbor indices.

        Returns:
            float: The computed conditional probability.
        """
        hk_values = [self._compute_hk3D(z, k_prime, theta, neighbors) for k_prime in range(theta.shape[0])]
        if np.max(hk_values) > 0:
            hk_values = hk_values / np.max(hk_values) * self.phi
        exp_values = np.exp(hk_values)
        res = exp_values / np.sum(exp_values)
        return res

    def get_custom_neighbors3D(self, 
                               i: int, 
                               j: int, 
                               k: int, 
                               neighborhood: list) -> list:
        """
        Get the custom neighbors for a pixel (i, j, k) in a 3-D grid based on the specified neighborhood.

        Args:
            i (int): Row (X) index of the current pixel.
            j (int): Column (Y) index of the current pixel.
            k (int): Z index of the current pixel.
            neighborhood (list): List of neighbor offsets relative to the current pixel.

        Returns:
            list: List of custom neighbor indices.
        """
        neighbors = []
        for offset in neighborhood:
            ni = i + offset[0]
            nj = j + offset[1]
            nk = k + offset[2]
            if 0 <= ni < self.shape[0] and 0 <= nj < self.shape[1] and 0 <= nk < self.shape[2]:
                if self.grid[ni, nj, nk] >= 0: ## Skip the empty cells
                    neighbors.append((ni, nj, nk)) ## Add the neighbor
        return neighbors
    
    def perturbation3D(self, 
                       step: float | int) -> None:
        """
        Perturb the coordinates of the grid. It adds Gaussian noise to the coordinates of the grid, so the cells are 
        randomly displaced instead of being perfectly aligned to the grid.

        Args:
            step (float | int): The standard deviation of the Gaussian noise to be added to the coordinates.
        Raises:
            ValueError: If the step is not a positive number.
        Examples:
            >>> sim.perturbation3D(step = 0.2)
        """
        if not isinstance(step, (float, int)):
            raise ValueError("Step must be a float or int.")
        if step < 0:
            raise ValueError("Step must be a positive number.")
        
        np.random.seed(self.seed)
        self._wide_to_long3D()
        self.meta['x'] = self.meta['x'] + np.random.normal(0, step, self.meta.shape[0])
        self.meta['y'] = self.meta['y'] + np.random.normal(0, step, self.meta.shape[0])
        self.meta['z'] = self.meta['z'] + np.random.normal(0, step, self.meta.shape[0])
        self.meta.reset_index(drop=True, inplace=True)

    ############################ Spatial Analysis      #########################
    def moran_I(self):
        """
        Calculate the global Moran's I. Results are stored in self.moran_I_value.
        """
        self.moran_I_value = spatial.integrate_morans_I(self.meta["state"], self.meta[['row', 'col']], self.meta["state"].unique())

    ############################ Expression Simulation #########################
    def create_omics(
            self, 
            bg_ratio = 0.2, 
            n_genes = 1000, 
            bg_param: typing.Tuple[float, float] = (1, 1),
            marker_param: typing.Tuple[float, float] = (5, 2),
            lr_ratio: float = 0.5,
            spatial=True,
            k_neighors: int = 20,
            spatial_effect: float = 3,
            se_threshold: float = 1.5,) -> None:
        """
        Create the reference-free omics data using Gamma-Poisson distribution.

        Args:
            bg_ratio (float): The ratio of background genes to total genes. Defaults to 0.2.
            n_genes (int): The total number of genes to simulate. Defaults to 1000.
            bg_param (tuple): The parameters for the background gene distribution (shape, scale). Defaults to (1, 1).
            marker_param (tuple): The parameters for the marker gene distribution (shape, scale). Defaults to (5, 2).
            lr_ratio (float): The ratio of ligand-receptor pairs to total genes. Defaults to 0.5.
            spatial (bool): Whether to simulate spatial omics data. Defaults to True.
            k_neighors (int): The number of neighbors to consider for spatial omics.
            spatial_effect (float): The spatial effect parameter for spatial omics. Defaults to 3.
            se_threshold (float): The threshold for spatial effect. Defaults to 1.5.
        Raises:
            ValueError: If the bg_ratio is not between 0 and 1, or if the n_genes is not a positive integer.
        Notes:
            This function simulates omics data based on the provided parameters.
            It uses the `omics` module to generate the omics data and metadata.
            The simulated omics data can be used for further analysis or visualization.
        Examples:
            >>> sim.create_omics(bg_ratio=0.2, n_genes=500)
        """
        self.gene_meta=omics.simOmicsMeta(
            meta=self.meta, 
            bg_ratio=bg_ratio,
            n_genes=n_genes, 
            bg_param=bg_param,
            marker_param=marker_param,
            lr_ratio=lr_ratio,
            seed=self.seed,
            )
        if bg_ratio < 0 or bg_ratio > 1:
            raise ValueError("bg_ratio must be between 0 and 1.")
        if not isinstance(n_genes, int) or n_genes <= 0:
            raise ValueError("n_genes must be a positive integer.")
    
        omics_res = omics.simOmics(
            omics_meta=self.gene_meta, 
            meta=self.meta, 
            seed=self.seed,
            )
        if spatial:
            self.omics = omics.simSpatialOmics(
                gene_data=omics_res, 
                gene_meta=self.gene_meta, 
                cell_meta=self.meta, 
                k_neighors=k_neighors, 
                spatial_effect=spatial_effect,
                se_threshold=se_threshold,
                seed=1)
        else:
            self.omics = omics_res
    
    def fit_scdesign(
            self,
            ref_count_path,
            ref_meta_path,
            group_col,
            spatial_x,
            spatial_y,
            seed = 0,
            isreturn=False,
            ):
        """
        Fit the scdesign model using the reference dataset. 

        Args:
            ref_count_path (str): Path to the reference count matrix.
            ref_meta_path (str): Path to the reference metadata.
            group_col (str): Column name in the metadata for grouping in the simulation.
            spatial_x (str): Column name in the metadata for the x-coordinate of spatial data.
            spatial_y (str): Column name in the metadata for the y-coordinate of spatial data.
            seed (int): Random seed for reproducibility. Defaults to 0.
            isreturn (bool): If True, return the simulated count matrix. Defaults to False.

        Raises:
            ValueError: If the reference dataset is too small or if scdesign_fit returns None for sim_meta or sim_count.
        
        Notes: 
            This function uses the `omics` module to fit the scdesign model based on the provided reference dataset.
            It ranks the cell types based on their frequency in the reference metadata and maps them to the simulation metadata.
            The simulated count matrix is stored in `self.omics` and can be returned if `isreturn` is set to True.

        Example:
            >>> sim.fit_scdesign(
            ... ref_count_path='path/to/ref_count.csv',
            ... ref_meta_path='path/to/ref_meta.csv',
            ... group_col='celltype',   # Should match the column in ref_meta
            ... spatial_x='x_coord',    # Column in ref_meta for x-coordinate
            ... spatial_y='y_coord',    # Column in ref_meta for y-coordinate
            ... seed=42)
        """
        ref_meta = pd.read_csv(ref_meta_path, index_col=0)
        # Rank the cell types based on their frequency
        ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(ref_meta[group_col].value_counts().index, 1)}
        # Map the 'celltype' to their frequency rank
        ref_meta['state_rank'] = ref_meta[group_col].map(ranked_cell_types)
        # Create a map for ref_meta['state_rank'] and ref_meta['celltype']
        state_rank_to_celltype = ref_meta.set_index('state_rank')[group_col].to_dict()
        # Use this map to create a 'celltype' column in self.meta
        self.meta['fitted_celltype'] = self.meta['state_rank'].map(state_rank_to_celltype)

        sim_count, _ = omics.scdesign_fit(
            ref_count_path,
            ref_meta_path,
            group_col,
            spatial_x,
            spatial_y,
            self.meta,
            seed=seed,
        )
        sim_count = sim_count.T
        
        self.omics = sim_count.copy()
        if isreturn:
            return sim_count
    
    def fit_srtsim(
            self,
            ref_count_path,
            ref_meta_path,
            group_col,
            spatial_x,
            spatial_y,
            seed = 0,
            isreturn=False,
            ):
        """
        Fit the scdesign model using the reference dataset.

        Args:
            ref_count_path (str): Path to the reference count matrix.
            ref_meta_path (str): Path to the reference metadata.
            group_col (str): Column name in the metadata for grouping in the simulation.
            spatial_x (str): Column name in the metadata for the x-coordinate of spatial data.
            spatial_y (str): Column name in the metadata for the y-coordinate of spatial data.
            seed (int): Random seed for reproducibility. Defaults to 0.
            isreturn (bool): If True, return the simulated count matrix. Defaults to False.
        
        Raises:
            ValueError: If the reference dataset is too small or if scdesign_fit returns None for sim_meta or sim_count.
        
        Notes: 
            This function uses the `omics` module to fit the scdesign model based on the provided reference dataset.
            It ranks the cell types based on their frequency in the reference metadata and maps them to their corresponding state ranks.
        """
        ref_meta = pd.read_csv(ref_meta_path, index_col=0)

        ranked_cell_types = {cell_type: rank for rank, cell_type in enumerate(ref_meta[group_col].value_counts().index, 1)}
        ref_meta['state_rank'] = ref_meta[group_col].map(ranked_cell_types)
        ref_n_cells = ref_meta.shape[0]
        n_cells_ratio = ref_n_cells / self.ncells
        if n_cells_ratio < 0.3:
            raise ValueError("The number of cells in the reference dataset is too small. Please provide a larger reference dataset.")
        elif n_cells_ratio < 1:
            n_rep = int(1 / n_cells_ratio) + 2
        else:
            n_rep = 2

        sim_count, sim_meta = omics.srtsim_fit(
            ref_count_path,
            ref_meta_path,
            group_col,
            spatial_x,
            spatial_y,
            seed=seed,
            n_rep=n_rep,
        )
        sim_meta.rename(columns={'label': group_col}, inplace=True)
        sim_meta['state_rank'] = sim_meta[group_col].map(ranked_cell_types)
        sim_count = sim_count.T

        if sim_meta is None or sim_count is None:
            raise ValueError("scdesign_fit returned None for sim_meta or sim_count. Please check the input data and parameters.")
        
        sim_count_new = pd.DataFrame(index=self.meta.index, columns=sim_count.columns)
        for state_rank in self.meta['state_rank'].unique():
            # Filter cells with the same state_rank in sim_meta
            matching_cells = sim_meta[sim_meta['state_rank'] == state_rank].index
            # Randomly sample rows from sim_count corresponding to the matching cells
            num_samples = self.meta[self.meta['state_rank'] == state_rank].shape[0]
            sampled_rows = sim_count.loc[matching_cells].sample(
                n=num_samples,
                replace=True,
                random_state=self.seed
            ).reset_index(drop=True)
            # Ensure lengths match before assignment
            if len(sampled_rows) == num_samples:
                sim_count_new.loc[self.meta[self.meta['state_rank'] == state_rank].index] = sampled_rows.values
            else:
                raise ValueError(f"Mismatch in lengths: {len(sampled_rows)} vs {num_samples}")
        sim_count_new.reset_index(drop=True, inplace=True)
        
        self.omics = sim_count_new.copy()
        if isreturn:
            return sim_count_new

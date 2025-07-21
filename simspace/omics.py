import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import typing
import subprocess
import os
from scipy.spatial import KDTree

def simOmicsMeta(
        meta: pd.DataFrame,
        n_genes: int = 1000,
        bg_ratio: float = 0.5,
        bg_param: typing.Tuple[float, float] = (1, 1),
        marker_param: typing.Tuple[float, float] = (5, 2),
        lr_ratio: float = 0.5,
        seed: int = 0,
        ) -> pd.DataFrame:
    """
    Simulate the metadata of omics data
    Args:
        meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        n_genes (int): The number of genes to simulate.
        bg_ratio (float): The ratio of background genes (non-marker genes) to the total number of genes.
        bg_param (tuple): The parameters for the gamma distribution to simulate the mean expression level of background genes.
        marker_param (tuple): The parameters for the gamma distribution to  simulate the mean expression level of marker genes.
        lr_ratio (float): The ratio of ligand-receptor pairs to the total number of marker genes.
        seed (int): The random seed for reproducibility.
    Returns:
        pd.DataFrame: A DataFrame containing the simulated metadata of the omics data, with columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
            - Type_{cell_type}: The mean expression level of the gene in the corresponding cell type.
    
    Raises:
        ValueError: If the background gene ratio or ligand-receptor pair ratio is not between 0 and 1, or if the metadata does not contain the 'state' column.
    
    Examples:
        >>> meta = pd.DataFrame({'state': ['A', 'B', 'C']})
        >>> omics_meta = simOmicsMeta(meta, n_genes=100, bg_ratio=0.3, lr_ratio=0.2, seed=42)
        >>> print(omics_meta.head())
    """

    np.random.seed(seed)

    if bg_ratio < 0 or bg_ratio > 1:
        raise ValueError('The background gene ratio should be between 0 and 1.')
    if lr_ratio < 0 or lr_ratio > 1:
        raise ValueError('The ligand-receptor pair ratio should be between 0 and 1.')
    if 'state' not in meta.columns:
        raise ValueError('The metadata does not contain the "state" column. Please check the column names.')
    
    type_arr = meta['state'].unique()
    type_arr = sorted(type_arr)
    n_types = len(type_arr)

    # Randomly assign a cell type to each gene
    marker_arr = np.random.choice(type_arr, n_genes)

    # Randomly set a part of genes as background genes (non-marker genes)
    n_bg = int(n_genes * bg_ratio)
    bg_indices = np.random.choice(n_genes, n_bg, replace=False)
    marker_arr[bg_indices] = -1

    # Randomly generates the ligand-receptor pair among genes that are not background genes
    ligand_receptor_arr = np.full(n_genes, -1)
    marker_indices = np.where(marker_arr != -1)[0]
    n_pairs = int(len(marker_indices) * lr_ratio)
    paired_indices = np.random.choice(marker_indices, n_pairs, replace=False)

    for i in paired_indices:
        if ligand_receptor_arr[i] != -1:
            continue
        possible_pairs = np.where((marker_arr != -1) & (ligand_receptor_arr == -1))[0]
        possible_pairs = possible_pairs[possible_pairs != i]
        if len(possible_pairs) > 0:
            tmp = np.random.choice(possible_pairs)
            ligand_receptor_arr[i] = tmp
            ligand_receptor_arr[tmp] = i

    # Create the omics metadata
    omics_meta = pd.DataFrame({
        "GeneID": np.arange(n_genes),
        "Marker": marker_arr,
        "LRindex": ligand_receptor_arr,
        })
    
    # Simulate the mean expression level of each gene
    for celltype in type_arr:
        omics_meta[f'Type_{celltype}'] = np.zeros(n_genes)

    for i in range(n_genes):
        gene_lam = np.random.gamma(bg_param[0], bg_param[1], size=n_types)
        omics_meta.iloc[i, 3:] = gene_lam

        if omics_meta.loc[i, 'Marker'] != -1:
            omics_meta.loc[i, f'Type_{omics_meta.loc[i, "Marker"]}'] = np.random.gamma(marker_param[0], marker_param[1], 1)

    return omics_meta

def simOmics(omics_meta: pd.DataFrame, 
             meta: pd.DataFrame,
             seed: int = 0,
             ) -> pd.DataFrame:
    """
    Simulate the omics data based on the metadata and cell metadata. 

    Args:
        omics_meta (pd.DataFrame): The metadata of the omics data, which should contain the following columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
            - Type_{cell_type}: The mean expression level of the gene in the corresponding cell type.
        meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        seed (int): The random seed for reproducibility.

    Returns:
        pd.DataFrame: A DataFrame containing the simulated omics data, with columns:
            - Gene_{gene_id}: The expression level of the gene in each cell.

    Raises:
        ValueError: If the metadata does not contain the 'state' column.
        TypeError: If the input data is not a DataFrame.

    Examples:
        >>> omics_meta = pd.DataFrame({
        ...     'GeneID': [0, 1, 2],
        ...     'Marker': ['A', 'B', -1],
        ...     'LRindex': [-1, 0, -1],
        ...     'Type_A': [10, 20, 0],
        ...     'Type_B': [0, 30, 0]
        ... })
        >>> meta = pd.DataFrame({'state': ['A', 'B']})
        >>> omics_data = simOmics(omics_meta, meta, seed=42)
        >>> print(omics_data.head())
    """
    if 'state' not in meta.columns:
        raise ValueError('The metadata does not contain the "state" column. Please check the column names.')
    if not isinstance(omics_meta, pd.DataFrame) or not isinstance(meta, pd.DataFrame):
        raise TypeError('The input data must be a pandas DataFrame.')
    
    np.random.seed(seed)
    
    n_genes = omics_meta.shape[0]
    n_cells = len(meta)

    # Simulate the expression level of each gene in each cell
    # The expression level is generated from a poisson distribution
    omics_data = np.zeros((n_cells, n_genes))
    for i in range(n_cells):
        cell_type = meta.loc[i, 'state']
        for j in range(n_genes):
            expr = np.random.poisson(omics_meta.loc[j, f'Type_{cell_type}'])
            omics_data[i, j] = expr

    omics_data = pd.DataFrame(omics_data, columns=[f'Gene_{gene_id}' for gene_id in omics_meta['GeneID']], index=meta.index)

    return omics_data

def simSpatialOmics(gene_data: pd.DataFrame,
                    gene_meta: pd.DataFrame,
                    cell_meta: pd.DataFrame,
                    k_neighors: int = 10,
                    spatial_effect: float = 1.0,
                    se_threshold: float = 1.5,
                    seed: int = 0,
                    ) -> pd.DataFrame:
    """
    Simulate the spatial omics data

    Args:
        gene_data (pd.DataFrame): The gene expression data, with cells as rows and genes as columns.
        gene_meta (pd.DataFrame): The metadata of the genes, which should contain the following columns:
            - GeneID: The index of the gene.
            - Marker: The cell type of the gene, -1 for background genes.
            - LRindex: The index of the ligand-receptor pair, -1 if the gene is not a ligand or receptor.
        cell_meta (pd.DataFrame): The metadata of the cells, which should contain a column named 'state' representing the cell types.
        k_neighors (int): The number of nearest neighbors to consider for spatial effects.
        spatial_effect (float): The factor by which to increase or decrease the expression level based on spatial effects.
        se_threshold (float): The threshold for spatial effect application.
        seed (int): The random seed for reproducibility.

    Returns:
        pd.DataFrame: A DataFrame containing the simulated spatial omics data, with cells as rows and genes as columns.

    Raises:
        ValueError: If the spatial effect is not greater than 1, or if the cell metadata does not contain the coordinates or cell types.
    
    Examples:
        >>> gene_data = pd.DataFrame({
        ...     'Gene_0': [10, 20, 30],
        ...     'Gene_1': [5, 15, 25],
        ...     'Gene_2': [0, 10, 20]
        ... }, index=['cell_1', 'cell_2', 'cell_3'])
        >>> gene_meta = pd.DataFrame({
        ...     'GeneID': [0, 1, 2],
        ...     'Marker': ['A', 'B', -1],
        ...     'LRindex': [-1, 0, -1]
        ... })
        >>> cell_meta = pd.DataFrame({
        ...     'state': ['A', 'B', 'A'],
        ...     'row': [0, 1, 0],
        ...     'col': [0, 1, 2]
        ... })
        >>> spatial_omics = simSpatialOmics(gene_data, gene_meta, cell_meta, k_neighors=2, spatial_effect=2, seed=42)
        >>> print(spatial_omics.head())
    """

    assert spatial_effect > 0, 'The spatial effect should be greater than 0.'

    if 'col' not in cell_meta.columns or 'row' not in cell_meta.columns:
        raise ValueError('The cell metadata does not contain the coordinates. Please check the column names.')
    
    if 'state' not in cell_meta.columns:
        raise ValueError('The cell metadata does not contain the cell types. Please check the column names.')
    
    output_data = gene_data.copy()

    np.random.seed(seed)
    n_gene = gene_data.shape[1]
    n_cell = gene_data.shape[0]

    # Extract coordinates
    coordinates = cell_meta[['row', 'col']].values

    # Build a KDTree for efficient neighbor search
    tree = KDTree(coordinates)

    # Find the top 10 nearest neighbors for each cell
    _, neighbors = tree.query(coordinates, k=k_neighors+1)  # k+1 because the closest neighbor is the point itself

    # Create a dictionary to store neighbors for each cell, excluding the cell itself
    cell_neighbors = {k: neighbors[k][1:] for k in range(n_cell)}

    for i in range(n_gene):
        LRindex = gene_meta.loc[i, 'LRindex']
        Marker = gene_meta.loc[i, 'Marker']
        if Marker == -1:
            continue
        if Marker == LRindex:
            continue
        if LRindex != -1:
            mean_expr = gene_data.iloc[:, LRindex].mean()
            for j in range(n_cell):
                if cell_meta.state[j] == gene_meta.loc[i, 'Marker']:
                    neighbor_indices = cell_neighbors[j]
                    neighbor_expr = gene_data.iloc[neighbor_indices, LRindex].values.sum()
                    neighbor_ratio = neighbor_expr / mean_expr / k_neighors
                    if neighbor_ratio > se_threshold:
                        output_data.iloc[j, i] *= spatial_effect * neighbor_ratio
                    else:
                        output_data.iloc[j, i] /= spatial_effect
                
    return output_data

def run_splatter(new_meta, 
                 ngene = 1000, 
                 r_script_path=None):
    """
    Run the splatter simulation to generate synthetic single-cell RNA-seq data.

    Args:
        new_meta (pd.DataFrame): A DataFrame containing simulated spatial metadata for omics simulation, which is derived from the .meta of the simspace object.
        ngene (int): The number of genes to simulate.
        r_script_path (str): The path to the R script that performs the splatter simulation. Default is None, which uses the script of simspace package.
    
    Returns:
        tuple: A tuple containing two DataFrames:
            - splatter_data: The simulated gene expression data.
            - splatter_meta: The metadata of the simulated cells.
    
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    
    Examples:
        >>> splatter_data, splatter_meta = run_splatter(sim.meta, ngene=1000) # sim is simulated simspace object
        >>> print(splatter_data.head())
        >>> print(splatter_meta.head())
    """

    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/splatter.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")
    new_meta.to_csv("./tmp/cells_meta.csv", index=False)
        
    result = subprocess.run(
        ["Rscript", r_script_path, "./tmp/cells_meta.csv", ngene],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("Splatter simulation complete.")
        # Load the simulated data
        splatter_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        splatter_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/cells_meta.csv")
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")

        return splatter_data, splatter_meta

def splatter_fit(count_path, 
                 meta_path, 
                 group_col, 
                 n_cells = 2000,
                 r_script_path=None):
    """
    Fit the splatter model to the given reference data and metadata to simulate new omics data.

    Args:
        count_path (str): The path to the count matrix file.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information.
        n_cells (int): The number of cells to simulate. Should match the number of cells in the spatial simulation results.
        r_script_path (str): The path to the R script that performs the splatter fitting. Default is None, which uses the script of simspace package.
    
    Returns:
        tuple: A tuple containing two DataFrames:
            - splatter_data: The simulated gene expression data.
            - splatter_meta: The metadata of the simulated cells.
    
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    
    Examples:
        >>> splatter_data, splatter_meta = splatter_fit("path/to/count.csv",
        ...                                              "path/to/meta.csv",
        ...                                              "group_column",
        ...                                              n_cells=2000)
        >>> print(splatter_data.head())
        >>> print(splatter_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/splatter_fit.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")

    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    result = subprocess.run(
        ["Rscript", r_script_path, meta_path, count_path, group_col, str(n_cells)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("Splatter fit complete.")
        splatter_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        splatter_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return splatter_data, splatter_meta

def scdesign_fit(count_path, 
                 meta_path, 
                 group_col, 
                 spatial_x,
                 spatial_y,
                 new_meta,
                 seed = 0,
                 r_script_path=None):
    """
    Fit the scdesign model to the given reference data and metadata to simulate new omics data.
    Args:
        count_path (str): The path to the count matrix file.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information.
        spatial_x (str): The column name in the metadata that contains the x-coordinate of the spatial location.
        spatial_y (str): The column name in the metadata that contains the y-coordinate of the spatial location.
        new_meta (pd.DataFrame): A DataFrame containing simulated spatial metadata for omics simulation, which is derived from the .meta of the simspace object.
        seed (int): The random seed for reproducibility.
        r_script_path (str): The path to the R script that performs the scDesign fitting. Default is None, which uses the script of simspace package.
    Returns:
        tuple: A tuple containing two DataFrames:
            - sim_data: The simulated gene expression data.
            - sim_meta: The metadata of the simulated cells.
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    Examples:
        >>> sim_data, sim_meta = scdesign_fit("path/to/count.csv",
        ...                              "path/to/meta.csv",
        ...                              "group_column",
        ...                              "x_coordinate",
        ...                              "y_coordinate",
        ...                              new_meta=sim.meta, # sim is simulated simspace object
        ...                              seed=42)
        >>> print(sim_data.head())
        >>> print(sim_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/scdesign.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    new_meta.to_csv("./tmp/new_meta.csv", index=False)

    result = subprocess.run(
        ["Rscript", r_script_path, meta_path, count_path, group_col, spatial_x, spatial_y, str(seed)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stdout)
        print(result.stderr)
        return None, None
    else:
        print("scDesgin fit complete.")
        sim_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        sim_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        os.remove("./tmp/new_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return sim_data, sim_meta
    

def srtsim_fit(count_path, 
               meta_path, 
               group_col = 'state', 
               spatial_x = 'x',
               spatial_y = 'y',
               n_rep = 1,
               seed = 0,
               r_script_path=None):
    """
    Fit the SRTsim model to the given reference data and metadata to simulate new omics data.
    Args:
        count_path (str): The path to the count matrix file.
        meta_path (str): The path to the metadata file.
        group_col (str): The column name in the metadata that contains the grouping information. Default is 'state'.
        spatial_x (str): The column name in the metadata that contains the x-coordinate of the spatial location.
        spatial_y (str): The column name in the metadata that contains the y-coordinate of the spatial location.
        n_rep (int): The number of replicates to simulate. Default is 1. Since SRTsim can only simulate exact same number of cells as the reference, this parameter is used when the number of cells in the reference is less than the number of cells in the spatial simulation results.
        seed (int): The random seed for reproducibility. Default is 0.
        r_script_path (str): The path to the R script that performs the SRTsim fitting. Default is None, which uses the script of simspace package.
    Returns:
        tuple: A tuple containing two DataFrames:
            - sim_data: The simulated gene expression data.
            - sim_meta: The metadata of the simulated cells.
    Raises:
        FileNotFoundError: If the R script file does not exist.
        Exception: If the R script fails to execute or returns an error.
    Examples:
        >>> sim_data, sim_meta = srtsim_fit("path/to/count.csv",
        ...                              "path/to/meta.csv",
        ...                              group_col='state',
        ...                              spatial_x='x',
        ...                              spatial_y='y',
        ...                              n_rep=1,
        ...                              seed=42)
        >>> print(sim_data.head())
        >>> print(sim_meta.head())
    """
    if r_script_path is None:
        r_script_path = os.path.join(os.path.dirname(__file__), "R/srtsim.R")
    if not os.path.exists(r_script_path):
        raise FileNotFoundError(f"The R script {r_script_path} does not exist. Please provide a valid path.")
    if not os.path.exists("./tmp"):
        os.makedirs("./tmp")
        print("Temporary directory created.")

    result = subprocess.run(
        ["Rscript", r_script_path, meta_path, count_path, group_col, spatial_x, spatial_y, str(n_rep), str(seed)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("R script failed:")
        print(result.stderr)
        return None, None
    else:
        print("SRTsim fit complete.")
        sim_data = pd.read_csv('./tmp/simulated_data.csv', sep=",", header=0, index_col=0)
        sim_meta = pd.read_csv('./tmp/simulated_meta.csv', sep=",", header=0, index_col=0)
        # Clean up temporary files
        os.remove("./tmp/simulated_data.csv")
        os.remove("./tmp/simulated_meta.csv")
        if os.path.exists("./tmp"):
            os.rmdir("./tmp")
        return sim_data, sim_meta
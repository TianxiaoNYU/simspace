import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def spatial_pie(
        SimSpace, 
        spot_meta: pd.DataFrame,
        kernel: tuple = (5, 5),
        figure_size: tuple = (5, 5),
        dpi: int = 300,
        save_path: str = None,
        ) -> None:
    """
    Plot the spatial pie chart of the convolved SimSpace dataset.

    Args:
        SimSpace: The SimSpace object containing the spatial data.
        spot_meta: DataFrame containing the metadata for each spot, including state proportions. Should have at least three columns: 'col' and 'row' as the first two columns, and state proportions.
        kernel: The kernel used for convolution, which determines the size of the pie chart.
        figure_size: Tuple specifying the size of the figure (width, height).
        dpi: Dots per inch for the figure resolution.
        save_path: Path to save the figure. If None, the figure will be displayed instead.

    Returns:
        None: Displays the pie chart or saves it to the specified path.

    Raises:
        ValueError: If the spot_meta DataFrame does not contain the expected columns. 
        ValueError: If spot_meta does not contain 'col' and 'row' columns

    """
    if not isinstance(spot_meta, pd.DataFrame):
        raise ValueError("spot_meta must be a pandas DataFrame.")
    if spot_meta.shape[1] < 3:
        raise ValueError("spot_meta must contain at least three columns: 'col', 'row', and state proportions.")
    if 'col' not in spot_meta.columns or 'row' not in spot_meta.columns:
        raise ValueError("spot_meta must contain 'col' and 'row' columns for spatial coordinates.")

    # Plot the outcome of mixing
    cmap = sns.color_palette('tab20', n_colors=SimSpace.num_states)
    state_names = spot_meta.columns[2:]
    state_name_mapping = {i: name for i, name in enumerate(state_names)}
    state_colors = {state_name_mapping[i]: cmap[i] for i in range(len(state_names))}
    # print(state_colors[0])
    fig, ax = plt.subplots()
    fig.set_size_inches(figure_size)
    fig.set_dpi(dpi)
    ax.set_aspect('equal')
    ax.set_xlim(-(SimSpace.size[0]/100*2), SimSpace.size[0] + (SimSpace.size[0]/100*2))
    ax.set_ylim(-(SimSpace.size[1]/100*2), SimSpace.size[1] + (SimSpace.size[1]/100*2))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Convolved SimSpace Dataset')
    for i in range(len(spot_meta)):
        centroid_x = spot_meta.iloc[i]['col']
        centroid_y = spot_meta.iloc[i]['row']
        state_proportions = spot_meta.iloc[i][2:]
        state_proportions = state_proportions[state_proportions > 0]
        state_proportions = state_proportions / state_proportions.sum()
        # for j, state in enumerate(state_proportions.index):
            # ax.add_patch(plt.Circle((centroid_x, centroid_y), state_proportions[state] * 3, color=state_colors[state], alpha=0.5))
        _, _ = ax.pie(state_proportions, 
                      colors=[state_colors[i] for i in state_proportions.index], 
                      startangle=90, 
                      radius=kernel[0]/3, 
                      center=(centroid_x, centroid_y), 
                      frame=True,
                      )
    # ax.invert_yaxis()
    if save_path is not None:
        plt.savefig(save_path)
    else:
        plt.show()

def plot_gene(
        coords: pd.DataFrame,
        feature: pd.Series,
        size=10,
        save_path=None,
        figsize=(6, 6),
        dpi=200,
        cmap=None,
        title=None
        ):
    """
    Plot the gene expression level on the spatial coordinates.

    Args:
        coords: DataFrame containing the spatial coordinates with columns 'col' and 'row'.
        feature: Series containing the gene expression levels, indexed by the same index as coords.
        size: Size of the scatter points.
        save_path: Path to save the figure. If None, the figure will be displayed instead.
        figsize: Tuple specifying the size of the figure (width, height).
        dpi: Dots per inch for the figure resolution.
        cmap: Colormap for the scatter plot. If None, a default colormap will be used.
        title: Title of the plot. If None, the name of the feature will be used.
    
    Returns:
        None: Displays the scatter plot or saves it to the specified path.

    Raises:
        ValueError: If coords does not contain 'col' and 'row' columns, or if feature is not a Series indexed by coords.
        TypeError: If coords or feature are not of the expected types.
    """
    if not isinstance(coords, pd.DataFrame):
        raise TypeError("coords must be a pandas DataFrame.")
    if not isinstance(feature, pd.Series):
        raise TypeError("feature must be a pandas Series.")
    if 'col' not in coords.columns or 'row' not in coords.columns:
        raise ValueError("coords must contain 'col' and 'row' columns.")
    if coords.shape[0] != feature.shape[0]:
        raise ValueError("coords and feature must have the same number of rows.")
    
    feature_tmp = feature.copy()
    coords_tmp = coords.copy()
    feature_tmp.index = coords_tmp.index
    df = pd.concat([coords_tmp, feature_tmp], axis=1)
    if cmap is None:
        cmap = sns.color_palette('flare', as_cmap=True)

    fig, ax = plt.subplots()
    fig.set_size_inches(figsize)
    fig.set_dpi(dpi)
    ax.set_aspect('equal')
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(f'{feature.name}')
    scatter = ax.scatter(df['col'], df['row'], c=df[feature.name], s=size, cmap=cmap, edgecolor='none')
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(feature.name)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, figsize=figsize, dpi=dpi)
    else:
        plt.show()


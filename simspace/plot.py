import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def spatial_pie(
        SimSpace, 
        spot_meta,
        kernel,
        figure_size=(5, 5),
        dpi=300,
        save_path=None,
        ):
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
    """
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


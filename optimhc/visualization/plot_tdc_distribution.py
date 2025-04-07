import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import logging
from optimhc.psm_container import PsmContainer
from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)

def visualize_target_decoy_features(psms: PsmContainer, num_cols=5, save_path=None, **kwargs):
    """
    Visualize the distribution of features in a DataFrame using kernel density estimation plots.

    Parameters:
        psms (PsmContainer): A PsmContainer object containing the features to visualize.
        num_cols (int, optional): The number of columns in the plot grid. Defaults to 5.
        save_path (str, optional): The file path to save the plot. If not provided, the plot is displayed.
        **kwargs: Additional plotting parameters such as `figsize` and `dpi`, etc.
    """
    rescoring_features = [item for sublist in psms.rescoring_features.values() for item in sublist if item != psms.hit_rank_column]
    num_features = len(rescoring_features)
    num_rows = (num_features + num_cols - 1) // num_cols
    
    figsize = kwargs.get('figsize', (15, num_rows * 15 / num_cols))
    dpi = kwargs.get('dpi', 300)
    
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize, dpi=dpi)
    axes = axes.flatten()
    
    psms_top_hits = psms.psms[psms.psms[psms.hit_rank_column] == 1].copy()
    # report the the number of true hits and decoys
    num_true_hits = len(psms_top_hits[psms_top_hits[psms.label_column] == True])
    num_decoys = len(psms_top_hits[psms_top_hits[psms.label_column] == False])
    logger.debug(f"Number of true hits: {num_true_hits}")
    logger.debug(f"Number of decoys: {num_decoys}")
    psms_top_hits[psms.label_column] = psms_top_hits[psms.label_column].map({True: 'Target', False: 'Decoy'})
    

    for i, feature in enumerate(rescoring_features):
        ax = axes[i]
        sns.histplot(
            data=psms_top_hits,
            x=feature,
            hue=psms.label_column,
            ax=ax,
            bins='auto',  
            common_bins=True,
            multiple='dodge',
            fill=True,
            alpha=0.3,
            stat='frequency', 
            kde=True,
            linewidth=0
        )
        ax.set_title(feature)
        ax.set_xlabel("")
        ax.set_ylabel("")

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    save_or_show_plot(save_path, logger)


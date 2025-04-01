import os
import matplotlib.pyplot as plt
from optimhc.visualization.save_or_show_plot import save_or_show_plot
import logging

logger = logging.getLogger(__name__)

def plot_qvalues(results, save_path=None, dpi=300, figsize=(15, 10), threshold=0.05, colors=None, **kwargs):
    """
    Plot q-values for the given results.

    Parameters:
        results: A list of results objects or a single result object.
                 Each result object should have a method `plot_qvalues`.
        save_path: Optional. If provided, saves the plot to the specified path.
        dpi: The resolution of the plot.
        figsize: The size of the figure.
        threshold: The q-value threshold for plotting.
        colors: A list of colors for the plots.
        **kwargs: Additional plotting parameters.
    
    """
    if not isinstance(results, list):
        results = [results]

    if colors is None:
        colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
            "#bcbd22", "#17becf"
        ]

    fig, axs = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
    
    for i, result in enumerate(results):
        for ax, level in zip(axs, ['psms', 'peptides']):
            result.plot_qvalues(
                level=level,
                c=colors[i % len(colors)],
                ax=ax,
                threshold=threshold,
                label=f'Result {i+1}' if len(results) > 1 else 'Results',
                linewidth=1,
                **kwargs
            )
            ax.legend(frameon=False)
            ax.set_title(f'{level}')

    plt.tight_layout()
    return save_or_show_plot(save_path, logger)

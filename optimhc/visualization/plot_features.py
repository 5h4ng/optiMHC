import numpy as np
import matplotlib.pyplot as plt
import logging
from itertools import cycle
from matplotlib.patches import Patch
from optimhc.psm_container import PsmContainer
import seaborn as sns
from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)


def plot_feature_importance(
    models, rescoring_features, save_path=None, sort=False, error=False, **kwargs
):
    """
    Unified function to plot average feature importance across multiple models.

    This function supports:
      - Linear models (e.g., Linear SVR) which provide an 'estimator' attribute with a 'coef_'.
        The absolute value of the coefficients is used for importance, and hatch patterns are applied
        to differentiate between positive and negative coefficients.
      - XGBoost models which provide a 'feature_importances_' attribute. Since these values are
        always positive, no hatch patterns are applied.

    Parameters
    ----------
    models : list
        A list of model objects.
        For linear models, each model should have an 'estimator' with 'coef_'.
        For XGBoost models, each model should have a 'feature_importances_' attribute.
    rescoring_features : dict
        A dictionary where keys are sources and values are lists of features.
    save_path : str, optional
        If provided, saves the plot to the specified path.
    sort : bool, optional
        If True, sorts the features by their importance in descending order.
        Default is False.
    error : bool, optional
        If True, adds error bars to the plot. Default is False.
    **kwargs : dict
        Additional plotting parameters such as 'figsize' and 'dpi'.

    Notes
    -----
    The function automatically detects the model type based on the presence of the corresponding attribute.
    For linear models, it uses hatch patterns to differentiate between positive and negative coefficients.
    For XGBoost models, it uses solid bars since the importances are always positive.
    """
    # Determine the model type based on the first model in the list.
    if hasattr(models[0].estimator, "coef_"):
        model_type = "linear"
    elif hasattr(models[0].estimator, "feature_importances_"):
        model_type = "xgb"
    else:
        raise ValueError(
            "Model type not recognized. Model must have 'estimator.coef_' for linear models or "
            "'estimator.feature_importances_' for XGBoost models."
        )

    if model_type == "linear":
        feature_importances = []
        for model in models:
            coefficients = model.estimator.coef_
            feature_importances.append(np.abs(coefficients).mean(axis=0))
            logger.debug(f"Model coefficients shape: {coefficients.shape}")

        average_feature_importance = np.mean(feature_importances, axis=0)
        std_feature_importance = np.std(feature_importances, axis=0)
        feature_signs = np.mean(
            [model.estimator.coef_.mean(axis=0) for model in models], axis=0
        )

    elif model_type == "xgb":
        feature_importances = []
        for model in models:
            # Use the XGBoost feature importances directly as they are always positive
            imp = model.estimator.feature_importances_
            feature_importances.append(imp)
            logger.debug(f"Model feature importances shape: {imp.shape}")

        average_feature_importance = np.mean(feature_importances, axis=0)
        std_feature_importance = np.std(feature_importances, axis=0)
        feature_signs = np.ones_like(average_feature_importance)

    logger.debug(
        f"Total rescoring features: {len(sum(rescoring_features.values(), []))}"
    )
    logger.debug(
        f"Average feature importance length: {len(average_feature_importance)}"
    )
    logger.debug(f"Features: {sum(rescoring_features.values(), [])}")

    figsize = kwargs.get("figsize", (15, 10))
    dpi = kwargs.get("dpi", 300)

    all_features = []
    all_importances = []
    all_errors = []
    all_colors = []
    all_hatches = []  # Hatch patterns will be applied only for linear models.

    color_cycle = cycle(plt.cm.tab10.colors)

    for source, features in rescoring_features.items():
        color = next(color_cycle)
        indices = [
            i
            for i, name in enumerate(sum(rescoring_features.values(), []))
            if name in features
        ]
        source_importances = average_feature_importance[indices]
        source_std = std_feature_importance[indices]

        if model_type == "linear":
            source_signs = feature_signs[indices]

        if sort:
            sorted_indices = np.argsort(-source_importances)
        else:
            sorted_indices = np.arange(len(features))

        sorted_features = [features[i] for i in sorted_indices]
        sorted_importances = source_importances[sorted_indices]
        sorted_std = source_std[sorted_indices]

        all_features.extend(sorted_features)
        all_importances.extend(sorted_importances)
        all_errors.extend(sorted_std)
        all_colors.extend([color] * len(sorted_features))

        if model_type == "linear":
            # For linear models, use hatch patterns to differentiate positive and negative coefficients.
            # An empty hatch ('') for positive and '\\' for negative coefficients.
            sorted_signs = source_signs[sorted_indices]
            all_hatches.extend(["" if sign >= 0 else "\\\\" for sign in sorted_signs])
        else:
            all_hatches.extend([""] * len(sorted_features))

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    if error:
        bars = ax.barh(
            all_features, all_importances, xerr=all_errors, color=all_colors, capsize=5
        )
    else:
        bars = ax.barh(all_features, all_importances, color=all_colors)

    if model_type == "linear":
        for bar, hatch in zip(bars, all_hatches):
            bar.set_hatch(hatch)
        legend_hatches = [
            Patch(facecolor="white", edgecolor="black", hatch="", label="Positive"),
            Patch(facecolor="white", edgecolor="black", hatch="\\\\", label="Negative"),
        ]
        legend_colors = [
            Patch(facecolor=color, edgecolor="black", label=source)
            for color, source in zip(plt.cm.tab10.colors, rescoring_features.keys())
        ]
        ax.legend(handles=legend_hatches + legend_colors, loc="best")
    else:
        legend_colors = [
            Patch(facecolor=color, edgecolor="black", label=source)
            for color, source in zip(plt.cm.tab10.colors, rescoring_features.keys())
        ]
        ax.legend(handles=legend_colors, loc="best")

    ax.set_xlabel("Average Feature Importance")
    ax.set_ylabel("Feature")

    save_or_show_plot(save_path, logger)


def visualize_feature_correlation(psms: PsmContainer, save_path=None, **kwargs):
    """
    Visualize the correlation between features in a DataFrame using a heatmap.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing the features to visualize.
    save_path : str, optional
        The file path to save the plot. If not provided, the plot is displayed.
    **kwargs : dict
        Additional plotting parameters such as `figsize` and `dpi`, etc.

    Notes
    -----
    This function:
    1. Extracts all rescoring features from the PsmContainer
    2. Calculates the correlation matrix between features
    3. Creates a heatmap visualization of the correlations
    4. Uses a coolwarm colormap to show positive and negative correlations
    """
    figsize = kwargs.get("figsize", (40, 36))
    dpi = kwargs.get("dpi", 300)

    rescoring_features = [
        item for sublist in psms.rescoring_features.values() for item in sublist
    ]
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    corr = psms.psms[rescoring_features].corr()
    # sns.heatmap(corr, annot=True, fmt=".2f", cmap='coolwarm', ax=ax)
    sns.heatmap(corr, cmap="coolwarm", ax=ax)
    ax.set_title("Feature Correlation Heatmap")

    save_or_show_plot(save_path, logger)

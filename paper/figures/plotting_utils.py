"""
Plotting utilities for paper figures.

This module contains reusable plotting components that preserve the exact
styling and formatting from the original notebook.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import os
import sys
from itertools import cycle
from typing import Optional, List, Dict, Tuple

# Add parent directory to path to import settings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import settings


def setup_plotting_theme():
    """
    Set up the seaborn theme for all plots.

    Original logic from notebook Cell 0.
    """
    sns.set_theme(style="whitegrid")


def format_method_names(data: pd.DataFrame, method_col: str = 'method') -> pd.DataFrame:
    """
    Apply method abbreviations and format method names.

    Original logic from notebook cells 7, 9, 11, etc.

    Args:
        data: DataFrame with method column
        method_col: Name of method column

    Returns:
        DataFrame with formatted method names
    """
    data[method_col] = data[method_col].map(
        lambda x: settings.METHOD_ABBREVIATIONS[x] if x in settings.METHOD_ABBREVIATIONS else x
    )
    data[method_col] = data[method_col].map(lambda x: x.replace('_', ' '))
    return data


def create_method_order(found_methods: set, filter_network_constraint: bool = False) -> List[str]:
    """
    Create ordered list of methods for seaborn plots.

    Maintains the order from settings.METHODS, applies abbreviations,
    and optionally filters out network constraint methods.

    Original logic from notebook cells 7, 9, 11, 13, etc.

    Args:
        found_methods: Set of methods actually found in data
        filter_network_constraint: If True, exclude network constraint methods

    Returns:
        Ordered list of formatted method names
    """
    # Step 1: Filter to found methods, maintaining order
    order = [method for method in settings.METHODS if method in found_methods]

    # Step 2: Apply abbreviations
    order = [settings.METHOD_ABBREVIATIONS[x] if x in settings.METHOD_ABBREVIATIONS else x
             for x in order]

    # Step 3: Replace underscores with spaces
    order = [x.replace('_', ' ') for x in order]

    # Step 4: Filter network constraint methods if requested
    if filter_network_constraint:
        order = [x for x in order if x not in settings.METHODS_NETWORK_CONSTRAINT]

    return order


def create_method_palette(order: List[str]) -> List[str]:
    """
    Create color palette for ordered methods.

    Original logic from notebook cells 13, 15, 20, 22.

    Args:
        order: Ordered list of method names

    Returns:
        List of colors corresponding to methods
    """
    palette = [settings.METHOD_PALETTE[x] for x in order if x in settings.METHOD_PALETTE]
    return palette


def create_simple_barplot(data: pd.DataFrame,
                          x_col: str,
                          y_col: str,
                          order: List[str],
                          ax: plt.Axes,
                          title: Optional[str] = None,
                          ylim: Tuple[float, float] = (0, 1),
                          capsize: float = 0.2,
                          palette: Optional[Dict] = None,
                          hide_x_labels: bool = False) -> plt.Axes:
    """
    Create a simple barplot without hue.

    Original logic from notebook cells 7, 9, 11.

    Args:
        data: DataFrame with plotting data
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        order: Order of x-axis categories
        ax: Matplotlib axis to plot on
        title: Plot title
        ylim: Y-axis limits
        capsize: Error bar cap size
        palette: Color palette (defaults to settings.METHOD_PALETTE)
        hide_x_labels: If True, hide x-axis tick labels

    Returns:
        Configured seaborn axis
    """
    if palette is None:
        palette = settings.METHOD_PALETTE

    sub_fig = sns.barplot(
        data=data,
        y=y_col,
        x=x_col,
        estimator=np.mean,
        errorbar=(lambda x: (min(x), max(x))),
        capsize=capsize,
        ax=ax,
        order=order,
        palette=palette
    )

    # Axis formatting
    if title:
        ax.set_title(title, fontsize=18, y=1.0, pad=5)

    ax.tick_params(axis='x', rotation=90)
    ax.set(ylim=ylim)
    sub_fig.set(xlabel=None, ylabel=None)
    sns.despine(bottom=True, left=True, right=True, top=True)

    # Hide x-axis labels if requested
    if hide_x_labels:
        ax.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False
        )

    return sub_fig


def create_hue_barplot(data: pd.DataFrame,
                       x_col: str,
                       y_col: str,
                       hue_col: str,
                       order: List[str],
                       hue_order: List[str],
                       ax: plt.Axes,
                       ylim: Tuple[float, float] = (0, 1),
                       capsize: float = 0.05,
                       legend_loc: str = 'upper right') -> plt.Axes:
    """
    Create a barplot with hue (e.g., TCGA vs METABRIC).

    Uses gray palette initially - colors applied via hatches later.

    Original logic from notebook cells 13, 15, 20, 22.

    Args:
        data: DataFrame with plotting data
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        hue_col: Column name for hue grouping
        order: Order of x-axis categories
        hue_order: Order of hue categories
        ax: Matplotlib axis to plot on
        ylim: Y-axis limits
        capsize: Error bar cap size
        legend_loc: Legend location

    Returns:
        Configured seaborn axis
    """
    # Gray palette - colors applied later with hatches
    palette = ['grey', 'grey']

    sub_fig = sns.barplot(
        palette=palette,
        hue_order=hue_order,
        data=data,
        x=x_col,
        y=y_col,
        hue=hue_col,
        ax=ax,
        estimator=np.mean,
        errorbar=(lambda x: (min(x), max(x))),
        capsize=capsize,
        order=order
    )

    # Legend and axis formatting
    plt.legend(loc=legend_loc)
    ax.tick_params(axis='x', rotation=90)
    ax.set(ylim=ylim)

    return sub_fig


def apply_hatch_pattern(sub_fig: plt.Axes,
                        palette: List[str],
                        hatches: List[str] = ['', '//']) -> None:
    """
    Apply color and hatch patterns to distinguish datasets.

    Original logic from notebook cells 13, 15, 20, 22.

    Args:
        sub_fig: Seaborn barplot axis
        palette: List of colors for each bar group
        hatches: Hatch patterns for each hue level
    """
    for bars, hatch, legend_handle in zip(sub_fig.containers, hatches, sub_fig.legend_.legendHandles):
        for bar, color in zip(bars, palette):
            bar.set_facecolor(color)
            bar.set_hatch(hatch)
        # Update the existing legend, use twice the hatching pattern to make it denser
        legend_handle.set_hatch(hatch + hatch)


def apply_violin_hatch_pattern(sub_fig: plt.Axes,
                                order: List[str]) -> None:
    """
    Apply hatch patterns to violin plot.

    Original logic from notebook cell 24.

    Args:
        sub_fig: Seaborn violin plot axis
        order: Order of methods
    """
    hatches = cycle(['', '//'])
    palette = iter([settings.METHOD_PALETTE[x] for x in order if x in settings.METHOD_PALETTE])
    prev_i = None
    color = next(palette)

    for i, patch in enumerate(sub_fig.get_children()):
        if not isinstance(patch, matplotlib.collections.PolyCollection):
            continue
        # Boxes from left to right
        hatch = next(hatches)
        patch.set_hatch(hatch)
        if prev_i is not None and i - prev_i == 1:
            color = next(palette)
        patch.set_facecolor(color)
        patch.set_edgecolor('white')
        prev_i = i

    # Update legend
    for i, hatch in enumerate(['', '//']):
        sub_fig.legend_.legend_handles[i].set_hatch(hatch + hatch)
        sub_fig.legend_.legend_handles[i].set_edgecolor('white')


def create_violinplot(data: pd.DataFrame,
                      x_col: str,
                      y_col: str,
                      hue_col: str,
                      order: List[str],
                      ax: plt.Axes,
                      ylim: Tuple[float, float] = (-0.2, 1),
                      palette: Optional[List] = None) -> plt.Axes:
    """
    Create a violin plot for method comparison.

    Original logic from notebook cell 24.

    Args:
        data: DataFrame with plotting data
        x_col: Column name for x-axis
        y_col: Column name for y-axis
        hue_col: Column name for hue grouping
        order: Order of x-axis categories
        ax: Matplotlib axis to plot on
        ylim: Y-axis limits
        palette: Color palette

    Returns:
        Configured seaborn axis
    """
    if palette is None:
        palette = ['grey', 'grey']

    sub_fig = sns.violinplot(
        data=data,
        gap=.1,
        split=True,
        fill=False,
        palette=palette,
        x=x_col,
        y=y_col,
        hue=hue_col,
        ax=ax,
        scale="width",
        cut=0,
        order=order
    )

    # Axis formatting
    ax.tick_params(axis='x', rotation=90)
    ax.set(ylim=ylim)
    sub_fig.set(xlabel=None, ylabel=None)
    sns.despine(bottom=False, left=True, right=True, top=False)

    return sub_fig


def create_figure_legend(fig: plt.Figure,
                         palette_dict: Dict[str, str],
                         title: str = "",
                         fontsize: str = 'medium',
                         loc: int = 1,
                         bbox_to_anchor: Tuple[float, float] = (0.25, 0.87)):
    """
    Create custom legend from palette dictionary.

    Original logic from notebook cells 7, 9, 11.

    Args:
        fig: Matplotlib figure
        palette_dict: Dictionary mapping labels to colors
        title: Legend title
        fontsize: Font size for legend
        loc: Legend location code
        bbox_to_anchor: Legend anchor position

    Returns:
        Legend object
    """
    leg = fig.legend(
        palette_dict,
        title=title,
        fontsize=fontsize,
        fancybox=True,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor
    )

    # Set colors for legend handles
    for i, (key, color) in enumerate(palette_dict.items()):
        leg.legendHandles[i].set_color(color)

    return leg


def add_ylabel(fig: plt.Figure,
               text: str,
               x: float = 0.08,
               y: float = 0.5,
               fontsize: int = 18) -> None:
    """
    Add rotated y-axis label at figure level.

    Original logic from notebook cells 7, 9, 11, 13, etc.

    Args:
        fig: Matplotlib figure
        text: Label text
        x: X position (figure coordinates)
        y: Y position (figure coordinates)
        fontsize: Font size
    """
    fig.text(x, y, text, va='center', rotation='vertical', fontsize=fontsize)


def save_figure(fig: plt.Figure,
                filename: str,
                dpi: int = 300,
                save_svg: bool = False,
                output_dir: str = 'paper') -> None:
    """
    Save figure to file.

    Original logic from notebook (multiple cells).

    Args:
        fig: Matplotlib figure
        filename: Filename without extension
        dpi: Resolution for PNG
        save_svg: If True, also save SVG version
        output_dir: Output directory
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save PNG
    png_path = os.path.join(output_dir, f'{filename}.png')
    plt.savefig(png_path, dpi=dpi, bbox_inches='tight')

    # Optionally save SVG
    if save_svg:
        svg_path = os.path.join(output_dir, f'{filename}.svg')
        plt.savefig(svg_path)


def rename_tcga_to_brca(data: pd.DataFrame, dataset_col: str = 'Dataset') -> pd.DataFrame:
    """
    Rename 'TCGA' to 'TCGA-BRCA' for display purposes.

    Original logic from notebook cells 13, 15, 20, 22, 24.

    Args:
        data: DataFrame with dataset column
        dataset_col: Name of dataset column

    Returns:
        DataFrame with renamed values
    """
    data[dataset_col] = data[dataset_col].map(
        lambda x: 'TCGA-BRCA' if x == 'TCGA' else x
    )
    return data


def prepare_long_format_data(df_wide: pd.DataFrame,
                              method_name: str,
                              dataset: Optional[str] = None,
                              value_name: str = 'Jaccard index',
                              var_name: str = 'Cancer type') -> pd.DataFrame:
    """
    Convert wide format data to long format for seaborn.

    Original logic from notebook cells 13, 20, 22.

    Args:
        df_wide: Wide format DataFrame
        method_name: Method name to add as column
        dataset: Optional dataset filter
        value_name: Name for value column after melting
        var_name: Name for variable column after melting

    Returns:
        Long format DataFrame
    """
    # Filter to specific dataset if provided
    if dataset is not None and 'dataset' in df_wide.columns:
        df_wide = df_wide[df_wide['dataset'] == dataset]
        del df_wide['dataset']

    # Melt to long format
    df_long = df_wide.melt()

    # Add method column
    df_long['Method'] = method_name

    # Rename columns
    df_long = df_long.rename(columns={'value': value_name, 'variable': var_name})

    return df_long


def create_method_dataframe_for_violin(data_dict: Dict[str, List],
                                        dataset: str) -> pd.DataFrame:
    """
    Create DataFrame from method-keyed dictionary for violin plot.

    Original logic from notebook cell 24.

    Args:
        data_dict: Dictionary mapping method names to value lists
        dataset: Dataset name to add as column

    Returns:
        Long format DataFrame with 'variable', 'value', and 'dataset' columns
    """
    # Create series for each method
    series_list = []
    for method, values in data_dict.items():
        s = pd.Series(values, name=method)
        series_list.append(s)

    # Concatenate and melt
    df_plot = pd.concat(series_list, axis=1)
    df_plot = df_plot.melt()

    # Filter network constraint methods if they exist
    if hasattr(settings, 'METHODS_NETWORK_CONSTRAINT'):
        df_plot = df_plot[~df_plot['variable'].isin(settings.METHODS_NETWORK_CONSTRAINT)]

    # Add dataset label
    df_plot['dataset'] = dataset

    return df_plot

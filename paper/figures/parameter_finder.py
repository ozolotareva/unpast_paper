"""
Parameter optimization utilities for paper figures.

This module contains logic for finding best and default parameters
for each method, preserving the exact logic from the original notebook.
"""

import pandas as pd
import sys
import os

import settings


def calc_average(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group data by scenario, gsize, and parameters, then average.

    Original logic from notebook Cell 7.

    Args:
        df: DataFrame with scenario, gsize, parameters columns

    Returns:
        Averaged DataFrame with reset index
    """
    df = df.fillna(0)
    return df.groupby(['scenario', 'gsize', 'parameters']).mean().reset_index()


def calc_best_params(df: pd.DataFrame, metric_col: str) -> str:
    """
    Find the best parameter combination across all runs.

    Merges runs with different seeds, then finds the parameter
    combination with highest mean metric value.

    Original logic from notebook Cell 7.

    Args:
        df: DataFrame with parameters and metric columns
        metric_col: Name of metric column to optimize (e.g., 'performance', 'PAM50')

    Returns:
        Parameter string of best configuration
    """
    # Merge runs with different seeds if needed
    if 'run' in df.columns or not df["parameters"].is_unique:
        df = calc_average(df)

    # Find parameter combination with highest mean metric
    top_params = (df.groupby('parameters')
                    .mean()
                    .reset_index()
                    .sort_values(metric_col, ascending=False)
                    .iloc[0]['parameters'])

    return top_params


def find_default_params(df: pd.DataFrame, metric_col: str, method: str) -> str:
    """
    Find the parameter combination matching default parameters from settings.

    Searches through DataFrame to find row where all default parameters
    are present in the parameter string.

    Original logic from notebook Cell 9.

    Args:
        df: DataFrame with parameters column
        metric_col: Name of metric column (not used but kept for consistency)
        method: Method name to lookup defaults in settings

    Returns:
        Parameter string of default configuration

    Raises:
        Exception: If default parameters cannot be found
    """
    default_params = settings.DEFAULT_PARAMETERS[method]

    mode = False
    for i, row in df.iterrows():
        mode = True
        for key, value in default_params.items():
            substring = f'{key}={value}'
            if substring in row['parameters']:
                continue
            mode = False

        if mode is True:
            break

    if mode is False:
        print(f'Could not find default parameters for {method}')
        raise Exception('Could not find default parameters')

    return row['parameters']


def get_best_run_values(df: pd.DataFrame, metric_col: str, params: str, method: str) -> pd.DataFrame:
    """
    Extract all run values for a specific parameter configuration.

    Returns values from all runs (not averaged) to enable error bar
    calculation in plots.

    Original logic from notebook Cell 7.

    Args:
        df: DataFrame with parameters and metric columns
        metric_col: Name of metric column
        params: Parameter string to filter on
        method: Method name

    Returns:
        DataFrame in seaborn format with 'value' and 'method' columns
    """
    # Get all run values (needed for error calculation)
    top_values = df[df['parameters'] == params][metric_col].values

    # Format as seaborn input
    df_top = pd.DataFrame(top_values, columns=['value'])
    df_top['method'] = method

    return df_top


def get_default_run_values(df: pd.DataFrame, metric_col: str, params: str, method: str) -> pd.DataFrame:
    """
    Extract all run values for default parameter configuration.

    Identical to get_best_run_values but named differently for clarity.

    Original logic from notebook Cell 9.

    Args:
        df: DataFrame with parameters and metric columns
        metric_col: Name of metric column
        params: Parameter string to filter on
        method: Method name

    Returns:
        DataFrame in seaborn format with 'value' and 'method' columns
    """
    # Get all run values (needed for error calculation)
    top_values = df[df['parameters'] == params][metric_col].values

    # Format as seaborn input
    df_top = pd.DataFrame(top_values, columns=['value'])
    df_top['method'] = method

    return df_top


def find_best_average_rank_params(df_tcga: pd.DataFrame, df_metabric: pd.DataFrame, metric_col: str) -> str:
    """
    Find parameter combination with best average rank across TCGA and METABRIC.

    This is used for finding parameters that work well across both datasets,
    rather than optimizing for one specific dataset.

    Process:
    1. Rank parameters independently in each dataset
    2. Find shared parameter configurations
    3. Sum ranks across datasets
    4. Return parameter with minimum rank sum

    Original logic from notebook Cell 5 and used in Cell 15.

    Args:
        df_tcga: TCGA results DataFrame
        df_metabric: METABRIC results DataFrame
        metric_col: Metric column to rank by

    Returns:
        Parameter string with best average rank
    """
    # Rank in each dataset independently (higher metric = better = lower rank number)
    df_tcga_ranked = df_tcga.sort_values(metric_col, ascending=False).copy()
    df_tcga_ranked['rank'] = list(range(len(df_tcga_ranked.index)))

    df_metabric_ranked = df_metabric.sort_values(metric_col, ascending=False).copy()
    df_metabric_ranked['rank'] = list(range(len(df_metabric_ranked.index)))

    # Get shared parameter configurations
    index_shared = df_tcga_ranked.index[df_tcga_ranked.index.isin(df_metabric_ranked.index)]

    # Concatenate and sum ranks
    df_combined = pd.concat([df_tcga_ranked, df_metabric_ranked]).groupby('parameters').sum()

    # Find parameter with minimum rank sum
    best_params = df_combined.sort_values('rank', ascending=True).index[0]

    return best_params


def replace_string_in_list(arr: list, orig: str, rep: str) -> list:
    """
    Replace a string in a list with another string.

    Original logic from notebook Cell 7.

    Args:
        arr: List to modify
        orig: Original string to find
        rep: Replacement string

    Returns:
        Modified list
    """
    index = arr.index(orig)
    assert index > -1
    arr[index] = rep
    return arr


def calculate_performance_delta(df: pd.DataFrame,
                                 metric_col: str,
                                 best_params: str,
                                 default_params: str,
                                 method: str) -> pd.DataFrame:
    """
    Calculate performance improvement: best - default.

    Used for figures showing performance gain from parameter optimization.

    Original logic from notebook Cell 11.

    Args:
        df: DataFrame with all runs
        metric_col: Metric column name
        best_params: Best parameter string
        default_params: Default parameter string
        method: Method name

    Returns:
        DataFrame with delta values in seaborn format
    """
    # Get values for both parameter sets
    df_best = get_best_run_values(df, metric_col, best_params, method)
    df_default = get_default_run_values(df, metric_col, default_params, method)

    # Calculate delta
    df_best['value'] = df_best['value'] - df_default['value']

    return df_best

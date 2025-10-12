"""
Data loading and preprocessing utilities for paper figures.

This module centralizes all data loading logic from plot_paper_figures.ipynb,
preserving the exact transformation order to ensure figures remain identical.

CRITICAL: The order of transformations matters! Any changes to the sequence
may alter the final figures.
"""

import pandas as pd
import os
from typing import Optional, List

# Get the directory containing this file to construct data paths
_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def remove_substring_from_params(x: str, string: str, chars_right: int) -> str:
    """
    Remove a substring and characters to its right from parameter string.

    Original logic from notebook Cell 0.

    Args:
        x: Parameter string to modify
        string: Substring to find and remove
        chars_right: Number of characters to remove after the substring

    Returns:
        Modified parameter string
    """
    x = str(x)
    pos = x.find(string)
    if pos > -1:
        return x[:pos] + x[pos + len(string) + chars_right:]
    return x


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove metadata columns that are not needed for plotting.

    Original logic from notebook Cell 0.

    Args:
        df: DataFrame to clean

    Returns:
        DataFrame with metadata columns removed
    """
    columns_to_remove = ['seed', 'run', 'time', 'parameters', 'Unnamed: 0']
    for col in columns_to_remove:
        if col in df.columns:
            del df[col]
    return df


def load_raw_data(file_path: str) -> pd.DataFrame:
    """
    Load raw TSV file and perform initial transformations.

    Step 1 of the data loading pipeline.

    Args:
        file_path: Path to TSV file

    Returns:
        DataFrame with initial transformations applied

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    df = pd.read_csv(file_path, sep='\t')

    # Handle duplicate columns (both 'simulated' and 'performance_simulated')
    if 'simulated' in df.columns and 'performance_simulated' in df.columns:
        del df['simulated']

    return df


def extract_grandforest_seeds(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract seed information from grandforest parameter strings.

    Grandforest stores seeds in the parameter string. This function
    extracts them to separate columns.

    Original logic from notebook Cell 9.

    Args:
        df: DataFrame potentially containing grandforest data

    Returns:
        DataFrame with seed/run columns added if grandforest detected
    """
    # Check if this is grandforest data (seed not in columns but in parameter string)
    if 'seed' not in df.columns and len(df) > 0 and 'seeds=' in str(df.iloc[0]['parameters']):
        # Assumes that seeds= is last parameter (grandforest convention)
        df['seed'] = df['parameters'].map(lambda x: x.split('seeds=')[-1])
        df['run'] = df['parameters'].map(lambda x: x.split('seeds=')[-1])
        df['parameters'] = df['parameters'].map(lambda x: x.split('seeds=')[0])

    return df


def standardize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename columns to standard names used throughout analysis.

    Step 2 of the data loading pipeline.

    Original logic from notebook Cell 9.

    Args:
        df: DataFrame with potentially non-standard column names

    Returns:
        DataFrame with standardized column names
    """
    rename_map = {
        'param': 'parameters',
        'simulated': 'performance',
        'performance_simulated': 'performance'
    }

    df = df.rename(columns=rename_map)

    return df


def clean_parameter_strings_simulated(df: pd.DataFrame, method: Optional[str] = None) -> pd.DataFrame:
    """
    Clean parameter strings for simulated data.

    Removes:
    - random_state parameter
    - All scenario/gsize path combinations (A/5, A/50, etc.)
    - Special handling for mclust: k=autotuned → k=default

    Original logic from notebook Cell 9.

    CRITICAL: Order of operations must be preserved exactly.

    Args:
        df: DataFrame with parameters column
        method: Method name for special case handling (e.g., 'mclust')

    Returns:
        DataFrame with cleaned parameter strings
    """
    # Remove random_state parameter
    df['parameters'] = df['parameters'].map(
        lambda x: remove_substring_from_params(x, 'random_state=', 1)
    )

    # Remove all scenario/gsize combinations
    scenario_gsize_combinations = [
        '/A/5/', '/A/50/', '/A/500/',
        '/B/5/', '/B/50/', '/B/500/',
        '/C/5/', '/C/50/', '/C/500/'
    ]

    for combination in scenario_gsize_combinations:
        df['parameters'] = df['parameters'].map(
            lambda x: remove_substring_from_params(x, combination, 0)
        )

    # Special case for mclust: replace 'k=autotuned' with 'k=default'
    if method == 'mclust':
        df['parameters'] = df['parameters'].map(
            lambda x: x.replace('k=autotuned', 'k=default')
        )

    return df


def clean_parameter_strings_real(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean parameter strings for real (BRCA) data.

    Removes:
    - random_state parameter
    - Dataset names (TCGA, METABRIC)

    Original logic from notebook Cell 5.

    Args:
        df: DataFrame with parameters column

    Returns:
        DataFrame with cleaned parameter strings
    """
    # Remove random_state parameter
    df['parameters'] = df['parameters'].map(
        lambda x: remove_substring_from_params(x, 'random_state=', 1)
    )

    # Remove dataset names
    df['parameters'] = df['parameters'].map(
        lambda x: remove_substring_from_params(x, 'TCGA', 0)
    )
    df['parameters'] = df['parameters'].map(
        lambda x: remove_substring_from_params(x, 'METABRIC', 0)
    )

    return df


def aggregate_runs(df: pd.DataFrame, groupby_cols: List[str] = ['parameters']) -> pd.DataFrame:
    """
    Aggregate multiple runs by taking mean.

    Original logic from notebook Cell 5.

    Args:
        df: DataFrame with potentially multiple runs
        groupby_cols: Columns to group by

    Returns:
        DataFrame with runs aggregated by mean
    """
    return df.groupby(groupby_cols).mean()


def aggregate_runs_with_index(df: pd.DataFrame, groupby_cols: List[str]) -> pd.DataFrame:
    """
    Aggregate multiple runs by taking mean and reset index.

    Original logic from notebook Cell 7.

    Args:
        df: DataFrame with potentially multiple runs
        groupby_cols: Columns to group by

    Returns:
        DataFrame with runs aggregated by mean and index reset
    """
    df = df.fillna(0)
    return df.groupby(groupby_cols).mean().reset_index()


def load_simulated_data(method: str, base_path: str = None) -> pd.DataFrame:
    """
    Load and preprocess simulated data for a single method.

    This is the main entry point for loading simulated data. It applies
    all transformations in the exact order used in the original notebook.

    Transformation order (CRITICAL - DO NOT CHANGE):
    1. Load raw TSV
    2. Handle duplicate columns
    3. Standardize column names
    4. Extract grandforest seeds (if applicable)
    5. Clean parameter strings

    Args:
        method: Method name (e.g., 'unpast', 'fabia')
        base_path: Base directory containing data files (defaults to data/simulated_data in module directory)

    Returns:
        Preprocessed DataFrame

    Raises:
        FileNotFoundError: If data file doesn't exist
    """
    if base_path is None:
        base_path = os.path.join(_MODULE_DIR, 'data', 'simulated_data')

    file_path = os.path.join(base_path, f'{method}_ABC.tsv')

    # Step 1: Load raw data
    df = load_raw_data(file_path)

    # Step 2: Standardize column names
    df = standardize_column_names(df)

    # Step 3: Extract grandforest seeds if needed
    df = extract_grandforest_seeds(df)

    # Step 4: Clean parameter strings
    df = clean_parameter_strings_simulated(df, method)

    return df


def load_real_data(method: str, dataset: str, base_path: str = None) -> pd.DataFrame:
    """
    Load and preprocess real (BRCA) data for a single method and dataset.

    This is the main entry point for loading real data. It applies
    all transformations in the exact order used in the original notebook.

    Transformation order (CRITICAL - DO NOT CHANGE):
    1. Load raw TSV
    2. Handle duplicate columns
    3. Standardize column names
    4. Extract grandforest seeds (if applicable)
    5. Clean parameter strings (remove random_state and dataset names)

    Args:
        method: Method name (e.g., 'unpast', 'fabia')
        dataset: Dataset name ('TCGA' or 'METABRIC')
        base_path: Base directory containing data files (defaults to data/real_data in module directory)

    Returns:
        Preprocessed DataFrame

    Raises:
        FileNotFoundError: If data file doesn't exist
    """
    if base_path is None:
        base_path = os.path.join(_MODULE_DIR, 'data', 'real_data')

    file_path = os.path.join(base_path, f'{method}_{dataset}.tsv')

    # Step 1: Load raw data
    df = load_raw_data(file_path)

    # Step 2: Standardize column names
    df = standardize_column_names(df)

    # Step 3: Extract grandforest seeds if needed
    df = extract_grandforest_seeds(df)

    # Step 4: Clean parameter strings
    df = clean_parameter_strings_real(df)

    return df


def get_data_path(method: str, dataset: str, data_type: str = 'real') -> str:
    """
    Construct file path for method and dataset.

    Args:
        method: Method name
        dataset: Dataset name (for real: 'TCGA' or 'METABRIC', for simulated: 'ABC')
        data_type: 'real' or 'simulated'

    Returns:
        Full file path
    """
    if data_type == 'simulated':
        base_path = 'data/simulated_data'
        return os.path.join(base_path, f'{method}_ABC.tsv')
    else:  # real
        base_path = 'data/real_data'
        return os.path.join(base_path, f'{method}_{dataset}.tsv')

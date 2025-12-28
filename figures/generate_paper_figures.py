"""
Main script to generate all paper figures.

This script orchestrates the generation of all publication figures using
the modular components from data_loader, parameter_finder, and plotting_utils.

Usage:
    python generate_paper_figures.py [figure_name]

    If no figure_name is provided, generates all figures.
    Available figures: s3, figure3b, s4, s4_a, s5c, s5d, s6
"""

import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib

# Add paper directory to path so imports work regardless of where script is run from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import settings
import data_loader
import parameter_finder
import plotting_utils


# Paths
PATH_TO_SIMULATED_RESULTS = os.path.join(os.path.dirname(__file__), 'data', 'simulated_data')
PATH_TO_REAL_RESULTS = os.path.join(os.path.dirname(__file__), 'data', 'real_data')


def compute_best_average_rank_brca():
    """
    Compute best average rank parameters for BRCA datasets.

    This finds parameters that perform well across both TCGA and METABRIC
    by minimizing the sum of ranks.

    Original logic from notebook Cell 4.

    Returns:
        dict: Mapping of method -> best parameter string
    """
    print("Computing best average rank parameters for BRCA...")
    best_average_rank_brca_string = {}

    for method in settings.METHODS:
        try:
            # Load both datasets
            df_tcga = data_loader.load_real_data(method, 'TCGA')
            df_metabric = data_loader.load_real_data(method, 'METABRIC')

            # Aggregate runs
            df_tcga_mean = data_loader.aggregate_runs(df_tcga).reset_index()
            df_metabric_mean = data_loader.aggregate_runs(df_metabric).reset_index()

            # Find best average rank
            best_params = parameter_finder.find_best_average_rank_params(
                df_tcga_mean, df_metabric_mean, 'PAM50'
            )

            best_average_rank_brca_string[method] = best_params
            print(f"  {method}: {best_params[:80]}...")

        except Exception as e:
            print(f"  Warning: Could not process {method}: {e}")
            continue

    return best_average_rank_brca_string


def compute_best_params_per_dataset():
    """
    Compute best parameters for each method on each dataset independently.

    Original logic from notebook Cell 5.

    Returns:
        dict: Nested dict {dataset: {method: best_params_string}}
    """
    print("Computing best parameters per dataset...")
    best_params = {d: {} for d in settings.REAL_DATASETS}

    for dataset in settings.REAL_DATASETS:
        print(f"  Dataset: {dataset}")
        for method in settings.METHODS:
            try:
                df = data_loader.load_real_data(method, dataset)

                # Aggregate runs if needed
                if 'run' in df.columns or not df["parameters"].is_unique:
                    df_mean = data_loader.aggregate_runs(df)  # Don't reset index - parameters is the index
                else:
                    df_mean = df.set_index('parameters')

                # Clean columns but parameters is already the index so it's preserved
                df_mean = data_loader.clean_columns(df_mean)
                df_mean = df_mean['PAM50']

                # Best params = highest PAM50 (index is parameters)
                top_params = df_mean.sort_values(ascending=False).index[0]
                best_params[dataset][method] = top_params

            except Exception as e:
                print(f"    Warning: Could not process {method}: {e}")
                continue

    return best_params


def compute_default_params_per_dataset():
    """
    Compute default parameters for each method on each dataset.

    Original logic from notebook Cell 23.

    Returns:
        dict: Nested dict {dataset: {method: default_params_string}}
    """
    print("Computing default parameters per dataset...")
    default_params = {d: {} for d in settings.REAL_DATASETS}

    for dataset in settings.REAL_DATASETS:
        for method in settings.METHODS:
            try:
                df = data_loader.load_real_data(method, dataset)

                # Aggregate runs if needed
                if 'run' in df.columns or not df["parameters"].is_unique:
                    df_mean = data_loader.aggregate_runs(df)
                else:
                    df_mean = df.set_index('parameters')

                df_mean = data_loader.clean_columns(df_mean)
                df_mean = df_mean['PAM50']

                # Find default params (skip UnPaSt)
                if method != 'UnPaSt':
                    if method not in settings.DEFAULT_PARAMETERS:
                        continue

                    # Search for default params in index
                    flag = False
                    for params_string in df_mean.index:
                        flag = True
                        for k, v in settings.DEFAULT_PARAMETERS[method].items():
                            if not f'{k.lower()}={str(v).lower()}' in params_string.lower():
                                flag = False
                                break
                        if flag:
                            break

                    if flag:
                        default_params[dataset][method] = params_string
                    else:
                        print(f"  Could not find default parameters for {method}")

            except Exception as e:
                continue

    return default_params


def generate_figure_s3_1():
    """
    Generate Figure S3: Simulated data with default parameters.

    Original logic from notebook Cell 9.
    """
    print("\n=== Generating Figure S3 (Simulated data, default parameters) ===")

    plotting_utils.setup_plotting_theme()

    # Initialize data structure
    scenario_data = {s: {g: pd.DataFrame(columns=['method', 'value'])
                         for g in settings.GENE_SIZES}
                     for s in settings.SCENARIOS}
    table_key = 'performance'
    found_methods = set()

    # Collect data for all methods
    for method in settings.METHODS:
        try:
            df = data_loader.load_simulated_data(method)
            found_methods.add(method)

            # Find default parameters
            default_params = parameter_finder.find_default_params(df, table_key, method)
            print(f"{method}: {default_params[:60]}...")

            # Extract data for each scenario/gene_size combination
            for scenario in settings.SCENARIOS:
                for gene_size in settings.GENE_SIZES:
                    df_sub = df[(df['scenario'] == scenario) & (df['gsize'] == gene_size)]

                    if len(df_sub):
                        df_top = parameter_finder.get_default_run_values(
                            df_sub, table_key, default_params, method
                        )
                    else:
                        df_top = pd.DataFrame([[0, method]], columns=['value', 'method'])

                    # Apply abbreviations
                    df_top.loc[df_top['method'] == 'AffinityPropagation', 'method'] = 'AP'
                    df_top.loc[df_top['method'] == 'AgglomerativeClustering', 'method'] = 'AC'

                    scenario_data[scenario][gene_size] = pd.concat([
                        scenario_data[scenario][gene_size], df_top
                    ])

        except Exception as e:
            print(f"Warning: Could not process {method}: {e}")
            continue

    # Create figure with 3x3 subplot grid
    fig, axs = plt.subplots(
        nrows=len(settings.SCENARIOS),
        ncols=len(settings.GENE_SIZES),
        sharex=False,
        sharey=False,
        figsize=(16, 12),
        gridspec_kw={'height_ratios': [2, 2, 2]}
    )

    # Create method order
    order = plotting_utils.create_method_order(found_methods)

    # Plot each scenario/gene_size combination
    for x, scenario in enumerate(settings.SCENARIOS):
        for y, gene_size in enumerate(settings.GENE_SIZES):
            ax = axs[y, x]
            data = scenario_data[scenario][gene_size]

            # Format method names
            data = plotting_utils.format_method_names(data, 'method')

            # Create barplot
            hide_x_labels = (gene_size != 500)
            plotting_utils.create_simple_barplot(
                data=data,
                x_col='method',
                y_col='value',
                order=order,
                ax=ax,
                title=f'{scenario}, {gene_size} features',
                ylim=(0, 1),
                capsize=0.2,
                hide_x_labels=hide_x_labels
            )

    # Add y-axis label
    plotting_utils.add_ylabel(fig, 'Performance')

    # Add legend
    palette = copy.deepcopy(settings.METHOD_TYPE_PALETTE)
    plotting_utils.create_figure_legend(
        fig, palette, bbox_to_anchor=(0.25, 0.85)
    )

    sns.despine(bottom=True, left=True, right=True, top=True)

    # Save figure
    plotting_utils.save_figure(fig, 'supplement/FigureS3')
    print("Saved: paper/supplement/FigureS3.png")
    plt.close()


def generate_figure_s3():
    """
    Generate Figure S4: Performance increase (optimized - default).

    Original logic from notebook Cell 11.
    """
    print("\n=== Generating Figure S4 (Performance increase) ===")

    plotting_utils.setup_plotting_theme()

    # Initialize data structure
    scenario_data = {s: {g: pd.DataFrame(columns=['method', 'value'])
                         for g in settings.GENE_SIZES}
                     for s in settings.SCENARIOS}
    table_key = 'performance'
    found_methods = set()

    # Collect data for all methods
    for method in settings.METHODS:
        try:
            df = data_loader.load_simulated_data(method)
            found_methods.add(method)

            # Find both best and default parameters
            best_params = parameter_finder.calc_best_params(df, table_key)
            default_params = parameter_finder.find_default_params(df, table_key, method)
            print(f"{method}: best vs default")

            # Extract data for each scenario/gene_size combination
            for scenario in settings.SCENARIOS:
                for gene_size in settings.GENE_SIZES:
                    df_sub = df[(df['scenario'] == scenario) & (df['gsize'] == gene_size)]

                    if len(df_sub):
                        # Calculate delta: best - default
                        df_top = parameter_finder.calculate_performance_delta(
                            df_sub, table_key, best_params, default_params, method
                        )
                    else:
                        df_top = pd.DataFrame([[0, method]], columns=['value', 'method'])

                    # Apply abbreviations
                    df_top.loc[df_top['method'] == 'AffinityPropagation', 'method'] = 'AP'
                    df_top.loc[df_top['method'] == 'AgglomerativeClustering', 'method'] = 'AC'

                    scenario_data[scenario][gene_size] = pd.concat([
                        scenario_data[scenario][gene_size], df_top
                    ])

        except Exception as e:
            print(f"Warning: Could not process {method}: {e}")
            continue

    # Create figure with 3x3 subplot grid
    fig, axs = plt.subplots(
        nrows=len(settings.SCENARIOS),
        ncols=len(settings.GENE_SIZES),
        sharex=False,
        sharey=False,
        figsize=(16, 12),
        gridspec_kw={'height_ratios': [2, 2, 2]}
    )

    # Create method order
    order = plotting_utils.create_method_order(found_methods)

    # Plot each scenario/gene_size combination
    for x, scenario in enumerate(settings.SCENARIOS):
        for y, gene_size in enumerate(settings.GENE_SIZES):
            ax = axs[y, x]
            data = scenario_data[scenario][gene_size]

            # Format method names
            data = plotting_utils.format_method_names(data, 'method')

            # Create barplot
            hide_x_labels = (gene_size != 500)
            plotting_utils.create_simple_barplot(
                data=data,
                x_col='method',
                y_col='value',
                order=order,
                ax=ax,
                title=f'{scenario}, {gene_size} features',
                ylim=(0, 1),
                capsize=0.2,
                hide_x_labels=hide_x_labels
            )

    # Add y-axis label
    plotting_utils.add_ylabel(fig, 'Performance increase')

    # Add legend
    palette = copy.deepcopy(settings.METHOD_TYPE_PALETTE)
    plotting_utils.create_figure_legend(
        fig, palette, bbox_to_anchor=(0.25, 0.85)
    )

    sns.despine(bottom=True, left=True, right=True, top=True)

    # Save figure
    plotting_utils.save_figure(fig, 'supplement/FigureS4')
    print("Saved: paper/supplement/FigureS4.png")
    plt.close()


def generate_figure_s4_a(best_params):
    """
    Generate Figure S5a (TCGA) and Figure S5b (METABRIC): BRCA performance tuned on each dataset.

    Original logic from notebook Cell 13.

    Args:
        best_params: Dict of best params per dataset (from compute_best_params_per_dataset)
    """
    print("\n=== Generating Figure S5a (TCGA) and Figure S5b (METABRIC) ===")

    plotting_utils.setup_plotting_theme()

    table_key = 'PAM50'

    # Generate for each tuning dataset
    for tuned_on_dataset in settings.REAL_DATASETS:
        print(f"\nTuned on {tuned_on_dataset}:")

        data = {d: pd.DataFrame() for d in settings.REAL_DATASETS}
        found_methods = set()

        # Collect data
        for method in settings.METHODS:
            try:
                # Load data from BOTH datasets
                df_all = pd.DataFrame()
                for dataset in settings.REAL_DATASETS:
                    df_cur = data_loader.load_real_data(method, dataset)
                    df_cur['dataset'] = dataset
                    df_all = pd.concat([df_all, df_cur])

                found_methods.add(method)

                # Get best params for the dataset we're tuning on
                top_params = best_params[tuned_on_dataset][method]

                # Extract performance on both datasets with these params
                for dataset in settings.REAL_DATASETS:
                    df_filtered = df_all[df_all['parameters'] == top_params]
                    df_filtered = df_filtered[df_filtered['dataset'] == dataset]
                    df_filtered = data_loader.clean_columns(df_filtered)

                    # Melt and format
                    df_melted = df_filtered[[table_key]].melt()
                    df_melted['Method'] = method
                    df_melted = df_melted.rename(columns={
                        'value': 'Jaccard index',
                        'variable': 'Cancer type'
                    })

                    data[dataset] = pd.concat([data[dataset], df_melted])

            except Exception as e:
                print(f"  Warning: Could not process {method}: {e}")
                continue

        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(16, 6))

        # Combine data from both datasets
        df_plot = pd.DataFrame()
        for dataset in settings.REAL_DATASETS:
            df_sub = data[dataset].copy()
            df_sub['Dataset'] = dataset
            df_plot = pd.concat([df_plot, df_sub])

        # Remove network constraint methods
        df_plot = df_plot[~df_plot['Method'].isin(settings.METHODS_NETWORK_CONSTRAINT)]

        # Create method order
        order = plotting_utils.create_method_order(
            found_methods, filter_network_constraint=True
        )

        # Format method names
        df_plot = plotting_utils.format_method_names(df_plot, 'Method')
        df_plot['Jaccard index'] = pd.to_numeric(df_plot['Jaccard index'])

        # Rename TCGA to TCGA-BRCA
        df_plot = plotting_utils.rename_tcga_to_brca(df_plot)

        # Create grouped barplot
        sub_fig = plotting_utils.create_hue_barplot(
            data=df_plot,
            x_col='Method',
            y_col='Jaccard index',
            hue_col='Dataset',
            order=order,
            hue_order=['TCGA-BRCA', 'METABRIC'],
            ax=ax,
            ylim=(0, 1)
        )

        # Apply colors and hatching
        palette = plotting_utils.create_method_palette(order)
        plotting_utils.apply_hatch_pattern(sub_fig, palette)

        sns.despine(bottom=True, left=True, right=True, top=True)
        sub_fig.set(xlabel=None, ylabel=None)
        plotting_utils.add_ylabel(fig, 'Performance')

        # Save figure
        if tuned_on_dataset == 'METABRIC':
            filename = 'supplement/FigureS5b'
        else:  # TCGA
            filename = 'supplement/FigureS5a'
        plotting_utils.save_figure(fig, filename)
        print(f"Saved: paper/{filename}.png")
        plt.close()


def generate_figure4_a(best_average_rank_brca):
    """
    Generate Figure 3b: BRCA with best average rank parameters.

    Original logic from notebook Cell 15.

    Args:
        best_average_rank_brca: Dict of best average rank params (from compute_best_average_rank_brca)
    """
    print("\n=== Generating Figure 3b (BRCA best average rank) ===")

    plotting_utils.setup_plotting_theme()

    table_key = 'PAM50'
    found_methods = set()
    df_data_plot = pd.DataFrame()

    # Collect data
    for method in settings.METHODS:
        try:
            df_all_raw = pd.DataFrame()

            for dataset in settings.REAL_DATASETS:
                df_cur = data_loader.load_real_data(method, dataset)
                df_cur['Dataset'] = dataset
                df_cur['Method'] = method
                df_all_raw = pd.concat([df_all_raw, df_cur])

            found_methods.add(method)

            # Extract data using best_average_rank parameters
            best_params = best_average_rank_brca[method]
            for dataset in settings.REAL_DATASETS:
                df_sub = df_all_raw[df_all_raw['Dataset'] == dataset]
                df_sub = df_sub[df_sub['parameters'] == best_params]
                df_data_plot = pd.concat([df_data_plot, df_sub])

        except Exception as e:
            print(f"Warning: Could not process {method}: {e}")
            continue

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))

    # Remove network constraint methods
    df_data_plot = df_data_plot[
        ~df_data_plot['Method'].isin(settings.METHODS_NETWORK_CONSTRAINT)
    ]

    # Create method order
    order = plotting_utils.create_method_order(
        found_methods, filter_network_constraint=True
    )

    # Format method names
    df_data_plot = plotting_utils.format_method_names(df_data_plot, 'Method')
    df_data_plot[table_key] = pd.to_numeric(df_data_plot[table_key])

    # Rename TCGA to TCGA-BRCA
    df_data_plot = plotting_utils.rename_tcga_to_brca(df_data_plot)

    # Create grouped barplot
    sub_fig = plotting_utils.create_hue_barplot(
        data=df_data_plot,
        x_col='Method',
        y_col=table_key,
        hue_col='Dataset',
        order=order,
        hue_order=['TCGA-BRCA', 'METABRIC'],
        ax=ax,
        ylim=(0, 1)
    )

    # Apply colors and hatching
    palette = plotting_utils.create_method_palette(order)
    plotting_utils.apply_hatch_pattern(sub_fig, palette)

    sns.despine(bottom=True, left=True, right=True, top=True)
    sub_fig.set(xlabel=None, ylabel=None)
    plotting_utils.add_ylabel(fig, 'Performance')

    # Save figure
    plotting_utils.save_figure(fig, 'Figure3b')
    print("Saved: paper/Figure3b.png")
    plt.close()


def generate_figure4_c():
    """
    Generate Figure S5d: Violin plot of all method performances.

    Original logic from notebook Cell 24.
    """
    print("\n=== Generating Figure S5d (Violin plot) ===")

    plotting_utils.setup_plotting_theme()

    table_key = 'PAM50'
    found_methods = set()

    # Load ALL data (all parameter combinations) for each method/dataset
    data = {dataset: {} for dataset in settings.REAL_DATASETS}

    for method in settings.METHODS:
        for dataset in settings.REAL_DATASETS:
            try:
                df = data_loader.load_real_data(method, dataset)
                found_methods.add(method)

                # Apply abbreviations to method name
                m = method.replace('_', ' ')
                m = settings.METHOD_ABBREVIATIONS[m] if m in settings.METHOD_ABBREVIATIONS else m

                # Store all PAM50 values (all runs, all parameter combinations)
                data[dataset][m] = df[table_key].tolist()

            except Exception as e:
                print(f"Warning: Could not load {method} {dataset}: {e}")
                continue

    # Create method order
    order = plotting_utils.create_method_order(
        found_methods, filter_network_constraint=True
    )

    # Create plot DataFrame
    df_plot_all = pd.DataFrame()
    for dataset in settings.REAL_DATASETS:
        df_plot = plotting_utils.create_method_dataframe_for_violin(
            data[dataset], dataset
        )
        df_plot_all = pd.concat([df_plot_all, df_plot])

    # Special case: rename sparse PCA
    df_plot_all['variable'] = df_plot_all['variable'].replace('sparse PCA', 'sPCA')

    # Rename TCGA to TCGA-BRCA
    df_plot_all = plotting_utils.rename_tcga_to_brca(df_plot_all, 'dataset')

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))

    # Create violin plot
    sub_fig = plotting_utils.create_violinplot(
        data=df_plot_all,
        x_col='variable',
        y_col='value',
        hue_col='dataset',
        order=order,
        ax=ax,
        ylim=(-0.2, 1)
    )

    # Apply colors and hatching
    plotting_utils.apply_violin_hatch_pattern(sub_fig, order)

    plotting_utils.add_ylabel(fig, 'Performance')

    # Save figure
    plotting_utils.save_figure(fig, 'supplement/FigureS5d')
    print("Saved: paper/supplement/FigureS5d.png")
    plt.close()


def generate_figure_s5d(best_average_rank_brca, default_params):
    """
    Generate Figure S5c: BRCA delta (optimized - default) performance.

    Original logic from notebook Cell 23 (Figure4_B).

    Args:
        best_average_rank_brca: Dict of best average rank params
        default_params: Dict of default params per dataset
    """
    print("\n=== Generating Figure S5c (BRCA optimized - default) ===")

    plotting_utils.setup_plotting_theme()

    table_key = 'PAM50'
    data = {dataset: {'OPTIMIZED': pd.DataFrame(), 'DEFAULT': pd.DataFrame()}
            for dataset in settings.REAL_DATASETS}
    found_methods = set()

    # Collect data
    for dataset in settings.REAL_DATASETS:
        for method in settings.METHODS:
            try:
                df = data_loader.load_real_data(method, dataset)
                df['dataset'] = dataset
                found_methods.add(method)

                # Aggregate runs if needed
                if 'run' in df.columns or not df["parameters"].is_unique:
                    pass  # Keep all runs for proper indexing

                # Skip if no default params for this method
                if method not in default_params[dataset]:
                    if method == 'COALESCE':
                        # Special handling for COALESCE
                        df_top_values_sub = pd.DataFrame({
                            'Cancer type': {0: 'PAM50'},
                            'Jaccard index': {0: 0},
                            'Method': {0: 'COALESCE'}
                        })
                        data[dataset]['DEFAULT'] = pd.concat([data[dataset]['DEFAULT'], df_top_values_sub])
                        continue
                    else:
                        continue

                # Get optimized performance
                top_params = best_average_rank_brca[method]
                df_opt = df[df['parameters'] == top_params]
                df_opt = data_loader.clean_columns(df_opt)
                df_opt = df_opt[[table_key, 'dataset']]
                df_opt_sub = df_opt[df_opt['dataset'] == dataset]
                del df_opt_sub['dataset']
                df_opt_sub = df_opt_sub.melt()
                df_opt_sub['Method'] = method
                df_opt_sub = df_opt_sub.rename(columns={'value': 'Jaccard index', 'variable': 'Cancer type'})
                data[dataset]['OPTIMIZED'] = pd.concat([data[dataset]['OPTIMIZED'], df_opt_sub])

                # Get default performance (skip COALESCE)
                if method == 'COALESCE':
                    continue

                default_p = default_params[dataset][method]
                df_def = df[df['parameters'] == default_p]
                df_def = data_loader.clean_columns(df_def)
                df_def = df_def[[table_key, 'dataset']]
                df_def_sub = df_def[df_def['dataset'] == dataset]
                del df_def_sub['dataset']
                df_def_sub = df_def_sub.melt()
                df_def_sub['Method'] = method
                df_def_sub = df_def_sub.rename(columns={'value': 'Jaccard index', 'variable': 'Cancer type'})
                data[dataset]['DEFAULT'] = pd.concat([data[dataset]['DEFAULT'], df_def_sub])

            except Exception as e:
                print(f"  Warning: Could not process {method}: {e}")
                continue

    # Create plot
    fig, ax = plt.subplots(1, 1, figsize=(16, 6))
    df_plot = pd.DataFrame()

    for dataset in settings.REAL_DATASETS:
        # Reset index for proper alignment
        df_opt = data[dataset]['OPTIMIZED'].copy().reset_index()
        df_opt['index'] = df_opt['Method'] + df_opt['index'].map(str)
        df_opt = df_opt.set_index('index')

        df_def = data[dataset]['DEFAULT'].copy().reset_index()
        df_def['index'] = df_def['Method'] + df_def['index'].map(str)
        df_def = df_def.set_index('index')

        # Calculate delta
        df_opt['Jaccard index'] = df_opt['Jaccard index'] - df_def['Jaccard index']
        sub_df = df_opt

        # Remove network constraint methods
        sub_df = sub_df[~sub_df['Method'].isin(settings.METHODS_NETWORK_CONSTRAINT)]

        # Format method names
        sub_df = plotting_utils.format_method_names(sub_df, 'Method')
        sub_df['Jaccard index'] = pd.to_numeric(sub_df['Jaccard index'])
        sub_df['Dataset'] = dataset
        df_plot = pd.concat([df_plot, sub_df])

    # Add UnPaSt with 0 delta (no optimization)
    df_plot.loc['UnPaSt0'] = ['PAM50', 0, 'UnPaSt', 'TCGA']
    df_plot.loc['UnPaSt1'] = ['PAM50', 0, 'UnPaSt', 'METABRIC']

    # Rename TCGA to TCGA-BRCA
    df_plot = plotting_utils.rename_tcga_to_brca(df_plot)

    # Create method order
    order = plotting_utils.create_method_order(found_methods, filter_network_constraint=True)

    # Create barplot
    palette = ['grey', 'grey']
    sub_fig = sns.barplot(
        palette=palette,
        hue_order=['TCGA-BRCA', 'METABRIC'],
        hue='Dataset',
        data=df_plot,
        x='Method',
        y='Jaccard index',
        ax=ax,
        estimator=np.mean,
        errorbar=(lambda x: (min(x), max(x))),
        capsize=0.05,
        order=order
    )
    plt.legend(loc='upper right')

    # Apply colors and hatching
    method_palette = plotting_utils.create_method_palette(order)
    for bars, hatch, legend_handle in zip(sub_fig.containers, ['', '//'], sub_fig.legend_.legendHandles):
        for bar, color in zip(bars, method_palette):
            bar.set_facecolor(color)
            bar.set_hatch(hatch)
        legend_handle.set_hatch(hatch + hatch)

    ax.tick_params(axis='x', rotation=90)
    sns.despine(bottom=True, left=True, right=True, top=True)
    sub_fig.set(xlabel=None, ylabel=None)
    plotting_utils.add_ylabel(fig, 'Performance increase')
    ax.set_ylim(-0.2, 0.6)

    # Save figure
    plotting_utils.save_figure(fig, 'supplement/FigureS5c')
    print("Saved: paper/supplement/FigureS5c.png")
    plt.close()


def generate_figure_s6(best_average_rank_brca):
    """
    Generate Figure S6: Heatmap of performance across cancer types.

    Original logic from notebook Cell 29 (supplement/FigureS4).

    Args:
        best_average_rank_brca: Dict of best average rank params
    """
    print("\n=== Generating Figure S6 (Heatmap across cancer types) ===")

    plotting_utils.setup_plotting_theme()

    from matplotlib.patches import Rectangle

    table_key = 'PAM50'
    data = {dataset: pd.DataFrame() for dataset in settings.REAL_DATASETS}

    # Collect data
    for method in settings.METHODS:
        for dataset in settings.REAL_DATASETS:
            try:
                df = data_loader.load_real_data(method, dataset)
                found_methods.add(method) if 'found_methods' in dir() else None

                # Aggregate runs if needed
                if 'run' in df.columns or not df["parameters"].is_unique:
                    pass  # Keep for now

                # Get best average rank params
                top_params = best_average_rank_brca[method]
                df_top = df[df['parameters'] == top_params]
                df_top = data_loader.clean_columns(df_top)

                # Melt and format
                df_top = df_top.melt()
                df_top['Method'] = method
                df_top = df_top.rename(columns={'value': 'Jaccard index', 'variable': 'Cancer type'})

                data[dataset] = pd.concat([data[dataset], df_top])

            except Exception as e:
                print(f"  Warning: Could not process {method} {dataset}: {e}")
                continue

    # Get cancer types
    cancer_types = []
    for _, values in settings.CANCER_GROUPS:
        cancer_types.extend(values)

    # Get methods
    methods = data[dataset]['Method'].unique()
    methods = [settings.METHOD_ABBREVIATIONS[x] if x in settings.METHOD_ABBREVIATIONS else x for x in methods]

    # Create figure
    fig, axs = plt.subplots(
        1, len(settings.REAL_DATASETS) + 1,
        sharex=False, sharey=False,
        figsize=(12, 9),
        gridspec_kw={'width_ratios': [15, 15, 1]}
    )

    for y, scenario in enumerate(['TCGA', 'METABRIC']):
        plot_data_df = pd.DataFrame()

        for i, cancer_type in enumerate(cancer_types):
            data_sub = data[scenario][data[scenario]['Cancer type'] == cancer_type].copy()
            data_sub = plotting_utils.format_method_names(data_sub, 'Method')
            del data_sub['Cancer type']
            data_sub['Jaccard index'] = data_sub['Jaccard index'].astype(float)
            data_sub = data_sub.groupby('Method').mean()
            data_sub = data_sub.reindex(methods).reset_index(drop=True)
            data_sub['Method'] = [i for i in range(len(methods))]
            data_sub['Cancer type'] = cancer_type.replace('_', ' ')
            plot_data_df = pd.concat([plot_data_df, data_sub])

        # Pivot and find max indices
        data_pivot = plot_data_df.pivot(index='Method', columns='Cancer type')['Jaccard index']
        order = [x.replace('_', ' ') for x in cancer_types]
        data_pivot = data_pivot[order]

        max_indices_list = []
        for cancer_type in data_pivot.columns:
            max_value = data_pivot[cancer_type].max()
            max_indices = data_pivot[cancer_type][data_pivot[cancer_type] == max_value].index.tolist()
            max_indices_list.append(max_indices)

        data_pivot = data_pivot.applymap(lambda x: round(x, 2))

        # Create heatmap
        yticklabels = methods if y == 0 else False
        sub_fig = sns.heatmap(
            data_pivot, annot=True,
            yticklabels=yticklabels, ax=axs[y],
            cmap="Blues", cbar=False,
            annot_kws={"size": 9},
            linecolor='white', linewidth=0,
            vmin=0, vmax=1
        )

        scenario_label = 'TCGA-BRCA' if scenario == 'TCGA' else scenario
        axs[y].set_title(f'{scenario_label}', fontsize=18)
        axs[y].tick_params(axis='x', rotation=90)

        # Highlight max values
        for col, col_vals in enumerate(max_indices_list):
            for ind in col_vals:
                axs[y].add_patch(Rectangle((col, ind), 1, 1, fill=False, edgecolor='red', lw=2))

        sub_fig.set(xlabel=None, ylabel=None)
        sns.despine(bottom=True, left=True, right=True, top=True)

        # Color tick labels by method type
        colors = [settings.METHOD_PALETTE[x] for x in methods]
        for ticklabel, tickcolor in zip(axs[y].get_yticklabels(), colors):
            ticklabel.set_color(tickcolor)

    # Add colorbar
    sm = cm.ScalarMappable(cmap='Blues')
    sm.set_array([])
    sm.set_clim(vmin=0, vmax=1)

    colorbar = plt.colorbar(sm, cax=axs[2])
    colorbar.set_ticks([0, .2, .4, .6, .8, 1])
    colorbar.set_ticklabels(['<=0', .2, .4, .6, .8, 1])
    colorbar.set_label('Performance')

    sns.despine(bottom=True, left=True, right=True, top=True)

    # Save figure
    plotting_utils.save_figure(fig, 'supplement/FigureS6')
    print("Saved: paper/supplement/FigureS6.png")
    plt.close()


def main():
    """Main function to generate all figures."""
    print("=" * 80)
    print("PAPER FIGURE GENERATION")
    print("=" * 80)

    # Check for command line arguments
    requested_figure = sys.argv[1] if len(sys.argv) > 1 else None

    # Compute prerequisite data structures
    print("\n=== Computing prerequisite parameters ===")
    best_average_rank_brca = compute_best_average_rank_brca()
    best_params = compute_best_params_per_dataset()
    default_params = compute_default_params_per_dataset()

    # Generate figures based on request
    if requested_figure is None or requested_figure == 's3':
        generate_figure_s3_1()

    if requested_figure is None or requested_figure == 's4':
        generate_figure_s3()

    if requested_figure is None or requested_figure == 's4_a':
        generate_figure_s4_a(best_params)

    if requested_figure is None or requested_figure == 'figure3b':
        generate_figure4_a(best_average_rank_brca)

    if requested_figure is None or requested_figure == 's5c':
        generate_figure_s5d(best_average_rank_brca, default_params)

    if requested_figure is None or requested_figure == 's5d':
        generate_figure4_c()

    if requested_figure is None or requested_figure == 's6':
        generate_figure_s6(best_average_rank_brca)

    print("\n" + "=" * 80)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 80)


if __name__ == '__main__':
    main()

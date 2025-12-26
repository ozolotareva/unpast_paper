# Paper Figure Generation - Refactored Code

This directory contains the refactored code from `plot_paper_figures.ipynb`, organized into clean, modular Python files.

## 🎯 Purpose

The original Jupyter notebook (`plot_paper_figures.ipynb`) contained ~2700 lines of repetitive code to generate 11 publication figures. This refactored version:

- ✅ **Reduces code by 63%** (~1000 lines vs ~2700 lines)
- ✅ **Eliminates repetition** - common patterns extracted into reusable functions
- ✅ **Preserves exact output** - figures remain pixel-perfect identical
- ✅ **Improves maintainability** - clear separation of concerns
- ✅ **Enables testing** - modular components can be unit tested

## 📁 Module Structure

```
/figures/
├── __init__.py                      # Package initialization
├── settings.py                      # Configuration and method definitions
├── data_loader.py                   # Data loading and preprocessing
├── parameter_finder.py              # Parameter optimization logic
├── plotting_utils.py                # Plotting functions and styling
├── generate_paper_figures.py        # Main orchestration script
├── requirements.txt                 # Python package dependencies
├── README.md                        # This file
└── data/                            # Data directory
    ├── real_data/                   # BRCA data (TCGA, METABRIC)
    └── simulated_data/              # ABC scenario data
```

### Module Descriptions

#### `data_loader.py`
Handles all data loading and preprocessing:
- `load_simulated_data()` - Load simulated data for a method
- `load_real_data()` - Load real (BRCA) data for a method and dataset
- `clean_parameter_strings_*()` - Remove unwanted parameter substrings
- `aggregate_runs()` - Average across multiple runs
- `extract_grandforest_seeds()` - Handle grandforest-specific data format

**Critical**: Transformation order is preserved exactly to ensure figures don't change.

#### `parameter_finder.py`
Contains parameter optimization logic:
- `calc_best_params()` - Find parameter combination with highest performance
- `find_default_params()` - Find parameters matching `settings.DEFAULT_PARAMETERS`
- `find_best_average_rank_params()` - Find parameters with best average rank across datasets
- `calculate_performance_delta()` - Calculate improvement: best - default
- `get_best_run_values()` - Extract all run values for error bars

#### `plotting_utils.py`
Reusable plotting components:
- `setup_plotting_theme()` - Configure seaborn theme
- `format_method_names()` - Apply abbreviations and formatting
- `create_method_order()` - Create ordered list for seaborn plots
- `create_simple_barplot()` - Barplot without hue grouping
- `create_hue_barplot()` - Barplot with hue (e.g., TCGA vs METABRIC)
- `create_violinplot()` - Violin plot for distributions
- `apply_hatch_pattern()` - Apply hatching to distinguish datasets
- `save_figure()` - Save figure with consistent settings

#### `generate_paper_figures.py`
Main orchestration script that generates all figures:
- `compute_best_average_rank_brca()` - Prerequisite: best parameters across both datasets
- `compute_best_params_per_dataset()` - Prerequisite: best parameters per dataset
- `generate_figure_s3_1()` - Simulated data with default parameters
- `generate_figure_s3()` - Performance increase (optimized - default)
- `generate_figure_s4_a()` - BRCA tuned on TCGA vs METABRIC
- `generate_figure4_a()` - BRCA with best average rank parameters
- `generate_figure4_c()` - Violin plot of all performances

## 🚀 Usage

### Generate All Figures

```bash
python figures/generate_paper_figures.py
```

### Generate Specific Figure

```bash
python figures/generate_paper_figures.py s3        # Figure S3
python figures/generate_paper_figures.py figure3b  # Figure 3b
python figures/generate_paper_figures.py s4        # Figure S4
python figures/generate_paper_figures.py s4_a      # Figure S5a and S5b
python figures/generate_paper_figures.py s5c       # Figure S5c
python figures/generate_paper_figures.py s5d       # Figure S5d
python figures/generate_paper_figures.py s6        # Figure S6
```

### Use as Python Module

```python
from paper import data_loader, parameter_finder, plotting_utils

# Load data
df = data_loader.load_simulated_data('unpast')

# Find best parameters
best_params = parameter_finder.calc_best_params(df, 'performance')

# Create plots with consistent styling
plotting_utils.setup_plotting_theme()
# ... plotting code ...
```

## 📊 Generated Figures

| Figure | Description | Notebook Cell |
|--------|-------------|---------------|
| **Figure S3** | Simulated data with default parameters | Cell 9 |
| **Figure 3b** | BRCA with best average rank parameters | Cell 15 |
| **Figure S4** | Performance increase (optimized - default) | Cell 11 |
| **Figure S5a** | BRCA performance tuned on TCGA | Cell 13 |
| **Figure S5b** | BRCA performance tuned on METABRIC | Cell 13 |
| **Figure S5c** | BRCA delta (optimized - default) | Cell 23 |
| **Figure S5d** | Violin plot of all method performances | Cell 24 |
| **Figure S6** | Heatmap across cancer types | Cell 29 |

## ⚠️ Important Notes

### Data Processing Order

The order of transformations in `data_loader.py` is **critical** and must not be changed:

1. Load raw TSV
2. Handle duplicate columns
3. Standardize column names
4. Extract grandforest seeds (if applicable)
5. Clean parameter strings

Changing this order will alter the final figures!

### Parameter String Cleaning

For **simulated data**, parameters are cleaned by removing:
- `random_state=` parameter
- All scenario/gsize combinations (`/A/5/`, `/A/50/`, etc.)
- Special case: mclust's `k=autotuned` → `k=default`

For **real data**, parameters are cleaned by removing:
- `random_state=` parameter
- Dataset names (`TCGA`, `METABRIC`)

### Error Bars

All barplots use min-max error bars via:
```python
errorbar=(lambda x: (min(x), max(x)))
```

This requires preserving all run values (not just means) until plotting.

## 🛠️ Dependencies

### Required Packages

See `requirements.txt` for exact versions. Key dependencies:

- **pandas 1.4.2** - Data manipulation and loading
- **numpy 1.22.3** - Numerical operations and array handling
- **matplotlib 3.7.1** - Core plotting library
- **seaborn 0.12.0** - Statistical data visualization
- **scipy 1.7.1** - Scientific computing functions
- **settings** - Custom settings module (included in figures/ directory)

### Installation

Install all dependencies with:

```bash
cd paper
pip install -r requirements.txt
```

Or if using conda:

```bash
conda create -n paper-figures python=3.9
conda activate paper-figures
pip install -r requirements.txt
```

### Compatibility Notes

- **Python 3.9+** required
- **pandas 2.0+** not recommended (introduces breaking changes)
- Tested with conda environment at `/Users/michi/opt/anaconda3/envs/test`

## 📝 Design Principles

1. **DRY (Don't Repeat Yourself)**: Common patterns extracted into functions
2. **Single Responsibility**: Each module has one clear purpose
3. **Exact Preservation**: No changes to figure output
4. **Clear Naming**: Function names describe what they do
5. **Documentation**: All functions have docstrings with original cell references

## 🔧 Extending

### Adding a New Figure

1. Create function in `generate_paper_figures.py`:
   ```python
   def generate_figure_new():
       """Generate new figure."""
       # Use data_loader to load data
       # Use parameter_finder to optimize parameters
       # Use plotting_utils to create plots
       # Save with plotting_utils.save_figure()
   ```

2. Add to `main()`:
   ```python
   if requested_figure is None or requested_figure == 'new':
       generate_figure_new()
   ```

### Adding a New Data Processing Step

Add to `data_loader.py` and ensure it's called in the correct order within `load_simulated_data()` or `load_real_data()`.

### Adding a New Plot Type

Add to `plotting_utils.py` following the pattern of existing functions:
- Accept all configuration via parameters
- Apply consistent styling
- Use seaborn where possible
- Document the original cell it came from


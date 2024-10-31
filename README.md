# qpcr-analyzer
A comprehensive Python package for analyzing qPCR data with support for multiple genes, samples, and flexible experimental designs.


## Features

- 📊 Handles both CSV and Excel input files
- 🧬 Supports multiple target genes and reference genes
- 📈 Automatic ΔCT and fold change calculations
- 📉 Publication-quality visualizations
- 📋 Comprehensive statistical analysis
- 🔍 Flexible group comparisons
- 📁 Easy-to-use data export

  ## Quick Start

```python
from qpcr_analyzer import QPCRAnalyzer

# Initialize analyzer
analyzer = QPCRAnalyzer(
    file_path='your_data.csv',
    target_genes=['Gene1', 'Gene2', 'Gene3'],
    reference_gene='GAPDH'
)

# Get statistics
stats = analyzer.get_statistics()

# Create visualizations
analyzer.plot_results()

# Save results
analyzer.save_results()
```

## Input Data Format

Your input file (CSV or Excel) should have the following columns:
- Sample: Sample identifiers
- Group: Group/condition labels
- Target gene columns: CT values for each target gene
- Reference gene column: CT values for reference gene

Example:
```csv
Sample,Group,Gene1,Gene2,Gene3,GAPDH
Sample1_A,Treatment1,25.3,27.5,23.4,20.0
Sample1_B,Treatment1,25.2,27.4,23.5,20.1
...
```

## Features in Detail

### 1. Data Analysis
- Automatic handling of technical replicates
- ΔCT calculation
- ΔΔCT and fold change calculations
- Statistical testing (t-tests)
- Multiple group comparisons
- Handles missing/undetermined values

### 2. Visualizations
- Bar plots with error bars
- Heatmaps showing expression patterns
- Violin plots with individual data points
- Statistical significance indicators
- Publication-ready figure formatting

### 3. Results Export
- Comprehensive Excel workbook with multiple sheets
- Raw data and calculated values
- Statistical analysis results
- Significant changes summary

## Advanced Usage

### Excluding Groups
```python
analyzer = QPCRAnalyzer(
    file_path='data.csv',
    target_genes=['Gene1', 'Gene2'],
    reference_gene='GAPDH',
    exclude_groups=['UnwantedGroup']
)
```

### Custom Output Directory
```python
# Specify output directory for plots
analyzer.plot_results(save_dir='path/to/output')

# Custom filename for results
analyzer.save_results('custom_results.xlsx')
```

## Sample Data

The repository includes a sample dataset (`sample_qpcr_data.csv`) demonstrating the required format:
- 3 target genes (Gene1, Gene2, Gene3)
- 1 reference gene (GAPDH)
- 4 groups with 3 replicates each
- 12 total samples

## Requirements

- Python 3.7+
- pandas
- numpy
- matplotlib
- seaborn
- scipy

## Citation

If you use this package in your research, please cite:
```
QPCRAnalyzer (2024)
https://github.com/alinamvarr/qpcr-analyzer
```

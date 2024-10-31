# File: qpcr_analyzer.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import List, Dict, Union, Tuple
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class QPCRAnalyzer:
    """
    A comprehensive qPCR data analysis tool.
    
    Features:
    - Handles both CSV and Excel files
    - Supports multiple genes and samples
    - Flexible reference gene selection
    - Publication-quality visualizations
    - Comprehensive statistical analysis
    """
    
    def __init__(self, 
                 file_path: str,
                 target_genes: List[str],
                 reference_gene: str,
                 sample_column: str = 'Sample',
                 group_column: str = 'Group',
                 control_group: str = 'Control',
                 exclude_groups: List[str] = None):
        """
        Initialize the QPCRAnalyzer.
        
        Args:
            file_path: Path to CSV or Excel file
            target_genes: List of target gene names
            reference_gene: Name of reference gene
            sample_column: Name of sample column
            group_column: Name of group column
            control_group: Name of control group
            exclude_groups: List of groups to exclude from analysis
        """
        self.file_path = Path(file_path)
        self.target_genes = target_genes
        self.reference_gene = reference_gene
        self.sample_column = sample_column
        self.group_column = group_column
        self.control_group = control_group
        self.exclude_groups = exclude_groups or []
        
        # Load and process data
        self.load_data()
        self.process_data()
        
    def load_data(self):
        """Load data from CSV or Excel file"""
        if self.file_path.suffix.lower() == '.csv':
            self.raw_data = pd.read_csv(self.file_path)
        else:
            self.raw_data = pd.read_excel(self.file_path)
            
        # Remove excluded groups
        if self.exclude_groups:
            self.raw_data = self.raw_data[~self.raw_data[self.group_column].isin(self.exclude_groups)]
            
        # Replace 'Undetermined' with NaN
        self.raw_data = self.raw_data.replace('Undetermined', np.nan)
        
        # Convert CT values to float
        for gene in self.target_genes + [self.reference_gene]:
            self.raw_data[gene] = pd.to_numeric(self.raw_data[gene], errors='coerce')
    
    def process_data(self):
        """Calculate ΔCT, ΔΔCT, and fold changes"""
        self.results = pd.DataFrame()
        
        # Calculate ΔCT for each target gene
        for gene in self.target_genes:
            dct = self.raw_data[gene] - self.raw_data[self.reference_gene]
            control_dct = dct[self.raw_data[self.group_column] == self.control_group].mean()
            ddct = dct - control_dct
            fold_change = np.power(2, -ddct)
            
            # Store results
            temp_df = pd.DataFrame({
                'Sample': self.raw_data[self.sample_column],
                'Group': self.raw_data[self.group_column],
                'Gene': gene,
                'CT': self.raw_data[gene],
                'DCT': dct,
                'DDCT': ddct,
                'Fold_Change': fold_change
            })
            self.results = pd.concat([self.results, temp_df], ignore_index=True)
    
    def get_statistics(self) -> pd.DataFrame:
        """Calculate comprehensive statistics"""
        stats_list = []
        
        for gene in self.target_genes:
            for group in self.results['Group'].unique():
                group_data = self.results[
                    (self.results['Gene'] == gene) & 
                    (self.results['Group'] == group)
                ]
                
                control_data = self.results[
                    (self.results['Gene'] == gene) & 
                    (self.results['Group'] == self.control_group)
                ]['Fold_Change']
                
                # Calculate statistics
                _, p_value = stats.ttest_ind(
                    group_data['Fold_Change'].dropna(),
                    control_data.dropna()
                )
                
                stats_list.append({
                    'Gene': gene,
                    'Group': group,
                    'N': len(group_data),
                    'CT_Mean': group_data['CT'].mean(),
                    'CT_SEM': group_data['CT'].sem(),
                    'DCT_Mean': group_data['DCT'].mean(),
                    'DCT_SEM': group_data['DCT'].sem(),
                    'Fold_Change_Mean': group_data['Fold_Change'].mean(),
                    'Fold_Change_SEM': group_data['Fold_Change'].sem(),
                    'P_Value': p_value,
                    'Significant': p_value < 0.05
                })
        
        return pd.DataFrame(stats_list)
    
    def plot_results(self, save_dir: str = None):
        """Create publication-quality visualizations"""
        save_dir = Path(save_dir) if save_dir else self.file_path.parent
        
        # 1. Bar plot with error bars
        plt.figure(figsize=(12, 6))
        for i, gene in enumerate(self.target_genes):
            plt.subplot(1, len(self.target_genes), i+1)
            
            stats_data = self.get_statistics()
            gene_data = stats_data[stats_data['Gene'] == gene]
            
            x = range(len(gene_data))
            plt.bar(x, gene_data['Fold_Change_Mean'], 
                   yerr=gene_data['Fold_Change_SEM'],
                   capsize=5)
            
            plt.xticks(x, gene_data['Group'], rotation=45, ha='right')
            plt.ylabel('Fold Change (2^-ΔΔCT)')
            plt.title(f'{gene} Expression')
            
            # Add significance stars
            for j, row in gene_data.iterrows():
                if row['Significant'] and row['Group'] != self.control_group:
                    plt.text(j, row['Fold_Change_Mean'] + row['Fold_Change_SEM'],
                           '*', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(save_dir / 'fold_changes.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # 2. Heatmap
        plt.figure(figsize=(10, 8))
        pivot_data = self.results.pivot_table(
            values='Fold_Change',
            index='Group',
            columns='Gene',
            aggfunc='mean'
        )
        
        sns.heatmap(pivot_data, annot=True, cmap='RdBu_r', center=1,
                   fmt='.2f', cbar_kws={'label': 'Fold Change'})
        plt.title('Expression Heatmap')
        plt.tight_layout()
        plt.savefig(save_dir / 'heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # 3. Violin plots
        plt.figure(figsize=(12, 6))
        for i, gene in enumerate(self.target_genes):
            plt.subplot(1, len(self.target_genes), i+1)
            
            gene_data = self.results[self.results['Gene'] == gene]
            sns.violinplot(data=gene_data, x='Group', y='Fold_Change')
            sns.swarmplot(data=gene_data, x='Group', y='Fold_Change', 
                         color='black', size=4, alpha=0.6)
            
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('Fold Change (2^-ΔΔCT)')
            plt.title(f'{gene} Expression Distribution')
        
        plt.tight_layout()
        plt.savefig(save_dir / 'distribution.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    def save_results(self, output_file: str = None):
        """Save comprehensive results to Excel file"""
        if output_file is None:
            output_file = self.file_path.with_name(
                self.file_path.stem + '_analysis_results.xlsx'
            )
        
        with pd.ExcelWriter(output_file) as writer:
            # Raw results
            self.results.to_excel(writer, sheet_name='Raw_Data', index=False)
            
            # Statistics
            self.get_statistics().to_excel(writer, sheet_name='Statistics', index=False)
            
            # Summary of significant changes
            significant_changes = self.get_statistics()[
                (self.get_statistics()['Significant']) & 
                (self.get_statistics()['Group'] != self.control_group)
            ]
            significant_changes.to_excel(writer, sheet_name='Significant_Changes', index=False)

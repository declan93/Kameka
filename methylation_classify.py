#!/usr/bin/env python3
"""
Methylation State Classification Tool
Classifies genomic regions based on methylation patterns and integrates with genomic annotations.
"""

import pandas as pd
import numpy as np
import argparse
import logging
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

class MethylationClassifier:
    """
    Comprehensive methylation state classifier with multiple classification methods
    """
    
    def __init__(self, hypo_threshold: float = 0.3, hyper_threshold: float = 0.7):
        self.hypo_threshold = hypo_threshold
        self.hyper_threshold = hyper_threshold
        self.logger = self._setup_logging()
        
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('methylation_classify.log'),
                logging.StreamHandler()
            ]
        )
        return logging.getLogger(__name__)
    
    def load_methylation_data(self, filepath: str, format_type: str = 'bedgraph') -> pd.DataFrame:
        """
        Load methylation data from various formats
        
        Args:
            filepath: Path to methylation data file
            format_type: Format type ('bedgraph', 'bismark', 'custom')
            
        Returns:
            DataFrame with columns: chr, start, end, methylation_level, coverage
        """
        self.logger.info(f"Loading methylation data from {filepath}")
        
        if format_type == 'bedgraph':
            df = pd.read_csv(filepath, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'methylation_level'])
            df['coverage'] = np.nan  # Not available in bedgraph
            
        elif format_type == 'bismark':
            df = pd.read_csv(filepath, sep='\t')
            # Assume standard Bismark format
            df = df.rename(columns={
                'chr': 'chr', 'pos': 'start', 'methylation_percentage': 'methylation_level',
                'count_methylated': 'meth_count', 'count_total': 'coverage'
            })
            df['end'] = df['start'] + 1
            df['methylation_level'] = df['methylation_level'] / 100.0
            
        elif format_type == 'custom':
            df = pd.read_csv(filepath, sep='\t')
            # Assume custom format with required columns
            
        else:
            raise ValueError(f"Unsupported format: {format_type}")
            
        # Validate and clean data
        df = df.dropna(subset=['methylation_level'])
        df = df[(df['methylation_level'] >= 0) & (df['methylation_level'] <= 1)]
        
        self.logger.info(f"Loaded {len(df)} methylation sites")
        return df
    
    def classify_simple_threshold(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Simple threshold-based classification
        """
        df = df.copy()
        conditions = [
            df['methylation_level'] < self.hypo_threshold,
            (df['methylation_level'] >= self.hypo_threshold) & 
            (df['methylation_level'] <= self.hyper_threshold),
            df['methylation_level'] > self.hyper_threshold
        ]
        choices = ['hypomethylated', 'intermediate', 'hypermethylated']
        df['methylation_class'] = np.select(conditions, choices)
        
        return df
    
    def classify_adaptive_threshold(self, df: pd.DataFrame, 
                                  percentile_low: float = 25, 
                                  percentile_high: float = 75) -> pd.DataFrame:
        """
        Adaptive threshold based on data distribution
        """
        df = df.copy()
        low_thresh = np.percentile(df['methylation_level'], percentile_low)
        high_thresh = np.percentile(df['methylation_level'], percentile_high)
        
        conditions = [
            df['methylation_level'] < low_thresh,
            (df['methylation_level'] >= low_thresh) & (df['methylation_level'] <= high_thresh),
            df['methylation_level'] > high_thresh
        ]
        choices = ['hypomethylated', 'intermediate', 'hypermethylated']
        df['methylation_class'] = np.select(conditions, choices)
        
        self.logger.info(f"Adaptive thresholds: {low_thresh:.3f} - {high_thresh:.3f}")
        return df
    
    def classify_kmeans(self, df: pd.DataFrame, n_clusters: int = 3) -> pd.DataFrame:
        """
        K-means clustering-based classification
        """
        df = df.copy()
        
        # Prepare features (can be extended with more features)
        features = df[['methylation_level']].values
        if 'coverage' in df.columns and not df['coverage'].isna().all():
            # Include coverage if available
            scaler = StandardScaler()
            features = scaler.fit_transform(df[['methylation_level', 'coverage']].fillna(0))
        
        # K-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(features)
        
        # Map clusters to methylation states based on mean methylation levels
        cluster_means = []
        for i in range(n_clusters):
            cluster_mean = df.loc[clusters == i, 'methylation_level'].mean()
            cluster_means.append((i, cluster_mean))
        
        # Sort by mean methylation level
        cluster_means.sort(key=lambda x: x[1])
        
        # Create mapping
        if n_clusters == 3:
            state_mapping = {
                cluster_means[0][0]: 'hypomethylated',
                cluster_means[1][0]: 'intermediate', 
                cluster_means[2][0]: 'hypermethylated'
            }
        else:
            # For other numbers of clusters, use percentile-based naming
            state_mapping = {}
            for idx, (cluster_id, _) in enumerate(cluster_means):
                percentile = (idx + 1) * 100 // n_clusters
                state_mapping[cluster_id] = f'cluster_{percentile}th_percentile'
        
        df['methylation_class'] = [state_mapping[c] for c in clusters]
        
        return df
    
    def add_genomic_context(self, df: pd.DataFrame, annotation_files: Dict[str, str]) -> pd.DataFrame:
        """
        Add genomic context annotations (CpG islands, promoters, etc.)
        
        Args:
            df: Methylation DataFrame
            annotation_files: Dict mapping annotation type to file path
                              e.g., {'cpg_islands': 'cpg.bed', 'promoters': 'promoters.bed'}
        """
        df = df.copy()
        
        # Convert to BedTool
        bed_cols = ['chr', 'start', 'end']
        methylation_bed = pybedtools.BedTool.from_dataframe(df[bed_cols])
        
        for annotation_type, annotation_file in annotation_files.items():
            self.logger.info(f"Adding {annotation_type} annotations")
            
            try:
                annotation_bed = pybedtools.BedTool(annotation_file)
                
                # Intersect with annotation
                intersected = methylation_bed.intersect(annotation_bed, wa=True, wb=True)
                
                # Convert back to DataFrame
                if len(intersected) > 0:
                    intersect_df = intersected.to_dataframe()
                    
                    # Create boolean column for this annotation
                    df[f'in_{annotation_type}'] = False
                    
                    # Mark intersecting regions
                    for _, row in intersect_df.iterrows():
                        mask = ((df['chr'] == row['chrom']) & 
                               (df['start'] >= row['start']) & 
                               (df['end'] <= row['end']))
                        df.loc[mask, f'in_{annotation_type}'] = True
                else:
                    df[f'in_{annotation_type}'] = False
                    
            except Exception as e:
                self.logger.warning(f"Could not process {annotation_type}: {e}")
                df[f'in_{annotation_type}'] = False
        
        return df
    
    def calculate_regional_methylation(self, df: pd.DataFrame, 
                                     window_size: int = 1000) -> pd.DataFrame:
        """
        Calculate regional methylation averages in sliding windows
        """
        regional_data = []
        
        for chrom in df['chr'].unique():
            chrom_data = df[df['chr'] == chrom].sort_values('start')
            
            if len(chrom_data) == 0:
                continue
                
            min_pos = chrom_data['start'].min()
            max_pos = chrom_data['end'].max()
            
            for window_start in range(min_pos, max_pos, window_size):
                window_end = window_start + window_size
                
                # Get sites in this window
                window_sites = chrom_data[
                    (chrom_data['start'] >= window_start) & 
                    (chrom_data['end'] <= window_end)
                ]
                
                if len(window_sites) > 0:
                    # Calculate regional statistics
                    regional_meth = window_sites['methylation_level'].mean()
                    meth_std = window_sites['methylation_level'].std()
                    site_count = len(window_sites)
                    
                    # Coverage-weighted average if coverage is available
                    if 'coverage' in window_sites.columns and not window_sites['coverage'].isna().all():
                        weights = window_sites['coverage'].fillna(1)
                        weighted_meth = np.average(window_sites['methylation_level'], weights=weights)
                    else:
                        weighted_meth = regional_meth
                    
                    regional_data.append({
                        'chr': chrom,
                        'start': window_start,
                        'end': window_end,
                        'methylation_level': regional_meth,
                        'weighted_methylation': weighted_meth,
                        'methylation_std': meth_std,
                        'site_count': site_count
                    })
        
        regional_df = pd.DataFrame(regional_data)
        
        if len(regional_df) > 0:
            # Classify regional methylation
            regional_df = self.classify_simple_threshold(regional_df)
        
        return regional_df
    
    def generate_summary_stats(self, df: pd.DataFrame) -> Dict:
        """
        Generate comprehensive summary statistics
        """
        stats_dict = {}
        
        # Overall distribution
        stats_dict['total_sites'] = len(df)
        stats_dict['methylation_mean'] = df['methylation_level'].mean()
        stats_dict['methylation_median'] = df['methylation_level'].median()
        stats_dict['methylation_std'] = df['methylation_level'].std()
        
        # Class distribution
        if 'methylation_class' in df.columns:
            class_counts = df['methylation_class'].value_counts()
            for class_name, count in class_counts.items():
                stats_dict[f'{class_name}_count'] = count
                stats_dict[f'{class_name}_fraction'] = count / len(df)
        
        # Genomic context distribution
        context_cols = [col for col in df.columns if col.startswith('in_')]
        for col in context_cols:
            context_name = col.replace('in_', '')
            if df[col].dtype == bool:
                stats_dict[f'{context_name}_sites'] = df[col].sum()
                stats_dict[f'{context_name}_fraction'] = df[col].mean()
        
        # Per-chromosome statistics
        chrom_stats = df.groupby('chr').agg({
            'methylation_level': ['count', 'mean', 'std']
        }).round(3)
        stats_dict['per_chromosome'] = chrom_stats.to_dict()
        
        return stats_dict
    
    def plot_methylation_distribution(self, df: pd.DataFrame, output_path: str):
        """
        Create comprehensive methylation distribution plots
        """
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Overall distribution
        axes[0,0].hist(df['methylation_level'], bins=50, alpha=0.7, edgecolor='black')
        axes[0,0].axvline(self.hypo_threshold, color='red', linestyle='--', label=f'Hypo threshold ({self.hypo_threshold})')
        axes[0,0].axvline(self.hyper_threshold, color='blue', linestyle='--', label=f'Hyper threshold ({self.hyper_threshold})')
        axes[0,0].set_xlabel('Methylation Level')
        axes[0,0].set_ylabel('Frequency')
        axes[0,0].set_title('Overall Methylation Distribution')
        axes[0,0].legend()
        
        # Classification pie chart
        if 'methylation_class' in df.columns:
            class_counts = df['methylation_class'].value_counts()
            axes[0,1].pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%')
            axes[0,1].set_title('Methylation Class Distribution')
        
        # Box plot by class
        if 'methylation_class' in df.columns:
            df.boxplot(column='methylation_level', by='methylation_class', ax=axes[1,0])
            axes[1,0].set_title('Methylation Level by Class')
            axes[1,0].set_xlabel('Methylation Class')
            axes[1,0].set_ylabel('Methylation Level')
        
        # Genomic context comparison
        context_cols = [col for col in df.columns if col.startswith('in_')]
        if context_cols:
            context_data = []
            for col in context_cols:
                context_name = col.replace('in_', '')
                in_context = df[df[col] == True]['methylation_level']
                out_context = df[df[col] == False]['methylation_level']
                
                context_data.extend([
                    {'Context': f'In {context_name}', 'Methylation': val} for val in in_context
                ])
                context_data.extend([
                    {'Context': f'Out {context_name}', 'Methylation': val} for val in out_context
                ])
            
            if context_data:
                context_df = pd.DataFrame(context_data)
                sns.boxplot(data=context_df, x='Context', y='Methylation', ax=axes[1,1])
                axes[1,1].set_title('Methylation by Genomic Context')
                axes[1,1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Methylation distribution plots saved to {output_path}")
    
    def save_results(self, df: pd.DataFrame, output_path: str, 
                    summary_stats: Dict, format_type: str = 'bed'):
        """
        Save classified methylation data and summary statistics
        """
        # Save main results
        if format_type == 'bed':
            # BED format output
            bed_cols = ['chr', 'start', 'end', 'methylation_level', 'methylation_class']
            available_cols = [col for col in bed_cols if col in df.columns]
            df[available_cols].to_csv(output_path, sep='\t', index=False, header=False)
        else:
            # Tab-separated format with header
            df.to_csv(output_path, sep='\t', index=False)
        
        # Save summary statistics
        summary_path = output_path.replace('.bed', '_summary.txt').replace('.tsv', '_summary.txt')
        with open(summary_path, 'w') as f:
            f.write("Methylation Classification Summary\n")
            f.write("=" * 50 + "\n\n")
            
            for key, value in summary_stats.items():
                if isinstance(value, dict):
                    f.write(f"{key}:\n")
                    for subkey, subvalue in value.items():
                        f.write(f"  {subkey}: {subvalue}\n")
                else:
                    f.write(f"{key}: {value}\n")
        
        self.logger.info(f"Results saved to {output_path}")
        self.logger.info(f"Summary statistics saved to {summary_path}")

def main():
    parser = argparse.ArgumentParser(description='Classify genomic regions by methylation state')
    parser.add_argument('--input', required=True, help='Input methylation file')
    parser.add_argument('--output', required=True, help='Output classified regions file')
    parser.add_argument('--format', default='bedgraph', 
                       choices=['bedgraph', 'bismark', 'custom'],
                       help='Input format type')
    parser.add_argument('--method', default='threshold',
                       choices=['threshold', 'adaptive', 'kmeans'],
                       help='Classification method')
    parser.add_argument('--hypo-threshold', type=float, default=0.3,
                       help='Hypomethylation threshold')
    parser.add_argument('--hyper-threshold', type=float, default=0.7,
                       help='Hypermethylation threshold')
    parser.add_argument('--cpg-islands', help='CpG islands BED file')
    parser.add_argument('--promoters', help='Promoters BED file')
    parser.add_argument('--enhancers', help='Enhancers BED file')
    parser.add_argument('--window-size', type=int, default=1000,
                       help='Window size for regional analysis')
    parser.add_argument('--plot-output', help='Output path for plots')
    parser.add_argument('--regional-analysis', action='store_true',
                       help='Perform regional methylation analysis')
    
    args = parser.parse_args()
    
    # Initialize classifier
    classifier = MethylationClassifier(
        hypo_threshold=args.hypo_threshold,
        hyper_threshold=args.hyper_threshold
    )
    
    try:
        # Load data
        df = classifier.load_methylation_data(args.input, args.format)
        
        # Apply classification method
        if args.method == 'threshold':
            df = classifier.classify_simple_threshold(df)
        elif args.method == 'adaptive':
            df = classifier.classify_adaptive_threshold(df)
        elif args.method == 'kmeans':
            df = classifier.classify_kmeans(df)
        
        # Add genomic context annotations
        annotation_files = {}
        if args.cpg_islands:
            annotation_files['cpg_islands'] = args.cpg_islands
        if args.promoters:
            annotation_files['promoters'] = args.promoters
        if args.enhancers:
            annotation_files['enhancers'] = args.enhancers
        
        if annotation_files:
            df = classifier.add_genomic_context(df, annotation_files)
        
        # Regional analysis if requested
        if args.regional_analysis:
            regional_df = classifier.calculate_regional_methylation(df, args.window_size)
            regional_output = args.output.replace('.bed', '_regional.bed').replace('.tsv', '_regional.tsv')
            regional_stats = classifier.generate_summary_stats(regional_df)
            classifier.save_results(regional_df, regional_output, regional_stats)
        
        # Generate summary statistics
        summary_stats = classifier.generate_summary_stats(df)
        
        # Create plots if requested
        if args.plot_output:
            classifier.plot_methylation_distribution(df, args.plot_output)
        
        # Save results
        classifier.save_results(df, args.output, summary_stats)
        
        print("Methylation classification completed successfully!")
        
    except Exception as e:
        logging.error(f"Error in methylation classification: {e}")
        raise

if __name__ == "__main__":
    main()

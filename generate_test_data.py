#!/usr/bin/env python3
"""
Test Data Generator for Methylation-Selection Analysis Pipeline
Generate synthetic datasets for testing pipeline functionality.
"""

import pandas as pd
import numpy as np
import argparse
import os
import random
from pathlib import Path
from datetime import datetime
import yaml

class TestDataGenerator:
    """
    Generate synthetic test data for methylation-selection analysis
    """
    
    def __init__(self, output_dir: str = "test_data", seed: int = 42):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.seed = seed
        np.random.seed(seed)
        random.seed(seed)
        
        # Define chromosomes and their approximate lengths (simplified)
        self.chromosomes = {
            'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559,
            'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
            'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717,
            'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
            'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
            'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
            'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983,
            'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
        }
        
        # Gene symbols for realistic annotation
        self.gene_symbols = [
            'TP53', 'KRAS', 'PIK3CA', 'APC', 'PTEN', 'BRCA1', 'BRCA2',
            'EGFR', 'MYC', 'RAS', 'BRAF', 'ATM', 'RB1', 'VHL', 'MLH1',
            'MSH2', 'MSH6', 'PMS2', 'CDKN2A', 'TERT', 'IDH1', 'IDH2',
            'NOTCH1', 'JAK2', 'KIT', 'PDGFRA', 'ALK', 'RET', 'MET',
            'FGFR1', 'FGFR2', 'FGFR3', 'ERBB2', 'ESR1', 'AR', 'CTNNB1'
        ] * 10  # Extend list for more genes
        
    def generate_methylation_data(self, n_sites: int = 100000) -> str:
        """Generate synthetic methylation data in BedGraph format"""
        print(f"Generating {n_sites:,} methylation sites...")
        
        methylation_data = []
        sites_per_chr = n_sites // len(self.chromosomes)
        
        for chrom, chrom_length in self.chromosomes.items():
            # Generate random genomic positions
            positions = sorted(np.random.randint(0, chrom_length-1000, sites_per_chr))
            
            for pos in positions:
                # Create realistic methylation patterns
                # - CpG islands: mostly hypomethylated
                # - Gene bodies: intermediate methylation
                # - Intergenic regions: hypermethylated
                
                region_type = np.random.choice(['cpg_island', 'gene_body', 'intergenic'], 
                                             p=[0.1, 0.4, 0.5])
                
                if region_type == 'cpg_island':
                    # CpG islands: low methylation with some variability
                    meth_level = np.clip(np.random.beta(1, 4), 0, 1)
                elif region_type == 'gene_body':
                    # Gene bodies: intermediate methylation
                    meth_level = np.clip(np.random.beta(2, 2), 0, 1)
                else:
                    # Intergenic: higher methylation
                    meth_level = np.clip(np.random.beta(4, 1), 0, 1)
                
                methylation_data.append({
                    'chr': chrom,
                    'start': pos,
                    'end': pos + 1,
                    'methylation_level': round(meth_level, 3)
                })
        
        # Save as BedGraph
        df = pd.DataFrame(methylation_data)
        output_file = self.output_dir / "synthetic_methylation.bedgraph"
        df.to_csv(output_file, sep='\t', index=False, header=False)
        
        print(f"Methylation data saved to: {output_file}")
        return str(output_file)
    
    def generate_mutations(self, mutation_type: str = 'germline', 
                          n_mutations: int = 10000) -> str:
        """Generate synthetic mutation data in VCF format"""
        print(f"Generating {n_mutations:,} {mutation_type} mutations...")
        
        vcf_header = [
            "##fileformat=VCFv4.2",
            f"##fileDate={datetime.now().strftime('%Y%m%d')}",
            "##source=synthetic_data_generator",
            "##reference=hg38",
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        ]
        
        mutations = []
        mutations_per_chr = n_mutations // len(self.chromosomes)
        
        # Define mutation patterns based on type
        if mutation_type == 'germline':
            # Germline mutations: more conservative, higher quality
            quality_dist = lambda: np.random.normal(100, 20)
            vaf_dist = lambda: np.random.choice([0.5, 1.0], p=[0.8, 0.2])  # Heterozygous/Homozygous
        else:
            # Somatic mutations: more variable quality and VAF
            quality_dist = lambda: np.random.exponential(50)
            vaf_dist = lambda: np.random.beta(2, 8)  # Lower VAF distribution
        
        for chrom, chrom_length in self.chromosomes.items():
            # Generate random positions
            positions = sorted(np.random.randint(1000, chrom_length-1000, mutations_per_chr))
            
            for i, pos in enumerate(positions):
                # Generate realistic mutation
                ref_base = np.random.choice(['A', 'T', 'G', 'C'])
                alt_base = np.random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref_base])
                
                # Quality score
                qual = max(1, int(quality_dist()))
                
                # Variant allele frequency and depth
                vaf = np.clip(vaf_dist(), 0.01, 0.99)
                depth = np.random.poisson(50)
                alt_depth = int(depth * vaf)
                ref_depth = depth - alt_depth
                
                # Genotype
                if vaf > 0.8:
                    genotype = "1/1"  # Homozygous alt
                elif vaf > 0.3:
                    genotype = "0/1"  # Heterozygous
                else:
                    genotype = "0/1"  # Low frequency het (somatic)
                
                # Create VCF record
                info = f"AC={alt_depth};AF={vaf:.3f}"
                format_field = "GT:AD:DP"
                sample_field = f"{genotype}:{ref_depth},{alt_depth}:{depth}"
                
                vcf_record = f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t{qual}\tPASS\t{info}\t{format_field}\t{sample_field}"
                mutations.append(vcf_record)
        
        # Write VCF file
        output_file = self.output_dir / f"synthetic_{mutation_type}_mutations.vcf"
        with open(output_file, 'w') as f:
            for header_line in vcf_header:
                f.write(header_line + '\n')
            for mutation in mutations:
                f.write(mutation + '\n')
        
        print(f"{mutation_type.title()} mutations saved to: {output_file}")
        return str(output_file)
    
    def generate_cpg_islands(self, n_islands: int = 1000) -> str:
        """Generate synthetic CpG islands annotation"""
        print(f"Generating {n_islands:,} CpG islands...")
        
        cpg_islands = []
        islands_per_chr = n_islands // len(self.chromosomes)
        
        for chrom, chrom_length in self.chromosomes.items():
            # Generate CpG island positions (typically 500-2000bp)
            for _ in range(islands_per_chr):
                start = np.random.randint(1000, chrom_length - 5000)
                length = np.random.randint(500, 2000)
                end = start + length
                
                cpg_islands.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'name': f"CpG_{len(cpg_islands)}",
                    'score': 1000,
                    'strand': '+'
                })
        
        # Save as BED file
        df = pd.DataFrame(cpg_islands)
        output_file = self.output_dir / "synthetic_cpg_islands.bed"
        df[['chr', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            output_file, sep='\t', index=False, header=False
        )
        
        print(f"CpG islands saved to: {output_file}")
        return str(output_file)
    
    def generate_promoters(self, n_promoters: int = 5000) -> str:
        """Generate synthetic promoter regions"""
        print(f"Generating {n_promoters:,} promoters...")
        
        promoters = []
        promoters_per_chr = n_promoters // len(self.chromosomes)
        
        for chrom, chrom_length in self.chromosomes.items():
            for i in range(promoters_per_chr):
                # Promoters are typically 2kb upstream of TSS
                tss = np.random.randint(2000, chrom_length - 2000)
                start = tss - 1000
                end = tss + 1000
                
                # Assign gene symbol
                gene_symbol = self.gene_symbols[len(promoters) % len(self.gene_symbols)]
                
                promoters.append({
                    'chr': chrom,
                    'start': start,
                    'end': end,
                    'name': f"{gene_symbol}_promoter",
                    'score': 1000,
                    'strand': np.random.choice(['+', '-'])
                })
        
        # Save as BED file
        df = pd.DataFrame(promoters)
        output_file = self.output_dir / "synthetic_promoters.bed"
        df[['chr', 'start', 'end', 'name', 'score', 'strand']].to_csv(
            output_file, sep='\t', index=False, header=False
        )
        
        print(f"Promoters saved to: {output_file}")
        return str(output_file)
    
    def generate_cadd_scores(self, mutations_file: str) -> str:
        """Generate synthetic CADD scores for mutations"""
        print("Generating CADD scores...")
        
        # Read mutations from VCF
        mutations = []
        with open(mutations_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    mutations.append({
                        'Chrom': fields[0],
                        'Pos': int(fields[1]),
                        'Ref': fields[3],
                        'Alt': fields[4]
                    })
        
        # Generate CADD scores (realistic distribution)
        cadd_data = []
        for mutation in mutations:
            # CADD scores follow a roughly exponential distribution
            # Higher scores = more deleterious
            raw_score = np.random.exponential(0.5)
            phred_score = max(0.1, raw_score * 10)
            
            cadd_data.append({
                'Chrom': mutation['Chrom'],
                'Pos': mutation['Pos'],
                'Ref': mutation['Ref'],
                'Alt': mutation['Alt'],
                'RawScore': round(raw_score, 3),
                'PHRED': round(phred_score, 1)
            })
        
        # Save CADD scores
        df = pd.DataFrame(cadd_data)
        output_file = self.output_dir / "synthetic_cadd_scores.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        
        print(f"CADD scores saved to: {output_file}")
        return str(output_file)
    
    def generate_gene_lengths(self, n_genes: int = 1000) -> str:
        """Generate synthetic gene length annotations"""
        print(f"Generating gene lengths for {n_genes} genes...")
        
        gene_lengths = []
        
        for i in range(n_genes):
            gene_symbol = self.gene_symbols[i % len(self.gene_symbols)]
            if i >= len(self.gene_symbols):
                gene_symbol += f"_{i // len(self.gene_symbols)}"
            
            # Realistic gene length distribution
            total_length = int(np.random.lognormal(8, 1.5))  # ~3kb median
            
            # Approximate synonymous/non-synonymous sites
            # Roughly 25% of positions are synonymous
            coding_length = int(total_length * 0.3)  # ~30% coding
            synonymous_sites = int(coding_length * 0.25)
            nonsynonymous_sites = coding_length - synonymous_sites
            
            gene_lengths.append({
                'gene_symbol': gene_symbol,
                'total_length': total_length,
                'coding_length': coding_length,
                'synonymous_sites': synonymous_sites,
                'nonsynonymous_sites': nonsynonymous_sites
            })
        
        # Save gene lengths
        df = pd.DataFrame(gene_lengths)
        output_file = self.output_dir / "synthetic_gene_lengths.tsv"
        df.to_csv(output_file, sep='\t', index=False)
        
        print(f"Gene lengths saved to: {output_file}")
        return str(output_file)
    
    def generate_test_config(self, generated_files: dict) -> str:
        """Generate test configuration file with generated data paths"""
        print("Creating test configuration file...")
        
        config = {
            'input_files': {
                'methylation_data': generated_files['methylation_data'],
                'germline_mutations': generated_files.get('germline_mutations'),
                'somatic_mutations': generated_files.get('somatic_mutations'),
                'cpg_islands': generated_files.get('cpg_islands'),
                'promoters': generated_files.get('promoters'),
                'cadd_scores': generated_files.get('cadd_scores'),
                'gene_lengths': generated_files.get('gene_lengths')
            },
            'parameters': {
                'hypo_threshold': 0.3,
                'hyper_threshold': 0.7,
                'classification_method': 'threshold',
                'methylation_format': 'bedgraph',
                'germline_format': 'vcf',
                'somatic_format': 'vcf',
                'regional_analysis': True,
                'window_size': 1000,
                'calculate_burdens': True,
                'calculate_selection': True,
                'min_mutations_per_gene': 3,  # Lower for test data
                'bootstrap_iterations': 100,  # Lower for speed
                'permutation_tests': 1000,    # Lower for speed
                'fdr_method': 'BH',
                'effect_size_threshold': 0.1,
                'significance_threshold': 0.05,
                'confidence_level': 0.95,
                'cores': 2,
                'continue_on_failure': False
            },
            'output_settings': {
                'cleanup_intermediates': False,
                'generate_html_report': True,
                'create_plots': True
            }
        }
        
        # Remove None values
        config['input_files'] = {k: v for k, v in config['input_files'].items() if v is not None}
        
        # Save configuration
        config_file = self.output_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=2)
        
        print(f"Test configuration saved to: {config_file}")
        return str(config_file)
    
    def generate_all_test_data(self, methylation_sites: int = 50000,
                              germline_mutations: int = 5000,
                              somatic_mutations: int = 8000) -> dict:
        """Generate complete test dataset"""
        print("="*60)
        print("GENERATING COMPLETE TEST DATASET")
        print("="*60)
        
        generated_files = {}
        
        # Generate core data files
        generated_files['methylation_data'] = self.generate_methylation_data(methylation_sites)
        generated_files['germline_mutations'] = self.generate_mutations('germline', germline_mutations)
        generated_files['somatic_mutations'] = self.generate_mutations('somatic', somatic_mutations)
        
        # Generate annotation files
        generated_files['cpg_islands'] = self.generate_cpg_islands(500)
        generated_files['promoters'] = self.generate_promoters(2000)
        generated_files['gene_lengths'] = self.generate_gene_lengths(500)
        
        # Generate functional annotation (CADD scores for germline mutations)
        generated_files['cadd_scores'] = self.generate_cadd_scores(
            generated_files['germline_mutations']
        )
        
        # Generate test configuration
        config_file = self.generate_test_config(generated_files)
        generated_files['config_file'] = config_file
        
        print("\n" + "="*60)
        print("TEST DATASET GENERATION COMPLETE")
        print("="*60)
        print(f"Output directory: {self.output_dir}")
        print(f"Configuration file: {config_file}")
        print("\nTo run the pipeline on test data:")
        print(f"python methylation_selection_pipeline.py --config {config_file} --output test_results/")
        print("\n" + "="*60)
        
        return generated_files

def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic test data for methylation-selection analysis pipeline'
    )
    parser.add_argument('--output', default='test_data', 
                       help='Output directory for test data')
    parser.add_argument('--methylation-sites', type=int, default=50000,
                       help='Number of methylation sites to generate')
    parser.add_argument('--germline-mutations', type=int, default=5000,
                       help='Number of germline mutations to generate')
    parser.add_argument('--somatic-mutations', type=int, default=8000,
                       help='Number of somatic mutations to generate')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducible data')
    parser.add_argument('--quick', action='store_true',
                       help='Generate smaller dataset for quick testing')
    
    args = parser.parse_args()
    
    # Adjust sizes for quick testing
    if args.quick:
        args.methylation_sites = 10000
        args.germline_mutations = 1000
        args.somatic_mutations = 1500
        print("Quick mode: generating smaller test dataset")
    
    # Initialize generator
    generator = TestDataGenerator(args.output, args.seed)
    
    # Generate complete test dataset
    generated_files = generator.generate_all_test_data(
        args.methylation_sites,
        args.germline_mutations, 
        args.somatic_mutations
    )
    
    # Summary
    total_size = sum(os.path.getsize(f) for f in generated_files.values() if os.path.exists(f))
    print(f"Total test data size: {total_size / 1024 / 1024:.1f} MB")
    
    return 0

if __name__ == "__main__":
    exit(main())

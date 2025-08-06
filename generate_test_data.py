#!/usr/bin/env python3
"""
Enhanced Test Data Generator for Methylation-Selection Analysis Pipeline
Generate realistic synthetic datasets with proper gene annotations.
"""

import pandas as pd
import numpy as np
import argparse
import os
import random
from pathlib import Path
from datetime import datetime
import yaml


class EnhancedTestDataGenerator:
    """
    Generate realistic synthetic test data with proper gene annotations
    """

    def __init__(self, output_dir: str = "test_data", seed: int = 42):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.seed = seed
        np.random.seed(seed)
        random.seed(seed)

        # Define chromosomes and their approximate lengths
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

        # Comprehensive list of real gene symbols
        self.gene_symbols = [
            # Tumor suppressors
            'TP53', 'RB1', 'APC', 'PTEN', 'VHL', 'NF1', 'NF2', 'BRCA1', 'BRCA2',
            'ATM', 'CHEK2', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'MUTYH', 'CDH1',
            'STK11', 'SMAD4', 'DPC4', 'CDKN2A', 'CDKN2B', 'CDKN1A', 'CDKN1B',

            # Oncogenes
            'KRAS', 'NRAS', 'HRAS', 'BRAF', 'PIK3CA', 'AKT1', 'AKT2', 'MYC',
            'MYCN', 'MYCL', 'EGFR', 'ERBB2', 'ERBB3', 'ERBB4', 'MET', 'ALK',
            'RET', 'ROS1', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'PDGFRA',
            'PDGFRB', 'KIT', 'FLT3', 'JAK1', 'JAK2', 'JAK3', 'STAT3', 'STAT5A',

            # DNA repair genes
            'BRIP1', 'PALB2', 'RAD51C', 'RAD51D', 'BARD1', 'NBN', 'MRE11A',
            'RAD50', 'BLM', 'WRN', 'RECQL4', 'FANCA', 'FANCC', 'FANCD2',
            'FANCF', 'FANCG', 'PARP1', 'PARP2', 'XRCC1', 'XRCC2', 'XRCC3',

            # Cell cycle genes
            'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CCNE2', 'CCNA1', 'CCNA2',
            'CCNB1', 'CCNB2', 'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CDC25A',
            'CDC25B', 'CDC25C', 'WEE1', 'CHEK1', 'CHEK2',

            # Apoptosis genes
            'BCL2', 'BCL2L1', 'BCL2L2', 'BCLAF1', 'BAX', 'BAK1', 'BAD',
            'BID', 'BIM', 'PUMA', 'NOXA', 'APAF1', 'CASP3', 'CASP7',
            'CASP8', 'CASP9', 'XIAP', 'SURVIVIN',

            # Transcription factors
            'MYB', 'ETS1', 'ETS2', 'FOS', 'JUN', 'SP1', 'SP3', 'E2F1',
            'E2F3', 'RB1', 'FOXO1', 'FOXO3', 'FOXP1', 'FOXP2', 'FOXP3',
            'RUNX1', 'RUNX2', 'RUNX3', 'GATA1', 'GATA2', 'GATA3',

            # Chromatin modifiers
            'DNMT1', 'DNMT3A', 'DNMT3B', 'TET1', 'TET2', 'TET3', 'IDH1',
            'IDH2', 'EZH2', 'SUZ12', 'EED', 'KMT2A', 'KMT2D', 'KDM5A',
            'KDM5B', 'KDM6A', 'CREBBP', 'EP300', 'HDAC1', 'HDAC2',

            # Metabolic genes
            'LDHA', 'LDHB', 'PKM', 'PFKL', 'PFKM', 'PFKP', 'G6PD',
            'TALDO1', 'TKT', 'GAPDH', 'ENO1', 'ENO2', 'PGK1', 'ALDOA',
            'FH', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'SUCLG1',

            # Immune genes
            'CD19', 'CD20', 'CD22', 'CD274', 'PDCD1', 'CTLA4', 'LAG3',
            'TIM3', 'TIGIT', 'ICOS', 'CD28', 'CD40', 'CD40L', 'TNF',
            'IFNG', 'IL2', 'IL10', 'IL12A', 'IL12B', 'IL6', 'IL1B',

            # Hormone signaling
            'ESR1', 'ESR2', 'AR', 'PGR', 'GR', 'TR', 'VDR', 'RXRA',
            'RXRB', 'PPARA', 'PPARG', 'RARA', 'RARB', 'RARG',

            # Development genes
            'WNT1', 'WNT2', 'WNT3A', 'WNT5A', 'WNT10B', 'CTNNB1', 'APC2',
            'GSK3B', 'AXIN1', 'AXIN2', 'TCF7L2', 'LEF1', 'DKK1', 'SFRP1',
            'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'JAG1', 'JAG2', 'DLL1',

            # Housekeeping genes
            'ACTB', 'GAPDH', 'HPRT1', 'TBP', 'RPLP0', 'RPL13A', 'RPS18',
            'YWHAZ', 'PGK1', 'LDHA', 'TFRC', 'GUSB', 'HMBS', 'UBC',

            # Additional cancer genes
            'SMARCA4', 'ARID1A', 'ARID1B', 'ARID2', 'PBRM1', 'BAP1',
            'SETD2', 'KDM5C', 'UTX', 'ASXL1', 'TET2', 'DNMT3A', 'NPM1',
            'FLT3', 'CEBPA', 'GATA2', 'RUNX1', 'ETV6', 'CBFB', 'MYH11',
            'PML', 'RARA', 'BCR', 'ABL1', 'EWSR1', 'FLI1', 'DDIT3',
            'CHOP', 'MDM2', 'MDM4', 'CDKN2A', 'CDKN2B'
        ]

        # Create gene database with realistic properties
        self.gene_database = self._create_gene_database()

    def _create_gene_database(self):
        """Create a database of genes with realistic genomic coordinates"""
        gene_db = []

        for i, gene_symbol in enumerate(self.gene_symbols):
            # Randomly assign to chromosome
            chrom = np.random.choice(list(self.chromosomes.keys()))
            chrom_length = self.chromosomes[chrom]

            # Generate realistic gene coordinates
            gene_length = int(np.random.lognormal(8.5, 1.5))  # ~5kb median, up to ~500kb
            gene_length = min(gene_length, 500000)  # Cap at 500kb

            start = np.random.randint(10000, chrom_length - gene_length - 10000)
            end = start + gene_length

            # Generate exon structure (simplified)
            n_exons = np.random.poisson(8) + 1  # Average ~8 exons
            coding_length = int(gene_length * np.random.uniform(0.05, 0.4))  # 5-40% coding

            # Synonymous/non-synonymous sites approximation
            synonymous_sites = int(coding_length * 0.25)
            nonsynonymous_sites = coding_length - synonymous_sites

            gene_db.append({
                'gene_symbol': gene_symbol,
                'chromosome': chrom,
                'start': start,
                'end': end,
                'strand': np.random.choice(['+', '-']),
                'gene_length': gene_length,
                'coding_length': coding_length,
                'n_exons': n_exons,
                'synonymous_sites': synonymous_sites,
                'nonsynonymous_sites': nonsynonymous_sites,
                'gene_type': self._assign_gene_type(gene_symbol)
            })

        return pd.DataFrame(gene_db)

    def _assign_gene_type(self, gene_symbol):
        """Assign gene type based on gene symbol"""
        tumor_suppressors = ['TP53', 'RB1', 'APC', 'PTEN', 'BRCA1', 'BRCA2', 'VHL', 'NF1', 'ATM']
        oncogenes = ['KRAS', 'NRAS', 'BRAF', 'PIK3CA', 'MYC', 'EGFR', 'ERBB2']
        dna_repair = ['MLH1', 'MSH2', 'MSH6', 'PMS2', 'BRIP1', 'PALB2', 'FANCA']
        housekeeping = ['ACTB', 'GAPDH', 'HPRT1', 'TBP', 'RPLP0']

        if gene_symbol in tumor_suppressors:
            return 'tumor_suppressor'
        elif gene_symbol in oncogenes:
            return 'oncogene'
        elif gene_symbol in dna_repair:
            return 'dna_repair'
        elif gene_symbol in housekeeping:
            return 'housekeeping'
        else:
            return 'other'

    def generate_annotated_mutations(self, mutation_type: str = 'germline',
                                     n_mutations: int = 10000) -> str:
        """Generate synthetic mutations with proper gene annotations"""
        print(f"Generating {n_mutations:,} {mutation_type} mutations with gene annotations...")

        mutations = []

        # Define mutation patterns based on selection pressures
        if mutation_type == 'germline':
            # Germline: stronger purifying selection
            functional_weights = {
                'tumor_suppressor': 0.3,  # Fewer deleterious mutations
                'oncogene': 0.7,  # Moderate constraint
                'dna_repair': 0.4,  # High constraint
                'housekeeping': 0.6,  # Moderate constraint
                'other': 1.0  # Neutral
            }
            quality_mean = 100
        else:
            # Somatic: more relaxed selection, cancer context
            functional_weights = {
                'tumor_suppressor': 1.5,  # More mutations in cancer
                'oncogene': 1.3,  # Activated in cancer
                'dna_repair': 1.2,  # Defective in cancer
                'housekeeping': 0.8,  # Still under selection
                'other': 1.0  # Neutral
            }
            quality_mean = 60

        # Generate mutations proportional to gene lengths and selection
        gene_weights = []
        for _, gene in self.gene_database.iterrows():
            base_weight = gene['gene_length'] / 1000  # Length-based
            func_weight = functional_weights.get(gene['gene_type'], 1.0)
            final_weight = base_weight * func_weight
            gene_weights.append(final_weight)

        gene_weights = np.array(gene_weights)
        gene_probs = gene_weights / gene_weights.sum()

        # Generate mutations
        selected_genes = np.random.choice(
            len(self.gene_database),
            size=n_mutations,
            p=gene_probs
        )

        for gene_idx in selected_genes:
            gene = self.gene_database.iloc[gene_idx]

            # Generate position within gene
            pos = np.random.randint(gene['start'], gene['end'])

            # Determine mutation type based on functional category
            if np.random.random() < 0.3:  # 30% coding
                if np.random.random() < 0.75:  # 75% missense
                    mut_consequence = 'missense_variant'
                    mut_impact = 'MODERATE'
                    mut_category = 'moderate_impact'
                elif np.random.random() < 0.9:  # Most remaining are synonymous
                    mut_consequence = 'synonymous_variant'
                    mut_impact = 'LOW'
                    mut_category = 'low_impact'
                else:  # Few are nonsense
                    mut_consequence = 'stop_gained'
                    mut_impact = 'HIGH'
                    mut_category = 'high_impact'
            else:  # 70% non-coding
                if np.random.random() < 0.4:  # Intronic
                    mut_consequence = 'intron_variant'
                    mut_impact = 'MODIFIER'
                    mut_category = 'intronic'
                else:  # Regulatory
                    mut_consequence = 'regulatory_region_variant'
                    mut_impact = 'MODIFIER'
                    mut_category = 'regulatory'

            # Generate alleles
            ref_base = np.random.choice(['A', 'T', 'G', 'C'])
            alt_base = np.random.choice([b for b in ['A', 'T', 'G', 'C'] if b != ref_base])

            # Quality and depth
            qual = max(10, int(np.random.normal(quality_mean, 20)))
            depth = np.random.poisson(40) + 10

            if mutation_type == 'germline':
                vaf = np.random.choice([0.5, 1.0], p=[0.8, 0.2])  # Het/Hom
                genotype = "0/1" if vaf == 0.5 else "1/1"
            else:
                vaf = np.random.beta(2, 8)  # Lower VAF for somatic
                vaf = max(0.05, vaf)  # Minimum detectable VAF
                genotype = "0/1"

            alt_depth = int(depth * vaf)
            ref_depth = depth - alt_depth

            mutations.append({
                'chr': gene['chromosome'],
                'pos': pos,
                'ref': ref_base,
                'alt': alt_base,
                'gene_symbol': gene['gene_symbol'],
                'consequence': mut_consequence,
                'impact': mut_impact,
                'mutation_type_category': mut_category,
                'gene_type': gene['gene_type'],
                'qual': qual,
                'depth': depth,
                'vaf': vaf,
                'genotype': genotype,
                'ref_depth': ref_depth,
                'alt_depth': alt_depth,
                'variant_id': f"{gene['chromosome']}:{pos}:{ref_base}:{alt_base}",
                'mutation_origin': mutation_type
            })

        # Save as TSV (easier for annotation pipeline)
        df = pd.DataFrame(mutations)
        output_file = self.output_dir / f"synthetic_{mutation_type}_mutations_annotated.tsv"
        df.to_csv(output_file, sep='\t', index=False)

        print(f"{mutation_type.title()} mutations saved to: {output_file}")

        # Also save as VCF for compatibility
        vcf_file = self._create_vcf_from_mutations(df, mutation_type)

        return str(output_file)

    def _create_vcf_from_mutations(self, mutations_df, mutation_type):
        """Create VCF file from annotated mutations DataFrame"""
        vcf_header = [
            "##fileformat=VCFv4.2",
            f"##fileDate={datetime.now().strftime('%Y%m%d')}",
            "##source=enhanced_synthetic_data_generator",
            "##reference=hg38",
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">',
            '##INFO=<ID=CSQ,Number=1,Type=String,Description="Consequence">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
            '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
        ]

        vcf_records = []
        for _, mut in mutations_df.iterrows():
            info = f"AC={mut['alt_depth']};AF={mut['vaf']:.3f};GENE={mut['gene_symbol']};CSQ={mut['consequence']}"
            format_field = "GT:AD:DP:VAF"
            sample_field = f"{mut['genotype']}:{mut['ref_depth']},{mut['alt_depth']}:{mut['depth']}:{mut['vaf']:.3f}"

            vcf_record = f"{mut['chr']}\t{mut['pos']}\t.\t{mut['ref']}\t{mut['alt']}\t{mut['qual']}\tPASS\t{info}\t{format_field}\t{sample_field}"
            vcf_records.append(vcf_record)

        # Write VCF file
        vcf_file = self.output_dir / f"synthetic_{mutation_type}_mutations.vcf"
        with open(vcf_file, 'w') as f:
            for header_line in vcf_header:
                f.write(header_line + '\n')
            for record in vcf_records:
                f.write(record + '\n')

        return str(vcf_file)

    def generate_methylation_data(self, n_sites: int = 100000) -> str:
        """Generate realistic methylation data with gene-aware patterns"""
        print(f"Generating {n_sites:,} methylation sites with gene-aware patterns...")

        methylation_data = []
        sites_per_chr = n_sites // len(self.chromosomes)

        for chrom, chrom_length in self.chromosomes.items():
            # Get genes on this chromosome
            chr_genes = self.gene_database[self.gene_database['chromosome'] == chrom]

            # Generate background methylation
            n_background = int(sites_per_chr * 0.7)  # 70% background
            background_positions = sorted(np.random.randint(1000, chrom_length - 1000, n_background))

            for pos in background_positions:
                # Background methylation varies by genomic context
                meth_level = np.random.beta(3, 2)  # Skewed toward higher methylation

                methylation_data.append({
                    'chr': chrom,
                    'start': pos,
                    'end': pos + 1,
                    'methylation_level': round(meth_level, 3)
                })

            # Generate gene-associated methylation
            n_gene_associated = sites_per_chr - n_background

            if len(chr_genes) > 0:
                for _ in range(n_gene_associated):
                    # Select random gene
                    gene = chr_genes.sample(1).iloc[0]

                    # Generate position relative to gene
                    region_type = np.random.choice(['promoter', 'gene_body', 'downstream'],
                                                   p=[0.3, 0.5, 0.2])

                    if region_type == 'promoter':
                        # Promoter region (-2kb to TSS)
                        pos = np.random.randint(max(1, gene['start'] - 2000), gene['start'])
                        # Promoters often hypomethylated
                        if gene['gene_type'] == 'housekeeping':
                            meth_level = np.random.beta(1, 5)  # Very low
                        else:
                            meth_level = np.random.beta(2, 3)  # Low-moderate

                    elif region_type == 'gene_body':
                        # Within gene body
                        pos = np.random.randint(gene['start'], gene['end'])
                        # Gene bodies often intermediately methylated
                        meth_level = np.random.beta(3, 3)  # Intermediate

                        # Tumor suppressors might have different patterns
                        if gene['gene_type'] == 'tumor_suppressor' and np.random.random() < 0.3:
                            meth_level = np.random.beta(5, 2)  # Higher (silencing)

                    else:
                        # Downstream region
                        pos = np.random.randint(gene['end'], min(chrom_length, gene['end'] + 5000))
                        meth_level = np.random.beta(4, 2)  # Higher methylation

                    if 1 <= pos < chrom_length:
                        methylation_data.append({
                            'chr': chrom,
                            'start': pos,
                            'end': pos + 1,
                            'methylation_level': round(np.clip(meth_level, 0, 1), 3)
                        })

        # Save methylation data
        df = pd.DataFrame(methylation_data)
        df = df.sort_values(['chr', 'start'])  # Sort for bedtools

        output_file = self.output_dir / "synthetic_methylation_with_genes.bedgraph"
        df.to_csv(output_file, sep='\t', index=False, header=False)

        print(f"Gene-aware methylation data saved to: {output_file}")
        return str(output_file)

    def generate_realistic_cadd_scores(self, mutations_file: str) -> str:
        """Generate realistic CADD scores based on mutation consequences"""
        print("Generating realistic CADD scores based on functional consequences...")

        # Read mutations
        mutations_df = pd.read_csv(mutations_file, sep='\t')

        cadd_data = []

        for _, mut in mutations_df.iterrows():
            # Base CADD score on consequence type
            consequence = mut.get('consequence', 'unknown')
            gene_type = mut.get('gene_type', 'other')

            if consequence == 'stop_gained':
                # Very high CADD scores for nonsense
                raw_score = np.random.gamma(2, 2) + 3
                phred_score = min(99, raw_score * 3 + np.random.normal(25, 5))
            elif consequence == 'missense_variant':
                # Moderate to high CADD scores for missense
                if gene_type in ['tumor_suppressor', 'dna_repair']:
                    raw_score = np.random.gamma(1.5, 1) + 1
                    phred_score = raw_score * 2 + np.random.normal(15, 8)
                else:
                    raw_score = np.random.gamma(1, 1) + 0.5
                    phred_score = raw_score * 2 + np.random.normal(10, 6)
            elif consequence == 'synonymous_variant':
                # Low CADD scores for synonymous
                raw_score = np.random.exponential(0.3)
                phred_score = raw_score + np.random.normal(2, 3)
            else:
                # Moderate scores for other variants
                raw_score = np.random.exponential(0.5)
                phred_score = raw_score * 1.5 + np.random.normal(5, 4)

            # Ensure realistic ranges
            phred_score = max(0.1, min(99, phred_score))
            raw_score = max(0.001, raw_score)

            cadd_data.append({
                'Chrom': mut['chr'],
                'Pos': mut['pos'],
                'Ref': mut['ref'],
                'Alt': mut['alt'],
                'RawScore': round(raw_score, 3),
                'PHRED': round(phred_score, 1)
            })

        # Save CADD scores
        df = pd.DataFrame(cadd_data)
        output_file = self.output_dir / "realistic_cadd_scores.tsv"
        df.to_csv(output_file, sep='\t', index=False)

        print(f"Realistic CADD scores saved to: {output_file}")
        return str(output_file)

    def save_gene_database(self) -> str:
        """Save the gene database for reference"""
        output_file = self.output_dir / "synthetic_gene_database.tsv"
        self.gene_database.to_csv(output_file, sep='\t', index=False)
        print(f"Gene database saved to: {output_file}")
        return str(output_file)

    def generate_comprehensive_test_data(self,
                                         methylation_sites: int = 50000,
                                         germline_mutations: int = 5000,
                                         somatic_mutations: int = 8000) -> dict:
        """Generate comprehensive test dataset with real gene annotations"""
        print("=" * 70)
        print("GENERATING COMPREHENSIVE TEST DATASET WITH REAL GENE ANNOTATIONS")
        print("=" * 70)

        generated_files = {}

        # Save gene database first
        generated_files['gene_database'] = self.save_gene_database()

        # Generate methylation data with gene-aware patterns
        generated_files['methylation_data'] = self.generate_methylation_data(methylation_sites)

        # Generate annotated mutations
        generated_files['germline_mutations'] = self.generate_annotated_mutations('germline', germline_mutations)
        generated_files['somatic_mutations'] = self.generate_annotated_mutations('somatic', somatic_mutations)

        # Generate realistic CADD scores
        generated_files['cadd_scores'] = self.generate_realistic_cadd_scores(
            generated_files['germline_mutations']
        )

        # Create analysis configuration
        config = self._create_analysis_config(generated_files)
        config_file = self.output_dir / "comprehensive_test_config.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, indent=2)

        generated_files['config_file'] = str(config_file)

        # Generate summary statistics
        self._generate_data_summary(generated_files)

        print("\n" + "=" * 70)
        print("COMPREHENSIVE TEST DATASET GENERATION COMPLETE")
        print("=" * 70)
        print(f"Output directory: {self.output_dir}")
        print(f"Configuration file: {config_file}")
        print(f"\nUnique genes in dataset: {len(self.gene_database)}")
        print(f"Gene types: {self.gene_database['gene_type'].value_counts().to_dict()}")

        print("\nTo run the pipeline:")
        print(f"# First, run methylation classification:")
        print(
            f"python methylation_classify.py --input {generated_files['methylation_data']} --output test_methylation_classified.tsv")
        print(f"\n# Then run mutation annotation:")
        print(
            f"python annotate_mutations.py --germline {generated_files['germline_mutations']} --somatic {generated_files['somatic_mutations']} --methylation-regions test_methylation_classified.tsv --output test_annotated_mutations.tsv")
        print(f"\n# Finally run selection analysis:")
        print(f"Rscript calculate_selection.r --mutations test_annotated_mutations.tsv --output test_selection_results")

        return generated_files

    def _create_analysis_config(self, generated_files):
        """Create analysis configuration"""
        return {
            'input_files': generated_files,
            'analysis_parameters': {
                'methylation_thresholds': {
                    'hypomethylated': 0.3,
                    'hypermethylated': 0.7
                },
                'selection_analysis': {
                    'min_mutations_per_gene': 3,
                    'bootstrap_iterations': 500,
                    'significance_threshold': 0.05
                },
                'quality_filters': {
                    'min_depth': 10,
                    'min_quality': 20,
                    'min_vaf': 0.05
                }
            },
            'gene_categories': {
                'high_constraint': ['tumor_suppressor', 'dna_repair', 'housekeeping'],
                'low_constraint': ['oncogene', 'other'],
                'methylation_sensitive': ['tumor_suppressor']
            }
        }

    def _generate_data_summary(self, generated_files):
        """Generate summary statistics for the test data"""
        summary_file = self.output_dir / "test_data_summary.txt"

        with open(summary_file, 'w') as f:
            f.write("TEST DATA SUMMARY\n")
            f.write("================\n\n")
            f.write(f"Generation date: {datetime.now()}\n")
            f.write(f"Random seed: {self.seed}\n\n")

            f.write("GENE DATABASE:\n")
            f.write(f"Total genes: {len(self.gene_database)}\n")
            for gene_type, count in self.gene_database['gene_type'].value_counts().items():
                f.write(f"  {gene_type}: {count}\n")

            f.write(f"\nAverage gene length: {self.gene_database['gene_length'].mean():.0f} bp\n")
            f.write(f"Average coding length: {self.gene_database['coding_length'].mean():.0f} bp\n\n")

            # Mutation summaries
            for mut_type in ['germline', 'somatic']:
                mut_file = generated_files.get(f'{mut_type}_mutations')
                if mut_file and os.path.exists(mut_file):
                    df = pd.read_csv(mut_file, sep='\t')
                    f.write(f"{mut_type.upper()} MUTATIONS:\n")
                    f.write(f"Total: {len(df)}\n")
                    f.write("By consequence:\n")
                    for cons, count in df['consequence'].value_counts().items():
                        f.write(f"  {cons}: {count}\n")
                    f.write("By gene type:\n")
                    for gtype, count in df['gene_type'].value_counts().items():
                        f.write(f"  {gtype}: {count}\n")
                    f.write(f"Average VAF: {df['vaf'].mean():.3f}\n\n")

        print(f"Data summary saved to: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate comprehensive test data with real gene annotations'
    )
    parser.add_argument('--output', default='realistic_test_data',
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

    # Initialize enhanced generator
    generator = EnhancedTestDataGenerator(args.output, args.seed)

    # Generate comprehensive test dataset
    generated_files = generator.generate_comprehensive_test_data(
        args.methylation_sites,
        args.germline_mutations,
        args.somatic_mutations
    )

    # Summary
    total_size = sum(
        os.path.getsize(f) for f in generated_files.values()
        if isinstance(f, str) and os.path.exists(f)
    )
    print(f"\nTotal test data size: {total_size / 1024 / 1024:.1f} MB")
    print(f"Files generated: {len(generated_files)}")

    return 0


if __name__ == "__main__":
    exit(main())
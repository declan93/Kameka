#!/usr/bin/env python3
"""
Mutation Annotation Tool
Annotates germline and somatic mutations with methylation states and functional consequences.
"""

import pandas as pd
import numpy as np
import argparse
import logging
from typing import Dict, List, Tuple, Optional, Set
import pysam
import pybedtools
from collections import defaultdict
import vcf
import warnings

warnings.filterwarnings('ignore')


class MutationAnnotator:
    """
    Comprehensive mutation annotator with methylation state integration
    """

    def __init__(self):
        self.logger = self._setup_logging()
        self.mutation_types = {
            'synonymous': 0,
            'missense': 1,
            'nonsense': 2,
            'splice_site': 3,
            'frameshift': 4,
            'inframe_indel': 5,
            'upstream': 6,
            'downstream': 7,
            'intronic': 8,
            'intergenic': 9
        }

    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('annotate_mutations.log'),
                logging.StreamHandler()
            ]
        )
        return logging.getLogger(__name__)

    def load_vcf_mutations(self, vcf_path: str, mutation_origin: str) -> pd.DataFrame:
        """
        Load mutations from VCF file

        Args:
            vcf_path: Path to VCF file
            mutation_origin: 'germline' or 'somatic'

        Returns:
            DataFrame with mutation information
        """
        self.logger.info(f"Loading {mutation_origin} mutations from {vcf_path}")

        mutations = []
        vcf_reader = vcf.Reader(open(vcf_path, 'r'))

        for record in vcf_reader:
            # Basic mutation info
            mutation_info = {
                'chr': str(record.CHROM),
                'pos': record.POS,
                'start': record.POS - 1,  # Convert to 0-based for BED compatibility
                'end': record.POS,
                'ref': record.REF,
                'alt': str(record.ALT[0]) if record.ALT else '.',
                'mutation_origin': mutation_origin,
                'variant_id': f"{record.CHROM}:{record.POS}:{record.REF}:{record.ALT[0] if record.ALT else '.'}",
                'qual': record.QUAL if record.QUAL else 0,
                'filter_status': ';'.join(record.FILTER) if record.FILTER else 'PASS'
            }

            # Add INFO fields
            for info_key, info_value in record.INFO.items():
                mutation_info[f'info_{info_key}'] = info_value

            # Process samples if available
            if record.samples:
                for i, sample in enumerate(record.samples):
                    sample_prefix = f'sample_{i}_' if len(record.samples) > 1 else ''

                    # Genotype information
                    if hasattr(sample.data, 'GT') and sample.data.GT:
                        mutation_info[f'{sample_prefix}genotype'] = sample.data.GT

                        # Calculate allele frequency for somatic mutations
                        if mutation_origin == 'somatic':
                            if hasattr(sample.data, 'AD') and sample.data.AD:
                                ref_count, alt_count = sample.data.AD[:2]
                                total_depth = ref_count + alt_count
                                if total_depth > 0:
                                    mutation_info[f'{sample_prefix}vaf'] = alt_count / total_depth
                                    mutation_info[f'{sample_prefix}depth'] = total_depth

                    # Add other sample-level annotations
                    for field in ['DP', 'GQ', 'PL']:
                        if hasattr(sample.data, field):
                            value = getattr(sample.data, field)
                            if value is not None:
                                mutation_info[f'{sample_prefix}{field.lower()}'] = value

            mutations.append(mutation_info)

        df = pd.DataFrame(mutations)
        self.logger.info(f"Loaded {len(df)} {mutation_origin} mutations")

        return df

    def load_mutation_table(self, table_path: str, mutation_origin: str) -> pd.DataFrame:
        """
        Load mutations from tab-separated table

        Expected columns: chr, pos, ref, alt, [additional annotations]
        """
        self.logger.info(f"Loading {mutation_origin} mutations from table {table_path}")

        df = pd.read_csv(table_path, sep='\t', low_memory=False)

        # Standardize column names
        column_mapping = {
            'chromosome': 'chr', 'chrom': 'chr',
            'position': 'pos', 'start_pos': 'pos',
            'reference': 'ref', 'ref_allele': 'ref',
            'alternate': 'alt', 'alt_allele': 'alt'
        }

        for old_col, new_col in column_mapping.items():
            if old_col in df.columns and new_col not in df.columns:
                df = df.rename(columns={old_col: new_col})

        # Ensure required columns exist
        required_cols = ['chr', 'pos', 'ref', 'alt']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Add standard columns
        df['start'] = df['pos'] - 1  # 0-based start
        df['end'] = df['pos']
        df['mutation_origin'] = mutation_origin
        df['variant_id'] = df.apply(
            lambda row: f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}", axis=1
        )

        self.logger.info(f"Loaded {len(df)} {mutation_origin} mutations")
        return df

    def annotate_with_methylation(self, mutations_df: pd.DataFrame,
                                  methylation_regions: str) -> pd.DataFrame:
        """
        Annotate mutations with methylation states from classified regions
        """
        self.logger.info("Annotating mutations with methylation states")

        # Load methylation regions
        try:
            meth_df = pd.read_csv(methylation_regions, sep='\t')
        except Exception as e:
            self.logger.error(f"Failed to load methylation regions file: {e}")
            raise

        self.logger.info(f"Loaded methylation file with columns: {list(meth_df.columns)}")
        self.logger.info(f"Methylation file shape: {meth_df.shape}")

        # Check if file has no header (first row used as column names)
        # This happens when column names look like data (chromosome names, positions, etc.)
        def looks_like_data_not_header(columns):
            """Check if column names look like data rather than headers"""
            if len(columns) < 4:
                return False

            first_col = str(columns[0])
            second_col = str(columns[1])

            # Check if first column looks like chromosome (chr1, chr2, 1, 2, etc.)
            chr_like = (first_col.startswith('chr') or first_col.isdigit() or first_col in ['X', 'Y', 'M', 'MT'])

            # Check if second column looks like genomic position (numeric)
            pos_like = second_col.replace('.', '').isdigit()

            return chr_like and pos_like

        # Handle different formats and column naming
        if looks_like_data_not_header(meth_df.columns):
            self.logger.info("Detected file without header - first row used as column names")
            # Re-read without treating first row as header
            try:
                meth_df = pd.read_csv(methylation_regions, sep='\t', header=None)
                if len(meth_df.columns) == 4:
                    meth_df.columns = ['chr', 'start', 'end', 'methylation_level']
                    self.logger.info("Assigned 4-column BED format column names")
                elif len(meth_df.columns) == 5:
                    meth_df.columns = ['chr', 'start', 'end', 'methylation_level', 'methylation_class']
                    self.logger.info("Assigned 5-column BED format column names")
                else:
                    raise ValueError(f"Unexpected number of columns: {len(meth_df.columns)}")
            except Exception as e:
                self.logger.error(f"Failed to re-read methylation file without header: {e}")
                raise

        elif len(meth_df.columns) == 4 and meth_df.columns[0] == 0:  # BED format without header (numeric columns)
            meth_df.columns = ['chr', 'start', 'end', 'methylation_level']
            self.logger.info("Detected 4-column BED format, assigned standard column names")
        elif len(meth_df.columns) >= 4 and all(col in meth_df.columns for col in ['chr', 'start', 'end']):
            # Has proper column names, check for methylation_level column
            if 'methylation_level' not in meth_df.columns:
                # Try to find methylation level column by position or alternative names
                possible_meth_cols = [col for col in meth_df.columns if
                                      any(term in col.lower() for term in ['methylation', 'level', 'beta', 'value'])]
                if possible_meth_cols:
                    meth_df['methylation_level'] = meth_df[possible_meth_cols[0]]
                    self.logger.info(f"Using column '{possible_meth_cols[0]}' as methylation_level")
                elif len(meth_df.columns) >= 4:
                    # Assume 4th column is methylation level
                    fourth_col = meth_df.columns[3]
                    meth_df['methylation_level'] = meth_df[fourth_col]
                    self.logger.info(f"Using column '{fourth_col}' as methylation_level")
                else:
                    raise ValueError("Could not identify methylation level column")
        else:
            raise ValueError(f"Unexpected methylation file format. Columns: {list(meth_df.columns)}")

        # Ensure chromosome column is string type
        meth_df['chr'] = meth_df['chr'].astype(str)

        # Ensure start and end are numeric
        meth_df['start'] = pd.to_numeric(meth_df['start'], errors='coerce')
        meth_df['end'] = pd.to_numeric(meth_df['end'], errors='coerce')

        # Ensure methylation_level is numeric
        meth_df['methylation_level'] = pd.to_numeric(meth_df['methylation_level'], errors='coerce')

        # Remove rows with missing essential data
        before_cleanup = len(meth_df)
        meth_df = meth_df.dropna(subset=['chr', 'start', 'end', 'methylation_level'])
        after_cleanup = len(meth_df)
        if before_cleanup != after_cleanup:
            self.logger.info(f"Removed {before_cleanup - after_cleanup} rows with missing data")

        # Classify methylation if not already done
        if 'methylation_class' not in meth_df.columns:
            self.logger.info("Classifying methylation states on-the-fly")
            conditions = [
                meth_df['methylation_level'] < 0.3,
                (meth_df['methylation_level'] >= 0.3) & (meth_df['methylation_level'] <= 0.7),
                meth_df['methylation_level'] > 0.7
            ]
            choices = ['hypomethylated', 'intermediate', 'hypermethylated']
            meth_df['methylation_class'] = np.select(conditions, choices)

        self.logger.info(f"Final methylation file shape: {meth_df.shape}")
        self.logger.info(f"Methylation classes: {meth_df['methylation_class'].value_counts().to_dict()}")

        def sort_chromosomes(df):
            """Sort dataframe by chromosomes in proper genomic order"""
            # Create a copy to avoid modifying the original
            df_copy = df.copy()

            # Create a sorting key for chromosomes
            def chr_sort_key(chr_name):
                chr_str = str(chr_name).upper()
                if chr_str.startswith('CHR'):
                    chr_str = chr_str[3:]  # Remove 'CHR' prefix

                # Handle numeric chromosomes
                if chr_str.isdigit():
                    return (0, int(chr_str))
                elif chr_str == 'X':
                    return (1, 0)
                elif chr_str == 'Y':
                    return (1, 1)
                elif chr_str in ['M', 'MT']:
                    return (1, 2)
                else:
                    # Other chromosomes sorted alphabetically
                    return (2, chr_str)

            df_copy['chr_sort_key'] = df_copy['chr'].apply(chr_sort_key)
            df_sorted = df_copy.sort_values(['chr_sort_key', 'start', 'end']).drop('chr_sort_key', axis=1)
            return df_sorted

        # Sort dataframes for bedtools compatibility
        mutations_sorted = mutations_df[['chr', 'start', 'end']].assign(name=mutations_df.index)
        mutations_sorted = sort_chromosomes(mutations_sorted)

        meth_sorted = meth_df[['chr', 'start', 'end', 'methylation_level', 'methylation_class']]
        meth_sorted = sort_chromosomes(meth_sorted)

        # Create BedTool objects for intersection
        mutations_bed = pybedtools.BedTool.from_dataframe(mutations_sorted)
        methylation_bed = pybedtools.BedTool.from_dataframe(meth_sorted)

        # Perform intersection with error handling
        try:
            intersected = mutations_bed.intersect(methylation_bed, wa=True, wb=True)
        except Exception as e:
            self.logger.warning(f"Bedtools intersect with sorted input failed: {e}")
            self.logger.info("Retrying with sorted=False")
            intersected = mutations_bed.intersect(methylation_bed, wa=True, wb=True, sorted=False)

        # Initialize methylation columns
        mutations_df['methylation_level'] = np.nan
        mutations_df['methylation_class'] = 'unknown'
        mutations_df['distance_to_meth_site'] = np.inf

        # Process intersections
        if len(intersected) > 0:
            try:
                intersect_df = intersected.to_dataframe()

                # Map intersected data back to mutations
                for _, row in intersect_df.iterrows():
                    mutation_idx = int(row['name'])  # Get original mutation index
                    if mutation_idx < len(mutations_df):
                        # The methylation data starts at column index 4 (after chr, start, end, name)
                        meth_level_col_idx = 7  # Position of methylation_level in intersected output
                        meth_class_col_idx = 8  # Position of methylation_class in intersected output

                        if len(row) > meth_level_col_idx:
                            mutations_df.loc[mutation_idx, 'methylation_level'] = row.iloc[meth_level_col_idx]
                        if len(row) > meth_class_col_idx:
                            mutations_df.loc[mutation_idx, 'methylation_class'] = row.iloc[meth_class_col_idx]
                        mutations_df.loc[mutation_idx, 'distance_to_meth_site'] = 0

            except Exception as e:
                self.logger.warning(f"Error processing intersections: {e}")
                # Fall back to simple overlap checking
                for idx, mut_row in mutations_df.iterrows():
                    overlapping_sites = meth_df[
                        (meth_df['chr'] == mut_row['chr']) &
                        (meth_df['start'] <= mut_row['end']) &
                        (meth_df['end'] >= mut_row['start'])
                        ]

                    if len(overlapping_sites) > 0:
                        closest_site = overlapping_sites.iloc[0]  # Take first overlapping site
                        mutations_df.loc[idx, 'methylation_level'] = closest_site['methylation_level']
                        mutations_df.loc[idx, 'methylation_class'] = closest_site['methylation_class']
                        mutations_df.loc[idx, 'distance_to_meth_site'] = 0

        # For mutations without direct overlap, find nearest methylation site
        unassigned_mask = mutations_df['methylation_class'] == 'unknown'
        if unassigned_mask.sum() > 0:
            self.logger.info(f"Finding nearest methylation sites for {unassigned_mask.sum()} mutations")

            # Create sorted BED file for unassigned mutations
            unassigned_data = mutations_df[unassigned_mask][['chr', 'start', 'end']].assign(
                name=mutations_df[unassigned_mask].index
            )
            unassigned_data = sort_chromosomes(unassigned_data)

            unassigned_bed = pybedtools.BedTool.from_dataframe(unassigned_data)

            # Find closest methylation sites (use sorted=False as backup)
            try:
                closest = unassigned_bed.closest(methylation_bed, d=True)
            except Exception as e:
                self.logger.warning(f"Bedtools closest with sorted input failed: {e}")
                self.logger.info("Retrying with sorted=False")
                closest = unassigned_bed.closest(methylation_bed, d=True, sorted=False)

            if len(closest) > 0:
                closest_df = closest.to_dataframe()
                for _, row in closest_df.iterrows():
                    mutation_idx = row['name']
                    if mutation_idx in mutations_df.index:
                        distance = row.iloc[-1]  # Distance is the last column
                        if distance <= 10000:  # Within 10kb
                            # Get methylation info from the correct columns
                            meth_level_idx = len(row) - 3  # methylation_level should be 3rd from end
                            meth_class_idx = len(row) - 2  # methylation_class should be 2nd from end

                            mutations_df.loc[mutation_idx, 'methylation_level'] = row.iloc[meth_level_idx]
                            mutations_df.loc[mutation_idx, 'methylation_class'] = row.iloc[meth_class_idx]
                            mutations_df.loc[mutation_idx, 'distance_to_meth_site'] = distance

        # Classify remaining unknown regions as intermediate
        unknown_mask = mutations_df['methylation_class'] == 'unknown'
        mutations_df.loc[unknown_mask, 'methylation_class'] = 'intermediate'

        self.logger.info(f"Methylation annotation complete. Distribution:")
        class_counts = mutations_df['methylation_class'].value_counts()
        for class_name, count in class_counts.items():
            self.logger.info(f"  {class_name}: {count} ({count / len(mutations_df) * 100:.1f}%)")

        return mutations_df

    def annotate_functional_consequences(self, mutations_df: pd.DataFrame,
                                         annotation_sources: Dict[str, str]) -> pd.DataFrame:
        """
        Annotate mutations with functional consequences from multiple sources

        Args:
            mutations_df: Mutations DataFrame
            annotation_sources: Dict with annotation type and file paths
                               e.g., {'vep': 'vep_annotations.txt',
                                     'annovar': 'annovar_annotations.txt',
                                     'cadd': 'cadd_scores.txt'}
        """
        self.logger.info("Annotating functional consequences")

        # Initialize functional annotation columns
        mutations_df['gene_symbol'] = 'unknown'
        mutations_df['consequence'] = 'unknown'
        mutations_df['impact'] = 'unknown'
        mutations_df['cadd_score'] = np.nan
        mutations_df['sift_score'] = np.nan
        mutations_df['polyphen_score'] = np.nan
        mutations_df['gerp_score'] = np.nan
        mutations_df['phylop_score'] = np.nan
        mutations_df['is_coding'] = False
        mutations_df['is_splice_site'] = False

        # Process VEP annotations
        if 'vep' in annotation_sources:
            mutations_df = self._annotate_vep(mutations_df, annotation_sources['vep'])

        # Process ANNOVAR annotations
        if 'annovar' in annotation_sources:
            mutations_df = self._annotate_annovar(mutations_df, annotation_sources['annovar'])

        # Process CADD scores
        if 'cadd' in annotation_sources:
            mutations_df = self._annotate_cadd(mutations_df, annotation_sources['cadd'])

        # Process conservation scores
        if 'conservation' in annotation_sources:
            mutations_df = self._annotate_conservation(mutations_df, annotation_sources['conservation'])

        # Assign mutation type categories
        mutations_df['mutation_type_category'] = mutations_df.apply(
            self._categorize_mutation, axis=1
        )

        return mutations_df

    def _annotate_vep(self, mutations_df: pd.DataFrame, vep_file: str) -> pd.DataFrame:
        """Annotate with VEP consequences"""
        try:
            vep_df = pd.read_csv(vep_file, sep='\t', comment='#', low_memory=False)

            # Create variant ID for matching
            if 'Location' in vep_df.columns:
                # Parse VEP location format
                vep_df['variant_id'] = vep_df.apply(
                    lambda row: self._parse_vep_location(row), axis=1
                )
            else:
                # Construct from individual columns if available
                required_cols = ['#Uploaded_variation', 'Location', 'Allele']
                if all(col in vep_df.columns for col in required_cols):
                    vep_df['variant_id'] = vep_df['#Uploaded_variation']

            # Merge with mutations
            vep_cols = ['variant_id', 'SYMBOL', 'Consequence', 'IMPACT',
                        'SIFT', 'PolyPhen', 'BIOTYPE']
            available_vep_cols = [col for col in vep_cols if col in vep_df.columns]

            merged = mutations_df.merge(
                vep_df[available_vep_cols],
                on='variant_id',
                how='left'
            )

            # Update annotations
            if 'SYMBOL' in merged.columns:
                mutations_df['gene_symbol'] = merged['SYMBOL'].fillna(mutations_df['gene_symbol'])
            if 'Consequence' in merged.columns:
                mutations_df['consequence'] = merged['Consequence'].fillna(mutations_df['consequence'])
            if 'IMPACT' in merged.columns:
                mutations_df['impact'] = merged['IMPACT'].fillna(mutations_df['impact'])

            self.logger.info("VEP annotations added")

        except Exception as e:
            self.logger.warning(f"Could not process VEP annotations: {e}")

        return mutations_df

    def _annotate_annovar(self, mutations_df: pd.DataFrame, annovar_file: str) -> pd.DataFrame:
        """Annotate with ANNOVAR consequences"""
        try:
            annovar_df = pd.read_csv(annovar_file, sep='\t', low_memory=False)

            # ANNOVAR typically has columns: Chr, Start, End, Ref, Alt, Func.refGene, Gene.refGene, etc.
            if all(col in annovar_df.columns for col in ['Chr', 'Start', 'Ref', 'Alt']):
                annovar_df['variant_id'] = annovar_df.apply(
                    lambda row: f"{row['Chr']}:{row['Start']}:{row['Ref']}:{row['Alt']}", axis=1
                )

                # Map ANNOVAR columns to standard names
                annovar_mapping = {
                    'Gene.refGene': 'gene_symbol',
                    'Func.refGene': 'consequence',
                    'ExonicFunc.refGene': 'exonic_function',
                    'CADD_raw': 'cadd_score',
                    'SIFT_score': 'sift_score',
                    'Polyphen2_HDIV_score': 'polyphen_score',
                    'GERP++_RS': 'gerp_score',
                    'phyloP30way_mammalian': 'phylop_score'
                }

                for annovar_col, std_col in annovar_mapping.items():
                    if annovar_col in annovar_df.columns:
                        merged = mutations_df.merge(
                            annovar_df[['variant_id', annovar_col]],
                            on='variant_id',
                            how='left'
                        )
                        mutations_df[std_col] = merged[annovar_col].fillna(mutations_df[std_col])

                self.logger.info("ANNOVAR annotations added")

        except Exception as e:
            self.logger.warning(f"Could not process ANNOVAR annotations: {e}")

        return mutations_df

    def _annotate_cadd(self, mutations_df: pd.DataFrame, cadd_file: str) -> pd.DataFrame:
        """Annotate with CADD scores"""
        try:
            cadd_df = pd.read_csv(cadd_file, sep='\t', comment='#', low_memory=False)

            # CADD format: Chrom, Pos, Ref, Alt, RawScore, PHRED
            if all(col in cadd_df.columns for col in ['Chrom', 'Pos', 'Ref', 'Alt']):
                cadd_df['variant_id'] = cadd_df.apply(
                    lambda row: f"{row['Chrom']}:{row['Pos']}:{row['Ref']}:{row['Alt']}", axis=1
                )

                score_col = 'PHRED' if 'PHRED' in cadd_df.columns else 'RawScore'

                merged = mutations_df.merge(
                    cadd_df[['variant_id', score_col]],
                    on='variant_id',
                    how='left'
                )

                mutations_df['cadd_score'] = pd.to_numeric(
                    merged[score_col], errors='coerce'
                ).fillna(mutations_df['cadd_score'])

                self.logger.info(f"CADD scores added (median: {mutations_df['cadd_score'].median():.2f})")

        except Exception as e:
            self.logger.warning(f"Could not process CADD scores: {e}")

        return mutations_df

    def _annotate_conservation(self, mutations_df: pd.DataFrame, conservation_file: str) -> pd.DataFrame:
        """Annotate with evolutionary conservation scores"""
        try:
            cons_df = pd.read_csv(conservation_file, sep='\t', low_memory=False)

            # Assume format: chr, pos, gerp_score, phylop_score
            if all(col in cons_df.columns for col in ['chr', 'pos']):
                cons_df['variant_key'] = cons_df['chr'].astype(str) + ':' + cons_df['pos'].astype(str)
                mutations_df['variant_key'] = mutations_df['chr'].astype(str) + ':' + mutations_df['pos'].astype(str)

                merged = mutations_df.merge(cons_df, on='variant_key', how='left', suffixes=('', '_cons'))

                if 'gerp_score' in merged.columns:
                    mutations_df['gerp_score'] = merged['gerp_score'].fillna(mutations_df['gerp_score'])
                if 'phylop_score' in merged.columns:
                    mutations_df['phylop_score'] = merged['phylop_score'].fillna(mutations_df['phylop_score'])

                self.logger.info("Conservation scores added")

        except Exception as e:
            self.logger.warning(f"Could not process conservation scores: {e}")

        return mutations_df

    def _categorize_mutation(self, row: pd.Series) -> str:
        """Categorize mutation based on consequence annotations"""
        consequence = str(row.get('consequence', '')).lower()
        impact = str(row.get('impact', '')).lower()

        # High impact mutations
        if any(term in consequence for term in ['stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor']):
            return 'high_impact'

        # Moderate impact mutations
        elif any(term in consequence for term in ['missense', 'inframe_deletion', 'inframe_insertion']):
            return 'moderate_impact'

        # Low impact mutations
        elif any(term in consequence for term in ['synonymous', 'stop_retained']):
            return 'low_impact'

        # Regulatory mutations
        elif any(term in consequence for term in ['regulatory', 'upstream', 'downstream', '5_prime', '3_prime']):
            return 'regulatory'

        # Intronic mutations
        elif 'intron' in consequence:
            return 'intronic'

        else:
            return 'other'

    def _parse_vep_location(self, row: pd.Series) -> str:
        """Parse VEP location format to create variant ID"""
        try:
            location = row['Location']
            allele = row['Allele']

            # Parse location (format: chr:pos or chr:start-end)
            if ':' in location:
                chrom, pos_info = location.split(':', 1)
                if '-' in pos_info:
                    start, end = pos_info.split('-')
                    pos = start
                else:
                    pos = pos_info

                # For VEP, we need to determine ref/alt
                # This is simplified - in practice, you'd need the original reference
                uploaded_var = row.get('#Uploaded_variation', '')
                if '_' in uploaded_var:
                    parts = uploaded_var.split('_')
                    if len(parts) >= 4:
                        return f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}"

                return f"{chrom}:{pos}:N:{allele}"  # Placeholder

        except Exception:
            pass

        return row.get('#Uploaded_variation', 'unknown')

    def calculate_mutation_burdens(self, mutations_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate various mutation burden metrics
        """
        self.logger.info("Calculating mutation burdens")

        # Overall burden metrics
        burden_stats = {}

        # By mutation origin
        for origin in mutations_df['mutation_origin'].unique():
            origin_muts = mutations_df[mutations_df['mutation_origin'] == origin]
            burden_stats[f'{origin}_total'] = len(origin_muts)

            # By methylation class
            for meth_class in origin_muts['methylation_class'].unique():
                class_muts = origin_muts[origin_muts['methylation_class'] == meth_class]
                burden_stats[f'{origin}_{meth_class}'] = len(class_muts)

            # By functional category
            for func_cat in origin_muts['mutation_type_category'].unique():
                func_muts = origin_muts[origin_muts['mutation_type_category'] == func_cat]
                burden_stats[f'{origin}_{func_cat}'] = len(func_muts)

        # Calculate ratios
        total_germline = burden_stats.get('germline_total', 0)
        total_somatic = burden_stats.get('somatic_total', 0)

        if total_germline > 0 and total_somatic > 0:
            burden_stats['somatic_germline_ratio'] = total_somatic / total_germline

        # Methylation-specific ratios
        for meth_class in ['hypomethylated', 'hypermethylated', 'intermediate']:
            germline_class = burden_stats.get(f'germline_{meth_class}', 0)
            somatic_class = burden_stats.get(f'somatic_{meth_class}', 0)

            if germline_class > 0:
                burden_stats[f'{meth_class}_somatic_germline_ratio'] = somatic_class / germline_class

        # Add burden stats as metadata
        mutations_df.attrs['burden_stats'] = burden_stats

        return mutations_df

    def calculate_selection_metrics(self, mutations_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate preliminary selection metrics for downstream analysis
        """
        self.logger.info("Calculating preliminary selection metrics")

        # Group by gene and methylation class
        gene_groups = mutations_df.groupby(['gene_symbol', 'methylation_class', 'mutation_origin'])

        selection_data = []

        for (gene, meth_class, origin), group in gene_groups:
            if len(group) < 2:  # Need minimum mutations for meaningful analysis
                continue

            # Count mutation types
            synonymous_count = len(group[group['mutation_type_category'] == 'low_impact'])
            nonsynonymous_count = len(group[group['mutation_type_category'].isin(['moderate_impact', 'high_impact'])])

            # Basic dN/dS approximation (proper calculation needs codon context)
            if synonymous_count > 0:
                dn_ds_ratio = nonsynonymous_count / synonymous_count
            else:
                dn_ds_ratio = np.inf if nonsynonymous_count > 0 else np.nan

            # Average pathogenicity scores
            avg_cadd = group['cadd_score'].mean()
            avg_gerp = group['gerp_score'].mean()

            # Mutation density (mutations per base pair, approximated)
            # This would need gene length information for accuracy
            mutation_density = len(group)  # Placeholder

            selection_data.append({
                'gene_symbol': gene,
                'methylation_class': meth_class,
                'mutation_origin': origin,
                'total_mutations': len(group),
                'synonymous_count': synonymous_count,
                'nonsynonymous_count': nonsynonymous_count,
                'dn_ds_ratio': dn_ds_ratio,
                'avg_cadd_score': avg_cadd,
                'avg_gerp_score': avg_gerp,
                'mutation_density': mutation_density,
                'high_impact_count': len(group[group['mutation_type_category'] == 'high_impact']),
                'moderate_impact_count': len(group[group['mutation_type_category'] == 'moderate_impact']),
                'low_impact_count': len(group[group['mutation_type_category'] == 'low_impact'])
            })

        selection_df = pd.DataFrame(selection_data)

        # Add selection metrics to mutation dataframe as gene-level annotations
        if len(selection_df) > 0:
            gene_metrics = selection_df.groupby(['gene_symbol', 'methylation_class']).agg({
                'dn_ds_ratio': 'mean',
                'avg_cadd_score': 'mean',
                'avg_gerp_score': 'mean',
                'total_mutations': 'sum'
            }).reset_index()

            mutations_df = mutations_df.merge(
                gene_metrics,
                on=['gene_symbol', 'methylation_class'],
                how='left',
                suffixes=('', '_gene')
            )

        # Save selection metrics separately
        mutations_df.attrs['selection_metrics'] = selection_df

        return mutations_df

    def save_annotated_mutations(self, mutations_df: pd.DataFrame, output_path: str):
        """
        Save annotated mutations with all metadata
        """
        # Main mutations file
        mutations_df.to_csv(output_path, sep='\t', index=False)

        # Save burden statistics
        if hasattr(mutations_df, 'attrs') and 'burden_stats' in mutations_df.attrs:
            burden_path = output_path.replace('.tsv', '_burden_stats.txt')
            with open(burden_path, 'w') as f:
                f.write("Mutation Burden Statistics\n")
                f.write("=" * 50 + "\n\n")
                for key, value in mutations_df.attrs['burden_stats'].items():
                    f.write(f"{key}: {value}\n")

        # Save selection metrics
        if hasattr(mutations_df, 'attrs') and 'selection_metrics' in mutations_df.attrs:
            selection_path = output_path.replace('.tsv', '_selection_metrics.tsv')
            mutations_df.attrs['selection_metrics'].to_csv(selection_path, sep='\t', index=False)

        self.logger.info(f"Annotated mutations saved to {output_path}")

        # Print summary
        self._print_annotation_summary(mutations_df)

    def _print_annotation_summary(self, mutations_df: pd.DataFrame):
        """Print summary of annotation results"""
        print("\nMutation Annotation Summary")
        print("=" * 50)

        print(f"Total mutations: {len(mutations_df)}")

        # By origin
        origin_counts = mutations_df['mutation_origin'].value_counts()
        print(f"\nBy origin:")
        for origin, count in origin_counts.items():
            print(f"  {origin}: {count} ({count / len(mutations_df) * 100:.1f}%)")

        # By methylation class
        meth_counts = mutations_df['methylation_class'].value_counts()
        print(f"\nBy methylation class:")
        for meth_class, count in meth_counts.items():
            print(f"  {meth_class}: {count} ({count / len(mutations_df) * 100:.1f}%)")

        # By functional category
        func_counts = mutations_df['mutation_type_category'].value_counts()
        print(f"\nBy functional impact:")
        for func_cat, count in func_counts.items():
            print(f"  {func_cat}: {count} ({count / len(mutations_df) * 100:.1f}%)")

        # Cross-tabulation
        print(f"\nMutation origin vs Methylation class:")
        crosstab = pd.crosstab(mutations_df['mutation_origin'], mutations_df['methylation_class'])
        print(crosstab)


def main():
    parser = argparse.ArgumentParser(
        description='Annotate mutations with methylation states and functional consequences')
    parser.add_argument('--germline', help='Germline mutations file (VCF or TSV)')
    parser.add_argument('--somatic', help='Somatic mutations file (VCF or TSV)')
    parser.add_argument('--germline-format', default='vcf', choices=['vcf', 'table'],
                        help='Germline mutations file format')
    parser.add_argument('--somatic-format', default='vcf', choices=['vcf', 'table'],
                        help='Somatic mutations file format')
    parser.add_argument('--methylation-regions', required=True,
                        help='Classified methylation regions file from methylation_classify.py')
    parser.add_argument('--output', required=True, help='Output annotated mutations file')

    # Functional annotation sources
    parser.add_argument('--vep-annotations', help='VEP annotations file')
    parser.add_argument('--annovar-annotations', help='ANNOVAR annotations file')
    parser.add_argument('--cadd-scores', help='CADD scores file')
    parser.add_argument('--conservation-scores', help='Conservation scores file')

    # Analysis options
    parser.add_argument('--calculate-burdens', action='store_true',
                        help='Calculate mutation burden metrics')
    parser.add_argument('--calculate-selection', action='store_true',
                        help='Calculate preliminary selection metrics')

    args = parser.parse_args()

    if not args.germline and not args.somatic:
        raise ValueError("Must provide either germline or somatic mutations (or both)")

    # Initialize annotator
    annotator = MutationAnnotator()

    try:
        mutations_list = []

        # Load germline mutations
        if args.germline:
            if args.germline_format == 'vcf':
                germline_df = annotator.load_vcf_mutations(args.germline, 'germline')
            else:
                germline_df = annotator.load_mutation_table(args.germline, 'germline')
            mutations_list.append(germline_df)

        # Load somatic mutations
        if args.somatic:
            if args.somatic_format == 'vcf':
                somatic_df = annotator.load_vcf_mutations(args.somatic, 'somatic')
            else:
                somatic_df = annotator.load_mutation_table(args.somatic, 'somatic')
            mutations_list.append(somatic_df)

        # Combine mutations
        mutations_df = pd.concat(mutations_list, ignore_index=True)

        # Annotate with methylation states
        mutations_df = annotator.annotate_with_methylation(mutations_df, args.methylation_regions)

        # Prepare functional annotation sources
        annotation_sources = {}
        if args.vep_annotations:
            annotation_sources['vep'] = args.vep_annotations
        if args.annovar_annotations:
            annotation_sources['annovar'] = args.annovar_annotations
        if args.cadd_scores:
            annotation_sources['cadd'] = args.cadd_scores
        if args.conservation_scores:
            annotation_sources['conservation'] = args.conservation_scores

        # Annotate functional consequences
        if annotation_sources:
            mutations_df = annotator.annotate_functional_consequences(mutations_df, annotation_sources)

        # Calculate mutation burdens if requested
        if args.calculate_burdens:
            mutations_df = annotator.calculate_mutation_burdens(mutations_df)

        # Calculate preliminary selection metrics if requested
        if args.calculate_selection:
            mutations_df = annotator.calculate_selection_metrics(mutations_df)

        # Save results
        annotator.save_annotated_mutations(mutations_df, args.output)

        print("Mutation annotation completed successfully!")

    except Exception as e:
        logging.error(f"Error in mutation annotation: {e}")
        raise


if __name__ == "__main__":
    main()

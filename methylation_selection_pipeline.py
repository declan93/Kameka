#!/usr/bin/env python3
"""
Updated Methylation-Selection Analysis Pipeline
With flexible configuration handling and better error management.
"""

import os
import sys
import argparse
import subprocess
import json
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
from datetime import datetime
import shutil


class FlexibleMethylationSelectionPipeline:
    """
    Enhanced pipeline with flexible configuration handling
    """

    def __init__(self, config_file: str, output_dir: str = "methylation_selection_analysis"):
        self.config_file = config_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Set up logging
        self.log_file = self.output_dir / "pipeline.log"
        self._setup_logging()

        # Load and normalize configuration
        self.config = self._load_and_normalize_config()

        # Track pipeline state
        self.pipeline_state = {
            'start_time': datetime.now().isoformat(),
            'steps_completed': [],
            'steps_failed': [],
            'intermediate_files': {},
            'final_results': {}
        }

        self.logger.info(f"Pipeline initialized. Output directory: {self.output_dir}")

    def _setup_logging(self):
        """Set up comprehensive logging"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def _load_and_normalize_config(self) -> Dict:
        """Load and normalize configuration from different formats"""
        self.logger.info(f"Loading configuration from {self.config_file}")

        with open(self.config_file, 'r') as f:
            if self.config_file.endswith('.yaml') or self.config_file.endswith('.yml'):
                raw_config = yaml.safe_load(f)
            else:
                raw_config = json.load(f)

        # Normalize configuration structure
        config = self._normalize_config_structure(raw_config)

        self.logger.info("Configuration loaded and normalized successfully")
        return config

    def _normalize_config_structure(self, raw_config: Dict) -> Dict:
        """Normalize different configuration formats to expected structure"""

        # Create normalized config with defaults
        config = {
            'input_files': {},
            'parameters': {},
            'output_settings': {
                'cleanup_intermediates': False,
                'generate_html_report': True,
                'create_plots': True
            }
        }

        # Handle input_files section
        if 'input_files' in raw_config:
            config['input_files'] = raw_config['input_files']

        # Handle parameters - check different possible section names
        parameters_section = None
        for section_name in ['parameters', 'analysis_parameters', 'config']:
            if section_name in raw_config:
                parameters_section = raw_config[section_name]
                break

        if parameters_section:
            # Flatten nested parameter structures if needed
            config['parameters'] = self._flatten_parameters(parameters_section)

        # Set parameter defaults
        default_params = {
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
            'min_mutations_per_gene': 5,
            'bootstrap_iterations': 1000,
            'permutation_tests': 10000,
            'fdr_method': 'BH',
            'effect_size_threshold': 0.1,
            'significance_threshold': 0.05,
            'confidence_level': 0.95,
            'cores': 4,
            'continue_on_failure': False
        }

        # Merge defaults with provided parameters
        for key, default_value in default_params.items():
            if key not in config['parameters']:
                config['parameters'][key] = default_value

        # Handle output_settings
        if 'output_settings' in raw_config:
            config['output_settings'].update(raw_config['output_settings'])

        # Store additional sections for reference
        for key, value in raw_config.items():
            if key not in ['input_files', 'parameters', 'analysis_parameters', 'output_settings']:
                config[key] = value

        return config

    def _flatten_parameters(self, parameters: Dict) -> Dict:
        """Flatten nested parameter structures"""
        flattened = {}

        for key, value in parameters.items():
            if isinstance(value, dict):
                # Handle nested dictionaries
                if key == 'methylation_thresholds':
                    flattened['hypo_threshold'] = value.get('hypomethylated', 0.3)
                    flattened['hyper_threshold'] = value.get('hypermethylated', 0.7)
                elif key == 'selection_analysis':
                    for subkey, subvalue in value.items():
                        # Convert camelCase/snake_case variations
                        normalized_key = subkey.replace('-', '_')
                        flattened[normalized_key] = subvalue
                elif key == 'quality_filters':
                    for subkey, subvalue in value.items():
                        flattened[f'min_{subkey}'] = subvalue
                else:
                    # Generic nested handling
                    for subkey, subvalue in value.items():
                        flattened[f'{key}_{subkey}'] = subvalue
            else:
                # Direct parameter
                flattened[key] = value

        return flattened

    def _validate_inputs(self) -> bool:
        """Validate all input files exist and are readable"""
        self.logger.info("Validating input files")

        input_files = self.config['input_files']
        missing_files = []

        # Check required files
        if 'methylation_data' not in input_files:
            missing_files.append("Missing required input: methylation_data")
        elif not os.path.exists(input_files['methylation_data']):
            missing_files.append(f"File not found: {input_files['methylation_data']}")

        # Check at least one mutation file exists
        mutation_files = ['germline_mutations', 'somatic_mutations']
        available_mutation_files = [key for key in mutation_files
                                    if key in input_files and os.path.exists(input_files[key])]

        if not available_mutation_files:
            missing_files.append("At least one mutation file (germline or somatic) is required")
        else:
            self.logger.info(f"Found mutation files: {available_mutation_files}")

        if missing_files:
            self.logger.error("Input validation failed:")
            for error in missing_files:
                self.logger.error(f"  - {error}")
            return False

        # Log available input files
        self.logger.info("Available input files:")
        for file_key, file_path in input_files.items():
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path) / 1024 / 1024  # MB
                self.logger.info(f"  ✓ {file_key}: {file_path} ({file_size:.1f} MB)")
            else:
                self.logger.warning(f"  - {file_key}: {file_path} (optional, not found)")

        return True

    def _run_command(self, command: List[str], step_name: str,
                     capture_output: bool = True, timeout: int = 3600) -> Tuple[bool, str]:
        """Execute a command with enhanced error handling and timeout"""
        self.logger.info(f"Starting step: {step_name}")
        self.logger.info(f"Command: {' '.join(command)}")

        try:
            if capture_output:
                result = subprocess.run(
                    command,
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    check=True
                )
                output = result.stdout
                if result.stderr:
                    self.logger.warning(f"Stderr: {result.stderr}")
            else:
                result = subprocess.run(command, timeout=timeout, check=True)
                output = "Command completed successfully"

            self.pipeline_state['steps_completed'].append(step_name)
            self.logger.info(f"Step completed successfully: {step_name}")
            return True, output

        except subprocess.TimeoutExpired as e:
            error_msg = f"Step timed out after {timeout}s: {step_name}"
            self.pipeline_state['steps_failed'].append(step_name)
            self.logger.error(error_msg)
            return False, error_msg

        except subprocess.CalledProcessError as e:
            error_msg = f"Step failed: {step_name}. Return code: {e.returncode}"
            if hasattr(e, 'stdout') and e.stdout:
                error_msg += f"\nStdout: {e.stdout}"
            if hasattr(e, 'stderr') and e.stderr:
                error_msg += f"\nStderr: {e.stderr}"

            self.pipeline_state['steps_failed'].append(step_name)
            self.logger.error(error_msg)
            return False, error_msg

        except FileNotFoundError as e:
            error_msg = f"Command not found for step {step_name}: {command[0]}. Error: {e}"
            self.pipeline_state['steps_failed'].append(step_name)
            self.logger.error(error_msg)
            return False, error_msg

    def step1_classify_methylation(self) -> bool:
        """Step 1: Classify methylation states"""
        self.logger.info("=" * 50)
        self.logger.info("STEP 1: METHYLATION CLASSIFICATION")
        self.logger.info("=" * 50)

        # Prepare output file
        output_file = self.output_dir / "methylation_classified.tsv"
        plot_file = self.output_dir / "methylation_distribution.png"

        # Build command
        command = [
            'python3', 'methylation_classify.py',
            '--input', self.config['input_files']['methylation_data'],
            '--output', str(output_file),
            '--format', self.config['parameters'].get('methylation_format', 'bedgraph'),
            '--method', self.config['parameters'].get('classification_method', 'threshold'),
            '--hypo-threshold', str(self.config['parameters'].get('hypo_threshold', 0.3)),
            '--hyper-threshold', str(self.config['parameters'].get('hyper_threshold', 0.7)),
            '--plot-output', str(plot_file)
        ]

        # Add optional annotation files
        input_files = self.config['input_files']
        optional_annotations = ['cpg_islands', 'promoters', 'enhancers']
        for annotation in optional_annotations:
            if annotation in input_files and os.path.exists(input_files[annotation]):
                command.extend([f'--{annotation.replace("_", "-")}', input_files[annotation]])

        if self.config['parameters'].get('regional_analysis', False):
            command.append('--regional-analysis')
            command.extend(['--window-size',
                            str(self.config['parameters'].get('window_size', 1000))])

        success, output = self._run_command(command, "methylation_classification")

        if success:
            self.pipeline_state['intermediate_files']['methylation_regions'] = str(output_file)

            # Validate output
            if output_file.exists() and output_file.stat().st_size > 0:
                self.logger.info(f"Methylation classification completed. Output: {output_file}")
                return True
            else:
                self.logger.error("Methylation classification output file is empty or missing")
                return False
        else:
            return False

    def step2_annotate_mutations(self) -> bool:
        """Step 2: Annotate mutations with methylation states and functional consequences"""
        self.logger.info("=" * 50)
        self.logger.info("STEP 2: MUTATION ANNOTATION")
        self.logger.info("=" * 50)

        # Check methylation regions file exists from step 1
        if 'methylation_regions' not in self.pipeline_state['intermediate_files']:
            self.logger.error("Methylation regions file not found. Step 1 may have failed.")
            return False

        # Prepare output file
        output_file = self.output_dir / "mutations_annotated.tsv"

        # Build command
        command = [
            'python3', 'annotate_mutations.py',
            '--methylation-regions', self.pipeline_state['intermediate_files']['methylation_regions'],
            '--output', str(output_file)
        ]

        # Add mutation files
        input_files = self.config['input_files']
        if 'germline_mutations' in input_files and os.path.exists(input_files['germline_mutations']):
            command.extend(['--germline', input_files['germline_mutations']])
            command.extend(['--germline-format',
                            self.config['parameters'].get('germline_format', 'vcf')])

        if 'somatic_mutations' in input_files and os.path.exists(input_files['somatic_mutations']):
            command.extend(['--somatic', input_files['somatic_mutations']])
            command.extend(['--somatic-format',
                            self.config['parameters'].get('somatic_format', 'vcf')])

        # Add functional annotation sources
        annotation_mappings = {
            'vep_annotations': '--vep-annotations',
            'annovar_annotations': '--annovar-annotations',
            'cadd_scores': '--cadd-scores',
            'conservation_scores': '--conservation-scores'
        }

        for ann_key, ann_flag in annotation_mappings.items():
            if ann_key in input_files and os.path.exists(input_files[ann_key]):
                command.extend([ann_flag, input_files[ann_key]])

        # Add analysis options
        if self.config['parameters'].get('calculate_burdens', True):
            command.append('--calculate-burdens')

        if self.config['parameters'].get('calculate_selection', True):
            command.append('--calculate-selection')

        success, output = self._run_command(command, "mutation_annotation", timeout=7200)  # 2 hours

        if success:
            self.pipeline_state['intermediate_files']['annotated_mutations'] = str(output_file)

            # Validate output
            if output_file.exists() and output_file.stat().st_size > 0:
                # Quick validation - check file structure
                try:
                    test_df = pd.read_csv(output_file, sep='\t', nrows=5)
                    self.logger.info(f"Mutation annotation completed. Output: {output_file}")
                    self.logger.info(f"Columns: {len(test_df.columns)}, sample rows: {len(test_df)}")
                    return True
                except Exception as e:
                    self.logger.error(f"Invalid annotation output file: {e}")
                    return False
            else:
                self.logger.error("Mutation annotation output file is empty or missing")
                return False
        else:
            return False

    def step3_calculate_selection(self) -> bool:
        """Step 3: Calculate comprehensive selection metrics"""
        self.logger.info("=" * 50)
        self.logger.info("STEP 3: SELECTION CALCULATION")
        self.logger.info("=" * 50)

        # Check annotated mutations file exists from step 2
        if 'annotated_mutations' not in self.pipeline_state['intermediate_files']:
            self.logger.error("Annotated mutations file not found. Step 2 may have failed.")
            return False

        # Prepare output prefix
        output_prefix = str(self.output_dir / "selection_analysis")

        # Build command
        command = [
            'Rscript', 'calculate_selection.r',  # Note: lowercase .r extension
            '--mutations', self.pipeline_state['intermediate_files']['annotated_mutations'],
            '--output', output_prefix,
            '--min-mutations', str(self.config['parameters'].get('min_mutations_per_gene', 5)),
            '--bootstrap-n', str(self.config['parameters'].get('bootstrap_iterations', 1000)),
            '--cores', str(self.config['parameters'].get('cores', 4)),
            '--alpha', str(self.config['parameters'].get('significance_threshold', 0.05))
        ]

        # Add optional files
        input_files = self.config['input_files']
        if 'gene_lengths' in input_files and os.path.exists(input_files['gene_lengths']):
            command.extend(['--gene-lengths', input_files['gene_lengths']])

        if 'codon_usage' in input_files and os.path.exists(input_files['codon_usage']):
            command.extend(['--codon-usage', input_files['codon_usage']])

        success, output = self._run_command(command, "selection_calculation",
                                            capture_output=False, timeout=7200)  # 2 hours

        if success:
            # Check for expected output files
            selection_metrics_file = Path(f"{output_prefix}_selection_metrics.tsv")
            rdata_file = Path(f"{output_prefix}_selection_analysis.RData")

            if selection_metrics_file.exists():
                self.pipeline_state['intermediate_files']['selection_metrics'] = str(selection_metrics_file)
                self.pipeline_state['final_results']['selection_analysis_data'] = str(rdata_file)
                self.logger.info(f"Selection calculation completed. Output: {selection_metrics_file}")
                return True
            else:
                self.logger.error("Selection calculation output file not found")
                # List what files were actually created
                self.logger.info("Files in output directory:")
                for f in self.output_dir.glob("*selection*"):
                    self.logger.info(f"  {f}")
                return False
        else:
            return False

    def step4_generate_summary(self) -> bool:
        """Step 4: Generate comprehensive summary and quality control"""
        self.logger.info("=" * 50)
        self.logger.info("STEP 4: SUMMARY AND QUALITY CONTROL")
        self.logger.info("=" * 50)

        summary_file = self.output_dir / "analysis_summary.txt"
        qc_file = self.output_dir / "quality_control.txt"

        try:
            # Generate analysis summary
            with open(summary_file, 'w') as f:
                f.write("METHYLATION-SELECTION ANALYSIS SUMMARY\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Configuration: {self.config_file}\n")
                f.write(f"Output Directory: {self.output_dir}\n\n")

                # Input summary
                f.write("INPUT FILES\n")
                f.write("-" * 20 + "\n")
                for file_key, file_path in self.config['input_files'].items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path) / 1024 / 1024
                        f.write(f"✓ {file_key}: {os.path.basename(file_path)} ({file_size:.1f} MB)\n")
                    else:
                        f.write(f"- {file_key}: not provided\n")

                # Pipeline status
                f.write(f"\nPIPELINE STATUS\n")
                f.write("-" * 20 + "\n")
                f.write(f"Steps completed: {len(self.pipeline_state['steps_completed'])}\n")
                f.write(f"Steps failed: {len(self.pipeline_state['steps_failed'])}\n")

                if self.pipeline_state['steps_completed']:
                    f.write("\nCompleted steps:\n")
                    for step in self.pipeline_state['steps_completed']:
                        f.write(f"  ✓ {step}\n")

                if self.pipeline_state['steps_failed']:
                    f.write("\nFailed steps:\n")
                    for step in self.pipeline_state['steps_failed']:
                        f.write(f"  ✗ {step}\n")

                # Results summary
                f.write(f"\nOUTPUT FILES\n")
                f.write("-" * 20 + "\n")

                all_outputs = {**self.pipeline_state['intermediate_files'],
                               **self.pipeline_state['final_results']}

                for file_key, file_path in all_outputs.items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path) / 1024 / 1024
                        rel_path = os.path.relpath(file_path, self.output_dir)
                        f.write(f"✓ {file_key}: {rel_path} ({file_size:.1f} MB)\n")

                # Quick data summary if annotated mutations exist
                if 'annotated_mutations' in self.pipeline_state['intermediate_files']:
                    ann_file = self.pipeline_state['intermediate_files']['annotated_mutations']
                    try:
                        df = pd.read_csv(ann_file, sep='\t')
                        f.write(f"\nDATA SUMMARY\n")
                        f.write("-" * 20 + "\n")
                        f.write(f"Total mutations: {len(df):,}\n")

                        if 'mutation_origin' in df.columns:
                            f.write("By origin:\n")
                            for origin, count in df['mutation_origin'].value_counts().items():
                                f.write(f"  {origin}: {count:,}\n")

                        if 'methylation_class' in df.columns:
                            f.write("By methylation class:\n")
                            for mclass, count in df['methylation_class'].value_counts().items():
                                f.write(f"  {mclass}: {count:,}\n")

                        if 'gene_symbol' in df.columns:
                            unique_genes = df['gene_symbol'].nunique()
                            f.write(f"Unique genes: {unique_genes:,}\n")

                    except Exception as e:
                        f.write(f"\nCould not read mutation data for summary: {e}\n")

            # Generate quality control report
            with open(qc_file, 'w') as f:
                f.write("QUALITY CONTROL REPORT\n")
                f.write("=" * 30 + "\n\n")

                # Check file integrity
                f.write("FILE INTEGRITY CHECK\n")
                f.write("-" * 25 + "\n")

                all_files = {**self.config['input_files'],
                             **self.pipeline_state['intermediate_files'],
                             **self.pipeline_state['final_results']}

                for file_key, file_path in all_files.items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path)
                        if file_size > 0:
                            f.write(f"✓ {file_key}: OK ({file_size:,} bytes)\n")
                        else:
                            f.write(f"⚠ {file_key}: Empty file\n")
                    else:
                        f.write(f"✗ {file_key}: Missing\n")

                # Runtime information
                start_time = datetime.fromisoformat(self.pipeline_state['start_time'])
                end_time = datetime.now()
                runtime = end_time - start_time

                f.write(f"\nRUNTIME INFORMATION\n")
                f.write("-" * 25 + "\n")
                f.write(f"Start: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"End: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total runtime: {runtime}\n")
                f.write(f"Log file: {self.log_file}\n")

            self.logger.info(f"Summary generated: {summary_file}")
            self.logger.info(f"Quality control report: {qc_file}")

            return True

        except Exception as e:
            self.logger.error(f"Failed to generate summary: {e}")
            return False

    def run_pipeline(self) -> bool:
        """Execute the complete pipeline with flexible steps"""
        self.logger.info("Starting flexible methylation-selection analysis pipeline")

        # Validate inputs
        if not self._validate_inputs():
            self.logger.error("Input validation failed. Pipeline aborted.")
            return False

        # Execute pipeline steps
        steps = [
            ("methylation_classification", self.step1_classify_methylation),
            ("mutation_annotation", self.step2_annotate_mutations),
            ("selection_calculation", self.step3_calculate_selection),
            ("summary_generation", self.step4_generate_summary)
        ]

        pipeline_success = True
        continue_on_failure = self.config['parameters'].get('continue_on_failure', False)

        for step_name, step_function in steps:
            self.logger.info(f"\n{'=' * 60}")
            self.logger.info(f"EXECUTING: {step_name.upper().replace('_', ' ')}")
            self.logger.info(f"{'=' * 60}")

            if not step_function():
                self.logger.error(f"Step failed: {step_name}")
                pipeline_success = False

                if not continue_on_failure:
                    self.logger.error("Pipeline aborted due to step failure")
                    break
                else:
                    self.logger.warning("Continuing despite step failure")
            else:
                self.logger.info(f"Step completed successfully: {step_name}")

        # Always generate final summary (even on partial failure)
        self.step4_generate_summary()

        # Final status
        if pipeline_success:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info("=" * 60)
            self.logger.info(f"Results available in: {self.output_dir}")
            self.logger.info("Key output files:")
            for key, path in self.pipeline_state['intermediate_files'].items():
                self.logger.info(f"  {key}: {os.path.relpath(path, self.output_dir)}")
        else:
            self.logger.error("\n" + "=" * 60)
            self.logger.error("PIPELINE COMPLETED WITH ERRORS")
            self.logger.error("=" * 60)
            self.logger.error(f"Check logs and partial results in: {self.output_dir}")
            self.logger.error(f"Log file: {self.log_file}")

        return pipeline_success


def main():
    parser = argparse.ArgumentParser(
        description='Flexible Methylation-Selection Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run pipeline with test data
  python flexible_pipeline.py --config comprehensive_test_config.yaml --output test_results/

  # Run pipeline with timeout and continue on failure
  python flexible_pipeline.py --config config.yaml --output results/ --continue-on-failure
        """
    )

    parser.add_argument('--config', type=str, required=True,
                        help='Configuration file (YAML or JSON)')
    parser.add_argument('--output', type=str, default='methylation_selection_results',
                        help='Output directory')
    parser.add_argument('--validate-only', action='store_true',
                        help='Only validate inputs without running pipeline')
    parser.add_argument('--continue-on-failure', action='store_true',
                        help='Continue pipeline even if some steps fail')

    args = parser.parse_args()

    try:
        # Initialize pipeline
        pipeline = FlexibleMethylationSelectionPipeline(args.config, args.output)

        # Override continue_on_failure if specified
        if args.continue_on_failure:
            pipeline.config['parameters']['continue_on_failure'] = True

        if args.validate_only:
            success = pipeline._validate_inputs()
            if success:
                print("✓ Input validation passed")
                print("\nConfiguration summary:")
                print(f"Methylation data: {pipeline.config['input_files']['methylation_data']}")
                for key in ['germline_mutations', 'somatic_mutations']:
                    if key in pipeline.config['input_files']:
                        print(f"{key}: {pipeline.config['input_files'][key]}")
                sys.exit(0)
            else:
                print("✗ Input validation failed")
                sys.exit(1)

        # Run complete pipeline
        success = pipeline.run_pipeline()

        if success:
            print(f"\n✓ Pipeline completed successfully!")
            print(f"Results available in: {pipeline.output_dir}")
            print(f"Summary: {pipeline.output_dir}/analysis_summary.txt")
            sys.exit(0)
        else:
            print(f"\n⚠ Pipeline completed with some errors")
            print(f"Check results and logs in: {pipeline.output_dir}")
            print(f"Log file: {pipeline.output_dir}/pipeline.log")
            sys.exit(1)

    except Exception as e:
        print(f"Pipeline error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
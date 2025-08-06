#!/usr/bin/env python3
"""
Methylation-Selection Analysis Pipeline
Complete pipeline orchestration for methylation-dependent selection analysis.
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

class MethylationSelectionPipeline:
    """
    Complete pipeline orchestration for methylation-dependent selection analysis
    """
    
    def __init__(self, config_file: str, output_dir: str = "methylation_selection_analysis"):
        self.config_file = config_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logging
        self.log_file = self.output_dir / "pipeline.log"
        self._setup_logging()
        
        # Load configuration
        self.config = self._load_config()
        
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
        
    def _load_config(self) -> Dict:
        """Load and validate configuration file"""
        self.logger.info(f"Loading configuration from {self.config_file}")
        
        with open(self.config_file, 'r') as f:
            if self.config_file.endswith('.yaml') or self.config_file.endswith('.yml'):
                config = yaml.safe_load(f)
            else:
                config = json.load(f)
        
        # Validate required sections
        required_sections = ['input_files', 'parameters', 'output_settings']
        missing_sections = [section for section in required_sections if section not in config]
        if missing_sections:
            raise ValueError(f"Missing required config sections: {missing_sections}")
        
        self.logger.info("Configuration loaded and validated successfully")
        return config
    
    def _validate_inputs(self) -> bool:
        """Validate all input files exist and are readable"""
        self.logger.info("Validating input files")
        
        input_files = self.config['input_files']
        missing_files = []
        
        # Check required files
        required_files = ['methylation_data']
        for file_key in required_files:
            if file_key not in input_files:
                missing_files.append(f"Missing required input: {file_key}")
            elif not os.path.exists(input_files[file_key]):
                missing_files.append(f"File not found: {input_files[file_key]}")
        
        # Check at least one mutation file exists
        mutation_files = ['germline_mutations', 'somatic_mutations']
        if not any(key in input_files and os.path.exists(input_files[key]) 
                  for key in mutation_files):
            missing_files.append("At least one mutation file (germline or somatic) is required")
        
        if missing_files:
            self.logger.error("Input validation failed:")
            for error in missing_files:
                self.logger.error(f"  - {error}")
            return False
        
        self.logger.info("All input files validated successfully")
        return True
    
    def _run_command(self, command: List[str], step_name: str, 
                    capture_output: bool = True) -> Tuple[bool, str]:
        """Execute a command and handle logging"""
        self.logger.info(f"Starting step: {step_name}")
        self.logger.info(f"Command: {' '.join(command)}")
        
        try:
            if capture_output:
                result = subprocess.run(
                    command, 
                    capture_output=True, 
                    text=True, 
                    check=True
                )
                output = result.stdout
                if result.stderr:
                    self.logger.warning(f"Stderr: {result.stderr}")
            else:
                result = subprocess.run(command, check=True)
                output = "Command completed successfully"
            
            self.pipeline_state['steps_completed'].append(step_name)
            self.logger.info(f"Step completed successfully: {step_name}")
            return True, output
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Step failed: {step_name}. Error: {e}"
            if hasattr(e, 'stderr') and e.stderr:
                error_msg += f"\nStderr: {e.stderr}"
            
            self.pipeline_state['steps_failed'].append(step_name)
            self.logger.error(error_msg)
            return False, error_msg
    
    def step1_classify_methylation(self) -> bool:
        """Step 1: Classify methylation states"""
        self.logger.info("="*50)
        self.logger.info("STEP 1: METHYLATION CLASSIFICATION")
        self.logger.info("="*50)
        
        # Prepare output file
        output_file = self.output_dir / "methylation_classified.bed"
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
        if 'cpg_islands' in input_files:
            command.extend(['--cpg-islands', input_files['cpg_islands']])
        if 'promoters' in input_files:
            command.extend(['--promoters', input_files['promoters']])
        if 'enhancers' in input_files:
            command.extend(['--enhancers', input_files['enhancers']])
        
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
        self.logger.info("="*50)
        self.logger.info("STEP 2: MUTATION ANNOTATION")
        self.logger.info("="*50)
        
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
        if 'germline_mutations' in input_files:
            command.extend(['--germline', input_files['germline_mutations']])
            command.extend(['--germline-format', 
                          self.config['parameters'].get('germline_format', 'vcf')])
        
        if 'somatic_mutations' in input_files:
            command.extend(['--somatic', input_files['somatic_mutations']])
            command.extend(['--somatic-format', 
                          self.config['parameters'].get('somatic_format', 'vcf')])
        
        # Add functional annotation sources
        annotation_files = ['vep_annotations', 'annovar_annotations', 
                          'cadd_scores', 'conservation_scores']
        for ann_file in annotation_files:
            if ann_file in input_files:
                command.extend([f'--{ann_file.replace("_", "-")}', input_files[ann_file]])
        
        # Add analysis options
        if self.config['parameters'].get('calculate_burdens', True):
            command.append('--calculate-burdens')
        
        if self.config['parameters'].get('calculate_selection', True):
            command.append('--calculate-selection')
        
        success, output = self._run_command(command, "mutation_annotation")
        
        if success:
            self.pipeline_state['intermediate_files']['annotated_mutations'] = str(output_file)
            
            # Validate output
            if output_file.exists() and output_file.stat().st_size > 0:
                self.logger.info(f"Mutation annotation completed. Output: {output_file}")
                return True
            else:
                self.logger.error("Mutation annotation output file is empty or missing")
                return False
        else:
            return False
    
    def step3_calculate_selection(self) -> bool:
        """Step 3: Calculate comprehensive selection metrics"""
        self.logger.info("="*50)
        self.logger.info("STEP 3: SELECTION CALCULATION")
        self.logger.info("="*50)
        
        # Check annotated mutations file exists from step 2
        if 'annotated_mutations' not in self.pipeline_state['intermediate_files']:
            self.logger.error("Annotated mutations file not found. Step 2 may have failed.")
            return False
        
        # Prepare output prefix
        output_prefix = str(self.output_dir / "selection_analysis")
        
        # Build command
        command = [
            'Rscript', 'calculate_selection.R',
            '--mutations', self.pipeline_state['intermediate_files']['annotated_mutations'],
            '--output', output_prefix,
            '--min-mutations', str(self.config['parameters'].get('min_mutations_per_gene', 5)),
            '--bootstrap-n', str(self.config['parameters'].get('bootstrap_iterations', 1000)),
            '--cores', str(self.config['parameters'].get('cores', 4)),
            '--alpha', str(self.config['parameters'].get('significance_threshold', 0.05))
        ]
        
        # Add optional files
        input_files = self.config['input_files']
        if 'gene_lengths' in input_files:
            command.extend(['--gene-lengths', input_files['gene_lengths']])
        
        if 'codon_usage' in input_files:
            command.extend(['--codon-usage', input_files['codon_usage']])
        
        success, output = self._run_command(command, "selection_calculation", capture_output=False)
        
        if success:
            # Check for expected output files
            selection_metrics_file = Path(f"{output_prefix}_selection_metrics.tsv")
            if selection_metrics_file.exists():
                self.pipeline_state['intermediate_files']['selection_metrics'] = str(selection_metrics_file)
                self.logger.info(f"Selection calculation completed. Output: {selection_metrics_file}")
                return True
            else:
                self.logger.error("Selection calculation output file not found")
                return False
        else:
            return False
    
    def step4_compare_selection(self) -> bool:
        """Step 4: Comprehensive selection comparison and statistical analysis"""
        self.logger.info("="*50)
        self.logger.info("STEP 4: SELECTION COMPARISON")
        self.logger.info("="*50)
        
        # Check selection metrics file exists from step 3
        if 'selection_metrics' not in self.pipeline_state['intermediate_files']:
            self.logger.error("Selection metrics file not found. Step 3 may have failed.")
            return False
        
        # Prepare output prefix
        output_prefix = str(self.output_dir / "comprehensive_analysis")
        
        # Build command
        command = [
            'Rscript', 'compare_selection.R',
            '--selection-scores', self.pipeline_state['intermediate_files']['selection_metrics'],
            '--output', output_prefix,
            '--permutations', str(self.config['parameters'].get('permutation_tests', 10000)),
            '--fdr-method', self.config['parameters'].get('fdr_method', 'BH'),
            '--cores', str(self.config['parameters'].get('cores', 4)),
            '--effect-size-threshold', str(self.config['parameters'].get('effect_size_threshold', 0.1)),
            '--confidence-level', str(self.config['parameters'].get('confidence_level', 0.95))
        ]
        
        # Add mutations file for enhanced analysis
        if 'annotated_mutations' in self.pipeline_state['intermediate_files']:
            command.extend(['--mutations', self.pipeline_state['intermediate_files']['annotated_mutations']])
        
        # Add HTML report generation if requested
        if self.config['output_settings'].get('generate_html_report', True):
            command.append('--generate-report')
        
        success, output = self._run_command(command, "selection_comparison", capture_output=False)
        
        if success:
            # Check for expected output files
            comprehensive_results = Path(f"{output_prefix}_comprehensive_results.RData")
            plots_file = Path(f"{output_prefix}_comprehensive_analysis.png")
            
            if comprehensive_results.exists():
                self.pipeline_state['final_results']['comprehensive_analysis'] = str(comprehensive_results)
                
                if plots_file.exists():
                    self.pipeline_state['final_results']['comprehensive_plots'] = str(plots_file)
                
                self.logger.info(f"Selection comparison completed. Results: {comprehensive_results}")
                return True
            else:
                self.logger.error("Selection comparison output files not found")
                return False
        else:
            return False
    
    def step5_quality_control(self) -> bool:
        """Step 5: Generate quality control report and validation"""
        self.logger.info("="*50)
        self.logger.info("STEP 5: QUALITY CONTROL AND VALIDATION")
        self.logger.info("="*50)
        
        qc_report = self.output_dir / "quality_control_report.txt"
        
        try:
            with open(qc_report, 'w') as f:
                f.write("METHYLATION-SELECTION ANALYSIS QUALITY CONTROL REPORT\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Pipeline Version: 1.0.0\n")
                f.write(f"Output Directory: {self.output_dir}\n\n")
                
                # Input validation
                f.write("INPUT VALIDATION\n")
                f.write("-" * 20 + "\n")
                
                input_files = self.config['input_files']
                for file_key, file_path in input_files.items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path)
                        f.write(f"✓ {file_key}: {file_path} ({file_size:,} bytes)\n")
                    else:
                        f.write(f"✗ {file_key}: {file_path} (NOT FOUND)\n")
                
                f.write("\nPIPELINE EXECUTION\n")
                f.write("-" * 20 + "\n")
                
                total_steps = len(self.pipeline_state['steps_completed']) + len(self.pipeline_state['steps_failed'])
                completed_steps = len(self.pipeline_state['steps_completed'])
                
                f.write(f"Total steps: {total_steps}\n")
                f.write(f"Completed steps: {completed_steps}\n")
                f.write(f"Failed steps: {len(self.pipeline_state['steps_failed'])}\n")
                
                if self.pipeline_state['steps_completed']:
                    f.write("\nCompleted steps:\n")
                    for step in self.pipeline_state['steps_completed']:
                        f.write(f"  ✓ {step}\n")
                
                if self.pipeline_state['steps_failed']:
                    f.write("\nFailed steps:\n")
                    for step in self.pipeline_state['steps_failed']:
                        f.write(f"  ✗ {step}\n")
                
                # Output validation
                f.write("\nOUTPUT VALIDATION\n")
                f.write("-" * 20 + "\n")
                
                # Check intermediate files
                for file_key, file_path in self.pipeline_state['intermediate_files'].items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path)
                        f.write(f"✓ {file_key}: {file_path} ({file_size:,} bytes)\n")
                        
                        # Additional validation for specific file types
                        if file_key == 'annotated_mutations':
                            try:
                                df = pd.read_csv(file_path, sep='\t', nrows=5)
                                f.write(f"    Columns: {len(df.columns)}, Sample rows: {len(df)}\n")
                            except Exception as e:
                                f.write(f"    ✗ Error reading file: {e}\n")
                    else:
                        f.write(f"✗ {file_key}: {file_path} (NOT FOUND)\n")
                
                # Check final results
                for file_key, file_path in self.pipeline_state['final_results'].items():
                    if os.path.exists(file_path):
                        file_size = os.path.getsize(file_path)
                        f.write(f"✓ {file_key}: {file_path} ({file_size:,} bytes)\n")
                    else:
                        f.write(f"✗ {file_key}: {file_path} (NOT FOUND)\n")
                
                # Performance metrics
                f.write("\nPERFORMANCE METRICS\n")
                f.write("-" * 20 + "\n")
                
                start_time = datetime.fromisoformat(self.pipeline_state['start_time'])
                end_time = datetime.now()
                total_runtime = end_time - start_time
                
                f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Total runtime: {total_runtime}\n")
                
                # Log file info
                if os.path.exists(self.log_file):
                    log_size = os.path.getsize(self.log_file)
                    f.write(f"Log file size: {log_size:,} bytes\n")
            
            self.logger.info(f"Quality control report generated: {qc_report}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to generate quality control report: {e}")
            return False
    
    def cleanup_intermediate_files(self):
        """Clean up intermediate files if requested"""
        if self.config['output_settings'].get('cleanup_intermediates', False):
            self.logger.info("Cleaning up intermediate files")
            
            # Keep only final results and key intermediate files
            files_to_keep = {
                'annotated_mutations',  # Keep for potential reanalysis
                'selection_metrics'     # Keep for potential reanalysis
            }
            
            for file_key, file_path in self.pipeline_state['intermediate_files'].items():
                if file_key not in files_to_keep:
                    try:
                        if os.path.exists(file_path):
                            os.remove(file_path)
                            self.logger.info(f"Removed intermediate file: {file_path}")
                    except Exception as e:
                        self.logger.warning(f"Could not remove {file_path}: {e}")
    
    def generate_summary_report(self):
        """Generate final summary report"""
        summary_file = self.output_dir / "ANALYSIS_SUMMARY.txt"
        
        with open(summary_file, 'w') as f:
            f.write("METHYLATION-SELECTION ANALYSIS SUMMARY\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Configuration file: {self.config_file}\n")
            f.write(f"Output directory: {self.output_dir}\n\n")
            
            # Success/failure summary
            total_steps = 5  # We have 5 main steps
            completed_steps = len(self.pipeline_state['steps_completed'])
            success_rate = completed_steps / total_steps * 100
            
            f.write(f"PIPELINE STATUS\n")
            f.write(f"Steps completed: {completed_steps}/{total_steps} ({success_rate:.1f}%)\n")
            
            if len(self.pipeline_state['steps_failed']) == 0:
                f.write("✓ All steps completed successfully\n\n")
            else:
                f.write(f"✗ {len(self.pipeline_state['steps_failed'])} steps failed\n\n")
            
            # Key results
            f.write("KEY OUTPUT FILES\n")
            f.write("-" * 20 + "\n")
            
            # Final results
            for result_type, file_path in self.pipeline_state['final_results'].items():
                rel_path = os.path.relpath(file_path, self.output_dir)
                f.write(f"• {result_type}: {rel_path}\n")
            
            # Important intermediate files
            important_intermediates = ['annotated_mutations', 'selection_metrics']
            for file_key in important_intermediates:
                if file_key in self.pipeline_state['intermediate_files']:
                    file_path = self.pipeline_state['intermediate_files'][file_key]
                    rel_path = os.path.relpath(file_path, self.output_dir)
                    f.write(f"• {file_key}: {rel_path}\n")
            
            f.write("\nFor detailed results, see:\n")
            f.write("• comprehensive_analysis_report.html (if generated)\n")
            f.write("• quality_control_report.txt\n")
            f.write("• pipeline.log\n")
        
        self.logger.info(f"Summary report generated: {summary_file}")
    
    def run_pipeline(self) -> bool:
        """Execute the complete pipeline"""
        self.logger.info("Starting methylation-selection analysis pipeline")
        
        # Validate inputs
        if not self._validate_inputs():
            self.logger.error("Input validation failed. Pipeline aborted.")
            return False
        
        # Execute pipeline steps
        steps = [
            ("methylation_classification", self.step1_classify_methylation),
            ("mutation_annotation", self.step2_annotate_mutations),
            ("selection_calculation", self.step3_calculate_selection),
            ("selection_comparison", self.step4_compare_selection),
            ("quality_control", self.step5_quality_control)
        ]
        
        pipeline_success = True
        
        for step_name, step_function in steps:
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"EXECUTING: {step_name.upper().replace('_', ' ')}")
            self.logger.info(f"{'='*60}")
            
            if not step_function():
                self.logger.error(f"Step failed: {step_name}")
                pipeline_success = False
                
                # Check if we should continue on failure
                if not self.config['parameters'].get('continue_on_failure', False):
                    self.logger.error("Pipeline aborted due to step failure")
                    break
            else:
                self.logger.info(f"Step completed successfully: {step_name}")
        
        # Cleanup if requested
        if self.config['output_settings'].get('cleanup_intermediates', False):
            self.cleanup_intermediate_files()
        
        # Generate final summary
        self.generate_summary_report()
        
        # Final status
        if pipeline_success:
            self.logger.info("\n" + "="*60)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info("="*60)
            self.logger.info(f"Results available in: {self.output_dir}")
        else:
            self.logger.error("\n" + "="*60)
            self.logger.error("PIPELINE COMPLETED WITH ERRORS")
            self.logger.error("="*60)
            self.logger.error(f"Check logs and partial results in: {self.output_dir}")
        
        return pipeline_success

def create_example_config(output_path: str = "example_config.yaml"):
    """Create an example configuration file"""
    example_config = {
        'input_files': {
            'methylation_data': '/path/to/methylation_data.bedgraph',
            'germline_mutations': '/path/to/germline_variants.vcf',
            'somatic_mutations': '/path/to/somatic_variants.vcf',
            'cpg_islands': '/path/to/cpg_islands.bed',
            'promoters': '/path/to/promoters.bed',
            'enhancers': '/path/to/enhancers.bed',
            'vep_annotations': '/path/to/vep_annotations.txt',
            'cadd_scores': '/path/to/cadd_scores.tsv',
            'gene_lengths': '/path/to/gene_lengths.tsv'
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
            'min_mutations_per_gene': 5,
            'bootstrap_iterations': 1000,
            'permutation_tests': 10000,
            'fdr_method': 'BH',
            'effect_size_threshold': 0.1,
            'significance_threshold': 0.05,
            'confidence_level': 0.95,
            'cores': 4,
            'continue_on_failure': False
        },
        'output_settings': {
            'cleanup_intermediates': False,
            'generate_html_report': True,
            'create_plots': True
        }
    }
    
    with open(output_path, 'w') as f:
        yaml.dump(example_config, f, default_flow_style=False, indent=2)
    
    print(f"Example configuration file created: {output_path}")
    print("Please edit the file paths and parameters as needed for your analysis.")

def main():
    parser = argparse.ArgumentParser(
        description='Methylation-Selection Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create example configuration file
  python methylation_selection_pipeline.py --create-config
  
  # Run complete pipeline
  python methylation_selection_pipeline.py --config config.yaml --output results/
  
  # Run pipeline with custom output directory
  python methylation_selection_pipeline.py --config my_config.yaml --output /path/to/results/
        """
    )
    
    parser.add_argument('--config', type=str, 
                       help='Configuration file (YAML or JSON)')
    parser.add_argument('--output', type=str, default='methylation_selection_analysis',
                       help='Output directory')
    parser.add_argument('--create-config', action='store_true',
                       help='Create example configuration file')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate inputs without running pipeline')
    
    args = parser.parse_args()
    
    if args.create_config:
        create_example_config()
        return
    
    if not args.config:
        parser.error("Configuration file is required (use --config or --create-config)")
    
    try:
        # Initialize pipeline
        pipeline = MethylationSelectionPipeline(args.config, args.output)
        
        if args.validate_only:
            success = pipeline._validate_inputs()
            if success:
                print("✓ Input validation passed")
                sys.exit(0)
            else:
                print("✗ Input validation failed")
                sys.exit(1)
        
        # Run complete pipeline
        success = pipeline.run_pipeline()
        
        if success:
            print(f"\n✓ Pipeline completed successfully!")
            print(f"Results available in: {pipeline.output_dir}")
            sys.exit(0)
        else:
            print(f"\n✗ Pipeline completed with errors")
            print(f"Check logs in: {pipeline.output_dir}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Pipeline error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

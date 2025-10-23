#!/usr/bin/env python3
"""
Unified Experiment Analyzer - Main Entry Point
Comprehensive analysis system for ID3 experimental results
"""

import argparse
import json
import logging
from pathlib import Path
from datetime import datetime
import sys
import os

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from id3.experiments.analysis.result_parser import ExperimentResultParser
from id3.experiments.analysis.visualization import ExperimentVisualizer
from id3.experiments.analysis.performance_tables import PerformanceTableGenerator

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class UnifiedExperimentAnalyzer:
    """Main analyzer for unified experiments"""
    
    def __init__(self, input_dir: str, output_dir: str = None):
        """
        Initialize analyzer
        
        Args:
            input_dir: Directory containing experiment result JSON files
            output_dir: Directory to save analysis outputs (default: input_dir/analysis)
        """
        self.input_dir = Path(input_dir)
        if not self.input_dir.exists():
            raise ValueError(f"Input directory does not exist: {input_dir}")
        
        # Set output directory
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = self.input_dir / 'analysis'
        
        # Create output subdirectories
        self.figures_dir = self.output_dir / 'figures'
        self.tables_dir = self.output_dir / 'tables'
        
        for dir_path in [self.output_dir, self.figures_dir, self.tables_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.parser = ExperimentResultParser(self.input_dir)
        self.visualizer = ExperimentVisualizer(self.figures_dir)
        self.table_gen = PerformanceTableGenerator(self.tables_dir)
        
        # Analysis results
        self.results = []
        self.analysis_summary = {}
        
    def run_analysis(self, plots_only: bool = False, tables_only: bool = False):
        """
        Run complete analysis pipeline
        
        Args:
            plots_only: Only generate plots
            tables_only: Only generate tables
        """
        logger.info("Starting unified experiment analysis...")
        
        # Load all results
        logger.info("Loading experiment results...")
        self.results = self.parser.load_all_results()
        logger.info(f"Loaded {len(self.results)} experiment results")
        
        # Generate analysis summary
        self.analysis_summary = {
            'timestamp': datetime.now().isoformat(),
            'input_dir': str(self.input_dir),
            'output_dir': str(self.output_dir),
            'total_experiments': len(self.results),
            'statistics': self.parser.get_summary_statistics(),
            'best_configurations': self.parser.find_best_configurations(),
            'cai_comparison': self.parser.compare_cai_effect(),
        }
        
        # Generate visualizations
        if not tables_only:
            logger.info("Generating visualizations...")
            self._generate_all_plots()
        
        # Generate tables
        if not plots_only:
            logger.info("Generating performance tables...")
            self._generate_all_tables()
        
        # Generate report
        logger.info("Generating analysis report...")
        self._generate_markdown_report()
        
        # Save analysis summary
        self._save_analysis_summary()
        
        logger.info(f"Analysis complete! Results saved to: {self.output_dir}")
        
    def _generate_all_plots(self):
        """Generate all visualization plots"""
        
        # 1. Convergence curves
        logger.info("  Creating convergence curves...")
        self.visualizer.plot_convergence_curves(
            self.results,
            save_path=self.figures_dir / 'convergence_curves.png'
        )
        
        # 2. Performance comparison
        logger.info("  Creating performance comparison...")
        self.visualizer.plot_performance_comparison(
            self.results,
            save_path=self.figures_dir / 'performance_comparison.png'
        )
        
        # 3. CAI effect analysis
        logger.info("  Creating CAI effect analysis...")
        self.visualizer.plot_cai_effect_analysis(
            self.results,
            save_path=self.figures_dir / 'cai_effect_analysis.png'
        )
        
        # 4. Protein comparison
        logger.info("  Creating protein comparison...")
        self.visualizer.plot_protein_comparison(
            self.results,
            save_path=self.figures_dir / 'protein_comparison.png'
        )
        
        # 5. Performance heatmap
        logger.info("  Creating performance heatmap...")
        matrix, proteins, combinations = self.parser.get_performance_matrix()
        self.visualizer.plot_performance_heatmap(
            matrix, proteins, combinations,
            save_path=self.figures_dir / 'performance_heatmap.png'
        )
        
        # 6. Constraint-variant analysis
        logger.info("  Creating constraint-variant analysis...")
        self.visualizer.plot_constraint_variant_analysis(
            self.results,
            save_path=self.figures_dir / 'constraint_variant_analysis.png'
        )
        
        logger.info(f"  Saved {6} plots to {self.figures_dir}")
        
    def _generate_all_tables(self):
        """Generate all performance tables"""
        
        logger.info("  Generating tables in multiple formats...")
        
        # Save all tables in CSV, LaTeX, and Markdown formats
        self.table_gen.save_tables(self.results, format='all')
        
        # Print summary to console
        self.table_gen.print_summary(self.results)
        
        logger.info(f"  Saved tables to {self.tables_dir}")
        
    def _generate_markdown_report(self):
        """Generate comprehensive Markdown report"""
        
        report_path = self.output_dir / 'report.md'
        
        with open(report_path, 'w') as f:
            # Header
            f.write("# Unified Experiment Analysis Report\n\n")
            f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"**Input Directory**: `{self.input_dir}`\n\n")
            f.write(f"**Total Experiments**: {len(self.results)}\n\n")
            
            # Summary Statistics
            f.write("## Summary Statistics\n\n")
            stats = self.analysis_summary['statistics']
            f.write(f"- **Total Experiments**: {stats.get('total_experiments', 0)}\n")
            f.write(f"- **Unique Proteins**: {stats.get('unique_proteins', 0)}\n")
            f.write(f"- **Unique Constraints**: {stats.get('unique_constraints', 0)}\n")
            f.write(f"- **Unique Variants**: {stats.get('unique_variants', 0)}\n")
            f.write(f"- **Average Improvement**: {stats.get('avg_improvement', 0):.4f} ± {stats.get('std_improvement', 0):.4f}\n")
            f.write(f"- **Average Final Accessibility**: {stats.get('avg_final_accessibility', 0):.4f} ± {stats.get('std_final_accessibility', 0):.4f}\n")
            
            if stats.get('cai_experiments'):
                f.write(f"\n### CAI Statistics\n")
                f.write(f"- **CAI-Enabled Experiments**: {stats.get('cai_experiments', 0)}\n")
                f.write(f"- **Average Final ECAI**: {stats.get('avg_final_ecai', 0):.4f} ± {stats.get('std_final_ecai', 0):.4f}\n")
            
            # Best Configurations
            f.write("\n## Best Configurations by Protein\n\n")
            f.write("| Protein | Constraint | Variant | Final Accessibility | Improvement | CAI Enabled |\n")
            f.write("|---------|------------|---------|-------------------|-------------|-------------|\n")
            
            best_configs = self.analysis_summary['best_configurations']
            for protein, config in sorted(best_configs.items(), 
                                         key=lambda x: x[1]['final_accessibility']):
                f.write(f"| {protein} | {config['constraint_type']} | {config['variant']} | "
                       f"{config['final_accessibility']:.4f} | {config['improvement']:.4f} | "
                       f"{'Yes' if config['cai_enabled'] else 'No'} |\n")
            
            # CAI Comparison
            f.write("\n## CAI Effect Analysis\n\n")
            cai_comp = self.analysis_summary['cai_comparison']
            if cai_comp.get('cai_enabled_count'):
                f.write(f"- **CAI-Enabled Count**: {cai_comp['cai_enabled_count']}\n")
                f.write(f"- **CAI-Disabled Count**: {cai_comp.get('cai_disabled_count', 0)}\n")
                
                if 'cai_avg_improvement' in cai_comp:
                    f.write(f"- **CAI Average Improvement**: {cai_comp['cai_avg_improvement']:.4f}\n")
                if 'no_cai_avg_improvement' in cai_comp:
                    f.write(f"- **No-CAI Average Improvement**: {cai_comp['no_cai_avg_improvement']:.4f}\n")
                if 'cai_avg_final' in cai_comp:
                    f.write(f"- **CAI Average Final**: {cai_comp['cai_avg_final']:.4f}\n")
                if 'no_cai_avg_final' in cai_comp:
                    f.write(f"- **No-CAI Average Final**: {cai_comp['no_cai_avg_final']:.4f}\n")
            
            # Figures Section
            f.write("\n## Visualization Results\n\n")
            
            figure_descriptions = [
                ("convergence_curves.png", "Convergence Curves", 
                 "Shows the optimization progress over iterations for different configurations."),
                ("performance_comparison.png", "Performance Comparison", 
                 "Compares average performance across constraint types and variants."),
                ("cai_effect_analysis.png", "CAI Effect Analysis", 
                 "Analyzes the impact of CAI optimization on results."),
                ("protein_comparison.png", "Protein Comparison", 
                 "Compares optimization performance across different proteins."),
                ("performance_heatmap.png", "Performance Heatmap", 
                 "Visual matrix of all experiment results."),
                ("constraint_variant_analysis.png", "Constraint-Variant Analysis", 
                 "Detailed analysis of constraint types and variant interactions."),
            ]
            
            for filename, title, description in figure_descriptions:
                f.write(f"### {title}\n\n")
                f.write(f"{description}\n\n")
                f.write(f"![{title}](figures/{filename})\n\n")
            
            # Tables Section
            f.write("## Performance Tables\n\n")
            f.write("The following tables have been generated in CSV, LaTeX, and Markdown formats:\n\n")
            f.write("1. **comprehensive_performance**: All experiments with key metrics\n")
            f.write("2. **best_configurations**: Best configuration for each protein\n")
            f.write("3. **cai_comparison**: Comparison of CAI vs non-CAI results\n")
            f.write("4. **statistical_summary**: Statistical summary by category\n")
            f.write("5. **variant_matrix**: Performance matrix across all variants\n\n")
            f.write("Tables are available in the `tables/` directory.\n\n")
            
            # Footer
            f.write("---\n\n")
            f.write("*This report was automatically generated by the Unified Experiment Analyzer.*\n")
        
        logger.info(f"  Saved report to {report_path}")
    
    def _save_analysis_summary(self):
        """Save analysis summary as JSON"""
        
        summary_path = self.output_dir / 'summary.json'
        
        # Convert non-serializable objects
        def convert_for_json(obj):
            if hasattr(obj, 'tolist'):  # numpy arrays
                return obj.tolist()
            elif hasattr(obj, '__dict__'):  # custom objects
                return obj.__dict__
            else:
                return str(obj)
        
        with open(summary_path, 'w') as f:
            json.dump(self.analysis_summary, f, indent=2, default=convert_for_json)
        
        logger.info(f"  Saved analysis summary to {summary_path}")


def main():
    """Main entry point for command-line usage"""
    
    parser = argparse.ArgumentParser(
        description='Analyze unified experiment results from ID3 framework'
    )
    
    parser.add_argument(
        '--input-dir', '-i',
        type=str,
        required=True,
        help='Directory containing experiment result JSON files'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        default=None,
        help='Directory to save analysis outputs (default: input_dir/analysis)'
    )
    
    parser.add_argument(
        '--plots-only',
        action='store_true',
        help='Only generate visualization plots'
    )
    
    parser.add_argument(
        '--tables-only',
        action='store_true',
        help='Only generate performance tables'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Create analyzer and run analysis
        analyzer = UnifiedExperimentAnalyzer(
            input_dir=args.input_dir,
            output_dir=args.output_dir
        )
        
        analyzer.run_analysis(
            plots_only=args.plots_only,
            tables_only=args.tables_only
        )
        
        print(f"\n✅ Analysis complete! Results saved to: {analyzer.output_dir}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
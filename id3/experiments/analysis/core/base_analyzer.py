"""
Base Analyzer Class

Provides common functionality for all analysis modules.
"""

import json
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Any, Optional
import pandas as pd

logger = logging.getLogger(__name__)


class BaseAnalyzer(ABC):
    """
    Abstract base class for all analyzers.
    
    Provides common functionality for loading experiment results,
    filtering data, and saving outputs.
    """
    
    def __init__(self, config):
        """
        Initialize base analyzer.
        
        Args:
            config: AnalysisConfig instance
        """
        self.config = config
        self.input_dir = Path(config.input_dir)
        self.output_dir = config.get_output_dir() / self.get_analyzer_name()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results = []
        self.filtered_results = []
        self.analysis_results = {}
        
    @abstractmethod
    def get_analyzer_name(self) -> str:
        """Get the name of this analyzer for output directory."""
        pass
    
    @abstractmethod
    def _perform_analysis(self) -> Dict[str, Any]:
        """
        Perform the actual analysis.
        Must be implemented by subclasses.
        
        Returns:
            Dictionary containing analysis results
        """
        pass
    
    @abstractmethod
    def _generate_visualizations(self) -> Dict[str, Path]:
        """
        Generate visualizations.
        Must be implemented by subclasses.
        
        Returns:
            Dictionary mapping visualization names to file paths
        """
        pass
    
    def analyze(self) -> Dict[str, Any]:
        """
        Main analysis pipeline.
        
        Returns:
            Dictionary containing all analysis results and output files
        """
        # Load experiment results
        self._load_results()
        
        # Apply filters
        self._apply_filters()
        
        # Perform analysis
        self.analysis_results = self._perform_analysis()
        
        # Generate visualizations if requested
        visualization_files = {}
        if self.config.visualize:
            visualization_files = self._generate_visualizations()
        
        # Save results
        output_files = self._save_results()
        
        # Combine all outputs
        return {
            'summary': self.analysis_results.get('summary', {}),
            'detailed_results': self.analysis_results,
            'output_files': {**output_files, **visualization_files},
            'num_experiments_analyzed': len(self.filtered_results),
        }
    
    def _load_results(self) -> None:
        """Load all experiment result files from input directory."""
        result_files = list(self.input_dir.glob("*.json"))
        
        # Exclude config and summary files
        result_files = [
            f for f in result_files 
            if f.name not in ['config.json', 'summary.json']
        ]
        
        if not result_files:
            raise FileNotFoundError(f"No experiment result files found in {self.input_dir}")
        
        logger.info(f"Loading {len(result_files)} experiment results...")
        
        for file_path in result_files:
            try:
                with open(file_path) as f:
                    data = json.load(f)
                    # Add filename for reference
                    data['_filename'] = file_path.name
                    self.results.append(data)
            except Exception as e:
                logger.warning(f"Failed to load {file_path}: {e}")
        
        logger.info(f"Successfully loaded {len(self.results)} results")
    
    def _apply_filters(self) -> None:
        """Apply configured filters to results."""
        self.filtered_results = self.results.copy()
        
        # Filter by protein
        if self.config.filter_proteins:
            self.filtered_results = [
                r for r in self.filtered_results
                if r.get('protein_name') in self.config.filter_proteins
            ]
            logger.info(f"Filtered to {len(self.filtered_results)} results by protein")
        
        # Filter by constraint type
        if self.config.filter_constraints:
            self.filtered_results = [
                r for r in self.filtered_results
                if r.get('constraint_type') in self.config.filter_constraints
            ]
            logger.info(f"Filtered to {len(self.filtered_results)} results by constraint")
        
        # Filter by variant
        if self.config.filter_variants:
            self.filtered_results = [
                r for r in self.filtered_results
                if r.get('variant') in self.config.filter_variants
            ]
            logger.info(f"Filtered to {len(self.filtered_results)} results by variant")
        
        if not self.filtered_results:
            raise ValueError("No results remain after applying filters")
        
        logger.info(f"Analyzing {len(self.filtered_results)} experiments after filtering")
    
    def _save_results(self) -> Dict[str, Path]:
        """
        Save analysis results to files.
        
        Returns:
            Dictionary mapping output types to file paths
        """
        output_files = {}
        
        # Save JSON if requested
        if self.config.output_format in ['json', 'both']:
            json_file = self.output_dir / f"{self.get_analyzer_name()}_results.json"
            with open(json_file, 'w') as f:
                json.dump(self.analysis_results, f, indent=2, default=str)
            output_files['json'] = json_file
            logger.info(f"Saved JSON results to {json_file}")
        
        # Save CSV if requested and applicable
        if self.config.output_format in ['csv', 'both']:
            if 'dataframe' in self.analysis_results:
                csv_file = self.output_dir / f"{self.get_analyzer_name()}_data.csv"
                df = self.analysis_results['dataframe']
                if isinstance(df, pd.DataFrame):
                    df.to_csv(csv_file, index=False)
                    output_files['csv'] = csv_file
                    logger.info(f"Saved CSV data to {csv_file}")
        
        # Save summary text
        summary_file = self.output_dir / f"{self.get_analyzer_name()}_summary.txt"
        self._write_summary(summary_file)
        output_files['summary'] = summary_file
        
        return output_files
    
    def _write_summary(self, filepath: Path) -> None:
        """
        Write human-readable summary to text file.
        
        Args:
            filepath: Path to output file
        """
        with open(filepath, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"{self.get_analyzer_name().upper()} ANALYSIS SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Input Directory: {self.input_dir}\n")
            f.write(f"Total Experiments Analyzed: {len(self.filtered_results)}\n")
            f.write(f"Analysis Time: {pd.Timestamp.now()}\n\n")
            
            # Subclasses should override this to add specific summary content
            self._write_specific_summary(f)
    
    def _write_specific_summary(self, file_handle) -> None:
        """
        Write analyzer-specific summary content.
        Should be overridden by subclasses.
        
        Args:
            file_handle: Open file handle for writing
        """
        pass
    
    def get_results_dataframe(self) -> pd.DataFrame:
        """
        Convert filtered results to pandas DataFrame.
        
        Returns:
            DataFrame with experiment results
        """
        if not self.filtered_results:
            return pd.DataFrame()
        
        # Extract relevant fields
        data = []
        for result in self.filtered_results:
            data.append({
                'protein': result.get('protein_name'),
                'constraint_type': result.get('constraint_type'),
                'variant': result.get('variant'),
                'seed': result.get('seed'),
                'best_accessibility': result.get('best_accessibility'),
                'final_accessibility': result.get('final_accessibility'),
                'initial_accessibility': result.get('initial_accessibility'),
                'improvement': result.get('improvement'),
                'iterations': result.get('iterations'),
                'optimization_time': result.get('optimization_time'),
                'amino_acids_match': result.get('amino_acids_match'),
                'filename': result.get('_filename'),
            })
        
        return pd.DataFrame(data)
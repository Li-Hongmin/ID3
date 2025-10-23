#!/usr/bin/env python3
"""
Results Management Module

Responsible for saving, loading, and analyzing experiment results
"""

import json
import time
from pathlib import Path
from typing import Dict, Any, Optional
import torch
import numpy as np


class ResultsManager:
    """Experiment results manager"""

    def __init__(self, base_dir: str = "experiments/results"):
        """
        Initialize results manager

        Args:
            base_dir: Base directory for saving results
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
    
    def create_experiment_dir(self, experiment_name: str) -> Path:
        """
        Create experiment directory

        Args:
            experiment_name: Experiment name

        Returns:
            Experiment directory path
        """
        exp_dir = self.base_dir / experiment_name
        exp_dir.mkdir(parents=True, exist_ok=True)
        return exp_dir
    
    def save_results(
        self,
        results: Dict[str, Any],
        experiment_name: str,
        protein_id: str,
        suffix: Optional[str] = None
    ) -> Path:
        """
        Save experiment results

        Args:
            results: Result data
            experiment_name: Experiment name
            protein_id: Protein ID
            suffix: Filename suffix

        Returns:
            Path to saved file
        """
        # Create experiment directory
        exp_dir = self.create_experiment_dir(experiment_name)

        # Build filename
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        if suffix:
            filename = f"{protein_id}_{suffix}_{timestamp}.json"
        else:
            filename = f"{protein_id}_{timestamp}.json"
        
        filepath = exp_dir / filename

        # Convert data types
        results_json = self._prepare_for_json(results)

        # Save
        with open(filepath, 'w') as f:
            json.dump(results_json, f, indent=2)

        print(f"Results saved: {filepath}")
        return filepath

    def load_results(self, filepath: str) -> Dict[str, Any]:
        """
        Load experiment results

        Args:
            filepath: Result file path

        Returns:
            Result data
        """
        with open(filepath, 'r') as f:
            return json.load(f)
    
    def _prepare_for_json(self, obj: Any) -> Any:
        """
        Prepare data for JSON serialization

        Args:
            obj: Object to serialize

        Returns:
            Serializable object
        """
        if isinstance(obj, dict):
            return {k: self._prepare_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._prepare_for_json(v) for v in obj]
        elif isinstance(obj, torch.Tensor):
            return obj.cpu().numpy().tolist()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif hasattr(obj, 'item'):
            return obj.item()
        else:
            return obj
    
    def compare_results(
        self,
        result1: Dict[str, Any],
        result2: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Compare two experiment results

        Args:
            result1: First result
            result2: Second result

        Returns:
            Comparison result
        """
        comparison = {}

        # Compare key metrics
        metrics = [
            'initial_accessibility',
            'final_accessibility',
            'improvement',
            'optimization_time',
            'amino_acids_correct'
        ]
        
        for metric in metrics:
            if metric in result1 and metric in result2:
                val1 = result1[metric]
                val2 = result2[metric]
                
                if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                    comparison[metric] = {
                        'result1': val1,
                        'result2': val2,
                        'difference': val1 - val2,
                        'relative_diff': (val1 - val2) / val2 if val2 != 0 else None
                    }
                else:
                    comparison[metric] = {
                        'result1': val1,
                        'result2': val2,
                        'equal': val1 == val2
                    }
        
        return comparison
"""
Progress tracking for experiments.

This module handles experiment progress tracking, saving, and resumption
capabilities for long-running experiment batches.
"""

import json
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Optional
from tqdm import tqdm

logger = logging.getLogger(__name__)


class ExperimentProgressTracker:
    """Tracks and manages experiment progress."""

    def __init__(self, output_dir: Path):
        """
        Initialize progress tracker.

        Args:
            output_dir: Directory for saving progress files
        """
        self.output_dir = output_dir
        self.progress_file = output_dir / 'progress.json'
        self.progress_data = None
        self.pbar = None
        self.completed_experiments = self._load_existing_progress()

    def _load_existing_progress(self) -> Set[str]:
        """
        Load existing progress from progress.json.

        Returns:
            Set of completed experiment identifiers
        """
        completed = set()
        if self.progress_file.exists():
            try:
                with open(self.progress_file, 'r') as f:
                    data = json.load(f)
                    completed_list = data.get('completed_experiments', [])
                    for exp in completed_list:
                        if isinstance(exp, dict):
                            # Create unique identifier for the experiment
                            exp_id = self._create_experiment_id(exp)
                            completed.add(exp_id)
                    logger.info(f"Loaded {len(completed)} completed experiments from progress file")
            except Exception as e:
                logger.warning(f"Failed to load progress file: {e}")
        return completed

    def _create_experiment_id(self, exp: Dict) -> str:
        """
        Create unique identifier for an experiment.

        Args:
            exp: Experiment configuration dictionary

        Returns:
            Unique string identifier
        """
        return f"{exp.get('protein_id', 'unknown')}_{exp.get('constraint_type', 'unknown')}_{exp.get('variant', 'unknown')}"

    def initialize(self, total_experiments: int, already_completed: int = 0):
        """
        Initialize progress tracking for a batch of experiments.

        Args:
            total_experiments: Total number of experiments to run
            already_completed: Number of experiments already completed
        """
        self.progress_data = {
            'total_experiments': total_experiments,
            'completed_experiments': [],
            'failed_experiments': [],
            'start_time': datetime.now().isoformat(),
            'status': 'running',
            'already_completed': already_completed
        }

        # Initialize progress bar
        self.pbar = tqdm(
            total=total_experiments,
            initial=already_completed,
            desc="Running experiments",
            unit="exp"
        )

        # Save initial progress
        self._save_progress()

    def update(self, status: str, experiment: Dict, index: int, result: Optional[Dict] = None):
        """
        Update progress for a single experiment.

        Args:
            status: Status of the experiment ('completed', 'failed', 'skipped')
            experiment: Experiment configuration
            index: Index of the experiment
            result: Result dictionary (if completed)
        """
        if self.progress_data is None:
            return

        if status == 'completed' and result:
            # Add to completed list
            exp_record = {
                **experiment,
                'index': index,
                'completion_time': datetime.now().isoformat()
            }
            self.progress_data['completed_experiments'].append(exp_record)

            # Update completed set
            exp_id = self._create_experiment_id(experiment)
            self.completed_experiments.add(exp_id)

        elif status == 'failed':
            # Add to failed list
            exp_record = {
                **experiment,
                'index': index,
                'failure_time': datetime.now().isoformat()
            }
            self.progress_data['failed_experiments'].append(exp_record)

        # Update progress bar
        if self.pbar and status in ['completed', 'failed']:
            self.pbar.update(1)

        # Save progress periodically (every 10 experiments)
        if len(self.progress_data['completed_experiments']) % 10 == 0:
            self._save_progress()

    def finalize(self):
        """Finalize progress tracking and save final state."""
        if self.progress_data:
            self.progress_data['end_time'] = datetime.now().isoformat()
            self.progress_data['status'] = 'completed'
            self._save_progress()

        if self.pbar:
            self.pbar.close()

    def _save_progress(self):
        """Save current progress to file."""
        if self.progress_data:
            try:
                with open(self.progress_file, 'w') as f:
                    json.dump(self.progress_data, f, indent=2)
            except Exception as e:
                logger.error(f"Failed to save progress: {e}")

    def is_experiment_completed(self, experiment: Dict) -> bool:
        """
        Check if an experiment has already been completed.

        Args:
            experiment: Experiment configuration

        Returns:
            True if experiment is already completed
        """
        exp_id = self._create_experiment_id(experiment)
        return exp_id in self.completed_experiments

    def get_stats(self) -> Dict:
        """
        Get current progress statistics.

        Returns:
            Dictionary with progress statistics
        """
        if self.progress_data:
            return {
                'total': self.progress_data['total_experiments'],
                'completed': len(self.progress_data['completed_experiments']),
                'failed': len(self.progress_data['failed_experiments']),
                'already_completed': self.progress_data.get('already_completed', 0)
            }
        return {'total': 0, 'completed': 0, 'failed': 0, 'already_completed': 0}
"""
Base class for experiment runners with common functionality.

This module provides the foundation for experiment execution,
separating concerns for better maintainability.
"""

import torch
import logging
from typing import Dict, Optional
from pathlib import Path
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
from id3.experiments.utils.data_loader import ProteinDataLoader

logger = logging.getLogger(__name__)


class ExperimentBase:
    """Base class with common experiment functionality."""

    def __init__(self, config: Dict):
        """
        Initialize base experiment runner.

        Args:
            config: Experiment configuration dictionary
        """
        self.config = config
        self.device = torch.device(config.get('device', 'cuda'))
        self.data_loader = ProteinDataLoader()
        self.deepraccess = DeepRaccessID3Wrapper()

        # Setup output directory
        if config.get('output_dir'):
            self.output_dir = Path(config['output_dir'])
        else:
            self.output_dir = self._get_default_output_dir()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _get_default_output_dir(self) -> Path:
        """Get default output directory based on config."""
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return Path(f"experiments/results/{timestamp}")

    def cleanup(self):
        """Clean up resources."""
        if hasattr(self, 'deepraccess'):
            del self.deepraccess
        torch.cuda.empty_cache()
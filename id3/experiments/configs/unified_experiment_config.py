#!/usr/bin/env python3
"""
Unified Experiment Configuration Module

Centralized configuration management using constants from id3.utils.constants.
No hardcoding allowed - all values come from constants or command-line arguments.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
import argparse
import logging
from pathlib import Path
from datetime import datetime

logger = logging.getLogger(__name__)

# Import all constants from the centralized location
from id3.utils.constants import (
    NUCLEOTIDES,
    NUCLEOTIDE_MAP,
    amino_acids_to_codons,
    num_codons,
    num_amino_acids
)


@dataclass
class UnifiedExperimentConfig:
    """
    Configuration for unified experiments.

    This class encapsulates all experiment parameters,
    using constants from id3.utils.constants instead of hardcoding.
    """

    # Protein selection
    proteins: List[str] = field(default_factory=lambda: ['O15263'])

    # Constraint configuration
    constraints: List[str] = field(default_factory=lambda: ['lagrangian'])
    variants: List[str] = field(default_factory=lambda: ['11'])

    # Optimization parameters
    iterations: int = 1000
    learning_rate: float = 0.05
    batch_size: int = 1
    seeds: int = 1
    base_seed: int = 42  # Starting seed value
    device: str = 'cuda'

    # CAI parameters
    enable_cai: bool = False
    cai_target: float = 0.8
    lambda_cai: float = 1.0
    species: str = 'ecoli_bl21de3'  # Can be overridden, not hardcoded in logic
    disable_constraint_penalty: bool = False  # Disable constraint penalty (CAI no penalty mode)
    
    # Adaptive lambda_cai parameters
    adaptive_lambda_cai: bool = False
    lambda_cai_lr: float = 0.1
    lambda_cai_max: float = 2.0
    lambda_cai_min: float = 0.01
    cai_tolerance: float = 0.05
    smoothing_factor: float = 0.9

    # Execution parameters
    # Removed parallel parameter, always use serial execution for optimal performance
    verbose: bool = False
    mixed_precision: bool = False  # Mixed precision training (fully enabled or fully disabled)
    gradient_clip: float = 1.0  # Gradient clipping threshold
    disable_inner_tqdm: bool = False  # Disable inner tqdm progress bar to reduce logging

    # Output parameters
    output_dir: Optional[str] = None
    save_trajectories: bool = True

    # Constants from id3.utils.constants (read-only)
    nucleotides: List[str] = field(default_factory=lambda: NUCLEOTIDES, init=False)
    nucleotide_map: Dict[str, int] = field(default_factory=lambda: NUCLEOTIDE_MAP.copy(), init=False)
    genetic_code: Dict[str, List[str]] = field(default_factory=lambda: amino_acids_to_codons.copy(), init=False)

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'UnifiedExperimentConfig':
        """
        Create configuration from command-line arguments.

        Preserves default configuration values and only overrides when arguments are provided.

        Args:
            args: Parsed command-line arguments

        Returns:
            UnifiedExperimentConfig instance
        """
        # Create base config with defaults
        config = cls()

        # Override only if arguments are provided
        if hasattr(args, 'proteins') and args.proteins:
            config.proteins = args.proteins.split(',')

        if hasattr(args, 'constraints') and args.constraints:
            config.constraints = args.constraints.split(',')

        if hasattr(args, 'variants') and args.variants:
            config.variants = args.variants.split(',')

        if hasattr(args, 'iterations') and args.iterations is not None:
            config.iterations = args.iterations

        if hasattr(args, 'learning_rate') and args.learning_rate is not None:
            config.learning_rate = args.learning_rate

        if hasattr(args, 'batch_size') and args.batch_size is not None:
            config.batch_size = args.batch_size

        if hasattr(args, 'seeds') and args.seeds is not None:
            config.seeds = args.seeds

        if hasattr(args, 'base_seed') and args.base_seed is not None:
            config.base_seed = args.base_seed

        if hasattr(args, 'device') and args.device:
            config.device = args.device

        if hasattr(args, 'enable_cai'):
            config.enable_cai = args.enable_cai

        if hasattr(args, 'cai_target') and args.cai_target is not None:
            config.cai_target = args.cai_target

        if hasattr(args, 'lambda_cai') and args.lambda_cai is not None:
            config.lambda_cai = args.lambda_cai

        if hasattr(args, 'species') and args.species:
            config.species = args.species
        
        if hasattr(args, 'disable_constraint_penalty'):
            config.disable_constraint_penalty = args.disable_constraint_penalty

        # Adaptive lambda_cai parameters
        if hasattr(args, 'adaptive_lambda_cai'):
            config.adaptive_lambda_cai = args.adaptive_lambda_cai

        if hasattr(args, 'lambda_cai_lr') and args.lambda_cai_lr is not None:
            config.lambda_cai_lr = args.lambda_cai_lr

        if hasattr(args, 'lambda_cai_max') and args.lambda_cai_max is not None:
            config.lambda_cai_max = args.lambda_cai_max

        if hasattr(args, 'lambda_cai_min') and args.lambda_cai_min is not None:
            config.lambda_cai_min = args.lambda_cai_min

        if hasattr(args, 'cai_tolerance') and args.cai_tolerance is not None:
            config.cai_tolerance = args.cai_tolerance

        if hasattr(args, 'smoothing_factor') and args.smoothing_factor is not None:
            config.smoothing_factor = args.smoothing_factor

        # Removed parallel parameter handling, always use serial execution

        if hasattr(args, 'verbose'):
            config.verbose = args.verbose

        if hasattr(args, 'disable_inner_tqdm'):
            config.disable_inner_tqdm = args.disable_inner_tqdm

        if hasattr(args, 'output_dir') and args.output_dir:
            config.output_dir = args.output_dir

        if hasattr(args, 'save_trajectories') and args.save_trajectories is not None:
            config.save_trajectories = args.save_trajectories

        return config

    def _get_mode_string(self) -> str:
        """Get mode description string"""
        if not self.enable_cai:
            return 'unified_accessibility_optimization'
        elif self.disable_constraint_penalty:
            return 'unified_cai_no_penalty'
        else:
            return 'unified_cai_with_penalty'

    def _get_theory_string(self) -> str:
        """Get theory formula string"""
        if not self.enable_cai:
            return 'L_Access + constraint_penalties'
        elif self.disable_constraint_penalty:
            return 'L_Access + λ_CAI * L_CAI (no constraint penalty)'
        else:
            return 'L_Access + λ_CAI * L_CAI + constraint_penalties'

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to dictionary for serialization.

        Returns:
            Dictionary representation of configuration
        """
        return {
            'proteins': self.proteins,
            'constraints': self.constraints,
            'variants': self.variants,
            'iterations': self.iterations,
            'learning_rate': self.learning_rate,
            'batch_size': self.batch_size,
            'seeds': self.seeds,
            'base_seed': self.base_seed,
            'device': self.device,
            'enable_cai': self.enable_cai,
            'cai_target': self.cai_target if self.enable_cai else None,
            'lambda_cai': self.lambda_cai if self.enable_cai else None,
            'species': self.species if self.enable_cai else None,
            'disable_constraint_penalty': self.disable_constraint_penalty if self.enable_cai else False,
            'verbose': self.verbose,
            'disable_inner_tqdm': self.disable_inner_tqdm,
            'output_dir': self.output_dir,
            'save_trajectories': self.save_trajectories,
            'mode': self._get_mode_string(),
            'theory': self._get_theory_string()
        }

    def get_output_dir(self) -> Path:
        """
        Get or create output directory.

        Returns:
            Path to output directory
        """
        if self.output_dir:
            output_dir = Path(self.output_dir)
        else:
            # If timestamp hasn't been set yet, create a fixed timestamp
            if not hasattr(self, '_output_timestamp'):
                self._output_timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

            mode_str = 'unified_cai' if self.enable_cai else 'unified_access'
            output_dir = Path(f'results/{self._output_timestamp}_{mode_str}_experiments')

        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def generate_experiments(self) -> List[Dict[str, Any]]:
        """
        Generate list of experiment specifications based on configuration.

        Returns:
            List of experiment dictionaries
        """
        experiments = []
        for protein in self.proteins:
            for constraint in self.constraints:
                for variant in self.variants:
                    for seed_offset in range(self.seeds):
                        experiments.append({
                            'protein_name': protein,
                            'constraint_type': constraint,
                            'variant': variant,
                            'seed': self.base_seed + seed_offset
                        })
        return experiments

    def validate(self) -> None:
        """
        Validate configuration parameters.

        Raises:
            ValueError: If configuration is invalid
        """
        # Validate constraint types
        valid_constraints = ['lagrangian', 'ams', 'cpc']
        for constraint in self.constraints:
            if constraint.lower() not in valid_constraints:
                raise ValueError(f"Invalid constraint type: {constraint}. "
                               f"Must be one of {valid_constraints}")

        # Validate variants
        for variant in self.variants:
            if len(variant) != 2 or not all(c in '01' for c in variant):
                raise ValueError(f"Invalid variant: {variant}. "
                               "Must be 2 digits of 0 or 1 (e.g., '00', '01', '10', '11')")

        # Validate CAI parameters
        if self.enable_cai:
            if not 0.0 <= self.cai_target <= 1.0:
                raise ValueError(f"CAI target must be between 0 and 1, got {self.cai_target}")
            if self.lambda_cai < 0:
                raise ValueError(f"Lambda CAI must be non-negative, got {self.lambda_cai}")

        # Validate optimization parameters
        if self.iterations <= 0:
            raise ValueError(f"Iterations must be positive, got {self.iterations}")
        if self.learning_rate <= 0:
            raise ValueError(f"Learning rate must be positive, got {self.learning_rate}")
        if self.seeds <= 0:
            raise ValueError(f"Seeds must be positive, got {self.seeds}")

    def print_summary(self) -> None:
        """Print configuration summary."""
        logger.info("🚀 Unified ID3-DeepRaccess Experiment Configuration")
        logger.info("=" * 60)

        # Display mode information
        if not self.enable_cai:
            logger.info(f"Mode: Accessibility-only Optimization")
            logger.info(f"Theory: L_total = L_Access (+ constraint penalties)")
        elif self.disable_constraint_penalty:
            logger.info(f"Mode: CAI No Penalty Optimization")
            logger.info(f"Theory: L_total = L_Access + λ_CAI * L_CAI (no constraint penalty)")
            logger.info(f"CAI Target: {self.cai_target}")
            logger.info(f"λ_CAI Coefficient: {self.lambda_cai}")
            logger.info(f"Species: {self.species}")
        else:
            logger.info(f"Mode: CAI With Penalty Optimization")
            logger.info(f"Theory: L_total = L_Access + λ_CAI * L_CAI + constraint_penalties")
            logger.info(f"CAI Target: {self.cai_target}")
            logger.info(f"λ_CAI Coefficient: {self.lambda_cai}")
            logger.info(f"Species: {self.species}")
        logger.info(f"Proteins: {len(self.proteins)} - {self.proteins[:3]}{'...' if len(self.proteins) > 3 else ''}")
        logger.info(f"Constraints: {self.constraints}")
        logger.info(f"Variants: {self.variants}")
        logger.info(f"Iterations: {self.iterations}")
        logger.info(f"Learning Rate: {self.learning_rate}")
        logger.info(f"Seeds: {self.seeds} (starting from {self.base_seed})")
        logger.info(f"Device: {self.device}")
        logger.info(f"Output: {self.get_output_dir()}")
        logger.info("")


# Preset configurations for common experiments
class ExperimentPresets:
    """Predefined experiment configurations for common use cases."""

    @staticmethod
    def quick_test() -> UnifiedExperimentConfig:
        """Quick test configuration for debugging - O15263 with 12 variants (no CAI)."""
        return UnifiedExperimentConfig(
            proteins=['O15263'],  # Use O15263 protein for quick test
            constraints=['lagrangian', 'ams', 'cpc'],  # 3 constraints
            variants=['00', '01', '10', '11'],  # 4 variants
            iterations=2,  # 2 iterations for quick test
            seeds=1,
            enable_cai=False,  # Do not enable CAI by default
            verbose=True
        )

    @staticmethod
    def quick_test_cai_penalty() -> UnifiedExperimentConfig:
        """Quick test configuration for debugging - O15263 with 12 variants (with CAI and constraint penalty)."""
        return UnifiedExperimentConfig(
            proteins=['O15263'],  # Use O15263 protein for quick test
            constraints=['lagrangian', 'ams', 'cpc'],  # 3 constraints
            variants=['00', '01', '10', '11'],  # 4 variants
            iterations=2,  # 2 iterations for quick test
            seeds=1,
            enable_cai=True,  # Enable CAI optimization
            cai_target=0.8,
            lambda_cai=1.0,
            species='ecoli_bl21de3',
            disable_constraint_penalty=False,  # Include constraint penalty
            verbose=True
        )
    
    @staticmethod
    def quick_test_cai_no_penalty() -> UnifiedExperimentConfig:
        """Quick test configuration for debugging - O15263 with 12 variants (CAI no penalty)."""
        return UnifiedExperimentConfig(
            proteins=['O15263'],  # Use O15263 protein for quick test
            constraints=['lagrangian', 'ams', 'cpc'],  # 3 constraints
            variants=['00', '01', '10', '11'],  # 4 variants
            iterations=2,  # 2 iterations for quick test
            seeds=1,
            enable_cai=True,  # Enable CAI optimization
            cai_target=0.8,
            lambda_cai=1.0,
            species='ecoli_bl21de3',
            disable_constraint_penalty=True,  # Disable constraint penalty
            verbose=True
        )

    @staticmethod
    def quick_test_both() -> List[UnifiedExperimentConfig]:
        """Quick test configuration for both CAI and non-CAI optimization."""
        return [
            ExperimentPresets.quick_test(),      # Configuration without CAI
            ExperimentPresets.quick_test_cai_penalty()   # Configuration with CAI (with constraint penalty)
        ]
    
    @staticmethod
    def quick_test_cai_comparison() -> List[UnifiedExperimentConfig]:
        """Quick test configuration comparing CAI with and without constraint penalty."""
        return [
            ExperimentPresets.quick_test_cai_penalty(),   # CAI with constraint penalty
            ExperimentPresets.quick_test_cai_no_penalty() # CAI without constraint penalty
        ]

    @staticmethod
    def full_12x12() -> UnifiedExperimentConfig:
        """Full 12x12 experiment matrix (without CAI)."""
        # Use actually available proteins
        return UnifiedExperimentConfig(
            proteins=['O15263', 'P00004', 'P01308', 'P01825',
                     'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',
                     'P31417', 'P42212', 'P61626', 'P99999'],
            constraints=['lagrangian', 'ams', 'cpc'],
            variants=['00', '01', '10', '11'],
            iterations=1000,
            seeds=1,
            enable_cai=False  # Pure accessibility optimization
        )

    @staticmethod
    def full_12x12_cai_penalty() -> UnifiedExperimentConfig:
        """Full 12x12 experiment matrix with CAI optimization (with constraint penalty)."""
        # Use actually available proteins
        return UnifiedExperimentConfig(
            proteins=['O15263', 'P00004', 'P01308', 'P01825',
                     'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',
                     'P31417', 'P42212', 'P61626', 'P99999'],
            constraints=['lagrangian', 'ams', 'cpc'],
            variants=['00', '01', '10', '11'],
            iterations=1000,
            seeds=1,
            enable_cai=True,  # Enable CAI optimization
            cai_target=0.8,
            lambda_cai=1.0,
            species='ecoli_bl21de3',
            disable_constraint_penalty=False  # Include constraint penalty
        )
    
    @staticmethod
    def full_12x12_cai_no_penalty() -> UnifiedExperimentConfig:
        """Full 12x12 experiment matrix with CAI optimization (no constraint penalty)."""
        # Use actually available proteins
        return UnifiedExperimentConfig(
            proteins=['O15263', 'P00004', 'P01308', 'P01825',
                     'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',
                     'P31417', 'P42212', 'P61626', 'P99999'],
            constraints=['lagrangian', 'ams', 'cpc'],
            variants=['00', '01', '10', '11'],
            iterations=1000,
            seeds=1,
            enable_cai=True,  # Enable CAI optimization
            cai_target=0.8,
            lambda_cai=1.0,
            species='ecoli_bl21de3',
            disable_constraint_penalty=True  # Disable constraint penalty
        )

    @staticmethod
    def adaptive_lambda_cai_test() -> UnifiedExperimentConfig:
        """Quick test with adaptive lambda_cai optimization."""
        return UnifiedExperimentConfig(
            proteins=['P42212'],  # Best protein
            constraints=['lagrangian'],  # Best constraint type
            variants=['11'],  # Best variant
            iterations=100,
            seeds=1,
            enable_cai=True,
            cai_target=0.8,
            lambda_cai=0.1,  # Initial value
            species='ecoli_bl21de3',
            # Enable adaptive lambda_cai optimization
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.1,
            lambda_cai_max=2.0,
            lambda_cai_min=0.01,
            cai_tolerance=0.05,
            smoothing_factor=0.9,
            verbose=True
        )

    @staticmethod
    def full_12x12_both() -> List[UnifiedExperimentConfig]:
        """Full 12x12 experiment matrix for both CAI and non-CAI optimization."""
        return [
            ExperimentPresets.full_12x12_cai_penalty(),   # Full experiment with CAI enabled (with constraint penalty)
            ExperimentPresets.full_12x12(),      # Full experiment without CAI
        ]
    
    @staticmethod
    def full_12x12_cai_comparison() -> List[UnifiedExperimentConfig]:
        """Full 12x12 experiment matrix comparing CAI with and without constraint penalty."""
        return [
            ExperimentPresets.full_12x12_cai_penalty(),   # CAI with constraint penalty
            ExperimentPresets.full_12x12_cai_no_penalty() # CAI without constraint penalty
        ]

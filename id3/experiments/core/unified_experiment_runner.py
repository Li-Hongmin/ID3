#!/usr/bin/env python3
"""
Unified Experiment Runner Module

Implements the paper's theoretical approach: L_unified = L_Access + λ_CAI * L_CAI
This module handles the actual optimization logic, separated from the entry point.
"""

import time
import torch
import gc
import json
import os
import numpy as np
from typing import Dict, Optional, Tuple, List
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore', category=UserWarning)
import logging
from tqdm import tqdm

# Configure logger
logger = logging.getLogger(__name__)

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
from id3.experiments.utils.data_loader import ProteinDataLoader
from id3.utils.sequence_utils import rna_to_amino_acids
from id3.utils.constants import NUCLEOTIDE_MAP, NUCLEOTIDES



class UnifiedExperimentRunner:
    """
    Unified experiment runner that implements L_unified = L_Access + λ_CAI * L_CAI.

    This class encapsulates all optimization logic, keeping the entry point clean.
    """

    def __init__(self, config: Dict):
        """
        Initialize the experiment runner with configuration.

        Args:
            config: Experiment configuration dictionary
        """
        self.config = config
        # Fix: Get device from dictionary config instead of using getattr
        self.device = torch.device(config.get('device', 'cuda'))
        self.data_loader = ProteinDataLoader()
        self.deepraccess = DeepRaccessID3Wrapper()

        # Set output directory for incremental saving
        # Fix: Get output_dir from dictionary config
        if config.get('output_dir'):
            self.output_dir = Path(config['output_dir'])
        else:
            self.output_dir = self._get_output_dir()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load existing progress (if any)
        self.completed_experiments = self._load_existing_progress()

    def string_to_one_hot_tensor(self, sequence: str) -> torch.Tensor:
        """
        Convert RNA sequence string to one-hot tensor for DeepRaccess.
        Uses constants from id3.utils.constants instead of hardcoding.
        """
        seq_len = len(sequence)
        one_hot = torch.zeros(1, seq_len, 4, device=self.device)

        for i, nucleotide in enumerate(sequence):
            if nucleotide in NUCLEOTIDE_MAP:
                one_hot[0, i, NUCLEOTIDE_MAP[nucleotide]] = 1.0

        return one_hot

    def create_constraint(self,
                         constraint_type: str,
                         amino_acid_sequence: str,
                         variant: str) -> torch.nn.Module:
        """
        Create constraint mechanism with unified CAI integration.

        Args:
            constraint_type: Type of constraint ('lagrangian', 'ams', 'cpc')
            amino_acid_sequence: Target amino acid sequence
            variant: Variant specification (00, 01, 10, 11)

        Returns:
            Constraint instance with unified CAI support
        """
        # Parse variant: first digit = deterministic/gumbel, second digit = soft/ste
        if len(variant) != 2:
            raise ValueError(f"Variant must be 2 digits, got: {variant}")

        # Extract CAI parameters from config (Fix: use dictionary access)
        enable_cai = self.config.get('enable_cai', False)
        cai_target = self.config.get('cai_target', 0.8)
        lambda_cai = self.config.get('lambda_cai', 1.0)
        batch_size = self.config.get('batch_size', 1)
        
        # Extract adaptive lambda_cai parameters
        adaptive_lambda_cai = self.config.get('adaptive_lambda_cai', False)
        lambda_cai_lr = self.config.get('lambda_cai_lr', 0.1)
        lambda_cai_max = self.config.get('lambda_cai_max', 2.0)
        lambda_cai_min = self.config.get('lambda_cai_min', 0.01)
        cai_tolerance = self.config.get('cai_tolerance', 0.05)
        smoothing_factor = self.config.get('smoothing_factor', 0.9)

        if constraint_type.lower() == 'lagrangian':
            return LagrangianConstraint(
                amino_acid_sequence=amino_acid_sequence,
                batch_size=batch_size,
                device=self.device,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=lambda_cai,
                # Adaptive lambda_cai parameters
                adaptive_lambda_cai=adaptive_lambda_cai,
                lambda_cai_lr=lambda_cai_lr,
                lambda_cai_max=lambda_cai_max,
                cai_tolerance=cai_tolerance,
                smoothing_factor=smoothing_factor
            )

        elif constraint_type.lower() == 'ams':
            constraint = AminoMatchingSoftmax(
                amino_acid_sequence=amino_acid_sequence,
                batch_size=batch_size,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=lambda_cai,
                device=str(self.device),
                # Adaptive lambda_cai parameters
                adaptive_lambda_cai=adaptive_lambda_cai,
                lambda_cai_lr=lambda_cai_lr,
                lambda_cai_max=lambda_cai_max,
                cai_tolerance=cai_tolerance,
                smoothing_factor=smoothing_factor
            )
            return constraint.to(self.device)

        elif constraint_type.lower() == 'cpc':
            constraint = CodonProfileConstraint(
                amino_acid_sequence=amino_acid_sequence,
                batch_size=batch_size,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=lambda_cai,
                device=str(self.device),
                verbose=self.config.get('verbose', False),  # Pass verbose parameter
                # Adaptive lambda_cai parameters
                adaptive_lambda_cai=adaptive_lambda_cai,
                lambda_cai_lr=lambda_cai_lr,
                lambda_cai_max=lambda_cai_max,
                cai_tolerance=cai_tolerance,
                smoothing_factor=smoothing_factor
            )
            return constraint.to(self.device)

        else:
            raise ValueError(f"Unknown constraint type: {constraint_type}")

    def run_single_experiment(self,
                            protein_name: str,
                            constraint_type: str,
                            variant: str,
                            seed: int = 42,
                            show_progress: bool = False,
                            enable_deferred_validation: bool = False) -> Dict:
        """
        Run a single unified experiment with integrated CAI optimization.

        This implements the paper's unified loss: L_total = L_Access + λ_CAI * L_CAI
        
        Args:
            protein_name: Name of the protein
            constraint_type: Type of constraint (lagrangian, ams, cpc)
            variant: Variant code (00, 01, 10, 11)
            seed: Random seed
            show_progress: Whether to show optimization progress bar
        """
        start_time = time.time()
        torch.manual_seed(seed)

        try:
            # Load protein information
            protein_info = self.data_loader.get_protein_info(protein_name)
            amino_acid_sequence = protein_info['sequence']

            # Create constraint with unified CAI support
            constraint = self.create_constraint(
                constraint_type=constraint_type,
                amino_acid_sequence=amino_acid_sequence,
                variant=variant
            )

            # Optimization parameters from config
            iterations = self.config.get('iterations', 1000)
            learning_rate = self.config.get('learning_rate', 0.01)
            optimizer = torch.optim.AdamW(constraint.parameters(), lr=learning_rate)

            # Mixed precision setup (fully enabled or fully disabled)
            use_amp = self.config.get('mixed_precision', False) and self.device.type == 'cuda'
            scaler = torch.cuda.amp.GradScaler() if use_amp else None
            if use_amp:
                logger.info("⚡ Mixed precision training enabled (Full AMP mode)")
            gradient_clip = self.config.get('gradient_clip', 1.0)

            # Track optimization progress
            trajectory = {
                'iterations': [],
                'timestamps': [],
                'accessibility': [],
                'unified_loss': [],
                'cai_loss': [],
                'ecai_values': [],
                'discrete_cai_values': [],  # New: Actual CAI values for discrete sequences
                # New: Detailed trajectory data
                'rna_sequences': [],       # Probability distribution at each iteration
                'discrete_sequences': [],  # Discrete sequence at each iteration
                'accessibility_values': [], # Accessibility value at each iteration
                'loss_values': [],          # Loss value at each iteration
                'amp_enabled': use_amp,    # Record whether mixed precision is used
                # Deferred validation cache
                'deferred_sequences_cache': [] if enable_deferred_validation else None,
                'amp_scale_growth': []     # Record GradScaler scale changes
            }

            # Parse variant for alpha and beta values (fixed values, no annealing)
            alpha = 0.1 if variant[0] == '1' else 0.0  # Gumbel noise
            beta = 1.0 if variant[1] == '1' else 0.0   # Fixed beta value: 0=soft probability, 1=STE

            best_accessibility = float('inf')
            best_sequence = None
            best_seq_design = None  # Complete best design information

            # Precompute UTR tensors once (optimization)
            utr5_tensor = self.string_to_one_hot_tensor(protein_info['utr5'])
            utr3_tensor = self.string_to_one_hot_tensor(protein_info['utr3'])

            # Create inner progress bar (if needed)
            if show_progress:
                pbar = tqdm(total=iterations, 
                           desc=f"Optimizing", 
                           position=1, 
                           leave=False,
                           ncols=100,
                           bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')

            # Optimization loop
            for iteration in range(iterations):
                optimizer.zero_grad()

                # Beta value is fixed, no annealing (correcting erroneous design)

                if use_amp:
                    # Fully enabled mode - All computations use mixed precision (including beta=1)
                    with torch.cuda.amp.autocast():
                        # Forward pass through constraint (with precomputed UTR tensors)
                        result = self._forward_pass(
                            constraint, constraint_type, amino_acid_sequence,
                            alpha, beta, protein_info, utr5_tensor, utr3_tensor,
                            enable_discrete_monitoring=not enable_deferred_validation
                        )
                        total_loss = result['total_loss']
                    
                    # Backward pass with gradient scaling
                    scaler.scale(total_loss).backward()
                    
                    # Gradient clipping if enabled
                    if gradient_clip > 0:
                        scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(constraint.parameters(), gradient_clip)
                    
                    # Optimizer step with scaler
                    scaler.step(optimizer)
                    scaler.update()
                else:
                    # Original FP32 mode
                    # Forward pass through constraint (with precomputed UTR tensors)
                    result = self._forward_pass(
                        constraint, constraint_type, amino_acid_sequence,
                        alpha, beta, protein_info, utr5_tensor, utr3_tensor,
                        enable_discrete_monitoring=not enable_deferred_validation
                    )

                    # Backward pass
                    total_loss = result['total_loss']
                    total_loss.backward()

                    # Gradient clipping if enabled (also for FP32)
                    if gradient_clip > 0:
                        torch.nn.utils.clip_grad_norm_(constraint.parameters(), gradient_clip)
                    
                    optimizer.step()

                # Update Lagrangian multiplier if applicable
                if hasattr(constraint, 'update_lambda') and 'constraint_penalty' in result:
                    constraint.update_lambda(result['constraint_penalty'].item())

                # Track progress - save every iteration
                self._record_trajectory(
                    trajectory, iteration, start_time,
                    result, self.config.get('verbose', False)
                )
                
                # Record AMP scale growth every 100 iterations
                if use_amp and iteration % 100 == 0:
                    trajectory['amp_scale_growth'].append({
                        'iteration': iteration,
                        'scale': scaler.get_scale()
                    })

                # Track best result - Only sequences satisfying constraints can be recorded as best
                # Note: Result during iteration doesn't have amino_acids_match field, needs real-time validation
                from id3.utils.sequence_utils import rna_to_amino_acids
                discrete_seq = result.get('discrete_sequence', '')
                if discrete_seq and len(discrete_seq) == len(amino_acid_sequence) * 3:
                    translated = rna_to_amino_acids(discrete_seq)
                    amino_acids_match = (translated == amino_acid_sequence)
                else:
                    amino_acids_match = False
                
                if result['accessibility'] < best_accessibility and amino_acids_match:
                    best_accessibility = result['accessibility']
                    best_sequence = result['discrete_sequence']

                    # Save complete best design information
                    best_seq_design = {
                        'accessibility': result['accessibility'],
                        'discrete_sequence': result['discrete_sequence'],
                        'iteration': iteration,
                        'timestamp': time.time() - start_time,
                    }

                    # If CAI is enabled, add CAI-related information
                    if self.config.get('enable_cai', False) and 'loss_components' in result:
                        loss_components = result['loss_components']

                        # ECAI value (continuous optimization target)
                        best_seq_design['ecai'] = loss_components.get('ecai_value', None)
                        if best_seq_design['ecai'] and isinstance(best_seq_design['ecai'], torch.Tensor):
                            best_seq_design['ecai'] = best_seq_design['ecai'].item()

                        # Discrete CAI value (CAI of actual sequence)
                        if 'eval_cai' in loss_components:
                            best_seq_design['discrete_cai'] = loss_components['eval_cai']
                        elif 'discrete_cai' in loss_components:
                            best_seq_design['discrete_cai'] = loss_components['discrete_cai']
                        else:
                            best_seq_design['discrete_cai'] = None
                        
                        # CAI loss
                        best_seq_design['cai_loss'] = loss_components.get('cai_loss', None)
                        if best_seq_design['cai_loss'] and isinstance(best_seq_design['cai_loss'], torch.Tensor):
                            best_seq_design['cai_loss'] = best_seq_design['cai_loss'].item()

                    # Add other useful information
                    best_seq_design['total_loss'] = result['total_loss'].item()

                # Update inner progress bar
                if show_progress:
                    postfix = {
                        'access': f"{result['accessibility']:.4f}",
                        'loss': f"{result['total_loss'].item():.4f}"
                    }

                    # If CAI is enabled, add CAI information
                    if self.config.get('enable_cai', False) and 'loss_components' in result:
                        loss_components = result['loss_components']
                        if 'ecai_value' in loss_components:
                            postfix['ecai'] = f"{loss_components['ecai_value'].item():.4f}"
                    
                    pbar.set_postfix(postfix)
                    pbar.update(1)

            # Close inner progress bar
            if show_progress:
                pbar.close()

            # 🚀 Deferred validation batch processing
            if enable_deferred_validation:
                # Check if it's STE mode, STE mode doesn't need deferred validation
                if beta == 0.0:  # Only soft probability mode needs deferred validation
                    self._process_deferred_discrete_validation(trajectory, protein_info)
                else:
                    logger.debug(f"Skip deferred validation: STE mode (beta={beta}) already has correct monitoring values")

            # Final evaluation
            return self._prepare_final_result(
                constraint, protein_info, amino_acid_sequence,
                trajectory, best_accessibility, best_sequence,
                best_seq_design,  # New parameter
                protein_name, constraint_type, variant, seed,
                iterations, learning_rate, start_time
            )

        except Exception as e:
            return {
                'protein_name': protein_name,
                'constraint_type': constraint_type,
                'variant': variant,
                'seed': seed,
                'error': str(e),
                'status': 'failed',
                'optimization_time': time.time() - start_time
            }

    def _forward_pass(self, constraint, constraint_type, amino_acid_sequence,
                     alpha, beta, protein_info, utr5_tensor=None, utr3_tensor=None, 
                     enable_discrete_monitoring=True) -> Dict:
        """Execute forward pass and compute losses."""
        # Execute constraint forward pass
        if constraint_type.lower() == 'lagrangian':
            result = constraint.forward(alpha=alpha, beta=beta, tau=1.0, compute_penalty=True)
        else:  # AMS or CPC
            result = constraint.forward(alpha=alpha, beta=beta, tau=1.0)

        # Get discrete sequence for validation
        # Optimization: Get discrete sequence directly from result to avoid redundant calls
        discrete_sequence = result.get('discrete_sequence', '')

        # If discrete_sequence is not in result (compatible with old constraint classes)
        if not discrete_sequence or len(discrete_sequence) != len(amino_acid_sequence) * 3:
            # Fallback: convert soft probabilities to discrete
            rna_sequence = result['rna_sequence']
            if isinstance(rna_sequence, torch.Tensor):
                if rna_sequence.dim() == 3:
                    rna_sequence = rna_sequence[0]
                indices = torch.argmax(rna_sequence, dim=-1)
                discrete_sequence = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
            else:
                discrete_sequence = str(rna_sequence)

        # Get soft probabilities for gradient computation (Continuous Path from paper)
        rna_probs = result.get('rna_sequence')  # [length, 4] or [batch, length, 4]

        if rna_probs is not None:
            # Ensure batch dimension
            if rna_probs.dim() == 2:
                rna_probs = rna_probs.unsqueeze(0)  # [seq_len, 4] -> [1, seq_len, 4]

            # Use precomputed UTR tensors if provided, otherwise compute them
            if utr5_tensor is None:
                utr5_tensor = self.string_to_one_hot_tensor(protein_info['utr5'])  # [1, utr5_len, 4]
            if utr3_tensor is None:
                utr3_tensor = self.string_to_one_hot_tensor(protein_info['utr3'])  # [1, utr3_len, 4]

            # Concatenate full sequence: UTR5 + CDS (soft probs) + UTR3
            full_rna_probs = torch.cat([utr5_tensor, rna_probs, utr3_tensor], dim=1)

            # Compute accessibility loss using soft probabilities (Continuous Path - preserves gradients)
            accessibility_loss = self.deepraccess.compute_atg_window_accessibility(
                full_rna_probs,
                atg_position=len(protein_info['utr5']),
                discrete=False  # Use continuous mode for gradient flow
            )
        else:
            # Fallback: use discrete sequence if no soft probabilities available
            full_rna = protein_info['utr5'] + discrete_sequence + protein_info['utr3']
            full_rna_tensor = self.string_to_one_hot_tensor(full_rna)
            accessibility_loss = self.deepraccess.compute_atg_window_accessibility(
                full_rna_tensor,
                atg_position=len(protein_info['utr5']),
                discrete=True  # Discrete mode when no soft probs available
            )

        # 🔥 Optimization: Conditional discrete monitoring (deferred validation optimization)
        if enable_discrete_monitoring:
            # Check if it's STE mode (beta=1)
            if beta == 1.0:
                # STE mode: Both paths should be numerically consistent, use continuous path result directly
                accessibility = accessibility_loss.item()
                logger.debug(f"STE mode: Using continuous path result {accessibility:.6f}")
            else:
                # Soft probability mode: Needs real discrete validation
                with torch.no_grad():
                    full_rna_discrete = protein_info['utr5'] + discrete_sequence + protein_info['utr3']
                    full_rna_discrete_tensor = self.string_to_one_hot_tensor(full_rna_discrete)
                    accessibility_discrete = self.deepraccess.compute_atg_window_accessibility(
                        full_rna_discrete_tensor,
                        atg_position=len(protein_info['utr5']),
                        discrete=True  # Discrete mode for evaluation
                    )
                    accessibility = accessibility_discrete.item() if isinstance(accessibility_discrete, torch.Tensor) else accessibility_discrete
        else:
            # Optimization mode: Deferred inference, temporarily use continuous value
            accessibility = accessibility_loss.item()  # Temporary value

        # Compute unified loss
        # Check if constraint penalty is disabled (CAI no penalty mode)
        disable_constraint_penalty = self.config.get('disable_constraint_penalty', False)
        
        if constraint_type.lower() == 'lagrangian':
            # If constraint penalty is disabled, pass zero penalty
            constraint_penalty_to_use = torch.zeros_like(result['constraint_penalty']) if disable_constraint_penalty else result['constraint_penalty']
            
            loss_components = constraint.compute_total_loss(
                accessibility_loss,
                constraint_penalty_to_use,
                probabilities=result['probabilities'],
                enhanced_sequence=result.get('enhanced_sequence'),  # Dual-path architecture support
                cai_metadata=result.get('cai_metadata')  # 🚀 FIX: Pass CAI metadata
            )
        else:
            # For CPC, need to pass enhanced distribution
            if constraint_type.lower() == 'cpc':
                loss_components = constraint.compute_total_loss(
                    accessibility_loss,
                    codon_probs=result.get('codon_probs'),
                    enhanced_codon_dist=result.get('enhanced_codon_dist')
                )
            else:
                loss_components = constraint.compute_total_loss(
                    accessibility_loss,
                    codon_probs=result.get('codon_probs')
                )

        return {
            'total_loss': loss_components['total_loss'],
            'accessibility': accessibility,
            'discrete_sequence': discrete_sequence,
            'constraint_penalty': result.get('constraint_penalty'),
            'loss_components': loss_components,
            # New: Ensure rna_sequence is returned for trajectory saving
            'rna_sequence': result.get('rna_sequence'),  # Soft probability distribution
            # Deferred validation optimization information
            'deferred_validation_enabled': not enable_discrete_monitoring,
            'full_rna_discrete': protein_info['utr5'] + discrete_sequence + protein_info['utr3'] if not enable_discrete_monitoring else None,
        }

    def _record_trajectory(self, trajectory, iteration, start_time, result, verbose):
        """Record optimization trajectory."""
        trajectory['iterations'].append(iteration)
        trajectory['timestamps'].append(time.time() - start_time)
        trajectory['accessibility'].append(result['accessibility'])
        trajectory['unified_loss'].append(result['total_loss'].item())

        # New: Save detailed trajectory data
        # 1. Save probability distribution (rna_sequence)
        if 'rna_sequence' in result:
            rna_seq = result['rna_sequence']
            if isinstance(rna_seq, torch.Tensor):
                # Convert to numpy and save as list
                rna_seq_np = rna_seq.detach().cpu().numpy()
                if rna_seq_np.ndim > 2:  # If there's a batch dimension, take the first one
                    rna_seq_np = rna_seq_np[0]
                trajectory['rna_sequences'].append(rna_seq_np.tolist())
            else:
                trajectory['rna_sequences'].append(rna_seq)
        else:
            trajectory['rna_sequences'].append(None)

        # 2. Save discrete sequence
        if 'discrete_sequence' in result:
            trajectory['discrete_sequences'].append(result['discrete_sequence'])
        else:
            trajectory['discrete_sequences'].append(None)

        # 3. accessibility is already recorded at line 542, no need to repeat

        # 4. loss is already recorded as unified_loss at line 542, no need to repeat

        # 5. Deferred validation cache processing
        if trajectory.get('deferred_sequences_cache') is not None and result.get('full_rna_discrete'):
            trajectory['deferred_sequences_cache'].append(result['full_rna_discrete'])

        loss_components = result['loss_components']
        if self.config.get('enable_cai', False) and 'cai_loss' in loss_components:
            trajectory['cai_loss'].append(loss_components['cai_loss'].item())
            trajectory['ecai_values'].append(loss_components['ecai_value'].item())

            # Record actual CAI values for discrete sequences
            if 'eval_cai' in loss_components:
                # eval_cai returned by CPC/AMS constraint types
                trajectory['discrete_cai_values'].append(loss_components['eval_cai'])
            elif 'discrete_cai' in loss_components:
                # Other possible naming
                trajectory['discrete_cai_values'].append(loss_components['discrete_cai'])
            else:
                # If not provided, record as None or compute
                trajectory['discrete_cai_values'].append(None)
        else:
            trajectory['cai_loss'].append(0.0)
            trajectory['ecai_values'].append(0.0)
            trajectory['discrete_cai_values'].append(None)

        if verbose and iteration % max(1, self.config.get('iterations', 1000) // 5) == 0:
            cai_info = ""
            if self.config.get('enable_cai', False) and 'ecai_value' in loss_components:
                cai_info = f", ECAI: {loss_components['ecai_value'].item():.4f}"
            logger.debug(f"   Iter {iteration:4d}: Access={result['accessibility']:.4f}, "
                  f"Loss={result['total_loss'].item():.4f}{cai_info}")

    def _process_deferred_discrete_validation(self, trajectory, protein_info):
        """🚀 Batch process deferred discrete validation"""

        sequences_cache = trajectory.get('deferred_sequences_cache')
        if not sequences_cache:
            logger.warning("⚠️ Deferred validation cache is empty, skipping batch processing")
            return

        logger.info(f"🔄 Starting batch inference for {len(sequences_cache)} discrete sequences...")
        
        batch_size = 32  # Configurable
        total_sequences = len(sequences_cache)
        batch_accessibilities = []

        with torch.no_grad():
            for i in range(0, total_sequences, batch_size):
                batch_end = min(i + batch_size, total_sequences)
                batch_sequences = sequences_cache[i:batch_end]

                # Build batch input
                batch_tensors = []
                for full_rna_seq in batch_sequences:
                    tensor = self.string_to_one_hot_tensor(full_rna_seq)
                    batch_tensors.append(tensor)

                # Batch inference
                if batch_tensors:
                    batch_input = torch.cat(batch_tensors, dim=0)
                    batch_results = self.deepraccess.compute_atg_window_accessibility(
                        batch_input,
                        atg_position=len(protein_info['utr5']),
                        discrete=True
                    )

                    # Convert to scalar list
                    for result in batch_results:
                        acc_value = result.item() if isinstance(result, torch.Tensor) else result
                        batch_accessibilities.append(acc_value)

        # Update accessibility values in trajectory
        if len(batch_accessibilities) == len(trajectory['accessibility']):
            trajectory['accessibility'] = batch_accessibilities
            logger.info(f"✅ Batch inference completed, updated {len(batch_accessibilities)} accessibility values")
        else:
            logger.warning(f"⚠️ Batch inference count mismatch: {len(batch_accessibilities)} vs {len(trajectory['accessibility'])}")

    def _prepare_final_result(self, constraint, protein_info, amino_acid_sequence,
                             trajectory, best_accessibility, best_sequence,
                             best_seq_design,  # New parameter
                             protein_name, constraint_type, variant, seed,
                             iterations, learning_rate, start_time) -> Dict:
        """Prepare final experiment result."""
        optimization_time = time.time() - start_time

        # Final evaluation using discrete sequence
        final_sequence = constraint.get_discrete_sequence()
        final_amino_acids = rna_to_amino_acids(final_sequence)
        amino_acids_match = final_amino_acids == amino_acid_sequence

        # Final accessibility using discrete sequence
        final_full_rna = protein_info['utr5'] + final_sequence + protein_info['utr3']
        final_full_rna_tensor = self.string_to_one_hot_tensor(final_full_rna)
        final_accessibility = self.deepraccess.compute_atg_window_accessibility(
            final_full_rna_tensor, atg_position=len(protein_info['utr5'])
        )
        # Convert tensor output to scalar
        final_accessibility = final_accessibility.item() if isinstance(final_accessibility, torch.Tensor) else final_accessibility

        # Prepare result
        result_dict = {
            'protein_name': protein_name,
            'constraint_type': constraint_type,
            'variant': variant,
            'seed': seed,
            'iterations': iterations,
            'learning_rate': learning_rate,
            'optimization_time': optimization_time,
            'iterations_per_second': iterations / optimization_time,

            'initial_accessibility': trajectory['accessibility'][0] if trajectory['accessibility'] else 0.0,
            'final_accessibility': final_accessibility,
            'improvement': trajectory['accessibility'][0] - final_accessibility if trajectory['accessibility'] else 0.0,
            'best_accessibility': best_accessibility,
            'best_seq_design': best_seq_design,  # Complete best design information

            'amino_acids_match': amino_acids_match,
            'amino_acids_correct': 100.0 if amino_acids_match else 0.0,
            'expected_amino_acids': amino_acid_sequence,
            'actual_amino_acids': final_amino_acids,
            'final_sequence': final_sequence,

            'unified_optimization': True,
            'cai_enabled': self.config.get('enable_cai', False),
            'trajectory': trajectory,
            'status': 'completed'
        }

        # Add CAI-specific results
        if self.config.get('enable_cai', False) and trajectory['ecai_values']:
            # Get final discrete CAI value to determine if target is achieved
            final_discrete_cai = None
            if trajectory['discrete_cai_values'] and trajectory['discrete_cai_values'][-1] is not None:
                final_discrete_cai = trajectory['discrete_cai_values'][-1]
            
            result_dict.update({
                'cai_target': self.config.get('cai_target', 0.8),
                'lambda_cai': self.config.get('lambda_cai', 0.1),
                'final_ecai': trajectory['ecai_values'][-1],
                'initial_ecai': trajectory['ecai_values'][0],
                'ecai_improvement': trajectory['ecai_values'][-1] - trajectory['ecai_values'][0],
                'final_discrete_cai': final_discrete_cai,
                # Use actual discrete CAI value to determine if target is achieved, not ECAI value
                'cai_target_achieved': final_discrete_cai >= self.config.get('cai_target', 0.8) if final_discrete_cai is not None else False
            })

        return result_dict

    def _load_existing_progress(self) -> set:
        """
        Load the set of completed experiment IDs.

        Returns:
            Set of completed experiment IDs
        """
        completed = set()
        progress_file = self.output_dir / 'progress.json'
        
        if progress_file.exists():
            try:
                with open(progress_file, 'r') as f:
                    progress = json.load(f)

                # Extract completed experiment IDs
                for exp in progress.get('experiments', []):
                    if exp.get('status') == 'completed':
                        completed.add(exp['id'])

                if completed:
                    logger.info(f"📂 Loaded progress from {self.output_dir.name}: {len(completed)} experiments completed")

            except Exception as e:
                logger.warning(f"⚠️ Unable to load progress file: {e}")
        
        return completed
    
    def _get_output_dir(self) -> Path:
        """
        Get the output directory path.

        Returns:
            Path to output directory
        """
        # Fix: Get output_dir from dictionary config
        if self.config.get('output_dir'):
            return Path(self.config['output_dir'])

        # Fix: Get enable_cai from dictionary config and generate correct mode name
        enable_cai = self.config.get('enable_cai', False)
        mode = 'unified_cai_experiments' if enable_cai else 'unified_access_experiments'
        if hasattr(self, '_output_timestamp'):
            timestamp = self._output_timestamp
        else:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            self._output_timestamp = timestamp

        # Timestamp at the beginning
        return Path(f'results/{timestamp}_{mode}')

    def _init_or_update_progress_tracker(self, total_experiments: int, already_completed: int = 0) -> None:
        """
        Initialize or update progress tracking file.

        Args:
            total_experiments: Total number of experiments
            already_completed: Number of already completed experiments
        """
        progress_file = self.output_dir / 'progress.json'

        # If file exists, preserve existing data
        if progress_file.exists():
            try:
                with open(progress_file, 'r') as f:
                    progress_data = json.load(f)
                # Update total (user may have changed experiment configuration)
                progress_data['total_experiments'] = total_experiments
                progress_data['last_update'] = datetime.now().isoformat()
                logger.info(f"📊 Updated progress tracking: {progress_file}")
            except Exception as e:
                logger.warning(f"⚠️ Unable to read progress file, creating new one: {e}")
                progress_data = self._create_new_progress_data(total_experiments, already_completed)
        else:
            progress_data = self._create_new_progress_data(total_experiments, already_completed)
            logger.info(f"📊 Initialized progress tracking: {progress_file}")
        
        try:
            with open(progress_file, 'w') as f:
                json.dump(progress_data, f, indent=2)
                f.flush()
                os.fsync(f.fileno())
        except Exception as e:
            logger.warning(f"⚠️ Unable to create progress file: {str(e)}")

    def _create_new_progress_data(self, total_experiments: int, already_completed: int = 0) -> dict:
        """
        Create new progress data structure.
        """
        return {
            'start_time': datetime.now().isoformat(),
            'total_experiments': total_experiments,
            'completed_experiments': already_completed,
            'failed_experiments': 0,
            'current_experiment': None,
            'experiments': [],
            'last_update': datetime.now().isoformat(),
            'version': '1.1'  # Add version number for future upgrades
        }

    def _init_progress_tracker(self, total_experiments: int) -> None:
        """
        Initialize progress tracking file (for compatibility).

        Args:
            total_experiments: Total number of experiments
        """
        self._init_or_update_progress_tracker(total_experiments, 0)
        
        try:
            with open(progress_file, 'w') as f:
                json.dump(progress_data, f, indent=2)
                f.flush()
                os.fsync(f.fileno())
            logger.info(f"📊 Initialized progress tracking: {progress_file}")
        except Exception as e:
            logger.warning(f"⚠️ Unable to create progress file: {str(e)}")

    def _update_progress(self, status: str, exp: dict, index: int,
                        saved_path: Path = None, error: str = None) -> None:
        """
        Update progress tracking file.

        Args:
            status: Experiment status ('running', 'completed', 'failed')
            exp: Experiment configuration dictionary
            index: Current experiment index
            saved_path: Saved file path (completed status only)
            error: Error message (failed status only)
        """
        progress_file = self.output_dir / 'progress.json'
        
        try:
            # Read existing progress data
            if progress_file.exists():
                with open(progress_file, 'r') as f:
                    progress_data = json.load(f)
            else:
                # If file doesn't exist, create new one
                progress_data = {
                    'start_time': datetime.now().isoformat(),
                    'total_experiments': 0,
                    'completed_experiments': 0,
                    'failed_experiments': 0,
                    'current_experiment': None,
                    'experiments': [],
                    'last_update': datetime.now().isoformat()
                }
            
            exp_id = f"{exp['protein_name']}-{exp['constraint_type']}-{exp['variant']}-{exp['seed']}"
            
            if status == 'running':
                progress_data['current_experiment'] = {
                    'id': exp_id,
                    'index': index,
                    'start_time': datetime.now().isoformat()
                }
            elif status == 'completed':
                progress_data['completed_experiments'] += 1
                progress_data['current_experiment'] = None
                progress_data['experiments'].append({
                    'id': exp_id,
                    'index': index,
                    'status': 'completed',
                    'file': saved_path.name if saved_path else None,
                    'timestamp': datetime.now().isoformat()
                })
            elif status == 'failed':
                progress_data['failed_experiments'] += 1
                progress_data['current_experiment'] = None
                progress_data['experiments'].append({
                    'id': exp_id,
                    'index': index,
                    'status': 'failed',
                    'error': error,
                    'timestamp': datetime.now().isoformat()
                })

            progress_data['last_update'] = datetime.now().isoformat()

            # Calculate progress percentage
            total = progress_data.get('total_experiments', 0)
            if total > 0:
                completed = progress_data['completed_experiments']
                failed = progress_data['failed_experiments']
                progress_pct = ((completed + failed) / total) * 100
                progress_data['progress_percentage'] = round(progress_pct, 2)

            # Atomic write (write to temp file then rename)
            temp_file = progress_file.with_suffix('.tmp')
            with open(temp_file, 'w') as f:
                json.dump(progress_data, f, indent=2)
                f.flush()
                os.fsync(f.fileno())

            # Atomic replacement
            temp_file.replace(progress_file)

            if self.config.get('verbose', False):
                logger.debug(f"📊 Updated progress: {exp_id} - {status}")

        except Exception as e:
            logger.warning(f"⚠️ Unable to update progress file: {str(e)}")

    def _save_experiment_result(self, result: dict) -> Path:
        """
        Incrementally save individual experiment result, ensuring immediate disk write.

        Args:
            result: Experiment result dictionary

        Returns:
            Path to saved file
        """
        # Generate filename: timestamp_protein_constraint_variant_seed.json
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        filename = (f"{timestamp}_{result['protein_name']}_{result['constraint_type']}_"
                   f"{result['variant']}_seed{result['seed']}.json")

        file_path = self.output_dir / filename

        # Serialize numpy arrays
        def serialize_numpy(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, (np.float32, np.float64)):
                return float(obj)
            elif isinstance(obj, (np.int32, np.int64)):
                return int(obj)
            raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

        # Save result and immediately flush to disk
        try:
            with open(file_path, 'w') as f:
                json.dump(result, f, indent=2, default=serialize_numpy)
                f.flush()  # Immediately flush buffer
                os.fsync(f.fileno())  # Force write to disk

            # Provide more detailed save information
            exp_id = f"{result['protein_name']}-{result['constraint_type']}-{result['variant']}-{result['seed']}"
            logger.info(f"✅ Saved experiment result: {filename} | Accessibility: {result.get('final_accessibility', 'N/A'):.4f}")

        except Exception as e:
            logger.error(f"❌ Failed to save experiment result: {filename} - {str(e)}")
            raise

        return file_path

    def run_batch(self, experiments: list, num_workers: int = None) -> list:
        """
        Run a batch of experiments with serial execution (optimal performance).

        Based on performance test results, serial execution is 2-5x faster than parallel, so always use serial mode.

        Args:
            experiments: List of experiment specifications
            num_workers: Ignored parameter, always use serial execution

        Returns:
            List of experiment results
        """
        logger.info(f"🚀 Serial execution of {len(experiments)} experiments (optimal performance mode)")
        logger.info(f"💡 Performance tests show: Serial is 2-5x faster than parallel")

        # Directly call serial execution method
        return self._run_batch_sequential(experiments)

    def _run_batch_sequential(self, experiments: list) -> list:
        """
        Run experiments sequentially with progress bar and incremental saving.
        Supports resuming from previous runs.
        """
        results = []

        # Filter completed experiments
        experiments_to_run = []
        skipped_count = 0

        for exp in experiments:
            exp_id = f"{exp['protein_name']}-{exp['constraint_type']}-{exp['variant']}-{exp['seed']}"
            if exp_id in self.completed_experiments:
                skipped_count += 1
                logger.debug(f"Skipping completed experiment: {exp_id}")
            else:
                experiments_to_run.append(exp)

        if skipped_count > 0:
            logger.info(f"⏭️  Skipped {skipped_count} completed experiments")
            logger.info(f"🚀 Will run {len(experiments_to_run)} remaining experiments")

        if not experiments_to_run:
            logger.info("✅ All experiments are already completed!")
            return results

        # Initialize or update progress tracking
        self._init_or_update_progress_tracker(len(experiments), len(self.completed_experiments))

        # Outer progress bar - experiment progress
        # Calculate total number of experiments (including completed ones)
        total_experiments = len(experiments)
        
        with tqdm(total=len(experiments_to_run), 
                  desc="Overall Experiments", 
                  position=0, 
                  leave=True,
                  ncols=100,
                  bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]',
                  initial=0) as pbar_experiments:
            
            for i, exp in enumerate(experiments_to_run, 1):
                exp_id = f"{exp['protein_name']}-{exp['constraint_type']}-{exp['variant']}-{exp['seed']}"

                # Calculate actual experiment index (including skipped ones)
                actual_index = i + skipped_count

                # Update progress bar description
                pbar_experiments.set_description(f"Exp: {exp_id} [{actual_index}/{total_experiments}]")

                # Update progress: starting experiment
                self._update_progress('running', exp, actual_index)

                # Run single experiment (pass show_progress parameter to enable inner progress bar)
                result = self.run_single_experiment(
                    protein_name=exp['protein_name'],
                    constraint_type=exp['constraint_type'],
                    variant=exp['variant'],
                    seed=exp['seed'],
                    show_progress=not (self.config.get('disable_inner_tqdm', False) if isinstance(self.config, dict) else getattr(self.config, 'disable_inner_tqdm', False))  # Determine whether to enable inner progress bar based on configuration
                )

                # Incrementally save experiment result
                if result.get('status') == 'completed':
                    saved_path = self._save_experiment_result(result)

                    # Update progress: experiment completed
                    self._update_progress('completed', exp, actual_index, saved_path=saved_path)

                    # Update outer progress bar suffix information
                    postfix = {
                        'access': f"{result['final_accessibility']:.4f}",
                        'AA': f"{result['amino_acids_correct']:.0f}%"
                    }

                    # If CAI is enabled, add CAI information
                    if self.config.get('enable_cai', False) and 'final_ecai' in result:
                        postfix['CAI'] = f"{result['final_ecai']:.4f}"
                        if result.get('cai_target_achieved', False):
                            postfix['CAI'] += "✓"

                    pbar_experiments.set_postfix(postfix)
                else:
                    # Update progress: experiment failed
                    self._update_progress('failed', exp, actual_index, error=result.get('error', 'Unknown error'))

                    # Display error when experiment fails
                    pbar_experiments.set_postfix({'status': 'FAILED'})
                    tqdm.write(f"❌ Failed: {exp_id} - {result.get('error', 'Unknown error')}")

                results.append(result)
                pbar_experiments.update(1)

        logger.info(f"💾 All experiment results saved to: {self.output_dir}")
        logger.info(f"📊 Progress tracking file: {self.output_dir / 'progress.json'}")

        return results

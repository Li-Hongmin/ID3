#!/usr/bin/env python3
"""
Ray Worker for ID3-DeepRaccess Experiments





"""

import ray
import torch
import logging
import traceback
from typing import Dict, Any, Optional
from pathlib import Path
import sys
import time
import numpy as np

logger = logging.getLogger(__name__)


@ray.remote
class ExperimentWorker:
    """

    

    """
    
    def __init__(self, num_gpus: float = 0.25):
        """

        

        """

        project_root = Path(__file__).parent.parent.parent.parent
        sys.path.insert(0, str(project_root))
        

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - Worker - %(levelname)s - %(message)s'
        )
        

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.num_gpus = num_gpus
        
        logger.info(f"Worker initialized on {self.device} with {num_gpus} GPU fraction")
        

        from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
        

        self.deepraccess = DeepRaccessID3Wrapper()
        logger.info("DeepRaccess model loaded successfully")
        

        self._import_modules()
    
    def _import_modules(self):


        global UnifiedExperimentConfig, UnifiedExperimentRunner
        global LagrangianConstraint, AminoMatchingSoftmax, CodonProfileConstraint
        global DiscreteCAISearcher, ProteinDataLoader
        
        from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
        from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
        from id3.constraints.lagrangian import LagrangianConstraint
        from id3.constraints.amino_matching import AminoMatchingSoftmax
        from id3.constraints.codon_profile import CodonProfileConstraint
        from id3.cai.discrete_binary_search import DiscreteCAISearcher
        from id3.experiments.utils.data_loader import ProteinDataLoader
    
    def run_experiment(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        运行单个实验
        
        Args:
            config: 实验配置字典，包含：
                - protein: 蛋白质名称
                - constraint: 约束类型
                - variant: 变体
                - iterations: 迭代次数
                - seed: 随机种子
                - enable_cai: 是否启用CAI
                - cai_target: CAI目标值
                - lambda_cai: CAI权重
                
        Returns:
            实验结果字典
        """
        try:
            start_time = time.time()
            

            protein_name = config['protein']
            constraint_type = config['constraint']
            variant = config['variant']
            iterations = config.get('iterations', 1000)
            seed = config.get('seed', 42)

            enable_cai = config.get('enable_cai', False)
            
            logger.info(f"Starting experiment: {protein_name}_{constraint_type}_{variant}")
            

            torch.manual_seed(seed)
            np.random.seed(seed)
            

            data_loader = ProteinDataLoader()
            protein_info = data_loader.get_protein_info(protein_name)
            
            if protein_info is None:
                raise ValueError(f"Failed to load protein: {protein_name}")
            
            amino_acid_sequence = protein_info['sequence']
            

            constraint = self._create_constraint(
                constraint_type=constraint_type,
                amino_acid_sequence=amino_acid_sequence,
                variant=variant,
                enable_cai=enable_cai,
                cai_target=config.get('cai_target', 0.8),
                cai_weight=config.get('lambda_cai', 0.1),
                cai_species=config.get('cai_species', 'ecoli_bl21de3')
            )
            
            # Import needed for conversion
            from id3.utils.constants import NUCLEOTIDES, NUCLEOTIDE_MAP
            from id3.utils.sequence_utils import rna_to_amino_acids
            from id3.config.utr_loader import get_default_utrs
            

            optimizer = torch.optim.Adam(constraint.parameters(), lr=0.01)
            

            utrs = get_default_utrs()
            utr5 = protein_info.get('utr5', utrs['utr5'])
            utr3 = protein_info.get('utr3', utrs['utr3'])
            

            utr5_tensor = self._string_to_one_hot_tensor(utr5)
            utr3_tensor = self._string_to_one_hot_tensor(utr3)
            if utr5_tensor.dim() == 2:
                utr5_tensor = utr5_tensor.unsqueeze(0)
            if utr3_tensor.dim() == 2:
                utr3_tensor = utr3_tensor.unsqueeze(0)
            



            
            best_accessibility = float('inf')
            best_sequence = None
            best_iteration = 0
            best_seq_design = None
            

            use_amp = torch.cuda.is_available() and config.get('use_amp', True)
            scaler = torch.cuda.amp.GradScaler() if use_amp else None
            gradient_clip = config.get('gradient_clip', 1.0)
            

            trajectory = {
                'iterations': [],
                'accessibility': [],
                'loss': [],
                'discrete_sequences': []
            }
            
            for iteration in range(iterations):
                optimizer.zero_grad()
                
                if use_amp:
                    with torch.cuda.amp.autocast():

                        if constraint_type == 'lagrangian':
                            result = constraint.forward(alpha=alpha, beta=beta, tau=1.0, compute_penalty=True)
                        else:
                            result = constraint.forward(alpha=alpha, beta=beta, tau=1.0)
                        

                        rna_probs = result.get('rna_sequence')
                        if rna_probs.dim() == 2:
                            rna_probs = rna_probs.unsqueeze(0)  # [1, length, 4]
                        

                        full_rna_probs = torch.cat([
                            utr5_tensor,
                            rna_probs,
                            utr3_tensor
                        ], dim=1)
                        

                        accessibility_loss = self.deepraccess.compute_atg_window_accessibility(
                            full_rna_probs,
                            atg_position=len(utr5),

                        )
                        

                        if constraint_type == 'lagrangian':
                            loss_components = constraint.compute_total_loss(
                                accessibility_loss,
                                result.get('constraint_penalty'),
                                probabilities=result.get('probabilities'),
                                enhanced_sequence=result.get('enhanced_sequence'),
                                cai_metadata=result.get('cai_metadata')
                            )
                        else:
                            loss_components = constraint.compute_total_loss(
                                accessibility_loss,
                                codon_probs=result.get('codon_probs')
                            )
                        
                        total_loss = loss_components['total_loss']
                    

                    scaler.scale(total_loss).backward()
                    

                    if gradient_clip > 0:
                        scaler.unscale_(optimizer)
                        torch.nn.utils.clip_grad_norm_(constraint.parameters(), gradient_clip)
                    

                    scaler.step(optimizer)
                    scaler.update()
                else:

                    if constraint_type == 'lagrangian':
                        result = constraint.forward(alpha=alpha, beta=beta, tau=1.0, compute_penalty=True)
                    else:
                        result = constraint.forward(alpha=alpha, beta=beta, tau=1.0)
                    

                    rna_probs = result.get('rna_sequence')
                    if rna_probs.dim() == 2:
                        rna_probs = rna_probs.unsqueeze(0)
                    

                    full_rna_probs = torch.cat([
                        utr5_tensor,
                        rna_probs,
                        utr3_tensor
                    ], dim=1)
                    

                    accessibility_loss = self.deepraccess.compute_atg_window_accessibility(
                        full_rna_probs,
                        atg_position=len(utr5),
                        discrete=False
                    )
                    

                    if constraint_type == 'lagrangian':
                        loss_components = constraint.compute_total_loss(
                            accessibility_loss,
                            result.get('constraint_penalty'),
                            probabilities=result.get('probabilities'),
                            enhanced_sequence=result.get('enhanced_sequence'),
                            cai_metadata=result.get('cai_metadata')
                        )
                    else:
                        loss_components = constraint.compute_total_loss(
                            accessibility_loss,
                            codon_probs=result.get('codon_probs')
                        )
                    
                    total_loss = loss_components['total_loss']
                    total_loss.backward()
                    

                    if gradient_clip > 0:
                        torch.nn.utils.clip_grad_norm_(constraint.parameters(), gradient_clip)
                    
                    optimizer.step()
                

                if hasattr(constraint, 'update_lambda') and 'constraint_penalty' in result:
                    constraint.update_lambda(result['constraint_penalty'].item())
                

                discrete_sequence = result.get('discrete_sequence', '')
                if discrete_sequence and len(discrete_sequence) == len(amino_acid_sequence) * 3:
                    translated = rna_to_amino_acids(discrete_sequence)
                    amino_acids_match = (translated == amino_acid_sequence)
                    

                    full_rna_discrete = utr5 + discrete_sequence + utr3
                    full_rna_discrete_tensor = self._string_to_one_hot_tensor(full_rna_discrete)

                    if full_rna_discrete_tensor.dim() == 2:
                        full_rna_discrete_tensor = full_rna_discrete_tensor.unsqueeze(0)
                    accessibility_discrete = self.deepraccess.compute_atg_window_accessibility(
                        full_rna_discrete_tensor,
                        atg_position=len(utr5),
                        discrete=True
                    )
                    current_accessibility = accessibility_discrete.item() if isinstance(accessibility_discrete, torch.Tensor) else accessibility_discrete
                else:
                    amino_acids_match = False
                    current_accessibility = float('inf')
                

                trajectory['iterations'].append(iteration)
                trajectory['accessibility'].append(current_accessibility if amino_acids_match else float('inf'))
                trajectory['loss'].append(total_loss.item())
                trajectory['discrete_sequences'].append(discrete_sequence if amino_acids_match else '')
                

                if current_accessibility < best_accessibility and amino_acids_match:
                    best_accessibility = current_accessibility
                    best_sequence = discrete_sequence
                    best_iteration = iteration
                    

                    best_seq_design = {
                        'accessibility': current_accessibility,
                        'discrete_sequence': discrete_sequence,
                        'iteration': iteration,
                        'total_loss': total_loss.item()
                    }
                    

                    if enable_cai and 'loss_components' in locals():
                        if 'ecai_value' in loss_components:
                            best_seq_design['ecai'] = loss_components['ecai_value'].item() if hasattr(loss_components['ecai_value'], 'item') else loss_components['ecai_value']
                        if 'discrete_cai' in loss_components:
                            best_seq_design['discrete_cai'] = loss_components['discrete_cai']
            

            final_accessibility = best_accessibility if best_sequence else float('inf')
            final_amino_acids_match = best_sequence is not None
            

            optimization_time = time.time() - start_time
            

            initial_accessibility = trajectory['accessibility'][0] if trajectory.get('accessibility') else 0.0
            

            final_amino_acids = rna_to_amino_acids(best_sequence) if best_sequence else ''
            
            result = {
                'protein_name': protein_name,
                'constraint_type': constraint_type,
                'variant': variant,
                'seed': seed,
                'iterations': iterations,
                'learning_rate': learning_rate,
                'optimization_time': optimization_time,
                'iterations_per_second': iterations / optimization_time if optimization_time > 0 else 0.0,
                
                # Accessibility metrics
                'initial_accessibility': initial_accessibility,
                'final_accessibility': final_accessibility,
                'improvement': initial_accessibility - final_accessibility if initial_accessibility != 0.0 else 0.0,
                'best_accessibility': best_accessibility,

                
                # Amino acid validation
                'amino_acids_match': final_amino_acids_match,
                'amino_acids_correct': 100.0 if final_amino_acids_match else 0.0,
                'expected_amino_acids': amino_acid_sequence,
                'actual_amino_acids': final_amino_acids,
                'final_sequence': best_sequence if best_sequence else '',
                
                # Optimization metadata
                'unified_optimization': True,

                'trajectory': trajectory,

            }
            

            if best_seq_design:
                result.update({
                    'best_seq_design': best_seq_design,
                    'ecai': best_seq_design.get('ecai'),
                    'discrete_cai': best_seq_design.get('discrete_cai')
                })
            

            if enable_cai and trajectory.get('ecai_values'):
                result.update({
                    'cai_target': config.get('cai_target', 0.8),
                    'lambda_cai': config.get('lambda_cai', 0.1),
                    'final_ecai': trajectory['ecai_values'][-1] if trajectory['ecai_values'] else 0.0,
                    'initial_ecai': trajectory['ecai_values'][0] if trajectory['ecai_values'] else 0.0,
                    'ecai_improvement': (trajectory['ecai_values'][-1] - trajectory['ecai_values'][0]) if trajectory['ecai_values'] else 0.0,
                    'cai_target_achieved': best_seq_design.get('discrete_cai', 0.0) >= config.get('cai_target', 0.8) if best_seq_design else False
                })
            elif enable_cai:

                result.update({
                    'cai_target': config.get('cai_target', 0.8),
                    'lambda_cai': config.get('lambda_cai', 0.1),
                    'cai_species': config.get('cai_species', 'ecoli_bl21de3')
                })
            
            logger.info(f"Completed: {protein_name}_{constraint_type}_{variant}, "
                       f"Best: {best_accessibility:.4f} kcal/mol")
            
            return result
            
        except Exception as e:
            logger.error(f"Experiment failed: {str(e)}")
            logger.error(traceback.format_exc())
            

            return {
                'protein_name': config.get('protein', 'unknown'),
                'constraint_type': config.get('constraint', 'unknown'),
                'variant': config.get('variant', 'unknown'),
                'seed': config.get('seed', 42),
                'iterations': config.get('iterations', 0),
                'learning_rate': config.get('learning_rate', 0.05),
                'optimization_time': time.time() - start_time if 'start_time' in locals() else 0.0,
                'iterations_per_second': 0.0,
                

                'initial_accessibility': float('inf'),
                'final_accessibility': float('inf'),
                'improvement': 0.0,
                'best_accessibility': float('inf'),
                

                'amino_acids_match': False,
                'amino_acids_correct': 0.0,
                'expected_amino_acids': '',
                'actual_amino_acids': '',
                'final_sequence': '',
                
                # Optimization metadata
                'unified_optimization': True,
                'cai_enabled': config.get('enable_cai', False),
                'trajectory': {},
                'status': 'failed',
                'error': str(e),
                'traceback': traceback.format_exc()
            }
    
    def _string_to_one_hot_tensor(self, sequence: str) -> torch.Tensor:
        """将RNA序列转换为one-hot张量"""
        from id3.utils.constants import NUCLEOTIDE_MAP
        
        one_hot = torch.zeros(len(sequence), 4, device=self.device)
        for i, nuc in enumerate(sequence):
            if nuc in NUCLEOTIDE_MAP:
                one_hot[i, NUCLEOTIDE_MAP[nuc]] = 1.0
        return one_hot
    
    def _create_constraint(self, constraint_type: str, amino_acid_sequence: str,
                          variant: str, enable_cai: bool = False,
                          cai_target: float = 0.8, cai_weight: float = 0.1,
                          cai_species: str = 'ecoli_bl21de3'):

        constraint_type = constraint_type.lower()
        

        cai_enhancement_operator = None
        if enable_cai:
            from id3.cai.discrete_binary_search import DiscreteCAISearcher
            cai_enhancement_operator = DiscreteCAISearcher(
                species=cai_species,
                target_cai=cai_target,
                device=self.device
            )
        
        if constraint_type == 'lagrangian':
            return LagrangianConstraint(
                amino_acid_sequence=amino_acid_sequence,
                variant=variant,
                device=self.device,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=cai_weight,
                cai_enhancement_operator=cai_enhancement_operator
            )
        elif constraint_type == 'ams':
            return AminoMatchingSoftmax(
                amino_acid_sequence=amino_acid_sequence,
                variant=variant,
                device=self.device,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=cai_weight,
                cai_enhancement_operator=cai_enhancement_operator
            )
        elif constraint_type == 'cpc':
            return CodonProfileConstraint(
                amino_acid_sequence=amino_acid_sequence,
                variant=variant,
                device=self.device,
                enable_cai=enable_cai,
                cai_target=cai_target,
                cai_weight=cai_weight,
                cai_enhancement_operator=cai_enhancement_operator
            )
        else:
            raise ValueError(f"Unknown constraint type: {constraint_type}")
    
    def test_connection(self) -> Dict[str, Any]:
        """测试Worker连接和GPU状态"""
        return {
            'status': 'connected',
            'device': str(self.device),
            'gpu_available': torch.cuda.is_available(),
            'gpu_count': torch.cuda.device_count() if torch.cuda.is_available() else 0,
            'num_gpus_allocated': self.num_gpus
        }
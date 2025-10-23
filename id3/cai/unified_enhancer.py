"""
Unified CAI Enhancer (KISS Simplified Version)

Simple CAI enhancement for codon probability distributions.
Removed complex binary search optimizer that broke amino acid constraints.
"""

import torch
import torch.nn.functional as F
import numpy as np
from typing import Dict, Optional, Tuple

# Import existing conversion functions
from id3.utils.constants import codon_to_rna_matrix, amino_acids_to_codons, codons
from id3.cai.probability import rna_to_codon_probabilities
from id3.constraints.codon_profile import CodonReconstruction
from id3.cai.discrete_binary_search import DiscreteCAISearcher
from id3.cai.incremental_cai_calculator import IncrementalCAICalculator, CachedIncrementalCalculator


class UnifiedCAIEnhancer:
    """
    Simplified CAI enhancement using probability reweighting.
    
    KISS Principle: Instead of complex binary search, simply reweight
    codon probabilities by their CAI values and renormalize.
    
    This preserves amino acid constraints while improving CAI scores.
    """
    
    def __init__(self,
                 cai_target: float = 0.8,
                 cai_weight: float = 0.1,  # Ignored in KISS version
                 device: str = 'cuda',
                 species: str = 'ecoli_bl21de3',
                 enable_cai: bool = True):
        """Initialize simplified CAI enhancer."""
        
        self.cai_target = cai_target
        self.cai_weight = cai_weight  # Kept for compatibility but ignored
        self.species = species
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')
        self.enable_cai = enable_cai
        
        # CAI enhancement using simple probability reweighting (KISS approach)
        if self.enable_cai:
            self.cai_weights = self._load_cai_weights()
        else:
            self.cai_weights = None
            
        # Cache conversion matrices
        self.codon_to_rna = codon_to_rna_matrix.to(self.device)
        self.amino_to_codons = amino_acids_to_codons
        
        # Initialize discrete binary searcher for more efficient CAI optimization
        # This is now the DEFAULT method due to 22x performance improvement
        if self.enable_cai and self.cai_weights:
            self.discrete_searcher = DiscreteCAISearcher(
                self.cai_weights, 
                str(self.device),
                use_incremental=True  # Enable incremental by default
            )
    
    def _load_cai_weights(self):
        """Load CAI weights in correct codon order."""
        import json
        import os
        from id3.utils.constants import codons
        
        # Load E. coli BL21(DE3) CAI weights
        cai_weights_file = os.path.join(
            os.path.dirname(__file__), '..', '..', 'data', 'codon_references',
            'ecoli_bl21de3_wi_weights_comparison.json'
        )
        
        if os.path.exists(cai_weights_file):
            with open(cai_weights_file, 'r') as f:
                data = json.load(f)
                wi_table = data['wi_table']
                

            weights = []
            for codon in codons:
                dna_codon = codon.replace('U', 'T')  # RNA -> DNA
                weight = wi_table.get(dna_codon, 0.1)
                weights.append(weight)
                
            return weights
        else:
            # Fallback: uniform weights (no CAI optimization)
            return [1.0] * len(codons)
        
    def enhance_codon_probs(self, 
                           codon_probs: torch.Tensor,
                           amino_sequence: str,
                           valid_codon_mask: Optional[torch.Tensor] = None,
                           codon_indices: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, Dict]:
        """

        """
        metadata = {
            'cai_enabled': self.enable_cai,
            'method': 'training_passthrough',
            'note': 'CAI只用于离散化选择，不影响训练'
        }
        return codon_probs, metadata
    
    def enhance_rna_probs(self,
                         rna_probs: torch.Tensor,
                         amino_sequence: str) -> Tuple[torch.Tensor, Dict]:
        """
        Enhance RNA probabilities with CAI optimization using simple approach (for Lagrangian)
        
        Args:
            rna_probs: RNA probability distributions [batch, seq_len, 4] or [seq_len, 4]
            amino_sequence: Target amino acid sequence
            
        Returns:
            enhanced_rna: CAI-weighted RNA probabilities  
            metadata: Optimization details
        """
        if not self.enable_cai or self.cai_weights is None:
            return rna_probs, {'cai_enabled': False}
            
        # For Lagrangian constraints, pass through without modification
        # Lagrangian uses penalty-based approach, so simple reweighting may not be appropriate
        metadata = {
            'cai_enabled': True,
            'method': 'passthrough_lagrangian',
            'note': 'Lagrangian uses penalty-based CAI constraints'
        }
        
        return rna_probs, metadata
    
    def apply_cai_to_discrete_selection(self,
                                       probs: torch.Tensor,
                                       amino_sequence: str,
                                       valid_codon_mask: Optional[torch.Tensor] = None,
                                       codon_indices: Optional[torch.Tensor] = None,
                                       cai_factor: float = 1.0) -> torch.Tensor:
        """

        
        Args:





            
        Returns:
            interpolated_probs: cai_factor*w + (1-cai_factor)*π
        """
        if not self.enable_cai or self.cai_weights is None:
            return probs
            
        with torch.no_grad():
            if valid_codon_mask is not None and codon_indices is not None:

                pi_probs = probs.clone()
                

                w_probs = torch.zeros_like(probs)
                device = probs.device
                valid_codon_mask = valid_codon_mask.to(device)
                codon_indices = codon_indices.to(device)
                cai_weights_tensor = torch.tensor(self.cai_weights, dtype=torch.float32, device=device)
                
                num_positions = probs.shape[0]
                for pos in range(num_positions):

                    valid_mask = valid_codon_mask[pos]
                    if not valid_mask.any():
                        continue
                        

                    for codon_slot in range(probs.shape[1]):
                        if valid_mask[codon_slot]:
                            global_codon_idx = codon_indices[pos, codon_slot].item()
                            if global_codon_idx < len(self.cai_weights):
                                cai_weight = cai_weights_tensor[global_codon_idx]
                                w_probs[pos, codon_slot] = cai_weight
                    


                    w_logits = torch.zeros_like(w_probs[pos])
                    

                    for codon_slot in range(probs.shape[1]):
                        if valid_mask[codon_slot]:
                            weight = w_probs[pos, codon_slot].item()
                            if weight > 0:


                                w_logits[codon_slot] = torch.log(torch.tensor(weight))
                            else:

                                w_logits[codon_slot] = -10.0
                    

                    w_logits[~valid_mask] = -float('inf')
                    

                    w_probs[pos] = torch.softmax(w_logits, dim=-1)
                



                interpolated_probs = cai_factor * w_probs + (1 - cai_factor) * pi_probs
                

                for pos in range(num_positions):
                    valid_mask = valid_codon_mask[pos]
                    if valid_mask.sum() > 0:

                        valid_sum = interpolated_probs[pos, valid_mask].sum()
                        if abs(valid_sum - 1.0) > 1e-6:
                            interpolated_probs[pos, valid_mask] /= valid_sum
                            
                return interpolated_probs
            else:

                return probs
    
    def binary_search_optimal_cai_factor(self,
                                        pi_probs: torch.Tensor,
                                        amino_sequence: str,
                                        valid_codon_mask: torch.Tensor,
                                        codon_indices: torch.Tensor,
                                        target_cai: float = 0.8,
                                        max_iterations: int = 50,
                                        verbose: bool = True,
                                        use_incremental: bool = True) -> float:
        """


        

        
        Args:







            
        Returns:

        """
        if not self.enable_cai or self.cai_weights is None:
            return 0.0
            
        with torch.no_grad():
            if verbose:
                print(f"🔍 CAI二分搜索开始 {'(增量计算)' if use_incremental else ''} - 序列:{amino_sequence}, 目标CAI≥{target_cai}")
            
            # Initialize incremental calculator if enabled
            incremental_calc = None
            if use_incremental:
                incremental_calc = CachedIncrementalCalculator(self.cai_weights, str(self.device))
            


            

            pi_probs = self._apply_symmetry_breaking_perturbation(
                pi_probs, amino_sequence, valid_codon_mask, codon_indices, verbose
            )
                


            for pos in range(pi_probs.shape[0]):
                valid_mask = valid_codon_mask[pos]
                if valid_mask.any():

                    valid_sum = pi_probs[pos, valid_mask].sum()
                    if abs(valid_sum - 1.0) > 1e-5:
                        if verbose:
                            print(f"⚠️ 位置{pos}的概率和={valid_sum:.6f}，需要归一化")
                        pi_probs[pos, valid_mask] = pi_probs[pos, valid_mask] / valid_sum
                
            # Import CAI validator if needed  
            from id3.cai.validator import compute_cai_from_sequence
            
            # Pre-compute pi and w best codons if using incremental
            pi_best_codons = {}
            w_best_codons = {}
            initial_sequence = []
            
            if use_incremental:
                num_positions = pi_probs.shape[0]
                for pos in range(num_positions):
                    valid_mask = valid_codon_mask[pos]
                    if not valid_mask.any():
                        continue
                    
                    # Get π-optimal codon
                    valid_probs = pi_probs[pos, valid_mask]
                    pi_best_local = torch.argmax(valid_probs).item()
                    valid_slots = torch.where(valid_mask)[0]
                    pi_best_global = codon_indices[pos, valid_slots[pi_best_local]].item()
                    pi_best_codons[pos] = pi_best_global
                    initial_sequence.append(pi_best_global)
                    
                    # Get w-optimal codon
                    w_best_local = -1
                    max_weight = -1
                    for local_idx, slot in enumerate(valid_slots.cpu().numpy()):
                        global_idx = codon_indices[pos, slot].item()
                        if global_idx < len(self.cai_weights):
                            if self.cai_weights[global_idx] > max_weight:
                                max_weight = self.cai_weights[global_idx]
                                w_best_local = local_idx
                    
                    if w_best_local >= 0:
                        w_best_global = codon_indices[pos, valid_slots[w_best_local]].item()
                        w_best_codons[pos] = w_best_global
                    else:
                        w_best_codons[pos] = pi_best_global
                
                # Initialize incremental calculator
                incremental_calc.initialize(initial_sequence)
            

            cai_factor_low = 0.0
            if use_incremental:
                # Compute switching events for incremental update
                from id3.cai.discrete_binary_search import DiscreteCAISearcher
                searcher = DiscreteCAISearcher(self.cai_weights, str(self.device), use_incremental=False)
                events = searcher.compute_switching_events(pi_probs, valid_codon_mask, codon_indices, amino_sequence)
                cai_low, _ = incremental_calc.get_cai_at_cai_factor(cai_factor_low, events, pi_best_codons, w_best_codons)
            else:
                probs_low = self.apply_cai_to_discrete_selection(
                    pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=cai_factor_low
                )

                discrete_seq_low = self._get_discrete_sequence_from_probs(probs_low, valid_codon_mask, codon_indices)
                from id3.cai.validator import compute_cai_from_sequence
                cai_low = compute_cai_from_sequence(discrete_seq_low)
            
            if verbose:
                print(f"📊 边界检查: CAI插值因子=0.0 (纯accessibility) -> CAI={cai_low:.4f}")
            

            if cai_low >= target_cai:
                if verbose:
                    print(f"✅ 纯accessibility已满足要求CAI={cai_low:.4f}≥{target_cai}")
                return cai_factor_low
            

            cai_factor_high = 1.0
            if use_incremental:
                cai_high, _ = incremental_calc.get_cai_at_cai_factor(cai_factor_high, events, pi_best_codons, w_best_codons)
            else:
                probs_high = self.apply_cai_to_discrete_selection(
                    pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=cai_factor_high
                )

                discrete_seq_high = self._get_discrete_sequence_from_probs(probs_high, valid_codon_mask, codon_indices)
                from id3.cai.validator import compute_cai_from_sequence
                cai_high = compute_cai_from_sequence(discrete_seq_high)
            
            if verbose:
                print(f"📊 边界检查: CAI插值因子=1.0 (纯CAI) -> CAI={cai_high:.4f}")
            

            if cai_high < target_cai:
                if verbose:
                    print(f"⚠️ 警告: 该序列理论最大CAI={cai_high:.4f} < 目标{target_cai}")
                    print(f"         无法满足CAI≥{target_cai}约束，返回CAI插值因子=1.0（最佳努力）")

                return cai_factor_high
            
            if verbose:
                print(f"🎯 开始二分搜索: [{cai_factor_low:.6f}, {cai_factor_high:.6f}]")
            

            best_cai_factor = cai_factor_high
            
            for iteration in range(max_iterations):
                cai_factor_mid = (cai_factor_low + cai_factor_high) / 2.0
                
                if use_incremental:
                    # Use incremental calculation with caching
                    cai_mid, cache_hit = incremental_calc.get_cai_at_cai_factor(
                        cai_factor_mid, events, pi_best_codons, w_best_codons
                    )
                else:

                    probs_mid = self.apply_cai_to_discrete_selection(
                        pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=cai_factor_mid
                    )
                    

                    discrete_sequence = self._get_discrete_sequence_from_probs(
                        probs_mid, valid_codon_mask, codon_indices
                    )
                    

                    from id3.cai.validator import compute_cai_from_sequence
                    cai_mid = compute_cai_from_sequence(discrete_sequence)
                
                if verbose:
                    status = "✅" if cai_mid >= target_cai else "❌"
                    action = "缩小CAI插值因子(保持accessibility)" if cai_mid >= target_cai else "增大CAI插值因子(提升CAI)"
                    print(f"  Iter{iteration:2d}: CAI插值因子={cai_factor_mid:.6f} -> CAI={cai_mid:.4f} {status} [{action}]")
                
                if cai_mid >= target_cai:

                    best_cai_factor = cai_factor_mid


                    if cai_mid >= target_cai + 0.05:
                        cai_factor_high = cai_factor_mid
                    else:

                        break
                else:

                    cai_factor_low = cai_factor_mid
                

                if abs(cai_factor_high - cai_factor_low) < 1e-6:

                    if best_cai_factor == cai_factor_mid and cai_mid < target_cai:

                        best_cai_factor = cai_factor_high
                    if verbose:
                        print(f"🎯 收敛: 区间宽度 < 1e-6, 最终CAI插值因子={best_cai_factor:.6f}")
                    break
            
            if verbose:
                print(f"🏁 搜索完成: 最优CAI插值因子={best_cai_factor:.6f}")
                
                # Print incremental calculator statistics if used
                if use_incremental and incremental_calc:
                    stats = incremental_calc.get_statistics()
                    cache_stats = incremental_calc.get_cache_statistics()
                    print(f"📊 增量计算统计:")
                    print(f"   效率比: {stats['efficiency_ratio']}")
                    print(f"   缓存命中率: {cache_stats['hit_rate']}")
                

                final_probs = self.apply_cai_to_discrete_selection(
                    pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=best_cai_factor
                )
                final_discrete_seq = self._get_discrete_sequence_from_probs(final_probs, valid_codon_mask, codon_indices)
                final_cai = compute_cai_from_sequence(final_discrete_seq)
                print(f"✅ 最终验证: CAI插值因子={best_cai_factor:.6f} -> 离散CAI={final_cai:.6f} (目标≥{target_cai})")
                

                if final_cai < target_cai:
                    print(f"⚠️ 警告：最终CAI {final_cai:.6f} < 目标 {target_cai}，可能是边界情况")
            

            final_probs = self.apply_cai_to_discrete_selection(
                pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=best_cai_factor
            )
            final_discrete_seq = self._get_discrete_sequence_from_probs(final_probs, valid_codon_mask, codon_indices)
            final_cai = compute_cai_from_sequence(final_discrete_seq)
            

            if final_cai < target_cai:
                if verbose:
                    print(f"⚠️ 需要调整: CAI={final_cai:.6f} < {target_cai}")
                

                for adjustment in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]:
                    test_cai_factor = min(best_cai_factor + adjustment, 1.0)
                    test_probs = self.apply_cai_to_discrete_selection(
                        pi_probs, amino_sequence, valid_codon_mask, codon_indices, cai_factor=test_cai_factor
                    )
                    test_seq = self._get_discrete_sequence_from_probs(test_probs, valid_codon_mask, codon_indices)
                    test_cai = compute_cai_from_sequence(test_seq)
                    
                    if test_cai >= target_cai:
                        best_cai_factor = test_cai_factor
                        if verbose:
                            print(f"✅ 调整成功: CAI插值因子={test_cai_factor:.6f} -> CAI={test_cai:.6f}")
                        break
                    elif test_cai_factor >= 1.0:

                        if verbose:
                            print(f"⚠️ 已达到CAI插值因子=1.0上限，CAI={test_cai:.6f}")
                        best_cai_factor = 1.0
                        break
            
            return best_cai_factor
    
    def _apply_symmetry_breaking_perturbation(self,
                                              pi_probs: torch.Tensor,
                                              amino_sequence: str,
                                              valid_codon_mask: torch.Tensor,
                                              codon_indices: torch.Tensor,
                                              verbose: bool = False) -> torch.Tensor:
        """

        


        

        
        Args:





            
        Returns:

        """
        with torch.no_grad():
            device = pi_probs.device
            num_positions = pi_probs.shape[0]
            

            switching_points = []
            for pos in range(num_positions):
                valid_mask = valid_codon_mask[pos].to(device)
                if not valid_mask.any():
                    switching_points.append(float('inf'))
                    continue
                

                pi_pos_probs = pi_probs[pos, valid_mask]
                pi_best_idx_local = torch.argmax(pi_pos_probs).item()
                

                w_best_idx_local = -1
                max_cai_weight = -1
                for local_idx, slot in enumerate(torch.where(valid_mask)[0].cpu().numpy()):
                    global_codon_idx = codon_indices[pos, slot].item()
                    if global_codon_idx < len(self.cai_weights):
                        if self.cai_weights[global_codon_idx] > max_cai_weight:
                            max_cai_weight = self.cai_weights[global_codon_idx]
                            w_best_idx_local = local_idx
                
                if pi_best_idx_local != w_best_idx_local and w_best_idx_local >= 0:


                    p_pi = pi_pos_probs[pi_best_idx_local].item()
                    p_w = pi_pos_probs[w_best_idx_local].item() if w_best_idx_local < len(pi_pos_probs) else 0
                    

                    alpha_star = p_pi - p_w if p_pi > p_w else 0
                    switching_points.append(alpha_star)
                else:

                    switching_points.append(float('inf'))
            

            from collections import defaultdict
            switch_groups = defaultdict(list)
            for pos, alpha_star in enumerate(switching_points):
                if alpha_star != float('inf'):

                    key = round(alpha_star, 8)
                    switch_groups[key].append(pos)
            

            perturbed_pi_probs = pi_probs.clone()
            epsilon_base = 1e-10
            
            for alpha_key, positions in switch_groups.items():
                if len(positions) > 1:
                    if verbose:
                        print(f"🔧 检测到{len(positions)}个位置在α*≈{alpha_key:.6f}处同时切换")
                        print(f"   位置: {positions[:5]}{'...' if len(positions) > 5 else ''}")
                    

                    for idx, pos in enumerate(positions):
                        valid_mask = valid_codon_mask[pos].to(device)
                        if not valid_mask.any():
                            continue
                        

                        epsilon = epsilon_base * (1 + idx)
                        

                        valid_probs = perturbed_pi_probs[pos, valid_mask]
                        

                        perturbation = torch.randn_like(valid_probs) * epsilon
                        perturbed_probs = valid_probs + perturbation
                        

                        perturbed_probs = torch.clamp(perturbed_probs, min=1e-12)
                        perturbed_probs = perturbed_probs / perturbed_probs.sum()
                        
                        perturbed_pi_probs[pos, valid_mask] = perturbed_probs
            
            if verbose and len(switch_groups) > 0:
                total_perturbed = sum(len(positions) for positions in switch_groups.values() if len(positions) > 1)
                if total_perturbed > 0:
                    print(f"✅ 已对{total_perturbed}个位置添加对称性破坏扰动")
            
            return perturbed_pi_probs
    
    def _compute_expected_cai(self, 
                             probs: torch.Tensor,
                             valid_codon_mask: torch.Tensor,
                             codon_indices: torch.Tensor) -> float:
        """

        
        Args:



            
        Returns:

        """
        if not self.enable_cai or self.cai_weights is None:
            return 0.0
            
        device = probs.device
        cai_weights_tensor = torch.tensor(self.cai_weights, dtype=torch.float32, device=device)
        
        total_cai_log = 0.0
        valid_positions = 0
        
        num_positions = probs.shape[0]
        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos].to(device)
            if not valid_mask.any():
                continue
                

            expected_weight = 0.0
            for slot in range(probs.shape[1]):
                if valid_mask[slot]:
                    codon_idx = codon_indices[pos, slot].item()
                    if 0 <= codon_idx < len(self.cai_weights):
                        prob = probs[pos, slot].item()
                        weight = cai_weights_tensor[codon_idx].item()
                        expected_weight += prob * weight
            
            if expected_weight > 1e-8:
                total_cai_log += torch.log(torch.tensor(expected_weight)).item()
                valid_positions += 1
        
        if valid_positions == 0:
            return 0.0
            

        expected_cai = torch.exp(torch.tensor(total_cai_log / valid_positions)).item()
        return expected_cai
    
    def _get_discrete_sequence_from_probs(self, 
                                         probs: torch.Tensor,
                                         valid_codon_mask: torch.Tensor,
                                         codon_indices: torch.Tensor) -> str:
        """

        
        Args:



            
        Returns:

        """
        from id3.utils.constants import codon_to_rna_matrix
        
        device = probs.device
        num_positions = probs.shape[0]
        selected_codons = []
        

        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos].to(device)
            if not valid_mask.any():
                continue
            

            valid_probs = probs[pos, valid_mask]
            if valid_probs.sum() < 1e-8:

                selected_idx = 0
            else:

                selected_idx = torch.argmax(valid_probs).item()
            

            valid_indices = torch.nonzero(valid_mask).flatten()
            actual_codon_slot = valid_indices[selected_idx].item()
            

            global_codon_idx = codon_indices[pos, actual_codon_slot].item()
            selected_codons.append(global_codon_idx)
        

        rna_sequence_parts = []
        for codon_idx in selected_codons:
            codon_rna = codon_to_rna_matrix[codon_idx]  # [3, 4] one-hot
            for pos in range(3):
                nucleotide_idx = torch.argmax(codon_rna[pos]).item()
                nucleotides = ['A', 'C', 'G', 'U']
                rna_sequence_parts.append(nucleotides[nucleotide_idx])
        
        return ''.join(rna_sequence_parts)
    
    def binary_search_discrete(self,
                              pi_probs: torch.Tensor,
                              amino_sequence: str,
                              valid_codon_mask: torch.Tensor,
                              codon_indices: torch.Tensor,
                              target_cai: float = 0.8,
                              verbose: bool = True) -> Tuple[float, Dict]:
        """
        Discrete binary search for CAI optimization.
        
        This is a more efficient alternative to binary_search_optimal_cai_factor()
        that operates in discrete switching event space rather than continuous

        
        Args:
            pi_probs: π probability distribution (accessibility-optimal)
            amino_sequence: Amino acid sequence
            valid_codon_mask: Valid codon mask
            codon_indices: Codon indices
            target_cai: Target CAI value
            verbose: Whether to print debug info
            
        Returns:
            optimal_alpha: Optimal α value
            metadata: Search metadata including performance metrics
        """
        if not self.enable_cai or self.cai_weights is None:
            return 0.0, {'method': 'disabled'}
        
        if not hasattr(self, 'discrete_searcher'):
            # Fallback to continuous search if discrete searcher not available
            if verbose:
                print("⚠️ Discrete searcher not initialized, using continuous search")
            return self.binary_search_optimal_cai_factor(
                pi_probs, amino_sequence, valid_codon_mask, 
                codon_indices, target_cai, verbose=verbose
            ), {'method': 'continuous_fallback'}
        
        # Use discrete binary search
        optimal_alpha, metadata = self.discrete_searcher.discrete_binary_search(
            pi_probs, valid_codon_mask, codon_indices,
            amino_sequence, target_cai, verbose
        )
        
        metadata['method'] = 'discrete'
        return optimal_alpha, metadata
    
    def binary_search_optimal_alpha(self, *args, **kwargs):
        """
        Alias for binary_search_optimal_cai_factor for compatibility.
        This method maintains backward compatibility with existing constraint code.
        """
        return self.binary_search_optimal_cai_factor(*args, **kwargs)
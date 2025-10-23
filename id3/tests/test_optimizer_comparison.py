"""



"""

import sys
import torch
import numpy as np
import time
import logging
from typing import Dict, List, Tuple


sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class OptimizerComparison:

    
    def __init__(self, test_sequence: str = None):
        """
        初始化测试
        
        Args:
            test_sequence: 测试用的氨基酸序列
        """

        self.sequence = test_sequence or "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
        self.target_cai = 0.8
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        logger.info(f"Device: {self.device}")
        logger.info(f"Sequence length: {len(self.sequence)}")
    
    def _create_test_distribution(self) -> Tuple[torch.Tensor, torch.Tensor]:
        """创建测试用的概率分布和掩码"""
        seq_len = len(self.sequence)
        max_codons = 6
        

        distribution = torch.rand(seq_len, max_codons, device=self.device)
        

        valid_mask = torch.zeros(seq_len, max_codons, device=self.device)
        
        for pos, aa in enumerate(self.sequence):
            if aa in amino_acids_to_codons:
                num_codons = len(amino_acids_to_codons[aa])
                valid_mask[pos, :num_codons] = 1.0

                distribution[pos] = distribution[pos] * valid_mask[pos]
                if distribution[pos].sum() > 0:
                    distribution[pos] = distribution[pos] / distribution[pos].sum()
        
        return distribution, valid_mask
    
    def test_single_optimization(self, method: str) -> Dict:
        """

        
        Args:

            
        Returns:

        """
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing method: {method}")
        logger.info(f"{'='*50}")
        

        operator = CAIEnhancementOperator(
            method=method,
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=self.sequence
        )
        

        distribution, valid_mask = self._create_test_distribution()
        codon_indices = torch.zeros_like(distribution, dtype=torch.long)
        

        start_time = time.time()
        result_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=distribution,
            amino_acid_sequence=self.sequence,
            valid_codon_mask=valid_mask,
            codon_indices=codon_indices,
            target_cai=self.target_cai
        )
        single_time = time.time() - start_time
        

        num_iterations = 10
        start_time = time.time()
        
        for i in range(num_iterations):

            test_dist = distribution + torch.rand_like(distribution) * 0.1
            test_dist = test_dist * valid_mask
            test_dist = test_dist / (test_dist.sum(dim=-1, keepdim=True) + 1e-10)
            
            result_dist, metadata = operator.apply_cai_enhancement(
                pi_accessibility=test_dist,
                amino_acid_sequence=self.sequence,
                valid_codon_mask=valid_mask,
                codon_indices=codon_indices,
                target_cai=self.target_cai
            )
        
        total_time = time.time() - start_time
        avg_time = total_time / num_iterations
        

        stats = operator.get_statistics()
        

        result = {
            'method': method,
            'single_optimization_time': single_time * 1000,
            'average_time_per_iteration': avg_time * 1000,
            'total_time': total_time,
            'num_iterations': num_iterations,
            'final_cai': metadata.get('final_cai', 0.0),
            'constraint_satisfied': metadata.get('constraint_satisfied', False),
            'statistics': stats
        }
        

        logger.info(f"Single optimization: {result['single_optimization_time']:.2f} ms")
        logger.info(f"Average per iteration: {result['average_time_per_iteration']:.2f} ms")
        logger.info(f"Final CAI: {result['final_cai']:.4f}")
        logger.info(f"Constraint satisfied: {result['constraint_satisfied']}")
        

        if 'diversity' in stats:
            diversity = stats['diversity']
            logger.info(f"Unique sequences: {diversity.get('num_unique', 0)}/{diversity.get('num_sequences', 0)}")
            logger.info(f"Repetition rate: {diversity.get('repetition_rate', 0):.1%}")
        
        return result
    
    def compare_methods(self) -> Dict:
        """

        
        Returns:

        """
        logger.info("\n" + "="*60)
        logger.info("CAI OPTIMIZER PERFORMANCE COMPARISON")
        logger.info("="*60)
        
        methods = CAIEnhancementOperator.available_methods()
        logger.info(f"Available methods: {methods}")
        
        results = {}
        

        for method in methods:
            try:
                result = self.test_single_optimization(method)
                results[method] = result
            except Exception as e:
                logger.error(f"Error testing {method}: {e}")
                results[method] = {'error': str(e)}
        

        if 'binary_search' in results and 'sado' in results:
            binary_time = results['binary_search']['average_time_per_iteration']
            sado_time = results['sado']['average_time_per_iteration']
            
            if sado_time > 0:
                speedup = binary_time / sado_time
                logger.info(f"\n{'='*60}")
                logger.info("PERFORMANCE COMPARISON")
                logger.info(f"{'='*60}")
                logger.info(f"Binary Search: {binary_time:.2f} ms/iteration")
                logger.info(f"SADO: {sado_time:.2f} ms/iteration")
                logger.info(f"SPEEDUP: {speedup:.1f}x faster")
                
                results['comparison'] = {
                    'speedup': speedup,
                    'binary_search_time': binary_time,
                    'sado_time': sado_time
                }
        
        return results
    
    def test_repetition_rate(self, method: str = 'sado', num_iterations: int = 100) -> float:
        """

        
        Args:


            
        Returns:

        """
        logger.info(f"\n{'='*50}")
        logger.info(f"Testing repetition rate for {method}")
        logger.info(f"{'='*50}")
        
        operator = CAIEnhancementOperator(
            method=method,
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=self.sequence
        )
        
        distribution, valid_mask = self._create_test_distribution()
        codon_indices = torch.zeros_like(distribution, dtype=torch.long)
        

        for i in range(num_iterations):

            test_dist = distribution + torch.rand_like(distribution) * 0.05
            test_dist = test_dist * valid_mask
            test_dist = test_dist / (test_dist.sum(dim=-1, keepdim=True) + 1e-10)
            
            operator.apply_cai_enhancement(
                pi_accessibility=test_dist,
                amino_acid_sequence=self.sequence,
                valid_codon_mask=valid_mask,
                codon_indices=codon_indices,
                target_cai=self.target_cai
            )
            
            if (i + 1) % 20 == 0:
                logger.info(f"Completed {i + 1}/{num_iterations} iterations")
        

        stats = operator.get_statistics()
        
        if method == 'sado' and 'diversity' in stats:
            diversity = stats['diversity']
            repetition_rate = diversity.get('repetition_rate', 0.0)
            logger.info(f"\nResults after {num_iterations} iterations:")
            logger.info(f"Unique sequences: {diversity.get('num_unique', 0)}")
            logger.info(f"Total sequences: {diversity.get('num_sequences', 0)}")
            logger.info(f"Repetition rate: {repetition_rate:.2%}")
            return repetition_rate
        
        return 0.0


def main():


    tester = OptimizerComparison()
    

    logger.info("\n" + "="*60)
    logger.info("STARTING PERFORMANCE COMPARISON TEST")
    logger.info("="*60)
    
    comparison_results = tester.compare_methods()
    

    logger.info("\n" + "="*60)
    logger.info("TESTING SADO REPETITION RATE")
    logger.info("="*60)
    
    repetition_rate = tester.test_repetition_rate('sado', num_iterations=100)
    

    logger.info("\n" + "="*60)
    logger.info("TEST SUMMARY")
    logger.info("="*60)
    
    if 'comparison' in comparison_results:
        comp = comparison_results['comparison']
        logger.info(f"✅ SADO is {comp['speedup']:.1f}x faster than Binary Search")
        logger.info(f"✅ SADO repetition rate: {repetition_rate:.2%}")
        

        if comp['speedup'] >= 50:
            logger.info("🎉 SADO achieves expected 50x+ speedup!")
        else:
            logger.info(f"⚠️  SADO speedup ({comp['speedup']:.1f}x) is less than expected (66x)")
        

            logger.info("🎉 SADO achieves near-zero repetition rate!")
        else:
            logger.info(f"⚠️  SADO repetition rate ({repetition_rate:.2%}) is higher than expected")
    
    logger.info("\n" + "="*60)
    logger.info("ALL TESTS COMPLETED")
    logger.info("="*60)


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""







"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import json
import torch
import numpy as np
from typing import Dict, List, Any

from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging

logger = setup_logging(level='INFO', name='sado_quick')

class QuickSADOTester:
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.sequences = self.load_test_sequences()
        logger.info(f"设备: {self.device}")
    
    def load_test_sequences(self) -> Dict[str, str]:

        try:
            with open('data/codon_references/ecoli_bl21de3_reference_sequences.json', 'r') as f:
                data = json.load(f)
            

            sequences = data.get('sequences', {})
            

            codon_table = {
                'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
                'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
                'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
            }
            
            amino_sequences = {}
            for name, dna_seq in sequences.items():
                rna_seq = dna_seq.replace('T', 'U')
                amino_seq = ""
                
                for i in range(0, len(rna_seq) - 2, 3):
                    codon = rna_seq[i:i+3]
                    if codon in codon_table:
                        aa = codon_table[codon]

                            break
                        amino_seq += aa
                

                    amino_sequences[name] = amino_seq
            

            return amino_sequences
            
        except Exception as e:

            return {}
    
    def compute_sequence_probability(self, indices, pi_distribution, optimizer):
        """计算序列的概率"""
        prob = 1.0
        
        for pos, idx in enumerate(indices):
            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):
                choice = optimizer.codon_choices[pos][idx]
                orig_idx = choice.get('original_local_index', idx)
                
                if orig_idx < pi_distribution.shape[1]:
                    p = pi_distribution[pos, orig_idx].item()
                    prob *= p
        
        return prob
    
    def test_speed_benchmark(self):

        logger.info("=" * 60)

        logger.info("=" * 60)
        
        test_cases = [
            ('short', 30),
            ('medium', 100),
            ('long', 200)
        ]
        
        results = {}
        
        for case_name, seq_length in test_cases:

            logger.info("-" * 40)
            

            amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
            test_seq = ''.join(np.random.choice(amino_acids, seq_length))
            

            optimizer = SADOOptimizer(
                species='ecoli_bl21de3',
                device=self.device,
                amino_acid_sequence=test_seq
            )
            

            num_codons = 6
            valid_mask = torch.zeros(seq_length, num_codons, dtype=torch.bool, device=self.device)
            
            for pos, aa in enumerate(test_seq):
                if aa in amino_acids_to_codons:
                    num_valid = len(amino_acids_to_codons[aa])
                    valid_mask[pos, :min(num_valid, num_codons)] = True
            

            pi_test = torch.rand(seq_length, num_codons, device=self.device)
            pi_test = pi_test * valid_mask.float()
            for pos in range(seq_length):
                if valid_mask[pos].any():
                    pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
            

            times = []
            cais = []
            probs = []
            satisfied = []
            
            for i in range(3):
                optimizer.reset()
                start_time = time.time()
                
                _, metadata = optimizer.optimize(
                    pi_accessibility=pi_test,
                    target_cai=0.8,
                    amino_acid_sequence=test_seq
                )
                
                end_time = time.time()
                execution_time = (end_time - start_time) * 1000  # ms
                
                times.append(execution_time)
                cais.append(metadata['final_cai'])
                satisfied.append(metadata['constraint_satisfied'])
                

                if optimizer.last_indices is not None:
                    prob = self.compute_sequence_probability(
                        optimizer.last_indices, pi_test, optimizer
                    )
                    probs.append(np.log(prob) if prob > 0 else float('-inf'))
                else:
                    probs.append(float('-inf'))
            

            avg_time = np.mean(times)
            std_time = np.std(times)
            avg_cai = np.mean(cais)
            avg_prob = np.mean(probs)
            satisfaction_rate = f"{sum(satisfied)}/{len(satisfied)}"
            




            
            results[case_name] = {
                'length': seq_length,
                'avg_time_ms': avg_time,
                'std_time_ms': std_time,
                'avg_cai': avg_cai,
                'avg_log_prob': avg_prob,
                'satisfaction_rate': satisfaction_rate
            }
        
        return results
    
    def test_real_sequences(self):
        """真实序列测试"""
        logger.info("=" * 60)
        logger.info("真实序列性能测试")
        logger.info("=" * 60)
        

        test_seqs = {}
        for name, seq in self.sequences.items():
            if 50 <= len(seq) <= 200:
                test_seqs[name] = seq
                if len(test_seqs) >= 3:
                    break
        
        results = {}
        
        for name, amino_seq in test_seqs.items():
            logger.info(f"\n测试序列: {name} (长度: {len(amino_seq)})")
            logger.info("-" * 40)
            

            optimizer = SADOOptimizer(
                species='ecoli_bl21de3',
                device=self.device,
                amino_acid_sequence=amino_seq
            )
            

            seq_len = len(amino_seq)
            num_codons = 6
            valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=self.device)
            
            for pos, aa in enumerate(amino_seq):
                if aa in amino_acids_to_codons:
                    num_valid = len(amino_acids_to_codons[aa])
                    valid_mask[pos, :min(num_valid, num_codons)] = True
            

            torch.manual_seed(42)
            pi_test = torch.rand(seq_len, num_codons, device=self.device)
            pi_test = pi_test * valid_mask.float()
            for pos in range(seq_len):
                if valid_mask[pos].any():
                    pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
            

            start_time = time.time()
            _, metadata = optimizer.optimize(
                pi_accessibility=pi_test,
                target_cai=0.8,
                amino_acid_sequence=amino_seq
            )
            end_time = time.time()
            
            execution_time = (end_time - start_time) * 1000  # ms
            

            log_prob = float('-inf')
            if optimizer.last_indices is not None:
                prob = self.compute_sequence_probability(
                    optimizer.last_indices, pi_test, optimizer
                )
                log_prob = np.log(prob) if prob > 0 else float('-inf')
            
            logger.info(f"执行时间: {execution_time:.2f} ms")
            logger.info(f"最终CAI: {metadata['final_cai']:.4f}")
            logger.info(f"约束满足: {'✅' if metadata['constraint_satisfied'] else '❌'}")
            logger.info(f"log(P): {log_prob:.2f}")
            logger.info(f"唯一序列数: {metadata['unique_sequences']}")
            
            results[name] = {
                'length': len(amino_seq),
                'execution_time_ms': execution_time,
                'final_cai': metadata['final_cai'],
                'constraint_satisfied': metadata['constraint_satisfied'],
                'log_prob': log_prob,
                'unique_sequences': metadata['unique_sequences']
            }
        
        return results
    
    def run_all_tests(self):




        

        benchmark_results = self.test_speed_benchmark()
        

        real_seq_results = self.test_real_sequences()
        

        logger.info("=" * 60)

        logger.info("=" * 60)
        

        for case, result in benchmark_results.items():
            logger.info(f"  {case}: {result['avg_time_ms']:.1f}ms, "
                       f"CAI={result['avg_cai']:.3f}, "

        

        for name, result in real_seq_results.items():
            logger.info(f"  {name}: {result['execution_time_ms']:.1f}ms, "
                       f"CAI={result['final_cai']:.3f}, "
                       f"{'✅' if result['constraint_satisfied'] else '❌'}")
        


        short_time = benchmark_results['short']['avg_time_ms']
        medium_time = benchmark_results['medium']['avg_time_ms']
        long_time = benchmark_results['long']['avg_time_ms']
        



        

        total_satisfied = sum(1 for r in real_seq_results.values() if r['constraint_satisfied'])
        total_tests = len(real_seq_results)
        cai_rate = total_satisfied / total_tests if total_tests > 0 else 0
        


        
        return {
            'benchmark_results': benchmark_results,
            'real_seq_results': real_seq_results,
            'summary': {
                'cai_satisfaction_rate': cai_rate,
                'avg_short_time': short_time,
                'avg_medium_time': medium_time,
                'avg_long_time': long_time
            }
        }

def main():
    """主函数"""
    tester = QuickSADOTester()
    results = tester.run_all_tests()
    logger.info("\n✅ 测试完成!")
    
    return results

if __name__ == "__main__":
    main()
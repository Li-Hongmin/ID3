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
import pandas as pd
from typing import Dict, List, Any
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
from collections import defaultdict

from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='sado_performance')

@dataclass
class PerformanceResult:

    optimizer_name: str
    sequence_name: str
    sequence_length: int
    execution_time_ms: float
    final_cai: float
    constraint_satisfied: bool
    probability_score: float  # log(P)
    unique_sequences: int
    collision_detected: bool
    iterations: int = 0
    memory_usage_mb: float = 0.0

class SADOPerformanceTester:
    """SADO性能测试器"""
    
    def __init__(self, device: str = 'auto'):
        self.device = self._setup_device(device)
        self.results: List[PerformanceResult] = []
        self.sequence_data = self._load_sequence_data()
        
    def _setup_device(self, device: str) -> torch.device:

        if device == 'auto':
            return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        return torch.device(device)
    
    def _load_sequence_data(self) -> Dict[str, str]:
        """加载E.coli序列数据"""
        data_path = Path(__file__).parent.parent.parent / "data" / "codon_references" / "ecoli_bl21de3_reference_sequences.json"
        
        try:
            with open(data_path, 'r') as f:
                data = json.load(f)
                sequences = data.get('sequences', {})
                

            amino_sequences = {}
            genetic_code = {
                'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
                'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
                'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
            }
            
            for name, dna_seq in sequences.items():

                rna_seq = dna_seq.replace('T', 'U')
                amino_seq = ""
                
                for i in range(0, len(rna_seq) - 2, 3):
                    codon = rna_seq[i:i+3].replace('U', 'T')
                    if len(codon) == 3 and codon in genetic_code:
                        aa = genetic_code[codon]
                        if aa == '*':
                            break
                        amino_seq += aa
                
                if amino_seq and len(amino_seq) >= 20:
                    amino_sequences[name] = amino_seq
                    
            logger.info(f"加载了{len(amino_sequences)}个有效的氨基酸序列")
            return amino_sequences
            
        except Exception as e:
            logger.warning(f"无法加载序列数据: {e}")

            return {
                'test_short': 'MSKGEELFTGVVPILVELD',
                'test_medium': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKLICTTGKLPVPWPTLVTTLGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYISHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK',
                'test_long': 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKLICTTGKLPVPWPTLVTTLGYGLQCFARYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYISHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKGSHHHHHH' * 2
            }
    
    def _create_test_distribution(self, seq_length: int, num_codons: int = 6, seed: int = 42) -> torch.Tensor:

        torch.manual_seed(seed)
        

        valid_mask = torch.zeros(seq_length, num_codons, dtype=torch.bool, device=self.device)
        


        

        pi = torch.rand(seq_length, num_codons, device=self.device)
        pi = pi * valid_mask.float()
        

        for pos in range(seq_length):
            if valid_mask[pos].any():
                pi[pos] = pi[pos] / pi[pos].sum()
        
        return pi
    
    def _compute_probability_score(self, indices: np.ndarray, pi: torch.Tensor, optimizer) -> float:
        """计算序列概率分数"""
        log_prob = 0.0
        
        for pos, idx in enumerate(indices):
            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):
                choice = optimizer.codon_choices[pos][idx]
                orig_idx = choice.get('original_local_index', idx)
                
                if orig_idx < pi.shape[1]:
                    prob = pi[pos, orig_idx].item()
                    if prob > 0:
                        log_prob += np.log(prob)
        
        return log_prob
    
    def benchmark_speed(self, num_runs: int = 5) -> None:

        logger.info("=" * 80)

        logger.info("=" * 80)
        
        test_cases = [
            ('short', 30),
            ('medium', 100),
            ('long', 200)
        ]
        
        for case_name, seq_length in test_cases:

            logger.info("-" * 60)
            


            pi_test = self._create_test_distribution(seq_length)
            

            optimizer = SADOOptimizer(species='ecoli_bl21de3', device=self.device)
            
            times = []
            cais = []
            probs = []
            
            for run in range(num_runs):
                optimizer.reset()
                
                start_time = time.perf_counter()
                
                result, metadata = optimizer.optimize(
                    pi_accessibility=pi_test,
                    target_cai=0.8,
                    amino_acid_sequence=test_seq
                )
                
                end_time = time.perf_counter()

                
                times.append(execution_time)
                cais.append(metadata['final_cai'])
                

                if optimizer.last_indices is not None:
                    prob_score = self._compute_probability_score(
                        optimizer.last_indices, pi_test, optimizer
                    )
                    probs.append(prob_score)
                

                self.results.append(PerformanceResult(
                    optimizer_name='SADO',
                    sequence_name=case_name,
                    sequence_length=seq_length,
                    execution_time_ms=execution_time,
                    final_cai=metadata['final_cai'],
                    constraint_satisfied=metadata['constraint_satisfied'],
                    probability_score=prob_score if probs else 0.0,
                    unique_sequences=metadata['unique_sequences'],
                    collision_detected=metadata['collision_detected']
                ))
            

            avg_time = np.mean(times)
            std_time = np.std(times)
            avg_cai = np.mean(cais)
            avg_prob = np.mean(probs) if probs else 0.0
            




    
    def test_real_sequences(self) -> None:
        """使用真实E.coli序列进行测试"""
        logger.info("=" * 80)
        logger.info("真实序列性能测试")
        logger.info("=" * 80)
        

        test_sequences = {}
        sequences_by_length = sorted(
            [(name, seq, len(seq)) for name, seq in self.sequence_data.items()],
            key=lambda x: x[2]
        )
        

        length_ranges = [
            (20, 50, "short"),
            (50, 100, "medium"), 
            (100, 200, "long"),
            (200, 500, "very_long")
        ]
        
        for min_len, max_len, category in length_ranges:
            candidates = [(name, seq, length) for name, seq, length in sequences_by_length 
                         if min_len <= length <= max_len]
            if candidates:

                mid_idx = len(candidates) // 2
                name, seq, length = candidates[mid_idx]
                test_sequences[f"{category}_{name}"] = seq
        
        logger.info(f"选择了{len(test_sequences)}个代表序列进行测试")
        
        for seq_name, amino_seq in test_sequences.items():
            logger.info(f"\n测试序列: {seq_name} (长度: {len(amino_seq)})")
            logger.info("-" * 60)
            

            pi_test = self._create_test_distribution(len(amino_seq))
            

            optimizer = SADOOptimizer(species='ecoli_bl21de3', device=self.device)
            

            start_time = time.perf_counter()
            
            try:
                result, metadata = optimizer.optimize(
                    pi_accessibility=pi_test,
                    target_cai=0.8,
                    amino_acid_sequence=amino_seq
                )
                
                end_time = time.perf_counter()
                execution_time = (end_time - start_time) * 1000
                

                prob_score = 0.0
                if optimizer.last_indices is not None:
                    prob_score = self._compute_probability_score(
                        optimizer.last_indices, pi_test, optimizer
                    )
                

                self.results.append(PerformanceResult(
                    optimizer_name='SADO_real',
                    sequence_name=seq_name,
                    sequence_length=len(amino_seq),
                    execution_time_ms=execution_time,
                    final_cai=metadata['final_cai'],
                    constraint_satisfied=metadata['constraint_satisfied'],
                    probability_score=prob_score,
                    unique_sequences=metadata['unique_sequences'],
                    collision_detected=metadata['collision_detected']
                ))
                
                logger.info(f"执行时间: {execution_time:.2f} ms")
                logger.info(f"最终CAI: {metadata['final_cai']:.4f}")
                logger.info(f"约束满足: {'✅' if metadata['constraint_satisfied'] else '❌'}")
                logger.info(f"log(P): {prob_score:.2f}")
                logger.info(f"唯一序列数: {metadata['unique_sequences']}")
                
            except Exception as e:
                logger.error(f"序列{seq_name}优化失败: {e}")
    
    def test_scalability(self) -> None:

        logger.info("=" * 80)

        logger.info("=" * 80)
        

        lengths = [20, 50, 100, 200, 300, 500, 1000]
        
        for length in lengths:

            logger.info("-" * 40)
            

            test_seq = 'M' + 'A' * (length - 1)
            pi_test = self._create_test_distribution(length)
            

            optimizer = SADOOptimizer(species='ecoli_bl21de3', device=self.device)
            
            try:
                start_time = time.perf_counter()
                
                result, metadata = optimizer.optimize(
                    pi_accessibility=pi_test,
                    target_cai=0.8,
                    amino_acid_sequence=test_seq
                )
                
                end_time = time.perf_counter()
                execution_time = (end_time - start_time) * 1000
                

                time_per_aa = execution_time / length
                


                logger.info(f"CAI: {metadata['final_cai']:.4f}")

                

                self.results.append(PerformanceResult(
                    optimizer_name='SADO_scale',
                    sequence_name=f'length_{length}',
                    sequence_length=length,
                    execution_time_ms=execution_time,
                    final_cai=metadata['final_cai'],
                    constraint_satisfied=metadata['constraint_satisfied'],

                    unique_sequences=metadata['unique_sequences'],
                    collision_detected=metadata['collision_detected']
                ))
                
            except Exception as e:

    
    def test_stability(self, num_runs: int = 10) -> None:
        """稳定性测试"""
        logger.info("=" * 80)
        logger.info("稳定性测试")
        logger.info("=" * 80)
        
        test_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
        pi_test = self._create_test_distribution(len(test_seq))
        
        results = []
        times = []
        
        for run in range(num_runs):
            optimizer = SADOOptimizer(species='ecoli_bl21de3', device=self.device)
            
            start_time = time.perf_counter()
            result, metadata = optimizer.optimize(
                pi_accessibility=pi_test,
                target_cai=0.8,
                amino_acid_sequence=test_seq
            )
            end_time = time.perf_counter()
            
            execution_time = (end_time - start_time) * 1000
            times.append(execution_time)
            results.append(metadata['final_cai'])
        

        avg_time = np.mean(times)
        std_time = np.std(times)
        avg_cai = np.mean(results)
        std_cai = np.std(results)
        min_cai = np.min(results)
        max_cai = np.max(results)
        
        logger.info(f"运行次数: {num_runs}")
        logger.info(f"平均执行时间: {avg_time:.2f} ± {std_time:.2f} ms")
        logger.info(f"时间变异系数: {std_time/avg_time*100:.1f}%")
        logger.info(f"平均CAI: {avg_cai:.4f} ± {std_cai:.4f}")
        logger.info(f"CAI范围: [{min_cai:.4f}, {max_cai:.4f}]")
        logger.info(f"CAI变异系数: {std_cai/avg_cai*100:.1f}%")
        logger.info(f"约束满足率: {sum(c >= 0.8 for c in results)}/{num_runs}")
    
    def generate_report(self) -> None:

        logger.info("=" * 80)

        logger.info("=" * 80)
        
        if not self.results:

            return
        

        grouped = defaultdict(list)
        for result in self.results:
            grouped[result.optimizer_name].append(result)
        

        for optimizer_name, results in grouped.items():

            logger.info("-" * 50)
            
            times = [r.execution_time_ms for r in results]
            cais = [r.final_cai for r in results]
            lengths = [r.sequence_length for r in results]
            satisfied = [r.constraint_satisfied for r in results]
            





            
            if lengths and times:

                time_per_aa = [t/l for t, l in zip(times, lengths)]

        

        self._save_csv_report()
        


    
    def _save_csv_report(self) -> None:
        """保存CSV报告"""
        if not self.results:
            return
            

        data = []
        for result in self.results:
            data.append({
                'optimizer': result.optimizer_name,
                'sequence': result.sequence_name,
                'length': result.sequence_length,
                'time_ms': result.execution_time_ms,
                'cai': result.final_cai,
                'constraint_satisfied': result.constraint_satisfied,
                'log_probability': result.probability_score,
                'unique_sequences': result.unique_sequences,
                'collision_detected': result.collision_detected,
                'time_per_aa_ms': result.execution_time_ms / result.sequence_length
            })
        
        df = pd.DataFrame(data)
        

        output_file = Path(__file__).parent / "sado_performance_results.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"结果保存到: {output_file}")
    
    def run_comprehensive_test(self) -> None:




        
        try:

            self.benchmark_speed(num_runs=3)
            

            self.test_real_sequences()
            

            self.test_scalability()
            

            self.test_stability(num_runs=5)
            

            self.generate_report()
            
        except Exception as e:

            raise


def main():
    """主函数"""
    tester = SADOPerformanceTester()
    tester.run_comprehensive_test()


if __name__ == "__main__":
    main()
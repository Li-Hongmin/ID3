#!/usr/bin/env python3
"""




"""

import torch
import numpy as np
from pathlib import Path
import json
from typing import Dict, List, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


torch.manual_seed(42)
np.random.seed(42)


import sys
sys.path.append(str(Path(__file__).parent.parent.parent))

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.experiments.utils.data_loader import ProteinDataLoader
from id3.utils.constants import amino_acids_to_codons


class CAIReachabilityAnalyzer:

    
    def __init__(self, species: str = 'ecoli_bl21de3'):
        """
        初始化分析器
        
        Args:
            species: 物种（用于CAI权重）
        """
        self.species = species
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.cai_operator = CAIEnhancementOperator(species=species, device=self.device)
        self.data_loader = ProteinDataLoader()
        

        self.wi_table = self.cai_operator.wi_table
        

        self.results = {}
        
    def get_test_proteins(self) -> List[str]:
        """获取测试蛋白质列表"""

        proteins = [
            "P0CG48",  # Ubiquitin
            "P04637",  # p53 tumor suppressor
            "P0DTC2",  # SARS-CoV-2 spike protein
            "P01308",  # Insulin
            "P00004",  # Cytochrome c oxidase
            "P42212",  # Green fluorescent protein (GFP)
            "P0DTC9",  # SARS-CoV-2 nucleocapsid
            "O15263",
            "P31417",  # CRHBP
            "P01825",
            "Q9BWJ5",
            "P63165"
        ]
        

        available_proteins = []
        for protein_id in proteins:
            try:
                seq = self.data_loader.load_protein_sequence(protein_id)
                if seq:
                    available_proteins.append(protein_id)
            except:
                print(f"警告: 无法加载蛋白质 {protein_id}")
                
        return available_proteins
    
    def compute_theoretical_max_cai(self, amino_acid_sequence: str) -> Tuple[float, Dict]:
        """

        
        Args:

            
        Returns:
            (max_cai, details)
        """
        log_weights = []
        codon_choices = []
        aa_weights = {}
        
        for aa in amino_acid_sequence:
            if aa not in amino_acids_to_codons:
                continue
                
            codons = amino_acids_to_codons[aa]
            weights = []
            

            for codon in codons:

                rna_codon = codon
                if rna_codon in self.wi_table:
                    weights.append(self.wi_table[rna_codon])
                else:
                    weights.append(0.1)
            

            max_weight = max(weights) if weights else 0.1
            max_idx = weights.index(max_weight) if weights else 0
            
            log_weights.append(np.log(max(max_weight, 1e-10)))
            codon_choices.append(codons[max_idx])
            

            if aa not in aa_weights:
                aa_weights[aa] = {
                    'codons': codons,
                    'weights': weights,
                    'max_weight': max_weight,
                    'best_codon': codons[max_idx]
                }
        

        if log_weights:
            mean_log_weight = np.mean(log_weights)
            max_cai = np.exp(mean_log_weight)
        else:
            max_cai = 0.0
        
        return max_cai, {
            'codon_choices': codon_choices,
            'aa_weights': aa_weights,
            'num_positions': len(amino_acid_sequence)
        }
    
    def analyze_amino_acid_composition(self, amino_acid_sequence: str) -> Dict:
        """

        
        Args:

            
        Returns:

        """
        aa_counts = Counter(amino_acid_sequence)
        total_aas = len(amino_acid_sequence)
        

        aa_analysis = {}
        problem_aas = []
        weight_threshold = 0.5
        
        for aa, count in aa_counts.items():
            if aa not in amino_acids_to_codons:
                continue
                
            codons = amino_acids_to_codons[aa]
            weights = []
            
            for codon in codons:
                rna_codon = codon
                if rna_codon in self.wi_table:
                    weights.append(self.wi_table[rna_codon])
                else:
                    weights.append(0.1)
            
            max_weight = max(weights) if weights else 0.1
            avg_weight = np.mean(weights) if weights else 0.1
            
            aa_analysis[aa] = {
                'count': count,
                'frequency': count / total_aas,
                'num_codons': len(codons),
                'max_weight': max_weight,
                'avg_weight': avg_weight,
                'weights': weights,
                'is_problem': max_weight < weight_threshold
            }
            
            if max_weight < weight_threshold:
                problem_aas.append(aa)
        

        problem_aa_count = sum(aa_analysis[aa]['count'] for aa in problem_aas if aa in aa_analysis)
        problem_aa_ratio = problem_aa_count / total_aas if total_aas > 0 else 0
        
        return {
            'aa_counts': dict(aa_counts),
            'aa_analysis': aa_analysis,
            'problem_aas': problem_aas,
            'problem_aa_ratio': problem_aa_ratio,
            'total_positions': total_aas
        }
    
    def test_actual_reachability(self, amino_acid_sequence: str, 
                                target_cai: float = 0.8,
                                num_iterations: int = 100) -> Dict:
        """

        
        Args:



            
        Returns:

        """

        self.cai_operator._load_or_compute_amino_acid_cache(amino_acid_sequence)
        

        max_achievable_cai = self.cai_operator.max_achievable_cai
        

        num_positions = len(amino_acid_sequence)
        max_codons = 6
        

        valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=self.device)
        codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=self.device)
        
        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                for i, codon in enumerate(codons):
                    if i < max_codons:
                        valid_codon_mask[pos, i] = True
                        codon_idx = self.cai_operator._codon_to_standard_index(codon)
                        codon_indices[pos, i] = codon_idx
        

        initial_dist = torch.zeros(num_positions, max_codons, device=self.device)
        for pos in range(num_positions):
            valid_count = valid_codon_mask[pos].sum()
            if valid_count > 0:
                initial_dist[pos][valid_codon_mask[pos]] = 1.0 / valid_count
        

        try:
            enhanced_dist, metadata = self.cai_operator.apply_cai_enhancement(
                initial_dist, amino_acid_sequence, valid_codon_mask, 
                codon_indices, target_cai
            )
            
            actual_cai = metadata['final_cai']
            success = metadata['constraint_satisfied']
            gamma = metadata['optimal_gamma']
            
        except Exception as e:
            print(f"优化失败: {e}")
            actual_cai = 0.0
            success = False
            gamma = 0.0
        
        return {
            'max_achievable_cai': max_achievable_cai,
            'target_cai': target_cai,
            'actual_cai': actual_cai,
            'success': success,
            'gamma': gamma,
            'gap': target_cai - actual_cai if actual_cai > 0 else target_cai
        }
    
    def analyze_species_specificity(self) -> Dict:
        """

        
        Returns:

        """

        species_analysis = {
            'species': self.species,
            'codon_weights': {},
            'aa_max_weights': {},
            'low_weight_codons': [],
            'high_weight_codons': []
        }
        

        for aa, codons in amino_acids_to_codons.items():
            weights = []
            codon_weight_map = {}
            
            for codon in codons:
                rna_codon = codon
                if rna_codon in self.wi_table:
                    weight = self.wi_table[rna_codon]
                else:
                    weight = 0.1
                    
                weights.append(weight)
                codon_weight_map[codon] = weight
                

                if weight >= 0.8:
                    species_analysis['high_weight_codons'].append((codon, aa, weight))
                elif weight <= 0.2:
                    species_analysis['low_weight_codons'].append((codon, aa, weight))
            
            max_weight = max(weights) if weights else 0.1
            species_analysis['aa_max_weights'][aa] = max_weight
            species_analysis['codon_weights'][aa] = codon_weight_map
        

        all_weights = [w for aa_weights in species_analysis['codon_weights'].values() 
                      for w in aa_weights.values()]
        species_analysis['weight_stats'] = {
            'mean': np.mean(all_weights),
            'std': np.std(all_weights),
            'min': min(all_weights),
            'max': max(all_weights),
            'median': np.median(all_weights)
        }
        
        return species_analysis
    
    def run_comprehensive_analysis(self, save_results: bool = True) -> Dict:
        """

        
        Args:

            
        Returns:

        """
        print("=" * 80)
        print("CAI目标可达性综合分析")
        print("=" * 80)
        

        proteins = self.get_test_proteins()
        print(f"\n找到 {len(proteins)} 个测试蛋白质")
        

        print("\n1. 分析物种特异性...")
        species_analysis = self.analyze_species_specificity()
        print(f"   - 物种: {species_analysis['species']}")
        print(f"   - 权重统计: mean={species_analysis['weight_stats']['mean']:.3f}, "
              f"std={species_analysis['weight_stats']['std']:.3f}")
        print(f"   - 高权重密码子(≥0.8): {len(species_analysis['high_weight_codons'])}个")
        print(f"   - 低权重密码子(≤0.2): {len(species_analysis['low_weight_codons'])}个")
        

        protein_results = {}
        
        for protein_id in proteins:
            print(f"\n2. 分析蛋白质 {protein_id}...")
            
            try:

                sequence = self.data_loader.load_protein_sequence(protein_id)
                print(f"   序列长度: {len(sequence)}")
                

                max_cai, max_cai_details = self.compute_theoretical_max_cai(sequence)
                print(f"   理论最大CAI: {max_cai:.4f}")
                

                aa_analysis = self.analyze_amino_acid_composition(sequence)
                print(f"   问题氨基酸比例: {aa_analysis['problem_aa_ratio']:.2%}")
                if aa_analysis['problem_aas']:
                    print(f"   问题氨基酸: {', '.join(aa_analysis['problem_aas'])}")
                

                reachability = self.test_actual_reachability(sequence, target_cai=0.8)
                print(f"   实际可达CAI: {reachability['actual_cai']:.4f}")
                print(f"   目标达成: {'✓' if reachability['success'] else '✗'}")
                

                protein_results[protein_id] = {
                    'sequence_length': len(sequence),
                    'theoretical_max_cai': max_cai,
                    'actual_max_cai': reachability['max_achievable_cai'],
                    'target_cai': 0.8,
                    'achieved_cai': reachability['actual_cai'],
                    'success': reachability['success'],
                    'gap': reachability['gap'],
                    'problem_aa_ratio': aa_analysis['problem_aa_ratio'],
                    'problem_aas': aa_analysis['problem_aas'],
                    'aa_analysis': aa_analysis['aa_analysis']
                }
                
            except Exception as e:
                print(f"   错误: {e}")
                protein_results[protein_id] = {'error': str(e)}
        

        print("\n" + "=" * 80)
        print("汇总分析")
        print("=" * 80)
        

        successful = sum(1 for r in protein_results.values() 
                        if 'success' in r and r['success'])
        total = sum(1 for r in protein_results.values() if 'success' in r)
        success_rate = successful / total if total > 0 else 0
        
        print(f"\n目标达成率: {success_rate:.1%} ({successful}/{total})")
        

        print("\n无法达到0.8目标的蛋白质分析:")
        for protein_id, result in protein_results.items():
            if 'success' in result and not result['success']:
                print(f"\n{protein_id}:")
                print(f"  - 理论最大: {result['theoretical_max_cai']:.4f}")
                print(f"  - 实际最大: {result['actual_max_cai']:.4f}")
                print(f"  - 差距: {result['gap']:.4f}")
                print(f"  - 问题氨基酸比例: {result['problem_aa_ratio']:.2%}")
        

        print("\n" + "=" * 80)
        print("建议")
        print("=" * 80)
        

        max_cais = [r['theoretical_max_cai'] for r in protein_results.values() 
                   if 'theoretical_max_cai' in r]
        if max_cais:
            percentiles = np.percentile(max_cais, [25, 50, 75, 90])
            print(f"\n理论最大CAI分布:")
            print(f"  - 25%分位数: {percentiles[0]:.4f}")
            print(f"  - 中位数: {percentiles[1]:.4f}")
            print(f"  - 75%分位数: {percentiles[2]:.4f}")
            print(f"  - 90%分位数: {percentiles[3]:.4f}")
            
            print(f"\n建议:")
            print(f"  1. 通用目标CAI = {percentiles[1]:.3f} (中位数)")
            print(f"  2. 宽松目标CAI = {percentiles[0]:.3f} (25%分位数)")
            print(f"  3. 严格目标CAI = {percentiles[2]:.3f} (75%分位数)")
            print(f"  4. 使用自适应目标: min(0.8, 0.95 * 理论最大CAI)")
        

        if save_results:
            results_dir = Path("id3/experiments/analysis")
            results_dir.mkdir(parents=True, exist_ok=True)
            

            results_file = results_dir / "cai_reachability_analysis.json"
            with open(results_file, 'w') as f:
                json.dump({
                    'species_analysis': species_analysis,
                    'protein_results': protein_results,
                    'summary': {
                        'success_rate': success_rate,
                        'num_proteins': len(proteins),
                        'percentiles': percentiles.tolist() if max_cais else []
                    }
                }, f, indent=2, default=str)
            print(f"\n结果已保存到: {results_file}")
            

            self.visualize_results(protein_results, species_analysis, results_dir)
        
        return {
            'species_analysis': species_analysis,
            'protein_results': protein_results,
            'summary': {
                'success_rate': success_rate,
                'num_proteins': len(proteins)
            }
        }
    
    def visualize_results(self, protein_results: Dict, species_analysis: Dict, 
                         save_dir: Path):
        """

        
        Args:



        """

        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_palette("husl")
        

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        

        ax = axes[0, 0]
        theoretical = []
        actual = []
        labels = []
        
        for protein_id, result in protein_results.items():
            if 'theoretical_max_cai' in result:
                theoretical.append(result['theoretical_max_cai'])
                actual.append(result['actual_max_cai'])
                labels.append(protein_id)
        
        ax.scatter(theoretical, actual, s=100, alpha=0.6)
        ax.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='y=x')
        ax.axhline(y=0.8, color='g', linestyle='--', alpha=0.5, label='Target=0.8')
        
        for i, label in enumerate(labels):
            ax.annotate(label, (theoretical[i], actual[i]), 
                       fontsize=8, alpha=0.7)
        
        ax.set_xlabel('理论最大CAI')
        ax.set_ylabel('实际最大CAI')
        ax.set_title('理论 vs 实际最大CAI')
        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[0, 1]
        problem_ratios = []
        gaps = []
        
        for protein_id, result in protein_results.items():
            if 'problem_aa_ratio' in result:
                problem_ratios.append(result['problem_aa_ratio'])
                gaps.append(result['gap'])
        
        ax.scatter(problem_ratios, gaps, s=100, alpha=0.6)
        ax.set_xlabel('问题氨基酸比例')
        ax.set_ylabel('CAI差距 (目标 - 实际)')
        ax.set_title('问题氨基酸对CAI的影响')
        ax.grid(True, alpha=0.3)
        

        ax = axes[1, 0]
        aa_weights = list(species_analysis['aa_max_weights'].values())
        ax.hist(aa_weights, bins=20, edgecolor='black', alpha=0.7)
        ax.axvline(x=0.5, color='r', linestyle='--', alpha=0.5, 
                  label='问题阈值=0.5')
        ax.set_xlabel('氨基酸最大权重')
        ax.set_ylabel('数量')
        ax.set_title('氨基酸最大权重分布')
        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[1, 1]
        proteins = []
        max_cais = []
        colors = []
        
        for protein_id, result in protein_results.items():
            if 'theoretical_max_cai' in result:
                proteins.append(protein_id)
                max_cais.append(result['theoretical_max_cai'])
                colors.append('green' if result['success'] else 'red')
        
        bars = ax.bar(range(len(proteins)), max_cais, color=colors, alpha=0.6)
        ax.axhline(y=0.8, color='b', linestyle='--', alpha=0.5, 
                  label='目标=0.8')
        ax.set_xticks(range(len(proteins)))
        ax.set_xticklabels(proteins, rotation=45, ha='right')
        ax.set_ylabel('理论最大CAI')
        ax.set_title('蛋白质CAI目标达成情况')
        ax.legend()
        ax.grid(True, alpha=0.3)
        

        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='green', alpha=0.6, label='可达目标'),
            Patch(facecolor='red', alpha=0.6, label='无法达标')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        

        plot_file = save_dir / "cai_reachability_visualization.png"
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"可视化已保存到: {plot_file}")
        plt.close()


def main():

    analyzer = CAIReachabilityAnalyzer(species='ecoli_bl21de3')
    results = analyzer.run_comprehensive_analysis(save_results=True)
    
    print("\n" + "=" * 80)

    print("=" * 80)
    
    return results


if __name__ == "__main__":
    main()
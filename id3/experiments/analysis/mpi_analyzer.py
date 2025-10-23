#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import os
import sys


sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from id3.experiments.analysis.fast_analyzer import FastExperimentAnalyzer

class MPIExperimentAnalyzer:

    
    def __init__(self):
        self.fast_analyzer = FastExperimentAnalyzer()
    
    def is_mpi_results_dir(self, dir_path: Path) -> bool:
        """检查是否是MPI结果目录"""
        if not dir_path.exists():
            return False
        

        seed_dirs = list(dir_path.glob("seed*"))
        if not seed_dirs:
            return False
        

        worker_files = list(dir_path.glob("*/worker_*_results.json"))
        return len(worker_files) > 0
    
    def collect_worker_results(self, mpi_dir: Path) -> List[Dict[str, Any]]:

        worker_results = []
        

        worker_files = list(mpi_dir.glob("*/worker_*_results.json"))
        

        
        for worker_file in worker_files:
            try:
                with open(worker_file, 'r') as f:
                    worker_data = json.load(f)
                

                worker_config = worker_data.get('worker_config', {})
                experiments = worker_data.get('experiments', [])
                summary = worker_data.get('summary', {})
                

                for exp in experiments:
                    exp_result = {
                        'worker_id': worker_config.get('worker_id'),
                        'seed': worker_config.get('seed'),
                        'seed_idx': worker_config.get('seed_idx'),
                        'variant': worker_config.get('variant'),
                        'constraint': worker_config.get('constraint'),
                        'protein': exp.get('protein'),
                        'status': exp.get('status'),
                        'execution_time': exp.get('execution_time'),
                        'output_dir': exp.get('output_dir'),
                        'worker_summary': summary
                    }
                    worker_results.append(exp_result)
                    
            except Exception as e:

        
        return worker_results
    
    def extract_experiment_data(self, worker_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """从每个实验目录提取详细数据"""
        experiment_data = []
        
        print(f"📈 提取 {len(worker_results)} 个实验的详细数据...")
        
        for i, worker_result in enumerate(worker_results):
            if i % 100 == 0:
                print(f"  进度: {i}/{len(worker_results)}")
            
            if worker_result['status'] != 'success':
                continue
            
            output_dir = Path(worker_result['output_dir'])
            if not output_dir.exists():
                continue
            

            json_files = list(output_dir.glob("*.json"))
            exp_files = [f for f in json_files 
                        if not f.name.startswith(('config', 'progress', 'summary'))]
            
            if not exp_files:
                continue
            

            exp_file = exp_files[0]
            exp_data = self.fast_analyzer.extract_key_fields(exp_file)
            
            if exp_data:

                exp_data.update({
                    'worker_id': worker_result['worker_id'],
                    'seed': worker_result['seed'],
                    'seed_idx': worker_result['seed_idx'],
                    'variant': worker_result['variant'],
                    'constraint_type': worker_result['constraint'],
                    'mpi_execution_time': worker_result['execution_time']
                })
                experiment_data.append(exp_data)
        
        print(f"✅ 成功提取 {len(experiment_data)} 个实验的数据")
        return experiment_data
    
    def analyze_mpi_experiments(self, mpi_dir: Path) -> Dict[str, Any]:


        

        worker_results = self.collect_worker_results(mpi_dir)
        if not worker_results:

        

        experiment_data = self.extract_experiment_data(worker_results)
        if not experiment_data:

        

        df = pd.DataFrame(experiment_data)
        

        analysis = {
            'basic_stats': {
                'total_workers': len(set(wr['worker_id'] for wr in worker_results)),
                'total_experiments': len(df),
                'successful_experiments': len(df[df.get('status') == 'success']) if 'status' in df else len(df),
                'unique_proteins': df['protein_name'].nunique() if 'protein_name' in df else 0,
                'unique_constraints': df['constraint_type'].nunique() if 'constraint_type' in df else 0,
                'unique_variants': df['variant'].nunique() if 'variant' in df else 0,
                'unique_seeds': df['seed'].nunique() if 'seed' in df else 0,
                'seeds_range': (df['seed'].min(), df['seed'].max()) if 'seed' in df else None
            }
        }
        

        if 'best_accessibility' in df:
            accessibility_stats = {
                'mean': df['best_accessibility'].mean(),
                'median': df['best_accessibility'].median(),
                'std': df['best_accessibility'].std(),
                'min': df['best_accessibility'].min(),
                'max': df['best_accessibility'].max(),
                'excellent_rate': (df['best_accessibility'] < 1.0).mean() * 100,
                'good_rate': (df['best_accessibility'] < 1.5).mean() * 100
            }
            analysis['accessibility_performance'] = accessibility_stats
        

        if 'constraint_type' in df and 'best_accessibility' in df:
            constraint_stats = {}
            for constraint in df['constraint_type'].unique():
                constraint_df = df[df['constraint_type'] == constraint]
                constraint_stats[constraint] = {
                    'count': len(constraint_df),
                    'mean': constraint_df['best_accessibility'].mean(),
                    'median': constraint_df['best_accessibility'].median(),
                    'min': constraint_df['best_accessibility'].min(),
                    'max': constraint_df['best_accessibility'].max(),
                    'std': constraint_df['best_accessibility'].std()
                }
            analysis['by_constraint'] = constraint_stats
        

        if 'variant' in df and 'best_accessibility' in df:
            variant_stats = {}
            for variant in df['variant'].unique():
                variant_df = df[df['variant'] == variant]
                variant_stats[variant] = {
                    'count': len(variant_df),
                    'mean': variant_df['best_accessibility'].mean(),
                    'median': variant_df['best_accessibility'].median(),
                    'min': variant_df['best_accessibility'].min()
                }
            analysis['by_variant'] = variant_stats
        

        if 'seed' in df and 'best_accessibility' in df:
            seed_stats = {}
            for seed in sorted(df['seed'].unique()):
                seed_df = df[df['seed'] == seed]
                seed_stats[f'seed{seed}'] = {
                    'count': len(seed_df),
                    'mean': seed_df['best_accessibility'].mean(),
                    'median': seed_df['best_accessibility'].median(),
                    'min': seed_df['best_accessibility'].min()
                }
            analysis['by_seed'] = seed_stats
        

        if 'amino_acids_match' in df:
            aa_match_rate = df['amino_acids_match'].mean() * 100
            analysis['amino_acid_match_rate'] = aa_match_rate
        
        if 'amino_acids_correct' in df:
            aa_correct_rate = df['amino_acids_correct'].mean()
            analysis['amino_acid_correct_rate'] = aa_correct_rate
        

        if 'final_ecai' in df and df['final_ecai'].notna().any():
            ecai_data = df['final_ecai'].dropna()
            cai_stats = {
                'mean': ecai_data.mean(),
                'median': ecai_data.median(),
                'std': ecai_data.std(),
                'min': ecai_data.min(),
                'max': ecai_data.max()
            }
            if 'cai_target_achieved' in df:
                cai_stats['target_achieved_rate'] = df['cai_target_achieved'].mean() * 100
            analysis['cai_performance'] = cai_stats
        

        if 'mpi_execution_time' in df:
            exec_time_stats = {
                'mean_seconds': df['mpi_execution_time'].mean(),
                'median_seconds': df['mpi_execution_time'].median(),
                'total_hours': df['mpi_execution_time'].sum() / 3600,
                'min_seconds': df['mpi_execution_time'].min(),
                'max_seconds': df['mpi_execution_time'].max()
            }
            analysis['execution_time'] = exec_time_stats
        

        analysis['dataframe'] = df
        
        return analysis
    
    def print_analysis_report(self, analysis: Dict[str, Any]) -> None:
        """打印分析报告"""
        if 'error' in analysis:
            print(f"❌ {analysis['error']}")
            return
        
        basic = analysis['basic_stats']
        
        print("\n" + "="*80)
        print("🎯 MPI实验分析报告")
        print("="*80)
        
        print(f"\n📊 基础统计:")
        print(f"  • 总Workers: {basic['total_workers']}")
        print(f"  • 总实验数: {basic['total_experiments']}")
        print(f"  • 成功实验: {basic['successful_experiments']}")
        print(f"  • 蛋白质数: {basic['unique_proteins']}")
        print(f"  • 约束类型: {basic['unique_constraints']}")
        print(f"  • 变体类型: {basic['unique_variants']}")
        print(f"  • Seeds数量: {basic['unique_seeds']}")
        if basic['seeds_range']:
            print(f"  • Seeds范围: {basic['seeds_range'][0]}-{basic['seeds_range'][1]}")
        

        if 'accessibility_performance' in analysis:
            acc = analysis['accessibility_performance']
            print(f"\n🎯 Accessibility性能:")
            print(f"  • 平均值: {acc['mean']:.4f} kcal/mol")
            print(f"  • 中位数: {acc['median']:.4f} kcal/mol")
            print(f"  • 最小值: {acc['min']:.4f} kcal/mol")
            print(f"  • 最大值: {acc['max']:.4f} kcal/mol")
            print(f"  • 标准差: {acc['std']:.4f}")
            print(f"  • 优秀率(<1.0): {acc['excellent_rate']:.1f}%")
            print(f"  • 良好率(<1.5): {acc['good_rate']:.1f}%")
        

        if 'by_constraint' in analysis:
            print(f"\n📊 按约束类型统计:")
            for constraint, stats in analysis['by_constraint'].items():
                print(f"  {constraint.upper()}:")
                print(f"    • 实验数: {stats['count']}")
                print(f"    • 平均: {stats['mean']:.4f} kcal/mol")
                print(f"    • 最佳: {stats['min']:.4f} kcal/mol")
        

        if 'by_variant' in analysis:
            print(f"\n🧬 按变体统计:")
            for variant, stats in analysis['by_variant'].items():
                print(f"  变体{variant}: 平均 {stats['mean']:.4f}, 最佳 {stats['min']:.4f}")
        

        if 'by_seed' in analysis:
            print(f"\n🌱 按Seed统计:")
            for seed, stats in analysis['by_seed'].items():
                print(f"  {seed}: 平均 {stats['mean']:.4f}, 最佳 {stats['min']:.4f} ({stats['count']}实验)")
        

        if 'amino_acid_match_rate' in analysis:
            print(f"\n🧬 氨基酸匹配:")
            print(f"  • 完全匹配率: {analysis['amino_acid_match_rate']:.1f}%")
        
        if 'amino_acid_correct_rate' in analysis:
            print(f"  • 平均正确率: {analysis['amino_acid_correct_rate']:.1f}%")
        

        if 'cai_performance' in analysis:
            cai = analysis['cai_performance']
            print(f"\n🎨 CAI优化结果:")
            print(f"  • 平均CAI: {cai['mean']:.4f}")
            print(f"  • 中位数CAI: {cai['median']:.4f}")
            print(f"  • 最小CAI: {cai['min']:.4f}")
            print(f"  • 最大CAI: {cai['max']:.4f}")
            if 'target_achieved_rate' in cai:
                print(f"  • 达到目标(0.8): {cai['target_achieved_rate']:.1f}%")
        

        if 'execution_time' in analysis:
            exec_time = analysis['execution_time']
            print(f"\n⏱️ 执行时间统计:")
            print(f"  • 平均时间: {exec_time['mean_seconds']:.1f}秒")
            print(f"  • 中位数时间: {exec_time['median_seconds']:.1f}秒")
            print(f"  • 总计时间: {exec_time['total_hours']:.1f}小时")
            print(f"  • 最快: {exec_time['min_seconds']:.1f}秒")
            print(f"  • 最慢: {exec_time['max_seconds']:.1f}秒")

    def compare_mpi_experiments(self, dir1: Path, dir2: Path) -> Dict[str, Any]:


        print("="*80)
        
        analysis1 = self.analyze_mpi_experiments(dir1)
        analysis2 = self.analyze_mpi_experiments(dir2)
        
        if 'error' in analysis1 or 'error' in analysis2:

        

        comparison = {
            'dir1_name': dir1.name,
            'dir2_name': dir2.name,
            'dir1_stats': analysis1.get('accessibility_performance'),
            'dir2_stats': analysis2.get('accessibility_performance')
        }
        
        if comparison['dir1_stats'] and comparison['dir2_stats']:
            stats1 = comparison['dir1_stats']
            stats2 = comparison['dir2_stats']
            
            comparison['differences'] = {
                'mean_diff': stats1['mean'] - stats2['mean'],
                'median_diff': stats1['median'] - stats2['median'],
                'min_diff': stats1['min'] - stats2['min'],
                'excellent_rate_diff': stats1['excellent_rate'] - stats2['excellent_rate']
            }
        
        return comparison
    
    def print_comparison_report(self, comparison: Dict[str, Any]) -> None:
        """打印对比报告"""
        if 'error' in comparison:
            print(f"❌ {comparison['error']}")
            return
        
        print(f"\n📊 性能对比:")
        print(f"  目录1 ({comparison['dir1_name']}):")
        if comparison['dir1_stats']:
            stats1 = comparison['dir1_stats']
            print(f"    • 平均: {stats1['mean']:.4f} kcal/mol")
            print(f"    • 中位数: {stats1['median']:.4f} kcal/mol")
            print(f"    • 最优: {stats1['min']:.4f} kcal/mol")
            print(f"    • 优秀率: {stats1['excellent_rate']:.1f}%")
        
        print(f"  目录2 ({comparison['dir2_name']}):")
        if comparison['dir2_stats']:
            stats2 = comparison['dir2_stats']
            print(f"    • 平均: {stats2['mean']:.4f} kcal/mol")
            print(f"    • 中位数: {stats2['median']:.4f} kcal/mol")
            print(f"    • 最优: {stats2['min']:.4f} kcal/mol")
            print(f"    • 优秀率: {stats2['excellent_rate']:.1f}%")
        
        if 'differences' in comparison:
            diff = comparison['differences']
            print(f"\n  差异 (目录1 - 目录2):")
            print(f"    • 平均差异: {diff['mean_diff']:+.4f} kcal/mol")
            print(f"    • 中位数差异: {diff['median_diff']:+.4f} kcal/mol")
            print(f"    • 最优差异: {diff['min_diff']:+.4f} kcal/mol")
            print(f"    • 优秀率差异: {diff['excellent_rate_diff']:+.1f}%")

def main():

    analyzer = MPIExperimentAnalyzer()
    

    mpi_dir = Path("external_storage/results/mpi_flexible")
    if mpi_dir.exists():
        analysis = analyzer.analyze_mpi_experiments(mpi_dir)
        analyzer.print_analysis_report(analysis)
    else:


if __name__ == "__main__":
    main()
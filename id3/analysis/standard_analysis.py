#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
from pathlib import Path
import statistics
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

class StandardAnalyzer:

    
    def __init__(self):
        """初始化分析器"""
        self.results = {}
        sns.set_style("whitegrid")
    
    def analyze_experiment_batch(self, dir_path: Path) -> Dict:


        print("=" * 80)
        

        basic_stats = self.get_basic_statistics(dir_path)
        

        constraint_stats = self.analyze_by_constraint(dir_path)
        

        constraint_validation = self.validate_constraints(dir_path)
        

        anomalies = self.detect_anomalies(constraint_stats)
        

        relative_performance = self.analyze_relative_performance(constraint_stats)
        
        return {
            'directory': dir_path.name,
            'basic_stats': basic_stats,
            'constraint_stats': constraint_stats,
            'constraint_validation': constraint_validation,
            'anomalies': anomalies,
            'relative_performance': relative_performance
        }
    
    def get_basic_statistics(self, dir_path: Path) -> Dict:
        """获取基础统计信息"""
        all_values = []
        file_count = 0
        
        for json_file in dir_path.glob("*.json"):
            if json_file.name.startswith(('config', 'progress', 'summary')):
                continue
            
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    if 'best_accessibility' in data:
                        all_values.append(data['best_accessibility'])
                        file_count += 1
            except:
                pass
        
        if all_values:
            return {
                'total_experiments': file_count,
                'mean': statistics.mean(all_values),
                'median': statistics.median(all_values),
                'std': statistics.stdev(all_values) if len(all_values) > 1 else 0,
                'min': min(all_values),
                'max': max(all_values),
                'excellent_rate': sum(1 for v in all_values if v < 1.0) / len(all_values) * 100,
                'good_rate': sum(1 for v in all_values if v < 1.5) / len(all_values) * 100
            }
        return {}
    
    def analyze_by_constraint(self, dir_path: Path) -> Dict:

        constraint_data = {'lagrangian': [], 'ams': [], 'cpc': []}
        
        for json_file in dir_path.glob("*.json"):
            if json_file.name.startswith(('config', 'progress', 'summary')):
                continue
            
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    constraint = data.get('constraint_type')
                    if constraint in constraint_data:
                        constraint_data[constraint].append({
                            'accessibility': data.get('best_accessibility', 0),
                            'aa_correct': data.get('amino_acids_correct', 0),
                            'aa_match': data.get('amino_acids_match', False),
                            'protein': data.get('protein_name', 'unknown'),
                            'variant': data.get('variant', 'unknown')
                        })
            except:
                pass
        
        results = {}
        for constraint, values in constraint_data.items():
            if values:
                accs = [v['accessibility'] for v in values]
                aa_corrects = [v['aa_correct'] for v in values]
                
                results[constraint] = {
                    'count': len(values),
                    'mean': statistics.mean(accs),
                    'median': statistics.median(accs),
                    'std': statistics.stdev(accs) if len(accs) > 1 else 0,
                    'min': min(accs),
                    'max': max(accs),
                    'aa_correct_avg': statistics.mean(aa_corrects),
                    'aa_match_rate': sum(1 for v in values if v['aa_match']) / len(values) * 100,
                    'excellent_rate': sum(1 for a in accs if a < 1.0) / len(accs) * 100,

                }
        
        return results
    
    def validate_constraints(self, dir_path: Path, sample_size: int = 10) -> Dict:
        """验证约束满足情况"""
        import random
        from id3.analysis.verify_constraint_satisfaction import check_experiment_file
        
        json_files = list(dir_path.glob("*.json"))
        experiment_files = [f for f in json_files 
                           if not f.name.startswith(('config', 'progress', 'summary'))]
        

        if len(experiment_files) > sample_size:
            sampled_files = random.sample(experiment_files, sample_size)
        else:
            sampled_files = experiment_files
        
        validation_results = {
            'lagrangian': {'satisfied': 0, 'total': 0},
            'ams': {'satisfied': 0, 'total': 0},
            'cpc': {'satisfied': 0, 'total': 0}
        }
        
        for file_path in sampled_files:
            result = check_experiment_file(file_path)
            if result:

                with open(file_path, 'r') as f:
                    data = json.load(f)
                    constraint = data.get('constraint_type')
                    
                if constraint in validation_results:
                    validation_results[constraint]['total'] += 1
                    if result['constraint_satisfied']:
                        validation_results[constraint]['satisfied'] += 1
        

        for constraint in validation_results:
            total = validation_results[constraint]['total']
            if total > 0:
                satisfied = validation_results[constraint]['satisfied']
                validation_results[constraint]['rate'] = satisfied / total * 100
            else:
                validation_results[constraint]['rate'] = 0
        
        return validation_results
    
    def detect_anomalies(self, constraint_stats: Dict) -> Dict:

        anomalies = {
            'warnings': [],
            'critical': [],
            'info': []
        }
        

        if 'lagrangian' in constraint_stats:
            lag_mean = constraint_stats['lagrangian']['mean']
            
            if lag_mean < 0.65:

            elif lag_mean > 0.90:

            else:

        

        if all(k in constraint_stats for k in ['lagrangian', 'ams', 'cpc']):
            lag_mean = constraint_stats['lagrangian']['mean']
            ams_mean = constraint_stats['ams']['mean']
            cpc_mean = constraint_stats['cpc']['mean']
            

            ams_diff = (ams_mean - lag_mean) / lag_mean * 100 if lag_mean > 0 else 0
            cpc_diff = (cpc_mean - lag_mean) / lag_mean * 100 if lag_mean > 0 else 0
            
            if ams_diff > 50:

            elif ams_diff < -30:

            
            if cpc_diff > 50:

            elif cpc_diff < -30:

        

        for constraint in constraint_stats:
            if constraint_stats[constraint]['aa_correct_avg'] < 100:
                anomalies['critical'].append(

                )
        
        return anomalies
    
    def analyze_relative_performance(self, constraint_stats: Dict) -> Dict:
        """分析相对性能"""
        if not all(k in constraint_stats for k in ['lagrangian', 'ams', 'cpc']):
            return {}
        
        lag_mean = constraint_stats['lagrangian']['mean']
        ams_mean = constraint_stats['ams']['mean']
        cpc_mean = constraint_stats['cpc']['mean']
        
        return {
            'ranking': sorted(
                [(k, constraint_stats[k]['mean']) for k in ['lagrangian', 'ams', 'cpc']],
                key=lambda x: x[1]
            ),
            'ams_vs_lag': ((ams_mean - lag_mean) / lag_mean * 100) if lag_mean > 0 else 0,
            'cpc_vs_lag': ((cpc_mean - lag_mean) / lag_mean * 100) if lag_mean > 0 else 0,
            'cpc_vs_ams': ((cpc_mean - ams_mean) / ams_mean * 100) if ams_mean > 0 else 0
        }
    
    def generate_comparison_plot(self, all_results: List[Dict], save_path: Path = None):

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        

        experiments = []
        lag_means = []
        ams_means = []
        cpc_means = []
        
        for result in all_results:
            exp_name = result['directory'][:20]
            experiments.append(exp_name)
            
            stats = result['constraint_stats']
            lag_means.append(stats.get('lagrangian', {}).get('mean', np.nan))
            ams_means.append(stats.get('ams', {}).get('mean', np.nan))
            cpc_means.append(stats.get('cpc', {}).get('mean', np.nan))
        

        ax1 = axes[0, 0]
        x = np.arange(len(experiments))
        width = 0.25
        
        ax1.bar(x - width, lag_means, width, label='Lagrangian', color='#2E86AB')
        ax1.bar(x, ams_means, width, label='AMS', color='#A23B72')
        ax1.bar(x + width, cpc_means, width, label='CPC', color='#F18F01')
        
        ax1.set_xlabel('Experiment')
        ax1.set_ylabel('Mean Accessibility (kcal/mol)')
        ax1.set_title('Constraint Type Performance Comparison')
        ax1.set_xticks(x)
        ax1.set_xticklabels(experiments, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        



        

        ax2 = axes[0, 1]
        relative_data = []
        for result in all_results:
            rel = result.get('relative_performance', {})
            relative_data.append([
                0,  # Lag vs Lag
                rel.get('ams_vs_lag', 0),
                rel.get('cpc_vs_lag', 0)
            ])
        
        if relative_data:
            im = ax2.imshow(relative_data, cmap='RdYlGn_r', aspect='auto', vmin=-50, vmax=50)
            ax2.set_xticks([0, 1, 2])
            ax2.set_xticklabels(['vs Lag', 'AMS vs Lag', 'CPC vs Lag'])
            ax2.set_yticks(range(len(experiments)))
            ax2.set_yticklabels(experiments)
            ax2.set_title('Relative Performance (%)')
            

            for i in range(len(experiments)):
                for j in range(3):
                    value = relative_data[i][j]
                    color = 'white' if abs(value) > 25 else 'black'
                    ax2.text(j, i, f'{value:.0f}%', ha='center', va='center', color=color)
            
            plt.colorbar(im, ax=ax2)
        

        ax3 = axes[1, 0]
        validation_data = []
        for result in all_results:
            val = result.get('constraint_validation', {})
            validation_data.append([
                val.get('lagrangian', {}).get('rate', 0),
                val.get('ams', {}).get('rate', 0),
                val.get('cpc', {}).get('rate', 0)
            ])
        
        if validation_data:
            x = np.arange(len(experiments))
            width = 0.25
            
            ax3.bar(x - width, [d[0] for d in validation_data], width, label='Lagrangian', color='#2E86AB')
            ax3.bar(x, [d[1] for d in validation_data], width, label='AMS', color='#A23B72')
            ax3.bar(x + width, [d[2] for d in validation_data], width, label='CPC', color='#F18F01')
            
            ax3.set_xlabel('Experiment')
            ax3.set_ylabel('Constraint Satisfaction Rate (%)')
            ax3.set_title('Constraint Validation')
            ax3.set_xticks(x)
            ax3.set_xticklabels(experiments, rotation=45, ha='right')
            ax3.legend()
            ax3.axhline(100, color='green', linestyle='--', alpha=0.5)
            ax3.set_ylim(0, 105)
        

        ax4 = axes[1, 1]
        ax4.axis('off')
        

        for i, result in enumerate(all_results):
            anomalies = result.get('anomalies', {})
            exp_name = result['directory'][:25]
            
            summary_text += f"{exp_name}:\n"
            if anomalies.get('critical'):
                summary_text += f"  ❌ {anomalies['critical'][0]}\n"
            elif anomalies.get('warnings'):
                summary_text += f"  ⚠️ {anomalies['warnings'][0]}\n"
            elif anomalies.get('info'):
                summary_text += f"  ✅ {anomalies['info'][0]}\n"
            summary_text += "\n"
        
        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, 
                fontsize=9, verticalalignment='top', fontfamily='monospace')
        
        plt.suptitle('Standard Analysis Report', fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def generate_report(self, all_results: List[Dict]) -> str:
        """生成文本报告"""
        report = []
        report.append("=" * 80)
        report.append("标准分析报告")
        report.append("=" * 80)
        
        for result in all_results:
            report.append(f"\n## {result['directory']}")
            report.append("-" * 40)
            

            basic = result.get('basic_stats', {})
            if basic:
                report.append(f"总实验数: {basic.get('total_experiments', 0)}")
                report.append(f"平均值: {basic.get('mean', 0):.4f} kcal/mol")
                report.append(f"优秀率: {basic.get('excellent_rate', 0):.1f}%")
            

            report.append("\n约束类型性能:")
            for constraint in ['lagrangian', 'ams', 'cpc']:
                stats = result['constraint_stats'].get(constraint, {})
                if stats:
                    report.append(f"  {constraint.upper():12} 平均:{stats['mean']:7.3f}  "
                                f"最佳:{stats['min']:7.3f}  约束满足:{stats['aa_correct_avg']:5.1f}%")
            

            rel = result.get('relative_performance', {})
            if rel:
                report.append(f"\n相对性能:")
                report.append(f"  AMS vs Lagrangian: {rel.get('ams_vs_lag', 0):+.1f}%")
                report.append(f"  CPC vs Lagrangian: {rel.get('cpc_vs_lag', 0):+.1f}%")
            

            anomalies = result.get('anomalies', {})
            if anomalies.get('critical'):
                report.append(f"\n❌ 严重问题: {', '.join(anomalies['critical'])}")
            if anomalies.get('warnings'):
                report.append(f"⚠️ 警告: {', '.join(anomalies['warnings'])}")
            if anomalies.get('info'):
                report.append(f"✅ 正常: {', '.join(anomalies['info'])}")
        
        return "\n".join(report)


def run_standard_analysis(directories: List[Path], output_dir: Path = None):

    analyzer = StandardAnalyzer()
    all_results = []
    

    print("=" * 80)
    

    for dir_path in directories:
        if dir_path.exists():
            result = analyzer.analyze_experiment_batch(dir_path)
            all_results.append(result)
        else:

    
    if not all_results:

        return
    

    print("\n" + "=" * 80)

    print("=" * 80)
    
    report = analyzer.generate_report(all_results)
    print(report)
    

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        

        plot_path = output_dir / 'standard_analysis_plot.png'
        analyzer.generate_comparison_plot(all_results, plot_path)

        

        report_path = output_dir / 'standard_analysis_report.txt'
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)

        

        json_path = output_dir / 'standard_analysis_data.json'
        with open(json_path, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)

    
    return all_results


if __name__ == "__main__":

    default_dirs = [
        Path('results_full/20250909_105804_unified_access_experiments'),
        Path('results_full/20250910_022355_unified_access_experiments'),
        Path('results/20250910_123342_unified_cai_experiments')
    ]
    
    run_standard_analysis(default_dirs, Path('analysis_output'))
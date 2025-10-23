#!/usr/bin/env python3
"""



"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
from datetime import datetime

from .fast_analyzer import FastExperimentAnalyzer


class ExperimentComparator:

    
    def __init__(self):
        """初始化对比器"""
        self.analyzer = FastExperimentAnalyzer()
    
    def load_experiment_data(self, dir_path: Path) -> Optional[pd.DataFrame]:
        """

        
        Args:

            
        Returns:

        """
        results = []
        json_files = list(dir_path.glob("*.json"))
        experiment_files = [f for f in json_files 
                           if not f.name.startswith(('config', 'progress', 'summary'))]
        
        for file_path in experiment_files:
            data = self.analyzer.extract_key_fields(file_path)
            if data:
                results.append(data)
        
        return pd.DataFrame(results) if results else None
    
    def compare_two_experiments(self, df1: pd.DataFrame, df2: pd.DataFrame, 
                               name1: str, name2: str) -> Dict:
        """

        
        Args:




            
        Returns:

        """
        result = {
            'name1': name1,
            'name2': name2,
            'accessibility_comparison': {},
            'constraint_comparison': {},
            'variant_comparison': {},
            'excellence_comparison': {}
        }
        

        result['accessibility_comparison'] = {
            f'{name1}_mean': df1['best_accessibility'].mean(),
            f'{name1}_median': df1['best_accessibility'].median(),
            f'{name1}_min': df1['best_accessibility'].min(),
            f'{name1}_std': df1['best_accessibility'].std(),
            f'{name2}_mean': df2['best_accessibility'].mean(),
            f'{name2}_median': df2['best_accessibility'].median(),
            f'{name2}_min': df2['best_accessibility'].min(),
            f'{name2}_std': df2['best_accessibility'].std(),
            'diff_mean': df1['best_accessibility'].mean() - df2['best_accessibility'].mean(),
            'diff_median': df1['best_accessibility'].median() - df2['best_accessibility'].median(),
            'diff_min': df1['best_accessibility'].min() - df2['best_accessibility'].min(),
            'diff_percent': (df1['best_accessibility'].mean() - df2['best_accessibility'].mean()) / df2['best_accessibility'].mean() * 100
        }
        

        for constraint in ['lagrangian', 'ams', 'cpc']:
            c1_data = df1[df1['constraint_type'] == constraint]['best_accessibility']
            c2_data = df2[df2['constraint_type'] == constraint]['best_accessibility']
            
            if len(c1_data) > 0 and len(c2_data) > 0:
                result['constraint_comparison'][constraint] = {
                    f'{name1}_mean': c1_data.mean(),
                    f'{name2}_mean': c2_data.mean(),
                    'diff': c1_data.mean() - c2_data.mean(),
                    'diff_percent': (c1_data.mean() - c2_data.mean()) / c2_data.mean() * 100
                }
        

        for variant in sorted(set(df1['variant'].unique()) | set(df2['variant'].unique())):
            v1_data = df1[df1['variant'] == variant]['best_accessibility']
            v2_data = df2[df2['variant'] == variant]['best_accessibility']
            
            if len(v1_data) > 0 and len(v2_data) > 0:
                result['variant_comparison'][variant] = {
                    f'{name1}_mean': v1_data.mean(),
                    f'{name2}_mean': v2_data.mean(),
                    'diff': v1_data.mean() - v2_data.mean()
                }
        

        excellent1 = (df1['best_accessibility'] < 1.0).sum()
        excellent2 = (df2['best_accessibility'] < 1.0).sum()
        result['excellence_comparison'] = {
            f'{name1}_excellent': excellent1,
            f'{name1}_total': len(df1),
            f'{name1}_rate': excellent1 / len(df1) * 100,
            f'{name2}_excellent': excellent2,
            f'{name2}_total': len(df2),
            f'{name2}_rate': excellent2 / len(df2) * 100
        }
        
        return result
    
    def compare_cai_vs_access(self, cai_dir: Path, access_dir: Path) -> Dict:
        """

        
        Args:


            
        Returns:

        """

        df_cai = self.load_experiment_data(cai_dir)
        df_access = self.load_experiment_data(access_dir)
        
        if df_cai is None or df_access is None:
            return {'error': '无法加载实验数据'}
        

        result = self.compare_two_experiments(df_cai, df_access, 'CAI', 'Access')
        

        if 'final_ecai' in df_cai.columns and df_cai['final_ecai'].notna().any():
            ecai_data = df_cai['final_ecai'].dropna()
            result['cai_metrics'] = {
                'mean_cai': ecai_data.mean(),
                'median_cai': ecai_data.median(),
                'min_cai': ecai_data.min(),
                'max_cai': ecai_data.max()
            }
            
            if 'cai_target_achieved' in df_cai.columns:
                achieved = df_cai['cai_target_achieved'].sum()
                result['cai_metrics']['target_achieved'] = achieved
                result['cai_metrics']['target_achieved_rate'] = achieved / len(df_cai) * 100
        

        diff_mean = result['accessibility_comparison']['diff_mean']
        if abs(diff_mean) < 0.01:
            result['impact_analysis'] = 'CAI优化对Accessibility几乎没有影响'
        elif diff_mean > 0:
            result['impact_analysis'] = f'CAI优化降低了Accessibility性能（增加{diff_mean:.4f} kcal/mol）'
        else:
            result['impact_analysis'] = f'CAI优化提升了Accessibility性能（减少{-diff_mean:.4f} kcal/mol）'
        
        return result
    
    def compare_access_timeline(self, experiment_dirs: List[Tuple[str, str, Path]]) -> Dict:
        """

        
        Args:

            
        Returns:

        """
        timeline_data = []
        
        for exp_name, exp_time, exp_dir in experiment_dirs:
            df = self.load_experiment_data(exp_dir)
            if df is not None:
                timeline_data.append({
                    'name': exp_name,
                    'time': exp_time,
                    'mean': df['best_accessibility'].mean(),
                    'median': df['best_accessibility'].median(),
                    'min': df['best_accessibility'].min(),
                    'std': df['best_accessibility'].std(),
                    'excellent_rate': (df['best_accessibility'] < 1.0).mean() * 100,
                    'data': df
                })
        
        if len(timeline_data) < 2:
            return {'error': '需要至少2个实验进行时间序列对比'}
        

        result = {
            'timeline': timeline_data,
            'performance_changes': []
        }
        
        for i in range(1, len(timeline_data)):
            prev = timeline_data[i-1]
            curr = timeline_data[i]
            
            change = {
                'from': prev['name'],
                'to': curr['name'],
                'time_gap': f"{prev['time']} → {curr['time']}",
                'mean_change': curr['mean'] - prev['mean'],
                'mean_change_percent': (curr['mean'] - prev['mean']) / prev['mean'] * 100,
                'median_change': curr['median'] - prev['median'],
                'min_change': curr['min'] - prev['min'],
                'excellence_change': curr['excellent_rate'] - prev['excellent_rate']
            }
            

            if change['mean_change'] > 0.1:
                change['trend'] = '性能恶化'
            elif change['mean_change'] < -0.1:
                change['trend'] = '性能改善'
            else:
                change['trend'] = '基本稳定'
            
            result['performance_changes'].append(change)
        

        result['constraint_trends'] = {}
        for constraint in ['lagrangian', 'ams', 'cpc']:
            constraint_timeline = []
            for exp_data in timeline_data:
                constraint_data = exp_data['data'][exp_data['data']['constraint_type'] == constraint]
                if len(constraint_data) > 0:
                    constraint_timeline.append({
                        'name': exp_data['name'],
                        'mean': constraint_data['best_accessibility'].mean()
                    })
            
            if len(constraint_timeline) > 1:
                result['constraint_trends'][constraint] = constraint_timeline
        
        return result
    
    def print_comparison_report(self, comparison: Dict):


        print("=" * 80)
        

        acc = comparison['accessibility_comparison']

        print(f"  {comparison['name1']}:")
        name1_mean_key = f"{comparison['name1']}_mean"
        name1_median_key = f"{comparison['name1']}_median"
        name1_min_key = f"{comparison['name1']}_min"




        print(f"\n  {comparison['name2']}:")
        name2_mean_key = f"{comparison['name2']}_mean"
        name2_median_key = f"{comparison['name2']}_median"
        name2_min_key = f"{comparison['name2']}_min"



        

        diff_mean = acc['diff_mean']


        

        if comparison['constraint_comparison']:

            for constraint, data in comparison['constraint_comparison'].items():
                name1_mean_key = f"{comparison['name1']}_mean"
                name2_mean_key = f"{comparison['name2']}_mean"
                print(f"  {constraint.upper()}: "
                      f"{comparison['name1']}={data[name1_mean_key]:.4f}, "
                      f"{comparison['name2']}={data[name2_mean_key]:.4f}, "

        

        exc = comparison['excellence_comparison']

        name1_excellent_key = f"{comparison['name1']}_excellent"
        name1_total_key = f"{comparison['name1']}_total"
        name1_rate_key = f"{comparison['name1']}_rate"
        name2_excellent_key = f"{comparison['name2']}_excellent"
        name2_total_key = f"{comparison['name2']}_total"
        name2_rate_key = f"{comparison['name2']}_rate"

        print(f"  • {comparison['name1']}: {exc[name1_excellent_key]}/{exc[name1_total_key]} "
              f"({exc[name1_rate_key]:.1f}%)")
        print(f"  • {comparison['name2']}: {exc[name2_excellent_key]}/{exc[name2_total_key]} "
              f"({exc[name2_rate_key]:.1f}%)")
        

        if 'cai_metrics' in comparison:

            metrics = comparison['cai_metrics']


            if 'target_achieved' in metrics:

        

        if 'impact_analysis' in comparison:

    
    def print_timeline_report(self, timeline_result: Dict):
        """打印时间序列对比报告"""
        if 'error' in timeline_result:
            print(f"❌ {timeline_result['error']}")
            return
        
        print("\n📈 性能时间序列分析")
        print("=" * 80)
        

        print("\n📊 各时间点性能:")
        print(f"{'实验':<25} {'时间':<20} {'平均值':>10} {'中位数':>10} {'最优值':>10} {'优秀率':>10}")
        print("-" * 85)
        
        for exp in timeline_result['timeline']:
            print(f"{exp['name']:<25} {exp['time']:<20} "
                  f"{exp['mean']:>10.4f} {exp['median']:>10.4f} "
                  f"{exp['min']:>10.4f} {exp['excellent_rate']:>9.1f}%")
        

        print("\n📉 性能变化分析:")
        for change in timeline_result['performance_changes']:
            print(f"\n  {change['time_gap']}:")
            print(f"    • 平均值变化: {change['mean_change']:+.4f} kcal/mol "
                  f"({change['mean_change_percent']:+.1f}%)")
            print(f"    • 优秀率变化: {change['excellence_change']:+.1f}%")
            print(f"    • 趋势: {change['trend']}")
        

        if 'constraint_trends' in timeline_result:
            print("\n📊 按约束类型的性能趋势:")
            for constraint, data in timeline_result['constraint_trends'].items():
                print(f"\n  {constraint.upper()}:")
                for exp in data:
                    print(f"    • {exp['name']}: {exp['mean']:.4f}")
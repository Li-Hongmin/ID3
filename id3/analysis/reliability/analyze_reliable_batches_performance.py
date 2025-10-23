#!/usr/bin/env python3
"""


"""

import json
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def analyze_batch_performance(batch_dir: Path) -> Dict:

    results = {
        'batch_name': batch_dir.name,
        'lagrangian': {'files': [], 'accessibility': [], 'final_accessibility': []},
        'ams': {'files': [], 'accessibility': [], 'final_accessibility': []},
        'cpc': {'files': [], 'accessibility': [], 'final_accessibility': []}
    }
    

    for method in ['lagrangian', 'ams', 'cpc']:
        method_files = list(batch_dir.glob(f"*{method}*.json"))
        results[method]['file_count'] = len(method_files)
        

        sample_files = method_files[:50] if len(method_files) > 50 else method_files
        
        for file_path in sample_files:
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                best_acc = data.get('best_accessibility', None)
                final_acc = data.get('final_accessibility', None)
                
                if best_acc is not None:
                    results[method]['accessibility'].append(best_acc)
                    results[method]['files'].append(file_path.name)
                
                if final_acc is not None:
                    results[method]['final_accessibility'].append(final_acc)
                    
            except Exception as e:

                continue
    
    return results

def calculate_statistics(values: List[float]) -> Dict:
    """计算统计数据"""
    if not values:
        return {
            'count': 0,
            'mean': None,
            'median': None,
            'std': None,
            'min': None,
            'max': None,
            'p99': None,
            'p99.9': None,
            'p99.99': None,
            'excellent_rate': None,  # < 1.0
            'good_rate': None,       # < 1.5
            'acceptable_rate': None  # < 2.0
        }
    
    arr = np.array(values)
    
    return {
        'count': len(values),
        'mean': np.mean(arr),
        'median': np.median(arr),
        'std': np.std(arr),
        'min': np.min(arr),
        'max': np.max(arr),
        'p99': np.percentile(arr, 99) if len(arr) >= 100 else np.max(arr),
        'p99.9': np.percentile(arr, 99.9) if len(arr) >= 1000 else np.max(arr),
        'p99.99': np.percentile(arr, 99.99) if len(arr) >= 10000 else np.max(arr),
        'excellent_rate': np.mean(arr < 1.0) * 100,  # < 1.0 kcal/mol
        'good_rate': np.mean(arr < 1.5) * 100,       # < 1.5 kcal/mol  
        'acceptable_rate': np.mean(arr < 2.0) * 100  # < 2.0 kcal/mol
    }

def generate_comparison_table(batch_results: List[Dict]) -> pd.DataFrame:

    comparison_data = []
    
    for batch in batch_results:
        batch_name = batch['batch_name']
        
        for method in ['lagrangian', 'ams', 'cpc']:
            method_data = batch[method]
            

            best_stats = calculate_statistics(method_data['accessibility'])
            

            final_stats = calculate_statistics(method_data['final_accessibility'])
            
            comparison_data.append({



                






                






            })
    
    return pd.DataFrame(comparison_data)

def check_parameter_consistency(batch_results: List[Dict]) -> Dict:
    """检查参数一致性和结果重现性"""
    consistency_report = {
        'parameter_analysis': {},
        'result_consistency': {},
        'recommendations': []
    }
    

    cai_batches = [b for b in batch_results if 'cai' in b['batch_name'].lower()]
    access_batches = [b for b in batch_results if 'access' in b['batch_name'].lower()]
    
    consistency_report['parameter_analysis'] = {
        'cai_batches': len(cai_batches),
        'access_batches': len(access_batches),
        'total_batches': len(batch_results)
    }
    

    if len(access_batches) >= 2:
        batch1 = access_batches[0]
        batch2 = access_batches[1]
        
        for method in ['lagrangian', 'ams', 'cpc']:
            stats1 = calculate_statistics(batch1[method]['accessibility'])
            stats2 = calculate_statistics(batch2[method]['accessibility'])
            
            if stats1['mean'] and stats2['mean']:
                diff = abs(stats1['mean'] - stats2['mean'])
                relative_diff = diff / stats1['mean'] * 100
                
                consistency_report['result_consistency'][f'{method}_mean_diff'] = {
                    'absolute_diff': diff,
                    'relative_diff_percent': relative_diff,
                    'consistent': relative_diff < 10.0
                }
    
    return consistency_report

def main():
    reliable_dir = Path("可靠实验结果")
    
    if not reliable_dir.exists():
        print("❌ 可靠实验结果目录不存在")
        return
    
    print("🔍 分析三个可信批次的详细性能")
    print("=" * 80)
    
    batch_results = []
    

    for batch_dir in sorted(reliable_dir.iterdir()):
        if batch_dir.is_dir():
            print(f"\n📊 分析批次: {batch_dir.name}")
            
            result = analyze_batch_performance(batch_dir)
            batch_results.append(result)
            

            for method in ['lagrangian', 'ams', 'cpc']:
                count = result[method]['file_count']
                if count > 0:
                    stats = calculate_statistics(result[method]['accessibility'])
                    print(f"   {method.upper()}: {count}个文件, 最佳={stats['min']:.3f}, 均值={stats['mean']:.3f}")
    

    print(f"\n📈 三个批次性能对比表格")
    print("=" * 80)
    
    comparison_df = generate_comparison_table(batch_results)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_colwidth', 15)
    print(comparison_df.to_string(index=False))
    

    print(f"\n🔍 参数一致性和结果重现性分析")
    print("=" * 80)
    
    consistency = check_parameter_consistency(batch_results)
    
    print(f"批次组成:")
    print(f"   CAI实验批次: {consistency['parameter_analysis']['cai_batches']} 个")
    print(f"   Access实验批次: {consistency['parameter_analysis']['access_batches']} 个")
    print(f"   总批次数: {consistency['parameter_analysis']['total_batches']} 个")
    
    if consistency['result_consistency']:
        print(f"\n两个Access批次结果一致性:")
        for method_key, consistency_data in consistency['result_consistency'].items():
            method_name = method_key.replace('_mean_diff', '').upper()
            is_consistent = "✅ 一致" if consistency_data['consistent'] else "❌ 不一致"
            print(f"   {method_name}: 相对差异 {consistency_data['relative_diff_percent']:.1f}% {is_consistent}")
    

    print(f"\n🎯 性能等级评估 (基于最佳值)")
    print("=" * 80)
    
    for batch in batch_results:
        print(f"\n批次: {batch['batch_name']}")
        for method in ['lagrangian', 'ams', 'cpc']:
            stats = calculate_statistics(batch[method]['accessibility'])
            if stats['min'] is not None:
                if stats['min'] < 1.0:
                    level = "🌟 优秀"
                elif stats['min'] < 1.5:
                    level = "✅ 良好"
                elif stats['min'] < 2.0:
                    level = "⚠️ 可接受"
                else:
                    level = "❌ 需改进"
                
                print(f"   {method.upper()}: {stats['min']:.3f} kcal/mol {level}")
    
    return batch_results

if __name__ == "__main__":
    main()
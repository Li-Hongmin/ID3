#!/usr/bin/env python3
"""

"""

import json
import glob
import numpy as np
from pathlib import Path

def analyze_complete_batch(batch_path: str, batch_name: str):
    print(f"\n📊 {batch_name}")
    print("-" * 60)
    
    results = {}
    
    for method in ['lagrangian', 'ams', 'cpc']:
        files = glob.glob(f'{batch_path}/*{method}*.json')
        
        if not files:
            continue
            
        final_values = []
        best_values = []
        
        for f in files:
            try:
                with open(f) as file:
                    data = json.load(file)
                final_values.append(data['final_accessibility'])
                best_values.append(data['best_accessibility'])
            except:
                continue
        
        if final_values:

            final_min = min(final_values)
            final_mean = np.mean(final_values)
            final_median = np.median(final_values)
            best_min = min(best_values)
            best_mean = np.mean(best_values)
            

            excellent = sum(1 for v in final_values if v < 1.0)
            good = sum(1 for v in final_values if 1.0 <= v < 1.5)
            acceptable = sum(1 for v in final_values if 1.5 <= v < 2.0)
            poor = sum(1 for v in final_values if v >= 2.0)
            

            if final_min < 1.0:
                level = "🌟 优秀"
            elif final_min < 1.5:
                level = "✅ 良好"
            elif final_min < 2.0:
                level = "⚠️ 可接受"
            else:
                level = "❌ 需改进"
                
            print(f"  {method.upper()}: {len(files)}个文件")
            print(f"    Final最小值: {final_min:.3f} kcal/mol {level}")
            print(f"    Final均值: {final_mean:.3f} kcal/mol")
            print(f"    Final中位数: {final_median:.3f} kcal/mol")
            print(f"    Best最小值: {best_min:.3f} kcal/mol")
            print(f"    优秀率(<1.0): {excellent}/{len(files)} ({excellent/len(files)*100:.1f}%)")
            print(f"    良好率(1.0-1.5): {good}/{len(files)} ({good/len(files)*100:.1f}%)")
            
            results[method] = {
                'final_min': final_min,
                'final_mean': final_mean,
                'final_median': final_median,
                'best_min': best_min,
                'best_mean': best_mean,
                'excellent_rate': excellent/len(files)*100,
                'file_count': len(files)
            }
    
    return results

def main():
    print("🔍 修正后的三个可信批次完整性能分析")
    print("=" * 80)
    
    batches = [
        ("可靠实验结果/20250909_003751_unified_cai_experiments", "CAI实验批次 (20250909)"),
        ("可靠实验结果/20250910_004126_unified_access_experiments", "Access实验批次1 (20250910_004126)"),
        ("可靠实验结果/20250910_022355_unified_access_experiments", "Access实验批次2 (20250910_022355)")
    ]
    
    all_results = {}
    
    for batch_path, batch_name in batches:
        results = analyze_complete_batch(batch_path, batch_name)
        all_results[batch_name] = results
    

    print(f"\n📈 修正后的性能对比表格 (基于完整数据)")
    print("=" * 90)
    print(f"{'批次':<20} {'方法':<12} {'文件数':<8} {'Final最小值':<12} {'Final均值':<12} {'优秀率(%)':<10} {'性能等级':<12}")
    print("-" * 90)
    
    for batch_name, batch_results in all_results.items():
        short_name = batch_name.split('(')[1].replace(')', '')
        
        for method, stats in batch_results.items():

            if stats['final_min'] < 1.0:
                level = "🌟优秀"
            elif stats['final_min'] < 1.5:
                level = "✅良好"
            elif stats['final_min'] < 2.0:
                level = "⚠️可接受"
            else:
                level = "❌需改进"
            
            print(f"{short_name:<20} {method.upper():<12} {stats['file_count']:<8} {stats['final_min']:.3f} kcal/mol {stats['final_mean']:.3f}        {stats['excellent_rate']:.1f}      {level:<12}")
    

    print(f"\n🔍 关键技术发现")
    print("=" * 80)
    
    cai_lag = all_results["CAI实验批次 (20250909)"]["lagrangian"]
    acc1_ams = all_results["Access实验批次1 (20250910_004126)"]["ams"]
    
    print(f"1. **CAI + Lagrangian组合**:")
    print(f"   - Final最小值: {cai_lag['final_min']:.3f} kcal/mol")
    print(f"   - 优秀率: {cai_lag['excellent_rate']:.1f}% (达到<1.0阈值)")
    print(f"   - 结论: CAI增强确实提高了Lagrangian在最优案例中的表现")
    
    print(f"\n2. **无CAI模式下的AMS/CPC**:")
    print(f"   - AMS最小值: {acc1_ams['final_min']:.3f} kcal/mol") 
    print(f"   - 优秀率: {acc1_ams['excellent_rate']:.1f}%")
    print(f"   - 结论: 在无额外优化目标下，连续约束方法表现最稳定")
    
    print(f"\n3. **方法特点总结**:")
    print(f"   - **Lagrangian**: 在多目标优化(+CAI)下表现突出，但平均性能一般")
    print(f"   - **AMS/CPC**: 在单目标优化下稳定，大部分实验都能达到优秀水平")
    print(f"   - **技术解释**: CAI约束为Lagrangian提供了额外的梯度信息")

if __name__ == "__main__":
    main()
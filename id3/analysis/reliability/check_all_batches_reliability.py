#!/usr/bin/env python3
"""


"""

import json
import sys
import os
from pathlib import Path
import glob

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def translate_rna_to_protein(rna_sequence: str) -> str:

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
    
    protein = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            protein += codon_table.get(codon, 'X')
    return protein

def check_batch_reliability(batch_dir):
    """检查单个批次的可信性"""
    batch_path = Path(batch_dir)
    batch_name = batch_path.name
    

    if 'INVALID' in batch_name or 'PROBLEMATIC' in batch_name:
        return {
            'batch': batch_name,
            'status': 'INVALID',
            'reason': 'Already marked as invalid',
            'lagrangian_files': 0,
            'violated_files': 0
        }
    
    print(f"\n🔍 检查批次: {batch_name}")
    

    lagrangian_files = list(batch_path.glob("*lagrangian*.json"))
    
    if not lagrangian_files:
        return {
            'batch': batch_name,
            'status': 'NO_LAGRANGIAN',
            'reason': 'No Lagrangian files found',
            'lagrangian_files': 0,
            'violated_files': 0
        }
    
    print(f"   发现 {len(lagrangian_files)} 个Lagrangian文件")
    

    violated_count = 0
    total_checked = 0
    sample_violations = []
    
    for file_path in lagrangian_files[:10]:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            

            best_accessibility = data.get('best_accessibility')
            expected_aa = data.get('expected_amino_acids', '')
            trajectory = data.get('trajectory', {})
            accessibility_values = trajectory.get('accessibility', [])
            sequences = trajectory.get('discrete_sequences', [])
            
            if not accessibility_values or not sequences or not expected_aa:
                continue
            

            best_positions = [i for i, v in enumerate(accessibility_values) if abs(v - best_accessibility) < 0.0001]
            
            if not best_positions or best_positions[0] >= len(sequences):
                continue
            

            best_iteration = best_positions[0]
            best_sequence = sequences[best_iteration]
            
            if best_sequence and len(best_sequence) == len(expected_aa) * 3:
                best_translated = translate_rna_to_protein(best_sequence)
                constraint_satisfied = (best_translated == expected_aa)
                
                total_checked += 1
                if not constraint_satisfied:
                    violated_count += 1
                    if len(sample_violations) < 3:
                        sample_violations.append({
                            'file': file_path.name,
                            'best_accessibility': best_accessibility,
                            'iteration': best_iteration,
                            'expected': expected_aa[:20] + '...' if len(expected_aa) > 20 else expected_aa,
                            'actual': best_translated[:20] + '...' if len(best_translated) > 20 else best_translated
                        })
                        
        except Exception as e:
            print(f"   ❌ 处理文件 {file_path.name} 时出错: {e}")
            continue
    
    if total_checked == 0:
        return {
            'batch': batch_name,
            'status': 'ERROR',
            'reason': 'Could not check any files',
            'lagrangian_files': len(lagrangian_files),
            'violated_files': 0
        }
    
    violation_rate = violated_count / total_checked
    

    if violation_rate == 0:
        status = 'RELIABLE'
        reason = 'All sampled best_accessibility values satisfy constraints'
    elif violation_rate < 0.3:
        status = 'PARTIALLY_RELIABLE'  
        reason = f'Low violation rate: {violation_rate:.1%}'
    else:
        status = 'UNRELIABLE'
        reason = f'High violation rate: {violation_rate:.1%}'
    
    result = {
        'batch': batch_name,
        'status': status,
        'reason': reason,
        'lagrangian_files': len(lagrangian_files),
        'checked_files': total_checked,
        'violated_files': violated_count,
        'violation_rate': violation_rate,
        'sample_violations': sample_violations
    }
    
    print(f"   结果: {status} ({violation_rate:.1%} violation rate)")
    
    return result

def main():
    results_dir = Path("results_full")
    
    if not results_dir.exists():
        print("❌ results_full目录不存在")
        return
    
    print("🔍 系统性检查所有实验批次的可信性")
    print("=" * 80)
    
    all_results = []
    

    for batch_dir in sorted(results_dir.iterdir()):
        if batch_dir.is_dir():
            result = check_batch_reliability(batch_dir)
            all_results.append(result)
    

    print(f"\n📊 可信性报告:")
    print("=" * 80)
    
    reliable_batches = []
    unreliable_batches = []
    partial_batches = []
    
    for result in all_results:
        status = result['status']
        batch = result['batch']
        
        if status == 'RELIABLE':
            reliable_batches.append(result)
            print(f"✅ {batch}: 可信 ({result['lagrangian_files']} Lagrangian文件)")
            
        elif status == 'UNRELIABLE':
            unreliable_batches.append(result)
            print(f"❌ {batch}: 不可信 ({result.get('violation_rate', 0):.1%} 违反率)")
            
        elif status == 'PARTIALLY_RELIABLE':
            partial_batches.append(result)
            print(f"⚠️ {batch}: 部分可信 ({result.get('violation_rate', 0):.1%} 违反率)")
            
        elif status == 'NO_LAGRANGIAN':
            print(f"ℹ️ {batch}: 无Lagrangian文件 (可能只有AMS/CPC)")
            
        elif status == 'INVALID':
            print(f"🚫 {batch}: 已标记无效")
        
        else:
            print(f"❓ {batch}: 检查出错")
    

    print(f"\n🎯 总结:")
    print(f"   可信批次: {len(reliable_batches)}")
    print(f"   不可信批次: {len(unreliable_batches)}")  
    print(f"   部分可信批次: {len(partial_batches)}")
    
    if reliable_batches:
        print(f"\n✅ 推荐使用的可信批次:")
        for batch in reliable_batches:
            print(f"   - {batch['batch']}")
    
    if unreliable_batches:
        print(f"\n❌ 需要隔离的不可信批次:")
        for batch in unreliable_batches:
            print(f"   - {batch['batch']} ({batch.get('violation_rate', 0):.1%} 违反率)")
    
    return all_results

if __name__ == "__main__":
    main()
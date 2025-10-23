#!/usr/bin/env python3
"""


"""

import json
import sys
import os
from pathlib import Path

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

def check_detailed_steps_file(file_path):
    """检查单个detailed_steps文件的可信性"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        return {
            'file': file_path.name,
            'status': 'ERROR',
            'error': f'Failed to load JSON: {e}'
        }
    

    if 'step_records' not in data:
        return {
            'file': file_path.name,
            'status': 'ERROR',
            'error': 'No step_records found'
        }
    
    step_records = data['step_records']
    metadata = data.get('metadata', {})
    
    if not step_records:
        return {
            'file': file_path.name,
            'status': 'ERROR', 
            'error': 'Empty step_records'
        }
    

    expected_aa = metadata.get('amino_acid_sequence', '')
    if not expected_aa:
        return {
            'file': file_path.name,
            'status': 'ERROR',
            'error': 'No amino_acid_sequence in metadata'
        }
    

    accessibility_values = []
    sequences = []
    
    for step in step_records:
        if 'true_accessibility' in step:
            accessibility_values.append(step['true_accessibility'])
        elif 'accessibility' in step:
            accessibility_values.append(step['accessibility'])
        else:
            continue
            

        if 'rna_sequence_discrete' in step:
            sequences.append(step['rna_sequence_discrete'])
        elif 'discrete_sequence' in step:
            sequences.append(step['discrete_sequence'])
        else:
            sequences.append(None)
    
    if not accessibility_values:
        return {
            'file': file_path.name,
            'status': 'ERROR',
            'error': 'No accessibility values found'
        }
    

    best_accessibility = min(accessibility_values)
    best_positions = [i for i, v in enumerate(accessibility_values) if abs(v - best_accessibility) < 0.0001]
    
    if not best_positions or best_positions[0] >= len(sequences):
        return {
            'file': file_path.name,
            'status': 'UNKNOWN',
            'reason': 'Cannot find corresponding sequence for best accessibility'
        }
    
    best_iteration = best_positions[0]
    best_sequence = sequences[best_iteration]
    
    if not best_sequence or len(best_sequence) != len(expected_aa) * 3:
        return {
            'file': file_path.name,
            'status': 'UNKNOWN',
            'reason': 'Invalid sequence length or missing sequence'
        }
    

    best_translated = translate_rna_to_protein(best_sequence)
    constraint_satisfied = (best_translated == expected_aa)
    

    final_sequence = sequences[-1] if sequences else None
    final_translated = translate_rna_to_protein(final_sequence) if final_sequence else ''
    final_constraint_satisfied = (final_translated == expected_aa)
    
    status = 'RELIABLE' if constraint_satisfied else 'UNRELIABLE'
    
    return {
        'file': file_path.name,
        'status': status,
        'best_accessibility': best_accessibility,
        'best_iteration': best_iteration,
        'total_iterations': len(accessibility_values),
        'best_position_pct': best_iteration / len(accessibility_values) * 100,
        'constraint_satisfied': constraint_satisfied,
        'final_constraint_satisfied': final_constraint_satisfied,
        'expected_aa_length': len(expected_aa),
        'best_sequence_length': len(best_sequence),
        'protein_name': file_path.name.split('_')[0],
    }

def main():
    detailed_steps_dir = Path("external_storage/archives/results/id3_deepraccess/detailed_steps")
    
    if not detailed_steps_dir.exists():
        print(f"❌ 目录不存在: {detailed_steps_dir}")
        return
    
    print("🔍 检查external_storage中Lagrangian详细步骤文件的可信性")
    print("=" * 80)
    

    lagrangian_files = list(detailed_steps_dir.glob("*ID3-L*.json"))
    
    if not lagrangian_files:
        print("❌ 未找到Lagrangian文件")
        return
    
    print(f"发现 {len(lagrangian_files)} 个Lagrangian详细步骤文件")
    

    sample_files = lagrangian_files[:20]
    print(f"检查前 {len(sample_files)} 个文件作为样本...")
    
    reliable_count = 0
    unreliable_count = 0
    error_count = 0
    unknown_count = 0
    
    results = []
    
    for file_path in sample_files:
        print(f"\n📄 检查: {file_path.name}")
        
        result = check_detailed_steps_file(file_path)
        results.append(result)
        
        status = result['status']
        if status == 'RELIABLE':
            reliable_count += 1
            print(f"   ✅ 可信: best={result['best_accessibility']:.4f}, 第{result['best_iteration']}次迭代")
        elif status == 'UNRELIABLE':
            unreliable_count += 1
            print(f"   ❌ 不可信: best={result['best_accessibility']:.4f}, 约束违反")
        elif status == 'ERROR':
            error_count += 1
            print(f"   💥 错误: {result.get('error', 'Unknown error')}")
        else:
            unknown_count += 1
            print(f"   ❓ 未知: {result.get('reason', 'Unknown reason')}")
    

    print(f"\n📊 检查结果总结:")
    print("=" * 80)
    print(f"   检查文件数: {len(sample_files)}")
    print(f"   可信文件: {reliable_count} ({reliable_count/len(sample_files)*100:.1f}%)")
    print(f"   不可信文件: {unreliable_count} ({unreliable_count/len(sample_files)*100:.1f}%)")
    print(f"   错误文件: {error_count} ({error_count/len(sample_files)*100:.1f}%)")
    print(f"   未知状态: {unknown_count} ({unknown_count/len(sample_files)*100:.1f}%)")
    

    protein_stats = {}
    for result in results:
        if result['status'] in ['RELIABLE', 'UNRELIABLE']:
            protein = result.get('protein_name', 'Unknown')
            if protein not in protein_stats:
                protein_stats[protein] = {'reliable': 0, 'unreliable': 0}
            if result['status'] == 'RELIABLE':
                protein_stats[protein]['reliable'] += 1
            else:
                protein_stats[protein]['unreliable'] += 1
    
    if protein_stats:
        print(f"\n🧬 按蛋白质分组统计:")
        print("-" * 40)
        for protein, stats in protein_stats.items():
            total = stats['reliable'] + stats['unreliable']
            reliable_pct = stats['reliable'] / total * 100 if total > 0 else 0
            print(f"   {protein}: {stats['reliable']}/{total} 可信 ({reliable_pct:.1f}%)")
    

    if len(sample_files) > 0:
        overall_reliable_rate = reliable_count / len(sample_files)
        if overall_reliable_rate >= 0.9:
            conclusion = "✅ 整体可信"
        elif overall_reliable_rate >= 0.7:
            conclusion = "⚠️ 大部分可信"
        elif overall_reliable_rate >= 0.3:
            conclusion = "❓ 部分可信"
        else:
            conclusion = "❌ 大部分不可信"
        
        print(f"\n🎯 整体结论: {conclusion} ({overall_reliable_rate:.1%} 可信率)")
    
    return results

if __name__ == "__main__":
    main()
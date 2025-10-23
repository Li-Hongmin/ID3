#!/usr/bin/env python3
"""

"""

import json
import sys
import os
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

def verify_best_sequence_constraint(file_path: str):
    """验证best时刻的序列约束满足情况"""
    with open(file_path, 'r') as f:
        data = json.load(f)
    

    best_accessibility = data['best_accessibility']
    expected_aa = data['expected_amino_acids']
    trajectory = data['trajectory']
    accessibility_values = trajectory['accessibility']
    sequences = trajectory.get('discrete_sequences', [])
    

    best_positions = [i for i, v in enumerate(accessibility_values) if abs(v - best_accessibility) < 0.0001]
    
    if not best_positions:
        return {
            'error': 'Could not find best_accessibility in trajectory'
        }
    
    best_iteration = best_positions[0]
    

    if best_iteration >= len(sequences):
        return {
            'error': f'No sequence available for best iteration {best_iteration}'
        }
    
    best_sequence = sequences[best_iteration]
    best_translated = translate_rna_to_protein(best_sequence)
    

    constraint_satisfied = (best_translated == expected_aa)
    

    final_sequence = data['final_sequence']
    final_translated = translate_rna_to_protein(final_sequence)
    final_constraint_satisfied = (final_translated == expected_aa)
    
    return {
        'file': file_path.split('/')[-1],
        'best_accessibility': best_accessibility,
        'best_iteration': best_iteration,
        'total_iterations': len(accessibility_values),
        'best_position_pct': best_iteration / len(accessibility_values) * 100,
        
        'best_sequence_length': len(best_sequence),
        'best_translated_aa': best_translated,
        'expected_aa': expected_aa,
        'best_constraint_satisfied': constraint_satisfied,
        
        'final_translated_aa': final_translated,
        'final_constraint_satisfied': final_constraint_satisfied,
        
        'sequences_match': best_sequence == final_sequence
    }

def main():
    import glob
    

    pattern = "results_full/20250909_105804_unified_access_experiments/*lagrangian*.json"
    test_files = sorted(glob.glob(pattern))[:3]
    
    print("🔍 验证best_accessibility时刻的约束满足情况")
    print("=" * 80)
    
    constraint_violated_count = 0
    total_files = 0
    
    for test_file in test_files:
        try:
            result = verify_best_sequence_constraint(test_file)
            
            if 'error' in result:
                print(f"❌ 错误 {result['file']}: {result['error']}")
                continue
            
            total_files += 1
            if not result['best_constraint_satisfied']:
                constraint_violated_count += 1
            
            print(f"\n📋 {result['file']}")
            print(f"🎯 Best: {result['best_accessibility']:.4f} (第{result['best_iteration']}次, {result['best_position_pct']:.1f}%)")
            print(f"✅ Best约束满足: {result['best_constraint_satisfied']}")
            print(f"✅ Final约束满足: {result['final_constraint_satisfied']}")
            print(f"📊 序列相同: {result['sequences_match']}")
                
        except Exception as e:
            print(f"❌ 分析失败 {test_file}: {e}")
    
    print(f"\n🎯 总结:")
    print(f"   检查文件数: {total_files}")
    print(f"   Best时刻违反约束: {constraint_violated_count}")
    print(f"   约束违反率: {constraint_violated_count/total_files*100:.1f}%")
    
    if constraint_violated_count > 0:
        print(f"🚨 结论: best_accessibility记录了违反约束的序列，这些值不可信！")

if __name__ == "__main__":
    main()
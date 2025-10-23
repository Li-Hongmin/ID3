#!/usr/bin/env python3
"""


"""

import json
from pathlib import Path
import random

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

def check_experiment_file(file_path: Path):
    """检查单个实验文件的约束满足情况"""
    with open(file_path, 'r') as f:

        content = f.read(100000)
        

        import re
        

        target_aa_match = re.search(r'"expected_amino_acids":\s*"([^"]+)"', content)
        actual_aa_match = re.search(r'"actual_amino_acids":\s*"([^"]+)"', content)
        final_seq_match = re.search(r'"final_sequence":\s*"([^"]+)"', content)
        best_acc_match = re.search(r'"best_accessibility":\s*([\d.]+)', content)
        aa_correct_match = re.search(r'"amino_acids_correct":\s*([\d.]+)', content)
        
        if not all([target_aa_match, best_acc_match]):
            return None
            
        target_aa = target_aa_match.group(1)
        best_acc = float(best_acc_match.group(1))
        actual_aa = actual_aa_match.group(1) if actual_aa_match else None
        final_seq = final_seq_match.group(1) if final_seq_match else None
        aa_correct = float(aa_correct_match.group(1)) if aa_correct_match else None
        

        translated_aa = None
        if final_seq:
            translated_aa = translate_rna_to_protein(final_seq)
        
        return {
            'file': file_path.name,
            'target_aa': target_aa[:20] + '...' if len(target_aa) > 20 else target_aa,
            'actual_aa': actual_aa[:20] + '...' if actual_aa and len(actual_aa) > 20 else actual_aa,
            'translated_aa': translated_aa[:20] + '...' if translated_aa and len(translated_aa) > 20 else translated_aa,
            'best_accessibility': best_acc,
            'aa_correct': aa_correct,
            'constraint_satisfied': target_aa == actual_aa if actual_aa else None,
            'translation_correct': target_aa == translated_aa if translated_aa else None
        }

def analyze_directory(dir_path: Path, sample_size: int = 5):


    print("=" * 80)
    
    json_files = list(dir_path.glob("*.json"))
    experiment_files = [f for f in json_files 
                        if not f.name.startswith(('config', 'progress', 'summary'))]
    

    

    if len(experiment_files) > sample_size:
        sampled_files = random.sample(experiment_files, sample_size)
    else:
        sampled_files = experiment_files
    
    results = []
    for file_path in sampled_files:
        result = check_experiment_file(file_path)
        if result:
            results.append(result)
    

    if results:

        print("-" * 80)
        
        for r in results:




            print(f"  Accessibility: {r['best_accessibility']:.4f} kcal/mol")



        

        print("\n" + "=" * 80)

        avg_acc = sum(r['best_accessibility'] for r in results) / len(results)

        
        constraint_satisfied = sum(1 for r in results if r['constraint_satisfied'])

        
        translation_correct = sum(1 for r in results if r['translation_correct'])


def main():
    """主函数"""
    print("=" * 80)
    print("🔬 验证约束满足情况")
    print("=" * 80)
    

    directories = [
        ("20250909_105804 (被认为是'好'的)", Path("results_full/20250909_105804_unified_access_experiments")),
        ("20250910_022355 (被认为是'差'的)", Path("results_full/20250910_022355_unified_access_experiments")),
    ]
    
    for desc, dir_path in directories:
        if dir_path.exists():
            print(f"\n{'='*80}")
            print(f"🔍 {desc}")
            analyze_directory(dir_path, sample_size=5)
        else:
            print(f"\n❌ 目录不存在: {dir_path}")
    
    print("\n" + "=" * 80)
    print("💡 结论")
    print("=" * 80)
    print("""



""")

if __name__ == "__main__":
    main()
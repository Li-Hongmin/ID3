#!/usr/bin/env python3
"""


"""

import json
from pathlib import Path
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

def analyze_best_validity(file_path: Path):
    """分析best_accessibility记录时是否满足约束"""
    with open(file_path, 'r') as f:
        data = json.load(f)
    

    result = {
        'file': file_path.name,
        'constraint_type': data.get('constraint_type'),
        'protein': data.get('protein_name'),
        'variant': data.get('variant'),
        'best_accessibility': data.get('best_accessibility'),
        'final_accessibility': data.get('final_accessibility'),
        'amino_acids_match': data.get('amino_acids_match'),
        'expected_aa': data.get('expected_amino_acids', ''),
        'actual_aa': data.get('actual_amino_acids', ''),
        'final_sequence': data.get('final_sequence', ''),
    }
    

    if result['final_sequence']:
        translated_aa = translate_rna_to_protein(result['final_sequence'])
        result['translated_aa'] = translated_aa
        result['translation_correct'] = (translated_aa == result['expected_aa'])
    

    if 'trajectory' in data and 'accessibility' in data['trajectory']:
        traj = data['trajectory']['accessibility']
        best_val = result['best_accessibility']
        

        best_positions = [i for i, v in enumerate(traj) if abs(v - best_val) < 0.0001]
        
        if best_positions:
            result['best_iteration'] = best_positions[0]
            result['best_position_pct'] = best_positions[0] / len(traj) * 100
            result['total_iterations'] = len(traj)
            

            final_vals = traj[-100:]
            result['final_mean'] = sum(final_vals) / len(final_vals)
            result['final_std'] = (sum((x - result['final_mean'])**2 for x in final_vals) / len(final_vals))**0.5
            result['best_vs_final_ratio'] = result['final_accessibility'] / best_val
            

            result['suspicious'] = (
                result['best_position_pct'] < 30 or
                result['best_vs_final_ratio'] > 1.5 or
                result['final_std'] > 0.1
            )
        else:
            result['best_iteration'] = None
            result['suspicious'] = True
    
    return result

def main():

    print("=" * 80)

    print("=" * 80)
    

    directories = [
        Path("results_full/20250909_105804_unified_access_experiments"),
        Path("results_full/20250910_022355_unified_access_experiments"),
    ]
    
    for dir_path in directories:
        if not dir_path.exists():

            continue
        

        print("=" * 80)
        

        constraint_results = {
            'lagrangian': [],
            'ams': [],
            'cpc': []
        }
        
        json_files = list(dir_path.glob("*.json"))
        experiment_files = [f for f in json_files 
                           if not f.name.startswith(('config', 'progress', 'summary'))]
        

        
        for file_path in experiment_files:
            try:
                result = analyze_best_validity(file_path)
                constraint_type = result['constraint_type']
                if constraint_type in constraint_results:
                    constraint_results[constraint_type].append(result)
            except Exception as e:

        


        print("-" * 80)
        
        for constraint_type in ['lagrangian', 'ams', 'cpc']:
            results = constraint_results[constraint_type]
            if not results:
                continue
            
            print(f"\n{constraint_type.upper()} (n={len(results)}):")
            

            suspicious_count = sum(1 for r in results if r.get('suspicious', False))
            constraint_satisfied = sum(1 for r in results if r.get('amino_acids_match', False))
            translation_correct = sum(1 for r in results if r.get('translation_correct', False))
            



            

            best_vals = [r['best_accessibility'] for r in results if r['best_accessibility']]
            final_vals = [r['final_accessibility'] for r in results if r['final_accessibility']]
            
            if best_vals and final_vals:



            

            suspicious_examples = [r for r in results if r.get('suspicious', False)]
            if suspicious_examples:

                for i, example in enumerate(suspicious_examples[:3], 1):
                    best_pos = example.get('best_position_pct', 0)
                    ratio = example.get('best_vs_final_ratio', 1)

        

        print("-" * 40)
        

        all_suspicious = []
        all_constraint_types = []
        for constraint_type, results in constraint_results.items():
            if results:
                suspicious_rate = sum(1 for r in results if r.get('suspicious', False)) / len(results)
                all_suspicious.append(suspicious_rate)
                all_constraint_types.append(constraint_type)

        
        if all_suspicious:
            most_suspicious = all_constraint_types[all_suspicious.index(max(all_suspicious))]


if __name__ == "__main__":
    main()
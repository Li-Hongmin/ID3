#!/usr/bin/env python3
"""


"""

import json
import numpy as np
from pathlib import Path
from Bio.Seq import Seq
import sys

# Add project path for imports
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')
from id3.cai.unified_calculator import UnifiedCAICalculator

def validate_sequence_constraint(rna_sequence: str, amino_sequence: str) -> tuple:

    if not rna_sequence or not amino_sequence:

    
    try:

        rna_seq = Seq(rna_sequence.replace('T', 'U'))
        

        translated = str(rna_seq.translate())
        

        if translated.endswith('*'):
            translated = translated[:-1]
        
        if translated == amino_sequence:

        else:

    except Exception as e:


def validate_cai_calculation(rna_sequence: str, expected_cai: float, 
                            species: str = 'ecoli_bl21de3', tolerance: float = 0.001) -> tuple:
    """验证CAI计算是否正确"""
    try:

        calculator = UnifiedCAICalculator(species=species)
        

        calculated_cai = calculator.compute_cai(rna_sequence, method='standard')
        

        diff = abs(calculated_cai - expected_cai)
        if diff < tolerance:
            return True, f"CAI正确: 计算值={calculated_cai:.6f}, 期望值={expected_cai:.6f}, 差异={diff:.6f}"
        else:
            return False, f"CAI不匹配: 计算值={calculated_cai:.6f}, 期望值={expected_cai:.6f}, 差异={diff:.6f}"
    except Exception as e:
        return False, f"CAI计算错误: {str(e)}"

def analyze_l11_experiments():

    exp_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/merged_cai_experiments")
    

    l11_files = list(exp_dir.glob("*lagrangian_11*.json"))
    
    print("="*80)

    print("="*80)

    
    results = []
    constraint_satisfied = 0
    cai_valid = 0
    
    for i, json_file in enumerate(l11_files, 1):

        
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            

            protein = data.get('protein_name', 'unknown')
            constraint = data.get('constraint_type', '')
            variant = data.get('variant', '')
            

            if constraint != 'lagrangian' or variant != '11':

                continue
            

            best_design = data.get('best_seq_design', {})
            discrete_sequence = best_design.get('discrete_sequence', '')
            discrete_cai = best_design.get('discrete_cai', None)
            accessibility = best_design.get('accessibility', None)
            expected_amino = data.get('expected_amino_acids', '')
            



            print(f"  Accessibility: {accessibility:.4f} kcal/mol")

            

            constraint_ok, constraint_msg = validate_sequence_constraint(discrete_sequence, expected_amino)
            if constraint_ok:
                constraint_satisfied += 1
                print(f"  ✅ {constraint_msg}")
            else:

            

            if discrete_cai is not None:
                cai_ok, cai_msg = validate_cai_calculation(discrete_sequence, discrete_cai)
                if cai_ok:
                    cai_valid += 1
                    print(f"  ✅ {cai_msg}")
                else:

            else:

                cai_ok = False
            

            results.append({
                'protein': protein,
                'file': json_file.name,
                'accessibility': accessibility,
                'cai': discrete_cai,
                'constraint_satisfied': constraint_ok,
                'cai_valid': cai_ok,
                'sequence_length': len(discrete_sequence),
                'amino_length': len(expected_amino)
            })
            
        except Exception as e:

        

    

    print("="*80)

    print("="*80)



    
    if results:
        accessibilities = [r['accessibility'] for r in results if r['accessibility'] is not None]
        cais = [r['cai'] for r in results if r['cai'] is not None]
        

        print(f"Accessibility:")




        
        print(f"CAI:")




        

        best_accessibility_idx = np.argmin(accessibilities)
        best_cai_idx = np.argmax(cais)
        



        print(f"  CAI: {results[best_accessibility_idx]['cai']:.6f}")
        


        print(f"  Accessibility: {results[best_cai_idx]['accessibility']:.4f} kcal/mol")
        


        sorted_results = sorted(results, key=lambda x: x['accessibility'] if x['accessibility'] is not None else float('inf'))
        for i, result in enumerate(sorted_results[:5], 1):
            print(f"  {i}. {result['protein']}: {result['accessibility']:.4f} kcal/mol, CAI: {result['cai']:.6f}")
    
    return results

def main():
    """主函数"""
    results = analyze_l11_experiments()
    

    output_path = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures/l11_validation_report.json')
    

    json_results = []
    for r in results:
        json_r = r.copy()

        for key, value in json_r.items():
            if isinstance(value, np.float64):
                json_r[key] = float(value)
            elif isinstance(value, np.int64):
                json_r[key] = int(value)
    
    with open(output_path, 'w') as f:
        json.dump({
            'summary': {
                'total_experiments': len(results),
                'constraint_satisfied_count': sum(1 for r in results if r['constraint_satisfied']),
                'cai_valid_count': sum(1 for r in results if r['cai_valid']),
                'constraint_satisfaction_rate': sum(1 for r in results if r['constraint_satisfied']) / len(results) * 100 if results else 0,
                'cai_validation_rate': sum(1 for r in results if r['cai_valid']) / len(results) * 100 if results else 0
            },
            'experiments': json_results
        }, f, indent=2)
    
    print(f"\n✅ 详细验证报告已保存到: {output_path}")
    print("\n✨ L11验证完成!")

if __name__ == "__main__":
    main()
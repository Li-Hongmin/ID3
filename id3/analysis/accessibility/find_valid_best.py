#!/usr/bin/env python3
"""

"""

import json
import sys
from pathlib import Path
import numpy as np

def find_valid_best_results(directory, top_n=10):

    experiment_files = list(Path(directory).glob("*.json"))
    experiment_files = [f for f in experiment_files 
                       if f.name not in ["config.json", "summary.json"]]
    
    valid_results = []
    invalid_results = []
    
    print("⚡ Scanning experiments for valid results...")
    
    for file_path in experiment_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                

                if 'best_accessibility' not in data:
                    continue
                

                aa_match = data.get('amino_acids_match', False)
                aa_correct = data.get('amino_acids_correct', 0)
                
                result = {
                    'file': file_path.name,
                    'protein': data.get('protein_name'),
                    'constraint': data.get('constraint_type'),
                    'variant': data.get('variant'),
                    'best_acc': data.get('best_accessibility'),
                    'final_acc': data.get('final_accessibility'),
                    'aa_match': aa_match,
                    'aa_correct': aa_correct
                }
                


                if aa_match or aa_correct == 100:
                    valid_results.append(result)
                else:
                    invalid_results.append(result)
                    
        except Exception as e:
            continue
    

    valid_results.sort(key=lambda x: x['best_acc'])
    invalid_results.sort(key=lambda x: x['best_acc'])
    
    print("\n" + "=" * 80)
    print("🏆 TOP VALID RESULTS (Final sequence satisfies constraints)")
    print("=" * 80)
    
    print(f"\n{'Rank':<5} {'Protein':<8} {'Constraint':<12} {'Variant':<8} {'Best Acc':<10} {'Final Acc':<10} {'AA Match'}")
    print("-" * 80)
    
    for i, r in enumerate(valid_results[:top_n], 1):
        print(f"{i:<5} {r['protein']:<8} {r['constraint']:<12} {r['variant']:<8} "
              f"{r['best_acc']:<10.4f} {r['final_acc']:<10.4f} "
              f"{'✅' if r['aa_match'] else '❌'}")
    

    by_constraint = {}
    for r in valid_results:
        constraint = r['constraint']
        if constraint not in by_constraint:
            by_constraint[constraint] = []
        by_constraint[constraint].append(r['best_acc'])
    
    print(f"\n📊 Valid Results by Constraint Type:")
    for constraint in sorted(by_constraint.keys()):
        values = by_constraint[constraint]
        print(f"  {constraint:12s}: Best={min(values):.4f}, Mean={np.mean(values):.4f}, Count={len(values)}")
    

    if invalid_results:
        print(f"\n⚠️  Invalid Results (for comparison):")
        print(f"{'Protein':<8} {'Constraint':<12} {'Variant':<8} {'Best Acc':<10} {'AA Correct %'}")
        print("-" * 60)
        for r in invalid_results[:5]:
            print(f"{r['protein']:<8} {r['constraint']:<12} {r['variant']:<8} "
                  f"{r['best_acc']:<10.4f} {r['aa_correct']:.1f}%")
    

    if valid_results:
        best_valid = valid_results[0]
        print(f"\n" + "=" * 80)
        print("✅ TRUE BEST VALID RESULT")
        print("=" * 80)
        print(f"File: {best_valid['file']}")
        print(f"Protein: {best_valid['protein']}")
        print(f"Constraint: {best_valid['constraint']}")
        print(f"Variant: {best_valid['variant']}")
        print(f"Best Accessibility: {best_valid['best_acc']:.4f}")
        print(f"Final Accessibility: {best_valid['final_acc']:.4f}")
        

        print(f"\n📈 Comparison:")
        print(f"  Invalid '0.2492' result: Not satisfying AA constraint (85.8% match)")
        print(f"  Valid best result: {best_valid['best_acc']:.4f} (100% AA match in final)")
        print(f"  Difference: {best_valid['best_acc'] - 0.2492:.4f}")
    
    return valid_results, invalid_results

if __name__ == "__main__":
    if len(sys.argv) < 2:
        directory = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_105804_unified_access_experiments"
    else:
        directory = sys.argv[1]
    
    valid, invalid = find_valid_best_results(directory)
    
    print(f"\n📊 Summary:")
    print(f"  Total Valid Results: {len(valid)}")
    print(f"  Total Invalid Results: {len(invalid)}")
    print(f"  Valid Rate: {len(valid)/(len(valid)+len(invalid))*100:.1f}%")
#!/usr/bin/env python3
"""

"""

import json
import sys
from pathlib import Path
import numpy as np

def find_best_result(directory):

    experiment_files = list(Path(directory).glob("*.json"))
    experiment_files = [f for f in experiment_files 
                       if f.name not in ["config.json", "summary.json"]]
    
    best_acc = float('inf')
    best_file = None
    best_data = None
    
    for file_path in experiment_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                if 'best_accessibility' in data:
                    if data['best_accessibility'] < best_acc:
                        best_acc = data['best_accessibility']
                        best_file = file_path
                        best_data = data
        except:
            continue
    
    if best_data:
        print("=" * 80)
        print("🔍 BEST RESULT VERIFICATION")
        print("=" * 80)
        
        print(f"\n📁 File: {best_file.name}")
        print(f"🧬 Protein: {best_data.get('protein_name')}")
        print(f"🔧 Constraint: {best_data.get('constraint_type')}")
        print(f"📊 Variant: {best_data.get('variant')}")
        print(f"🎯 Best Accessibility: {best_data.get('best_accessibility'):.6f}")
        print(f"📈 Final Accessibility: {best_data.get('final_accessibility'):.6f}")
        print(f"✅ AA Match Rate: {best_data.get('aa_match_rate', 0)*100:.1f}%")
        print(f"🔄 Iterations: {best_data.get('iterations')}")
        

        if 'best_design' in best_data:
            best_design = best_data['best_design']
            print(f"\n📝 Best Design Details:")
            print(f"   Iteration: {best_design.get('iteration')}")
            print(f"   Accessibility: {best_design.get('accessibility'):.6f}")
            

            if 'sequence' in best_design:
                seq = best_design['sequence']
                print(f"   Sequence Length: {len(seq)}")
                print(f"   First 60 nt: {seq[:60]}")
                

                valid_nucleotides = set('ATCG')
                is_valid = all(n in valid_nucleotides for n in seq)
                print(f"   Valid Nucleotides: {'✅ Yes' if is_valid else '❌ No'}")
                

                if 'amino_acid_match' in best_design:
                    aa_match = best_design['amino_acid_match']
                    print(f"   AA Constraint Satisfied: {'✅ Yes' if aa_match else '❌ No'}")
        

        if 'trajectory' in best_data and 'accessibility' in best_data['trajectory']:
            traj = best_data['trajectory']['accessibility']
            print(f"\n📈 Optimization Trajectory:")
            print(f"   Initial: {traj[0]:.4f}")
            print(f"   Best: {min(traj):.4f}")
            print(f"   Final: {traj[-1]:.4f}")
            print(f"   Best at iteration: {traj.index(min(traj))}")
            

            if abs(min(traj) - best_acc) > 0.0001:
                print(f"   ⚠️  WARNING: Trajectory min ({min(traj):.6f}) != reported best ({best_acc:.6f})")
        

        print(f"\n🔒 Constraint Verification:")
        constraint_type = best_data.get('constraint_type')
        variant = best_data.get('variant', '00')
        
        print(f"   Type: {constraint_type}")
        print(f"   Variant: {variant} (α={variant[0]}, β={variant[1]})")
        

        if 'final_constraint_satisfaction' in best_data:
            satisfaction = best_data['final_constraint_satisfaction']
            print(f"   Final Constraint Satisfaction: {satisfaction}")
        

        output_file = Path("best_result_verification.json")
        with open(output_file, 'w') as f:
            json.dump({
                'file': str(best_file),
                'best_accessibility': best_acc,
                'protein': best_data.get('protein_name'),
                'constraint': best_data.get('constraint_type'),
                'variant': best_data.get('variant'),
                'aa_match_rate': best_data.get('aa_match_rate'),
                'valid': True if best_data.get('aa_match_rate', 0) == 1.0 else False
            }, f, indent=2)
        
        print(f"\n💾 Verification saved to: {output_file}")
        
        return best_data
    else:
        print("❌ No valid results found")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python verify_best_result.py <experiment_dir>")
        sys.exit(1)
    
    directory = sys.argv[1]
    find_best_result(directory)
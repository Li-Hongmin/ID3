#!/usr/bin/env python3
"""

"""

import json
import sys
from pathlib import Path
import torch
import numpy as np

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.utils.constants import amino_acids_to_codons

def translate_sequence(dna_sequence):

    codon_to_aa = {}
    for aa, codons in amino_acids_to_codons.items():
        for codon in codons:
            codon_to_aa[codon.replace('U', 'T')] = aa
    

    aa_sequence = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if codon in codon_to_aa:
            aa_sequence.append(codon_to_aa[codon])
        else:
            aa_sequence.append('X')  # Unknown
    
    return ''.join(aa_sequence)

def verify_experiment(file_path):
    """深度验证实验结果"""
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    print("=" * 80)
    print("🔬 DEEP VERIFICATION REPORT")
    print("=" * 80)
    

    print(f"\n📁 File: {Path(file_path).name}")
    print(f"🧬 Protein: {data.get('protein_name')}")
    print(f"🔧 Constraint: {data.get('constraint_type')}")
    print(f"📊 Variant: {data.get('variant')}")
    

    print(f"\n📈 Performance Metrics:")
    print(f"   Initial Accessibility: {data.get('initial_accessibility'):.4f}")
    print(f"   Best Accessibility: {data.get('best_accessibility'):.4f}")
    print(f"   Final Accessibility: {data.get('final_accessibility'):.4f}")
    

    print(f"\n🔒 Constraint Satisfaction:")
    print(f"   AA Match (reported): {data.get('amino_acids_match')}")
    print(f"   AA Correct % (reported): {data.get('amino_acids_correct', 0):.1f}%")
    

    final_seq = data.get('final_sequence', '')
    if final_seq:

        final_dna = final_seq.replace('U', 'T')
        

        actual_aa = translate_sequence(final_dna)
        expected_aa = data.get('expected_amino_acids', '')
        
        print(f"\n🧬 Final Sequence Verification:")
        print(f"   Sequence Length: {len(final_dna)}")
        print(f"   Expected AA Length: {len(expected_aa)}")
        print(f"   Actual AA Length: {len(actual_aa)}")
        

        if expected_aa:
            matches = sum(1 for a, e in zip(actual_aa, expected_aa) if a == e)
            match_rate = matches / len(expected_aa) if expected_aa else 0
            print(f"   Calculated AA Match: {match_rate*100:.1f}%")
            
            if match_rate < 1.0:
                print(f"   ⚠️  WARNING: AA sequence doesn't match!")

                for i, (a, e) in enumerate(zip(actual_aa[:10], expected_aa[:10])):
                    if a != e:
                        print(f"      Position {i}: {e} → {a}")
    

    if 'trajectory' in data and 'accessibility' in data['trajectory']:
        traj = data['trajectory']['accessibility']
        best_idx = np.argmin(traj)
        best_value = traj[best_idx]
        
        print(f"\n📊 Trajectory Analysis:")
        print(f"   Best Value in Trajectory: {best_value:.4f}")
        print(f"   Best Iteration: {best_idx}")
        print(f"   Reported Best: {data.get('best_accessibility'):.4f}")
        

        if abs(best_value - data.get('best_accessibility', 0)) > 0.001:
            print(f"   ❌ ERROR: Trajectory best != reported best!")
        else:
            print(f"   ✅ Trajectory consistent with reported best")
        

        if 'discrete_sequences' in data['trajectory']:
            sequences = data['trajectory']['discrete_sequences']
            if best_idx < len(sequences):
                best_seq = sequences[best_idx]
                best_aa = translate_sequence(best_seq.replace('U', 'T'))
                
                print(f"\n🏆 Best Point Sequence Analysis:")
                print(f"   Best Iteration: {best_idx}")
                print(f"   Best Sequence Length: {len(best_seq)}")
                

                if expected_aa:
                    best_matches = sum(1 for a, e in zip(best_aa, expected_aa) if a == e)
                    best_match_rate = best_matches / len(expected_aa) if expected_aa else 0
                    print(f"   Best Sequence AA Match: {best_match_rate*100:.1f}%")
                    
                    if best_match_rate < 1.0:
                        print(f"   ⚠️  WARNING: Best sequence doesn't satisfy AA constraint!")
                        print(f"   This explains the low accessibility - it's not a valid solution!")
    

    variant = data.get('variant', '00')
    print(f"\n⚙️  Variant Configuration:")
    print(f"   Variant: {variant}")
    print(f"   α (exploration): {variant[0]}")
    print(f"   β (discretization): {variant[1]}")
    
    if variant[1] == '0':
        print(f"   ⚠️  Note: β=0 means continuous optimization (no discrete sequences)")
    

    print(f"\n" + "=" * 80)
    print("📋 SUMMARY")
    print("=" * 80)
    
    issues = []
    

    if 'trajectory' in data and 'discrete_sequences' in data['trajectory']:
        sequences = data['trajectory']['discrete_sequences']
        traj = data['trajectory']['accessibility']
        best_idx = np.argmin(traj)
        
        if best_idx < len(sequences):
            best_seq = sequences[best_idx]
            best_aa = translate_sequence(best_seq.replace('U', 'T'))
            
            if expected_aa:
                best_match_rate = sum(1 for a, e in zip(best_aa, expected_aa) if a == e) / len(expected_aa)
                if best_match_rate < 1.0:
                    issues.append(f"Best sequence (at iteration {best_idx}) doesn't satisfy AA constraint ({best_match_rate*100:.1f}% match)")
    
    if issues:
        print("❌ CRITICAL ISSUES FOUND:")
        for issue in issues:
            print(f"   • {issue}")
        print("\nThis result is INVALID because it doesn't satisfy the constraints!")
    else:
        print("✅ Result appears valid")
    
    return data

if __name__ == "__main__":
    if len(sys.argv) < 2:

        file_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_105804_unified_access_experiments/20250909_112021_P04637_lagrangian_01_seed42.json"
    else:
        file_path = sys.argv[1]
    
    verify_experiment(file_path)
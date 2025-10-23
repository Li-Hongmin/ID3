#!/usr/bin/env python3
"""
Simple verification of Lagrangian constraint fix.
Runs small experiments using the command-line interface.
"""

import subprocess
import json
from pathlib import Path
from datetime import datetime
import time

def run_experiment(protein, constraint, variant, iterations=200):
    """Run a single experiment using subprocess."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = f"results/{timestamp}_verify_{constraint}_{variant}"
    
    cmd = [
        "python", "run_unified_experiment.py",
        "--proteins", protein,
        "--constraints", constraint,
        "--variants", variant,
        "--iterations", str(iterations),
        "--enable-cai",
        "--cai-target", "0.8",
        "--lambda-cai", "0.1",
        "--output-dir", output_dir,
        "--seeds", "1",
        "--base-seed", "42",
    ]
    
    print(f"Running: {constraint}_{variant}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running {constraint}_{variant}:")
        print(result.stderr)
        return None
    
    # Parse results
    output_path = Path(output_dir)
    experiment_files = list(output_path.glob("*.json"))
    experiment_files = [f for f in experiment_files if f.name not in ["config.json", "summary.json"]]
    
    if experiment_files:
        with open(experiment_files[0], 'r') as f:
            exp_data = json.load(f)
            
        # Calculate uniqueness
        if 'trajectory' in exp_data and 'discrete_sequences' in exp_data['trajectory']:
            sequences = exp_data['trajectory']['discrete_sequences']
            unique_sequences = set(sequences)
            uniqueness_rate = len(unique_sequences) / len(sequences) if sequences else 0
            
            return {
                "uniqueness_rate": uniqueness_rate,
                "total_sequences": len(sequences),
                "unique_sequences": len(unique_sequences),
                "best_accessibility": exp_data.get("best_accessibility"),
            }
    
    return None


def main():
    print("=" * 80)
    print("VERIFYING LAGRANGIAN CONSTRAINT FIX")
    print("=" * 80)
    print("\nRunning small-scale experiments to verify uniqueness improvement...")
    print("Expected: Lagrangian should now have ~40-50% uniqueness (was 6.2%)\n")
    
    protein = "O15263"
    variant = "11"  # alpha=1, beta=1 for maximum stochasticity
    
    results = {}
    
    # Test all three constraint types
    for constraint in ["lagrangian", "ams", "cpc"]:
        print(f"\n{'='*60}")
        print(f"Testing {constraint.upper()} constraint (variant {variant})")
        print(f"{'='*60}")
        
        result = run_experiment(protein, constraint, variant, iterations=200)
        
        if result:
            results[constraint] = result
            print(f"✓ Uniqueness Rate: {result['uniqueness_rate']:.1%}")
            print(f"  Unique/Total: {result['unique_sequences']}/{result['total_sequences']}")
            print(f"  Best Accessibility: {result['best_accessibility']:.4f}")
        else:
            print(f"✗ Failed to run {constraint}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    if "lagrangian" in results:
        lag_rate = results["lagrangian"]["uniqueness_rate"]
        
        print(f"\nLagrangian uniqueness: {lag_rate:.1%}")
        
        if "ams" in results:
            print(f"AMS uniqueness: {results['ams']['uniqueness_rate']:.1%}")
        if "cpc" in results:
            print(f"CPC uniqueness: {results['cpc']['uniqueness_rate']:.1%}")
        
        if lag_rate > 0.25:  # Expecting significant improvement
            print(f"\n✅ SUCCESS! Lagrangian fix is working!")
            print(f"   New rate: {lag_rate:.1%} (was ~6.2% before fix)")
            print(f"   The constraint now properly adds Gumbel noise.")
        elif lag_rate > 0.15:
            print(f"\n⚠️  PARTIAL SUCCESS: Some improvement detected")
            print(f"   New rate: {lag_rate:.1%} (was ~6.2% before fix)")
            print(f"   Consider increasing alpha or iterations for better exploration.")
        else:
            print(f"\n❌ PROBLEM: Lagrangian still has low diversity")
            print(f"   Current rate: {lag_rate:.1%}")
            print(f"   Please verify the fix is properly applied.")
    else:
        print("\n❌ Could not run Lagrangian experiment")
    
    print("\n" + "=" * 80)
    
    # Save verification report
    report_path = Path(f"results/lagrangian_fix_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(report_path, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "results": results,
            "conclusion": "fixed" if results.get("lagrangian", {}).get("uniqueness_rate", 0) > 0.25 else "not_fixed"
        }, f, indent=2)
    
    print(f"Report saved to: {report_path}")


if __name__ == "__main__":
    main()
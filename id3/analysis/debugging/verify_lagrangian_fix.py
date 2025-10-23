#!/usr/bin/env python3
"""
Verification experiment for Lagrangian constraint Gumbel noise fix.

This script runs a small-scale experiment to verify that the fixed Lagrangian
constraint now achieves similar uniqueness rates to AMS/CPC constraints.
"""

import json
import logging
import sys
from pathlib import Path
from datetime import datetime
import torch
import numpy as np

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def calculate_uniqueness_rate(sequences):
    """Calculate the uniqueness rate of sequences."""
    if not sequences:
        return 0.0
    unique_sequences = set(sequences)
    return len(unique_sequences) / len(sequences)


def run_verification_experiment():
    """Run a small-scale experiment to verify the Lagrangian fix."""
    
    print("=" * 80)
    print("LAGRANGIAN CONSTRAINT FIX VERIFICATION")
    print("=" * 80)
    
    # Create output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"results/{timestamp}_lagrangian_fix_verification")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Test configurations
    test_configs = [
        {
            "name": "lagrangian_11",
            "constraint_type": "lagrangian",
            "alpha": 1.0,
            "beta": 1.0,
        },
        {
            "name": "ams_11",
            "constraint_type": "ams",
            "alpha": 1.0,
            "beta": 1.0,
        },
        {
            "name": "cpc_11",
            "constraint_type": "cpc",
            "alpha": 1.0,
            "beta": 1.0,
        },
    ]
    
    results = {}
    
    for test_config in test_configs:
        print(f"\n{'='*60}")
        print(f"Testing: {test_config['name']}")
        print(f"{'='*60}")
        
        # Create experiment configuration
        config = UnifiedExperimentConfig(
            protein_name="O15263",  # Use a real protein
            constraint_type=test_config["constraint_type"],
            alpha=test_config["alpha"],
            beta=test_config["beta"],
            num_iterations=200,  # Small scale for quick verification
            cai_mode=True,
            target_cai=0.8,
            lambda_cai=0.1,
            species="ecoli_bl21de3",
            random_seed=42,
            output_dir=str(output_dir / test_config["name"]),
            save_trajectory=True,
            save_final_state=True,
        )
        
        # Run experiment
        runner = UnifiedExperimentRunner(config)
        result = runner.run()
        
        # Extract discrete sequences
        if 'trajectory' in result and 'discrete_sequences' in result['trajectory']:
            sequences = result['trajectory']['discrete_sequences']
            uniqueness_rate = calculate_uniqueness_rate(sequences)
            
            results[test_config["name"]] = {
                "uniqueness_rate": uniqueness_rate,
                "total_sequences": len(sequences),
                "unique_sequences": len(set(sequences)),
                "best_accessibility": result.get("best_accessibility"),
                "final_accessibility": result.get("final_accessibility"),
            }
            
            print(f"\nResults for {test_config['name']}:")
            print(f"  Uniqueness Rate: {uniqueness_rate:.1%}")
            print(f"  Unique Sequences: {len(set(sequences))}/{len(sequences)}")
            print(f"  Best Accessibility: {result.get('best_accessibility'):.4f}")
        else:
            print(f"  ERROR: No discrete sequences found")
            results[test_config["name"]] = {"error": "No discrete sequences found"}
    
    # Print comparison
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    
    print("\nUniqueness Rates:")
    for name, data in results.items():
        if "uniqueness_rate" in data:
            print(f"  {name:20s}: {data['uniqueness_rate']:.1%}")
        else:
            print(f"  {name:20s}: ERROR")
    
    # Check if Lagrangian is fixed
    if "lagrangian_11" in results and "uniqueness_rate" in results["lagrangian_11"]:
        lag_rate = results["lagrangian_11"]["uniqueness_rate"]
        
        if lag_rate > 0.20:  # Expecting >20% after fix (was 6.2% before)
            print("\n✅ SUCCESS: Lagrangian constraint now shows good diversity!")
            print(f"   New uniqueness rate: {lag_rate:.1%} (was ~6.2% before fix)")
        else:
            print("\n❌ ISSUE: Lagrangian constraint still shows low diversity")
            print(f"   Uniqueness rate: {lag_rate:.1%}")
            print("   Please verify the fix is properly applied")
    
    # Save results
    results_file = output_dir / "verification_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to: {results_file}")
    print("=" * 80)
    
    return results


if __name__ == "__main__":
    try:
        results = run_verification_experiment()
        
        # Exit with appropriate code
        if "lagrangian_11" in results and "uniqueness_rate" in results["lagrangian_11"]:
            if results["lagrangian_11"]["uniqueness_rate"] > 0.20:
                sys.exit(0)  # Success
            else:
                sys.exit(1)  # Fix not working as expected
        else:
            sys.exit(2)  # Error in experiment
            
    except Exception as e:
        logger.error(f"Verification failed: {e}", exc_info=True)
        sys.exit(3)
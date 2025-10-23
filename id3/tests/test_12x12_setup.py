#!/usr/bin/env python3
"""
Test script for 12x12 experiment setup
Tests the protein loading, CAI integration, and basic functionality
"""

import sys
from pathlib import Path

# Add project path
sys.path.append(str(Path(__file__).parent))

from run_12x12_v3_complete import ProteinSequenceLoader, Complete12x12ExperimentRunner
from id3.cai.module import CAIModule

def test_protein_loading():
    """Test protein sequence loading from data/ folder"""
    print("🧪 Testing protein sequence loading...")
    
    loader = ProteinSequenceLoader()
    sequences = loader.get_sequences()
    
    print(f"📊 Loaded {len(sequences)} protein sequences:")
    for protein_id, data in sequences.items():
        print(f"   {protein_id}: {data['length']} amino acids")
        print(f"      First 50 chars: {data['sequence'][:50]}...")
        print(f"      Source: {Path(data['source_file']).name}")
    
    return len(sequences) >= 8  # Need at least 8 for a meaningful test

def test_cai_integration():
    """Test CAI integration with E.coli BL21(DE3) reference"""
    print("\n🧪 Testing CAI integration...")
    
    try:
        cai_module = CAIModule(reference_species='ecoli_bl21de3')
        
        # Test with a sample RNA sequence
        test_rna = "ATGAAACGCGAAGAAGGCACCGAACGT"  # Simple test sequence
        cai_score = cai_module.compute_cai_score(test_rna)
        cai_stats = cai_module.get_cai_statistics(test_rna)
        
        print(f"   ✅ CAI module initialized successfully")
        print(f"   Test CAI score: {cai_score:.4f}")
        print(f"   Rare codon ratio: {cai_stats['rare_codon_ratio']:.2%}")
        print(f"   Total codons: {cai_stats['total_codons']}")
        
        return True
    except Exception as e:
        print(f"   ❌ CAI integration failed: {e}")
        return False

def test_constraint_imports():
    """Test that v3 constraint implementations can be imported"""
    print("\n🧪 Testing v3 constraint imports...")
    
    try:
        # Test imports with correct class names
        from id3.constraints.cpc_v3_stable import CodonProfileConstraintV3Stable
        from id3.constraints.ams_v3_stable import AdaptiveMultilevelSamplingV3Stable
        from id3.constraints.lagrangian_v3_stable import LagrangianConstraintV3Stable
        
        print("   ✅ CPC v3 constraint imported successfully")
        print("   ✅ AMS v3 constraint imported successfully")
        print("   ✅ Lagrangian v3 constraint imported successfully")
        
        return True
    except ImportError as e:
        print(f"   ❌ Constraint import failed: {e}")
        return False

def test_experiment_runner_init():
    """Test experiment runner initialization"""
    print("\n🧪 Testing experiment runner initialization...")
    
    try:
        runner = Complete12x12ExperimentRunner(output_dir="test_output")
        
        available_proteins = runner.protein_loader.get_sequence_list()
        print(f"   ✅ Runner initialized successfully")
        print(f"   Available proteins: {len(available_proteins)}")
        print(f"   Constraint types: {len(runner.constraint_types)}")
        print(f"   Variants: {len(runner.variants)}")
        
        return True
    except Exception as e:
        print(f"   ❌ Runner initialization failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_codon_reference_data():
    """Test E.coli BL21(DE3) codon reference data loading"""
    print("\n🧪 Testing codon reference data...")
    
    try:
        import json
        ref_file = "data/codon_references/ecoli_bl21de3_complete_reference_sequences.json"
        
        if Path(ref_file).exists():
            with open(ref_file, 'r') as f:
                data = json.load(f)
            
            sequences = data.get('sequences', {})
            metadata = data.get('metadata', {})
            
            print(f"   ✅ Reference file loaded successfully")
            print(f"   Total sequences: {metadata.get('total_sequences', 'unknown')}")
            print(f"   Sequence count: {len(sequences)}")
            print(f"   Sample genes: {list(sequences.keys())[:5]}")
            
            return True
        else:
            print(f"   ❌ Reference file not found: {ref_file}")
            return False
            
    except Exception as e:
        print(f"   ❌ Codon reference test failed: {e}")
        return False

def run_all_tests():
    """Run all setup tests"""
    print("🚀 Running 12x12 Experiment Setup Tests")
    print("=" * 50)
    
    tests = [
        ("Protein Loading", test_protein_loading),
        ("CAI Integration", test_cai_integration),
        ("Constraint Imports", test_constraint_imports),
        ("Codon Reference Data", test_codon_reference_data),
        ("Experiment Runner", test_experiment_runner_init),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"   ❌ Test {test_name} crashed: {e}")
            results.append((test_name, False))
    
    print("\n📊 Test Summary:")
    print("-" * 30)
    
    passed = 0
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name:<25} {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{len(tests)} tests passed")
    
    if passed == len(tests):
        print("🎉 All tests passed! Ready for 12x12 experiment.")
        return True
    else:
        print("⚠️  Some tests failed. Please fix issues before running full experiment.")
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
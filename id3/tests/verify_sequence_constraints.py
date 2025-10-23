#!/usr/bin/env python3
"""




"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from id3.utils.sequence_utils import rna_to_amino_acids

class SequenceConstraintVerifier:

    
    def __init__(self):
        self.violation_patterns = {
            'total_violations': 0,


            'method_violations': {'lagrangian': 0, 'ams': 0, 'cpc': 0}
        }
    
    def verify_single_file(self, file_path: Path) -> Dict:
        """验证单个实验文件的序列约束"""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            

            discrete_sequence = data.get('final_sequence', '')
            expected_amino_acids = data.get('expected_amino_acids', '')
            actual_amino_acids = data.get('actual_amino_acids', '')
            

            amino_acids_match = data.get('amino_acids_match', False)
            amino_acids_correct = data.get('amino_acids_correct', 0.0)
            

            if discrete_sequence:
                computed_amino_acids = rna_to_amino_acids(discrete_sequence)
                recomputed_match = (computed_amino_acids == expected_amino_acids)
            else:
                computed_amino_acids = ""
                recomputed_match = False
            

            best_accessibility = data.get('best_accessibility', float('inf'))
            best_seq_design = data.get('best_seq_design', {})
            
            best_sequence_constraint_check = None
            if best_seq_design and 'discrete_sequence' in best_seq_design:
                best_discrete_seq = best_seq_design['discrete_sequence']
                if best_discrete_seq:
                    best_amino_acids = rna_to_amino_acids(best_discrete_seq)
                    best_sequence_constraint_check = (best_amino_acids == expected_amino_acids)
            
            result = {
                'file': file_path.name,
                'method': self._extract_method(file_path.name),
                'variant': self._extract_variant(file_path.name),
                'protein': self._extract_protein(file_path.name),
                
                # Final sequence verification
                'final_sequence_length': len(discrete_sequence),
                'expected_aa_length': len(expected_amino_acids),
                'official_match': amino_acids_match,
                'official_correct': amino_acids_correct,
                'recomputed_match': recomputed_match,
                'consistency_check': (amino_acids_match == recomputed_match),
                
                # Sequences
                'expected_amino_acids': expected_amino_acids,
                'official_actual_amino_acids': actual_amino_acids,
                'recomputed_amino_acids': computed_amino_acids,
                
                # Best accessibility sequence check
                'best_accessibility': best_accessibility,
                'best_sequence_has_constraint_check': best_sequence_constraint_check,
                
                # Error details
                'sequence_errors': []
            }
            

            if not recomputed_match and expected_amino_acids and computed_amino_acids:
                result['sequence_errors'] = self._analyze_sequence_differences(
                    expected_amino_acids, computed_amino_acids
                )
            
            return result
            
        except Exception as e:
            return {
                'file': file_path.name,
                'error': str(e),
                'verification_failed': True
            }
    
    def _extract_method(self, filename: str) -> str:

        if 'lagrangian' in filename.lower():
            return 'lagrangian'
        elif 'ams' in filename.lower():
            return 'ams'
        elif 'cpc' in filename.lower():
            return 'cpc'
        return 'unknown'
    
    def _extract_variant(self, filename: str) -> str:
        """从文件名提取变体"""
        for variant in ['00', '01', '10', '11']:
            if f'_{variant}_' in filename:
                return variant
        return 'unknown'
    
    def _extract_protein(self, filename: str) -> str:

        parts = filename.split('_')
        for part in parts:
            if part.startswith(('P', 'O')) and part[1:].isdigit():
                return part
        return 'unknown'
    
    def _analyze_sequence_differences(self, expected: str, actual: str) -> List[Dict]:
        """分析序列差异的具体位置"""
        errors = []
        min_len = min(len(expected), len(actual))
        
        for i in range(min_len):
            if expected[i] != actual[i]:
                errors.append({
                    'position': i,
                    'expected': expected[i],
                    'actual': actual[i]
                })
        

        if len(expected) != len(actual):
            errors.append({
                'position': 'length',
                'expected_length': len(expected),
                'actual_length': len(actual)
            })
        
        return errors
    
    def analyze_batch(self, batch_path: Path) -> Dict:


        
        json_files = list(batch_path.glob("*.json"))
        if not json_files:

            return {}
        
        results = {
            'batch_name': batch_path.name,
            'total_files': len(json_files),
            'verification_results': [],
            'summary': {
                'methods': {'lagrangian': 0, 'ams': 0, 'cpc': 0},
                'violations': {'lagrangian': 0, 'ams': 0, 'cpc': 0},
                'best_sequence_violations': 0,
                'consistency_errors': 0
            }
        }
        
        for file_path in json_files:
            if file_path.name in ['summary.json', 'progress.json']:
                continue
                
            result = self.verify_single_file(file_path)
            results['verification_results'].append(result)
            

            if 'error' not in result:
                method = result['method']
                if method not in results['summary']['methods']:
                    results['summary']['methods'][method] = 0
                    results['summary']['violations'][method] = 0
                results['summary']['methods'][method] += 1
                

                if not result['recomputed_match']:
                    results['summary']['violations'][method] += 1
                

                if result['best_sequence_has_constraint_check'] == False:
                    results['summary']['best_sequence_violations'] += 1
                

                if not result['consistency_check']:
                    results['summary']['consistency_errors'] += 1
        

        for method in ['lagrangian', 'ams', 'cpc']:
            total = results['summary']['methods'][method]
            violations = results['summary']['violations'][method]
            if total > 0:
                results['summary'][f'{method}_violation_rate'] = violations / total * 100
            else:
                results['summary'][f'{method}_violation_rate'] = 0
        
        return results
    
    def print_analysis_report(self, results: Dict):
        """打印详细分析报告"""
        print(f"\n📊 约束验证分析报告: {results['batch_name']}")
        print("=" * 60)
        
        summary = results['summary']
        
        print(f"总文件数: {results['total_files']}")
        print(f"验证文件数: {len(results['verification_results'])}")
        
        print(f"\n📈 方法分布:")
        for method in ['lagrangian', 'ams', 'cpc']:
            count = summary['methods'][method]
            violations = summary['violations'][method]
            rate = summary[f'{method}_violation_rate']
            
            if count > 0:
                status = "🔴" if rate > 10 else "🟡" if rate > 0 else "🟢"
                print(f"  {method.upper()}: {count}个文件, {violations}个违反 ({rate:.1f}%) {status}")
        
        print(f"\n🚨 关键问题:")
        print(f"  Best序列约束违反: {summary['best_sequence_violations']}")
        print(f"  记录一致性错误: {summary['consistency_errors']}")
        

        violation_examples = [r for r in results['verification_results'] 
                             if not r.get('recomputed_match', True) and 'error' not in r]
        
        if violation_examples:
            print(f"\n❌ 违反约束的案例 (前5个):")
            for example in violation_examples[:5]:
                print(f"  {example['file']}: {example['method'].upper()} "
                      f"expected={len(example['expected_amino_acids'])}aa, "
                      f"got={len(example['recomputed_amino_acids'])}aa, "
                      f"errors={len(example['sequence_errors'])}")
        

        best_violations = [r for r in results['verification_results'] 
                          if r.get('best_sequence_has_constraint_check') == False]
        
        if best_violations:
            print(f"\n⚠️ Best序列违反约束的案例:")
            for example in best_violations[:3]:
                print(f"  {example['file']}: best_accessibility={example['best_accessibility']:.3f}")

def main():

    verifier = SequenceConstraintVerifier()
    

    test_batch = Path("results/20250910_133850_unified_access_experiments")
    
    if test_batch.exists():

        results = verifier.analyze_batch(test_batch)
        verifier.print_analysis_report(results)
        

        output_file = Path("sequence_constraint_verification_report.json")
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

        
    else:

        


        for results_dir in results_dirs:
            if results_dir.exists():

                for batch_dir in results_dir.iterdir():
                    if batch_dir.is_dir():
                        json_count = len(list(batch_dir.glob("*.json")))


if __name__ == "__main__":
    main()
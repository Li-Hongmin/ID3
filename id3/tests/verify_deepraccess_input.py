#!/usr/bin/env python3
"""




"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import torch
import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
from id3.utils.sequence_utils import rna_to_amino_acids
from id3.experiments.utils.data_loader import ProteinDataLoader
from id3.utils.constants import NUCLEOTIDE_MAP, NUCLEOTIDES

class DeepRaccessInputVerifier:

    
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.deepraccess = DeepRaccessID3Wrapper()
        self.data_loader = ProteinDataLoader()
        self.sequence_issues = []
    
    def string_to_one_hot_tensor(self, sequence: str) -> torch.Tensor:
        """将RNA序列转换为one-hot张量（复制实验代码的实现）"""
        seq_len = len(sequence)
        one_hot = torch.zeros(1, seq_len, 4, device=self.device)
        
        for i, nucleotide in enumerate(sequence):
            if nucleotide in NUCLEOTIDE_MAP:
                one_hot[0, i, NUCLEOTIDE_MAP[nucleotide]] = 1.0
        
        return one_hot
    
    def verify_sequence_reconstruction(self, experiment_file: Path) -> Dict:

        try:
            with open(experiment_file, 'r') as f:
                data = json.load(f)
            

            protein_name = self._extract_protein(experiment_file.name)
            final_sequence = data.get('final_sequence', '')
            expected_amino_acids = data.get('expected_amino_acids', '')
            final_accessibility = data.get('final_accessibility', 0.0)
            

            protein_info = self.data_loader.get_protein_info(protein_name)
            

            cds_verification = self._verify_cds_sequence(
                final_sequence, expected_amino_acids
            )
            

            full_sequence = protein_info['utr5'] + final_sequence + protein_info['utr3']
            full_sequence_verification = self._verify_full_sequence(
                protein_info, final_sequence, full_sequence
            )
            

            deepraccess_input_verification = self._verify_deepraccess_input_format(
                full_sequence
            )
            

            accessibility_verification = self._verify_accessibility_calculation(
                full_sequence, protein_info, final_accessibility
            )
            
            result = {
                'file': experiment_file.name,
                'protein': protein_name,
                'method': self._extract_method(experiment_file.name),
                'variant': self._extract_variant(experiment_file.name),
                

                'cds_length': len(final_sequence),
                'full_sequence_length': len(full_sequence),
                'utr5_length': len(protein_info['utr5']),
                'utr3_length': len(protein_info['utr3']),
                

                'cds_verification': cds_verification,
                'full_sequence_verification': full_sequence_verification,
                'deepraccess_input_verification': deepraccess_input_verification,
                'accessibility_verification': accessibility_verification,
                

                'original_final_accessibility': final_accessibility,
                'sequences': {
                    'cds': final_sequence[:100] + ('...' if len(final_sequence) > 100 else ''),
                    'full': full_sequence[:200] + ('...' if len(full_sequence) > 200 else ''),
                }
            }
            
            return result
            
        except Exception as e:
            return {
                'file': experiment_file.name,
                'error': str(e),
                'verification_failed': True
            }
    
    def _verify_cds_sequence(self, cds_sequence: str, expected_amino_acids: str) -> Dict:
        """验证CDS序列的正确性"""
        verification = {
            'sequence_valid': True,
            'length_correct': True,
            'nucleotides_valid': True,
            'amino_acid_match': True,
            'issues': []
        }
        

        if len(cds_sequence) % 3 != 0:
            verification['length_correct'] = False
            verification['issues'].append(f"CDS长度不是3的倍数: {len(cds_sequence)}")
        

        valid_nucleotides = set(NUCLEOTIDES)
        invalid_nucleotides = set(cds_sequence) - valid_nucleotides
        if invalid_nucleotides:
            verification['nucleotides_valid'] = False
            verification['issues'].append(f"非法核苷酸: {invalid_nucleotides}")
        

        if cds_sequence:
            translated_amino_acids = rna_to_amino_acids(cds_sequence)
            if translated_amino_acids != expected_amino_acids:
                verification['amino_acid_match'] = False
                verification['issues'].append(
                    f"氨基酸不匹配: expected={len(expected_amino_acids)}, got={len(translated_amino_acids)}"
                )
        
        verification['sequence_valid'] = (
            verification['length_correct'] and 
            verification['nucleotides_valid'] and 
            verification['amino_acid_match']
        )
        
        return verification
    
    def _verify_full_sequence(self, protein_info: Dict, cds_sequence: str, full_sequence: str) -> Dict:

        verification = {
            'splicing_correct': True,
            'length_consistent': True,
            'issues': []
        }
        
        utr5 = protein_info['utr5']
        utr3 = protein_info['utr3']
        expected_full = utr5 + cds_sequence + utr3
        

        if full_sequence != expected_full:
            verification['splicing_correct'] = False

        

        expected_length = len(utr5) + len(cds_sequence) + len(utr3)
        if len(full_sequence) != expected_length:
            verification['length_consistent'] = False
            verification['issues'].append(

            )
        

        atg_position = len(utr5)
        if atg_position + 3 <= len(full_sequence):
            start_codon = full_sequence[atg_position:atg_position+3]
            if start_codon != 'AUG':

        
        return verification
    
    def _verify_deepraccess_input_format(self, full_sequence: str) -> Dict:
        """验证DeepRaccess输入格式的正确性"""
        verification = {
            'tensor_creation_success': True,
            'tensor_shape_correct': True,
            'tensor_device_correct': True,
            'encoding_correct': True,
            'issues': []
        }
        
        try:

            tensor = self.string_to_one_hot_tensor(full_sequence)
            

            expected_shape = (1, len(full_sequence), 4)
            if tensor.shape != expected_shape:
                verification['tensor_shape_correct'] = False
                verification['issues'].append(
                    f"张量形状错误: expected={expected_shape}, got={tensor.shape}"
                )
            

            if tensor.device != self.device:
                verification['tensor_device_correct'] = False
                verification['issues'].append(
                    f"设备不匹配: expected={self.device}, got={tensor.device}"
                )
            

            for i in range(min(10, len(full_sequence))):
                nucleotide = full_sequence[i]
                if nucleotide in NUCLEOTIDE_MAP:
                    expected_index = NUCLEOTIDE_MAP[nucleotide]
                    actual_hot = tensor[0, i, :]
                    max_index = torch.argmax(actual_hot).item()
                    
                    if max_index != expected_index or actual_hot[expected_index].item() != 1.0:
                        verification['encoding_correct'] = False
                        verification['issues'].append(
                            f"位置{i}编码错误: {nucleotide} -> expected={expected_index}, got={max_index}"
                        )
                        break
            
        except Exception as e:
            verification['tensor_creation_success'] = False
            verification['issues'].append(f"张量创建失败: {str(e)}")
        
        return verification
    
    def _verify_accessibility_calculation(self, full_sequence: str, protein_info: Dict, 
                                        original_accessibility: float) -> Dict:

        verification = {
            'calculation_success': True,
            'value_reasonable': True,
            'consistent_with_original': True,
            'issues': []
        }
        
        try:

            full_tensor = self.string_to_one_hot_tensor(full_sequence)
            atg_position = len(protein_info['utr5'])
            
            recalculated_accessibility = self.deepraccess.compute_atg_window_accessibility(
                full_tensor,
                atg_position=atg_position,
                discrete=True
            )
            

            if isinstance(recalculated_accessibility, torch.Tensor):
                recalculated_value = recalculated_accessibility.item()
            else:
                recalculated_value = recalculated_accessibility
            
            verification['recalculated_accessibility'] = recalculated_value
            

            if not (0.01 <= recalculated_value <= 15.0):
                verification['value_reasonable'] = False
                verification['issues'].append(

                )
            

            diff = abs(recalculated_value - original_accessibility)
            if diff > 1e-5:
                verification['consistent_with_original'] = False
                verification['issues'].append(

                    f"recalculated={recalculated_value:.6f}, diff={diff:.6f}"
                )
            
        except Exception as e:
            verification['calculation_success'] = False

        
        return verification
    
    def _extract_method(self, filename: str) -> str:
        """从文件名提取方法"""
        if 'lagrangian' in filename.lower():
            return 'lagrangian'
        elif 'ams' in filename.lower():
            return 'ams'
        elif 'cpc' in filename.lower():
            return 'cpc'
        return 'unknown'
    
    def _extract_variant(self, filename: str) -> str:

        for variant in ['00', '01', '10', '11']:
            if f'_{variant}_' in filename:
                return variant
        return 'unknown'
    
    def _extract_protein(self, filename: str) -> str:
        """从文件名提取蛋白质名称"""
        parts = filename.split('_')
        for part in parts:
            if part.startswith(('P', 'O')) and part[1:].isdigit():
                return part
        return 'unknown'
    
    def analyze_batch(self, batch_path: Path, max_files: int = 20) -> Dict:


        
        json_files = [f for f in batch_path.glob("*.json") 
                     if f.name not in ['summary.json', 'progress.json']]
        
        if not json_files:

            return {}
        

        if len(json_files) > max_files:
            json_files = json_files[:max_files]

        
        results = {
            'batch_name': batch_path.name,
            'analyzed_files': len(json_files),
            'verification_results': [],
            'summary': {
                'cds_issues': 0,
                'splicing_issues': 0,
                'tensor_issues': 0,
                'accessibility_issues': 0,
                'inconsistent_accessibility': 0
            }
        }
        
        for file_path in json_files:

            result = self.verify_sequence_reconstruction(file_path)
            results['verification_results'].append(result)
            

            if 'error' not in result:
                if not result['cds_verification']['sequence_valid']:
                    results['summary']['cds_issues'] += 1
                
                if not result['full_sequence_verification']['splicing_correct']:
                    results['summary']['splicing_issues'] += 1
                
                if not result['deepraccess_input_verification']['tensor_creation_success']:
                    results['summary']['tensor_issues'] += 1
                
                if not result['accessibility_verification']['calculation_success']:
                    results['summary']['accessibility_issues'] += 1
                
                if not result['accessibility_verification'].get('consistent_with_original', True):
                    results['summary']['inconsistent_accessibility'] += 1
        
        return results
    
    def print_analysis_report(self, results: Dict):
        """打印DeepRaccess输入验证报告"""
        print(f"\n📊 DeepRaccess输入验证报告: {results['batch_name']}")
        print("=" * 60)
        
        summary = results['summary']
        total = results['analyzed_files']
        
        print(f"分析文件数: {total}")
        
        print(f"\n🔍 发现的问题:")
        issues = [
            ('CDS序列问题', summary['cds_issues']),
            ('序列拼接问题', summary['splicing_issues']),
            ('张量创建问题', summary['tensor_issues']),
            ('Accessibility计算问题', summary['accessibility_issues']),
            ('Accessibility值不一致', summary['inconsistent_accessibility'])
        ]
        
        for issue_name, count in issues:
            if count > 0:
                rate = count / total * 100
                status = "🔴" if rate > 20 else "🟡" if rate > 0 else "🟢"
                print(f"  {issue_name}: {count}/{total} ({rate:.1f}%) {status}")
            else:
                print(f"  {issue_name}: 无问题 🟢")
        

        problem_cases = [r for r in results['verification_results'] 
                        if 'error' not in r and (
                            not r['cds_verification']['sequence_valid'] or
                            not r['full_sequence_verification']['splicing_correct'] or
                            not r['accessibility_verification'].get('consistent_with_original', True)
                        )]
        
        if problem_cases:
            print(f"\n❌ 问题案例 (前3个):")
            for case in problem_cases[:3]:
                issues = []
                if not case['cds_verification']['sequence_valid']:
                    issues.extend(case['cds_verification']['issues'])
                if not case['full_sequence_verification']['splicing_correct']:
                    issues.extend(case['full_sequence_verification']['issues'])
                if not case['accessibility_verification'].get('consistent_with_original', True):
                    issues.extend(case['accessibility_verification']['issues'])
                
                print(f"  {case['file']}: {case['method'].upper()}")
                for issue in issues[:2]:
                    print(f"    - {issue}")

def main():

    verifier = DeepRaccessInputVerifier()
    

    test_batch = Path("results/20250910_133850_unified_access_experiments")
    
    if test_batch.exists():

        results = verifier.analyze_batch(test_batch, max_files=15)
        verifier.print_analysis_report(results)
        

        output_file = Path("deepraccess_input_verification_report.json")
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False, default=str)

        
    else:


if __name__ == "__main__":
    main()
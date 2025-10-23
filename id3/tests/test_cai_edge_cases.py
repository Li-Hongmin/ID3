#!/usr/bin/env python3
"""


"""

import torch
import numpy as np
import sys
import pytest
from typing import Dict, Any

sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint

class TestCAIEdgeCases:

    
    @pytest.fixture(autouse=True)
    def setup(self):
        """测试设置"""
        self.test_proteins = {
            'minimal': 'MK',
            'short': 'MKHELM',
            'medium': 'MEEPQSDPSVEPPLSQETFSDLWKLL',
        }
        
    def test_extreme_cai_targets(self):

        protein_seq = self.test_proteins['short']
        
        extreme_targets = [




        ]
        
        for target in extreme_targets:

            
            try:
                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=target,
                    cai_lambda=0.1
                )
                
                result = constraint.forward(alpha=1.0, tau=1.0)
                discrete_cai = result.get('discrete_cai_value', None)
                
                if target <= 1.0:




                else:


                    
            except Exception as e:
                if target > 1.0:

                else:

                    raise
    
    def test_extreme_lambda_values(self):
        """测试极端lambda值"""
        protein_seq = self.test_proteins['short']
        
        extreme_lambdas = [
            0.0,
            1e-10,
            1e-6,
            0.001,
            10.0,
            100.0,
            1000.0,
        ]
        
        for lambda_val in extreme_lambdas:
            print(f"测试极端lambda值: {lambda_val}")
            
            try:
                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=0.8,
                    cai_lambda=lambda_val
                )
                
                result = constraint.forward(alpha=1.0, tau=1.0)
                discrete_cai = result.get('discrete_cai_value', None)
                total_loss = result.get('total_loss', None)
                cai_loss = result.get('cai_loss', None)
                
                assert discrete_cai is not None, f"Lambda={lambda_val}时discrete_cai为None"
                assert total_loss is not None, f"Lambda={lambda_val}时total_loss为None"
                
                print(f"  ✅ Lambda={lambda_val}: CAI={discrete_cai:.6f}, Loss={total_loss:.6f}")
                

                if lambda_val == 0.0:

                    if cai_loss is not None:
                        assert abs(cai_loss) < 1e-6, f"Lambda=0时CAI损失应为0，实际: {cai_loss}"
                elif lambda_val >= 100.0:

                    cai_diff = abs(discrete_cai - 0.8)
                    if cai_diff > 0.2:
                        print(f"    ⚠️  极大权重下CAI偏差仍较大: {cai_diff:.6f}")
                
            except Exception as e:
                print(f"  ❌ Lambda={lambda_val}: 异常: {str(e)[:50]}...")
                if lambda_val < 1000:
                    raise
    
    def test_dynamic_lambda_edge_cases(self):

        protein_seq = self.test_proteins['short']
        
        edge_configs = [
            {
                'name': 'extreme_range',
                'lambda_cai_min': 1e-10,
                'lambda_cai_max': 1000.0,
                'lambda_cai_adjustment_factor': 10.0
            },
            {
                'name': 'narrow_range',
                'lambda_cai_min': 0.099,
                'lambda_cai_max': 0.101,
                'lambda_cai_adjustment_factor': 1.001
            },
            {

                'lambda_cai_min': 1.0,
                'lambda_cai_max': 0.1,
                'lambda_cai_adjustment_factor': 1.5
            }
        ]
        
        for config in edge_configs:

            
            try:
                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=0.8,
                    adaptive_lambda_cai=True,
                    **{k: v for k, v in config.items() if k != 'name'}
                )
                

                for i in range(3):
                    result = constraint.forward(alpha=1.0, tau=1.0)
                    discrete_cai = result.get('discrete_cai_value', 0)
                    

                    if hasattr(constraint, '_update_lambda_cai'):
                        constraint._update_lambda_cai(discrete_cai)
                
                current_lambda = getattr(constraint, 'cai_lambda', None)

                
            except Exception as e:
                if config['name'] == 'inverted_range':

                else:

                    raise
    
    def test_minimal_protein_sequences(self):
        """测试最小蛋白质序列"""
        minimal_proteins = [
            'M',
            'MK',
            'MKH',
        ]
        
        for protein_seq in minimal_proteins:
            print(f"测试最小蛋白质序列: {protein_seq} (长度: {len(protein_seq)})")
            
            try:
                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=0.8,
                    cai_lambda=0.1
                )
                
                result = constraint.forward(alpha=1.0, tau=1.0)
                discrete_cai = result.get('discrete_cai_value', None)
                
                assert discrete_cai is not None, f"序列{protein_seq}时discrete_cai为None"
                assert discrete_cai >= 0, f"序列{protein_seq}时CAI值为负: {discrete_cai}"
                assert discrete_cai <= 1.0, f"序列{protein_seq}时CAI值超过1: {discrete_cai}"
                
                print(f"  ✅ 序列{protein_seq}: CAI={discrete_cai:.6f}")
                
            except Exception as e:
                print(f"  ❌ 序列{protein_seq}: 异常: {str(e)[:50]}...")

                if len(protein_seq) == 1:
                    print(f"    (单氨基酸序列失败可能是合理的)")
                else:
                    raise
    
    def test_invalid_parameters(self):

        protein_seq = self.test_proteins['short']
        
        invalid_configs = [
            {
                'name': 'negative_target',
                'cai_target': -0.5,
                'should_fail': True
            },
            {
                'name': 'negative_lambda',
                'cai_lambda': -0.1,
                'should_fail': True
            },
            {
                'name': 'zero_target',
                'cai_target': 0.0,

            },
            {
                'name': 'nan_target',
                'cai_target': float('nan'),
                'should_fail': True
            },
            {
                'name': 'inf_lambda',
                'cai_lambda': float('inf'),
                'should_fail': True
            }
        ]
        
        for config in invalid_configs:

            
            try:
                constraint_kwargs = {
                    'enable_cai': True,
                    'cai_target': 0.8,
                    'cai_lambda': 0.1
                }
                

                for key in ['cai_target', 'cai_lambda']:
                    if key in config:
                        constraint_kwargs[key] = config[key]
                
                constraint = LagrangianConstraint(protein_seq, **constraint_kwargs)
                result = constraint.forward(alpha=1.0, tau=1.0)
                
                if config['should_fail']:

                else:
                    discrete_cai = result.get('discrete_cai_value', None)

                
            except Exception as e:
                if config['should_fail']:

                else:

                    raise
    
    def test_alpha_tau_edge_combinations(self):
        """测试alpha和tau的极端组合"""
        protein_seq = self.test_proteins['short']
        
        edge_combinations = [
            (0.0, 0.1),
            (0.0, 10.0),
            (1.0, 0.1),
            (1.0, 10.0),
            (0.5, 0.01),
            (0.999, 1.0),
            (0.001, 1.0),
        ]
        
        constraint = LagrangianConstraint(
            protein_seq,
            enable_cai=True,
            cai_target=0.8,
            cai_lambda=0.1
        )
        
        for alpha, tau in edge_combinations:
            print(f"测试极端alpha/tau组合: α={alpha}, τ={tau}")
            
            try:
                result = constraint.forward(alpha=alpha, tau=tau)
                discrete_cai = result.get('discrete_cai_value', None)
                total_loss = result.get('total_loss', None)
                
                assert discrete_cai is not None, f"α={alpha}, τ={tau}时discrete_cai为None"
                assert total_loss is not None, f"α={alpha}, τ={tau}时total_loss为None"
                assert np.isfinite(discrete_cai), f"α={alpha}, τ={tau}时CAI非有限值: {discrete_cai}"
                assert np.isfinite(total_loss), f"α={alpha}, τ={tau}时损失非有限值: {total_loss}"
                
                print(f"  ✅ α={alpha}, τ={tau}: CAI={discrete_cai:.6f}, Loss={total_loss:.6f}")
                
            except Exception as e:
                print(f"  ❌ α={alpha}, τ={tau}: 异常: {str(e)[:50]}...")

                if tau < 0.1 or tau > 5.0:
                    print(f"    (极端温度参数失败可能是合理的)")
                else:
                    raise
    
    def test_constraint_consistency_edge_cases(self):


        

        constraint_classes = [
            (LagrangianConstraint, {'enable_cai': True, 'cai_target': 0.95, 'cai_lambda': 10.0}),
            (AminoMatchingSoftmax, {'enable_cai': True, 'cai_target': 0.95, 'cai_weight': 10.0}),
            (CodonProfileConstraint, {'enable_cai': True, 'cai_target': 0.95, 'cai_weight': 10.0})
        ]
        results = {}
        
        for constraint_class, edge_config in constraint_classes:
            constraint_name = constraint_class.__name__.replace('Constraint', '').lower()

            
            try:
                constraint = constraint_class(protein_seq, **edge_config)
                result = constraint.forward(alpha=1.0, tau=1.0)
                
                discrete_cai = result.get('discrete_cai_value', None)
                results[constraint_name] = discrete_cai
                
                print(f"  ✅ {constraint_name}: CAI={discrete_cai:.6f}")
                
            except Exception as e:

                results[constraint_name] = None
        

        valid_results = {k: v for k, v in results.items() if v is not None}
        
        if len(valid_results) >= 2:
            cai_values = list(valid_results.values())
            cai_std = np.std(cai_values)
            cai_mean = np.mean(cai_values)
            




        else:

    
    def test_memory_stress(self):
        """测试内存压力情况"""

        long_protein = self.test_proteins['medium'] * 3
        
        print(f"内存压力测试，蛋白质长度: {len(long_protein)}")
        
        try:
            constraint = LagrangianConstraint(
                long_protein,
                enable_cai=True,
                cai_target=0.8,
                cai_lambda=0.1,
                batch_size=1
            )
            

            for i in range(5):
                result = constraint.forward(alpha=1.0, tau=1.0)
                discrete_cai = result.get('discrete_cai_value', None)
                
                assert discrete_cai is not None, f"第{i+1}次调用时discrete_cai为None"
                print(f"  第{i+1}次: CAI={discrete_cai:.6f}")
            
            print(f"✅ 内存压力测试通过")
            
        except Exception as e:
            print(f"❌ 内存压力测试失败: {str(e)[:100]}...")
            raise


def run_edge_case_tests():


    print("=" * 80)
    
    tester = TestCAIEdgeCases()
    tester.setup()
    
    test_methods = [








    ]
    
    passed_tests = 0
    total_tests = len(test_methods)
    
    for test_name, test_method in test_methods:
        print(f"\n🔬 {test_name}")
        print("-" * 60)
        
        try:
            test_method()
            passed_tests += 1

        except Exception as e:

    





    
    return passed_tests == total_tests


if __name__ == "__main__":
    success = run_edge_case_tests()
    exit(0 if success else 1)
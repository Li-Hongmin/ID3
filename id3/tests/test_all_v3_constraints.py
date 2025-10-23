#!/usr/bin/env python3
"""



- CPC v3 Stable  
- AMS v3 Stable
- Lagrangian v3 Stable









"""

import sys
import os
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

import torch
import time
import logging
import numpy as np
from typing import Dict, List, Tuple, Any
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


try:
    from id3.constraints.cpc_v3_stable import CodonProfileConstraintV3Stable, create_cpc_v3_constraint
    from id3.constraints.ams_v3_stable import AdaptiveMultilevelSamplingV3Stable, create_ams_v3_constraint  
    from id3.constraints.lagrangian_v3_stable import LagrangianConstraintV3Stable, create_lagrangian_v3_constraint
    from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
except ImportError as e:
    logger.error(f"V3模块导入失败: {e}")
    sys.exit(1)


class UnifiedV3ConstraintTester:

    
    def __init__(self, device: str = 'cuda'):
        self.device = device
        self.test_sequences = {
            'short': "MKLLILTLSCLEIALVETTAIIKYSPKEEAE",  # 30 AA
            'medium': "MKLLILTLSCLEIALVETTAIIKYSPKEEAEGIAPAIPMNGAAPDVMDAYKKYVAYVDCPQSHKDQLSVLKKALEGYVKSDKKYHVPTIIVQAVKDLLETLAADMQTLPQGRRPEELPTVVVEKALRRLGDTGDVLIITDPSQKPGAFVAELRSDGT",  # 150 AA  
            'long': "MKLLILTLSCLEIALVETTAIIKYSPKEEAEGIAPAIPMNGAAPDVMDAYKKYVAYVDCPQSHKDQLSVLKKALEGYVKSDKKYHVPTIIVQAVKDLLETLAADMQTLPQGRRPEELPTVVVEKALRRLGDTGDVLIITDPSQKPGAFVAELRSDGTVAEVIQKDGVPSRGDVLIMSSLGLPYRGEYAAVSNDGPTAPGFQYALVDVNSHEVLRYEGHTYETYCSPDQPTMEELEELWRTFGALHAWHVPVPDIEEAWALLSLGPQEEEPLPEVRALPELKFLFTQIFQLQFMAAFSFLHVKDPIGIMVQGGGPQALSDQRRALMVRLGHAAGMLNPQVRRLLGQQQSTLLAAFRSISQVLTVLPQILSLVQVVQGLIQMIVFLVIIQIRQMVDTSLSQMVDQSLEAMLVFHCQASDDMDGGRTLSHGYGSYDNVCQVTTFLGMFLIKVLGQLLHRYKAADDHVLDRILSEFIEQDLEDDERAVSSGQLSATMIYAAGFQQDKKNAASGSADAAQDKAASEAEEAKKEKQEEEHEEEAAAAESEPPPPAPAEQAAAVSSGALSGDAHAQEAVADRLFVNQSLQAALLSSGIPQEAQSFASMGMVQSGVSYRGPLYTGGQMVHVDLGSVQESASMTGGVTQQGNMYMQRGQLGGLYRGLQPQVGHSYMFNEPSLPQVEAIAPVMVQRGMPSAQALMVQMDLDVMDVVDANALTQGDHVLTVQGSHVLMVHDSDFSCAQVVSGGLQFQVQSSQCVQ"  # 755 AA
        }
        
        self.constraint_configs = {
            'CPC': {
                'factory': create_cpc_v3_constraint,
                'class': CodonProfileConstraintV3Stable,
                'variants': ['standard']
            },
            'AMS': {
                'factory': create_ams_v3_constraint, 
                'class': AdaptiveMultilevelSamplingV3Stable,
                'variants': ['standard']
            },
            'Lagrangian': {
                'factory': create_lagrangian_v3_constraint,
                'class': LagrangianConstraintV3Stable, 
                'variants': ['L00', 'L01', 'L10', 'L11']
            }
        }
        
        self.test_results = {}
        
    def setup_deepraccess(self) -> Any:
        """初始化DeepRaccess模型"""
        logger.info("初始化DeepRaccess模型...")
        try:
            deepraccess_wrapper = DeepRaccessID3Wrapper(device=self.device)
            deepraccess_wrapper.deepraccess_model.eval()
            logger.info("DeepRaccess模型初始化完成")
            return deepraccess_wrapper
        except Exception as e:
            logger.error(f"DeepRaccess初始化失败: {e}")
            raise
    
    def test_constraint_creation(self, deepraccess_wrapper) -> Dict:


        
        results = {}
        sequence = self.test_sequences['short']
        
        for constraint_type, config in self.constraint_configs.items():

            results[constraint_type] = {}
            
            for variant in config['variants']:
                try:

                    if constraint_type == 'Lagrangian':
                        constraint = config['factory'](
                            amino_acid_sequence=sequence,
                            deepraccess_model=deepraccess_wrapper,
                            device=self.device,
                            variant=variant
                        )
                    else:
                        constraint = config['factory'](
                            amino_acid_sequence=sequence,
                            deepraccess_model=deepraccess_wrapper,
                            device=self.device
                        )
                    

                    result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                    
                    results[constraint_type][variant] = {
                        'creation_success': True,
                        'forward_success': True,
                        'total_loss': result['loss'].item(),
                        'accessibility_loss': result['accessibility_loss'].item(),
                        'cai_loss': result.get('cai_loss', torch.tensor(0.0)).item(),
                        'forward_time_ms': result.get('forward_time_ms', 0.0)
                    }
                    

                    
                except Exception as e:

                    results[constraint_type][variant] = {
                        'creation_success': False,
                        'error': str(e)
                    }
        
        return results
    
    def test_zero_return_stability(self, deepraccess_wrapper, iterations: int = 10) -> Dict:
        """测试DeepRaccess零值返回问题的解决情况"""
        logger.info(f"测试零值返回稳定性 (迭代{iterations}次)...")
        
        results = {}
        sequence = self.test_sequences['medium']
        
        for constraint_type, config in self.constraint_configs.items():
            logger.info(f"测试 {constraint_type} 零值返回问题...")
            results[constraint_type] = {}
            

            test_variant = config['variants'][0] if constraint_type != 'Lagrangian' else 'L11'
            
            try:

                if constraint_type == 'Lagrangian':
                    constraint = config['factory'](
                        amino_acid_sequence=sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device,
                        variant=test_variant
                    )
                else:
                    constraint = config['factory'](
                        amino_acid_sequence=sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device
                    )
                

                zero_count = 0
                total_calls = 0
                execution_times = []
                accessibility_values = []
                
                for i in range(iterations):
                    try:
                        start_time = time.time()
                        result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                        execution_time = time.time() - start_time
                        
                        accessibility_loss = result['accessibility_loss'].item()
                        total_calls += 1
                        execution_times.append(execution_time)
                        accessibility_values.append(accessibility_loss)
                        

                        if abs(accessibility_loss) < 1e-8:
                            zero_count += 1
                            logger.warning(f"{constraint_type}-{test_variant} 迭代 {i}: 零值 {accessibility_loss}")
                        
                    except Exception as e:
                        logger.error(f"{constraint_type}-{test_variant} 迭代 {i} 失败: {e}")
                
                success_rate = (total_calls - zero_count) / max(total_calls, 1)
                results[constraint_type][test_variant] = {
                    'total_calls': total_calls,
                    'zero_count': zero_count,
                    'success_rate': success_rate,
                    'avg_execution_time': np.mean(execution_times),
                    'std_execution_time': np.std(execution_times),
                    'min_accessibility': np.min(accessibility_values),
                    'max_accessibility': np.max(accessibility_values),
                    'mean_accessibility': np.mean(accessibility_values)
                }
                
                logger.info(f"{constraint_type}-{test_variant}: {success_rate:.2%} 成功率, 零值{zero_count}/{total_calls}")
                
            except Exception as e:
                logger.error(f"{constraint_type} 零值测试失败: {e}")
                results[constraint_type][test_variant] = {'error': str(e)}
        
        return results
    
    def test_cai_functionality(self, deepraccess_wrapper) -> Dict:


        
        results = {}
        sequence = self.test_sequences['short']
        
        for constraint_type, config in self.constraint_configs.items():

            results[constraint_type] = {}
            

            for enable_cai in [False, True]:
                test_key = f"CAI_{'enabled' if enable_cai else 'disabled'}"
                
                try:

                    if constraint_type == 'Lagrangian':
                        variant = 'L01' if enable_cai else 'L00'
                        constraint = config['factory'](
                            amino_acid_sequence=sequence,
                            deepraccess_model=deepraccess_wrapper,
                            device=self.device,
                            variant=variant
                        )
                    else:
                        constraint = config['factory'](
                            amino_acid_sequence=sequence,
                            deepraccess_model=deepraccess_wrapper,
                            enable_cai=enable_cai,
                            device=self.device
                        )
                    
                    result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                    
                    results[constraint_type][test_key] = {
                        'success': True,
                        'total_loss': result['loss'].item(),
                        'accessibility_loss': result['accessibility_loss'].item(),
                        'cai_loss': result.get('cai_loss', torch.tensor(0.0)).item(),
                        'cai_expected': enable_cai
                    }
                    

                    cai_loss = results[constraint_type][test_key]['cai_loss']
                    if enable_cai and cai_loss <= 1e-8:

                    elif not enable_cai and cai_loss > 1e-8:

                    

                    
                except Exception as e:

                    results[constraint_type][test_key] = {
                        'success': False,
                        'error': str(e)
                    }
        
        return results
    
    def test_performance_comparison(self, deepraccess_wrapper) -> Dict:
        """性能对比测试"""
        logger.info("执行性能对比测试...")
        
        results = {}
        sequence = self.test_sequences['medium']
        iterations = 5
        
        for constraint_type, config in self.constraint_configs.items():
            logger.info(f"测试 {constraint_type} 性能...")
            

            test_variant = config['variants'][0] if constraint_type != 'Lagrangian' else 'L11'
            
            try:

                if constraint_type == 'Lagrangian':
                    constraint = config['factory'](
                        amino_acid_sequence=sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device,
                        variant=test_variant
                    )
                else:
                    constraint = config['factory'](
                        amino_acid_sequence=sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device
                    )
                

                execution_times = []
                
                for i in range(iterations):
                    start_time = time.time()
                    result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                    execution_time = time.time() - start_time
                    execution_times.append(execution_time)
                

                performance_stats = {}
                if hasattr(constraint, 'get_performance_stats'):
                    performance_stats = constraint.get_performance_stats()
                
                results[constraint_type] = {
                    'variant': test_variant,
                    'avg_execution_time': np.mean(execution_times),
                    'std_execution_time': np.std(execution_times),
                    'min_execution_time': np.min(execution_times),
                    'max_execution_time': np.max(execution_times),
                    'iterations': iterations,
                    'performance_stats': performance_stats
                }
                
                logger.info(f"{constraint_type} 性能: 平均{results[constraint_type]['avg_execution_time']:.4f}s")
                
            except Exception as e:
                logger.error(f"{constraint_type} 性能测试失败: {e}")
                results[constraint_type] = {'error': str(e)}
        
        return results
    
    def test_long_sequence_stability(self, deepraccess_wrapper) -> Dict:


        
        results = {}
        long_sequence = self.test_sequences['long']  # 755 AA

        
        for constraint_type, config in self.constraint_configs.items():

            
            test_variant = config['variants'][0] if constraint_type != 'Lagrangian' else 'L11'
            
            try:

                if constraint_type == 'Lagrangian':
                    constraint = config['factory'](
                        amino_acid_sequence=long_sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device,
                        variant=test_variant
                    )
                else:
                    constraint = config['factory'](
                        amino_acid_sequence=long_sequence,
                        deepraccess_model=deepraccess_wrapper,
                        device=self.device
                    )
                

                execution_times = []
                zero_count = 0
                total_calls = 0
                
                for i in range(iterations):
                    try:
                        start_time = time.time()
                        result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                        execution_time = time.time() - start_time
                        
                        accessibility_loss = result['accessibility_loss'].item()
                        total_calls += 1
                        execution_times.append(execution_time)
                        
                        if abs(accessibility_loss) < 1e-8:
                            zero_count += 1
                            
                    except Exception as e:

                
                success_rate = (total_calls - zero_count) / max(total_calls, 1)
                results[constraint_type] = {
                    'sequence_length': len(long_sequence),
                    'total_calls': total_calls,
                    'zero_count': zero_count, 
                    'success_rate': success_rate,
                    'avg_execution_time': np.mean(execution_times) if execution_times else 0,
                    'max_execution_time': np.max(execution_times) if execution_times else 0
                }
                

                
            except Exception as e:

                results[constraint_type] = {'error': str(e)}
        
        return results
    
    def run_comprehensive_test(self) -> Dict:
        """运行全面的v3约束测试"""
        logger.info("开始V3约束全面测试套件...")
        
        comprehensive_results = {
            'test_timestamp': datetime.now().isoformat(),
            'device': str(self.device),
            'gpu_available': torch.cuda.is_available()
        }
        
        try:

            deepraccess_wrapper = self.setup_deepraccess()
            

            logger.info("=" * 60)
            comprehensive_results['constraint_creation'] = self.test_constraint_creation(deepraccess_wrapper)
            

            logger.info("=" * 60)
            comprehensive_results['zero_return_stability'] = self.test_zero_return_stability(deepraccess_wrapper, iterations=10)
            

            logger.info("=" * 60)
            comprehensive_results['cai_functionality'] = self.test_cai_functionality(deepraccess_wrapper)
            

            logger.info("=" * 60)
            comprehensive_results['performance_comparison'] = self.test_performance_comparison(deepraccess_wrapper)
            

            logger.info("=" * 60)
            comprehensive_results['long_sequence_stability'] = self.test_long_sequence_stability(deepraccess_wrapper)
            
            logger.info("V3约束全面测试完成！")
            
        except Exception as e:
            logger.error(f"测试过程中出现错误: {e}")
            comprehensive_results['error'] = str(e)
        
        return comprehensive_results
    
    def generate_comprehensive_report(self, results: Dict) -> str:

        report = []

        report.append("=" * 80)



        report.append("")
        

        if 'constraint_creation' in results:

            report.append("")
            creation = results['constraint_creation']
            
            for constraint_type in ['CPC', 'AMS', 'Lagrangian']:
                if constraint_type in creation:

                    for variant, result in creation[constraint_type].items():
                        if result.get('creation_success', False):





                        else:

                    report.append("")
        

        if 'zero_return_stability' in results:

            report.append("")
            stability = results['zero_return_stability']
            
            for constraint_type in ['CPC', 'AMS', 'Lagrangian']:
                if constraint_type in stability:
                    for variant, result in stability[constraint_type].items():
                        if 'error' not in result:
                            report.append(f"### {constraint_type}-{variant}")




                            report.append("")
                        else:

                            report.append("")
        

        if 'cai_functionality' in results:

            report.append("")
            cai_results = results['cai_functionality']
            
            for constraint_type in ['CPC', 'AMS', 'Lagrangian']:
                if constraint_type in cai_results:

                    for test_key, result in cai_results[constraint_type].items():
                        if result.get('success', False):




                        else:

                    report.append("")
        

        if 'performance_comparison' in results:

            report.append("")
            performance = results['performance_comparison']
            
            for constraint_type in ['CPC', 'AMS', 'Lagrangian']:
                if constraint_type in performance and 'error' not in performance[constraint_type]:
                    result = performance[constraint_type]





                    

                    if 'performance_stats' in result and result['performance_stats']:
                        stats = result['performance_stats']

                        if 'total_calls' in stats:

                        if 'success_rate' in stats:

                    report.append("")
        

        if 'long_sequence_stability' in results:

            report.append("")
            long_seq = results['long_sequence_stability']
            
            for constraint_type in ['CPC', 'AMS', 'Lagrangian']:
                if constraint_type in long_seq and 'error' not in long_seq[constraint_type]:
                    result = long_seq[constraint_type]






                    report.append("")
        


        report.append("")
        

        all_success = True
        constraint_count = 0
        successful_constraints = 0
        
        if 'constraint_creation' in results:
            for constraint_type, variants in results['constraint_creation'].items():
                for variant, result in variants.items():
                    constraint_count += 1
                    if result.get('creation_success', False):
                        successful_constraints += 1
                    else:
                        all_success = False
        
        if 'zero_return_stability' in results:
            for constraint_type, variants in results['zero_return_stability'].items():
                for variant, result in variants.items():
                    if 'error' in result or result.get('success_rate', 0) < 1.0:
                        all_success = False
        

        if constraint_count > 0:
            creation_rate = successful_constraints / constraint_count

        
        if all_success:

            report.append("")






        else:

        
        report.append("")

        
        return "\n".join(report)


def main():
    """主测试函数"""
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    logger.info(f"使用设备: {device}")
    logger.info("开始ID3 V3约束统一测试...")
    

    tester = UnifiedV3ConstraintTester(device=device)
    

    results = tester.run_comprehensive_test()
    

    report = tester.generate_comprehensive_report(results)
    

    print("\n" + report)
    

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"v3_constraints_test_report_{timestamp}.txt"
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    logger.info(f"详细测试报告已保存到: {report_file}")
    
    return results


if __name__ == "__main__":
    results = main()
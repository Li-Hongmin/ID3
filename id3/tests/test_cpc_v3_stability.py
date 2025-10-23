#!/usr/bin/env python3
"""








"""

import sys
import os
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

import torch
import time
import logging
import numpy as np
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import seaborn as sns


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


try:
    from id3.constraints.cpc_v3_stable import CodonProfileConstraintV3Stable, create_cpc_v3_constraint
    from id3.constraints.cpc_v2_efficient import CodonProfileConstraintV2Efficient
    from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
except ImportError as e:
    logger.error(f"模块导入失败: {e}")
    sys.exit(1)


class CPCStabilityTester:

    
    def __init__(self, device: str = 'cuda'):
        self.device = device
        self.test_sequences = [
            "MKLLILTLSCLEIALVETTAIIKYSPKEEAEEGIAPAIPMNGAAPDVMDAYKKYVAYVDCPQSHKDQLSVLKKALEGYVKSDKKYHVPTIIVQAVKDLLETLAADMQTLPQGRRPEELPTVVVEKALRRLGDTGDVLIITDPSQKPGAFVAELRSDGTVAEVIQKDGVPSRGDVLIMSSLGLPYRGEYAAVSNDGPTAPGFQYALVDVNSHEVLRYEGHTYETYCSPDQPTMEELEELWRTFGALHAWHVPVPDIEEAWALLSLGPQEEEPLPEVRALPELKFLFTQIFQLQFMAAFSFLHVKDPIGIMVQGGGPQALSDQRRALMVRLGHAAGMLNPQVRRLLGQQQSTLLAAFRSISQVLTVLPQILSLVQVVQGLIQMIVFLVIIQIRQMVDTSLSQMVDQSLEAMLVFHCQASDDMDGGRTLSHGYGSYDNVCQVTTFLGMFLIKVLGQLLHRYKAADDHVLDRILSEFIEQDLEDDERAVSSGQLSATMIYAAGFQQDKKNAASGSADAAQDKAASEAEEAKKEKQEEEHEEEAAAAESEPPPPAPAEQAAAVSSGALSGDAHAQEAVADRLFVNQSLQAALLSSGIPQEAQSFASMGMVQSGVSYRGPLYTGGQMVHVDLGSVQESASMTGGVTQQGNMYMQRGQLGGLYRGLQPQVGHSYMFNEPSLPQVEAIAPVMVQRGMPSAQALMVQMDLDVMDVVDANALTQGDHVLTVQGSHVLMVHDSDFSCAQVVSGGLQFQVQSSQCVQ"
        ]
        
        self.test_results = {
            'v2_results': [],
            'v3_results': [],
            'performance_metrics': {}
        }
    
    def setup_models(self) -> Tuple:
        """设置测试用的模型"""
        logger.info("初始化DeepRaccess模型...")
        
        try:

            deepraccess_wrapper = DeepRaccessID3Wrapper(device=self.device)
            deepraccess_wrapper.deepraccess_model.eval()
            
            logger.info("DeepRaccess模型初始化完成")
            return deepraccess_wrapper
            
        except Exception as e:
            logger.error(f"模型初始化失败: {e}")
            raise
    
    def test_zero_return_issue(self, constraint, test_name: str, iterations: int = 10) -> Dict:


        
        results = {
            'zero_count': 0,
            'total_calls': 0,
            'accessibility_values': [],
            'execution_times': [],
            'success_rate': 0.0,
            'errors': []
        }
        
        for i in range(iterations):
            try:
                start_time = time.time()
                

                result = constraint.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
                
                execution_time = time.time() - start_time
                results['execution_times'].append(execution_time)
                results['total_calls'] += 1
                

                accessibility_loss = result['accessibility_loss'].item()
                results['accessibility_values'].append(accessibility_loss)
                

                if abs(accessibility_loss) < 1e-8:
                    results['zero_count'] += 1

                else:

                
            except Exception as e:

                results['errors'].append(str(e))
        

        if results['total_calls'] > 0:
            results['success_rate'] = (results['total_calls'] - results['zero_count']) / results['total_calls']
        

        
        return results
    
    def test_performance_comparison(self, deepraccess_wrapper) -> Dict:
        """性能对比测试"""
        logger.info("开始性能对比测试...")
        
        sequence = self.test_sequences[0][:50]
        

        logger.info("测试v2版本性能...")
        try:
            constraint_v2 = CodonProfileConstraintV2Efficient(
                amino_acid_sequence=sequence,
                deepraccess_model=deepraccess_wrapper,
                device=self.device
            )
            v2_results = self.test_zero_return_issue(constraint_v2, "v2", iterations=5)
        except Exception as e:
            logger.error(f"v2版本测试失败: {e}")
            v2_results = {'error': str(e)}
        

        logger.info("测试v3版本性能...")
        try:
            constraint_v3 = create_cpc_v3_constraint(
                amino_acid_sequence=sequence,
                deepraccess_model=deepraccess_wrapper,
                device=self.device
            )
            v3_results = self.test_zero_return_issue(constraint_v3, "v3", iterations=5)
            

            v3_performance = constraint_v3.get_performance_stats()
            v3_results['detailed_performance'] = v3_performance
            
        except Exception as e:
            logger.error(f"v3版本测试失败: {e}")
            v3_results = {'error': str(e)}
        
        return {
            'v2_results': v2_results,
            'v3_results': v3_results
        }
    
    def test_cai_functionality(self, deepraccess_wrapper) -> Dict:


        
        sequence = self.test_sequences[0][:30]
        

        constraint_with_cai = create_cpc_v3_constraint(
            amino_acid_sequence=sequence,
            deepraccess_model=deepraccess_wrapper,
            enable_cai=True,
            cai_target=0.8,
            lambda_cai=0.1,
            device=self.device
        )
        

        constraint_without_cai = create_cpc_v3_constraint(
            amino_acid_sequence=sequence,
            deepraccess_model=deepraccess_wrapper,
            enable_cai=False,
            device=self.device
        )
        
        results = {}
        

        try:
            result_with_cai = constraint_with_cai.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
            results['with_cai'] = {
                'total_loss': result_with_cai['loss'].item(),
                'accessibility_loss': result_with_cai['accessibility_loss'].item(),
                'cai_loss': result_with_cai['cai_loss'].item(),
                'success': True
            }

        except Exception as e:
            results['with_cai'] = {'success': False, 'error': str(e)}

        

        try:
            result_without_cai = constraint_without_cai.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
            results['without_cai'] = {
                'total_loss': result_without_cai['loss'].item(),
                'accessibility_loss': result_without_cai['accessibility_loss'].item(),
                'cai_loss': result_without_cai['cai_loss'].item(),
                'success': True
            }

        except Exception as e:
            results['without_cai'] = {'success': False, 'error': str(e)}

        
        return results
    
    def test_memory_efficiency(self, deepraccess_wrapper) -> Dict:
        """测试内存使用效率"""
        logger.info("测试内存使用效率...")
        
        if not torch.cuda.is_available():
            logger.warning("CUDA不可用，跳过GPU内存测试")
            return {'gpu_available': False}
        
        sequence = self.test_sequences[0][:100]
        

        torch.cuda.empty_cache()
        initial_memory = torch.cuda.memory_allocated()
        

        constraint_v3 = create_cpc_v3_constraint(
            amino_acid_sequence=sequence,
            deepraccess_model=deepraccess_wrapper,
            device=self.device
        )
        
        creation_memory = torch.cuda.memory_allocated()
        

        max_memory = creation_memory
        for i in range(10):
            result = constraint_v3.forward_with_loss(alpha=0.0, tau=1.0, beta=0.0)
            current_memory = torch.cuda.memory_allocated()
            max_memory = max(max_memory, current_memory)
        
        final_memory = torch.cuda.memory_allocated()
        
        return {
            'gpu_available': True,
            'initial_memory_mb': initial_memory / 1024 / 1024,
            'creation_memory_mb': creation_memory / 1024 / 1024,
            'max_memory_mb': max_memory / 1024 / 1024,
            'final_memory_mb': final_memory / 1024 / 1024,
            'memory_increase_mb': (final_memory - initial_memory) / 1024 / 1024,
            'peak_increase_mb': (max_memory - initial_memory) / 1024 / 1024
        }
    
    def run_comprehensive_test(self) -> Dict:


        
        test_results = {}
        
        try:

            deepraccess_wrapper = self.setup_models()
            

            logger.info("=" * 50)
            performance_results = self.test_performance_comparison(deepraccess_wrapper)
            test_results['performance_comparison'] = performance_results
            

            logger.info("=" * 50)
            cai_results = self.test_cai_functionality(deepraccess_wrapper)
            test_results['cai_functionality'] = cai_results
            

            logger.info("=" * 50)
            memory_results = self.test_memory_efficiency(deepraccess_wrapper)
            test_results['memory_efficiency'] = memory_results
            

            logger.info("=" * 50)


            constraint_long = create_cpc_v3_constraint(
                amino_acid_sequence=long_sequence,
                deepraccess_model=deepraccess_wrapper,
                device=self.device
            )
            
            long_seq_results = self.test_zero_return_issue(constraint_long, "v3_long_sequence", iterations=3)
            test_results['long_sequence_stability'] = long_seq_results
            

            
        except Exception as e:

            test_results['error'] = str(e)
        
        return test_results
    
    def generate_report(self, results: Dict) -> str:
        """生成测试报告"""
        report = []
        report.append("CPC v3稳定性测试报告")
        report.append("=" * 50)
        report.append("")
        

        if 'performance_comparison' in results:
            perf = results['performance_comparison']
            report.append("
            report.append("")
            
            if 'error' not in perf['v2_results']:
                v2 = perf['v2_results']
                report.append(f"v2版本:")
                report.append(f"  - 成功率: {v2.get('success_rate', 0):.2%}")
                report.append(f"  - 零值返回: {v2.get('zero_count', 0)}/{v2.get('total_calls', 0)}")
                report.append(f"  - 平均执行时间: {np.mean(v2.get('execution_times', [0])):.4f}s")
            else:
                report.append(f"v2版本: 测试失败 - {perf['v2_results']['error']}")
            
            if 'error' not in perf['v3_results']:
                v3 = perf['v3_results']
                report.append(f"v3版本:")
                report.append(f"  - 成功率: {v3.get('success_rate', 0):.2%}")
                report.append(f"  - 零值返回: {v3.get('zero_count', 0)}/{v3.get('total_calls', 0)}")
                report.append(f"  - 平均执行时间: {np.mean(v3.get('execution_times', [0])):.4f}s")
                
                if 'detailed_performance' in v3:
                    stats = v3['detailed_performance']
                    report.append(f"  - DeepRaccess调用: {stats.get('total_calls', 0)}")
                    report.append(f"  - BatchNorm重置次数: {stats.get('reset_count', 0)}")
                    report.append(f"  - 平均推理时间: {stats.get('avg_inference_time_ms', 0):.2f}ms")
            else:
                report.append(f"v3版本: 测试失败 - {perf['v3_results']['error']}")
            
            report.append("")
        

        if 'cai_functionality' in results:
            cai = results['cai_functionality']
            report.append("
            report.append("")
            
            if cai.get('with_cai', {}).get('success'):
                with_cai = cai['with_cai']
                report.append(f"CAI启用:")
                report.append(f"  - 总损失: {with_cai['total_loss']:.6f}")
                report.append(f"  - 可及性损失: {with_cai['accessibility_loss']:.6f}")
                report.append(f"  - CAI损失: {with_cai['cai_loss']:.6f}")
            
            if cai.get('without_cai', {}).get('success'):
                without_cai = cai['without_cai']
                report.append(f"CAI禁用:")
                report.append(f"  - 总损失: {without_cai['total_loss']:.6f}")
                report.append(f"  - 可及性损失: {without_cai['accessibility_loss']:.6f}")
                report.append(f"  - CAI损失: {without_cai['cai_loss']:.6f}")
            
            report.append("")
        

        if 'memory_efficiency' in results:
            mem = results['memory_efficiency']
            report.append("
            report.append("")
            
            if mem.get('gpu_available'):
                report.append(f"GPU内存使用:")
                report.append(f"  - 初始内存: {mem['initial_memory_mb']:.2f} MB")
                report.append(f"  - 创建后内存: {mem['creation_memory_mb']:.2f} MB")
                report.append(f"  - 峰值内存: {mem['max_memory_mb']:.2f} MB")
                report.append(f"  - 最终内存: {mem['final_memory_mb']:.2f} MB")
                report.append(f"  - 内存增长: {mem['memory_increase_mb']:.2f} MB")
                report.append(f"  - 峰值增长: {mem['peak_increase_mb']:.2f} MB")
            else:
                report.append("GPU内存测试不可用")
            
            report.append("")
        

        if 'long_sequence_stability' in results:
            long_seq = results['long_sequence_stability']
            report.append("
            report.append("")
            report.append(f"  - 成功率: {long_seq.get('success_rate', 0):.2%}")
            report.append(f"  - 零值返回: {long_seq.get('zero_count', 0)}/{long_seq.get('total_calls', 0)}")
            report.append(f"  - 平均执行时间: {np.mean(long_seq.get('execution_times', [0])):.4f}s")
            report.append("")
        

        report.append("
        report.append("")
        

        all_passed = True
        if 'performance_comparison' in results:
            perf = results['performance_comparison']
            if 'error' in perf['v3_results'] or perf['v3_results'].get('success_rate', 0) < 1.0:
                all_passed = False
        
        if 'cai_functionality' in results:
            cai = results['cai_functionality']
            if not (cai.get('with_cai', {}).get('success') and cai.get('without_cai', {}).get('success')):
                all_passed = False
        
        if 'long_sequence_stability' in results:
            if results['long_sequence_stability'].get('success_rate', 0) < 1.0:
                all_passed = False
        
        if all_passed:
            report.append("✅ 所有测试通过！CPC v3版本稳定可用。")
        else:
            report.append("❌ 部分测试未通过，需要进一步优化。")
        
        report.append("")
        report.append("测试完成时间: " + time.strftime('%Y-%m-%d %H:%M:%S'))
        
        return "\n".join(report)


def main():


    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    

    tester = CPCStabilityTester(device=device)
    

    results = tester.run_comprehensive_test()
    

    report = tester.generate_report(results)
    

    print("\n" + report)
    

    timestamp = time.strftime('%Y%m%d_%H%M%S')
    report_file = f"cpc_v3_test_report_{timestamp}.txt"
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    

    

    return results


if __name__ == "__main__":
    results = main()
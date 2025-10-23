#!/usr/bin/env python3
"""
Comprehensive CAI Fixes Validation Framework







"""

import torch
import numpy as np
import json
import time
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass

sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('comprehensive_cai_validation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class TestResult:

    test_name: str
    success: bool
    cai_value: Optional[float]
    error_message: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
    performance_metrics: Optional[Dict[str, float]] = None

class CAIValidationFramework:
    """CAI修复验证框架"""
    
    def __init__(self):
        self.results: List[TestResult] = []
        self.test_proteins = {
            'short': 'MKHELM',
            'medium': 'MEEPQSDPSVEPPLSQETFSDLWKLL',
            'long': 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMMDDLMLSPDDIEQWFTEDPGPDEAP',
        }
        self.test_constraints = ['lagrangian', 'amino_matching', 'codon_profile']
        
    def log_test_start(self, test_name: str):


        print(f"\n{'='*80}")

        print(f"{'='*80}")
        
    def log_test_result(self, result: TestResult):
        """记录测试结果"""
        self.results.append(result)
        status = "✅ 成功" if result.success else "❌ 失败"
        logger.info(f"测试结果 [{result.test_name}]: {status}")
        if result.cai_value is not None:
            logger.info(f"  CAI值: {result.cai_value:.6f}")
        if result.error_message:
            logger.error(f"  错误: {result.error_message}")
            
    def test_lagrangian_cai_recovery(self) -> List[TestResult]:


        
        test_results = []
        
        for protein_name, protein_seq in self.test_proteins.items():

            
            try:

                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=0.8,
                    cai_lambda=0.1,
                    verbose=False
                )
                

                start_time = time.time()

                end_time = time.time()
                

                discrete_cai = result.get('discrete_cai_value', None)
                cai_metadata = result.get('cai_metadata', {})
                final_cai = cai_metadata.get('final_cai', None) if cai_metadata else None
                

                success = True
                error_msg = None
                
                if discrete_cai is None:
                    success = False


                    success = False

                elif final_cai is not None and abs(discrete_cai - final_cai) > 0.1:
                    success = False

                

                test_result = TestResult(
                    test_name=f"Lagrangian_CAI_Recovery_{protein_name}",
                    success=success,
                    cai_value=discrete_cai,
                    error_message=error_msg,
                    metadata={
                        'final_cai': final_cai,
                        'original_cai': cai_metadata.get('original_cai') if cai_metadata else None,
                        'protein_sequence': protein_seq,
                        'execution_time': end_time - start_time
                    }
                )
                
                test_results.append(test_result)
                self.log_test_result(test_result)
                

                if success:
                    print(f"  ✅ discrete_cai_value: {discrete_cai:.6f}")
                    print(f"  ✅ final_cai: {final_cai:.6f if final_cai else 'N/A'}")

                else:
                    print(f"  ❌ {error_msg}")
                    
            except Exception as e:
                error_msg = str(e)
                test_result = TestResult(
                    test_name=f"Lagrangian_CAI_Recovery_{protein_name}",
                    success=False,
                    cai_value=None,
                    error_message=error_msg
                )
                test_results.append(test_result)
                self.log_test_result(test_result)

                
        return test_results
        
    def test_dynamic_lambda_cai_adjustment(self) -> List[TestResult]:
        """测试2: 动态lambda_cai调整验证"""
        self.log_test_start("动态lambda_cai调整验证")
        
        test_results = []
        protein_seq = self.test_proteins['medium']
        
        print(f"测试蛋白质: {protein_seq}")
        

        test_configs = [
            {
                'name': 'adaptive_enabled',
                'adaptive_lambda_cai': True,
                'lambda_cai_min': 0.01,
                'lambda_cai_max': 1.0,
                'lambda_cai_adjustment_factor': 1.5
            },
            {
                'name': 'adaptive_disabled',
                'adaptive_lambda_cai': False,
                'cai_lambda': 0.1
            }
        ]
        
        for config in test_configs:
            print(f"\n--- 配置: {config['name']} ---")
            
            try:

                constraint_kwargs = {
                    'enable_cai': True,
                    'cai_target': 0.8,
                    'verbose': False
                }
                constraint_kwargs.update({k: v for k, v in config.items() if k != 'name'})
                
                constraint = LagrangianConstraint(protein_seq, **constraint_kwargs)
                

                lambda_history = []
                cai_history = []
                
                start_time = time.time()
                
                for iteration in range(5):
                    result = constraint.forward(alpha=1.0, tau=1.0)
                    
                    current_lambda = getattr(constraint, 'cai_lambda', None)
                    discrete_cai = result.get('discrete_cai_value', 0)
                    
                    lambda_history.append(current_lambda)
                    cai_history.append(discrete_cai)
                    
                    if iteration < 4:
                        if hasattr(constraint, '_update_lambda_cai'):
                            constraint._update_lambda_cai(discrete_cai)
                
                end_time = time.time()
                

                final_lambda = lambda_history[-1]
                final_cai = cai_history[-1]
                
                success = True
                error_msg = None
                
                if config['name'] == 'adaptive_enabled':

                    lambda_variance = np.var(lambda_history) if len(lambda_history) > 1 else 0
                    if lambda_variance < 1e-10:
                        success = False
                        error_msg = "自适应模式下lambda_cai未发生变化"
                elif config['name'] == 'adaptive_disabled':

                    lambda_variance = np.var(lambda_history) if len(lambda_history) > 1 else 0
                    if lambda_variance > 1e-6:
                        success = False
                        error_msg = "固定模式下lambda_cai发生了意外变化"
                

                test_result = TestResult(
                    test_name=f"Dynamic_Lambda_CAI_{config['name']}",
                    success=success,
                    cai_value=final_cai,
                    error_message=error_msg,
                    metadata={
                        'lambda_history': lambda_history,
                        'cai_history': cai_history,
                        'lambda_variance': lambda_variance,
                        'final_lambda': final_lambda,
                        'execution_time': end_time - start_time,
                        'config': config
                    }
                )
                
                test_results.append(test_result)
                self.log_test_result(test_result)
                

                print(f"  Lambda历史: {[f'{l:.4f}' if l else 'None' for l in lambda_history]}")
                print(f"  CAI历史: {[f'{c:.6f}' for c in cai_history]}")
                print(f"  Lambda方差: {lambda_variance:.2e}")
                print(f"  最终结果: {'✅' if success else '❌'}")
                if error_msg:
                    print(f"  错误: {error_msg}")
                    
            except Exception as e:
                error_msg = str(e)
                test_result = TestResult(
                    test_name=f"Dynamic_Lambda_CAI_{config['name']}",
                    success=False,
                    cai_value=None,
                    error_message=error_msg
                )
                test_results.append(test_result)
                self.log_test_result(test_result)
                print(f"  ❌ 异常: {error_msg}")
        
        return test_results
        
    def test_constraint_consistency(self) -> List[TestResult]:


        
        test_results = []

        

        

        constraint_configs = {
            'lagrangian': {
                'class': LagrangianConstraint,
                'kwargs': {
                    'enable_cai': True,
                    'cai_target': 0.8,
                    'cai_lambda': 0.1
                }
            },
            'amino_matching': {
                'class': AminoMatchingSoftmax,
                'kwargs': {
                    'enable_cai': True,
                    'cai_target': 0.8,
                    'cai_weight': 0.1  # Uses cai_weight instead of cai_lambda
                }
            },
            'codon_profile': {
                'class': CodonProfileConstraint,
                'kwargs': {
                    'enable_cai': True,
                    'cai_target': 0.8,
                    'cai_weight': 0.1  # Uses cai_weight instead of cai_lambda
                }
            }
        }
        
        cai_results = {}
        
        for constraint_name, config in constraint_configs.items():

            
            try:

                constraint = config['class'](protein_seq, **config['kwargs'])
                

                start_time = time.time()
                result = constraint.forward(alpha=1.0, tau=1.0)
                end_time = time.time()
                

                discrete_cai = result.get('discrete_cai_value', None)
                cai_metadata = result.get('cai_metadata', {})
                final_cai = cai_metadata.get('final_cai', None) if cai_metadata else None
                
                cai_results[constraint_name] = {
                    'discrete_cai': discrete_cai,
                    'final_cai': final_cai,
                    'execution_time': end_time - start_time
                }
                

                success = discrete_cai is not None and discrete_cai > 1e-6

                
                test_result = TestResult(
                    test_name=f"Constraint_Consistency_{constraint_name}",
                    success=success,
                    cai_value=discrete_cai,
                    error_message=error_msg,
                    metadata={
                        'constraint_type': constraint_name,
                        'final_cai': final_cai,
                        'execution_time': end_time - start_time
                    }
                )
                
                test_results.append(test_result)
                self.log_test_result(test_result)
                
                print(f"  discrete_cai: {discrete_cai:.6f if discrete_cai else 'None'}")
                print(f"  final_cai: {final_cai:.6f if final_cai else 'None'}")

                
            except Exception as e:
                error_msg = str(e)
                test_result = TestResult(
                    test_name=f"Constraint_Consistency_{constraint_name}",
                    success=False,
                    cai_value=None,
                    error_message=error_msg
                )
                test_results.append(test_result)
                self.log_test_result(test_result)

                cai_results[constraint_name] = {'error': error_msg}
        


        valid_cai_values = []
        for name, data in cai_results.items():
            if 'discrete_cai' in data and data['discrete_cai'] is not None:
                valid_cai_values.append((name, data['discrete_cai']))
        
        if len(valid_cai_values) >= 2:

            cai_vals = [val for _, val in valid_cai_values]
            cai_std = np.std(cai_vals)
            cai_mean = np.mean(cai_vals)
            

            



            

            consistency_result = TestResult(
                test_name="Constraint_CAI_Consistency",
                success=consistency_ok,
                cai_value=cai_mean,

                metadata={
                    'cai_std': cai_std,
                    'cai_mean': cai_mean,
                    'individual_results': dict(valid_cai_values)
                }
            )
            test_results.append(consistency_result)
            self.log_test_result(consistency_result)
        else:

        
        return test_results
        
    def test_regression_compatibility(self) -> List[TestResult]:
        """测试4: 回归测试和向后兼容性"""
        self.log_test_start("回归测试和向后兼容性")
        
        test_results = []
        protein_seq = self.test_proteins['medium']
        

        regression_tests = [
            {
                'name': 'default_behavior_no_cai',
                'kwargs': {'enable_cai': False},
                'expected': 'no_cai_values'
            },
            {
                'name': 'default_behavior_with_cai',
                'kwargs': {'enable_cai': True, 'cai_target': 0.8},
                'expected': 'valid_cai_values'
            },
            {
                'name': 'legacy_parameters',
                'kwargs': {
                    'enable_cai': True,
                    'cai_target': 0.7,
                    'cai_lambda': 0.2,
                    'adaptive_lambda_cai': False
                },
                'expected': 'valid_cai_values'
            }
        ]
        
        for test_config in regression_tests:
            print(f"\n--- 回归测试: {test_config['name']} ---")
            
            try:

                constraint = LagrangianConstraint(protein_seq, **test_config['kwargs'])
                
                start_time = time.time()
                result = constraint.forward(alpha=0.5, tau=1.0)
                end_time = time.time()
                

                success = True
                error_msg = None
                
                discrete_cai = result.get('discrete_cai_value', None)
                
                if test_config['expected'] == 'no_cai_values':

                    if 'discrete_cai_value' in result or 'cai_metadata' in result:
                        success = False
                        error_msg = "禁用CAI时仍返回CAI相关值"
                elif test_config['expected'] == 'valid_cai_values':

                    if discrete_cai is None or discrete_cai < 1e-6:
                        success = False
                        error_msg = f"期望有效CAI值，实际: {discrete_cai}"
                
                test_result = TestResult(
                    test_name=f"Regression_{test_config['name']}",
                    success=success,
                    cai_value=discrete_cai,
                    error_message=error_msg,
                    metadata={
                        'test_config': test_config,
                        'execution_time': end_time - start_time,
                        'result_keys': list(result.keys())
                    }
                )
                
                test_results.append(test_result)
                self.log_test_result(test_result)
                
                print(f"  结果键: {list(result.keys())}")
                print(f"  discrete_cai: {discrete_cai}")
                print(f"  执行时间: {end_time - start_time:.3f}秒")
                print(f"  状态: {'✅' if success else '❌'}")
                if error_msg:
                    print(f"  错误: {error_msg}")
                
            except Exception as e:
                error_msg = str(e)
                test_result = TestResult(
                    test_name=f"Regression_{test_config['name']}",
                    success=False,
                    cai_value=None,
                    error_message=error_msg
                )
                test_results.append(test_result)
                self.log_test_result(test_result)
                print(f"  ❌ 异常: {error_msg}")
        
        return test_results
        
    def test_performance_benchmarks(self) -> List[TestResult]:


        
        test_results = []
        

        perf_tests = [
            {
                'name': 'small_protein_performance',
                'protein': self.test_proteins['short'],
                'iterations': 10,

            },
            {
                'name': 'medium_protein_performance', 
                'protein': self.test_proteins['medium'],
                'iterations': 5,

            }
        ]
        
        for perf_test in perf_tests:

            
            try:
                protein_seq = perf_test['protein']
                iterations = perf_test['iterations']
                

                constraint = LagrangianConstraint(
                    protein_seq,
                    enable_cai=True,
                    cai_target=0.8,
                    cai_lambda=0.1,
                    verbose=False
                )
                

                execution_times = []
                cai_values = []
                
                total_start = time.time()
                
                for i in range(iterations):
                    iter_start = time.time()
                    result = constraint.forward(alpha=1.0, tau=1.0)
                    iter_end = time.time()
                    
                    execution_times.append(iter_end - iter_start)
                    cai_values.append(result.get('discrete_cai_value', 0))
                
                total_end = time.time()
                

                avg_time = np.mean(execution_times)
                std_time = np.std(execution_times)
                min_time = np.min(execution_times)
                max_time = np.max(execution_times)
                
                avg_cai = np.mean(cai_values)
                std_cai = np.std(cai_values)
                

                success = True
                error_msg = None
                
                if avg_time > perf_test['expected_time_per_iter']:
                    success = False


                    success = False


                    success = False

                
                test_result = TestResult(
                    test_name=f"Performance_{perf_test['name']}",
                    success=success,
                    cai_value=avg_cai,
                    error_message=error_msg,
                    performance_metrics={
                        'avg_execution_time': avg_time,
                        'std_execution_time': std_time,
                        'min_execution_time': min_time,
                        'max_execution_time': max_time,
                        'avg_cai': avg_cai,
                        'std_cai': std_cai,
                        'total_time': total_end - total_start,
                        'iterations': iterations
                    }
                )
                
                test_results.append(test_result)
                self.log_test_result(test_result)
                






                if error_msg:

                
            except Exception as e:
                error_msg = str(e)
                test_result = TestResult(
                    test_name=f"Performance_{perf_test['name']}",
                    success=False,
                    cai_value=None,
                    error_message=error_msg
                )
                test_results.append(test_result)
                self.log_test_result(test_result)

        
        return test_results
        
    def run_comprehensive_validation(self) -> Dict[str, Any]:
        """运行全面验证测试"""
        print("🚀 开始执行CAI修复全面验证测试")
        print(f"{'='*100}")
        
        start_time = time.time()
        

        test_suites = [
            ("Lagrangian CAI恢复验证", self.test_lagrangian_cai_recovery),
            ("动态lambda_cai调整验证", self.test_dynamic_lambda_cai_adjustment),
            ("约束一致性验证", self.test_constraint_consistency),
            ("回归测试和向后兼容性", self.test_regression_compatibility),
            ("性能基准和收敛分析", self.test_performance_benchmarks)
        ]
        
        all_results = []
        for suite_name, test_func in test_suites:
            try:
                suite_results = test_func()
                all_results.extend(suite_results)
            except Exception as e:
                logger.error(f"测试套件 {suite_name} 执行失败: {str(e)}")
                error_result = TestResult(
                    test_name=f"{suite_name}_SUITE_ERROR",
                    success=False,
                    cai_value=None,
                    error_message=str(e)
                )
                all_results.append(error_result)
        
        end_time = time.time()
        

        report = self.generate_validation_report(all_results, end_time - start_time)
        
        return report
        
    def generate_validation_report(self, results: List[TestResult], total_time: float) -> Dict[str, Any]:

        print(f"\n{'='*100}")

        print(f"{'='*100}")
        

        total_tests = len(results)
        successful_tests = sum(1 for r in results if r.success)
        failed_tests = total_tests - successful_tests
        success_rate = successful_tests / total_tests * 100 if total_tests > 0 else 0
        

        category_stats = {}
        for result in results:
            category = result.test_name.split('_')[0]
            if category not in category_stats:
                category_stats[category] = {'total': 0, 'success': 0}
            category_stats[category]['total'] += 1
            if result.success:
                category_stats[category]['success'] += 1
        

        valid_cai_values = [r.cai_value for r in results if r.cai_value is not None and r.success]
        







        

        for category, stats in category_stats.items():
            cat_success_rate = stats['success'] / stats['total'] * 100
            print(f"  {category:20}: {stats['success']}/{stats['total']} ({cat_success_rate:.1f}%)")
        
        if valid_cai_values:





        

        failed_results = [r for r in results if not r.success]
        if failed_results:
            for result in failed_results:
                print(f"  • {result.test_name}: {result.error_message}")
        else:

        


        

        lagrangian_results = [r for r in results if 'Lagrangian' in r.test_name and 'Recovery' in r.test_name]
        lagrangian_success_rate = sum(1 for r in lagrangian_results if r.success) / len(lagrangian_results) * 100 if lagrangian_results else 0
        
        if lagrangian_success_rate >= 80:

        else:

        

        dynamic_results = [r for r in results if 'Dynamic' in r.test_name]
        dynamic_success = all(r.success for r in dynamic_results)
        
        if dynamic_success and dynamic_results:

        elif dynamic_results:

        else:

        

        consistency_results = [r for r in results if 'Consistency' in r.test_name]
        consistency_ok = any(r.success and 'Constraint_CAI_Consistency' == r.test_name for r in consistency_results)
        
        if consistency_ok:

        else:

        

        regression_results = [r for r in results if 'Regression' in r.test_name]
        regression_success_rate = sum(1 for r in regression_results if r.success) / len(regression_results) * 100 if regression_results else 0
        
        if regression_success_rate >= 90:

        else:

        


        if success_rate >= 90:


        elif success_rate >= 75:


        elif success_rate >= 50:


        else:


        

        report_data = {
            'summary': {
                'total_tests': total_tests,
                'successful_tests': successful_tests,
                'failed_tests': failed_tests,
                'success_rate': success_rate,
                'total_time': total_time,
                'health_status': health_status
            },
            'category_stats': category_stats,
            'cai_analysis': {
                'valid_values_count': len(valid_cai_values),
                'mean': np.mean(valid_cai_values) if valid_cai_values else None,
                'std': np.std(valid_cai_values) if valid_cai_values else None,
                'min': np.min(valid_cai_values) if valid_cai_values else None,
                'max': np.max(valid_cai_values) if valid_cai_values else None,
            },
            'key_findings': {
                'lagrangian_cai_fixed': lagrangian_success_rate >= 80,
                'dynamic_lambda_working': dynamic_success,
                'constraint_consistency': consistency_ok,
                'backward_compatibility': regression_success_rate >= 90
            },
            'detailed_results': [
                {
                    'test_name': r.test_name,
                    'success': r.success,
                    'cai_value': r.cai_value,
                    'error_message': r.error_message,
                    'metadata': r.metadata,
                    'performance_metrics': r.performance_metrics
                }
                for r in results
            ]
        }
        
        return report_data


def main():
    """主函数"""

    validator = CAIValidationFramework()
    

    report = validator.run_comprehensive_validation()
    

    report_file = Path('comprehensive_cai_validation_report.json')
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2, default=str)
    
    logger.info(f"验证报告已保存到: {report_file}")
    

    key_findings = report['key_findings']
    all_key_tests_pass = all([
        key_findings['lagrangian_cai_fixed'],
        key_findings['dynamic_lambda_working'],
        key_findings['constraint_consistency'],
        key_findings['backward_compatibility']
    ])
    
    if all_key_tests_pass:
        print(f"\n🎉 所有关键修复已验证成功！")
        return 0
    else:
        print(f"\n⚠️  部分关键功能仍存在问题，请检查详细报告")
        return 1


if __name__ == "__main__":
    exit(main())
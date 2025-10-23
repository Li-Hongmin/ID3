#!/usr/bin/env python3
"""




"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any
import re


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

class DataFlowTracer:

    
    def __init__(self):
        self.traced_values = []
        self.inconsistencies = []
        
    def trace_accessibility_flow(self, experiment_file: Path) -> Dict:
        """跟踪单个实验文件中accessibility值的完整流程"""
        try:
            with open(experiment_file, 'r') as f:
                data = json.load(f)
            

            accessibility_values = self._extract_all_accessibility_values(data)
            

            consistency_analysis = self._analyze_consistency(accessibility_values)
            

            trajectory_analysis = self._analyze_trajectory(data.get('trajectory', {}))
            

            best_final_analysis = self._analyze_best_final_relationship(data)
            

            mode_analysis = self._identify_computation_mode(data)
            
            result = {
                'file': experiment_file.name,
                'method': self._extract_method(experiment_file.name),
                'variant': self._extract_variant(experiment_file.name),
                'protein': self._extract_protein(experiment_file.name),
                
                'accessibility_values': accessibility_values,
                'consistency_analysis': consistency_analysis,
                'trajectory_analysis': trajectory_analysis,
                'best_final_analysis': best_final_analysis,
                'mode_analysis': mode_analysis,
                
                'data_flow_issues': []
            }
            

            issues = []
            if consistency_analysis.get('has_inconsistencies', False):
                issues.extend(consistency_analysis.get('inconsistencies', []))
            
            if trajectory_analysis.get('has_anomalies', False):
                issues.extend(trajectory_analysis.get('anomalies', []))
                
            if best_final_analysis.get('suspicious_relationship', False):
                issues.extend(best_final_analysis.get('issues', []))
            
            result['data_flow_issues'] = issues
            result['total_issues'] = len(issues)
            
            return result
            
        except Exception as e:
            return {
                'file': experiment_file.name,
                'error': str(e),
                'tracing_failed': True
            }
    
    def _extract_all_accessibility_values(self, data: Dict) -> Dict:

        values = {

            'final_accessibility': data.get('final_accessibility'),
            'best_accessibility': data.get('best_accessibility'),
            'initial_accessibility': data.get('initial_accessibility'),
            

            'trajectory_accessibility': data.get('trajectory', {}).get('accessibility', []),
            

            'best_seq_design_accessibility': data.get('best_seq_design', {}).get('accessibility'),
            

            'improvement': data.get('improvement'),
            

            'trajectory_count': len(data.get('trajectory', {}).get('accessibility', [])),
            'trajectory_min': None,
            'trajectory_max': None,
            'trajectory_final': None
        }
        

        traj_acc = values['trajectory_accessibility']
        if traj_acc:
            values['trajectory_min'] = min(traj_acc)
            values['trajectory_max'] = max(traj_acc)
            values['trajectory_final'] = traj_acc[-1] if traj_acc else None
        
        return values
    
    def _analyze_consistency(self, values: Dict) -> Dict:
        """分析accessibility值之间的一致性"""
        analysis = {
            'has_inconsistencies': False,
            'inconsistencies': [],
            'consistency_checks': {}
        }
        

        if (values['final_accessibility'] is not None and 
            values['trajectory_final'] is not None):
            diff = abs(values['final_accessibility'] - values['trajectory_final'])
            analysis['consistency_checks']['final_vs_trajectory_final'] = {
                'difference': diff,
                'consistent': diff < 1e-5
            }
            
            if diff >= 1e-5:
                analysis['has_inconsistencies'] = True
                analysis['inconsistencies'].append(
                    f"Final accessibility不匹配轨迹最终值: {diff:.6f}"
                )
        

        if (values['best_accessibility'] is not None and 
            values['trajectory_min'] is not None):

            diff = abs(values['best_accessibility'] - values['trajectory_min'])
            analysis['consistency_checks']['best_vs_trajectory_min'] = {
                'difference': diff,
                'best_value': values['best_accessibility'],
                'trajectory_min': values['trajectory_min'],
                'consistent': values['best_accessibility'] <= values['trajectory_min'] + 1e-5
            }
            
            if values['best_accessibility'] > values['trajectory_min'] + 1e-5:
                analysis['has_inconsistencies'] = True
                analysis['inconsistencies'].append(
                    f"Best accessibility大于轨迹最小值: best={values['best_accessibility']:.6f}, "
                    f"traj_min={values['trajectory_min']:.6f}"
                )
        

        if (values['best_seq_design_accessibility'] is not None and 
            values['best_accessibility'] is not None):
            diff = abs(values['best_seq_design_accessibility'] - values['best_accessibility'])
            analysis['consistency_checks']['best_vs_best_seq_design'] = {
                'difference': diff,
                'consistent': diff < 1e-5
            }
            
            if diff >= 1e-5:
                analysis['has_inconsistencies'] = True
                analysis['inconsistencies'].append(
                    f"Best sequence design中的accessibility不匹配: {diff:.6f}"
                )
        

        if (values['improvement'] is not None and 
            values['initial_accessibility'] is not None and 
            values['final_accessibility'] is not None):
            expected_improvement = values['initial_accessibility'] - values['final_accessibility']
            diff = abs(values['improvement'] - expected_improvement)
            analysis['consistency_checks']['improvement_calculation'] = {
                'recorded_improvement': values['improvement'],
                'calculated_improvement': expected_improvement,
                'difference': diff,
                'consistent': diff < 1e-5
            }
            
            if diff >= 1e-5:
                analysis['has_inconsistencies'] = True
                analysis['inconsistencies'].append(
                    f"Improvement计算不一致: recorded={values['improvement']:.6f}, "
                    f"calculated={expected_improvement:.6f}"
                )
        
        return analysis
    
    def _analyze_trajectory(self, trajectory: Dict) -> Dict:

        analysis = {
            'has_anomalies': False,
            'anomalies': [],
            'trajectory_stats': {}
        }
        
        if not trajectory:
            return analysis
        
        accessibility = trajectory.get('accessibility', [])
        iterations = trajectory.get('iterations', [])
        
        if not accessibility:
            return analysis
        

        analysis['trajectory_stats'] = {
            'length': len(accessibility),
            'min_value': min(accessibility),
            'max_value': max(accessibility),
            'first_value': accessibility[0] if accessibility else None,
            'last_value': accessibility[-1] if accessibility else None,
            'decreasing_trend': accessibility[0] > accessibility[-1] if len(accessibility) >= 2 else None
        }
        

        for i in range(1, len(accessibility)):
            diff = abs(accessibility[i] - accessibility[i-1])

                analysis['has_anomalies'] = True
                analysis['anomalies'].append(

                )
        

        non_decreasing_count = 0
        for i in range(1, len(accessibility)):
            if accessibility[i] > accessibility[i-1]:
                non_decreasing_count += 1
        
        non_decreasing_rate = non_decreasing_count / (len(accessibility) - 1) if len(accessibility) > 1 else 0

            analysis['has_anomalies'] = True
            analysis['anomalies'].append(

            )
        

        if accessibility:
            min_val, max_val = min(accessibility), max(accessibility)
            if min_val < 0.01 or max_val > 20.0:
                analysis['has_anomalies'] = True
                analysis['anomalies'].append(

                )
        
        return analysis
    
    def _analyze_best_final_relationship(self, data: Dict) -> Dict:
        """分析best和final accessibility的关系"""
        analysis = {
            'suspicious_relationship': False,
            'issues': [],
            'relationship_stats': {}
        }
        
        best = data.get('best_accessibility')
        final = data.get('final_accessibility')
        
        if best is None or final is None:
            return analysis
        
        ratio = best / final if final != 0 else float('inf')
        difference = final - best
        
        analysis['relationship_stats'] = {
            'best_value': best,
            'final_value': final,
            'ratio': ratio,
            'difference': difference,
            'best_is_better': best < final
        }
        

        if best < final * 0.5:
            analysis['suspicious_relationship'] = True
            analysis['issues'].append(
                f"Best显著好于final: best={best:.6f}, final={final:.6f}, ratio={ratio:.3f}"
            )
        

        if best < 0.2 and final > 1.0:
            analysis['suspicious_relationship'] = True
            analysis['issues'].append(
                f"Best极端优秀但final一般: best={best:.6f}, final={final:.6f}"
            )
        

        improvement = data.get('improvement', 0)
        if improvement < 0 and best < final:
            analysis['suspicious_relationship'] = True
            analysis['issues'].append(
                f"负improvement但best<final: improvement={improvement:.6f}"
            )
        
        return analysis
    
    def _identify_computation_mode(self, data: Dict) -> Dict:

        analysis = {
            'identified_mode': 'unknown',
            'mode_indicators': {},
            'confidence': 0.0
        }
        
        indicators = {}
        

        opt_time = data.get('optimization_time')
        if opt_time is not None:
            if opt_time < 20:
                indicators['fast_optimization'] = True
            elif opt_time > 100:
                indicators['slow_optimization'] = True
        

        variant = self._extract_variant(data.get('trajectory', {}).get('iterations', [0])[0] if data.get('trajectory') else "")
        if variant:

            if variant in ['01', '11']:
                indicators['ste_mode_likely'] = True  # beta=1
            elif variant in ['00', '10']:
                indicators['soft_prob_mode_likely'] = True  # beta=0
        

        best = data.get('best_accessibility')
        final = data.get('final_accessibility')
        if best is not None and final is not None:
            if abs(best - final) < 1e-5:

            elif best < final * 0.8:

        

        trajectory = data.get('trajectory', {}).get('accessibility', [])
        if len(trajectory) > 10:

            diffs = [abs(trajectory[i] - trajectory[i-1]) for i in range(1, len(trajectory))]
            if diffs:
                smoothness = np.std(diffs)
                if smoothness < 0.01:
                    indicators['smooth_trajectory'] = True
                elif smoothness > 0.1:
                    indicators['noisy_trajectory'] = True
        
        analysis['mode_indicators'] = indicators
        

        ste_score = 0
        soft_score = 0
        
        if indicators.get('fast_optimization', False):
            ste_score += 0.3
        if indicators.get('slow_optimization', False):
            soft_score += 0.3
        if indicators.get('ste_mode_likely', False):
            ste_score += 0.4
        if indicators.get('soft_prob_mode_likely', False):
            soft_score += 0.4
        if indicators.get('values_identical', False):
            ste_score += 0.3
        if indicators.get('significant_difference', False):
            soft_score += 0.3
            
        if ste_score > soft_score and ste_score > 0.5:
            analysis['identified_mode'] = 'ste'
            analysis['confidence'] = min(ste_score, 1.0)
        elif soft_score > ste_score and soft_score > 0.5:
            analysis['identified_mode'] = 'soft_probability'
            analysis['confidence'] = min(soft_score, 1.0)
        
        return analysis
    
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
    
    def analyze_batch(self, batch_path: Path, max_files: int = 25) -> Dict:


        
        json_files = [f for f in batch_path.glob("*.json") 
                     if f.name not in ['summary.json', 'progress.json']]
        
        if not json_files:

            return {}
        

        if len(json_files) > max_files:
            json_files = json_files[:max_files]

        
        results = {
            'batch_name': batch_path.name,
            'analyzed_files': len(json_files),
            'tracing_results': [],
            'summary': {
                'total_issues': 0,
                'consistency_issues': 0,
                'trajectory_anomalies': 0,
                'suspicious_relationships': 0,
                'mode_distribution': {'ste': 0, 'soft_probability': 0, 'unknown': 0}
            }
        }
        
        for file_path in json_files:
            result = self.trace_accessibility_flow(file_path)
            results['tracing_results'].append(result)
            

            if 'error' not in result:
                results['summary']['total_issues'] += result.get('total_issues', 0)
                
                if result['consistency_analysis'].get('has_inconsistencies', False):
                    results['summary']['consistency_issues'] += 1
                
                if result['trajectory_analysis'].get('has_anomalies', False):
                    results['summary']['trajectory_anomalies'] += 1
                
                if result['best_final_analysis'].get('suspicious_relationship', False):
                    results['summary']['suspicious_relationships'] += 1
                

                identified_mode = result['mode_analysis'].get('identified_mode', 'unknown')
                results['summary']['mode_distribution'][identified_mode] += 1
        
        return results
    
    def print_analysis_report(self, results: Dict):
        """打印数据流跟踪报告"""
        print(f"\n📊 数据流跟踪报告: {results['batch_name']}")
        print("=" * 60)
        
        summary = results['summary']
        total = results['analyzed_files']
        
        print(f"分析文件数: {total}")
        print(f"发现问题总数: {summary['total_issues']}")
        
        print(f"\n🔍 问题分布:")
        issues = [
            ('数据一致性问题', summary['consistency_issues']),
            ('轨迹异常', summary['trajectory_anomalies']),  
            ('可疑best/final关系', summary['suspicious_relationships'])
        ]
        
        for issue_name, count in issues:
            if count > 0:
                rate = count / total * 100
                status = "🔴" if rate > 30 else "🟡" if rate > 0 else "🟢"
                print(f"  {issue_name}: {count}/{total} ({rate:.1f}%) {status}")
            else:
                print(f"  {issue_name}: 无问题 🟢")
        
        print(f"\n🎯 识别的计算模式:")
        mode_dist = summary['mode_distribution']
        for mode, count in mode_dist.items():
            if count > 0:
                rate = count / total * 100
                print(f"  {mode}: {count}/{total} ({rate:.1f}%)")
        

        problem_cases = [r for r in results['tracing_results'] 
                        if 'error' not in r and r.get('total_issues', 0) > 0]
        
        if problem_cases:
            print(f"\n❌ 典型问题案例:")
            for case in problem_cases[:3]:
                print(f"  {case['file']}: {case['method'].upper()} "
                      f"({case.get('total_issues', 0)}个问题)")
                

                for issue in case.get('data_flow_issues', [])[:2]:
                    print(f"    - {issue}")

def main():

    tracer = DataFlowTracer()
    

    test_batch = Path("results/20250910_133850_unified_access_experiments")
    
    if test_batch.exists():

        results = tracer.analyze_batch(test_batch, max_files=20)
        tracer.print_analysis_report(results)
        

        output_file = Path("data_flow_tracing_report.json")
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False, default=str)

        
    else:


if __name__ == "__main__":
    main()
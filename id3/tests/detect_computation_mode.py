#!/usr/bin/env python3
"""




"""

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any
import re
import numpy as np


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

class ComputationModeDetector:

    
    def __init__(self):
        self.mode_patterns = {
            'ste_mode': {
                'indicators': [
                    'values_nearly_identical',    # best ≈ final



                ],
                'weight': 0.0
            },
            'soft_probability_mode': {
                'indicators': [
                    'significant_best_final_diff',  # best << final



                ],
                'weight': 0.0
            },
            'deferred_validation_mode': {
                'indicators': [




                ],
                'weight': 0.0
            }
        }
    
    def detect_single_file(self, experiment_file: Path) -> Dict:
        """检测单个实验文件的计算模式"""
        try:
            with open(experiment_file, 'r') as f:
                data = json.load(f)
            

            basic_info = self._extract_basic_info(data, experiment_file.name)
            

            mode_evidence = self._collect_mode_evidence(data)
            

            detected_modes = self._determine_modes(mode_evidence)
            

            special_cases = self._detect_special_cases(data)
            

            deferred_validation_analysis = self._analyze_deferred_validation(data)
            
            result = {
                'file': experiment_file.name,
                'basic_info': basic_info,
                'mode_evidence': mode_evidence,
                'detected_modes': detected_modes,
                'special_cases': special_cases,
                'deferred_validation': deferred_validation_analysis,
                'confidence_scores': self._calculate_confidence_scores(mode_evidence),
                'issues': []
            }
            

            issues = []
            
            if detected_modes['primary_mode'] == 'ambiguous':
                issues.append("无法明确识别计算模式")
            
            if special_cases.get('suspicious_accessibility_pattern', False):
                issues.append("发现可疑的accessibility模式")
                
            if deferred_validation_analysis.get('likely_deferred', False):
                issues.append("可能使用了延迟验证模式")
                
            if mode_evidence.get('inconsistent_indicators', False):
                issues.append("模式指示器相互矛盾")
            
            result['issues'] = issues
            result['total_issues'] = len(issues)
            
            return result
            
        except Exception as e:
            return {
                'file': experiment_file.name,
                'error': str(e),
                'detection_failed': True
            }
    
    def _extract_basic_info(self, data: Dict, filename: str) -> Dict:

        return {
            'method': self._extract_method(filename),
            'variant': self._extract_variant(filename),
            'protein': self._extract_protein(filename),
            'optimization_time': data.get('optimization_time'),
            'iterations': data.get('iterations'),
            'final_accessibility': data.get('final_accessibility'),
            'best_accessibility': data.get('best_accessibility'),
            'improvement': data.get('improvement'),
            'amino_acids_correct': data.get('amino_acids_correct', 0.0)
        }
    
    def _collect_mode_evidence(self, data: Dict) -> Dict:
        """收集各种模式的证据"""
        evidence = {

            'ste_evidence': {},

            'soft_prob_evidence': {},

            'deferred_evidence': {},

            'inconsistent_indicators': False
        }
        
        basic_info = data
        best = data.get('best_accessibility')
        final = data.get('final_accessibility')
        opt_time = data.get('optimization_time')
        trajectory = data.get('trajectory', {})
        

        if best is not None and final is not None:
            diff = abs(best - final)
            if diff < 1e-5:
                evidence['ste_evidence']['values_identical'] = True
            elif diff < 0.001:
                evidence['ste_evidence']['values_nearly_identical'] = True
            else:
                evidence['ste_evidence']['significant_difference'] = diff
        
        if opt_time is not None:
            if opt_time < 30:
                evidence['ste_evidence']['fast_optimization'] = True
            elif opt_time > 100:
                evidence['soft_prob_evidence']['slow_optimization'] = True
        

        variant = self._extract_variant(data.get('file', ''))
        if variant in ['01', '11']:
            evidence['ste_evidence']['variant_suggests_ste'] = True
        elif variant in ['00', '10']:
            evidence['soft_prob_evidence']['variant_suggests_soft'] = True
        

        traj_acc = trajectory.get('accessibility', [])
        if traj_acc:

            if len(traj_acc) > 5:
                diffs = [abs(traj_acc[i] - traj_acc[i-1]) for i in range(1, len(traj_acc))]
                smoothness = np.std(diffs) if diffs else 0
                
                if smoothness < 0.01:
                    evidence['ste_evidence']['smooth_trajectory'] = True
                elif smoothness > 0.1:
                    evidence['soft_prob_evidence']['noisy_trajectory'] = True
            

            min_traj = min(traj_acc)
            if min_traj < 0.2:
                evidence['deferred_evidence']['suspicious_low_values'] = min_traj
            

            if traj_acc and final is not None:
                final_traj = traj_acc[-1]
                if abs(final_traj - final) > 1e-5:
                    evidence['deferred_evidence']['trajectory_final_mismatch'] = {
                        'trajectory_final': final_traj,
                        'recorded_final': final,
                        'difference': abs(final_traj - final)
                    }
        

        if best is not None and final is not None:
            if best < final * 0.8:
                evidence['soft_prob_evidence']['significant_best_final_diff'] = {
                    'best': best,
                    'final': final,
                    'ratio': best / final if final != 0 else float('inf')
                }
        


        iterations = trajectory.get('iterations', [])
        timestamps = trajectory.get('timestamps', [])
        
        if len(iterations) > 0 and len(timestamps) > 0:

            if len(timestamps) > 10:
                time_diffs = [timestamps[i] - timestamps[i-1] for i in range(1, len(timestamps))]

                large_gaps = [d for d in time_diffs if d > 5.0]
                if large_gaps:
                    evidence['deferred_evidence']['batch_processing_gaps'] = len(large_gaps)
        

        ste_count = len([k for k, v in evidence['ste_evidence'].items() if v])
        soft_count = len([k for k, v in evidence['soft_prob_evidence'].items() if v])
        
        if ste_count > 0 and soft_count > 0:
            evidence['inconsistent_indicators'] = True
        
        return evidence
    
    def _determine_modes(self, evidence: Dict) -> Dict:

        modes = {
            'primary_mode': 'unknown',
            'secondary_modes': [],
            'confidence': 0.0,
            'evidence_summary': {}
        }
        

        ste_score = len([k for k, v in evidence['ste_evidence'].items() if v and k != 'significant_difference'])
        soft_score = len([k for k, v in evidence['soft_prob_evidence'].items() if v])
        deferred_score = len([k for k, v in evidence['deferred_evidence'].items() if v])
        

        if 'significant_difference' in evidence['ste_evidence']:
            ste_score -= 1
        
        modes['evidence_summary'] = {
            'ste_score': ste_score,
            'soft_prob_score': soft_score,
            'deferred_score': deferred_score
        }
        

        max_score = max(ste_score, soft_score, deferred_score)
        
        if max_score == 0:
            modes['primary_mode'] = 'unknown'
            modes['confidence'] = 0.0
        elif max_score == ste_score and ste_score > 0:
            modes['primary_mode'] = 'ste'
            modes['confidence'] = min(ste_score / 4.0, 1.0)
        elif max_score == soft_score and soft_score > 0:
            modes['primary_mode'] = 'soft_probability'
            modes['confidence'] = min(soft_score / 4.0, 1.0)
        elif max_score == deferred_score and deferred_score > 0:
            modes['primary_mode'] = 'deferred_validation'
            modes['confidence'] = min(deferred_score / 4.0, 1.0)
        else:
            modes['primary_mode'] = 'ambiguous'
            modes['confidence'] = 0.5
        

        scores = [('ste', ste_score), ('soft_probability', soft_score), ('deferred_validation', deferred_score)]
        sorted_scores = sorted(scores, key=lambda x: x[1], reverse=True)
        
        for mode_name, score in sorted_scores[1:]:
            if score > 0:
                modes['secondary_modes'].append(mode_name)
        
        return modes
    
    def _detect_special_cases(self, data: Dict) -> Dict:
        """检测特殊情况和异常模式"""
        special = {
            'suspicious_accessibility_pattern': False,
            'extreme_values': False,
            'optimization_failure': False,
            'data_corruption': False,
            'special_notes': []
        }
        
        best = data.get('best_accessibility')
        final = data.get('final_accessibility')
        improvement = data.get('improvement', 0)
        

        if best is not None and best < 0.15:
            special['suspicious_accessibility_pattern'] = True
            special['special_notes'].append(f"Best accessibility异常低: {best:.6f}")
        
        if final is not None and final < 0.1:
            special['extreme_values'] = True
            special['special_notes'].append(f"Final accessibility异常低: {final:.6f}")
        

        if improvement < -1.0:
            special['optimization_failure'] = True
            special['special_notes'].append(f"负向优化: improvement={improvement:.3f}")
        

        amino_acids_correct = data.get('amino_acids_correct', 100.0)
        if amino_acids_correct < 100.0:
            special['data_corruption'] = True
            special['special_notes'].append(f"氨基酸约束违反: {amino_acids_correct:.1f}%")
        
        return special
    
    def _analyze_deferred_validation(self, data: Dict) -> Dict:

        analysis = {
            'likely_deferred': False,
            'evidence': {},
            'confidence': 0.0
        }
        
        trajectory = data.get('trajectory', {})
        accessibility = trajectory.get('accessibility', [])
        
        if not accessibility:
            return analysis
        
        evidence = {}
        

        min_traj = min(accessibility)
        final = data.get('final_accessibility')
        
        if min_traj < 0.3 and final and final > 1.0:
            evidence['low_trajectory_normal_final'] = {
                'min_trajectory': min_traj,
                'final': final,
                'ratio': final / min_traj
            }
        

        opt_time = data.get('optimization_time')
        iterations = data.get('iterations', 1000)
        
        if opt_time and iterations:
            iter_per_second = iterations / opt_time

                evidence['slow_iteration_rate'] = {
                    'iterations_per_second': iter_per_second,
                    'total_time': opt_time
                }
        

        if len(accessibility) != iterations:
            evidence['trajectory_length_mismatch'] = {
                'trajectory_length': len(accessibility),
                'expected_iterations': iterations
            }
        
        analysis['evidence'] = evidence
        

        confidence = 0.0
        if evidence.get('low_trajectory_normal_final'):
            confidence += 0.4
        if evidence.get('slow_iteration_rate'):
            confidence += 0.3
        if evidence.get('trajectory_length_mismatch'):
            confidence += 0.3
        
        analysis['confidence'] = min(confidence, 1.0)
        analysis['likely_deferred'] = confidence > 0.5
        
        return analysis
    
    def _calculate_confidence_scores(self, evidence: Dict) -> Dict:
        """计算各种模式判断的置信度"""
        scores = {
            'ste_confidence': 0.0,
            'soft_prob_confidence': 0.0,
            'deferred_confidence': 0.0
        }
        

        ste_indicators = evidence['ste_evidence']
        ste_score = 0.0
        
        if ste_indicators.get('values_identical'):
            ste_score += 0.4
        elif ste_indicators.get('values_nearly_identical'):
            ste_score += 0.3
        
        if ste_indicators.get('fast_optimization'):
            ste_score += 0.2
        if ste_indicators.get('variant_suggests_ste'):
            ste_score += 0.2
        if ste_indicators.get('smooth_trajectory'):
            ste_score += 0.2
        
        scores['ste_confidence'] = min(ste_score, 1.0)
        

        soft_indicators = evidence['soft_prob_evidence']
        soft_score = 0.0
        
        if soft_indicators.get('significant_best_final_diff'):
            soft_score += 0.4
        if soft_indicators.get('slow_optimization'):
            soft_score += 0.2
        if soft_indicators.get('variant_suggests_soft'):
            soft_score += 0.2
        if soft_indicators.get('noisy_trajectory'):
            soft_score += 0.2
        
        scores['soft_prob_confidence'] = min(soft_score, 1.0)
        

        deferred_indicators = evidence['deferred_evidence']
        deferred_score = 0.0
        
        if deferred_indicators.get('suspicious_low_values'):
            deferred_score += 0.3
        if deferred_indicators.get('trajectory_final_mismatch'):
            deferred_score += 0.3
        if deferred_indicators.get('batch_processing_gaps'):
            deferred_score += 0.4
        
        scores['deferred_confidence'] = min(deferred_score, 1.0)
        
        return scores
    
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
    
    def analyze_batch(self, batch_path: Path, max_files: int = 30) -> Dict:
        """分析批次中的计算模式"""
        print(f"🔍 检测计算模式: {batch_path.name}")
        
        json_files = [f for f in batch_path.glob("*.json") 
                     if f.name not in ['summary.json', 'progress.json']]
        
        if not json_files:
            print(f"❌ 未找到JSON文件在 {batch_path}")
            return {}
        

        if len(json_files) > max_files:
            json_files = json_files[:max_files]
            print(f"⚠️ 限制分析文件数量为 {max_files}")
        
        results = {
            'batch_name': batch_path.name,
            'analyzed_files': len(json_files),
            'detection_results': [],
            'summary': {
                'mode_distribution': {
                    'ste': 0,
                    'soft_probability': 0,
                    'deferred_validation': 0,
                    'ambiguous': 0,
                    'unknown': 0
                },
                'special_cases': {
                    'suspicious_patterns': 0,
                    'extreme_values': 0,
                    'optimization_failures': 0,
                    'data_corruption': 0
                },
                'high_confidence_detections': 0,
                'deferred_validation_likely': 0
            }
        }
        
        for file_path in json_files:
            result = self.detect_single_file(file_path)
            results['detection_results'].append(result)
            

            if 'error' not in result:
                primary_mode = result['detected_modes']['primary_mode']
                results['summary']['mode_distribution'][primary_mode] += 1
                
                confidence = result['detected_modes']['confidence']
                if confidence > 0.7:
                    results['summary']['high_confidence_detections'] += 1
                

                special = result['special_cases']
                if special['suspicious_accessibility_pattern']:
                    results['summary']['special_cases']['suspicious_patterns'] += 1
                if special['extreme_values']:
                    results['summary']['special_cases']['extreme_values'] += 1
                if special['optimization_failure']:
                    results['summary']['special_cases']['optimization_failures'] += 1
                if special['data_corruption']:
                    results['summary']['special_cases']['data_corruption'] += 1
                

                if result['deferred_validation']['likely_deferred']:
                    results['summary']['deferred_validation_likely'] += 1
        
        return results
    
    def print_analysis_report(self, results: Dict):


        print("=" * 60)
        
        summary = results['summary']
        total = results['analyzed_files']
        


        

        mode_dist = summary['mode_distribution']
        for mode, count in mode_dist.items():
            if count > 0:
                rate = count / total * 100
                emoji = "🟢" if mode in ['ste', 'soft_probability'] else "🟡" if mode == 'deferred_validation' else "🔴"
                print(f"  {mode}: {count}/{total} ({rate:.1f}%) {emoji}")
        

        special = summary['special_cases']
        for case_type, count in special.items():
            if count > 0:
                rate = count / total * 100
                status = "🔴" if rate > 20 else "🟡" if count > 0 else "🟢"
                print(f"  {case_type}: {count}/{total} ({rate:.1f}%) {status}")
        
        deferred_count = summary['deferred_validation_likely']
        if deferred_count > 0:
            deferred_rate = deferred_count / total * 100

        

        detection_results = results['detection_results']
        

        ste_cases = [r for r in detection_results if r.get('detected_modes', {}).get('primary_mode') == 'ste']
        if ste_cases:

            case = ste_cases[0]
            info = case['basic_info']

            print(f"    best={info['best_accessibility']:.6f}, final={info['final_accessibility']:.6f}")
        

        suspicious_cases = [r for r in detection_results 
                           if r.get('special_cases', {}).get('suspicious_accessibility_pattern', False)]
        if suspicious_cases:

            for case in suspicious_cases[:2]:
                info = case['basic_info']
                print(f"  {case['file']}: {info['method'].upper()}")
                print(f"    best={info['best_accessibility']:.6f}, final={info['final_accessibility']:.6f}")

def main():
    """主函数"""
    detector = ComputationModeDetector()
    

    test_batch = Path("results/20250910_133850_unified_access_experiments")
    
    if test_batch.exists():
        print("🎯 验证可能性3&5: 计算模式和延迟验证问题")
        results = detector.analyze_batch(test_batch, max_files=25)
        detector.print_analysis_report(results)
        

        output_file = Path("computation_mode_detection_report.json")
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, ensure_ascii=False, default=str)
        print(f"\n💾 详细报告已保存: {output_file}")
        
    else:
        print(f"❌ 批次目录不存在: {test_batch}")

if __name__ == "__main__":
    main()
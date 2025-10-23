#!/usr/bin/env python3
"""


"""

import json
import re
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging

logger = logging.getLogger(__name__)

class FastExperimentAnalyzer:

    
    @staticmethod
    def extract_key_fields(file_path: Path) -> Optional[Dict]:
        """
        快速提取JSON关键字段（不完全解析）
        只读取需要的顶层字段，避免加载大的trajectory数据
        """
        try:
            with open(file_path, 'r') as f:

                content = f.read()
                
                result = {'filename': file_path.stem}
                

                patterns = {
                    'protein_name': r'"protein_name":\s*"([^"]+)"',
                    'constraint_type': r'"constraint_type":\s*"([^"]+)"',
                    'variant': r'"variant":\s*"([^"]+)"',
                    'best_accessibility': r'"best_accessibility":\s*([\d.]+)',
                    'final_accessibility': r'"final_accessibility":\s*([\d.]+)',
                    'aa_match_rate': r'"aa_match_rate":\s*([\d.]+)',
                    'amino_acids_match': r'"amino_acids_match":\s*(true|false)',
                    'amino_acids_correct': r'"amino_acids_correct":\s*([\d.]+)',
                    'iterations': r'"iterations":\s*(\d+)',
                    'final_ecai': r'"final_ecai":\s*([\d.]+)',
                    'cai_target_achieved': r'"cai_target_achieved":\s*(true|false)',
                }
                
                for key, pattern in patterns.items():
                    match = re.search(pattern, content)
                    if match:
                        value = match.group(1)

                        if key in ['best_accessibility', 'final_accessibility', 'aa_match_rate', 'final_ecai', 'amino_acids_correct']:
                            result[key] = float(value)
                        elif key == 'iterations':
                            result[key] = int(value)
                        elif key in ['cai_target_achieved', 'amino_acids_match']:
                            result[key] = value == 'true'
                        else:
                            result[key] = value
                


                try:
                    data = json.loads(content)
                    if 'best_seq_design' in data:
                        result['best_seq_design'] = data['best_seq_design']
                except:

                
                return result
                
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
            return None
    
    @staticmethod
    def analyze_directory(directory: Path, need_sequences: bool = False) -> Dict[str, Any]:
        """分析实验目录"""
        experiment_files = list(directory.glob("*.json"))
        experiment_files = [f for f in experiment_files 
                           if f.name not in ["config.json", "summary.json"]]
        
        logger.info(f"Fast analyzing {len(experiment_files)} files...")
        
        results = []
        for i, file_path in enumerate(experiment_files):
            if need_sequences:

                try:
                    with open(file_path, 'r') as f:
                        data = json.load(f)
                        if 'trajectory' in data and 'discrete_sequences' in data['trajectory']:
                            sequences = data['trajectory']['discrete_sequences']
                            unique_sequences = len(set(sequences))
                            data['uniqueness_rate'] = unique_sequences / len(sequences) if sequences else 0
                        results.append(data)
                except:
                    continue
            else:

                result = FastExperimentAnalyzer.extract_key_fields(file_path)
                if result:
                    results.append(result)
            

            if (i + 1) % 30 == 0:
                logger.info(f"  Processed {i + 1}/{len(experiment_files)} files...")
        

        analysis = FastExperimentAnalyzer._compute_statistics(results)
        analysis['raw_results'] = results
        
        return analysis
    
    @staticmethod
    def _compute_statistics(results: List[Dict]) -> Dict[str, Any]:

        if not results:
            return {'error': 'No results to analyze'}
        

        best_accs = [r.get('best_accessibility', float('inf')) for r in results 
                    if 'best_accessibility' in r]
        

        by_constraint = {}
        for r in results:
            if 'best_accessibility' in r:
                constraint = r.get('constraint_type', 'unknown')
                if constraint not in by_constraint:
                    by_constraint[constraint] = []
                by_constraint[constraint].append(r['best_accessibility'])
        

        has_cai = any('final_ecai' in r for r in results)
        

        constraint_satisfied = []
        for r in results:

            if r.get('amino_acids_match', False) or r.get('amino_acids_correct', 0) == 100.0:
                constraint_satisfied.append(r)
        

        valid_best = float('inf')
        best_result = None
        for r in constraint_satisfied:
            if r.get('best_accessibility', float('inf')) < valid_best:
                valid_best = r['best_accessibility']
                best_result = r
        
        return {
            'total_experiments': len(results),
            'has_cai': has_cai,
            'best_overall': min(best_accs) if best_accs else None,
            'valid_best': valid_best if valid_best < float('inf') else None,
            'best_result': best_result,
            'mean_best': np.mean(best_accs) if best_accs else None,
            'std_best': np.std(best_accs) if best_accs else None,
            'by_constraint': by_constraint,
            'constraint_satisfied_count': len(constraint_satisfied),
            'constraint_satisfaction_rate': len(constraint_satisfied) / len(results) if results else 0
        }
    
    @staticmethod
    def compare_experiments(dir1: Path, dir2: Path) -> Dict[str, Any]:
        """比较两个实验"""
        analysis1 = FastExperimentAnalyzer.analyze_directory(dir1)
        analysis2 = FastExperimentAnalyzer.analyze_directory(dir2)
        
        comparison = {
            'experiment1': {
                'path': str(dir1),
                'type': 'CAI' if analysis1.get('has_cai') else 'No CAI',
                'best': analysis1.get('best_overall'),
                'valid_best': analysis1.get('valid_best'),
                'mean': analysis1.get('mean_best'),
                'constraint_satisfaction': analysis1.get('constraint_satisfaction_rate', 0)
            },
            'experiment2': {
                'path': str(dir2),
                'type': 'CAI' if analysis2.get('has_cai') else 'No CAI',
                'best': analysis2.get('best_overall'),
                'valid_best': analysis2.get('valid_best'),
                'mean': analysis2.get('mean_best'),
                'constraint_satisfaction': analysis2.get('constraint_satisfaction_rate', 0)
            }
        }
        

        if analysis1.get('valid_best') and analysis2.get('valid_best'):
            if analysis1['valid_best'] < analysis2['valid_best']:
                comparison['winner'] = 'experiment1'
                comparison['improvement'] = analysis2['valid_best'] - analysis1['valid_best']
            else:
                comparison['winner'] = 'experiment2'
                comparison['improvement'] = analysis1['valid_best'] - analysis2['valid_best']
        
        return comparison
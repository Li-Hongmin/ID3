#!/usr/bin/env python3
"""





"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count
from functools import partial
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def process_single_experiment(args):
    """

    
    Args:
        args: (filepath, is_cai, cai_threshold, amino_sequence)
        
    Returns:

    """
    filepath, is_cai, cai_threshold, target_amino = args
    
    try:

        with open(filepath, 'r') as f:
            data = json.load(f)
        

        result = {
            'protein': data.get('protein_name'),
            'constraint': data.get('constraint_type'),
            'variant': data.get('variant'),
            'amino_sequence': data.get('amino_sequence', '')
        }
        
        if is_cai:

            result.update(process_cai_experiment(data, cai_threshold))
        else:

            result.update(process_access_experiment(data))
        

        if result.get('sequence'):
            result['constraint_satisfied'] = check_constraint_fast(
                result['sequence'], 
                result['amino_sequence']
            )
        else:
            result['constraint_satisfied'] = False
        
        return result
        
    except Exception as e:
        logger.warning(f"Error processing {filepath}: {e}")
        return None


def check_constraint_fast(rna_sequence: str, amino_sequence: str) -> bool:
    """

    """
    if not rna_sequence or not amino_sequence:
        return False
    
    try:

        rna_seq = Seq(rna_sequence.replace('T', 'U'))
        translated = str(rna_seq.translate())
        

        if translated.endswith('*'):
            translated = translated[:-1]
        
        return translated == amino_sequence
    except:
        return False


def process_cai_experiment(data: Dict, cai_threshold: float) -> Dict:

    result = {}
    trajectory = data.get('trajectory', {})
    

    cai_values = trajectory.get('discrete_cai_values', [])
    access_values = trajectory.get('accessibility', [])
    sequences = trajectory.get('discrete_sequences', [])
    
    if not cai_values or not access_values:

        best = data.get('best_seq_design', {})
        result['accessibility'] = best.get('accessibility', np.nan)
        result['cai'] = best.get('discrete_cai', 0)
        result['sequence'] = best.get('discrete_sequence', '')
        result['corrected'] = False
        return result
    

    best_acc = float('inf')
    best_cai = 0
    best_seq = ''
    best_idx = -1
    
    for i in range(min(len(cai_values), len(access_values))):
        if cai_values[i] >= cai_threshold and access_values[i] < best_acc:
            best_acc = access_values[i]
            best_cai = cai_values[i]
            best_idx = i
            if i < len(sequences):
                best_seq = sequences[i]
    
    if best_idx >= 0:
        result['accessibility'] = best_acc
        result['cai'] = best_cai
        result['sequence'] = best_seq
        result['corrected'] = True
    else:

        if cai_values:
            max_cai_idx = np.argmax(cai_values)
            result['accessibility'] = access_values[max_cai_idx] if max_cai_idx < len(access_values) else np.nan
            result['cai'] = cai_values[max_cai_idx]
            result['sequence'] = sequences[max_cai_idx] if max_cai_idx < len(sequences) else ''
            result['corrected'] = False
        else:
            result['accessibility'] = np.nan
            result['cai'] = 0
            result['sequence'] = ''
            result['corrected'] = False
    
    return result


def process_access_experiment(data: Dict) -> Dict:
    """处理Access-only实验数据"""
    result = {}
    best = data.get('best_seq_design', {})
    
    result['accessibility'] = best.get('accessibility', np.nan)
    result['cai'] = None
    result['sequence'] = best.get('discrete_sequence', '')
    result['corrected'] = False
    
    return result


class ParallelPaperPerformanceTableGenerator:

    

    PROTEINS = [
        'O15263', 'P00004', 'P01308', 'P01825',
        'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',
        'P31417', 'P42212', 'P61626', 'P99999'
    ]
    

    CONSTRAINT_NAMES = {
        'lagrangian': 'Lagrangian Multiplier',
        'ams': 'Amino Matching Softmax',
        'cpc': 'Codon Profile Constraint'
    }
    

    VARIANTS = ['00', '01', '10', '11']
    
    def __init__(self, n_processes: int = 40):
        """
        初始化表格生成器
        
        Args:
            n_processes: 并行进程数（默认40）
        """
        self.n_processes = min(n_processes, cpu_count())
        self.cai_threshold = 0.8
        logger.info(f"Using {self.n_processes} parallel processes")
    
    def load_experiments_parallel(self, directory: Path, is_cai: bool = False) -> List[Dict]:
        """
        并行加载和处理实验数据
        
        Args:
            directory: 实验目录
            is_cai: 是否为CAI实验
            
        Returns:
            处理后的实验数据列表
        """

        json_files = list(directory.glob('*.json'))
        experiment_files = [
            f for f in json_files 
            if not f.name.startswith(('config', 'progress', 'summary'))
        ]
        
        logger.info(f"Processing {len(experiment_files)} files from {directory.name}")
        

        process_args = [
            (f, is_cai, self.cai_threshold, '') 
            for f in experiment_files
        ]
        

        start_time = time.time()
        with Pool(self.n_processes) as pool:
            results = pool.map(process_single_experiment, process_args)
        

        valid_results = [r for r in results if r is not None]
        
        elapsed = time.time() - start_time
        logger.info(f"Processed {len(valid_results)} files in {elapsed:.2f} seconds")
        
        return valid_results
    
    def generate_12x12_table(self, access_dir: Path, cai_dir: Path) -> pd.DataFrame:
        """
        生成12×12性能表格
        
        Args:
            access_dir: Access-Only实验目录
            cai_dir: CAI实验目录
            
        Returns:
            性能表格DataFrame
        """

        logger.info("Loading Access-only experiments...")
        access_data = self.load_experiments_parallel(access_dir, is_cai=False)
        
        logger.info("Loading CAI experiments...")
        cai_data = self.load_experiments_parallel(cai_dir, is_cai=True)
        

        access_dict = {}
        for item in access_data:
            key = (item['protein'], item['constraint'], item['variant'])
            access_dict[key] = item
        
        cai_dict = {}
        for item in cai_data:
            key = (item['protein'], item['constraint'], item['variant'])
            cai_dict[key] = item
        

        table_data = []
        
        for protein in self.PROTEINS:
            row = {'Protein': protein}
            
            for constraint in ['lagrangian', 'ams', 'cpc']:
                for variant in self.VARIANTS:
                    column_name = f"{constraint}_{variant}"
                    key = (protein, constraint, variant)
                    

                    access_result = access_dict.get(key, {})
                    access_acc = access_result.get('accessibility', float('inf'))
                    access_satisfied = access_result.get('constraint_satisfied', False)
                    

                    cai_result = cai_dict.get(key, {})
                    cai_acc = cai_result.get('accessibility', float('inf'))
                    cai_val = cai_result.get('cai', 0)
                    cai_satisfied = cai_result.get('constraint_satisfied', False)
                    cai_corrected = cai_result.get('corrected', False)
                    

                    if access_satisfied and cai_satisfied:

                        if access_acc <= cai_acc:
                            row[column_name] = access_acc
                            row[f"{column_name}_source"] = 'access'
                        else:
                            row[column_name] = cai_acc
                            row[f"{column_name}_source"] = 'cai'
                            row[f"{column_name}_cai"] = cai_val
                    elif access_satisfied:
                        row[column_name] = access_acc
                        row[f"{column_name}_source"] = 'access'
                    elif cai_satisfied:
                        row[column_name] = cai_acc
                        row[f"{column_name}_source"] = 'cai'
                        row[f"{column_name}_cai"] = cai_val
                    else:

                        if access_acc <= cai_acc:
                            row[column_name] = access_acc
                            row[f"{column_name}_source"] = 'access'
                        else:
                            row[column_name] = cai_acc
                            row[f"{column_name}_source"] = 'cai'
                            row[f"{column_name}_cai"] = cai_val
            
            table_data.append(row)
        

        df = pd.DataFrame(table_data)
        main_columns = ['Protein']
        for constraint in ['lagrangian', 'ams', 'cpc']:
            for variant in self.VARIANTS:
                main_columns.append(f"{constraint}_{variant}")
        
        return df[main_columns]
    
    def generate_latex_table(self, df: pd.DataFrame) -> str:
        """
        生成LaTeX格式的表格（论文格式）
        """
        latex = []
        latex.append("\\begin{table*}[!t]")
        latex.append("\\centering")
        latex.append("\\caption{RNA accessibility optimization performance (kcal/mol) across ID3 variants\\label{tab:performance_matrix}}")
        latex.append("\\tiny")
        

        headers = ["\\textbf{Protein}"]
        for constraint in ['LAG', 'AMS', 'CPC']:
            for variant in self.VARIANTS:
                headers.append(f"\\rotatebox{{90}}{{\\textbf{{{constraint}-{variant}}}}}")
        headers.append("\\textbf{Avg}")
        

        latex.append("\\begin{tabular}{l|" + "c" * 12 + "|c}")
        latex.append("\\toprule")
        latex.append(" & ".join(headers) + " \\\\")
        latex.append("\\midrule")
        

        for _, row in df.iterrows():
            values = [row['Protein']]
            row_values = []
            
            for constraint in ['lagrangian', 'ams', 'cpc']:
                for variant in self.VARIANTS:
                    col_name = f"{constraint}_{variant}"
                    val = row[col_name]
                    if pd.isna(val) or val == float('inf'):
                        values.append("---")
                    else:
                        values.append(f"{val:.3f}")
                        row_values.append(val)
            

            if row_values:
                avg = np.mean(row_values)
                values.append(f"{avg:.3f}")
            else:
                values.append("---")
            
            latex.append(" & ".join(values) + " \\\\")
        

        latex.append("\\midrule")
        avg_values = ["\\textbf{Average}"]
        
        all_means = []
        for constraint in ['lagrangian', 'ams', 'cpc']:
            for variant in self.VARIANTS:
                col_name = f"{constraint}_{variant}"
                col_values = df[col_name].replace([float('inf')], np.nan).dropna()
                if len(col_values) > 0:
                    mean_val = col_values.mean()
                    avg_values.append(f"{mean_val:.3f}")
                    all_means.append(mean_val)
                else:
                    avg_values.append("---")
        

        if all_means:
            total_avg = np.mean(all_means)
            avg_values.append(f"\\textbf{{{total_avg:.3f}}}")
        else:
            avg_values.append("---")
        
        latex.append(" & ".join(avg_values) + " \\\\")
        
        latex.append("\\bottomrule")
        latex.append("\\end{tabular}")
        latex.append("\\end{table*}")
        
        return "\n".join(latex)
    
    def print_summary(self, df: pd.DataFrame):
        """打印表格摘要"""
        print("\n" + "=" * 80)
        print("📊 性能摘要")
        print("=" * 80)
        

        for constraint in ['lagrangian', 'ams', 'cpc']:
            values = []
            for variant in self.VARIANTS:
                col_name = f"{constraint}_{variant}"
                col_values = df[col_name].replace([float('inf')], np.nan).dropna()
                values.extend(col_values.tolist())
            
            if values:
                constraint_name = self.CONSTRAINT_NAMES[constraint]
                print(f"  {constraint_name}: {np.mean(values):.4f} kcal/mol (min: {np.min(values):.4f})")
        

        best_value = float('inf')
        best_config = None
        best_protein = None
        
        for _, row in df.iterrows():
            for constraint in ['lagrangian', 'ams', 'cpc']:
                for variant in self.VARIANTS:
                    col_name = f"{constraint}_{variant}"
                    val = row[col_name]
                    if val < best_value:
                        best_value = val
                        best_config = f"{self.CONSTRAINT_NAMES[constraint]} {variant}"
                        best_protein = row['Protein']
        
        if best_config:
            print(f"\n  最佳结果: {best_value:.4f} kcal/mol")
            print(f"  配置: {best_config}")
            print(f"  蛋白质: {best_protein}")


def main():

    import argparse
    

    parser.add_argument('--access-dir',
                       default='paper_experiment_results/access_only',

    parser.add_argument('--cai-dir',
                       default='paper_experiment_results/cai_with_penalty',

    parser.add_argument('--output-dir',
                       default='paper_experiment_results/tables',

    parser.add_argument('--processes', type=int, default=40,

    
    args = parser.parse_args()
    

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    

    generator = ParallelPaperPerformanceTableGenerator(n_processes=args.processes)
    


    start_time = time.time()
    
    df = generator.generate_12x12_table(Path(args.access_dir), Path(args.cai_dir))
    

    csv_path = output_dir / 'performance_12x12_parallel.csv'
    df.to_csv(csv_path, index=False)

    

    latex_content = generator.generate_latex_table(df)
    latex_path = output_dir / 'performance_12x12_parallel.tex'
    with open(latex_path, 'w') as f:
        f.write(latex_content)

    

    generator.print_summary(df)
    
    elapsed = time.time() - start_time



if __name__ == "__main__":
    main()
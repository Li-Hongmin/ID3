#!/usr/bin/env python3
"""
Paper Table Generator for ID3 Framework


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

from id3.experiments.analysis.fast_analyzer import FastExperimentAnalyzer

class PaperTableGenerator:

    

    PROTEINS = ['O15263', 'P00004', 'P01308', 'P01825', 'P04637', 'P0CG48', 
                'P0DTC2', 'P0DTC9', 'P31417', 'P42212', 'P61626', 'P99999']
    

    ID3_VARIANTS = {
        # Lagrangian (L)
        'lagrangian_00': 'ID3-L00', 'lagrangian_01': 'ID3-L01', 
        'lagrangian_10': 'ID3-L10', 'lagrangian_11': 'ID3-L11',
        # Amino Matching Softmax (A)  
        'ams_00': 'ID3-A00', 'ams_01': 'ID3-A01',
        'ams_10': 'ID3-A10', 'ams_11': 'ID3-A11',
        # Codon Profile Constraint (C)
        'cpc_00': 'ID3-C00', 'cpc_01': 'ID3-C01',
        'cpc_10': 'ID3-C10', 'cpc_11': 'ID3-C11'
    }
    

    SHORT_VARIANTS = {
        'lagrangian_00': 'L00', 'lagrangian_01': 'L01',
        'lagrangian_10': 'L10', 'lagrangian_11': 'L11',
        'ams_00': 'A00', 'ams_01': 'A01',
        'ams_10': 'A10', 'ams_11': 'A11',
        'cpc_00': 'C00', 'cpc_01': 'C01',
        'cpc_10': 'C10', 'cpc_11': 'C11'
    }
    
    def __init__(self, num_workers: int = 40):
        """
        初始化表格生成器
        
        Args:
            num_workers: 并行工作进程数，默认40
        """
        self.num_workers = num_workers
        self.analyzer = FastExperimentAnalyzer()
    
    def process_single_file(self, file_path: Path) -> Optional[Dict]:
        """处理单个JSON文件"""
        try:
            return self.analyzer.extract_key_fields(file_path)
        except Exception as e:
            print(f"⚠️ 处理文件失败 {file_path.name}: {e}")
            return None
    
    def collect_experiment_data_parallel(self, base_dir: Path, 
                                        experiment_type: str = "") -> pd.DataFrame:
        """

        
        Args:


            
        Returns:

        """
        print(f"📂 分析目录: {base_dir.name}")
        

        json_files = [f for f in base_dir.glob("*.json") 
                     if not f.name.startswith(('config', 'progress', 'summary'))]
        
        print(f"   找到 {len(json_files)} 个实验文件")
        

        results = []
        with ProcessPoolExecutor(max_workers=self.num_workers) as executor:

            future_to_file = {
                executor.submit(self.process_single_file, f): f 
                for f in json_files
            }
            

            with tqdm(total=len(json_files), desc="   处理文件", leave=False) as pbar:
                for future in as_completed(future_to_file):
                    data = future.result()
                    if data:
                        results.append(data)
                    pbar.update(1)
        
        print(f"   成功处理 {len(results)}/{len(json_files)} 个文件")
        

        if results:
            df = pd.DataFrame(results)
            return df
        else:
            return pd.DataFrame()
    
    def generate_accessibility_table(self, access_dir: Path, 
                                    output_dir: Path,
                                    use_short_names: bool = False) -> pd.DataFrame:
        """

        
        Args:



            
        Returns:

        """
        print("\n📊 生成RNA Accessibility性能表格...")
        

        df = self.collect_experiment_data_parallel(access_dir, 'access_only')
        
        if df.empty:
            print("❌ 没有找到有效数据")
            return pd.DataFrame()
        

        variants = self.SHORT_VARIANTS if use_short_names else self.ID3_VARIANTS
        

        matrix_data = []
        for protein in self.PROTEINS:
            row = {'Protein': protein}
            protein_df = df[df['protein_name'] == protein]
            
            for old_name, new_name in variants.items():
                constraint, variant = old_name.rsplit('_', 1)
                matching = protein_df[
                    (protein_df['constraint_type'] == constraint) & 
                    (protein_df['variant'] == variant)
                ]
                
                if not matching.empty:

                    row[new_name] = matching.iloc[0]['best_accessibility']
                else:
                    row[new_name] = np.nan
            
            matrix_data.append(row)
        
        result_df = pd.DataFrame(matrix_data)
        result_df.set_index('Protein', inplace=True)
        

        csv_path = output_dir / 'table1_accessibility_performance.csv'
        result_df.to_csv(csv_path)
        print(f"   ✅ 保存CSV: {csv_path}")
        

        self.generate_latex_table(result_df, 
                                 output_dir / 'table1_accessibility_performance.tex',
                                 caption="RNA Accessibility Optimization Performance (kcal/mol)",
                                 label="tab:accessibility_performance")
        

        self.print_statistics(result_df)
        
        return result_df
    
    def generate_cai_comparison_table(self, cai_penalty_dir: Path, 
                                     cai_no_penalty_dir: Path,
                                     output_dir: Path) -> pd.DataFrame:
        """

        
        Args:



            
        Returns:

        """
        print("\n📊 生成CAI对比表格...")
        

        df_penalty = self.collect_experiment_data_parallel(cai_penalty_dir, 'cai_penalty')
        df_no_penalty = self.collect_experiment_data_parallel(cai_no_penalty_dir, 'cai_no_penalty')
        

        deterministic_variants = ['L00', 'L01', 'A00', 'A01', 'C00', 'C01']
        

        comparison_data = []
        for protein in self.PROTEINS:
            for variant_key in ['lagrangian_00', 'lagrangian_01', 
                               'ams_00', 'ams_01', 'cpc_00', 'cpc_01']:
                constraint, variant = variant_key.rsplit('_', 1)
                short_name = self.SHORT_VARIANTS[variant_key]
                

                penalty_match = df_penalty[
                    (df_penalty['protein_name'] == protein) &
                    (df_penalty['constraint_type'] == constraint) &
                    (df_penalty['variant'] == variant)
                ]
                

                no_penalty_match = df_no_penalty[
                    (df_no_penalty['protein_name'] == protein) &
                    (df_no_penalty['constraint_type'] == constraint) &
                    (df_no_penalty['variant'] == variant)
                ]
                

                penalty_cai = np.nan
                no_penalty_cai = np.nan
                
                if not penalty_match.empty:
                    penalty_data = penalty_match.iloc[0]
                    if 'best_seq_design' in penalty_data and isinstance(penalty_data['best_seq_design'], dict):
                        penalty_cai = penalty_data['best_seq_design'].get('discrete_cai', np.nan)
                
                if not no_penalty_match.empty:
                    no_penalty_data = no_penalty_match.iloc[0]
                    if 'best_seq_design' in no_penalty_data and isinstance(no_penalty_data['best_seq_design'], dict):
                        no_penalty_cai = no_penalty_data['best_seq_design'].get('discrete_cai', np.nan)
                
                row = {
                    'Protein': protein,
                    'Variant': short_name,
                    'With_Penalty_Access': penalty_match.iloc[0]['best_accessibility'] if not penalty_match.empty else np.nan,
                    'No_Penalty_Access': no_penalty_match.iloc[0]['best_accessibility'] if not no_penalty_match.empty else np.nan,
                    'With_Penalty_CAI': penalty_cai,
                    'No_Penalty_CAI': no_penalty_cai,
                }
                comparison_data.append(row)
        
        df_comparison = pd.DataFrame(comparison_data)
        

        df_pivot = df_comparison.pivot_table(
            index='Protein',
            columns='Variant', 
            values=['With_Penalty_Access', 'No_Penalty_Access', 
                   'With_Penalty_CAI', 'No_Penalty_CAI']
        )
        

        csv_path = output_dir / 'table2_cai_comparison.csv'
        df_pivot.to_csv(csv_path)
        print(f"   ✅ 保存CSV: {csv_path}")
        

        self.generate_latex_table(df_pivot,
                                 output_dir / 'table2_cai_comparison.tex',
                                 caption="CAI Integration Performance Comparison",
                                 label="tab:cai_comparison")
        
        return df_comparison
    
    def generate_combined_tables_from_files(self, 
                                           cai_penalty_path: Path,
                                           cai_no_penalty_path: Path,
                                           output_dir: Path):
        """

        
        Args:



        """
        print("\n📊 生成组合表格（Access + CAI）...")
        

        proteins = ['O15263', 'P00004', 'P01308', 'P01825', 'P04637', 'P0CG48', 
                   'P0DTC2', 'P0DTC9', 'P31417', 'P42212', 'P61626', 'P99999']
        

        variant_map = {
            ('lagrangian', '00'): 'L00', ('lagrangian', '01'): 'L01',
            ('lagrangian', '10'): 'L10', ('lagrangian', '11'): 'L11',
            ('ams', '00'): 'A00', ('ams', '01'): 'A01',
            ('ams', '10'): 'A10', ('ams', '11'): 'A11',
            ('cpc', '00'): 'C00', ('cpc', '01'): 'C01',
            ('cpc', '10'): 'C10', ('cpc', '11'): 'C11',
        }
        

        from id3.experiments.analysis.fast_analyzer import FastExperimentAnalyzer
        

        penalty_results = FastExperimentAnalyzer.analyze_directory(cai_penalty_path, need_sequences=False)
        penalty_data = penalty_results.get('raw_results', [])
        

        no_penalty_results = FastExperimentAnalyzer.analyze_directory(cai_no_penalty_path, need_sequences=False)
        no_penalty_data = no_penalty_results.get('raw_results', [])
        

        combined_penalty_data = []
        for protein in proteins:
            row = {'Protein': protein}
            for (method, variant), short_name in variant_map.items():

                found = False
                for result in penalty_data:
                    if (result.get('protein_id') == protein and 
                        result.get('constraint_type') == method and
                        result.get('discretization_variant') == variant):
                        access_val = result.get('best_accessibility', np.nan)
                        cai_val = np.nan
                        if 'best_seq_design' in result and isinstance(result['best_seq_design'], dict):
                            cai_val = result['best_seq_design'].get('discrete_cai', np.nan)
                        row[f'{short_name}_Access'] = access_val
                        row[f'{short_name}_CAI'] = cai_val
                        found = True
                        break
                
                if not found:
                    row[f'{short_name}_Access'] = np.nan
                    row[f'{short_name}_CAI'] = np.nan
            
            combined_penalty_data.append(row)
        

        mean_row = {'Protein': 'Mean'}
        for (method, variant), short_name in variant_map.items():
            access_vals = [row[f'{short_name}_Access'] for row in combined_penalty_data if not np.isnan(row.get(f'{short_name}_Access', np.nan))]
            cai_vals = [row[f'{short_name}_CAI'] for row in combined_penalty_data if not np.isnan(row.get(f'{short_name}_CAI', np.nan))]
            mean_row[f'{short_name}_Access'] = np.mean(access_vals) if access_vals else np.nan
            mean_row[f'{short_name}_CAI'] = np.mean(cai_vals) if cai_vals else np.nan
        combined_penalty_data.append(mean_row)
        

        df_penalty = pd.DataFrame(combined_penalty_data)
        csv_path_penalty = output_dir / 'table2_combined_with_penalty.csv'
        df_penalty.to_csv(csv_path_penalty, index=False)
        print(f"   ✅ 保存CSV: {csv_path_penalty}")
        

        combined_no_penalty_data = []
        for protein in proteins:
            row = {'Protein': protein}
            for (method, variant), short_name in variant_map.items():

                found = False
                for result in no_penalty_data:
                    if (result.get('protein_id') == protein and 
                        result.get('constraint_type') == method and
                        result.get('discretization_variant') == variant):
                        access_val = result.get('best_accessibility', np.nan)
                        cai_val = np.nan
                        if 'best_seq_design' in result and isinstance(result['best_seq_design'], dict):
                            cai_val = result['best_seq_design'].get('discrete_cai', np.nan)
                        row[f'{short_name}_Access'] = access_val
                        row[f'{short_name}_CAI'] = cai_val
                        found = True
                        break
                
                if not found:
                    row[f'{short_name}_Access'] = np.nan
                    row[f'{short_name}_CAI'] = np.nan
            
            combined_no_penalty_data.append(row)
        

        mean_row = {'Protein': 'Mean'}
        for (method, variant), short_name in variant_map.items():
            access_vals = [row[f'{short_name}_Access'] for row in combined_no_penalty_data if not np.isnan(row.get(f'{short_name}_Access', np.nan))]
            cai_vals = [row[f'{short_name}_CAI'] for row in combined_no_penalty_data if not np.isnan(row.get(f'{short_name}_CAI', np.nan))]
            mean_row[f'{short_name}_Access'] = np.mean(access_vals) if access_vals else np.nan
            mean_row[f'{short_name}_CAI'] = np.mean(cai_vals) if cai_vals else np.nan
        combined_no_penalty_data.append(mean_row)
        

        df_no_penalty = pd.DataFrame(combined_no_penalty_data)
        csv_path_no_penalty = output_dir / 'table2_combined_no_penalty.csv'
        df_no_penalty.to_csv(csv_path_no_penalty, index=False)
        print(f"   ✅ 保存CSV: {csv_path_no_penalty}")
        
        return df_penalty, df_no_penalty
    
    def generate_summary_statistics(self, access_df: pd.DataFrame,
                                   output_dir: Path) -> pd.DataFrame:
        """

        
        Args:


            
        Returns:

        """
        print("\n📊 生成统计汇总表...")
        
        stats = []
        

        mechanisms = [
            ('Lagrangian', ['L00', 'L01', 'L10', 'L11']),
            ('Amino Matching Softmax', ['A00', 'A01', 'A10', 'A11']),
            ('Codon Profile Constraint', ['C00', 'C01', 'C10', 'C11'])
        ]
        
        for name, variants in mechanisms:
            values = access_df[variants].values.flatten()
            values = values[~np.isnan(values)]
            
            if len(values) > 0:
                stats.append({
                    'Category': name,
                    'Mean': np.mean(values),
                    'Std': np.std(values),
                    'Min': np.min(values),
                    'Max': np.max(values),
                    'Median': np.median(values),
                    'Count': len(values),
                    'Excellent (<1.0)': (values < 1.0).sum(),
                    'Good (<1.5)': (values < 1.5).sum()
                })
        

        strategies = [
            ('Deterministic Soft (00)', ['L00', 'A00', 'C00']),
            ('Deterministic Hard (01)', ['L01', 'A01', 'C01']),
            ('Stochastic Soft (10)', ['L10', 'A10', 'C10']),
            ('Stochastic Hard (11)', ['L11', 'A11', 'C11'])
        ]
        
        for name, variants in strategies:
            values = access_df[variants].values.flatten()
            values = values[~np.isnan(values)]
            
            if len(values) > 0:
                stats.append({
                    'Category': name,
                    'Mean': np.mean(values),
                    'Std': np.std(values),
                    'Min': np.min(values),
                    'Max': np.max(values),
                    'Median': np.median(values),
                    'Count': len(values),
                    'Excellent (<1.0)': (values < 1.0).sum(),
                    'Good (<1.5)': (values < 1.5).sum()
                })
        
        df_stats = pd.DataFrame(stats)
        

        csv_path = output_dir / 'table3_summary_statistics.csv'
        df_stats.to_csv(csv_path, index=False)
        print(f"   ✅ 保存CSV: {csv_path}")
        

        self.generate_latex_table(df_stats.set_index('Category'),
                                 output_dir / 'table3_summary_statistics.tex',
                                 caption="Summary Statistics by Constraint Mechanism and Strategy",
                                 label="tab:summary_statistics")
        
        return df_stats
    
    def generate_latex_table(self, df: pd.DataFrame, output_path: Path,
                            caption: str, label: str):

        latex = df.to_latex(
            float_format="%.4f",
            na_rep='--',
            caption=caption,
            label=label,
            column_format='l' + 'r' * len(df.columns)
        )
        

        latex = latex.replace('\\toprule', '\\hline')
        latex = latex.replace('\\midrule', '\\hline')
        latex = latex.replace('\\bottomrule', '\\hline')
        
        with open(output_path, 'w') as f:
            f.write(latex)

    
    def print_statistics(self, df: pd.DataFrame):
        """打印统计信息"""
        print("\n📈 统计汇总:")
        

        lagrangian_cols = [c for c in df.columns if 'L0' in c or 'L1' in c]
        ams_cols = [c for c in df.columns if 'A0' in c or 'A1' in c]
        cpc_cols = [c for c in df.columns if 'C0' in c or 'C1' in c]
        
        for name, cols in [('Lagrangian', lagrangian_cols),
                          ('AMS', ams_cols),
                          ('CPC', cpc_cols)]:
            if cols:
                values = df[cols].values.flatten()
                values = values[~np.isnan(values)]
                if len(values) > 0:
                    print(f"   {name}: {np.mean(values):.4f} ± {np.std(values):.4f} kcal/mol")
                    print(f"      最佳: {np.min(values):.4f}, 最差: {np.max(values):.4f}")
    
    def run(self, access_dir: str, cai_penalty_dir: str, cai_no_penalty_dir: str,
            output_dir: str = 'paper_experiment_results/tables'):
        """

        
        Args:




        """
        print("=" * 60)
        print("🔬 ID3论文表格生成器")
        print(f"   并行进程数: {self.num_workers}")
        print("=" * 60)
        

        access_path = Path(access_dir)
        cai_penalty_path = Path(cai_penalty_dir)
        cai_no_penalty_path = Path(cai_no_penalty_dir)
        output_path = Path(output_dir)
        

        output_path.mkdir(parents=True, exist_ok=True)
        

        access_df = self.generate_accessibility_table(access_path, output_path, use_short_names=True)
        
        if not access_df.empty:
            cai_df = self.generate_cai_comparison_table(cai_penalty_path, cai_no_penalty_path, output_path)
            stats_df = self.generate_summary_statistics(access_df, output_path)
            

            self.generate_combined_tables_from_files(
                cai_penalty_path,
                cai_no_penalty_path,
                output_path
            )
            
            print("\n" + "=" * 60)
            print("✅ 所有表格生成完成！")
            print("=" * 60)
            

            print("\n🔍 数据完整性检查:")
            total_expected = 12 * 12  # 12 proteins × 12 variants
            access_count = access_df.notna().sum().sum()
            print(f"   Access实验: {access_count}/{total_expected} ({access_count/total_expected*100:.1f}%)")
            
            print("\n📁 输出文件:")
            for file in sorted(output_path.glob('*')):
                print(f"   - {file.name}")
        else:
            print("❌ 生成失败：没有找到有效数据")
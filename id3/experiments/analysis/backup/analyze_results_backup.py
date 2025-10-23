#!/usr/bin/env python3
"""


"""

import argparse
import sys
import os
import json
from pathlib import Path


sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from id3.experiments.analysis.fast_analyzer import FastExperimentAnalyzer
from id3.experiments.analysis.experiment_validator import ExperimentValidator
from id3.experiments.analysis.experiment_comparator import ExperimentComparator
from id3.experiments.analysis.generate_tab1_access_comparison import ParallelPaperPerformanceTableGenerator as PaperPerformanceTableGenerator
from id3.experiments.analysis.mpi_analyzer import MPIExperimentAnalyzer
import pandas as pd
import numpy as np

def analyze_single_dir(dir_path: Path):


    print("=" * 80)
    
    analyzer = FastExperimentAnalyzer()
    

    results = []
    json_files = list(dir_path.glob("*.json"))
    experiment_files = [f for f in json_files 
                        if not f.name.startswith(('config', 'progress', 'summary'))]
    

    
    for file_path in experiment_files:
        data = analyzer.extract_key_fields(file_path)
        if data:
            results.append(data)
    
    if not results:

        return None
    

    df = pd.DataFrame(results)
    






    

    if 'best_accessibility' in df:






        

        excellent = (df['best_accessibility'] < 1.0).sum()

    

    if 'amino_acids_match' in df:
        aa_match = df['amino_acids_match'].sum()


    if 'amino_acids_correct' in df:

    

    if 'final_ecai' in df and df['final_ecai'].notna().any():

        ecai_data = df['final_ecai'].dropna()




        
        if 'cai_target_achieved' in df:
            achieved = df['cai_target_achieved'].sum()

    

    if 'constraint_type' in df and 'best_accessibility' in df:

        for constraint in df['constraint_type'].unique():
            constraint_df = df[df['constraint_type'] == constraint]
            print(f"  {constraint.upper()}:")


    
    return df

def compare_dirs(dir1: Path, dir2: Path):
    """对比两个实验目录"""
    print(f"\n🔄 对比分析")
    print("=" * 80)
    
    df1 = analyze_single_dir(dir1)
    df2 = analyze_single_dir(dir2)
    
    if df1 is None or df2 is None:
        print("❌ 无法进行对比")
        return
    
    print(f"\n📊 性能对比:")
    print(f"  目录1 ({dir1.name}):")
    print(f"    • 平均: {df1['best_accessibility'].mean():.4f}")
    print(f"    • 中位数: {df1['best_accessibility'].median():.4f}")
    print(f"    • 最优: {df1['best_accessibility'].min():.4f}")
    
    print(f"  目录2 ({dir2.name}):")
    print(f"    • 平均: {df2['best_accessibility'].mean():.4f}")
    print(f"    • 中位数: {df2['best_accessibility'].median():.4f}")
    print(f"    • 最优: {df2['best_accessibility'].min():.4f}")
    

    diff_mean = df1['best_accessibility'].mean() - df2['best_accessibility'].mean()
    diff_median = df1['best_accessibility'].median() - df2['best_accessibility'].median()
    
    print(f"\n  差异 (目录1 - 目录2):")
    print(f"    • 平均差异: {diff_mean:+.4f} kcal/mol ({diff_mean/df2['best_accessibility'].mean()*100:+.1f}%)")
    print(f"    • 中位数差异: {diff_median:+.4f} kcal/mol")

def find_latest_results():

    results_dirs = []
    

    if Path("results_full").exists():
        results_dirs.extend(Path("results_full").glob("*"))
    

    if Path("results").exists():
        results_dirs.extend(Path("results").glob("*"))
    

    results_dirs = sorted(results_dirs, key=lambda x: x.stat().st_mtime, reverse=True)
    
    if results_dirs:
        return results_dirs[0]
    return None

def main():
    parser = argparse.ArgumentParser(

        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:

  python analyze_results.py
  

  python analyze_results.py results_full/20250910_004126_unified_access_experiments
  

  python analyze_results.py --compare dir1 dir2
  

  python analyze_results.py --all-complete
  

  python analyze_results.py --check-complete
  python analyze_results.py --check-complete results_full
  

  python analyze_results.py --check-cai results_full/20250910_004126_unified_access_experiments
  

  python analyze_results.py --compare-cai-access
  

  python analyze_results.py --verify-full
  

  python analyze_results.py --timeline-access


  python analyze_results.py --generate-paper-figures
  

  python analyze_results.py --mpi external_storage/results/mpi_flexible
  

  python analyze_results.py --compare-mpi external_storage/results/mpi_flexible_no_cai external_storage/results/mpi_flexible_cai_penalty
  

  python analyze_results.py --auto-detect-mpi
        """
    )
    
    parser.add_argument('directory', nargs='?',

    parser.add_argument('--compare', nargs=2, metavar=('DIR1', 'DIR2'),

    parser.add_argument('--all-complete', action='store_true',

    parser.add_argument('--check-complete', nargs='?', const='results', metavar='DIR',

    parser.add_argument('--check-cai', metavar='DIR',

    parser.add_argument('--compare-cai-access', action='store_true',

    parser.add_argument('--verify-full', action='store_true',

    parser.add_argument('--timeline-access', action='store_true',

    parser.add_argument('--verify-constraints', nargs='*', metavar='DIR',

    parser.add_argument('--visualize', metavar='DIR',

    parser.add_argument('--compare-visual', nargs='+', metavar='DIR',

    parser.add_argument('--standard', nargs='*', metavar='DIR',

    parser.add_argument('--paper-tables', nargs='+', metavar='DIR',

    parser.add_argument('--generate-all-tables', action='store_true',

    parser.add_argument('--generate-paper-figures', action='store_true',

    parser.add_argument('--mpi', metavar='DIR',

    parser.add_argument('--compare-mpi', nargs=2, metavar=('DIR1', 'DIR2'),

    parser.add_argument('--auto-detect-mpi', action='store_true',

    parser.add_argument('--mpi-summary', metavar='DIR',


    parser.add_argument('--verbose', action='store_true',

    
    args = parser.parse_args()
    

    validator = ExperimentValidator()
    comparator = ExperimentComparator()
    mpi_analyzer = MPIExperimentAnalyzer()


    if args.mpi_summary:
        mpi_dir = Path(args.mpi_summary)
        if not mpi_dir.exists():

            return 1

        output_csv = None
        if args.output:
            output_csv = Path(args.output)

        result_df = generate_mpi_summary_table(mpi_dir, output_csv)
        return 0 if result_df is not None else 1
    

    if args.check_complete is not None:
        stats = validator.check_experiment_completeness(args.check_complete)
        if 'error' in stats:
            print(f"❌ {stats['error']}")
            return 1
        validator.print_completeness_report(stats)
        return 0
    

    if args.check_cai:
        dir_path = Path(args.check_cai)
        if not dir_path.exists():

            return 1
        status = validator.check_cai_status(dir_path)
        validator.print_cai_status_report(dir_path, status)
        return 0
    

    if args.verify_full:
        stats = validator.verify_results_full()
        validator.print_results_full_report(stats)
        return 0
    

    if args.compare_cai_access:

        cai_dir = Path("results_full/20250909_003751_unified_cai_experiments")
        access_dir = Path("results_full/20250909_105804_unified_access_experiments")
        
        if not cai_dir.exists() or not access_dir.exists():

            return 1
        

        print("=" * 80)
        result = comparator.compare_cai_vs_access(cai_dir, access_dir)
        if 'error' not in result:
            comparator.print_comparison_report(result)
        else:
            print(f"❌ {result['error']}")
        return 0
    

    if args.verify_constraints is not None:
        from id3.analysis.verify_constraint_satisfaction import analyze_directory
        
        if len(args.verify_constraints) == 0:

            dirs_to_check = [


            ]
        else:
            dirs_to_check = [(Path(d).name, Path(d)) for d in args.verify_constraints]
        

        print("=" * 80)
        
        for desc, dir_path in dirs_to_check:
            if dir_path.exists():
                print(f"\n🔍 {desc}")
                analyze_directory(dir_path, sample_size=10)
            else:

        return 0
    

    if args.visualize:
        from id3.analysis.performance_visualizer import PerformanceVisualizer
        
        dir_path = Path(args.visualize)
        if not dir_path.exists():

            return 1
        

        visualizer = PerformanceVisualizer()
        df = visualizer.load_experiment_data(dir_path)
        
        if len(df) == 0:

            return 1
        

        output_dir = dir_path / 'analysis' / 'visualizations'
        output_dir.mkdir(parents=True, exist_ok=True)
        

        visualizer.plot_constraint_comparison(df, output_dir / 'constraint_comparison.png')
        

        visualizer.plot_performance_distribution(df, output_dir / 'performance_distribution.png')
        

        visualizer.plot_protein_performance(df, output_dir / 'protein_performance.png')
        

        report = visualizer.generate_performance_report(df)
        





        

        for constraint, stats in report['by_constraint'].items():

        

        with open(output_dir / 'performance_report.json', 'w') as f:
            json.dump(report, f, indent=2)
        

        return 0
    

    if args.compare_visual:
        from id3.analysis.performance_visualizer import PerformanceVisualizer
        
        experiments = []
        for dir_str in args.compare_visual:
            dir_path = Path(dir_str)
            if dir_path.exists():
                experiments.append((dir_path.name[:20], dir_path))
            else:

        
        if len(experiments) < 2:

            return 1
        

        visualizer = PerformanceVisualizer()
        

        output_path = Path('comparison_visualization.png')
        visualizer.compare_experiments(experiments, output_path)
        

        return 0
    

    if args.standard is not None:
        from id3.analysis.standard_analysis import run_standard_analysis
        
        if len(args.standard) == 0:

            dirs_to_analyze = [
                Path('results_full/20250909_105804_unified_access_experiments'),
                Path('results_full/20250910_004126_unified_access_experiments'),
                Path('results_full/20250910_022355_unified_access_experiments')
            ]

            if Path('results/20250910_123342_unified_cai_experiments').exists():
                dirs_to_analyze.append(Path('results/20250910_123342_unified_cai_experiments'))
        else:
            dirs_to_analyze = [Path(d) for d in args.standard]
        

        valid_dirs = [d for d in dirs_to_analyze if d.exists()]
        
        if not valid_dirs:

            return 1
        


        

        output_dir = Path('analysis_output') / 'standard_analysis'
        

        run_standard_analysis(valid_dirs, output_dir)
        


        return 0
    

    if args.generate_all_tables:
        from id3.experiments.analysis.paper_table_generator import PaperTableGenerator
        

        paper_results_dir = Path('paper_experiment_results')
        if not paper_results_dir.exists():

            return 1
        

        access_dir = paper_results_dir / 'access_only'
        cai_penalty_dir = paper_results_dir / 'cai_with_penalty'
        cai_no_penalty_dir = paper_results_dir / 'cai_no_penalty'
        
        if not all([access_dir.exists(), cai_penalty_dir.exists(), cai_no_penalty_dir.exists()]):


            return 1
        

        print("=" * 80)
        

        output_dir = paper_results_dir / 'tables'
        output_dir.mkdir(parents=True, exist_ok=True)
        

        table_gen = PaperTableGenerator(num_workers=40)
        


        table_gen.run(
            str(access_dir),
            str(cai_penalty_dir),
            str(cai_no_penalty_dir),
            str(output_dir)
        )
        







        return 0
    

    if args.paper_tables:
        from id3.experiments.analysis.paper_table_generator import PaperTableGenerator
        

        if len(args.paper_tables) == 3:

            access_dir = Path(args.paper_tables[0])
            cai_penalty_dir = Path(args.paper_tables[1])
            cai_no_penalty_dir = Path(args.paper_tables[2])
            
            if not all([access_dir.exists(), cai_penalty_dir.exists(), cai_no_penalty_dir.exists()]):

                return 1
            

            print("=" * 80)
            

            output_dir = Path('paper_experiment_results/tables')
            output_dir.mkdir(parents=True, exist_ok=True)
            

            table_gen = PaperTableGenerator(num_workers=40)
            

            table_gen.run(
                str(access_dir),
                str(cai_penalty_dir),
                str(cai_no_penalty_dir),
                str(output_dir)
            )
            
        elif len(args.paper_tables) == 2:

            access_dir = Path(args.paper_tables[0])
            cai_dir = Path(args.paper_tables[1])
            
            if not access_dir.exists() or not cai_dir.exists():

                return 1
            

            print("=" * 80)
            

            output_dir = Path('paper_experiment_results/tables')
            output_dir.mkdir(parents=True, exist_ok=True)
            

            table_gen = PaperPerformanceTableGenerator()
            


            csv_path = output_dir / 'performance_12x12.csv'
            table_gen.run(str(access_dir), str(cai_dir), str(csv_path))
            


        else:



            return 1
        
        return 0


    if args.generate_paper_figures:
        from id3.experiments.analysis.paper_figures import PaperFigureGenerator


        print("=" * 80)


        data_dir = Path('paper_experiment_results')
        output_dir = data_dir / 'figures'

        if not data_dir.exists():

            return 1


        figure_gen = PaperFigureGenerator(data_dir, output_dir)






        figure_gen.generate_all_figures()








        return 0


    if args.mpi:
        mpi_dir = Path(args.mpi)
        if not mpi_dir.exists():

            return 1
        
        if not mpi_analyzer.is_mpi_results_dir(mpi_dir):

            return 1
        
        analysis = mpi_analyzer.analyze_mpi_experiments(mpi_dir)
        mpi_analyzer.print_analysis_report(analysis)
        

        if args.output and 'dataframe' in analysis:
            df = analysis['dataframe']
            output_path = Path(args.output)
            if output_path.suffix == '.csv':
                df.to_csv(output_path, index=False)

            elif output_path.suffix == '.json':
                df.to_json(output_path, orient='records', indent=2)

        
        return 0
    

    if args.compare_mpi:
        dir1 = Path(args.compare_mpi[0])
        dir2 = Path(args.compare_mpi[1])
        
        if not dir1.exists() or not dir2.exists():

            return 1
        
        if not mpi_analyzer.is_mpi_results_dir(dir1) or not mpi_analyzer.is_mpi_results_dir(dir2):

            return 1
        
        comparison = mpi_analyzer.compare_mpi_experiments(dir1, dir2)
        mpi_analyzer.print_comparison_report(comparison)
        return 0
    

    if args.auto_detect_mpi:



        mpi_base = Path("external_storage/results")
        mpi_dirs = []

        if mpi_base.exists():
            for dir_path in mpi_base.glob("mpi_flexible*"):
                if mpi_analyzer.is_mpi_results_dir(dir_path):
                    mpi_dirs.append(dir_path)

        if not mpi_dirs:

            return 1


        for i, mpi_dir in enumerate(mpi_dirs, 1):
            print(f"\n{'='*60}")

            print('='*60)

            analysis = mpi_analyzer.analyze_mpi_experiments(mpi_dir)
            mpi_analyzer.print_analysis_report(analysis)

        return 0

def generate_mpi_summary_table(mpi_results_dir: Path, output_csv: Path = None):
    """
    汇总MPI实验结果，生成类似data_tab1_access_comparison.csv的表格
    计算不同seed的平均值和标准差
    """



    all_results = []


    seed_dirs = list(mpi_results_dir.glob("seed*"))


    for seed_dir in seed_dirs:
        seed_num = seed_dir.name.replace('seed', '')


        worker_dirs = list(seed_dir.glob("worker*"))

        for worker_dir in worker_dirs:

            worker_name = worker_dir.name
            parts = worker_name.split('_')
            if len(parts) >= 5:
                variant = parts[1]  # v00, v01, v10, v11
                constraint = parts[4][1:]  # clagrangian -> lagrangian


                protein_dirs = [d for d in worker_dir.iterdir() if d.is_dir()]

                for protein_dir in protein_dirs:

                    json_files = [f for f in protein_dir.glob("*.json")
                                if not f.name.startswith(('config', 'progress', 'summary'))]

                    if json_files:


                        try:
                            with open(result_file) as f:
                                data = json.load(f)


                            result = {
                                'protein': data.get('protein_name', 'unknown'),
                                'constraint': constraint,
                                'variant': variant,
                                'seed': int(seed_num),
                                'best_accessibility': data.get('best_accessibility', float('inf')),
                                'amino_acids_match': data.get('amino_acids_match', False),
                                'amino_acids_correct': data.get('amino_acids_correct', 0.0),
                                'final_accessibility': data.get('final_accessibility', float('inf')),
                                'improvement': data.get('improvement', 0.0)
                            }
                            all_results.append(result)

                        except Exception as e:


    if not all_results:

        return None


    df = pd.DataFrame(all_results)







    match_rate = df['amino_acids_match'].mean() * 100



    grouped = df.groupby(['protein', 'constraint', 'variant'])['best_accessibility'].agg(['mean', 'std', 'count']).reset_index()


    proteins = sorted(df['protein'].unique())
    constraints = sorted(df['constraint'].unique())
    variants = sorted(df['variant'].unique())


    columns = ['Protein']
    for constraint in constraints:
        for variant in variants:
            columns.append(f"{constraint}_{variant}")


    result_rows = []

    for protein in proteins:
        row = {'Protein': protein}

        for constraint in constraints:
            for variant in variants:

                key = f"{constraint}_{variant}"
                mask = (grouped['protein'] == protein) & \
                       (grouped['constraint'] == constraint) & \
                       (grouped['variant'] == variant)

                matching_rows = grouped[mask]
                if len(matching_rows) > 0:
                    mean_val = matching_rows.iloc[0]['mean']
                    std_val = matching_rows.iloc[0]['std']
                    count_val = matching_rows.iloc[0]['count']


                    if pd.notna(std_val) and count_val > 1:
                        row[key] = f"{mean_val:.4f}±{std_val:.4f}"
                    else:
                        row[key] = f"{mean_val:.4f}"
                else:
                    row[key] = "N/A"

        result_rows.append(row)


    result_df = pd.DataFrame(result_rows)


    if output_csv is None:
        output_csv = mpi_results_dir / f"mpi_summary_{mpi_results_dir.name}.csv"

    result_df.to_csv(output_csv, index=False)



    print("=" * 100)
    print(result_df.to_string(index=False))










    overall_mean = df['best_accessibility'].mean()
    overall_std = df['best_accessibility'].std()
    overall_min = df['best_accessibility'].min()
    overall_max = df['best_accessibility'].max()








    for constraint in constraints:
        constraint_data = df[df['constraint'] == constraint]['best_accessibility']
        print(f"  {constraint.upper()}:")






    return result_df

def main():
        dir1 = Path(args.compare[0])
        dir2 = Path(args.compare[1])
        if not dir1.exists() or not dir2.exists():

            return 1
        compare_dirs(dir1, dir2)
        return 0
    

    if args.all_complete:

        complete_dirs = []
        
        for base_dir in ['results_full', 'results']:
            if Path(base_dir).exists():
                for dir_path in Path(base_dir).glob("*"):
                    if dir_path.is_dir():
                        json_count = len(list(dir_path.glob("*.json")))

                        exp_count = json_count - 3 if json_count > 3 else json_count
                        if exp_count >= 144:
                            complete_dirs.append(dir_path)
        

        for dir_path in complete_dirs:
            analyze_single_dir(dir_path)
        return 0
    

    if args.directory:
        dir_path = Path(args.directory)
        if not dir_path.exists():

            return 1
    else:

        dir_path = find_latest_results()
        if not dir_path:


            return 1

    
    df = analyze_single_dir(dir_path)
    

    if args.output and df is not None:
        output_path = Path(args.output)
        if output_path.suffix == '.csv':
            df.to_csv(output_path, index=False)

        elif output_path.suffix == '.json':
            df.to_json(output_path, orient='records', indent=2)

    
    return 0

if __name__ == "__main__":
    sys.exit(main())
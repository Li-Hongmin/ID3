#!/usr/bin/env python3
"""


"""

import json
import sys
import os
from pathlib import Path
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

def process_experiment_file(file_info):
    """

    file_info: (result_file_path, protein, constraint, variant, seed)
    """
    result_file, protein, constraint, variant, seed = file_info

    try:
        with open(result_file) as f:
            data = json.load(f)


        return {
            'protein': data.get('protein_name', protein),
            'constraint': constraint,
            'variant': variant,
            'seed': seed,
            'best_accessibility': data.get('best_accessibility', float('inf')),
            'amino_acids_match': data.get('amino_acids_match', False),
            'amino_acids_correct': data.get('amino_acids_correct', 0.0),
            'final_accessibility': data.get('final_accessibility', float('inf')),
            'improvement': data.get('improvement', 0.0),
            'file_path': str(result_file)
        }
    except Exception as e:
        print(f"⚠️ 解析文件失败 {result_file}: {e}")
        return None

def collect_experiment_files(mpi_results_dir: Path):
    """


    """
    file_info_list = []


    seed_dirs = list(mpi_results_dir.glob("seed*"))

    for seed_dir in seed_dirs:
        seed_num = int(seed_dir.name.replace('seed', ''))


        worker_dirs = list(seed_dir.glob("worker*"))

        for worker_dir in worker_dirs:


            worker_name = worker_dir.name
            parts = worker_name.split('_')
            if len(parts) >= 4:
                variant = parts[1][1:]  # v00 -> 00
                constraint = parts[3][1:]  # clagrangian -> lagrangian


                protein_dirs = [d for d in worker_dir.iterdir() if d.is_dir()]

                for protein_dir in protein_dirs:

                    protein_name_from_dir = protein_dir.name.split('_')[0]


                    json_files = [f for f in protein_dir.glob("*.json")
                                if not f.name.startswith(('config', 'progress', 'summary'))]

                    if json_files:
                        result_file = json_files[0]
                        file_info_list.append((
                            result_file,
                            protein_name_from_dir,
                            constraint,
                            variant,
                            seed_num
                        ))

    return file_info_list

def generate_mpi_summary_table(mpi_results_dir: Path, output_csv: Path = None, num_workers: int = None):
    """


    """
    print(f"🔍 汇总MPI实验结果: {mpi_results_dir}")

    if num_workers is None:
        num_workers = min(cpu_count(), 32)

    print(f"使用 {num_workers} 个并行进程")


    print("📁 收集实验文件...")
    file_info_list = collect_experiment_files(mpi_results_dir)

    if not file_info_list:
        print("❌ 没有找到实验文件")
        return None

    print(f"找到 {len(file_info_list)} 个实验文件")


    print("⚡ 并行处理实验文件...")
    all_results = []
    result_dict = defaultdict(list)  # {(protein, constraint, variant): [values]}

    with ProcessPoolExecutor(max_workers=num_workers) as executor:

        future_to_file = {executor.submit(process_experiment_file, file_info): file_info for file_info in file_info_list}


        completed = 0
        for future in as_completed(future_to_file):
            completed += 1
            if completed % 50 == 0:
                print(f"  已处理: {completed}/{len(file_info_list)}")

            result = future.result()
            if result is not None:
                all_results.append(result)


                protein = result['protein']
                constraint = result['constraint']
                variant = result['variant']
                best_accessibility = result['best_accessibility']

                key = (protein, constraint, variant)
                result_dict[key].append(best_accessibility)

    print(f"✅ 并行处理完成，共处理 {len(all_results)} 个有效实验")

    if not all_results:
        print("❌ 没有找到有效的实验结果")
        return None

    print(f"收集到 {len(all_results)} 个实验结果")


    proteins = set(r['protein'] for r in all_results)
    constraints = set(r['constraint'] for r in all_results)
    variants = set(r['variant'] for r in all_results)
    seeds = set(r['seed'] for r in all_results)

    print(f"蛋白质数量: {len(proteins)}")
    print(f"约束类型: {sorted(list(constraints))}")
    print(f"变体类型: {sorted(list(variants))}")
    print(f"随机种子: {sorted(list(seeds))}")


    match_rate = sum(r['amino_acids_match'] for r in all_results) / len(all_results) * 100
    print(f"氨基酸完全匹配率: {match_rate:.1f}%")


    proteins_sorted = sorted(list(proteins))
    constraints_sorted = sorted(list(constraints))
    variants_sorted = sorted(list(variants))


    columns = ['Protein']
    for constraint in constraints_sorted:
        for variant in variants_sorted:
            columns.append(f"{constraint}_{variant}")


    result_rows = []

    for protein in proteins_sorted:
        row = [protein]

        for constraint in constraints_sorted:
            for variant in variants_sorted:
                key = (protein, constraint, variant)
                if key in result_dict and len(result_dict[key]) > 0:
                    values = result_dict[key]
                    mean_val = np.mean(values)


                    if len(values) > 1:
                        std_val = np.std(values, ddof=1)
                        row.append(f"{mean_val:.4f}±{std_val:.4f}")
                    else:
                        row.append(f"{mean_val:.4f}±0.0000")
                else:
                    row.append("N/A")

        result_rows.append(row)


    if output_csv is None:
        output_csv = mpi_results_dir / f"mpi_summary_{mpi_results_dir.name}.csv"

    with open(output_csv, 'w') as f:

        f.write(','.join(columns) + '\n')

        for row in result_rows:
            f.write(','.join(map(str, row)) + '\n')


    print(f"\n📊 MPI实验汇总表格预览:")
    print("=" * 100)


    print(','.join(columns))


    for row in result_rows:
        print(','.join(map(str, row)))


    print(f"\n📈 统计摘要:")
    print(f"  • 蛋白质数量: {len(proteins_sorted)}")
    print(f"  • 约束类型: {len(constraints_sorted)} ({', '.join(constraints_sorted)})")
    print(f"  • 变体类型: {len(variants_sorted)} ({', '.join(variants_sorted)})")
    print(f"  • 随机种子数: {len(seeds)}")
    print(f"  • 总实验次数: {len(all_results)}")


    all_values = [r['best_accessibility'] for r in all_results]
    overall_mean = np.mean(all_values)
    overall_std = np.std(all_values, ddof=1)
    overall_min = np.min(all_values)
    overall_max = np.max(all_values)

    print(f"\n🎯 整体性能 (best_accessibility):")
    print(f"  • 平均值: {overall_mean:.4f} ± {overall_std:.4f} kcal/mol")
    print(f"  • 最优值: {overall_min:.4f} kcal/mol")
    print(f"  • 最差值: {overall_max:.4f} kcal/mol")


    print(f"\n📊 按约束类型统计:")
    for constraint in constraints_sorted:
        constraint_values = [r['best_accessibility'] for r in all_results if r['constraint'] == constraint]
        if constraint_values:
            mean_val = np.mean(constraint_values)
            std_val = np.std(constraint_values, ddof=1) if len(constraint_values) > 1 else 0
            min_val = np.min(constraint_values)
            print(f"  {constraint.upper()}:")
            print(f"    • 平均: {mean_val:.4f} ± {std_val:.4f}")
            print(f"    • 最优: {min_val:.4f}")
            print(f"    • 实验数: {len(constraint_values)}")

    print(f"\n✅ 汇总表格已保存到: {output_csv}")

    return result_rows

def main():
    if len(sys.argv) < 2:
        print("用法: python mpi_summary.py <mpi_results_dir> [output.csv] [num_workers]")
        print("示例:")
        print("  python mpi_summary.py external_storage/results/mpi_flexible")
        print("  python mpi_summary.py external_storage/results/mpi_flexible results.csv")
        print("  python mpi_summary.py external_storage/results/mpi_flexible results.csv 16")
        return 1

    mpi_dir = Path(sys.argv[1])
    if not mpi_dir.exists():
        print(f"❌ MPI实验目录不存在: {mpi_dir}")
        return 1

    output_csv = None
    if len(sys.argv) > 2 and sys.argv[2] != 'auto':
        output_csv = Path(sys.argv[2])

    num_workers = None
    if len(sys.argv) > 3:
        try:
            num_workers = int(sys.argv[3])
        except ValueError:
            print(f"⚠️ 无效的进程数: {sys.argv[3]}, 使用默认值")

    result = generate_mpi_summary_table(mpi_dir, output_csv, num_workers)
    return 0 if result is not None else 1

if __name__ == "__main__":
    sys.exit(main())
#!/usr/bin/env python3
"""


"""

import json
from pathlib import Path
import pandas as pd

def analyze_experiment_gaps():




    target_seeds = {



    }

    current_seeds = {



    }


    missing_seeds = {}
    for exp_type in target_seeds:
        missing = [s for s in target_seeds[exp_type] if s not in current_seeds[exp_type]]
        missing_seeds[exp_type] = missing


    total_missing = 0
    for exp_type, seeds in missing_seeds.items():
        completed = len(current_seeds[exp_type])
        target = len(target_seeds[exp_type])
        missing_count = len(seeds)
        total_missing += missing_count

        print(f"  {exp_type}:")





    return missing_seeds

def estimate_computation_cost(missing_seeds):
    """估算计算成本"""
    print("\n💻 计算成本估算:")


    base_experiments_per_seed = 1150 // 8
    base_time_per_experiment = 2.5

    for exp_type, seeds in missing_seeds.items():
        num_seeds = len(seeds)

        if exp_type == 'access_only':

            experiments_per_seed = 144
        else:

            experiments_per_seed = 90

        total_experiments = num_seeds * experiments_per_seed
        estimated_hours = (total_experiments * base_time_per_experiment) / 3600

        print(f"  {exp_type}:")
        print(f"    • 缺失seeds: {num_seeds}")
        print(f"    • 估算实验数: {total_experiments}")
        print(f"    • 估算时间: {estimated_hours:.1f} 小时")


        if num_seeds >= 4:
            parallel_hours = estimated_hours / min(num_seeds, 8)
            print(f"    • 并行时间: {parallel_hours:.1f} 小时 (8进程)")

def generate_execution_plan(missing_seeds):



    phases = []


    if missing_seeds['access_only']:
        phases.append({

            'seeds': missing_seeds['access_only'],
            'priority': 'HIGH',


        })





    phases.append({

        'seeds': {'cai_penalty': cai_penalty_partial, 'cai_no_penalty': cai_no_penalty_partial},
        'priority': 'MEDIUM',


    })


    remaining_cai_penalty = missing_seeds['cai_penalty'][8:]
    remaining_cai_no_penalty = missing_seeds['cai_no_penalty'][8:]

    if remaining_cai_penalty or remaining_cai_no_penalty:
        phases.append({

            'seeds': {'cai_penalty': remaining_cai_penalty, 'cai_no_penalty': remaining_cai_no_penalty},
            'priority': 'LOW',


        })

    for i, phase in enumerate(phases, 1):



        if isinstance(phase['seeds'], dict):
            for exp_type, seeds in phase['seeds'].items():
                if seeds:
                    print(f"    • {exp_type}: {len(seeds)} seeds {seeds[:5]}{'...' if len(seeds) > 5 else ''}")
        else:
            print(f"    • seeds: {len(phase['seeds'])} {phase['seeds'][:5]}{'...' if len(phase['seeds']) > 5 else ''}")

    return phases

def create_safe_execution_commands(missing_seeds):
    """生成安全的执行命令（考虑进程管理）"""
    print("\n🛡️ 安全执行策略:")

    print("  进程管理考虑:")
    print("    • 单任务执行：避免多个MPI任务冲突")
    print("    • 资源监控：每个任务前检查CPU/内存使用")
    print("    • 进程限制：最多使用系统50%资源")
    print("    • 错误恢复：每个seed独立，失败不影响其他")

    commands = []


    if missing_seeds['access_only']:
        print(f"\n  Phase 1 执行命令：")
        for seed in missing_seeds['access_only'][:3]:
            cmd = f"""

python run_unified_experiment.py --seeds {seed} --seeds {seed} --preset full-12x12 \\
    --output-dir external_storage/results/mpi_flexible_extended \\
    --max-workers 8 --timeout 7200


python mpi_summary.py external_storage/results/mpi_flexible_extended \\
    paper_experiment_results/tables/progress_check_seed_{seed}.csv 4
            """
            commands.append(cmd.strip())
            print(f"    Seed {seed}: python run_unified_experiment.py --seeds {seed} --preset full-12x12")

    return commands

def main():
    print("=" * 80)
    print("🧬 ID3-DeepRaccess实验补充计划制定")
    print("=" * 80)


    missing_seeds = analyze_experiment_gaps()


    estimate_computation_cost(missing_seeds)


    phases = generate_execution_plan(missing_seeds)


    commands = create_safe_execution_commands(missing_seeds)

    print("\n" + "=" * 80)
    print("📝 下一步行动建议:")
    print("  1. 先执行3个Access Only seeds作为测试")
    print("  2. 监控资源使用和完成时间")
    print("  3. 根据测试结果调整并行度")
    print("  4. 逐步扩展到完整计划")
    print("=" * 80)

    return missing_seeds, phases, commands

if __name__ == "__main__":
    main()
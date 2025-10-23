#!/usr/bin/env python3
"""
Unified ID3-DeepRaccess Experiment Entry Point (Pure Router)

This is a clean entry point that only handles:
1. Command-line argument parsing
2. Configuration setup
3. Experiment runner invocation
4. Result reporting

All business logic is delegated to modules in id3.experiments.
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

# Add project path
sys.path.append(str(Path(__file__).parent))

# Import modularized components
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import (
    UnifiedExperimentConfig,
    ExperimentPresets
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def _generate_performance_table(successful_results: list, config: UnifiedExperimentConfig) -> dict:
    """
    Generate detailed performance table for experiment combinations.
    
    Args:
        successful_results: List of successful experiment results
        config: Experiment configuration
        
    Returns:
        Dictionary containing performance tables
    """
    # Group results by protein, constraint, variant
    table_data = {}
    
    # Create structure: protein -> constraint -> variant -> metrics
    for result in successful_results:
        protein = result['protein_name']
        constraint = result['constraint_type']
        variant = result['variant']
        
        if protein not in table_data:
            table_data[protein] = {}
        if constraint not in table_data[protein]:
            table_data[protein][constraint] = {}
        if variant not in table_data[protein][constraint]:
            table_data[protein][constraint][variant] = []
        
        # Store key metrics for this experiment
        metrics = {
            'seed': result['seed'],
            'final_accessibility': result['final_accessibility'],
            'best_accessibility': result['best_accessibility'],
            'improvement': result.get('improvement', 0),
            'amino_acids_correct': result.get('amino_acids_correct', 0),
            'optimization_time': result.get('optimization_time', 0)
        }
        
        # Add CAI metrics if available
        if config.enable_cai and 'final_ecai' in result:
            metrics.update({
                'final_ecai': result['final_ecai'],
                'cai_target_achieved': result.get('cai_target_achieved', False),
                'ecai_improvement': result.get('ecai_improvement', 0)
            })
        
        table_data[protein][constraint][variant].append(metrics)
    
    # Generate summary tables
    performance_table = {
        'detailed_results': table_data,
        'best_accessibility_table': {},
        'constraint_comparison': {},
        'variant_comparison': {}
    }
    
    # Generate best accessibility table (protein × constraint)
    best_table = {}
    for protein in table_data:
        best_table[protein] = {}
        for constraint in table_data[protein]:
            # Find best accessibility across all variants and seeds for this protein-constraint combo
            all_accessibilities = []
            for variant in table_data[protein][constraint]:
                for metrics in table_data[protein][constraint][variant]:
                    all_accessibilities.append(metrics['best_accessibility'])
            
            best_table[protein][constraint] = {
                'best_accessibility': min(all_accessibilities) if all_accessibilities else float('inf'),
                'variant_count': len(table_data[protein][constraint]),
                'total_experiments': sum(len(table_data[protein][constraint][v]) for v in table_data[protein][constraint])
            }
    
    performance_table['best_accessibility_table'] = best_table
    
    # Generate constraint comparison (average performance per constraint)
    constraint_stats = {}
    for constraint in config.constraints:
        constraint_accessibilities = []
        for protein in table_data:
            if constraint in table_data[protein]:
                for variant in table_data[protein][constraint]:
                    for metrics in table_data[protein][constraint][variant]:
                        constraint_accessibilities.append(metrics['best_accessibility'])
        
        if constraint_accessibilities:
            constraint_stats[constraint] = {
                'average_best': sum(constraint_accessibilities) / len(constraint_accessibilities),
                'best_overall': min(constraint_accessibilities),
                'worst_overall': max(constraint_accessibilities),
                'experiment_count': len(constraint_accessibilities)
            }
    
    performance_table['constraint_comparison'] = constraint_stats
    
    # Generate variant comparison (average performance per variant)
    variant_stats = {}
    for variant in config.variants:
        variant_accessibilities = []
        for protein in table_data:
            for constraint in table_data[protein]:
                if variant in table_data[protein][constraint]:
                    for metrics in table_data[protein][constraint][variant]:
                        variant_accessibilities.append(metrics['best_accessibility'])
        
        if variant_accessibilities:
            variant_stats[variant] = {
                'average_best': sum(variant_accessibilities) / len(variant_accessibilities),
                'best_overall': min(variant_accessibilities),
                'worst_overall': max(variant_accessibilities),
                'experiment_count': len(variant_accessibilities)
            }
    
    performance_table['variant_comparison'] = variant_stats
    
    return performance_table


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description='Unified ID3-DeepRaccess Experiment Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick test
  python run_unified_experiment.py --preset quick-test
  
  # Standard CAI optimization
  python run_unified_experiment.py --proteins P04637,P0DTC2 --enable-cai
  
  # Custom configuration
  python run_unified_experiment.py \\
      --proteins P04637,P0DTC9 \\
      --constraints lagrangian,cpc \\
      --variants 11 \\
      --iterations 1000 \\
      --enable-cai \\
      --cai-target 0.8 \\
      --lambda-cai 0.1
        """
    )
    
    # Preset configurations
    parser.add_argument('--preset', type=str, 
                       choices=['quick-test', 'quick-test-cai-penalty', 'quick-test-cai-no-penalty', 'quick-test-cai-comparison', 'quick-test-both', 
                                'full-12x12', 'full-12x12-cai-penalty', 'full-12x12-cai-no-penalty', 'full-12x12-cai-comparison', 'full-12x12-both', 
                                'adaptive-lambda-cai-test'],
                       help='Use a preset configuration (full-12x12: without CAI, full-12x12-cai-penalty: CAI with constraint penalty, full-12x12-cai-no-penalty: CAI without constraint penalty, *-comparison: compare CAI modes, *-both: CAI vs no-CAI)')
    
    # Basic parameters
    # 支持单数和复数形式
    parser.add_argument('--proteins', '--protein', type=str,
                       help='Protein(s), comma-separated for multiple')
    parser.add_argument('--constraints', '--constraint', type=str,
                       help='Constraint type(s), comma-separated (lagrangian,ams,cpc)')
    parser.add_argument('--variants', '--variant', type=str,
                       help='Variant(s), comma-separated (00,01,10,11)')
    
    # Optimization parameters
    parser.add_argument('--iterations', type=int,
                       help='Number of optimization iterations')
    parser.add_argument('--learning-rate', type=float,
                       help='Learning rate')
    parser.add_argument('--batch-size', type=int, default=1,
                       help='Batch size (default: 1)')
    parser.add_argument('--seeds', type=int,
                       help='Number of random seeds')
    parser.add_argument('--base-seed', '--seed', type=int,
                       help='Starting seed value (default: 42)')
    
    # CAI parameters
    parser.add_argument('--enable-cai', action='store_true',
                       help='Enable unified CAI optimization')
    parser.add_argument('--cai-target', type=float,
                       help='Target CAI value')
    parser.add_argument('--lambda-cai', '--lambda_cai', type=float,
                       help='Weight for CAI loss (λ_CAI coefficient)')
    parser.add_argument('--species', type=str,
                       help='Species for CAI weights')
    parser.add_argument('--disable-constraint-penalty', action='store_true',
                       help='Disable constraint penalty in CAI optimization (CAI no penalty mode)')
    
    # Adaptive lambda_cai parameters
    parser.add_argument('--adaptive-lambda-cai', action='store_true',
                       help='Enable adaptive lambda_cai optimization')
    parser.add_argument('--lambda-cai-lr', type=float, default=0.1,
                       help='Learning rate for lambda_cai adaptation (default: 0.1)')
    parser.add_argument('--lambda-cai-max', type=float, default=2.0,
                       help='Maximum allowed lambda_cai value (default: 2.0)')
    parser.add_argument('--lambda-cai-min', type=float, default=0.01,
                       help='Minimum allowed lambda_cai value (default: 0.01)')
    parser.add_argument('--cai-tolerance', type=float, default=0.05,
                       help='Tolerance for CAI target satisfaction (default: 0.05)')
    parser.add_argument('--smoothing-factor', type=float, default=0.9,
                       help='Exponential smoothing factor for CAI measurements (default: 0.9)')
    
    # Execution parameters
    parser.add_argument('--device', type=str,
                       help='Computation device (cuda/cpu)')
    parser.add_argument('--disable-inner-tqdm', action='store_true',
                       help='Disable inner tqdm progress bars to reduce log verbosity')
    # 移除parallel参数，始终使用串行执行以获得最佳性能
    parser.add_argument('--verbose', action='store_true',
                       help='Show detailed progress')
    parser.add_argument('--mixed-precision', action='store_true',
                       help='Enable full mixed precision training (all or nothing)')
    
    # Output parameters
    parser.add_argument('--output-dir', type=str,
                       help='Output directory')
    parser.add_argument('--no-trajectories', action='store_true',
                       help='Do not save optimization trajectories')
    parser.add_argument('--force', action='store_true',
                       help='Force run even if GPU has other Python processes')
    
    return parser.parse_args()


def load_configuration(args: argparse.Namespace):
    """
    Load configuration from arguments or preset.
    
    Args:
        args: Command-line arguments
        
    Returns:
        UnifiedExperimentConfig instance or list of configs for dual mode
    """
    def apply_arg_overrides(config: UnifiedExperimentConfig, args: argparse.Namespace):
        """Apply command-line argument overrides to a configuration."""
        if args.proteins:
            config.proteins = args.proteins.split(',')
        if args.constraints:
            config.constraints = args.constraints.split(',')
        if args.variants:
            config.variants = args.variants.split(',')
        if args.iterations is not None:
            config.iterations = args.iterations
        if args.learning_rate is not None:
            config.learning_rate = args.learning_rate
        if args.seeds is not None:
            config.seeds = args.seeds
        if args.base_seed is not None:
            config.base_seed = args.base_seed
            # 如果只指定了seed而没有指定seeds，则设置seeds=1
            if args.seeds is None:
                config.seeds = 1
        if args.enable_cai:
            config.enable_cai = True
        if args.cai_target is not None:
            config.cai_target = args.cai_target
        if args.lambda_cai is not None:
            config.lambda_cai = args.lambda_cai
        if args.disable_constraint_penalty:
            config.disable_constraint_penalty = True
        if args.device:
            config.device = args.device
        # 移除parallel设置，始终串行执行
        if args.verbose:
            config.verbose = True
        if args.mixed_precision:
            config.mixed_precision = True
        if args.output_dir:
            config.output_dir = args.output_dir
        if args.no_trajectories:
            config.save_trajectories = False
        return config

    # Load preset if specified
    if args.preset:
        if args.preset == 'quick-test':
            config = ExperimentPresets.quick_test()
            return apply_arg_overrides(config, args)
        elif args.preset == 'quick-test-cai-penalty':
            config = ExperimentPresets.quick_test_cai_penalty()
            return apply_arg_overrides(config, args)
        elif args.preset == 'quick-test-cai-no-penalty':
            config = ExperimentPresets.quick_test_cai_no_penalty()
            return apply_arg_overrides(config, args)
        elif args.preset == 'quick-test-cai-comparison':
            # Apply overrides to both configs in CAI comparison mode
            configs = ExperimentPresets.quick_test_cai_comparison()
            return [apply_arg_overrides(config, args) for config in configs]
        elif args.preset == 'quick-test-both':
            # Apply overrides to both configs in dual mode
            configs = ExperimentPresets.quick_test_both()
            return [apply_arg_overrides(config, args) for config in configs]
        elif args.preset == 'full-12x12':
            config = ExperimentPresets.full_12x12()
            return apply_arg_overrides(config, args)
        elif args.preset == 'full-12x12-cai-penalty':
            config = ExperimentPresets.full_12x12_cai_penalty()
            return apply_arg_overrides(config, args)
        elif args.preset == 'full-12x12-cai-no-penalty':
            config = ExperimentPresets.full_12x12_cai_no_penalty()
            return apply_arg_overrides(config, args)
        elif args.preset == 'full-12x12-cai-comparison':
            # Apply overrides to both configs in CAI comparison mode
            configs = ExperimentPresets.full_12x12_cai_comparison()
            return [apply_arg_overrides(config, args) for config in configs]
        elif args.preset == 'full-12x12-both':
            # Apply overrides to both configs in dual mode
            configs = ExperimentPresets.full_12x12_both()
            return [apply_arg_overrides(config, args) for config in configs]
        elif args.preset == 'adaptive-lambda-cai-test':
            config = ExperimentPresets.adaptive_lambda_cai_test()
            return apply_arg_overrides(config, args)
        else:
            raise ValueError(f"Unknown preset: {args.preset}")
    else:
        # Create configuration from arguments
        config = UnifiedExperimentConfig.from_args(args)
        return apply_arg_overrides(config, args)


# 注意：save_individual_experiment 函数已移至 UnifiedExperimentRunner 类中
# 作为 _save_experiment_result 方法进行增量保存


def save_summary_results(config: UnifiedExperimentConfig, results: list, 
                        total_time: float, experiment_files: list) -> Path:
    """
    Save experiment summary and statistics to summary.json.
    
    Args:
        config: Experiment configuration
        results: List of all experiment results
        total_time: Total execution time
        experiment_files: List of individual experiment file names
        
    Returns:
        Path to saved summary file
    """
    output_dir = config.get_output_dir()
    
    # Calculate statistics for successful experiments
    successful_results = [r for r in results if r.get('status') == 'completed']
    failed_results = [r for r in results if r.get('status') == 'failed']
    
    summary_data = {
        'experiment_overview': {
            'total_experiments': len(results),
            'successful': len(successful_results),
            'failed': len(failed_results),
            'total_time': total_time,
            'average_time_per_experiment': total_time / len(results) if results else 0
        },
        'configuration': config.to_dict(),
        'experiment_files': [f.name if isinstance(f, Path) else f for f in experiment_files],  # 只保存文件名
        'statistics': {},
        'performance_table': {}  # 添加性能表格
    }
    
    # Add accessibility statistics
    if successful_results:
        accessibilities = [r['final_accessibility'] for r in successful_results]
        summary_data['statistics']['accessibility'] = {
            'average': sum(accessibilities) / len(accessibilities),
            'best': min(accessibilities),
            'worst': max(accessibilities),
            'std': 0.0  # 可以后续添加标准差计算
        }
        
        # Add CAI statistics if enabled
        if config.enable_cai:
            cai_results = [r for r in successful_results if 'final_ecai' in r]
            if cai_results:
                ecai_values = [r['final_ecai'] for r in cai_results]
                summary_data['statistics']['cai'] = {
                    'target': config.cai_target,
                    'lambda_cai': config.lambda_cai,
                    'average_ecai': sum(ecai_values) / len(ecai_values),
                    'target_achieved_count': sum(1 for r in cai_results if r.get('cai_target_achieved', False)),
                    'target_achieved_rate': sum(1 for r in cai_results if r.get('cai_target_achieved', False)) / len(cai_results)
                }
    
    # Add failed experiments info if any
    if failed_results:
        summary_data['failed_experiments'] = [
            {
                'experiment_id': f"{r['protein_name']}-{r['constraint_type']}-{r['variant']}-{r['seed']}",
                'error': r.get('error', 'Unknown error')
            }
            for r in failed_results
        ]
    
    # Generate detailed performance table
    summary_data['performance_table'] = _generate_performance_table(successful_results, config)
    
    # Save summary file
    summary_file = output_dir / 'summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    return summary_file


def save_initial_config(config: UnifiedExperimentConfig) -> Path:
    """
    Save initial configuration immediately after directory creation.
    
    Args:
        config: Experiment configuration
        
    Returns:
        Path to saved config file
    """
    output_dir = config.get_output_dir()  # This creates the directory
    
    # Save configuration file immediately
    config_file = output_dir / 'config.json'
    with open(config_file, 'w') as f:
        json.dump(config.to_dict(), f, indent=2)
    
    logger.info(f"💾 Configuration saved to: {config_file}")
    return config_file


def save_results(config: UnifiedExperimentConfig, results: list, total_time: float) -> Path:
    """
    Save summary results after experiments are complete.
    Individual experiments are already saved incrementally during execution.
    
    Args:
        config: Experiment configuration
        results: List of experiment results
        total_time: Total execution time
        
    Returns:
        Path to saved summary file
    """
    output_dir = config.get_output_dir()
    
    # 收集已保存的实验文件（增量保存已完成）
    experiment_files = list(output_dir.glob('*.json'))
    # 排除config.json和summary.json
    experiment_files = [f for f in experiment_files 
                       if f.name not in ['config.json', 'summary.json']]
    
    # Save configuration file (skip if already exists from save_initial_config)
    config_file = output_dir / 'config.json'
    if not config_file.exists():
        with open(config_file, 'w') as f:
            json.dump(config.to_dict(), f, indent=2)
    
    # Save summary file
    summary_file = save_summary_results(config, results, total_time, experiment_files)
    
    return summary_file


def print_summary(config: UnifiedExperimentConfig, results: list, total_time: float, summary_file: Path) -> None:
    """
    Print experiment summary.
    
    Args:
        config: Experiment configuration
        results: List of experiment results
        total_time: Total execution time
        summary_file: Path to saved summary file
    """
    successful = sum(1 for r in results if r.get('status') == 'completed')
    failed = sum(1 for r in results if r.get('status') == 'failed')
    
    print(f"\n📊 Experiment Summary")
    print("=" * 60)
    print(f"Successful: {successful}/{len(results)}")
    print(f"Failed: {failed}/{len(results)}")
    
    if successful > 0:
        completed_results = [r for r in results if r.get('status') == 'completed']
        avg_accessibility = sum(r['final_accessibility'] for r in completed_results) / len(completed_results)
        best_accessibility = min(r['final_accessibility'] for r in completed_results)
        
        print(f"Average accessibility: {avg_accessibility:.4f}")
        print(f"Best accessibility: {best_accessibility:.4f}")
        
        if config.enable_cai:
            cai_results = [r for r in completed_results if 'final_ecai' in r]
            if cai_results:
                avg_ecai = sum(r['final_ecai'] for r in cai_results) / len(cai_results)
                target_achieved = sum(1 for r in cai_results if r.get('cai_target_achieved', False))
                print(f"Average ECAI: {avg_ecai:.4f}")
                print(f"CAI target achieved: {target_achieved}/{len(cai_results)} "
                      f"({target_achieved/len(cai_results)*100:.1f}%)")
    
    print(f"Total time: {total_time:.1f}s")
    print(f"Summary saved: {summary_file}")
    print(f"Individual experiments saved in: {summary_file.parent}/")
    
    print(f"\n🎉 Experiment complete!")
    if config.enable_cai:
        if config.disable_constraint_penalty:
            print("   • CAI No Penalty mode: L_total = L_Access + λ_CAI * L_CAI")
            print("   • Constraint penalty disabled for pure CAI optimization")
        else:
            print("   • CAI With Penalty mode: L_total = L_Access + λ_CAI * L_CAI + constraint_penalties")
            print("   • Full constraint enforcement with CAI optimization")
        print("   • Used ECAI formula for gradient-based CAI optimization")
        print("   • Applied CAI enhancement at β=1 discretization")
    else:
        print("   • Pure accessibility optimization with constraint penalties")


def check_gpu_processes():
    """
    检查GPU是否有其他Python进程在运行。
    如果有，则报错退出，避免多个实验同时运行导致GPU内存溢出。
    
    Returns:
        bool: True if GPU is available, False if occupied
    """
    try:
        # 检查nvidia-smi是否可用
        result = subprocess.run(
            ['nvidia-smi', '--query-compute-apps=pid,name', '--format=csv,noheader'],
            capture_output=True, text=True, timeout=5
        )
        
        if result.returncode != 0:
            logger.warning("⚠️  无法检查GPU状态，nvidia-smi命令失败")
            return True  # 如果无法检查，继续运行
        
        # 解析输出，查找Python进程
        lines = result.stdout.strip().split('\n')
        python_processes = []
        current_pid = os.getpid()
        
        for line in lines:
            if line.strip():
                parts = line.strip().split(',')
                if len(parts) >= 2:
                    pid = int(parts[0].strip())
                    process_name = parts[1].strip()
                    
                    # 排除当前进程
                    if pid != current_pid and 'python' in process_name.lower():
                        python_processes.append((pid, process_name))
        
        if python_processes:
            logger.error("❌ 检测到GPU上有其他Python进程正在运行！")
            logger.error("   为避免GPU内存溢出，请先停止其他实验。")
            logger.error("   检测到的进程：")
            for pid, name in python_processes:
                logger.error(f"     - PID {pid}: {name}")
            logger.error("\n   可以使用以下命令停止这些进程：")
            logger.error(f"     kill {' '.join(str(p[0]) for p in python_processes)}")
            return False
        
        logger.info("✅ GPU可用，没有其他Python进程占用")
        return True
        
    except subprocess.TimeoutExpired:
        logger.warning("⚠️  nvidia-smi命令超时，跳过GPU检查")
        return True
    except FileNotFoundError:
        logger.warning("⚠️  未找到nvidia-smi，可能没有安装NVIDIA驱动，跳过GPU检查")
        return True
    except Exception as e:
        logger.warning(f"⚠️  检查GPU状态时出错: {e}，跳过检查")
        return True


def main():
    """
    Main function - Pure router that delegates to modules.
    
    This function only:
    1. Checks GPU availability
    2. Parses arguments
    3. Loads configuration
    4. Creates and runs experiment runner
    5. Saves and reports results
    """
    # Parse arguments
    args = parse_arguments()
    
    # 添加 --force 参数跳过GPU检查
    if not getattr(args, 'force', False):
        if not check_gpu_processes():
            logger.error("\n💡 提示：如果确实要强制运行，可以使用 --force 参数跳过GPU检查")
            sys.exit(1)
    
    # Load configuration
    try:
        config_or_configs = load_configuration(args)
        
        # Handle dual mode (list of configs) or single mode
        if isinstance(config_or_configs, list):
            # Dual mode: run both CAI and non-CAI experiments
            configs = config_or_configs
            print(f"🔄 Running dual mode preset '{args.preset}' with {len(configs)} configuration sets")
            
            all_results = []
            total_start_time = time.time()
            
            for i, config in enumerate(configs):
                config.validate()
                mode_name = "CAI" if config.enable_cai else "NoCAI"
                print(f"\n🚀 Running configuration {i+1}/{len(configs)}: {mode_name} mode")
                config.print_summary()
                
                # Save configuration immediately after directory creation
                save_initial_config(config)
                
                # Generate and run experiments for this config
                experiments = config.generate_experiments()
                logger.info(f"🎯 Running {len(experiments)} experiments...")
                
                # 确保配置中包含输出目录路径
                runner_config = config.to_dict()
                runner_config['output_dir'] = str(config.get_output_dir())
                runner = UnifiedExperimentRunner(runner_config)
                start_time = time.time()
                results = runner.run_batch(experiments)
                config_time = time.time() - start_time
                
                # Add mode information to results
                for result in results:
                    result['mode'] = mode_name
                    result['enable_cai'] = config.enable_cai
                
                all_results.extend(results)
                
                # Save results for this config
                summary_file = save_results(config, results, config_time)
                print_summary(config, results, config_time, summary_file)
            
            total_time = time.time() - total_start_time
            print(f"\n🎉 All dual mode experiments complete! Total time: {total_time:.1f}s")
            print(f"   • Total experiments: {len(all_results)}")
            print(f"   • CAI experiments: {sum(1 for r in all_results if r.get('enable_cai'))}")
            print(f"   • Non-CAI experiments: {sum(1 for r in all_results if not r.get('enable_cai'))}")
            
        else:
            # Single mode: run one configuration
            config = config_or_configs
            config.validate()
            
            # Print configuration summary
            config.print_summary()
            
            # Save configuration immediately after directory creation
            save_initial_config(config)
            
            # Generate experiments
            experiments = config.generate_experiments()
            logger.info(f"🎯 Running {len(experiments)} experiments...")
            
            # Create runner and execute experiments
            # 确保配置中包含输出目录路径
            runner_config = config.to_dict()
            runner_config['output_dir'] = str(config.get_output_dir())
            runner = UnifiedExperimentRunner(runner_config)
            
            start_time = time.time()
            results = runner.run_batch(experiments)
            total_time = time.time() - start_time
            
            # Save results
            summary_file = save_results(config, results, total_time)
            
            # Print summary
            print_summary(config, results, total_time, summary_file)
        
    except ValueError as e:
        logger.error(f"❌ Configuration error: {e}")
        sys.exit(1)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
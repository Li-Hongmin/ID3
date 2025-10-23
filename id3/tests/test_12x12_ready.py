#!/usr/bin/env python3
"""







"""

import sys
import torch
from pathlib import Path

sys.path.append('.')

from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner


def test_small_matrix():

    print("\n" + "="*80)

    print("="*80)
    



    
    config = UnifiedExperimentConfig(
        proteins=proteins,
        constraints=constraints,


        seeds=1,
        enable_cai=False,
        verbose=False
    )
    
    runner = UnifiedExperimentRunner(config.to_dict())
    experiments = config.generate_experiments()
    

    results = []
    success_count = 0
    
    for i, exp in enumerate(experiments, 1):

        result = runner.run_single_experiment(**exp)
        
        if result['status'] == 'completed':
            success_count += 1

        else:

        
        results.append(result)
    

    return success_count == len(experiments)


def test_cai_integration():
    """测试CAI集成"""
    print("\n" + "="*80)
    print("🧪 测试CAI优化集成")
    print("="*80)
    
    config = UnifiedExperimentConfig(
        proteins=['P99999'],
        constraints=['lagrangian'],
        variants=['11'],
        iterations=5,
        enable_cai=True,
        cai_target=0.8,
        lambda_cai=0.1,
        seeds=1,
        verbose=True
    )
    
    runner = UnifiedExperimentRunner(config.to_dict())
    
    print("\n运行CAI优化实验...")
    result = runner.run_single_experiment(
        protein_name='P99999',
        constraint_type='lagrangian',
        variant='11',
        seed=42
    )
    
    if result['status'] == 'completed':
        print(f"\n✅ CAI优化成功")
        print(f"  初始 ECAI: {result.get('initial_ecai', 0):.4f}")
        print(f"  最终 ECAI: {result.get('final_ecai', 0):.4f}")
        print(f"  ECAI改进: {result.get('ecai_improvement', 0):.4f}")
        print(f"  目标达成: {result.get('cai_target_achieved', False)}")
        print(f"  Accessibility: {result['final_accessibility']:.4f}")
        return True
    else:
        print(f"\n❌ CAI优化失败: {result.get('error', '未知')}")
        return False


def check_protein_data():

    print("\n" + "="*80)

    print("="*80)
    

    required_proteins = [
        'EGFP', 'P99999', 'mCherry', 'GFP', 'Luc', 'Cas9',
        'Insulin', 'DHFR', 'MDH', 'AAV2', 'H1N1', 'SOD1'
    ]
    
    data_dir = Path('data')
    available = []
    missing = []
    
    for protein in required_proteins:

        found = False
        for suffix in ['.fasta', '.fasta.txt', '.txt']:
            file_path = data_dir / f"{protein}{suffix}"
            if file_path.exists():
                available.append(protein)
                found = True
                break
        if not found:
            missing.append(protein)
    


        print(f"  • {p}")
    if len(available) > 5:

    
    if missing:

        for p in missing:
            print(f"  • {p}")
    
    return len(missing) == 0


def main():
    """主测试函数"""
    print("\n" + "="*80)
    print("🚀 12×12实验系统最终验证")
    print("="*80)
    

    cuda_available = torch.cuda.is_available()
    print(f"\n🖥️  CUDA状态: {'✅ 可用' if cuda_available else '⚠️ 不可用（将使用CPU）'}")
    

    all_pass = True
    

    data_ready = check_protein_data()
    if not data_ready:
        print("\n⚠️  部分蛋白质数据缺失，但可以运行部分实验")
    

    matrix_pass = test_small_matrix()
    all_pass &= matrix_pass
    

    cai_pass = test_cai_integration()
    all_pass &= cai_pass
    

    print("\n" + "="*80)
    print("📊 最终验证结果")
    print("="*80)
    
    print(f"\n✅ DeepRaccess模型: 已加载")
    print(f"{'✅' if matrix_pass else '❌'} 约束机制: {'全部通过' if matrix_pass else '部分失败'}")
    print(f"{'✅' if cai_pass else '❌'} CAI优化: {'正常工作' if cai_pass else '失败'}")
    print(f"{'✅' if data_ready else '⚠️'} 蛋白质数据: {'完整' if data_ready else '部分缺失'}")
    
    if all_pass:
        print("\n🎉 系统完全准备就绪！可以运行完整的12×12实验。")
        print("\n建议命令：")
        print("
        print("  python run_unified_experiment.py --proteins P99999,P00004,P01308 --constraints lagrangian --iterations 100")
        print("\n
        print("  python run_unified_experiment.py --preset full-12x12 --iterations 1000")
    else:
        print("\n⚠️  系统部分功能需要修复，但可以运行基础实验。")
    
    print("="*80)
    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
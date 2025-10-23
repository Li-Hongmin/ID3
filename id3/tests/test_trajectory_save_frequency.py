#!/usr/bin/env python3


import json
import tempfile
from pathlib import Path
import sys


sys.path.append(str(Path(__file__).parent.parent.parent))

from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig


def test_trajectory_save_frequency():
    """测试轨迹是否每个迭代都保存"""
    

    config = UnifiedExperimentConfig(
        proteins=['O15263'],
        constraints=['lagrangian'],
        variants=['00'],
        iterations=10,
        seeds=1,
        base_seed=42,
        enable_cai=False,
        verbose=False,
        save_trajectories=True,
        output_dir=tempfile.mkdtemp()
    )
    

    runner_config = config.to_dict()
    runner_config['output_dir'] = str(config.get_output_dir())
    runner = UnifiedExperimentRunner(runner_config)
    

    experiments = config.generate_experiments()
    print(f"生成了 {len(experiments)} 个实验")
    

    results = runner.run_batch(experiments)
    

    if results and len(results) > 0:
        result = results[0]
        trajectory = result.get('trajectory', {})
        num_iterations = len(trajectory.get('iterations', []))
        expected_iterations = config.iterations
        
        print(f"配置的迭代次数: {expected_iterations}")
        print(f"实际保存的轨迹点数: {num_iterations}")
        print(f"轨迹迭代列表: {trajectory.get('iterations', [])[:20]}")
        
        if num_iterations == expected_iterations:
            print("✅ 测试通过: 每个迭代都保存了轨迹")
            return True
        else:
            print(f"❌ 测试失败: 期望 {expected_iterations} 个轨迹点，实际 {num_iterations} 个")
            return False
    else:
        print("❌ 测试失败: 没有返回结果")
        return False


if __name__ == "__main__":
    success = test_trajectory_save_frequency()
    sys.exit(0 if success else 1)
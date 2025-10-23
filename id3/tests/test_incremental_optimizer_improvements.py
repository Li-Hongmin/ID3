#!/usr/bin/env python3
"""








"""

import sys
import torch
import numpy as np
import logging
from pathlib import Path


project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from id3.optimizers.cai.incremental import IncrementalCAIOptimizer


logging.basicConfig(level=logging.DEBUG, format='%(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def test_constants():



    optimizer = IncrementalCAIOptimizer()


    constants = [
        'MIN_K_POSITIONS', 'MAX_K_POSITIONS', 'K_POSITION_DIVISOR',
        'CONSECUTIVE_FAILURE_THRESHOLD', 'RANDOM_REINIT_PROBABILITY',
        'NUMERICAL_EPSILON', 'CAI_EPSILON', 'PROBABILITY_EPSILON'
    ]

    for const in constants:

        value = getattr(optimizer, const)
        print(f"✓ {const} = {value}")




def test_input_validation():
    """测试输入验证"""
    print("\n=== 测试输入验证 ===")

    optimizer = IncrementalCAIOptimizer()


    try:
        optimizer._validate_inputs(
            torch.randn(10, 64),
            "ACGT",
            1.5
        )
        assert False, "应该抛出ValueError"
    except ValueError as e:
        print(f"✓ 正确捕获无效CAI值: {e}")


    try:
        optimizer._validate_inputs(
            torch.randn(10, 64),
            "",
            0.8
        )
        assert False, "应该抛出ValueError"
    except ValueError as e:
        print(f"✓ 正确捕获空序列: {e}")

    print("✅ 输入验证功能正常")


def test_numerical_stability():



    optimizer = IncrementalCAIOptimizer()




    probs = np.zeros(10)
    normalized = optimizer._normalize_probabilities(probs)




    probs = np.array([1.0, np.nan, 2.0, 3.0])
    normalized = optimizer._normalize_probabilities(probs)





    probs = np.array([1.0, np.inf, 2.0, 3.0])
    normalized = optimizer._normalize_probabilities(probs)





    probs = np.array([1e-100, 1e-99, 1e-98])
    normalized = optimizer._normalize_probabilities(probs)






def test_caching():
    """测试缓存功能"""
    print("\n=== 测试缓存功能 ===")


    optimizer_cached = IncrementalCAIOptimizer(
        amino_acid_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT",
        enable_caching=True
    )


    optimizer_no_cache = IncrementalCAIOptimizer(
        amino_acid_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT",
        enable_caching=False
    )


    test_indices = np.random.randint(0, 4, size=50)


    cai1 = optimizer_cached._compute_cai_from_indices(test_indices)
    assert optimizer_cached._cache_misses == 1, "应有1次缓存未命中"
    print(f"✓ 第一次计算: CAI={cai1:.4f}, 缓存未命中={optimizer_cached._cache_misses}")


    cai2 = optimizer_cached._compute_cai_from_indices(test_indices)
    assert optimizer_cached._cache_hits == 1, "应有1次缓存命中"
    assert cai1 == cai2, "缓存返回值应相同"
    print(f"✓ 第二次计算: CAI={cai2:.4f}, 缓存命中={optimizer_cached._cache_hits}")


    stats = optimizer_cached.get_statistics()
    assert 'cache_hit_rate' in stats, "统计信息应包含缓存命中率"
    print(f"✓ 缓存命中率: {stats['cache_hit_rate']:.2%}")


    optimizer_cached.clear_cache()
    assert len(optimizer_cached._cai_cache) == 0, "缓存应被清空"
    print("✓ 缓存清空成功")

    print("✅ 缓存功能测试通过")


def test_refactored_functions():



    optimizer = IncrementalCAIOptimizer(
        amino_acid_sequence="MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT"
    )


    k_pos = optimizer._calculate_k_positions(200)
    assert optimizer.MIN_K_POSITIONS <= k_pos <= optimizer.MAX_K_POSITIONS
    print(f"✓ _calculate_k_positions(200) = {k_pos}")

    k_prob = optimizer._calculate_k_probability_positions(200)
    assert optimizer.MIN_K_PROBABILITY <= k_prob <= optimizer.MAX_K_PROBABILITY
    print(f"✓ _calculate_k_probability_positions(200) = {k_prob}")

    perturb = optimizer._calculate_perturb_positions(100)
    assert 1 <= perturb <= optimizer.MAX_PERTURB_POSITIONS
    print(f"✓ _calculate_perturb_positions(100) = {perturb}")




def test_incremental_optimization():
    """测试增量优化主流程"""
    print("\n=== 测试增量优化流程 ===")


    seq_len = 50
    amino_acid_sequence = "M" * seq_len

    optimizer = IncrementalCAIOptimizer(
        amino_acid_sequence=amino_acid_sequence,
        random_seed=42,
        enable_caching=True
    )


    pi_accessibility = torch.rand(seq_len, 64)
    target_cai = 0.8


    valid_codon_mask = torch.ones(seq_len, 64)


    result1, metadata1 = optimizer.optimize(
        pi_accessibility,
        target_cai,
        amino_acid_sequence,
        valid_codon_mask=valid_codon_mask
    )

    assert metadata1['method'] == 'incremental_bs', "第一次应使用Binary Search"
    assert metadata1['iteration'] == 1
    print(f"✓ 第一次优化: 方法={metadata1['method']}, CAI={metadata1.get('final_cai', 'N/A')}")


    result2, metadata2 = optimizer.optimize(
        pi_accessibility,
        target_cai,
        amino_acid_sequence,
        valid_codon_mask=valid_codon_mask
    )

    assert metadata2['method'] == 'incremental_inc', "第二次应使用增量优化"
    assert metadata2['iteration'] == 2
    print(f"✓ 第二次优化: 方法={metadata2['method']}, 增量变化={metadata2.get('incremental_changes', 0)}")


    stats = optimizer.get_statistics()
    print(f"✓ 总迭代次数: {stats['iteration_count']}")
    print(f"✓ 唯一序列数: {stats['unique_sequences']}")
    if 'cache_hit_rate' in stats:
        print(f"✓ 缓存命中率: {stats['cache_hit_rate']:.2%}")

    print("✅ 增量优化流程测试通过")


def run_all_tests():

    print("=" * 60)

    print("=" * 60)

    test_functions = [
        test_constants,
        test_input_validation,
        test_numerical_stability,
        test_caching,
        test_refactored_functions,
        test_incremental_optimization
    ]

    failed = []
    for test_func in test_functions:
        try:
            test_func()
        except Exception as e:
            failed.append((test_func.__name__, str(e)))


    print("\n" + "=" * 60)
    if not failed:

    else:

        for name, error in failed:
            print(f"  - {name}: {error}")
    print("=" * 60)

    return len(failed) == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
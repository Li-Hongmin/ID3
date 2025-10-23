"""



"""

import sys
import torch
import time
from pathlib import Path


project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from id3.cai.unified_calculator import UnifiedCAICalculator, compute_cai
from id3.utils.logging_config import setup_logging, PerformanceLogger


logger = setup_logging(level='INFO', colored=True, name='verify')


def test_unified_cai_calculator():

    print("\n" + "="*60)

    print("="*60)
    

    test_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    test_rna = "AUGAGCAAGGGTGAAGAACTGTTCACCGGTGTTGTGCCGATCCTGGTTGAACTG"
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    

    calc = UnifiedCAICalculator(species='ecoli_bl21de3', device=device)

    

    cai = calc.compute_cai(test_rna, method='standard')

    

    result = calc.compute_cai(test_rna, method='standard', return_log=True)

    

    seq_len = len(test_sequence)
    num_codons = 6
    codon_probs = torch.rand(seq_len, num_codons, device=device)
    codon_probs = codon_probs / codon_probs.sum(dim=-1, keepdim=True)
    
    diff_cai = calc.compute_cai(
        codon_probs, 
        method='differentiable',
        amino_acid_sequence=test_sequence
    )

    

    batch_probs = torch.rand(4, seq_len, num_codons, device=device)
    batch_probs = batch_probs / batch_probs.sum(dim=-1, keepdim=True)
    batch_cai = calc.compute_cai(batch_probs, method='batch')

    

    max_cai = calc.get_max_achievable_cai(test_sequence)

    

    quick_cai = compute_cai(test_rna)

    
    return True


def test_logging_config():
    """测试日志配置"""
    print("\n" + "="*60)
    print("测试日志配置")
    print("="*60)
    

    test_logger = setup_logging(level='DEBUG', name='test_module')
    
    test_logger.debug("调试信息")
    test_logger.info("信息日志")
    test_logger.warning("警告信息")
    test_logger.error("错误信息")
    

    with PerformanceLogger(test_logger, "测试操作"):
        time.sleep(0.1)
    
    logger.info("✅ 日志配置测试完成")
    
    return True


def test_performance_comparison():

    print("\n" + "="*60)

    print("="*60)
    
    test_rna = "AUGAGCAAGGGTGAAGAACTGTTCACCGGTGTTGTGCCGATCCTGGTTGAACTG" * 10
    test_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG" * 10
    

    calc = UnifiedCAICalculator()
    

    start = time.time()
    for _ in range(100):
        cai = calc.compute_cai(test_rna, method='standard')
    new_time = time.time() - start

    

    try:
        from id3.cai.validator import compute_cai_from_sequence
        
        start = time.time()
        for _ in range(100):
            cai = compute_cai_from_sequence(test_rna)
        old_time = time.time() - start

        
        improvement = (old_time - new_time) / old_time * 100

    except ImportError:

    
    return True


def test_integration():
    """集成测试"""
    print("\n" + "="*60)
    print("集成测试")
    print("="*60)
    
    test_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    try:

        from id3.constraints.lagrangian import LagrangianConstraint
        
        constraint = LagrangianConstraint(
            amino_acid_sequence=test_sequence,
            protein_sequence=test_sequence,
            enable_cai=True,
            target_cai=0.8,
            device=device
        )
        

        result = constraint.forward(alpha=1.0, tau=1.0)
        
        if 'prob' in result:
            logger.info(f"✅ 约束输出形状: {result['prob'].shape}")
        

        calc = UnifiedCAICalculator(device=device)
        
        if 'enhanced_sequence' in result and result['enhanced_sequence'] is not None:
            logger.info("✅ 增强序列生成成功")
        
        logger.info("✅ 与约束集成测试通过")
    except Exception as e:
        logger.error(f"集成测试失败: {e}")
        return False
    
    return True


def main():

    print("\n" + "="*80)

    print("="*80)
    
    results = []
    

    tests = [




    ]
    
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:

            results.append((name, False))
    

    print("\n" + "="*80)

    print("="*80)
    
    total = len(results)
    passed = sum(1 for _, success in results if success)
    
    for name, success in results:

        print(f"{name:20} {status}")
    

    
    if passed == total:






    else:

    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
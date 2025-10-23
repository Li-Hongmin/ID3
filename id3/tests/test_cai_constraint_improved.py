"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.constraints.cai_constraint_improved import CAIConstraint, CAIEnhancedPsiFunction
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_cai_improved')


def test_cai_constraint_basic():

    logger.info("="*60)

    logger.info("="*60)
    

    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    constraint = CAIConstraint(
        amino_acid_sequence=amino_acid_sequence,
        target_cai=0.8,
        lambda_cai=0.1,
        species='ecoli_bl21de3',
        device=device,
        optimizer_type='binary_search'
    )
    

    seq_len = len(amino_acid_sequence)
    num_codons = 6
    batch_size = 1
    
    codon_probs = torch.rand(batch_size, seq_len, num_codons, device=device)

    codon_probs = codon_probs / codon_probs.sum(dim=-1, keepdim=True)
    


    result = constraint.forward(codon_probs, beta=1.0, compute_loss=True)
    
    assert 'discrete_sequence' in result
    assert 'metadata' in result
    assert 'cai_loss' in result
    assert 'total_loss' in result
    



    

    metadata = result['metadata']



    
    return True


def test_cai_psi_function():
    """测试CAI增强的Psi函数"""
    logger.info("\n2. 测试CAI增强的Psi函数")
    
    amino_acid_sequence = "MKAI"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    psi = CAIEnhancedPsiFunction(
        amino_acid_sequence=amino_acid_sequence,
        species='ecoli_bl21de3',
        device=device,
        optimizer_type='binary_search'
    )
    

    seq_len = len(amino_acid_sequence)
    num_codons = 6
    codon_probs = torch.rand(1, seq_len, num_codons, device=device)
    codon_probs = codon_probs / codon_probs.sum(dim=-1, keepdim=True)
    

    test_cases = [
        {'enable_cai': False, 'beta': 0.0, 'name': 'No CAI + Soft'},
        {'enable_cai': False, 'beta': 1.0, 'name': 'No CAI + Hard'},
        {'enable_cai': True, 'beta': 0.0, 'name': 'CAI + Soft'},
        {'enable_cai': True, 'beta': 1.0, 'name': 'CAI + Hard (STE)'},
    ]
    
    for case in test_cases:
        discrete_seq, metadata = psi.apply(
            codon_probs,
            target_cai=0.8,
            beta=case['beta'],
            enable_cai=case['enable_cai']
        )
        
        logger.info(f"   {case['name']}:")
        logger.info(f"     形状: {discrete_seq.shape}")
        logger.info(f"     CAI启用: {metadata['enable_cai']}")
        
        if case['enable_cai']:
            logger.info(f"     实际CAI: {metadata.get('actual_cai', 0):.3f}")
            logger.info(f"     满足约束: {metadata.get('constraint_satisfied', False)}")
    
    logger.info("   ✅ Psi函数测试通过")
    return True


def test_algorithm_consistency():


    
    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    torch.manual_seed(42)
    np.random.seed(42)
    

    seq_len = len(amino_acid_sequence)
    num_codons = 6
    codon_probs = torch.rand(1, seq_len, num_codons, device=device)
    codon_probs = codon_probs / codon_probs.sum(dim=-1, keepdim=True)
    codon_probs_copy = codon_probs.clone()
    

    constraint_improved = CAIConstraint(
        amino_acid_sequence=amino_acid_sequence,
        target_cai=0.8,
        lambda_cai=0.1,
        species='ecoli_bl21de3',
        device=device
    )
    
    result_improved = constraint_improved.forward(codon_probs, beta=1.0)
    

    try:
        from id3.constraints.cai_constraint import CAIConstraint as CAIConstraintOld
        

        torch.manual_seed(42)
        np.random.seed(42)
        
        constraint_old = CAIConstraintOld(
            amino_acid_sequence=amino_acid_sequence,
            target_cai=0.8,
            lambda_cai=0.1,
            species='ecoli_bl21de3',
            device=device
        )
        
        result_old = constraint_old.forward(codon_probs_copy, beta=1.0)
        

        if 'discrete_sequence' in result_old:
            diff = torch.abs(result_improved['discrete_sequence'] - result_old['discrete_sequence']).sum()

            
            if diff < 0.01:

            else:

        
    except ImportError:

    


    

    discrete_seq = result_improved['discrete_sequence']
    if discrete_seq.requires_grad:

    

    metadata = result_improved['metadata']
    if 'actual_cai' in metadata:

    

    if 'cai_loss' in result_improved:
        loss = result_improved['cai_loss']
        if loss >= 0:

    

    return True


def test_performance():
    """测试性能"""
    logger.info("\n4. 性能测试")
    
    import time
    
    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG" * 5
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    constraint = CAIConstraint(
        amino_acid_sequence=amino_acid_sequence,
        target_cai=0.8,
        device=device
    )
    

    seq_len = len(amino_acid_sequence)
    num_codons = 6
    codon_probs = torch.rand(1, seq_len, num_codons, device=device)
    codon_probs = codon_probs / codon_probs.sum(dim=-1, keepdim=True)
    

    for _ in range(5):
        _ = constraint.forward(codon_probs, beta=1.0)
    

    num_iterations = 20
    start_time = time.time()
    
    for _ in range(num_iterations):
        result = constraint.forward(codon_probs, beta=1.0)
    
    elapsed = time.time() - start_time
    avg_time = elapsed / num_iterations
    
    logger.info(f"   平均前向传播时间: {avg_time*1000:.2f} ms")
    logger.info(f"   序列长度: {seq_len}")
    logger.info(f"   吞吐量: {seq_len/avg_time:.0f} residues/s")
    
    if avg_time < 0.1:
        logger.info("   ✅ 性能满足要求")
    else:
        logger.warning("   ⚠️ 性能可能需要优化")
    
    return True


def main():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    tests = [




    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:

            results.append((name, False))
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:

        logger.info(f"{name:15} {status}")
    

    
    if passed == total:






    else:

    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
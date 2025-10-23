"""



"""

import pytest
import torch
import numpy as np
from id3.cai.unified_calculator import UnifiedCAICalculator, compute_cai
from id3.utils.logging_config import setup_logging, PerformanceLogger


logger = setup_logging(level='DEBUG', name='test_cai')


class TestUnifiedCAICalculator:

    
    def test_initialization(self, device):
        """测试初始化"""
        calc = UnifiedCAICalculator(species='ecoli_bl21de3', device=device)
        assert calc.species == 'ecoli_bl21de3'
        assert calc.device == device
        assert calc.wi_table is not None
        assert calc.weights_tensor is not None
        logger.info("✅ CAI calculator initialization successful")
    
    def test_standard_cai(self, test_rna_sequence):

        calc = UnifiedCAICalculator()
        

        cai = calc.compute_cai(test_rna_sequence, method='standard')
        

        assert isinstance(cai, float)
        assert 0.0 <= cai <= 1.0
        logger.info(f"✅ Standard CAI: {cai:.4f}")
    
    def test_standard_cai_with_log(self, test_rna_sequence):
        """测试带详细日志的标准CAI计算"""
        calc = UnifiedCAICalculator()
        

        result = calc.compute_cai(test_rna_sequence, method='standard', return_log=True)
        

        assert 'cai' in result
        assert 'num_codons' in result
        assert 'weights' in result
        assert 'codons' in result
        

        assert isinstance(result['cai'], float)
        assert result['num_codons'] == len(test_rna_sequence) // 3
        assert len(result['weights']) == result['num_codons']
        assert len(result['codons']) == result['num_codons']
        
        logger.info(f"✅ Standard CAI with log: {result['cai']:.4f}, {result['num_codons']} codons")
    
    def test_differentiable_cai(self, test_codon_probs, test_amino_sequence, device):

        calc = UnifiedCAICalculator(device=device)
        

        cai = calc.compute_cai(
            test_codon_probs, 
            method='differentiable',
            amino_acid_sequence=test_amino_sequence
        )
        

        assert isinstance(cai, torch.Tensor)

        assert 0.0 <= cai.item() <= 1.0
        
        logger.info(f"✅ Differentiable CAI: {cai.item():.4f}")
    
    def test_batch_cai(self, device):
        """测试批量CAI计算"""
        calc = UnifiedCAICalculator(device=device)
        

        batch_size = 4
        seq_len = 10
        num_codons = 6
        batch_probs = torch.rand(batch_size, seq_len, num_codons, device=device)
        batch_probs = batch_probs / batch_probs.sum(dim=-1, keepdim=True)
        

        cai_values = calc.compute_cai(batch_probs, method='batch')
        

        assert isinstance(cai_values, np.ndarray)
        assert len(cai_values) == batch_size
        assert all(0.0 <= v <= 1.0 for v in cai_values)
        
        logger.info(f"✅ Batch CAI: mean={cai_values.mean():.4f}, std={cai_values.std():.4f}")
    
    def test_max_achievable_cai(self, test_amino_sequence):

        calc = UnifiedCAICalculator()
        

        max_cai = calc.get_max_achievable_cai(test_amino_sequence)
        

        assert isinstance(max_cai, float)
        assert 0.0 <= max_cai <= 1.0
        
        logger.info(f"✅ Max achievable CAI: {max_cai:.4f}")
    
    def test_cai_validation(self, test_rna_sequence):
        """测试CAI约束验证"""
        calc = UnifiedCAICalculator()
        

        target_cais = [0.3, 0.5, 0.7]
        
        for target in target_cais:
            satisfied, actual = calc.validate_cai_constraint(
                test_rna_sequence, 
                target_cai=target,
                tolerance=0.01
            )
            
            logger.info(f"Target: {target:.2f}, Actual: {actual:.4f}, Satisfied: {satisfied}")
    
    def test_convenience_function(self, test_rna_sequence, test_codon_probs):


        cai1 = compute_cai(test_rna_sequence)
        assert isinstance(cai1, float)
        

        cai2 = compute_cai(test_codon_probs)
        assert isinstance(cai2, torch.Tensor)
        
        logger.info(f"✅ Convenience function: string={cai1:.4f}, tensor={cai2.item():.4f}")
    
    def test_cache_functionality(self, test_amino_sequence):
        """测试缓存功能"""
        calc = UnifiedCAICalculator()
        

        max_cai1 = calc.get_max_achievable_cai(test_amino_sequence)
        

        max_cai2 = calc.get_max_achievable_cai(test_amino_sequence)
        
        assert max_cai1 == max_cai2
        

        calc.clear_cache()
        logger.info("✅ Cache functionality working")
    
    def test_performance(self, test_rna_sequence, test_amino_sequence):

        calc = UnifiedCAICalculator()
        

        with PerformanceLogger(logger, "Standard CAI computation"):
            for _ in range(100):
                cai = calc.compute_cai(test_rna_sequence, method='standard')
        

        with PerformanceLogger(logger, "Max CAI computation"):
            for _ in range(100):
                max_cai = calc.get_max_achievable_cai(test_amino_sequence)
        
        logger.info("✅ Performance tests completed")


class TestIntegration:
    """集成测试"""
    
    def test_with_real_constraint(self, test_amino_sequence, device):

        from id3.constraints.lagrangian import LagrangianConstraint
        

        constraint = LagrangianConstraint(
            protein_sequence=test_amino_sequence,
            enable_cai=True,
            target_cai=0.8,
            device=device
        )
        

        calc = UnifiedCAICalculator(device=device)
        

        result = constraint.forward(alpha=1.0, tau=1.0)
        
        if 'enhanced_sequence' in result:

            enhanced_seq = result['enhanced_sequence']
            if isinstance(enhanced_seq, torch.Tensor):


                pass
        
        logger.info("✅ Integration with constraints successful")
    
    def test_comparison_with_old_methods(self, test_rna_sequence):
        """比较新旧方法的一致性"""

        calc = UnifiedCAICalculator()
        new_cai = calc.compute_cai(test_rna_sequence, method='standard')
        

        try:
            from id3.cai.validator import compute_cai_from_sequence
            old_cai = compute_cai_from_sequence(test_rna_sequence)
            

            assert abs(new_cai - old_cai) < 1e-6, f"CAI mismatch: new={new_cai:.6f}, old={old_cai:.6f}"
            logger.info(f"✅ Consistency check: new={new_cai:.6f}, old={old_cai:.6f}")
        except ImportError:
            logger.warning("Old CAI method not available for comparison")


@pytest.mark.parametrize("species", ['ecoli_bl21de3'])
def test_different_species(species):

    calc = UnifiedCAICalculator(species=species)
    

    test_seq = "AUGAAACCCUUUGGG"
    
    cai = calc.compute_cai(test_seq, method='standard')
    assert 0.0 <= cai <= 1.0
    
    logger.info(f"✅ Species {species}: CAI={cai:.4f}")


if __name__ == "__main__":

    pytest.main([__file__, "-v", "--tb=short"])
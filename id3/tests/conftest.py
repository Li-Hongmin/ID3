"""



"""

import pytest
import torch
import numpy as np
from pathlib import Path
import sys


project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))




@pytest.fixture(scope='session')
def device():

    return torch.device('cuda' if torch.cuda.is_available() else 'cpu')


@pytest.fixture(scope='session')
def cpu_device():
    """Force use of CPU device"""
    return torch.device('cpu')




@pytest.fixture
def test_amino_sequence():

    return "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"


@pytest.fixture
def test_rna_sequence():
    """Corresponding test RNA sequence"""

    return "AUGAGCAAGGGTGAAGAACTGTTCACCGGTGTTGTGCCGATCCTGGTTGAACTGGATGGTGATGTGAACGGTCACAAATTCAGCGTGAGCGGT"


@pytest.fixture
def short_amino_sequence():

    return "MKHELM"


@pytest.fixture
def long_amino_sequence():
    """Long test sequence (lysozyme fragment)"""
    return "KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQIN"




@pytest.fixture
def test_codon_probs(device):

    seq_len = 10
    num_codons = 6
    probs = torch.rand(seq_len, num_codons, device=device)

    probs = probs / probs.sum(dim=-1, keepdim=True)
    return probs


@pytest.fixture
def test_valid_mask(test_amino_sequence):
    """Generate valid codon mask"""
    from id3.utils.constants import amino_acids_to_codons

    seq_len = len(test_amino_sequence)
    max_codons = 6
    valid_mask = torch.zeros(seq_len, max_codons, dtype=torch.bool)

    for pos, aa in enumerate(test_amino_sequence):
        if aa in amino_acids_to_codons:
            num_codons = len(amino_acids_to_codons[aa])
            valid_mask[pos, :num_codons] = True

    return valid_mask




@pytest.fixture
def constraint_config():

    return {
        'alpha': 1.0,
        'beta': 0.0,
        'tau': 1.0,
        'target_cai': 0.8,
        'lambda_cai': 0.1,
        'species': 'ecoli_bl21de3'
    }


@pytest.fixture
def cai_config():
    """CAI-related configuration"""
    return {
        'species': 'ecoli_bl21de3',
        'target_cai': 0.8,
        'lambda_cai': 0.1,
        'enable_cai': True,
        'cai_method': 'sado'
    }




@pytest.fixture
def experiment_config():

    return {
        'protein': 'O15263',
        'constraint_type': 'lagrangian',
        'variant': '11',  # beta=1, tau=1
        'enable_cai': False,

        'learning_rate': 0.001,
        'batch_size': 1
    }




@pytest.fixture
def mock_deepraccess():
    """Mock DeepRaccess model (for tests that don't need the real model)"""
    class MockDeepRaccess:
        def __init__(self, device):
            self.device = device
            self.model = self

        def eval(self):
            return self

        def compute_atg_window_accessibility(self, rna_sequence, atg_position=0, discrete=False):

            if isinstance(rna_sequence, torch.Tensor):
                batch_size = rna_sequence.shape[0] if rna_sequence.dim() > 2 else 1
                if batch_size > 1:
                    return torch.rand(batch_size, device=self.device) * 10.0
                else:
                    return torch.rand(1, device=self.device).squeeze() * 10.0
            else:
                return torch.tensor(5.0, device=self.device)

    return MockDeepRaccess




@pytest.fixture
def temp_dir(tmp_path):

    return tmp_path


@pytest.fixture
def results_dir(temp_dir):
    """Test results directory"""
    results = temp_dir / "results"
    results.mkdir(exist_ok=True)
    return results




@pytest.fixture
def performance_thresholds():

    return {
        'cai_computation_time': 0.1,  # 100ms
        'constraint_forward_time': 0.5,  # 500ms
        'optimization_step_time': 1.0,  # 1s
        'memory_usage_mb': 1000,  # 1GB
    }




@pytest.fixture
def assert_tensor_close():
    """Provide assertion function for tensor approximate equality"""
    def _assert(actual, expected, rtol=1e-5, atol=1e-8, msg=""):
        """

        
        Args:





        """
        if isinstance(actual, torch.Tensor):
            actual = actual.detach().cpu()
        if isinstance(expected, torch.Tensor):
            expected = expected.detach().cpu()
        
        torch.testing.assert_close(actual, expected, rtol=rtol, atol=atol, msg=msg)
    
    return _assert


@pytest.fixture
def assert_cai_satisfied():

    def _assert(actual_cai, target_cai, tolerance=0.01):
        """
        Assert CAI constraint is satisfied

        Args:
            actual_cai: Actual CAI value
            target_cai: Target CAI value
            tolerance: Tolerance
        """
        assert actual_cai >= target_cai - tolerance, \
            f"CAI constraint not satisfied: {actual_cai:.4f} < {target_cai - tolerance:.4f}"

    return _assert




@pytest.fixture(autouse=True)
def reset_random_seeds():
    """Reset random seeds before each test to ensure reproducibility"""
    import random
    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(42)
        torch.cuda.manual_seed_all(42)


@pytest.fixture(autouse=True)
def cleanup_gpu_memory():

    yield
    if torch.cuda.is_available():
        torch.cuda.empty_cache()




def pytest_configure(config):
    """Add custom markers"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "gpu: marks tests that require GPU"
    )
    config.addinivalue_line(
        "markers", "integration: marks integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks unit tests"
    )
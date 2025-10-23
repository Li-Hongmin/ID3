# ID3 Framework for mRNA Sequence Design

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Research code for the paper: **"Gradient-based Optimization for mRNA Sequence Design with ID3 Framework"**

## Overview

This repository contains the complete implementation of the ID3 (Iterative Deep Learning-based Design) framework for optimizing mRNA sequences while maintaining biological constraints. The framework implements 12 optimization variants combining three constraint mechanisms with four optimization modes.

### Key Features

- **12 Optimization Variants**: 3 constraint mechanisms × 4 optimization modes
  - **Constraints**: Codon Profile Constraint (CPC), Amino Matching Softmax (AMS), Lagrangian Multiplier
  - **Modes**: Deterministic/Stochastic × Soft/Hard
- **DeepRaccess Integration**: RNA accessibility prediction for ribosome binding
- **CAI Optimization**: Codon Adaptation Index for translation efficiency
- **GPU Support**: CUDA acceleration for faster optimization

## Installation

### Prerequisites

- Python 3.8 or higher
- PyTorch 1.9 or higher
- CUDA-compatible GPU (recommended)

### Setup

**Quick Start (Recommended)**:
```bash
# Clone repository
git clone https://github.com/username/id3-framework.git
cd id3-framework

# Install dependencies
pip install -r requirements.txt

# Run demo - DeepRaccess will be set up automatically on first run
python demo.py
```

**That's it!** The first time you run the demo, it will automatically detect that DeepRaccess is missing and offer to set it up for you.

**Manual Setup**:
```bash
# Clone repository
git clone https://github.com/username/id3-framework.git
cd id3-framework

# Install dependencies
pip install -r requirements.txt

# Option 1: Automatic DeepRaccess setup
bash setup_deepraccess.sh

# Option 2: Manual DeepRaccess setup
git clone https://github.com/hmdlab/DeepRaccess.git
export DEEPRACCESS_PATH=$(pwd)/DeepRaccess
```

📖 **For detailed setup instructions, see [SETUP.md](SETUP.md)**

## Quick Start

### Simple Demo

```bash
# Interactive demo (default: 20 iterations)
python demo.py

# From FASTA file with more iterations
python demo.py --protein-file data/proteins/P04637.fasta --iterations 100

# Try different constraint mechanisms
python demo.py --constraint lagrangian --iterations 50       # Default, recommended
python demo.py --constraint amino_matching --iterations 50   # Softmax-based matching
python demo.py --constraint codon_profile --iterations 50    # Profile preservation

# Save optimized sequence
python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG --output result.fasta
```

### Production Experiments

For systematic experiments (research/paper reproduction):

```bash
# Quick test (5 iterations, 1 seed)
python run_unified_experiment.py --preset quick-test

# Full 12x12 experiments (1000 iterations, 12 seeds) - Accessibility only
python run_unified_experiment.py --preset full-12x12

# Full experiments with CAI optimization
python run_unified_experiment.py --preset full-12x12-cai-penalty

# Custom experiment
python run_unified_experiment.py \
    --proteins O15263,P04637 \
    --constraints lagrangian,ams,cpc \
    --variants 00,01,10,11 \
    --iterations 1000 \
    --seeds 12 \
    --enable-cai \
    --device cpu
```

Results saved to `results/` directory with detailed metrics and trajectories.

### What Each Tool Does

**`demo.py`** - Interactive demonstration
- ✅ Single protein optimization
- ✅ Visual progress feedback
- ✅ Good for testing and learning

**`run_unified_experiment.py`** - Systematic experiments
- ✅ Batch experiments (multiple proteins/constraints/variants)
- ✅ Multiple random seeds for statistical analysis
- ✅ Detailed result tracking and analysis
- ✅ Used for paper results

Both tools optimize:
- **Amino acid constraints** (3 mechanisms: Lagrangian, AMS, CPC)
- **CAI optimization** (Codon Adaptation Index)
- **RNA accessibility** (DeepRaccess prediction)

### Case Study Example

Run complete optimization with automatic visualization:

```bash
# Quick case study (20 iterations)
bash run_case_study.sh MSKGEELFT lagrangian 20

# Full case study (100 iterations)
bash run_case_study.sh MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG lagrangian 100

# Try different constraints
bash run_case_study.sh MSKGEELFT amino_matching 50
bash run_case_study.sh MSKGEELFT codon_profile 50
```

Outputs:
- `case_study_results/optimized_sequence.fasta` - Optimized mRNA
- `case_study_results/optimization_result.json` - Detailed data
- `case_study_results/*_convergence.png` - Convergence curves
- `case_study_results/*_sequence.png` - Sequence analysis

### Help
```bash
python demo.py --help                      # Demo options
python run_unified_experiment.py --help    # Experiment options
python visualize.py --help                 # Visualization options
```

## Repository Structure

```
id3-framework/
├── demo.py                      # CLI demo script
├── README.md                    # This file
├── LICENSE                      # CC BY-NC-SA 4.0 license
├── CITATION.cff                # Citation information
├── requirements.txt             # Python dependencies
│
├── id3/                        # Source code (49 files)
│   ├── constraints/            # Constraint mechanisms (12 files)
│   ├── optimizers/             # Optimization engines (10 files)
│   ├── cai/                    # CAI module (9 files)
│   └── utils/                  # Utility functions (18 files)
│
└── data/                        # Data files
    ├── proteins/               # Test protein sequences (12 files)
    ├── codon_references/       # CAI reference data (15 files)
    └── utr_templates/          # UTR templates (3 files)
```

## Usage

### Command Line (Recommended)

```bash
# Basic usage
python demo.py

# With CAI optimization
python demo.py --enable-cai

# From FASTA file
python demo.py --protein-file data/proteins/P04637.fasta --enable-cai

# Custom parameters
python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG \
               --constraint lagrangian \
               --enable-cai \
               --cai-target 0.9 \
               --iterations 10

# Save output
python demo.py --protein-file data/proteins/P0DTC9.fasta \
               --enable-cai \
               --output optimized.fasta
```

### Python API

```python
import sys
sys.path.insert(0, 'src')

from id3.constraints.lagrangian import LagrangianConstraint

# Create constraint (access-only)
constraint = LagrangianConstraint(
    protein_sequence,
    enable_cai=False
)

# Generate RNA sequence
result = constraint.forward(alpha=0.5, beta=0.5)
rna_seq = result['discrete_sequence']

# With CAI optimization
constraint_cai = LagrangianConstraint(
    protein_sequence,
    enable_cai=True,
    cai_target=0.8,
    cai_lambda=0.1
)

result = constraint_cai.forward(alpha=0.5, beta=0.5)
rna_seq = result['discrete_sequence']
cai_value = result['cai_metadata']['final_cai']
```

## Constraint Mechanisms

The ID3 framework provides 3 constraint mechanisms to ensure RNA sequences encode the correct amino acids. **All 3 mechanisms support joint optimization with DeepRaccess**.

### 1. Lagrangian Multiplier
- **Method**: Soft penalty-based optimization with adaptive λ
- **Formula**: `L = f_accessibility + λ·C_amino + λ_CAI·L_CAI`
- **Advantages**: Flexible penalty adjustment, stable optimization
- **Used in**: Both `demo.py` (default) and `demo_accessibility.py` (default)

### 2. Amino Matching Softmax (AMS)
- **Method**: Softmax-based amino acid probability matching
- **Advantages**: Differentiable, enforces constraints naturally
- **Used in**: Both `demo.py` and `demo_accessibility.py` (with `--constraint amino_matching`)

### 3. Codon Profile Constraint (CPC)
- **Method**: Maintains codon usage distribution from initial sequence
- **Advantages**: Preserves codon usage patterns
- **Used in**: Both `demo.py` and `demo_accessibility.py` (with `--constraint codon_profile`)

**Key Insight**: All constraint mechanisms output soft probability distributions (`rna_sequence`) that can be used for gradient-based optimization with DeepRaccess. The gradient flows through:
```
Constraint → Soft Probabilities → DeepRaccess → Accessibility Loss → Backprop
```

## Optimization Modes

- **det_soft**: Deterministic gradient descent with soft constraints
- **det_hard**: Deterministic gradient descent with hard constraints
- **sto_soft**: Stochastic sampling with soft constraints
- **sto_hard**: Stochastic sampling with hard constraints

## Demo Examples

### Example 1: Basic Access-Only
```bash
python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG
```

### Example 2: With CAI Optimization
```bash
python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG --enable-cai
```

### Example 3: From File with Custom Parameters
```bash
python demo.py --protein-file data/proteins/P04637.fasta \
               --enable-cai \
               --cai-target 0.9 \
               --cai-lambda 0.2 \
               --output result.fasta
```

### Example 4: Different Constraint Types
```bash
# Lagrangian Multiplier
python demo.py --constraint lagrangian

# Amino Matching Softmax
python demo.py --constraint amino_matching

# Codon Profile Constraint
python demo.py --constraint codon_profile
```

## Expected Performance

Production settings (1000 iterations on V100 GPU):

| Constraint | Mode | Accessibility Δ | Constraint Rate | Time |
|------------|------|-----------------|-----------------|------|
| Lagrangian | sto_hard | -0.65 kcal/mol | 100% | ~45 min |
| AMS | sto_soft | -0.85 kcal/mol | 100% | ~40 min |
| CPC | det_hard | -0.90 kcal/mol | 100% | ~35 min |

## Citation

If you use this code in your research, please cite:

```bibtex
@article{li2025id3,
  title={Gradient-based Optimization for mRNA Sequence Design with ID3 Framework},
  author={Li, Hongmin and Terai, Goro and Otagaki, Takumi and Asai, Kiyoshi},
  year={2025},
  note={In preparation}
}
```

*Note: DOI and journal information will be updated upon publication.*

## License

This work is licensed under [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](LICENSE).

**Academic Use**: ✅ Freely permitted
**Commercial Use**: ❌ Prohibited without permission
**Attribution**: ✅ Required in all publications

For commercial licensing inquiries: lihongmin@edu.k.u-tokyo.ac.jp

See [LICENSE-SUMMARY.md](LICENSE-SUMMARY.md) for detailed terms.

## Contact

- **Research Questions**: lihongmin@edu.k.u-tokyo.ac.jp
- **Bug Reports**: GitHub Issues
- **Commercial Licensing**: lihongmin@edu.k.u-tokyo.ac.jp

## Acknowledgments

- DeepRaccess model: [https://github.com/hmdlab/DeepRaccess](https://github.com/hmdlab/DeepRaccess)
- University of Tokyo

---

**Version**: 1.0.0
**Last Updated**: January 15, 2025
**Maintained by**: University of Tokyo

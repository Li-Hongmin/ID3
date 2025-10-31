# ID3 Framework for mRNA Sequence Design

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

Research code for the paper: **"Gradient-based Optimization for mRNA Sequence Design"**

**Publication**: Li Hongmin, Goro Terai, Takumi Otagaki, Kiyoshi Asai. *bioRxiv* 2025.10.22.683691; doi: https://doi.org/10.1101/2025.10.22.683691

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
- CUDA-compatible GPU (optional, CPU supported)

### Quick Setup

```bash
# Clone repository
git clone https://github.com/Li-Hongmin/ID3.git
cd ID3

# Install dependencies
pip install -r requirements.txt

# Run demo - DeepRaccess will be set up automatically
bash run_demo.sh
```

**That's it!** The demo automatically detects and installs DeepRaccess on first run.

## Quick Start

### One-Click Case Study Demo

```bash
# Default: O15263 protein, 1000 iterations
bash run_demo.sh

# Custom protein and iterations
bash run_demo.sh O15263 1000     # Protein ID, iterations
bash run_demo.sh P04637 500      # Different protein

# Results saved to examples/demo_<timestamp>/
# Includes: optimized sequence, trajectory data, and visualizations
```

The demo automatically:
- ✅ Checks and installs DeepRaccess if needed
- ✅ Runs mRNA optimization with AMS constraint
- ✅ Generates publication-quality evolution figures
- ✅ Saves all results to `examples/` directory

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

**`run_demo.sh`** - Quick case study demo
- ✅ One-click complete workflow
- ✅ Automatic DeepRaccess setup
- ✅ Single protein optimization with visualization
- ✅ Results saved to `examples/` directory
- ✅ Perfect for quick demonstrations

**`run_unified_experiment.py`** - Research-grade experiments
- ✅ Batch experiments (multiple proteins/constraints/variants)
- ✅ Multiple random seeds for statistical analysis
- ✅ 12 optimization variants (3 constraints × 4 modes)
- ✅ Detailed result tracking and analysis
- ✅ Used for paper results

Both tools optimize:
- **Amino acid constraints** (3 mechanisms: Lagrangian, AMS, CPC)
- **CAI optimization** (Codon Adaptation Index)
- **RNA accessibility** (DeepRaccess prediction)


## Repository Structure

```
ID3/
├── run_demo.sh                  # One-click case study demo
├── demo.py                      # Main demo script
├── run_unified_experiment.py    # Research experiment framework
├── README.md                    # This file
├── requirements.txt             # Python dependencies
│
├── scripts/                     # Auxiliary scripts
│   ├── evolution_figure.py     # Visualization generator
│   ├── setup_deepraccess.sh    # DeepRaccess installer
│   └── README.md               # Scripts documentation
│
├── id3/                        # Source code
│   ├── constraints/            # Constraint mechanisms (AMS, CPC, Lagrangian)
│   ├── optimizers/             # Optimization engines
│   ├── cai/                    # CAI module
│   └── utils/                  # Utility functions
│
├── data/                        # Data files
│   ├── proteins/               # Test protein sequences (.fasta.txt)
│   ├── codon_references/       # CAI reference data
│   └── utr_templates/          # UTR templates
│
└── examples/                    # Demo results (with visualizations)
    └── demo_20251031_233130/   # Example: 1000-iter optimization
```

## Usage

### Command Line (Recommended)

```bash
# Quick demo (default: 1000 iterations)
bash run_demo.sh

# Custom protein and iterations
bash run_demo.sh O15263 1000
bash run_demo.sh P04637 500

# Research-grade experiments
python run_unified_experiment.py --preset quick-test
python run_unified_experiment.py --preset full-12x12
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
- **Usage**: `demo.py --constraint lagrangian` (default)

### 2. Amino Matching Softmax (AMS)
- **Method**: Softmax-based amino acid probability matching
- **Advantages**: Differentiable, enforces constraints naturally
- **Usage**: `demo.py --constraint amino_matching` or `run_demo.sh` (default)

### 3. Codon Profile Constraint (CPC)
- **Method**: Maintains codon usage distribution from initial sequence
- **Advantages**: Preserves codon usage patterns
- **Usage**: `demo.py --constraint codon_profile`

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
python demo.py --protein-file data/proteins/P04637.fasta.txt \
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

## Citation

If you use this code in your research, please cite:

```bibtex
@article{li2025gradient,
  title={Gradient-based Optimization for mRNA Sequence Design},
  author={Li, Hongmin and Terai, Goro and Otagaki, Takumi and Asai, Kiyoshi},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.10.22.683691},
  url={https://doi.org/10.1101/2025.10.22.683691}
}
```

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

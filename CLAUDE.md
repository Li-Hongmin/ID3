# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ID3 (Iterative Deep Learning-based Design) framework for gradient-based mRNA sequence optimization. The framework optimizes RNA sequences for both accessibility (via DeepRaccess) and translation efficiency (via CAI - Codon Adaptation Index) while maintaining amino acid constraints.

**Key concept**: Optimizes mRNA sequences by balancing RNA accessibility (ribosome binding) and codon usage (translation efficiency) through differentiable constraint mechanisms.

## Development Commands

### Setup

**Automatic Setup (Recommended)**:
```bash
# Install Python dependencies
pip install -r requirements.txt

# Run demo - it will automatically detect and offer to set up DeepRaccess
python demo.py
```

The first time you run `demo.py`, it will:
1. Detect that DeepRaccess is not installed
2. Offer to automatically clone and configure it
3. Clone from https://github.com/hmdlab/DeepRaccess
4. Verify the installation

**Manual Setup**:
```bash
# Install dependencies
pip install -r requirements.txt

# Option 1: Use the setup script
bash setup_deepraccess.sh

# Option 2: Manual clone
git clone https://github.com/hmdlab/DeepRaccess.git
export DEEPRACCESS_PATH=$(pwd)/DeepRaccess
```

**Skip DeepRaccess Check**:
```bash
# For development/testing without DeepRaccess
export SKIP_DEEPRACCESS_CHECK=1
python demo.py
```

### Demo vs Production Experiments

**Interactive Demo** (`demo.py`):
```bash
# Quick testing and visualization
python demo.py --protein MSKGEELFT --iterations 50
python demo.py --constraint amino_matching --iterations 100
```

**Systematic Experiments** (`run_unified_experiment.py`):
```bash
# Quick test (5 iterations)
python run_unified_experiment.py --preset quick-test --device cpu

# Full 12x12 matrix (1000 iterations, 12 seeds)
python run_unified_experiment.py --preset full-12x12 --device cpu

# With CAI optimization
python run_unified_experiment.py --preset full-12x12-cai-penalty --device cpu

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

**Key differences**:
- `demo.py`: Single protein, interactive, visual feedback
- `run_unified_experiment.py`: Batch mode, multiple seeds, results in `results/` directory

### Common Parameters

**demo.py**:
- `--protein` / `--protein-file`: Protein sequence
- `--constraint`: lagrangian/amino_matching/codon_profile
- `--iterations`: Number of optimization steps (default: 20)
- `--cai-target`: Target CAI value (default: 0.8)
- `--cai-weight`: CAI weight in loss (default: 0.1)
- `--utr5-file` / `--utr3-file`: Custom UTR sequences
- `--output`: Save result to file

**run_unified_experiment.py**:
- `--proteins`: Comma-separated protein IDs
- `--constraints`: Comma-separated constraint types
- `--variants`: 00/01/10/11 (det/sto × soft/hard)
- `--iterations`: Default 1000 for production
- `--seeds`: Number of random seeds (default: 12)
- `--enable-cai`: Enable CAI optimization
- `--device`: cuda/cpu

## Architecture

### Core Components

**1. Constraint Mechanisms** (`id3/constraints/`)

The framework provides 3 constraint types to ensure RNA encodes correct amino acids:

**All 3 constraint types support gradient-based joint optimization with DeepRaccess**

**A. Lagrangian Multiplier** (`lagrangian.py`) - **Recommended**
- Soft penalty-based optimization: `L = f_accessibility + λ·C_amino + λ_CAI·L_CAI`
- Adaptive λ for constraint strength
- Most stable for long optimization runs
- **Used in**: Both demos (default)

**B. Amino Matching Softmax** (`amino_matching.py`)
- Softmax-based amino acid probability matching
- Differentiable constraint enforcement
- Natural probability distribution over valid codons
- **Used in**: Both demos (`--constraint amino_matching`)

**C. Codon Profile Constraint** (`codon_profile.py`)
- Maintains codon usage distribution from initial sequence
- Preserves codon frequency patterns
- Good for sequences with known optimal codon usage
- **Used in**: Both demos (`--constraint codon_profile`)

**Key mechanism**: All constraints output soft probability distributions (`rna_sequence`) that enable gradient flow:
```
Constraint → Soft Probs → DeepRaccess → Accessibility Loss → Backprop
```

**Base class**: `BaseConstraint` - Unified parent with CAI integration for all constraint types

**2. CAI Module** (`id3/cai/`)
- `unified_calculator.py`: Main CAI calculation engine
- `incremental_cai_calculator.py`: Efficient incremental CAI updates
- `discrete_binary_search.py`: Finds optimal discrete sequences
- `unified_enhancer.py`: CAI enhancement operations
- `validator.py`: CAI validation and verification

**3. Optimization Flow** (Lagrangian example)
```
Θ (parameters) → Π (transform) → P (probabilities) → Ψ (discretize) → S (sequence)
→ f_model (accessibility) + λC (constraint penalty) → L (loss) → ∇ (gradient) → Θ^(t+1)
```

**4. Utilities** (`id3/utils/`)
- `functions.py`: Codon/amino acid mappings, token conversions
- `sequence_utils.py`: RNA/protein sequence operations
- `global_cache.py`: Performance optimization via caching
- `constraint_satisfied_argmax.py`: Constraint-preserving argmax for discretization

### Key Design Patterns

**1. Unified CAI Integration**
All constraint classes inherit from `BaseConstraint` which includes `AdaptiveLambdaCAIMixin`. This provides:
- Unified CAI component initialization
- Common loss computation interface
- Shared precomputation and caching mechanisms

**2. Four Optimization Modes**
Controlled by `alpha` (stochasticity) and `beta` (output type):
- `det_soft`: alpha=0, beta=0 - Deterministic gradient, soft outputs
- `det_hard`: alpha=0, beta=1 - Deterministic gradient, hard outputs
- `sto_soft`: alpha=1, beta=0 - Stochastic sampling, soft outputs
- `sto_hard`: alpha=1, beta=1 - Stochastic sampling, hard outputs

**3. Constraint Satisfaction**
Each constraint type ensures the RNA sequence encodes the correct amino acid sequence:
- **Lagrangian**: Penalty-based approach with `ConstraintPenalty.compute_vectorized()`
- **AMS**: Softmax matching with `amino_acid_to_codon_matrix` lookup
- **CPC**: Profile preservation via codon distribution matching

**4. Gumbel Noise for Exploration**
When `alpha > 0`, Gumbel noise is added to logits before softmax:
```python
gumbel_noise = -torch.log(-torch.log(torch.rand_like(theta) + 1e-10) + 1e-10)
theta_noisy = theta + alpha * gumbel_noise
```

### Critical Implementation Details

**Codon Encoding**:
- RNA sequences represented as one-hot vectors (4 nucleotides × sequence length)
- Codons are 3-nucleotide groups → amino acids (genetic code mapping in `utils/functions.py`)
- Each amino acid position must select from valid codons (enforced by `valid_codon_mask`)

**CAI Calculation**:
- Uses species-specific codon usage tables (`data/codon_references/`)
- Differentiable approximation via weighted sums over codon probabilities
- Incremental updates for efficiency during optimization

**Gradient Flow**:
- Only Lagrangian supports iterative optimization via `theta` parameter
- AMS and CPC use single forward pass (no learnable parameters)
- CAI loss integrated via `result['cai_loss']` when `enable_cai=True`

## Working with This Codebase

### Adding New Constraint Mechanisms
1. Inherit from `BaseConstraint` in `id3/constraints/base.py`
2. Implement `forward()` method returning dict with:
   - `rna_sequence`: Soft probabilities
   - `discrete_sequence`: Hard sequence string
   - `constraint_loss`: Penalty term (if applicable)
   - `cai_metadata`: CAI info (handled by parent class)
3. Implement amino acid verification in `verify_amino_acid_constraint()`

### Modifying CAI Behavior
- CAI configuration in `id3/cai/config.py`
- Species-specific codon tables in `data/codon_references/`
- Core calculation logic in `unified_calculator.py`

### GPU/CPU Handling
- Device auto-detection in constraint classes
- All tensors automatically moved to correct device
- CUDA recommended for sequences >500 amino acids

### DeepRaccess Integration
- External dependency for RNA accessibility prediction
- **Automatic setup**: `demo.py` will offer to clone and configure on first run
- **Manual setup**: Run `bash setup_deepraccess.sh` or clone manually
- Optional environment variable: `DEEPRACCESS_PATH` (auto-detected if in project root)
- Model called during forward pass for accessibility scoring
- Wrapper location: `id3/utils/deepraccess_wrapper.py`
- Pre-trained models: Should be in `DeepRaccess/path/*.pth`

**Key Implementation Details**:
- Soft embedding for gradient flow: Probability-weighted embedding matrices
- Index mapping: ID3 [A,C,G,U] ↔ DeepRaccess [pad,mask,A,U,G,C]
- Sliding window: 440nt windows with 110nt overlap for long sequences
- ATG window: Extracts -19 to +15 position (35nt) for ribosome binding optimization

## Data Files

- `data/proteins/`: Test protein sequences (FASTA format)
- `data/codon_references/`: Species-specific codon usage tables for CAI
- `data/utr_templates/`: UTR (untranslated region) templates

## Known Considerations

- Production experiments typically use 1000 iterations (~35-45 min on V100 GPU)
- Demo defaults to 10 iterations for quick testing
- CAI optimization adds ~20-30% computational overhead
- Constraint satisfaction rate should be 100% (verify with `verify_amino_acid_constraint()`)
- /Users/lihongmin/Research/ID3_DeepRaccess_CAI_Paper 这个是原来的我们的项目 现在要整理出一个 github公开版
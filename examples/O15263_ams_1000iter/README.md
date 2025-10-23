# O15263 Case Study - AMS Constraint (1000 Iterations)

## Configuration

- **Protein**: O15263 (64 amino acids)
- **Constraint**: Amino Matching Softmax (AMS)
- **Iterations**: 1000
- **Device**: CPU

## Results

### Optimization Metrics

- **Initial accessibility**: 1.746 kcal/mol
- **Best accessibility**: 0.908 kcal/mol
- **Final accessibility**: 2.146 kcal/mol
- **Improvement**: 48.0% (initial → best)

### CAI Metrics

- **Final CAI**: 0.801
- **Target CAI**: 0.8
- **Status**: ✅ Target achieved

### Sequence

Final optimized mRNA sequence (192 nucleotides):
```
AUGCGUGUUCUGUAUCUGCUGUUUUCUUUUCUGUUCAUCUUUCUGAUGCCGCUGCCGGGC
GUUUUCGGUGGUAUCGGUGACCCAGUCACCUGCCUGAAAUCUGGUGCUAUUUGCCACCCG
GUAUUUUGCCCACGUCGUUACAAACAAAUCGGUACUUGUGGUCUGCCGGGUACCAAAUGUUGCAAAAAACCG
```

See: `optimized_sequence.fasta`

## Visualization

The 3-panel evolution figure shows:

### Top Panel: Nucleotide Probability Evolution
- RGB color-coded heatmap (A=Red, C=Green, G=Blue, U=Black)
- First 45 positions of the mRNA sequence
- Codon boundaries marked every 3 nucleotides
- Final sequence displayed on right axis

### Middle Panel: Convergence Curve
- Accessibility score evolution over 1000 iterations
- Best point marked at iteration 999 (0.908 kcal/mol)
- 48% improvement from initial state

### Bottom Panel: AU Content Evolution
- AU ratio changes during optimization
- Tracks compositional evolution

## Files

- `optimized_sequence.fasta` - Final mRNA sequence with metrics
- `O15263_amino_matching_figure.png` - Evolution figure (355KB, 300 DPI)
- `O15263_amino_matching_figure.pdf` - Vector format (57KB)
- `README.md` - This file

## Reproduction

To reproduce this case study:

```bash
bash run_case_study.sh O15263 ams 1000
```

Results will be saved to `case_study_results/` directory.

---

**Generated**: 2025-10-23
**Framework**: ID3 v1.0.0
**Constraint**: Amino Matching Softmax (AMS)
**Optimization time**: ~8 minutes on CPU

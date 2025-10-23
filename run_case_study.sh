#!/bin/bash
# ID3 Framework - Case Study Example
#
# This script runs a complete case study: optimization + visualization
# Demonstrates the full ID3 workflow with automatic result generation

set -e

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║      ID3 Framework - Case Study Example                       ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
PROTEIN="${1:-MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG}"
CONSTRAINT="${2:-lagrangian}"
ITERATIONS="${3:-100}"
OUTPUT_DIR="case_study_results"

mkdir -p "$OUTPUT_DIR"

echo "Configuration:"
echo "  Protein: $PROTEIN"
echo "  Constraint: $CONSTRAINT"
echo "  Iterations: $ITERATIONS"
echo "  Output: $OUTPUT_DIR/"
echo ""

# Step 1: Run optimization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 1: Running ID3 optimization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

python demo.py \
    --protein "$PROTEIN" \
    --constraint "$CONSTRAINT" \
    --iterations "$ITERATIONS" \
    --output "$OUTPUT_DIR/optimized_sequence.fasta" \
    --save-result "$OUTPUT_DIR/optimization_result.json"

echo ""

# Step 2: Generate visualization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 2: Generating visualization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

python visualize_results.py "$OUTPUT_DIR/optimization_result.json"

# Move generated figures to output directory
mv *.png "$OUTPUT_DIR/" 2>/dev/null || true

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║      ✅ Case Study Complete!                                   ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo "Results saved to: $OUTPUT_DIR/"
echo ""
echo "Files generated:"
echo "  📄 optimized_sequence.fasta  - Optimized mRNA sequence"
echo "  📊 optimization_result.json  - Detailed optimization data"
echo "  📈 *_convergence.png         - Convergence curves"
echo "  🧬 *_sequence.png            - Sequence analysis"
echo ""
echo "To view figures:"
echo "  open $OUTPUT_DIR/*.png"
echo ""

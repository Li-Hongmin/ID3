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
PROTEIN_ARG="${1:-MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG}"
CONSTRAINT="${2:-lagrangian}"
ITERATIONS="${3:-100}"
OUTPUT_DIR="case_study_results"

mkdir -p "$OUTPUT_DIR"

# Check if protein is a file or sequence
if [ -f "$PROTEIN_ARG" ]; then
    PROTEIN_FILE="$PROTEIN_ARG"
    PROTEIN_NAME=$(basename "$PROTEIN_FILE" .fasta.txt)
    echo "Configuration:"
    echo "  Protein file: $PROTEIN_FILE"
    echo "  Constraint: $CONSTRAINT"
    echo "  Iterations: $ITERATIONS"
    echo "  Output: $OUTPUT_DIR/"
elif [ -f "data/proteins/$PROTEIN_ARG.fasta.txt" ]; then
    PROTEIN_FILE="data/proteins/$PROTEIN_ARG.fasta.txt"
    PROTEIN_NAME="$PROTEIN_ARG"
    echo "Configuration:"
    echo "  Protein: $PROTEIN_NAME (from data/proteins/)"
    echo "  Constraint: $CONSTRAINT"
    echo "  Iterations: $ITERATIONS"
    echo "  Output: $OUTPUT_DIR/"
else
    PROTEIN_FILE=""
    PROTEIN_NAME="custom"
    echo "Configuration:"
    echo "  Protein sequence: $PROTEIN_ARG"
    echo "  Constraint: $CONSTRAINT"
    echo "  Iterations: $ITERATIONS"
    echo "  Output: $OUTPUT_DIR/"
fi

echo ""

# Step 1: Run optimization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 1: Running ID3 optimization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -n "$PROTEIN_FILE" ]; then
    python demo.py \
        --protein-file "$PROTEIN_FILE" \
        --constraint "$CONSTRAINT" \
        --iterations "$ITERATIONS" \
        --output "$OUTPUT_DIR/optimized_sequence.fasta" \
        --save-result "$OUTPUT_DIR/optimization_result.json"
else
    python demo.py \
        --protein "$PROTEIN_ARG" \
        --constraint "$CONSTRAINT" \
        --iterations "$ITERATIONS" \
        --output "$OUTPUT_DIR/optimized_sequence.fasta" \
        --save-result "$OUTPUT_DIR/optimization_result.json"
fi

echo ""

# Step 2: Generate visualization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 2: Generating visualization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Generate paper-quality 3-panel figure
python visualize_case_study.py --json "$OUTPUT_DIR/optimization_result.json" --output "$OUTPUT_DIR/${PROTEIN_NAME}_${CONSTRAINT}_figure.png"

# Move any remaining figures to output directory
mv *.png *.pdf "$OUTPUT_DIR/" 2>/dev/null || true

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

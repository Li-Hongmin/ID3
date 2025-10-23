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
CONSTRAINT_INPUT="${2:-lagrangian}"
ITERATIONS_ARG="${3}"  # Optional, will auto-calculate if not provided
OUTPUT_DIR="case_study_results"

# Map constraint abbreviations to full names
case "$CONSTRAINT_INPUT" in
    lag|lagrangian)
        CONSTRAINT="lagrangian"
        CONSTRAINT_SHORT="lag"
        ;;
    ams|amino_matching)
        CONSTRAINT="amino_matching"
        CONSTRAINT_SHORT="ams"
        ;;
    cpc|codon_profile)
        CONSTRAINT="codon_profile"
        CONSTRAINT_SHORT="cpc"
        ;;
    *)
        echo "❌ Error: Unknown constraint '$CONSTRAINT_INPUT'"
        echo "Valid options: lagrangian (lag), amino_matching (ams), codon_profile (cpc)"
        exit 1
        ;;
esac

mkdir -p "$OUTPUT_DIR"

# Check if protein is a file or sequence and get protein length
if [ -f "$PROTEIN_ARG" ]; then
    PROTEIN_FILE="$PROTEIN_ARG"
    PROTEIN_NAME=$(basename "$PROTEIN_FILE" .fasta.txt)
    PROTEIN_SEQ=$(grep -v "^>" "$PROTEIN_FILE" | tr -d '\n')
    PROTEIN_LENGTH=${#PROTEIN_SEQ}
elif [ -f "data/proteins/$PROTEIN_ARG.fasta.txt" ]; then
    PROTEIN_FILE="data/proteins/$PROTEIN_ARG.fasta.txt"
    PROTEIN_NAME="$PROTEIN_ARG"
    PROTEIN_SEQ=$(grep -v "^>" "$PROTEIN_FILE" | tr -d '\n')
    PROTEIN_LENGTH=${#PROTEIN_SEQ}
else
    PROTEIN_FILE=""
    PROTEIN_NAME="custom"
    PROTEIN_SEQ="$PROTEIN_ARG"
    PROTEIN_LENGTH=${#PROTEIN_SEQ}
fi

# Auto-calculate iterations if not provided
if [ -z "$ITERATIONS_ARG" ]; then
    # Adaptive iteration calculation based on protein length and constraint
    if [ $PROTEIN_LENGTH -le 20 ]; then
        BASE_ITERATIONS=50
    elif [ $PROTEIN_LENGTH -le 50 ]; then
        BASE_ITERATIONS=100
    elif [ $PROTEIN_LENGTH -le 100 ]; then
        BASE_ITERATIONS=200
    elif [ $PROTEIN_LENGTH -le 200 ]; then
        BASE_ITERATIONS=300
    else
        BASE_ITERATIONS=500
    fi

    # Adjust for constraint type (AMS/CPC typically converge faster)
    if [ "$CONSTRAINT" = "amino_matching" ] || [ "$CONSTRAINT" = "codon_profile" ]; then
        ITERATIONS=$((BASE_ITERATIONS * 3 / 4))
    else
        ITERATIONS=$BASE_ITERATIONS
    fi

    AUTO_ITERATIONS="(auto-calculated)"
else
    ITERATIONS=$ITERATIONS_ARG
    AUTO_ITERATIONS=""
fi

# Display configuration
echo "Configuration:"
if [ -n "$PROTEIN_FILE" ]; then
    echo "  Protein: $PROTEIN_NAME (${PROTEIN_LENGTH} amino acids)"
else
    echo "  Protein sequence: ${PROTEIN_SEQ:0:30}... (${PROTEIN_LENGTH} amino acids)"
fi
echo "  Constraint: $CONSTRAINT ($CONSTRAINT_SHORT)"
echo "  Iterations: $ITERATIONS $AUTO_ITERATIONS"
echo "  Output: $OUTPUT_DIR/"

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
python visualize.py --json "$OUTPUT_DIR/optimization_result.json" --output "$OUTPUT_DIR/${PROTEIN_NAME}_${CONSTRAINT}_figure.png"

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

#!/bin/bash
# ID3 Framework - One-Click Case Study Demo
# Automatically handles DeepRaccess setup and generates complete case study results

set -e

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║      ID3 Framework - Case Study Demo                          ║"
echo "║      One-Click RNA Optimization Example                       ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
PROTEIN="${1:-O15263}"           # Default: O15263 protein (MeCP2)
ITERATIONS="${2:-200}"           # Default: 200 iterations
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTPUT_DIR="$SCRIPT_DIR/examples/demo_$TIMESTAMP"

# Color codes
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Configuration:${NC}"
echo "  Protein: $PROTEIN"
echo "  Iterations: $ITERATIONS"
echo "  Output: $OUTPUT_DIR"
echo ""

# Step 0: Check and setup DeepRaccess if needed
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 0: Checking DeepRaccess dependency..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ ! -d "$SCRIPT_DIR/DeepRaccess" ]; then
    echo -e "${YELLOW}⚠️  DeepRaccess not found. Installing automatically...${NC}"
    echo ""
    bash "$SCRIPT_DIR/scripts/setup_deepraccess.sh"
    echo ""
else
    echo -e "${GREEN}✓ DeepRaccess found${NC}"
    echo ""
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Step 1: Run optimization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 1: Running ID3 optimization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Check if protein is a file or ID
if [ -f "data/proteins/${PROTEIN}.fasta.txt" ]; then
    PROTEIN_ARG="--protein-file data/proteins/${PROTEIN}.fasta.txt"
else
    PROTEIN_ARG="--protein $PROTEIN"
fi

python "$SCRIPT_DIR/demo.py" \
    $PROTEIN_ARG \
    --constraint amino_matching \
    --iterations "$ITERATIONS" \
    --alpha 0.5 \
    --beta 0.0 \
    --output "$OUTPUT_DIR/optimized_sequence.fasta" \
    --save-result "$OUTPUT_DIR/optimization_result.json"

echo ""

# Step 2: Generate visualization
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Step 2: Generating visualization..."
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

python "$SCRIPT_DIR/scripts/evolution_figure.py" \
    --json "$OUTPUT_DIR/optimization_result.json" \
    --output "$OUTPUT_DIR/${PROTEIN}_ams_figure.png"

echo ""

# Generate README for this case study
cat > "$OUTPUT_DIR/README.md" << EOF
# ID3 Case Study: $PROTEIN

Generated: $(date)

## Configuration
- Protein: $PROTEIN
- Constraint: Amino Acid Matching (AMS)
- Iterations: $ITERATIONS
- Device: CPU
- Alpha: 0.5 (soft probabilities)
- Beta: 0.0 (deterministic)

## Files
- \`optimized_sequence.fasta\` - Optimized mRNA sequence
- \`optimization_result.json\` - Complete optimization trajectory
- \`${PROTEIN}_ams_figure.png\` - Evolution visualization (PNG)
- \`${PROTEIN}_ams_figure.pdf\` - Evolution visualization (PDF)

## How to Use
View the figures:
\`\`\`bash
open ${PROTEIN}_ams_figure.png
\`\`\`

## Results
Check \`optimization_result.json\` for detailed metrics:
- Accessibility scores over iterations
- CAI (Codon Adaptation Index) evolution
- Nucleotide probability evolution (smooth gradients)
- Constraint satisfaction metrics
EOF

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                                                                ║"
echo "║      ✅ Case Study Complete!                                   ║"
echo "║                                                                ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""
echo -e "${GREEN}Results saved to:${NC} $OUTPUT_DIR/"
echo ""
echo "Files generated:"
echo "  📄 optimized_sequence.fasta  - Optimized mRNA sequence"
echo "  📊 optimization_result.json  - Detailed optimization data"
echo "  📈 ${PROTEIN}_ams_figure.png - Evolution visualization (PNG)"
echo "  📈 ${PROTEIN}_ams_figure.pdf - Evolution visualization (PDF)"
echo "  📖 README.md                - Case study documentation"
echo ""
echo -e "${GREEN}To view results:${NC}"
echo "  open $OUTPUT_DIR/${PROTEIN}_ams_figure.png"
echo ""
echo -e "${GREEN}To commit to git:${NC}"
echo "  git add $OUTPUT_DIR/"
echo "  git commit -m \"Add case study: $PROTEIN ($TIMESTAMP)\""
echo ""

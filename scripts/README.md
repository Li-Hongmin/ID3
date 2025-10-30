# ID3 Scripts Directory

This directory contains auxiliary scripts and tools for the ID3 framework.

## Scripts Overview


### `evolution_figure.py`
**Purpose**: Generate paper-quality 3-panel evolution figures

**Usage**:
```bash
python scripts/evolution_figure.py \
    --json optimization_result.json \
    --output figure.png
```

**Features**:
- **Panel 1**: Nucleotide probability evolution heatmap (RGB color mixing)
- **Panel 2**: Accessibility convergence curve
- **Panel 3**: AU content evolution

**Outputs**:
- PNG format (300 DPI, publication-ready)
- PDF format (vector graphics)

---

### `setup_deepraccess.sh`
**Purpose**: Automatic DeepRaccess dependency installation

**Usage**:
```bash
bash scripts/setup_deepraccess.sh
```

**Actions**:
1. Checks if DeepRaccess directory exists
2. If not, clones from https://github.com/hmdlab/DeepRaccess.git
3. Verifies essential files and models
4. Provides setup status and next steps

**Note**: This script is automatically called by `run_demo.sh` if needed.

---

## Integration with Main Scripts

### Called by `run_demo.sh`:
```bash
# Step 0: Setup DeepRaccess
bash scripts/setup_deepraccess.sh

# Step 1: Run optimization
python scripts/demo_case_study.py --protein O15263 --output examples/demo_...

# Step 2: Generate visualization
python scripts/evolution_figure.py --json ... --output ...
```

### Standalone Use:
All scripts can be used independently for custom workflows.

---

## Development Notes

### Adding New Scripts
1. Place script in `scripts/` directory
2. Add executable permissions: `chmod +x scripts/your_script.sh`
3. Update this README with usage instructions
4. Reference from main scripts if needed

### Dependencies
- Python 3.8+
- PyTorch (for optimization scripts)
- matplotlib, numpy (for visualization)
- DeepRaccess (for accessibility prediction)

See `requirements.txt` in root directory for complete list.

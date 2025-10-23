# ID3 Framework Setup Guide

This guide explains how to set up the ID3 framework, including automatic DeepRaccess configuration.

## Quick Start (Recommended)

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Run demo - DeepRaccess will be set up automatically
python demo.py
```

The first time you run `demo.py`, it will automatically:
- Detect that DeepRaccess is missing
- Ask if you want to set it up automatically
- Clone and configure DeepRaccess if you agree
- Verify the installation

## Manual Setup Options

### Option 1: Use the Automatic Setup Script

```bash
bash setup_deepraccess.sh
```

This script will:
- Check if DeepRaccess is already installed
- Clone from https://github.com/hmdlab/DeepRaccess
- Verify the installation
- Check for pre-trained models

### Option 2: Manual Clone

```bash
git clone https://github.com/hmdlab/DeepRaccess.git
```

The framework will automatically detect DeepRaccess in the project root directory.

## Verifying Installation

Run the test script to verify your setup:

```bash
python test_deepraccess_setup.py
```

This will check:
- ✓ DeepRaccess directory exists
- ✓ Required files are present (mymodel.py)
- ✓ Pre-trained models (if available)
- ✓ Import functionality

## DeepRaccess Pre-trained Models

The DeepRaccess repository may include pre-trained models in `DeepRaccess/path/`:
- `FCN_structured.pth` - Recommended model for structured RNA
- `FCN_uniform.pth` - Alternative model

If models are not included, the framework will use random initialization (for testing only).

## Skipping DeepRaccess (Development Mode)

For development or testing without DeepRaccess:

```bash
export SKIP_DEEPRACCESS_CHECK=1
python demo.py
```

**Note**: This will limit functionality to CAI optimization only.

## Troubleshooting

### DeepRaccess Clone Fails

If `git clone` fails:
1. Check your internet connection
2. Verify GitHub access
3. Try manual download from https://github.com/hmdlab/DeepRaccess

### Import Errors

If you see "No module named 'DeepRaccess'":
1. Verify DeepRaccess is in the project root directory
2. Check that `DeepRaccess/mymodel.py` exists
3. Restart Python session if you just installed DeepRaccess

### Missing Pre-trained Models

If no `.pth` files are found:
1. Check `DeepRaccess/path/` directory
2. Download models from DeepRaccess repository
3. Or use random initialization (for testing only)

## Advanced Configuration

### Custom DeepRaccess Location

If DeepRaccess is installed elsewhere:

```bash
export DEEPRACCESS_PATH=/path/to/DeepRaccess
```

### Custom Model Path

In Python code:

```python
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper

wrapper = DeepRaccessID3Wrapper(
    deepraccess_model_path='/path/to/model.pth',
    device='cuda'
)
```

## What Gets Installed

When you run automatic setup:

```
ID3-github/
├── DeepRaccess/              # Cloned repository
│   ├── mymodel.py           # DeepRaccess model definition
│   ├── path/                # Pre-trained models (if available)
│   │   ├── FCN_structured.pth
│   │   └── FCN_uniform.pth
│   └── ...
├── id3/                     # ID3 framework
│   └── utils/
│       └── deepraccess_wrapper.py  # Integration layer
└── demo.py                  # Entry point with auto-setup
```

## Dependencies

Required Python packages (from `requirements.txt`):
- numpy>=1.20.0
- torch>=1.9.0
- pyyaml>=5.4.0
- biopython>=1.79
- pandas>=1.3.0
- matplotlib>=3.4.0
- seaborn>=0.11.0
- tqdm>=4.62.0
- scikit-learn>=0.24.0
- scipy>=1.7.0

## Support

For issues with:
- **ID3 Framework**: Open issue on this repository
- **DeepRaccess**: Visit https://github.com/hmdlab/DeepRaccess

## License

- ID3 Framework: CC BY-NC-SA 4.0
- DeepRaccess: See DeepRaccess repository for license

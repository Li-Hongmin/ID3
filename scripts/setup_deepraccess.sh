#!/bin/bash
# DeepRaccess Automatic Setup Script for ID3 Framework
# This script automatically clones and configures DeepRaccess dependency

set -e  # Exit on error

echo "========================================"
echo "DeepRaccess Setup for ID3 Framework"
echo "========================================"
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments for non-interactive mode
NON_INTERACTIVE=false
if [[ "$1" == "-y" || "$1" == "--yes" ]]; then
    NON_INTERACTIVE=true
fi

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Install DeepRaccess to project root (parent of scripts directory)
PROJECT_ROOT="$SCRIPT_DIR/.."
DEEPRACCESS_DIR="$PROJECT_ROOT/DeepRaccess"

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if DeepRaccess already exists
if [ -d "$DEEPRACCESS_DIR" ]; then
    print_warning "DeepRaccess directory already exists at: $DEEPRACCESS_DIR"

    if [ "$NON_INTERACTIVE" = true ]; then
        # Non-interactive mode: verify and exit
        print_info "Non-interactive mode: Keeping existing installation"
        print_info "Verifying installation..."
    else
        # Interactive mode: ask user
        read -p "Do you want to remove and re-install? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_info "Removing existing DeepRaccess directory..."
            rm -rf "$DEEPRACCESS_DIR"
        else
            print_info "Keeping existing DeepRaccess installation."
            print_info "Verifying installation..."
        fi
    fi

    # Only verify if directory still exists (not removed)
    if [ -d "$DEEPRACCESS_DIR" ]; then
        # Verify key files exist
        if [ -f "$DEEPRACCESS_DIR/mymodel.py" ]; then
            print_info "✓ DeepRaccess model file found"
        else
            print_error "✗ DeepRaccess model file not found!"
            exit 1
        fi

        # Check for pre-trained models
        if [ -d "$DEEPRACCESS_DIR/path" ]; then
            MODEL_COUNT=$(ls -1 "$DEEPRACCESS_DIR/path"/*.pth 2>/dev/null | wc -l)
            if [ $MODEL_COUNT -gt 0 ]; then
                print_info "✓ Found $MODEL_COUNT pre-trained model(s)"
            else
                print_warning "No pre-trained models found in DeepRaccess/path/"
                print_warning "You may need to download models manually"
            fi
        fi

        print_info "DeepRaccess setup verification complete!"
        exit 0
    fi
fi

# Clone DeepRaccess repository
print_info "Cloning DeepRaccess repository..."
if git clone https://github.com/hmdlab/DeepRaccess.git "$DEEPRACCESS_DIR"; then
    print_info "✓ DeepRaccess cloned successfully"
else
    print_error "Failed to clone DeepRaccess repository"
    print_error "Please check your internet connection and try again"
    exit 1
fi

# Verify essential files
print_info "Verifying DeepRaccess installation..."

REQUIRED_FILES=(
    "mymodel.py"
    "__init__.py"
)

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$DEEPRACCESS_DIR/$file" ]; then
        print_info "✓ Found $file"
    else
        print_warning "✗ Missing $file (may not be critical)"
    fi
done

# Check for pre-trained models
print_info "Checking for pre-trained models..."
if [ -d "$DEEPRACCESS_DIR/path" ]; then
    MODEL_COUNT=$(ls -1 "$DEEPRACCESS_DIR/path"/*.pth 2>/dev/null | wc -l)
    if [ $MODEL_COUNT -gt 0 ]; then
        print_info "✓ Found $MODEL_COUNT pre-trained model(s)"
        ls -lh "$DEEPRACCESS_DIR/path"/*.pth 2>/dev/null
    else
        print_warning "No pre-trained models found in DeepRaccess/path/"
        print_warning "The models should be downloaded separately from the DeepRaccess repository"
        print_info ""
        print_info "You can either:"
        print_info "  1. Download models from DeepRaccess repository (recommended)"
        print_info "  2. Use the framework without pre-trained models (random initialization)"
    fi
else
    print_warning "Model directory 'path/' not found in DeepRaccess"
    print_info "Creating directory structure..."
    mkdir -p "$DEEPRACCESS_DIR/path"
fi

# Set environment variable suggestion
print_info ""
print_info "========================================"
print_info "Setup Complete!"
print_info "========================================"
print_info ""
print_info "DeepRaccess has been installed at: $DEEPRACCESS_DIR"
print_info ""
print_info "Optional: Add to your shell profile (~/.bashrc or ~/.zshrc):"
print_info "  export DEEPRACCESS_PATH=\"$DEEPRACCESS_DIR\""
print_info ""
print_info "Note: The ID3 framework will automatically detect DeepRaccess"
print_info "      in the project directory, so setting DEEPRACCESS_PATH is optional."
print_info ""

# Test import
print_info "Testing DeepRaccess import..."
if python3 -c "import sys; sys.path.insert(0, '$DEEPRACCESS_DIR'); from mymodel import FCN; print('✓ Import successful')" 2>/dev/null; then
    print_info "✓ DeepRaccess can be imported successfully"
else
    print_warning "Could not test import (this may be normal if dependencies are missing)"
    print_info "The framework will handle import errors gracefully"
fi

print_info ""
print_info "You can now run: python demo.py"
print_info ""

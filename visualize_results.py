#!/usr/bin/env python3
"""
ID3 Framework - Result Visualization

Visualizes optimization results including convergence curves, accessibility evolution,
and sequence analysis.
"""

import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300


def load_result(result_file):
    """Load optimization result from JSON file"""
    with open(result_file, 'r') as f:
        return json.load(f)


def plot_convergence(result, save_path=None):
    """
    Plot convergence curves for accessibility, CAI, and total loss

    Args:
        result: Result dictionary with trajectory data
        save_path: Path to save figure (optional)
    """
    trajectory = result.get('trajectory', {})

    if not trajectory:
        print("No trajectory data found in result")
        return

    # Extract data
    iterations = trajectory.get('iterations', [])
    accessibility = trajectory.get('accessibility', [])
    cai_values = trajectory.get('discrete_cai_values', trajectory.get('ecai_values', []))
    loss_values = trajectory.get('loss_values', trajectory.get('unified_loss', []))

    if not iterations:
        print("No iteration data found")
        return

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Accessibility convergence
    ax1 = axes[0, 0]
    if accessibility:
        ax1.plot(iterations, accessibility, 'b-', linewidth=2, label='Accessibility')
        ax1.axhline(y=min(accessibility), color='r', linestyle='--', alpha=0.5, label=f'Best: {min(accessibility):.3f}')
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('RNA Accessibility (kcal/mol)')
        ax1.set_title('RNA Accessibility Convergence')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    # Plot 2: CAI evolution
    ax2 = axes[0, 1]
    if cai_values:
        ax2.plot(iterations, cai_values, 'g-', linewidth=2, label='CAI')
        target_cai = result.get('configuration', {}).get('cai_target', 0.8)
        ax2.axhline(y=target_cai, color='r', linestyle='--', alpha=0.5, label=f'Target: {target_cai}')
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('CAI Value')
        ax2.set_title('CAI Evolution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    # Plot 3: Total loss
    ax3 = axes[1, 0]
    if loss_values:
        ax3.plot(iterations, loss_values, 'purple', linewidth=2, label='Total Loss')
        ax3.set_xlabel('Iteration')
        ax3.set_ylabel('Loss')
        ax3.set_title('Total Loss Convergence')
        ax3.set_yscale('log')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

    # Plot 4: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Display summary
    protein = result.get('protein_name', 'Unknown')
    constraint = result.get('constraint_type', 'Unknown')
    variant = result.get('variant', '??')
    seed = result.get('seed', 0)

    final_acc = result.get('final_accessibility', 'N/A')
    best_acc = result.get('best_accessibility', 'N/A')
    improvement = result.get('improvement', 0)

    # Format numbers safely
    final_acc_str = f"{final_acc:.4f}" if isinstance(final_acc, (int, float)) else str(final_acc)
    best_acc_str = f"{best_acc:.4f}" if isinstance(best_acc, (int, float)) else str(best_acc)
    improvement_str = f"{improvement:.4f}" if isinstance(improvement, (int, float)) else str(improvement)

    summary_text = f"""
    Optimization Summary
    ═══════════════════════════

    Protein: {protein}
    Constraint: {constraint.upper()}
    Variant: {variant}
    Seed: {seed}

    Results:
    ────────────────────────────
    Final Accessibility: {final_acc_str}
    Best Accessibility:  {best_acc_str}
    Improvement: {improvement_str}

    Iterations: {len(iterations) if iterations else 0}
    """

    ax4.text(0.1, 0.5, summary_text, fontsize=11, family='monospace',
            verticalalignment='center', transform=ax4.transAxes)

    # Overall title
    fig.suptitle(f'ID3 Optimization Results: {protein} - {constraint.upper()} {variant}',
                fontsize=16, fontweight='bold')

    plt.tight_layout()

    # Save or show
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"✓ Figure saved: {save_path}")
    else:
        plt.show()

    plt.close()


def plot_sequence_analysis(result, save_path=None):
    """
    Plot sequence-level analysis: nucleotide composition, codon usage

    Args:
        result: Result dictionary
        save_path: Path to save figure
    """
    # Get final sequence
    best_seq = result.get('best_seq_design', {}).get('discrete_sequence', '')

    if not best_seq:
        print("No sequence data found")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Nucleotide composition
    ax1 = axes[0]
    nucleotides = ['A', 'U', 'G', 'C']
    counts = [best_seq.count(nt) for nt in nucleotides]
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A']

    ax1.bar(nucleotides, counts, color=colors, alpha=0.7, edgecolor='black')
    ax1.set_xlabel('Nucleotide')
    ax1.set_ylabel('Count')
    ax1.set_title('Nucleotide Composition')
    ax1.grid(True, alpha=0.3, axis='y')

    # Add percentage labels
    total = sum(counts)
    for i, (nt, count) in enumerate(zip(nucleotides, counts)):
        pct = count / total * 100 if total > 0 else 0
        ax1.text(i, count, f'{pct:.1f}%', ha='center', va='bottom')

    # Plot 2: AU vs GC content
    ax2 = axes[1]
    au_count = counts[0] + counts[1]  # A + U
    gc_count = counts[2] + counts[3]  # G + C

    labels = ['AU', 'GC']
    sizes = [au_count, gc_count]
    colors_pie = ['#FF9999', '#66B2FF']

    ax2.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 12})
    ax2.set_title('AU vs GC Content')

    # Add CAI and sequence info
    protein = result.get('protein_name', 'Unknown')
    final_cai = result.get('final_ecai', result.get('final_cai', 'N/A'))
    final_cai_str = f"{final_cai:.4f}" if isinstance(final_cai, (int, float)) else str(final_cai)

    fig.suptitle(f'Sequence Analysis: {protein} | Final CAI: {final_cai_str}',
                fontsize=14, fontweight='bold')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"✓ Figure saved: {save_path}")
    else:
        plt.show()

    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Visualize ID3 optimization results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize result from demo.py output
  python visualize_results.py result.json

  # Visualize result from run_unified_experiment.py
  python visualize_results.py results/20251023_*/O15263_lagrangian_00_seed42.json

  # Save to specific location
  python visualize_results.py result.json --output figures/convergence.png
        """
    )

    parser.add_argument('result_file', type=str, help='Path to result JSON file')
    parser.add_argument('--output', '-o', type=str, help='Output file path for figure')
    parser.add_argument('--type', choices=['convergence', 'sequence', 'both'],
                       default='both', help='Type of plot to generate')

    args = parser.parse_args()

    # Load result
    print(f"Loading result from: {args.result_file}")
    result = load_result(args.result_file)

    protein = result.get('protein_name', 'Unknown')
    constraint = result.get('constraint_type', 'unknown')
    variant = result.get('variant', '??')

    print(f"\nProtein: {protein}")
    print(f"Constraint: {constraint}")
    print(f"Variant: {variant}")
    print(f"Iterations: {len(result.get('trajectory', {}).get('iterations', []))}")

    # Generate plots
    if args.type in ['convergence', 'both']:
        output_conv = args.output if args.output and args.type == 'convergence' else None
        if output_conv is None and args.type == 'both':
            output_conv = f"{protein}_{constraint}_{variant}_convergence.png"
        plot_convergence(result, output_conv)

    if args.type in ['sequence', 'both']:
        output_seq = args.output if args.output and args.type == 'sequence' else None
        if output_seq is None and args.type == 'both':
            output_seq = f"{protein}_{constraint}_{variant}_sequence.png"
        plot_sequence_analysis(result, output_seq)

    print("\n✅ Visualization complete!")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
ID3 Case Study Visualization Script
Enhanced evolution figure with improved visualization matching reference style.

Key features:
- Better color intensity and contrast
- Enhanced probability visualization
- Improved data sampling and processing
- Paper-quality 3-panel layout
"""

import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set matplotlib style for publication-quality figures
plt.style.use('default')
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.grid.axis'] = 'both'
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['figure.dpi'] = 300


def load_trajectory_data(json_file):
    """Load trajectory data from JSON file"""
    print(f"Loading trajectory data from: {json_file}")

    with open(json_file, 'r') as f:
        data = json.load(f)

    trajectory = data['trajectory']

    # Extract core trajectory information
    iterations = trajectory['iterations']
    accessibility = trajectory['accessibility']
    rna_sequences = trajectory.get('rna_sequences', [])
    discrete_sequences = trajectory['discrete_sequences']

    # Calculate basic metrics
    initial_score = accessibility[0]
    best_score = min(accessibility)
    best_step = iterations[accessibility.index(best_score)]
    improvement = (initial_score - best_score) / initial_score * 100

    protein_name = data.get('protein_name', 'Unknown')

    print(f"  Protein: {protein_name}")
    print(f"  Trajectory length: {len(iterations)} steps")
    print(f"  Initial accessibility: {initial_score:.3f}")
    print(f"  Best accessibility: {best_score:.3f}")
    print(f"  Improvement: {improvement:.1f}%")
    print(f"  Convergence step: {best_step}")

    return {
        'iterations': iterations,
        'accessibility': accessibility,
        'rna_sequences': rna_sequences,
        'discrete_sequences': discrete_sequences,
        'protein_name': protein_name,
        'initial_score': initial_score,
        'best_score': best_score,
        'best_step': best_step,
        'improvement': improvement,
        'constraint_type': data.get('constraint_type', 'unknown')
    }


def enhance_probability_contrast(prob_matrix, contrast_factor=2.0, min_intensity=0.3):
    """
    Enhance probability contrast to make visualization more clear.
    Similar to how logits->softmax transformation enhances contrast.

    Args:
        prob_matrix: [4, positions, steps] probability matrix
        contrast_factor: Enhancement strength (higher = more contrast)
        min_intensity: Minimum intensity for visualization

    Returns:
        Enhanced probability matrix with same shape
    """
    enhanced_matrix = np.copy(prob_matrix)

    for pos in range(prob_matrix.shape[1]):
        for step in range(prob_matrix.shape[2]):
            probs = prob_matrix[:, pos, step]

            if np.sum(probs) > 0:
                # Apply contrast enhancement
                enhanced_probs = np.power(probs, 1.0/contrast_factor)
                enhanced_probs = enhanced_probs / np.sum(enhanced_probs)

                # Ensure minimum intensity for visualization
                max_prob = np.max(enhanced_probs)
                if max_prob > 0:
                    enhanced_probs = enhanced_probs * (min_intensity + (1-min_intensity) * max_prob)

                enhanced_matrix[:, pos, step] = enhanced_probs

    return enhanced_matrix


def calculate_nucleotide_composition(discrete_sequences, iterations):
    """Calculate AU content evolution over iterations"""
    composition_data = []

    for i, seq in enumerate(discrete_sequences):
        if i % 25 == 0 or i == len(discrete_sequences) - 1:  # Sample every 25 steps + final
            total_count = len(seq)
            a_count = seq.count('A')
            u_count = seq.count('U')

            au_ratio = (a_count + u_count) / total_count if total_count > 0 else 0.5
            composition_data.append({
                'iteration': iterations[i],
                'au_ratio': au_ratio
            })

    return composition_data


def generate_visualization(json_file, output_file=None):
    """Generate the 3-panel visualization figure"""
    print("\nGenerating 3-panel visualization...")

    # Load data
    data = load_trajectory_data(json_file)

    # Extract components
    iterations = data['iterations']
    accessibility = data['accessibility']
    rna_sequences = data['rna_sequences']
    discrete_sequences = data['discrete_sequences']
    final_sequence = discrete_sequences[-1] if discrete_sequences else ""

    # Calculate AU composition
    composition_data = calculate_nucleotide_composition(discrete_sequences, iterations)

    # Create figure with 3 panels
    fig = plt.figure(figsize=(14, 10))

    # Create GridSpec with height ratios matching reference
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 0.6, 0.5], hspace=0.08)

    # ============= TOP PANEL: Nucleotide Probability Evolution =============
    ax_top = fig.add_subplot(gs[0])

    # Process probability data with better sampling
    num_positions = 45
    num_steps_to_show = 200  # Fixed 200 for optimal visualization

    # Sample steps evenly across trajectory
    total_steps = len(rna_sequences) if rna_sequences else len(discrete_sequences)
    step_indices = np.linspace(0, total_steps-1, num_steps_to_show, dtype=int)

    # Initialize probability matrix (4 nucleotides x positions x steps)
    prob_matrix = np.zeros((4, num_positions, len(step_indices)))

    if rna_sequences:
        # Use provided probability sequences
        for step_idx, step_num in enumerate(step_indices):
            if step_num < len(rna_sequences):
                prob_data = rna_sequences[step_num]
                if prob_data and len(prob_data) >= num_positions:
                    for pos in range(num_positions):
                        if pos < len(prob_data):
                            probs = prob_data[pos]  # [A, C, G, U] probabilities
                            if len(probs) >= 4:
                                prob_matrix[:, pos, step_idx] = probs[:4]
    else:
        # Generate dummy probabilities from discrete sequences
        nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
        for step_idx, step_num in enumerate(step_indices):
            if step_num < len(discrete_sequences):
                seq = discrete_sequences[step_num]
                for pos in range(min(num_positions, len(seq))):
                    nuc = seq[pos]
                    if nuc in nucleotide_map:
                        prob_matrix[nucleotide_map[nuc], pos, step_idx] = 1.0

    # Apply contrast enhancement
    enhanced_prob_matrix = enhance_probability_contrast(prob_matrix, contrast_factor=1.5, min_intensity=0.4)

    # Create nucleotide probability heatmap with RGB color mixing
    nucleotides = ['A', 'C', 'G', 'U']
    colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)]  # Red, Green, Blue, Black

    # Create RGB image
    rgb_image = np.zeros((num_positions, len(step_indices), 3))

    for pos in range(num_positions):
        for step_idx in range(len(step_indices)):
            probs = enhanced_prob_matrix[:, pos, step_idx]

            if np.sum(probs) > 0:
                # Weighted color mixing
                weighted_color = np.zeros(3)
                for nuc_idx, (prob, color) in enumerate(zip(probs, colors)):
                    weighted_color += prob * np.array(color)

                # Boost intensity for better visibility
                intensity = np.max(probs)
                rgb_image[pos, step_idx, :] = weighted_color * (0.5 + 0.5 * intensity)

    # Display RGB image (fixed extent for consistency)
    max_iter_display = max(iterations[-1], 1000) if iterations else 1000
    ax_top.imshow(rgb_image, aspect='auto', origin='lower',
                 extent=[0, max_iter_display, 0, num_positions], interpolation='nearest')

    # Add position lines every 3 bases (codon boundaries)
    for pos in range(3, num_positions, 3):
        ax_top.axhline(y=pos-0.5, color='black', linewidth=0.8, alpha=0.4)

    # Add vertical lines for major step intervals
    for step in range(0, max_iter_display + 1, 200):
        ax_top.axvline(x=step, color='gray', linewidth=0.3, alpha=0.3)

    # Formatting
    ax_top.set_ylabel('Top 45 Positions', fontsize=11, fontweight='bold')
    ax_top.set_xlim(0, max_iter_display)
    ax_top.set_ylim(0, num_positions)
    ax_top.set_title('Nucleotide Evolution', fontsize=12, fontweight='bold')

    # Add final sequence on right y-axis
    ax_top_right = ax_top.twinx()
    ax_top_right.set_ylim(0, num_positions)
    ax_top_right.set_ylabel('Final RNA Sequence', fontsize=11, fontweight='bold',
                           rotation=270, labelpad=20)

    # Add sequence letters
    ax_top_right.set_yticks(range(num_positions))
    ax_top_right.set_yticklabels([final_sequence[i] if i < len(final_sequence) else ''
                                  for i in range(num_positions)],
                                 fontsize=9, fontweight='bold')
    ax_top_right.tick_params(axis='y', colors='black', pad=2)

    # Remove x-axis labels for top panel
    ax_top.set_xticks([])
    ax_top.set_xticklabels([])

    # Add nucleotide legend
    legend_elements = [patches.Patch(facecolor=color, label=nuc, alpha=1.0)
                      for nuc, color in zip(nucleotides, ['red', 'green', 'blue', 'black'])]
    legend = ax_top.legend(handles=legend_elements, loc='upper left', ncol=2, fontsize=11,
                          frameon=True, fancybox=True, shadow=True,
                          bbox_to_anchor=(0.02, 0.98))
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(0.5)

    # ============= MIDDLE PANEL: Convergence Curve =============
    ax_middle = fig.add_subplot(gs[1])

    # Calculate normalized loss (0 to 1) for dual-axis display
    acc_losses = []
    for score in accessibility:
        # Normalize loss: best=0, initial=1
        synthetic_loss = max(0, (score - data['best_score']) / (data['initial_score'] - data['best_score']))
        acc_losses.append(synthetic_loss)

    # Plot accessibility score
    color1 = '#2E86AB'
    ax_middle.plot(iterations, accessibility, color=color1, linewidth=2, alpha=0.8,
                   label='Accessibility Score')
    ax_middle.set_ylabel('Score / Loss', fontsize=11, fontweight='bold')

    # Mark best point
    ax_middle.plot(data['best_step'], data['best_score'], 'r*', markersize=15,
                  label=f'Best Score ({data["best_score"]:.3f})', zorder=5)

    # Add reference lines
    ax_middle.axhline(y=data['initial_score'], color='gray', linestyle='--', alpha=0.5,
                      label=f'Initial ({data["initial_score"]:.3f})')
    ax_middle.axhline(y=data['best_score'], color='red', linestyle='--', alpha=0.5,
                      label=f'Optimal ({data["best_score"]:.3f})')

    # Plot normalized loss data
    color2 = '#F18F01'
    ax_middle.plot(iterations, acc_losses, color=color2, linewidth=2, alpha=0.8,
                  linestyle='--', label='Accessibility Loss')

    # Mark minimum loss point
    if acc_losses:
        min_loss_idx = acc_losses.index(min(acc_losses))
        min_loss_step = iterations[min_loss_idx]
        min_loss_value = acc_losses[min_loss_idx]
        ax_middle.plot(min_loss_step, min_loss_value, 'o', color=color2, markersize=8,
                      label=f'Min Loss ({min_loss_value:.3f})')

    # Legend and formatting
    ax_middle.legend(loc='upper center', fontsize=9)
    ax_middle.grid(True, alpha=0.3)
    ax_middle.set_xlim(0, max_iter_display)
    ax_middle.set_xticks([])
    ax_middle.set_xticklabels([])

    # ============= BOTTOM PANEL: AU Content Evolution =============
    ax_bottom = fig.add_subplot(gs[2])

    # Extract AU content data
    au_steps = [comp['iteration'] for comp in composition_data]
    au_ratios = [comp['au_ratio'] * 100 for comp in composition_data]  # Convert to percentage

    # Plot AU content evolution
    color3 = '#8B4513'  # Brown color for AU content
    ax_bottom.plot(au_steps, au_ratios, color=color3, linewidth=2, alpha=0.8,
                  label='AU Content Ratio')

    # Add reference line at 50% (balanced)
    ax_bottom.axhline(y=50, color='gray', linestyle=':', alpha=0.7, label='Balanced (50%)')

    # Mark initial and final AU content
    if au_ratios:
        initial_au = au_ratios[0]
        final_au = au_ratios[-1]
        ax_bottom.plot(au_steps[0], initial_au, 'o', color=color3, markersize=8,
                      label=f'Initial AU: {initial_au:.1f}%')
        # Force final point to appear at exactly max iteration for visual clarity
        ax_bottom.plot(max_iter_display, final_au, 's', color=color3, markersize=8,
                      label=f'Final AU: {final_au:.1f}%')

    # Formatting
    ax_bottom.set_xlabel('Optimization Step', fontsize=11, fontweight='bold')
    ax_bottom.set_ylabel('AU Content (%)', fontsize=11, fontweight='bold')
    ax_bottom.set_xlim(0, max_iter_display)
    ax_bottom.set_ylim(30, 70)
    ax_bottom.grid(True, alpha=0.3)
    ax_bottom.legend(loc='upper center', fontsize=9)

    # Overall figure title (removed for cleaner appearance)
    # constraint_name = data['constraint_type'].replace('_', ' ').title()
    # fig.suptitle(f'{data["protein_name"]} - {constraint_name} Constraint Optimization Analysis',
    #             fontsize=14, fontweight='bold', y=0.95)

    plt.tight_layout()

    # Save figure
    if output_file is None:
        output_file = Path(json_file).stem + '_visualization.png'

    output_path = Path(output_file)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    # Also save PDF version
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')

    plt.close()

    print(f"\n✅ Visualization complete!")
    print(f"   PNG: {output_path}")
    print(f"   PDF: {pdf_path}")

    return output_path


def main():
    """Main function with command-line interface"""
    parser = argparse.ArgumentParser(description='Generate ID3 evolution figures')
    parser.add_argument('--json', type=str, required=True, help='Path to optimization result JSON file')
    parser.add_argument('--output', type=str, default=None, help='Output PNG file path')

    args = parser.parse_args()

    try:
        generate_visualization(args.json, args.output)
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

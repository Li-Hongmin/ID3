#!/usr/bin/env python3
"""


"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import json
import pandas as pd
import seaborn as sns


plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 14

def load_cai_data():

    json_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables/data_tab2_access_cai_comparison.json"

    with open(json_path, 'r') as f:
        data = json.load(f)


    deterministic_variants = ['ams_00', 'ams_01', 'cpc_00', 'cpc_01']
    variant_labels = {
        'ams_00': 'A0',
        'ams_01': 'A1',
        'cpc_00': 'C0',
        'cpc_01': 'C1'
    }


    cai_data = {v: [] for v in deterministic_variants}


    for protein in data:
        if isinstance(data[protein], dict):
            for variant in deterministic_variants:
                if variant in data[protein]:
                    variant_data = data[protein][variant]
                    if 'discrete_cai' in variant_data:
                        cai_data[variant].append(variant_data['discrete_cai'])
                    elif 'cai' in variant_data:
                        cai_data[variant].append(variant_data['cai'])

    return cai_data, variant_labels, initial_cai

def create_figure():
    """创建Figure 8"""


    cai_data, variant_labels, initial_cai = load_cai_data()


    fig = plt.figure(figsize=(12, 5))
    gs = gridspec.GridSpec(1, 2, wspace=0.4)


    colors = {
        'ams_00': '
        'ams_01': '
        'cpc_00': '
        'cpc_01': '
    }

    # Panel A: Final CAI Distribution
    ax1 = fig.add_subplot(gs[0, 0])


    positions = []
    labels = []
    data_for_plot = []
    box_colors = []

    for i, (variant, label) in enumerate(variant_labels.items()):
        if variant in cai_data and len(cai_data[variant]) > 0:
            positions.append(i + 1)
            labels.append(label)
            data_for_plot.append(cai_data[variant])
            box_colors.append(colors[variant])


    bp = ax1.boxplot(data_for_plot, positions=positions, widths=0.6,
                     patch_artist=True, showmeans=True,
                     meanprops=dict(marker='o', markerfacecolor='white',
                                  markeredgecolor='black', markersize=7),
                     medianprops=dict(color='black', linewidth=1.5),
                     boxprops=dict(linewidth=1.2),
                     whiskerprops=dict(linewidth=1.2),
                     capprops=dict(linewidth=1.2))


    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)


    for i, (pos, data) in enumerate(zip(positions, data_for_plot)):
        mean_val = np.mean(data)
        ax1.text(pos, mean_val + 0.01, f'{mean_val:.3f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')


    ax1.axhline(y=0.8, color='red', linestyle='--', linewidth=2,
               alpha=0.7, label='Target CAI = 0.8')

    ax1.set_xlabel('Variant', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Final CAI', fontsize=11, fontweight='bold')
    ax1.set_title('Panel A: Final CAI Distribution', fontsize=12, fontweight='bold')
    ax1.set_xticks(positions)
    ax1.set_xticklabels(labels)
    ax1.set_ylim(0.75, 1.0)
    ax1.legend(loc='lower right', frameon=True, fancybox=True, shadow=True)
    ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Panel B: CAI Improvement
    ax2 = fig.add_subplot(gs[0, 1])


    improvements = []
    improvement_labels = []
    improvement_colors = []

    for variant, label in variant_labels.items():
        if variant in cai_data and len(cai_data[variant]) > 0:
            final_mean = np.mean(cai_data[variant])
            improvement = ((final_mean - initial_cai) / initial_cai) * 100
            improvements.append(improvement)
            improvement_labels.append(label)
            improvement_colors.append(colors[variant])


    x_pos = np.arange(len(improvements))
    bars = ax2.bar(x_pos, improvements, color=improvement_colors,
                   alpha=0.7, edgecolor='black', linewidth=1.5)


    for i, (bar, improvement) in enumerate(zip(bars, improvements)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'+{improvement:.0f}%',
                ha='center', va='bottom', fontsize=10, fontweight='bold')


    avg_improvement = np.mean(improvements)
    ax2.axhline(y=avg_improvement, color='green', linestyle='--',
               linewidth=2, alpha=0.7,
               label=f'Avg: +{avg_improvement:.0f}%')


    best_idx = np.argmax(improvements)
    best_improvement = improvements[best_idx]
    best_label = improvement_labels[best_idx]


    ax2.annotate(f'Best: {best_label}\n+{best_improvement:.0f}%',
                xy=(best_idx, best_improvement),
                xytext=(best_idx + 0.5, best_improvement + 5),
                arrowprops=dict(arrowstyle='->', color='red', lw=2),
                fontsize=10, color='red', fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))

    ax2.set_xlabel('Variant', fontsize=11, fontweight='bold')
    ax2.set_ylabel('CAI Improvement (%)', fontsize=11, fontweight='bold')
    ax2.set_title('Panel B: CAI Improvement', fontsize=12, fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(improvement_labels)
    ax2.set_ylim(0, max(improvements) + 15)
    ax2.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
    ax2.grid(True, alpha=0.3, axis='y', linestyle='-', linewidth=0.5)


    fig.suptitle('CAI Intelligent Retreat Algorithm Performance Validation',
                fontsize=14, fontweight='bold', y=1.02)


    all_cai = []
    for variant_data in cai_data.values():
        all_cai.extend(variant_data)

    if all_cai:
        success_rate = sum(1 for c in all_cai if c >= 0.8) / len(all_cai) * 100
        stats_text = (f"48 Experiments | Mean CAI: {np.mean(all_cai):.3f} ± {np.std(all_cai):.3f} | "
                     f"Success Rate: {success_rate:.0f}% | Avg Improvement: +{avg_improvement:.0f}%")
        fig.text(0.5, -0.02, stats_text, ha='center', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5))


    plt.tight_layout()


    output_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures"
    fig.savefig(f"{output_path}/figure_8_cai_validation.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_path}/figure_8_cai_validation.png", dpi=300, bbox_inches='tight')

    print(f"Figure 8 saved to {output_path}/figure_8_cai_validation.[pdf/png]")


    print("\n=== CAI Validation Statistics ===")
    for variant, label in variant_labels.items():
        if variant in cai_data and len(cai_data[variant]) > 0:
            mean_cai = np.mean(cai_data[variant])
            std_cai = np.std(cai_data[variant])
            improvement = ((mean_cai - initial_cai) / initial_cai) * 100
            print(f"{label}: {mean_cai:.3f} ± {std_cai:.3f} (Improvement: +{improvement:.0f}%)")

    plt.show()

if __name__ == "__main__":
    create_figure()
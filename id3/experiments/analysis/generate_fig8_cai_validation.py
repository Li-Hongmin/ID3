#!/usr/bin/env python3
"""


Panel A: Final CAI distribution by variant
Panel B: CAI improvement achieved by each variant
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import json
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


plt.style.use('seaborn-v0_8-whitegrid')

def load_cai_data(base_dir):



    json_path = os.path.join(base_dir, "tables", "data_tab2_access_cai_comparison.json")

    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            data = json.load(f)

        return data
    else:

        return None

def extract_cai_values(data):
    """从数据中提取CAI值"""

    cai_values = {}
    initial_cai_values = {}


    if isinstance(data, dict):

        for protein, protein_data in data.items():
            if isinstance(protein_data, dict):

                for variant_key, variant_data in protein_data.items():
                    if isinstance(variant_data, dict):

                        variant_map = {
                            'lagrangian': 'L',
                            'ams': 'A',
                            'cpc': 'C'
                        }


                        parts = variant_key.split('_')
                        if len(parts) == 2:
                            constraint_type = parts[0].lower()
                            strategy = parts[1]

                            if constraint_type in variant_map:
                                variant = variant_map[constraint_type] + strategy

                                if variant not in cai_values:
                                    cai_values[variant] = []
                                    initial_cai_values[variant] = []


                                if 'discrete_cai' in variant_data:
                                    cai_values[variant].append(variant_data['discrete_cai'])
                                elif 'cai' in variant_data:
                                    cai_values[variant].append(variant_data['cai'])


                                if 'initial_cai' in variant_data:
                                    initial_cai_values[variant].append(variant_data['initial_cai'])

    return cai_values, initial_cai_values

def create_cai_validation_figure(base_dir, output_path):



    data = load_cai_data(base_dir)

    if data is None:

        return


    cai_values, initial_cai_values = extract_cai_values(data)


    if not cai_values:

        csv_path = os.path.join(base_dir, "tables", "table2_cai_comparison.csv")
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)


            cai_columns = [col for col in df.columns if 'With_Penalty_CAI' in col]

            for col in cai_columns:
                variant = col.split(',')[-1] if ',' in col else col.split('_')[-1]

                    if variant not in cai_values:
                        cai_values[variant] = []


                    values = df[col].dropna()
                    for val in values:
                        try:

                            if isinstance(val, str):
                                val = val.replace('w', '').strip()
                            cai_values[variant].append(float(val))
                        except:
                            continue

    if not cai_values:

        return


    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(1, 2, wspace=0.3)


    constraint_colors = {



    }


    variant_order = ['L00', 'L01', 'L10', 'L11',
                    'A00', 'A01', 'A10', 'A11',
                    'C00', 'C01', 'C10', 'C11']


    available_variants = [v for v in variant_order if v in cai_values and len(cai_values[v]) > 0]

    # Panel A: Final CAI Distribution
    ax1 = fig.add_subplot(gs[0, 0])


    plot_data = []
    plot_labels = []
    plot_colors = []

    for variant in available_variants:
        if variant in cai_values and len(cai_values[variant]) > 0:
            plot_data.append(cai_values[variant])
            plot_labels.append(variant)
            plot_colors.append(constraint_colors.get(variant[0], 'gray'))

    if plot_data:

        bp = ax1.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                         showmeans=True, meanline=False,
                         medianprops=dict(color='black', linewidth=2),
                         meanprops=dict(marker='o', markerfacecolor='white',
                                      markeredgecolor='black', markersize=8))


        for patch, color in zip(bp['boxes'], plot_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)


        for i, (variant, data) in enumerate(zip(plot_labels, plot_data), 1):
            mean_val = np.mean(data)
            ax1.text(i, mean_val + 0.01, f'{mean_val:.3f}',
                    ha='center', va='bottom', fontsize=9, fontweight='bold')


        ax1.axhline(y=0.8, color='red', linestyle='--', linewidth=2,
                   label='Target CAI = 0.8', alpha=0.7)

        ax1.set_xlabel('ID3 Variant', fontsize=12)
        ax1.set_ylabel('Final CAI Value', fontsize=12)
        ax1.set_title('Panel A: Final CAI Distribution by Variant', fontsize=14, fontweight='bold')
        ax1.set_ylim(0.75, 1.0)
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='lower right')


        ax1.set_xticklabels(plot_labels, rotation=45, ha='right')

    # Panel B: CAI Improvement
    ax2 = fig.add_subplot(gs[0, 1])


    improvements = {}
    for variant in available_variants:
        if variant in cai_values and len(cai_values[variant]) > 0:
            final_mean = np.mean(cai_values[variant])


            initial_cai = 0.563
            if variant in initial_cai_values and len(initial_cai_values[variant]) > 0:
                initial_cai = np.mean(initial_cai_values[variant])

            improvement = ((final_mean - initial_cai) / initial_cai) * 100
            improvements[variant] = improvement

    if improvements:

        variants = list(improvements.keys())
        improvement_values = list(improvements.values())
        colors = [constraint_colors.get(v[0], 'gray') for v in variants]

        bars = ax2.bar(range(len(variants)), improvement_values, color=colors, alpha=0.7,
                      edgecolor='black', linewidth=1.5)


        for i, (variant, improvement) in enumerate(improvements.items()):
            ax2.text(i, improvement + 1, f'+{improvement:.1f}%',
                    ha='center', va='bottom', fontsize=9, fontweight='bold')


        avg_improvement = np.mean(improvement_values)
        ax2.axhline(y=avg_improvement, color='green', linestyle='--', linewidth=2,
                   label=f'Average Improvement: +{avg_improvement:.1f}%', alpha=0.7)

        ax2.set_xticks(range(len(variants)))
        ax2.set_xticklabels(variants, rotation=45, ha='right')
        ax2.set_xlabel('ID3 Variant', fontsize=12)
        ax2.set_ylabel('CAI Improvement (%)', fontsize=12)
        ax2.set_title('Panel B: CAI Improvement from Initial Sequences', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.legend(loc='upper right')


        best_variant = max(improvements, key=improvements.get)
        best_improvement = improvements[best_variant]
        best_idx = variants.index(best_variant)


        ax2.annotate(f'Best: {best_variant}\n+{best_improvement:.1f}%',
                    xy=(best_idx, best_improvement),
                    xytext=(best_idx + 1, best_improvement + 10),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2),
                    fontsize=10, color='red', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))


    fig.suptitle('CAI Intelligent Retreat Algorithm Performance Validation',
                fontsize=16, fontweight='bold')


    if cai_values:
        all_cai = []
        for variant_data in cai_values.values():
            all_cai.extend(variant_data)

        if all_cai:
            stats_text = (f"Total Experiments: {len(all_cai)} | "
                         f"Mean CAI: {np.mean(all_cai):.3f} | "
                         f"Success Rate (CAI ≥ 0.8): {sum(1 for c in all_cai if c >= 0.8)/len(all_cai)*100:.1f}%")
            fig.text(0.5, 0.02, stats_text, ha='center', fontsize=11,
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5))


    plt.tight_layout(rect=[0, 0.05, 1, 0.96])


    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    base_name = os.path.splitext(output_path)[0]
    plt.savefig(f"{base_name}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{base_name}.png", dpi=300, bbox_inches='tight')


    plt.show()



    print("-" * 50)
    for variant in sorted(cai_values.keys()):
        if len(cai_values[variant]) > 0:
            mean_cai = np.mean(cai_values[variant])
            std_cai = np.std(cai_values[variant])
            print(f"{variant}: {mean_cai:.3f} ± {std_cai:.3f} (n={len(cai_values[variant])})")

def main():
    """主函数"""
    base_dir = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results"
    output_path = os.path.join(base_dir, "figures", "fig8_cai_validation")

    create_cai_validation_figure(base_dir, output_path)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""


"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['figure.titlesize'] = 12

def load_trajectory_data():

    csv_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures/trajectory_access_only_1000steps.csv"
    df = pd.read_csv(csv_path)


    numeric_cols = ['step', 'true_accessibility']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df

def create_convergence_figure():
    """创建Figure 9"""


    df = load_trajectory_data()


    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3,
                          height_ratios=[1.2, 0.8])


    colors = {
        'A': '
        'C': '
        'L': '
    }

    variants = sorted(df['variant'].unique())


    gs_a = gridspec.GridSpecFromSubplotSpec(3, 4, subplot_spec=gs[0, 0],
                                           hspace=0.5, wspace=0.4)

    for idx, variant in enumerate(variants):
        ax = fig.add_subplot(gs_a[idx // 4, idx % 4])


        variant_data = df[df['variant'] == variant]


        grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean', 'std', 'count'])
        grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])


        constraint_type = variant.split('-')[1][0]
        color = colors.get(constraint_type, 'gray')


        ax.plot(grouped.index, grouped['mean'], color=color, linewidth=1.5, alpha=0.9)


        ax.fill_between(grouped.index,
                        grouped['mean'] - grouped['sem'],
                        grouped['mean'] + grouped['sem'],
                        color=color, alpha=0.2)


        ax.set_title(variant.replace('ID3-', ''), fontsize=8, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3, linewidth=0.3)


        if idx >= 8:
            ax.set_xlabel('Steps', fontsize=7)
        if idx % 4 == 0:
            ax.set_ylabel('Accessibility', fontsize=7)


    gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 1],
                                           hspace=0.4, wspace=0.4)

    strategies = ['00', '01', '10', '11']
    strategy_names = {
        '00': 'Det-Soft',
        '01': 'Det-Hard',
        '10': 'Stoch-Soft',
        '11': 'Stoch-Hard'
    }

    for idx, strategy in enumerate(strategies):
        ax = fig.add_subplot(gs_b[idx // 2, idx % 2])


        strategy_variants = [v for v in variants if v.endswith(strategy)]

        for variant in strategy_variants:
            variant_data = df[df['variant'] == variant]
            grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean'])

            constraint_type = variant.split('-')[1][0]
            color = colors.get(constraint_type, 'gray')
            label = {'L': 'Lagrangian', 'A': 'AMS', 'C': 'CPC'}[constraint_type]

            ax.plot(grouped.index, grouped['mean'], color=color,
                   linewidth=1.5, label=label, alpha=0.8)

        ax.set_title(strategy_names[strategy], fontsize=9, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)
        ax.set_xlabel('Steps', fontsize=8)
        ax.set_ylabel('Accessibility', fontsize=8)
        ax.legend(loc='upper right', fontsize=6, frameon=True)
        ax.grid(True, alpha=0.3, linewidth=0.3)
        ax.tick_params(labelsize=7)


    gs_c = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 0],
                                           hspace=0.4)

    constraint_names = {
        'L': 'Lagrangian Multiplier',
        'A': 'Amino Matching Softmax',
        'C': 'Codon Profile Constraint'
    }

    for idx, (prefix, name) in enumerate(constraint_names.items()):
        ax = fig.add_subplot(gs_c[idx, 0])


        constraint_variants = [v for v in variants if v.split('-')[1].startswith(prefix)]

        for variant in constraint_variants:
            variant_data = df[df['variant'] == variant]
            grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean'])

            strategy = variant[-2:]
            linestyle = {
                '00': '-',
                '01': '--',
                '10': '-.',
                '11': ':'
            }[strategy]

            ax.plot(grouped.index, grouped['mean'],
                   color=colors[prefix], linestyle=linestyle,
                   linewidth=1.5, label=strategy_names[strategy], alpha=0.8)

        ax.set_title(name, fontsize=9, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)

        if idx == 2:
            ax.set_xlabel('Steps', fontsize=8)
        ax.set_ylabel('Accessibility', fontsize=8)
        ax.legend(loc='upper right', fontsize=6, ncol=2, frameon=True)
        ax.grid(True, alpha=0.3, linewidth=0.3)
        ax.tick_params(labelsize=7)


    ax_d = fig.add_subplot(gs[1, 1])


    final_performance = []
    for variant in variants:
        variant_data = df[(df['variant'] == variant) & (df['step'] >= 900)]
        mean_val = variant_data['true_accessibility'].mean()
        std_val = variant_data['true_accessibility'].std()
        constraint_type = variant.split('-')[1][0]
        final_performance.append({
            'variant': variant,
            'mean': mean_val,
            'std': std_val,
            'constraint': constraint_type
        })

    final_df = pd.DataFrame(final_performance)
    final_df = final_df.sort_values('mean')


    x_pos = np.arange(len(final_df))
    bar_colors = [colors[row['constraint']] for _, row in final_df.iterrows()]

    bars = ax_d.bar(x_pos, final_df['mean'], color=bar_colors,
                   alpha=0.7, edgecolor='black', linewidth=1)


    ax_d.errorbar(x_pos, final_df['mean'], yerr=final_df['std'],
                 fmt='none', color='black', capsize=3, capthick=1)


    ax_d.set_xticks(x_pos)
    ax_d.set_xticklabels([v.replace('ID3-', '') for v in final_df['variant']],
                         rotation=45, ha='right', fontsize=7)

    ax_d.set_xlabel('Variant', fontsize=9, fontweight='bold')
    ax_d.set_ylabel('Final Accessibility (kcal/mol)', fontsize=9, fontweight='bold')
    ax_d.set_title('Panel D: Final Performance (Steps 900-1000)',
                  fontsize=10, fontweight='bold')


    ax_d.yaxis.grid(True, alpha=0.3, linewidth=0.3)
    ax_d.set_axisbelow(True)


    best_idx = 0
    best_variant = final_df.iloc[best_idx]['variant']
    best_mean = final_df.iloc[best_idx]['mean']
    best_std = final_df.iloc[best_idx]['std']

    ax_d.annotate(f'Best: {best_variant}\n{best_mean:.3f} ± {best_std:.3f}',
                 xy=(best_idx, best_mean),
                 xytext=(best_idx + 2, best_mean + 0.3),
                 arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                 fontsize=8, color='red', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))


    constraint_stats = []
    for prefix, name in [('L', 'Lagrangian'), ('A', 'AMS'), ('C', 'CPC')]:
        constraint_data = final_df[final_df['constraint'] == prefix]
        mean_perf = constraint_data['mean'].mean()
        constraint_stats.append(f"{name}: {mean_perf:.3f}")

    stats_text = ' | '.join(constraint_stats)
    ax_d.text(0.5, 0.02, stats_text, transform=ax_d.transAxes,
             fontsize=8, ha='center',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=0.5))


    fig.suptitle('ID3 Framework: Comprehensive Convergence Analysis',
                fontsize=14, fontweight='bold', y=0.98)


    fig.text(0.02, 0.95, 'Panel A: Individual Variants', fontsize=10, fontweight='bold')
    fig.text(0.52, 0.95, 'Panel B: By Strategy', fontsize=10, fontweight='bold')
    fig.text(0.02, 0.45, 'Panel C: By Constraint', fontsize=10, fontweight='bold')


    plt.tight_layout(rect=[0, 0.01, 1, 0.96])


    output_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures"
    fig.savefig(f"{output_path}/figure_9_convergence_analysis.pdf", dpi=300, bbox_inches='tight')
    fig.savefig(f"{output_path}/figure_9_convergence_analysis.png", dpi=300, bbox_inches='tight')

    print(f"Figure 9 saved to {output_path}/figure_9_convergence_analysis.[pdf/png]")


    print("\n=== Convergence Analysis Statistics ===")
    print(f"Best variant: {best_variant} ({best_mean:.3f} ± {best_std:.3f})")

    for prefix, name in [('L', 'Lagrangian'), ('A', 'AMS'), ('C', 'CPC')]:
        constraint_data = final_df[final_df['constraint'] == prefix]
        mean_perf = constraint_data['mean'].mean()
        std_perf = constraint_data['mean'].std()
        print(f"{name} average: {mean_perf:.3f} ± {std_perf:.3f}")

    plt.show()

if __name__ == "__main__":
    create_convergence_figure()
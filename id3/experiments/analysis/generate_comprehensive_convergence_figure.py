#!/usr/bin/env python3
"""






"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')


plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("tab20")

def load_trajectory_data(access_only_path, cai_path=None):




    df_access = pd.read_csv(access_only_path)
    df_access['experiment_type'] = 'Access-only'


    if cai_path and os.path.exists(cai_path):
        df_cai = pd.read_csv(cai_path)
        df_cai['experiment_type'] = 'Access+CAI'
        df = pd.concat([df_access, df_cai], ignore_index=True)
    else:
        df = df_access


    numeric_cols = ['step', 'accessibility_loss', 'true_accessibility',
                   'constraint_loss', 'temperature']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')





    return df

def create_comprehensive_figure(df, output_path='fig_convergence_comprehensive.pdf'):
    """创建综合收敛分析图"""


    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.25)


    variants = sorted(df['variant'].unique())


    ams_color = '
    cpc_color = '
    lag_color = '


    variant_colors = {
        'ID3-A00': ams_color, 'ID3-A01': ams_color,
        'ID3-A10': ams_color, 'ID3-A11': ams_color,
        'ID3-C00': cpc_color, 'ID3-C01': cpc_color,
        'ID3-C10': cpc_color, 'ID3-C11': cpc_color,
        'ID3-L00': lag_color, 'ID3-L01': lag_color,
        'ID3-L10': lag_color, 'ID3-L11': lag_color
    }


    variant_labels = {
        'ID3-L00': 'Lagrangian Deterministic Soft',
        'ID3-L01': 'Lagrangian Deterministic Hard',
        'ID3-L10': 'Lagrangian Stochastic Soft',
        'ID3-L11': 'Lagrangian Stochastic Hard',
        'ID3-A00': 'Amino Matching Deterministic Soft',
        'ID3-A01': 'Amino Matching Deterministic Hard',
        'ID3-A10': 'Amino Matching Stochastic Soft',
        'ID3-A11': 'Amino Matching Stochastic Hard',
        'ID3-C00': 'Codon Profile Deterministic Soft',
        'ID3-C01': 'Codon Profile Deterministic Hard',
        'ID3-C10': 'Codon Profile Stochastic Soft',
        'ID3-C11': 'Codon Profile Stochastic Hard'
    }


    gs_a = gridspec.GridSpecFromSubplotSpec(3, 4, subplot_spec=gs[0, 0],
                                           hspace=0.4, wspace=0.3)

    for idx, variant in enumerate(variants):
        ax = fig.add_subplot(gs_a[idx // 4, idx % 4])


        variant_data = df[df['variant'] == variant]


        grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean', 'std', 'count'])
        grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])


        ax.plot(grouped.index, grouped['mean'],
               color=variant_colors.get(variant, 'gray'),
               linewidth=2, label=variant)

        ax.fill_between(grouped.index,
                        grouped['mean'] - grouped['sem'],
                        grouped['mean'] + grouped['sem'],
                        color=variant_colors.get(variant, 'gray'),
                        alpha=0.2)


        ax.set_title(variant, fontsize=10, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)


        if idx >= 8:
            ax.set_xlabel('Optimization Steps', fontsize=8)
        if idx % 4 == 0:
            ax.set_ylabel('Accessibility (kcal/mol)', fontsize=8)

        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)


    fig.text(0.25, 0.95, 'A. Individual Convergence Curves for All 12 ID3 Variants',
            fontsize=14, fontweight='bold', ha='center')


    gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 1],
                                           hspace=0.3, wspace=0.3)

    strategies = ['00', '01', '10', '11']
    strategy_names = {
        '00': 'Deterministic Soft',
        '01': 'Deterministic Hard',
        '10': 'Stochastic Soft',
        '11': 'Stochastic Hard'
    }

    for idx, strategy in enumerate(strategies):
        ax = fig.add_subplot(gs_b[idx // 2, idx % 2])


        strategy_variants = [v for v in variants if v.endswith(strategy)]

        for variant in strategy_variants:
            variant_data = df[df['variant'] == variant]
            grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean', 'std', 'count'])
            grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])

            constraint_type = variant.split('-')[1][0]
            label = {
                'L': 'Lagrangian',
                'A': 'Amino Matching',
                'C': 'Codon Profile'
            }[constraint_type]

            ax.plot(grouped.index, grouped['mean'],
                   color=variant_colors.get(variant, 'gray'),
                   linewidth=2, label=label, alpha=0.8)

        ax.set_title(strategy_names[strategy], fontsize=11, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)
        ax.set_xlabel('Optimization Steps', fontsize=9)
        ax.set_ylabel('Accessibility (kcal/mol)', fontsize=9)
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)


    fig.text(0.75, 0.95, 'B. Performance by Base Strategy',
            fontsize=14, fontweight='bold', ha='center')


    gs_c = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 0],
                                           hspace=0.3)

    constraint_types = [('L', 'Lagrangian Multiplier', lag_color),
                       ('A', 'Amino Matching Softmax', ams_color),
                       ('C', 'Codon Profile Constraint', cpc_color)]

    for idx, (prefix, name, color) in enumerate(constraint_types):
        ax = fig.add_subplot(gs_c[idx, 0])


        constraint_variants = [v for v in variants if v.split('-')[1].startswith(prefix)]

        for variant in constraint_variants:
            variant_data = df[df['variant'] == variant]
            grouped = variant_data.groupby('step')['true_accessibility'].agg(['mean', 'std', 'count'])
            grouped['sem'] = grouped['std'] / np.sqrt(grouped['count'])


            strategy = variant[-2:]
            linestyle = {
                '00': '-',
                '01': '--',
                '10': '-.',
                '11': ':'
            }[strategy]

            ax.plot(grouped.index, grouped['mean'],
                   color=color, linestyle=linestyle,
                   linewidth=2, label=strategy_names[strategy], alpha=0.8)

        ax.set_title(name, fontsize=11, fontweight='bold')
        ax.set_xlim(0, 1000)
        ax.set_ylim(0.5, 2.5)

        if idx == 2:
            ax.set_xlabel('Optimization Steps', fontsize=10)
        ax.set_ylabel('Accessibility (kcal/mol)', fontsize=9)
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.3)


    fig.text(0.25, 0.48, 'C. Performance by Constraint Mechanism',
            fontsize=14, fontweight='bold', ha='center')


    ax_d = fig.add_subplot(gs[1, 1])


    final_performance = []
    for variant in variants:
        variant_data = df[(df['variant'] == variant) & (df['step'] >= 900)]
        mean_val = variant_data['true_accessibility'].mean()
        std_val = variant_data['true_accessibility'].std()
        final_performance.append({
            'variant': variant,
            'mean': mean_val,
            'std': std_val
        })

    final_df = pd.DataFrame(final_performance)
    final_df = final_df.sort_values('mean')


    bars = ax_d.bar(range(len(final_df)), final_df['mean'],
                    color=[variant_colors.get(v, 'gray') for v in final_df['variant']],
                    edgecolor='black', linewidth=1)


    ax_d.errorbar(range(len(final_df)), final_df['mean'], yerr=final_df['std'],
                 fmt='none', color='black', capsize=3, capthick=1)


    ax_d.set_xticks(range(len(final_df)))
    ax_d.set_xticklabels([v.replace('ID3-', '') for v in final_df['variant']],
                         rotation=45, ha='right')

    ax_d.set_xlabel('ID3 Variant', fontsize=11)
    ax_d.set_ylabel('Final Accessibility (kcal/mol)', fontsize=11)
    ax_d.set_title('D. Final Performance Comparison (Steps 900-1000)',
                  fontsize=14, fontweight='bold')


    ax_d.yaxis.grid(True, alpha=0.3)
    ax_d.set_axisbelow(True)


    best_variant = final_df.iloc[0]['variant']
    best_mean = final_df.iloc[0]['mean']
    best_std = final_df.iloc[0]['std']

    ax_d.annotate(f'Best: {best_variant}\n{best_mean:.3f} ± {best_std:.3f}',
                 xy=(0, best_mean), xytext=(2, best_mean + 0.3),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2),
                 fontsize=10, color='red', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.3))


    constraint_stats = []
    for prefix, name, color in constraint_types:
        constraint_variants = [v for v in variants if v.split('-')[1].startswith(prefix)]
        constraint_data = final_df[final_df['variant'].isin(constraint_variants)]
        mean_perf = constraint_data['mean'].mean()
        std_perf = constraint_data['mean'].std()
        constraint_stats.append(f"{name}: {mean_perf:.3f} ± {std_perf:.3f}")


    stats_text = ' | '.join(constraint_stats)
    fig.text(0.5, 0.02, stats_text, fontsize=11, ha='center',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5))


    fig.suptitle('ID3 Framework: Comprehensive Convergence Analysis',
                fontsize=16, fontweight='bold', y=0.98)


    plt.tight_layout(rect=[0, 0.03, 1, 0.96])


    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    base_name = os.path.splitext(output_path)[0]
    plt.savefig(f"{base_name}.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(f"{base_name}.png", dpi=300, bbox_inches='tight')
    print(f"图形已保存: {base_name}.pdf 和 {base_name}.png")

    plt.show()

    return final_df

def main():


    base_dir = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results"


    access_only_path = os.path.join(base_dir, "figures", "trajectory_access_only_1000steps.csv")
    cai_path = os.path.join(base_dir, "figures", "trajectory_cai_1000steps.csv")


    output_path = os.path.join(base_dir, "figures", "fig_convergence_comprehensive")


    df = load_trajectory_data(access_only_path, cai_path)


    df_access = df[df['experiment_type'] == 'Access-only'] if 'experiment_type' in df.columns else df

    final_performance = create_comprehensive_figure(
        df_access,
        output_path=os.path.join(base_dir, "figures", "fig9_convergence_comprehensive_access")
    )



    print(final_performance.to_string())


    if 'experiment_type' in df.columns and 'Access+CAI' in df['experiment_type'].unique():
        df_cai = df[df['experiment_type'] == 'Access+CAI']

        final_performance_cai = create_comprehensive_figure(
            df_cai,
            output_path=os.path.join(base_dir, "figures", "fig9_convergence_comprehensive_cai")
        )


        print(final_performance_cai.to_string())

if __name__ == "__main__":
    main()
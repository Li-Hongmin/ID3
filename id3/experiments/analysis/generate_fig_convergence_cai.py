#!/usr/bin/env python3
"""





"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("tab20")

print("🔄 加载Access+CAI轨迹数据...")
# Load the trajectory data from our parallel extraction
df = pd.read_csv('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures/trajectory_cai_1000steps.csv')

# Remove header row if present
df = df[df['protein'] != 'protein'].copy()

# Convert numeric columns
numeric_cols = ['step', 'true_accessibility', 'accessibility_loss']
for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

print(f"📊 数据加载完成: {len(df)} 记录")
print(f"🧬 变体: {sorted(df['variant'].unique())}")


final_step_data = df[df['step'] == 999]
if len(final_step_data) > 0:
    print(f"🎯 最终步骤性能范围: {final_step_data['true_accessibility'].min():.3f} - {final_step_data['true_accessibility'].max():.3f}")
    print(f"🎯 最终步骤平均性能: {final_step_data['true_accessibility'].mean():.3f}")

# Create the figure - 4 main panels
fig = plt.figure(figsize=(20, 16))
gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.25)

# Define color scheme
variants = sorted(df['variant'].unique())

# Color scheme for constraint types
ams_color = '#2E86AB'  # Blue for Amino Matching Softmax
cpc_color = '#A23B72'  # Purple for Codon Profile Constraint
lag_color = '#F18F01'  # Orange for Lagrangian Multiplier

# Variant colors by constraint type
variant_colors = {}
for variant in variants:
    if 'A' in variant:
        variant_colors[variant] = ams_color
    elif 'C' in variant:
        variant_colors[variant] = cpc_color
    elif 'L' in variant:
        variant_colors[variant] = lag_color

# ============= PANEL A: 12 Individual Variant Subplots =============
print("🎨 生成Panel A: 12个独立变体子图...")
gs_a = gridspec.GridSpecFromSubplotSpec(3, 4, subplot_spec=gs[0, 0], 
                                        hspace=0.4, wspace=0.3)

# Determine common y-axis range for consistency
all_means = []
for variant in variants:
    variant_data = df[df['variant'] == variant].copy()
    if len(variant_data) > 0:
        step_means = variant_data.groupby('step')['true_accessibility'].mean()
        all_means.extend(step_means.values)

y_min = min(all_means) * 0.95 if all_means else 0
y_max = max(all_means) * 1.05 if all_means else 3

print(f"🎯 Y轴范围: {y_min:.3f} - {y_max:.3f}")

# Create subplot for each variant
for idx, variant in enumerate(variants):
    row = idx // 4
    col = idx % 4
    ax = fig.add_subplot(gs_a[row, col])
    
    variant_data = df[df['variant'] == variant].copy()
    
    if len(variant_data) > 0:
        # Group by step and calculate statistics
        step_stats = variant_data.groupby('step')['true_accessibility'].agg(['mean', 'std', 'count']).reset_index()
        
        # Plot mean line
        ax.plot(step_stats['step'], step_stats['mean'], 
                color=variant_colors.get(variant, '#333333'), linewidth=1.5, alpha=0.9)
        
        # Add shaded area for std
        if 'std' in step_stats.columns and len(step_stats) > 0:
            ax.fill_between(step_stats['step'], 
                            step_stats['mean'] - step_stats['std'], 
                            step_stats['mean'] + step_stats['std'], 
                            color=variant_colors.get(variant, '#333333'), alpha=0.2)
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, 1000)
    ax.set_title(variant, fontsize=9, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Only add axis labels for edge subplots
    if row == 2:
        ax.set_xlabel('Step', fontsize=8)
    if col == 0:
        ax.set_ylabel('Accessibility (kcal/mol)', fontsize=8)
    
    ax.tick_params(labelsize=7)

# Add panel title
fig.text(0.25, 0.92, 'A) All 12 ID3 Variants - Access+CAI Optimization', fontsize=12, fontweight='bold', ha='center')

# ============= PANEL B: 4 Base Strategy Subplots =============
print("🎨 生成Panel B: 4个基础策略子图...")
gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0, 1], 
                                        hspace=0.3, wspace=0.3)

base_strategy_groups = {
    'ID3-00\\n(Det. Soft)': ['ID3-A00', 'ID3-C00', 'ID3-L00'],
    'ID3-01\\n(Stoch. Soft)': ['ID3-A01', 'ID3-C01', 'ID3-L01'],
    'ID3-10\\n(Det. Hard)': ['ID3-A10', 'ID3-C10', 'ID3-L10'],
    'ID3-11\\n(Stoch. Hard)': ['ID3-A11', 'ID3-C11', 'ID3-L11']
}

# Colors for each variant within base strategies
base_variant_colors = {
    'Amino Matching Softmax': '#2E86AB',
    'Codon Profile Constraint': '#A23B72', 
    'Lagrangian Multiplier': '#F18F01'
}

strategy_names = list(base_strategy_groups.keys())
for idx, (base_name, variant_list) in enumerate(base_strategy_groups.items()):
    row = idx // 2
    col = idx % 2
    ax = fig.add_subplot(gs_b[row, col])
    
    for variant in variant_list:
        if variant in df['variant'].values:
            variant_data = df[df['variant'] == variant].copy()
            if len(variant_data) > 0:
                step_means = variant_data.groupby('step')['true_accessibility'].mean()
                
                # Determine constraint type for color
                if 'A' in variant:
                    color = base_variant_colors['Amino Matching Softmax']
                    label = f"{variant} (AMS)"
                elif 'C' in variant:
                    color = base_variant_colors['Codon Profile Constraint']
                    label = f"{variant} (CPC)"
                else:
                    color = base_variant_colors['Lagrangian Multiplier']
                    label = f"{variant} (Lagrangian)"
                
                ax.plot(step_means.index, step_means.values, 
                       color=color, linewidth=1.5, label=label, alpha=0.8)
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, 1000)
    ax.set_title(base_name, fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
    
    # Only add axis labels for edge subplots
    if row == 1:
        ax.set_xlabel('Optimization Step', fontsize=9)
    if col == 0:
        ax.set_ylabel('True Accessibility (kcal/mol)', fontsize=9)
    
    ax.tick_params(labelsize=8)

# Add panel title
fig.text(0.75, 0.92, 'B) Convergence by Base Strategy', fontsize=12, fontweight='bold', ha='center')

# ============= PANEL C: 3 Constraint Mechanism Subplots =============
print("🎨 生成Panel C: 3个约束机制子图...")
gs_c = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1, 0], 
                                        hspace=0.3, wspace=0.3)

constraint_groups = {
    'Amino Matching\\nSoftmax (AMS)': ['ID3-A00', 'ID3-A01', 'ID3-A10', 'ID3-A11'],
    'Codon Profile\\nConstraint (CPC)': ['ID3-C00', 'ID3-C01', 'ID3-C10', 'ID3-C11'],
    'Lagrangian\\nMultiplier': ['ID3-L00', 'ID3-L01', 'ID3-L10', 'ID3-L11']
}

# Colors for each variant within constraints
strategy_colors = {
    '00': '#e74c3c',  # Red for Det. Soft
    '01': '#3498db',  # Blue for Stoch. Soft
    '10': '#2ecc71',  # Green for Det. Hard
    '11': '#f39c12'   # Orange for Stoch. Hard
}

for idx, (mechanism, variant_list) in enumerate(constraint_groups.items()):
    ax = fig.add_subplot(gs_c[idx, 0])
    
    for variant in variant_list:
        if variant in df['variant'].values:
            variant_data = df[df['variant'] == variant].copy()
            if len(variant_data) > 0:
                step_means = variant_data.groupby('step')['true_accessibility'].mean()
                
                # Determine strategy type for color and label
                suffix = variant[-2:]
                color = strategy_colors.get(suffix, '#333333')
                
                strategy_labels = {
                    '00': 'Det. Soft',
                    '01': 'Stoch. Soft',
                    '10': 'Det. Hard',
                    '11': 'Stoch. Hard'
                }
                
                ax.plot(step_means.index, step_means.values, 
                       color=color, linewidth=1.5, 
                       label=f"{variant} ({strategy_labels.get(suffix, '')})", alpha=0.8)
    
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(0, 1000)
    ax.set_title(mechanism, fontsize=10, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
    
    # Only add x-axis label for bottom subplot
    if idx == 2:
        ax.set_xlabel('Optimization Step', fontsize=9)
    
    ax.set_ylabel('True Accessibility (kcal/mol)', fontsize=9)
    
    ax.tick_params(labelsize=8)

# Add panel title
fig.text(0.25, 0.485, 'C) Convergence by Constraint Mechanism', fontsize=12, fontweight='bold', ha='center')

# ============= PANEL D: Final Performance Comparison =============
print("🎨 生成Panel D: 最终性能对比...")
ax4 = fig.add_subplot(gs[1, 1])

# Get final step data for each variant
final_performance = []
for variant in variants:
    variant_data = df[(df['variant'] == variant) & (df['step'] == 999)]
    if len(variant_data) > 0:
        final_values = variant_data['true_accessibility'].values
        final_performance.extend([(variant, val) for val in final_values])

if final_performance:
    final_df = pd.DataFrame(final_performance, columns=['Variant', 'Final_Accessibility'])
    
    # Create grouped bar plot
    variant_means = []
    variant_stds = []
    for v in variants:
        v_data = final_df[final_df['Variant'] == v]['Final_Accessibility']
        if len(v_data) > 0:
            variant_means.append(v_data.mean())
            variant_stds.append(v_data.std() if len(v_data) > 1 else 0)
        else:
            variant_means.append(0)
            variant_stds.append(0)
    
    x_pos = np.arange(len(variants))
    bars = ax4.bar(x_pos, variant_means, yerr=variant_stds, 
                   capsize=5, alpha=0.8, edgecolor='black', linewidth=1)
    
    # Color bars by variant
    for bar, variant in zip(bars, variants):
        bar.set_facecolor(variant_colors.get(variant, '#333333'))
    
    ax4.set_xlabel('ID3 Variant', fontweight='bold')
    ax4.set_ylabel('Final Accessibility Score (kcal/mol)', fontweight='bold')
    ax4.set_title('D) Final Performance Comparison', fontweight='bold', fontsize=12)
    ax4.set_xticks(x_pos)
    ax4.set_xticklabels(variants, rotation=45, ha='right')
    ax4.grid(True, alpha=0.3, axis='y')
    
    # Add horizontal line for average
    overall_mean = np.mean(variant_means)
    ax4.axhline(y=overall_mean, color='red', linestyle='--', alpha=0.5, label=f'Mean: {overall_mean:.3f}')
    ax4.legend(loc='upper right')

plt.tight_layout()

# Save the figure
output_path = '/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures/fig_convergence_analysis_cai.png'
print(f"💾 保存图表到: {output_path}")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.savefig(output_path.replace('.png', '.pdf'), dpi=300, bbox_inches='tight')
plt.close()

print("✅ Access+CAI收敛分析图表生成完成!")
print(f"🎯 输出文件: fig_convergence_analysis_cai.png/.pdf")

# Generate summary statistics
print("\n=== ACCESS+CAI收敛分析总结 ===")
print(f"📊 变体分析数量: {len(variants)}")
print("🧬 数据覆盖范围:")
for variant in variants:
    count = len(df[df['variant'] == variant])
    proteins = df[df['variant'] == variant]['protein'].nunique()
    print(f"  {variant}: {count:,} 数据点 ({proteins} 蛋白质)")

if final_performance:
    print("\n🏆 最终性能排名 (平均值 ± 标准差) - Access+CAI:")
    final_summary = final_df.groupby('Variant')['Final_Accessibility'].agg(['mean', 'std', 'count'])
    final_summary = final_summary.sort_values('mean')
    
    for variant, row in final_summary.iterrows():
        std_val = row['std'] if not pd.isna(row['std']) else 0.0
        print(f"  {variant}: {row['mean']:.3f} ± {std_val:.3f} (n={int(row['count'])})")

print("\n✨ 图表特色:")
print("- Panel A: 12个独立子图显示每个ID3变体的Access+CAI收敛轨迹")
print("- Panel B: 4个子图按基础策略分组")
print("- Panel C: 3个子图按约束机制分组")
print("- Panel D: 最终性能条形图")
print("- 颜色编码：蓝色(AMS), 紫色(CPC), 橙色(Lagrangian)")
print(f"- 性能范围: {y_min:.3f} - {y_max:.3f} kcal/mol")
print("- 反映CAI约束对accessibility性能的影响")
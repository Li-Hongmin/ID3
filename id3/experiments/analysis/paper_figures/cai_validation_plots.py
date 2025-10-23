#!/usr/bin/env python3
"""


"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

class CAIValidationPlotter:


    def __init__(self, data_dir: Path = None, output_dir: Path = None):
        """
        初始化CAI验证图表生成器

        Args:
            data_dir: 数据目录路径
            output_dir: 输出目录路径
        """
        self.data_dir = Path(data_dir) if data_dir else Path("paper_experiment_results")
        self.output_dir = Path(output_dir) if output_dir else self.data_dir / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)


        self.constraint_colors = {
            'Lagrangian': '#2E86AB',
            'Augmented Max Sampling': '#A23B72',
            'Constrained Path Continuation': '#F18F01',
        }

        self.variant_labels = {
            '00': 'Deterministic Soft',
            '01': 'Deterministic Hard',
            '10': 'Stochastic Soft',
            '11': 'Stochastic Hard',
        }


        self.full_labels = {
            'L00': 'Lagrangian Deterministic Soft',
            'L01': 'Lagrangian Deterministic Hard',
            'L10': 'Lagrangian Stochastic Soft',
            'L11': 'Lagrangian Stochastic Hard',
            'A00': 'Augmented Max Sampling Deterministic Soft',
            'A01': 'Augmented Max Sampling Deterministic Hard',
            'A10': 'Augmented Max Sampling Stochastic Soft',
            'A11': 'Augmented Max Sampling Stochastic Hard',
            'C00': 'Constrained Path Continuation Deterministic Soft',
            'C01': 'Constrained Path Continuation Deterministic Hard',
            'C10': 'Constrained Path Continuation Stochastic Soft',
            'C11': 'Constrained Path Continuation Stochastic Hard',
        }

    def load_cai_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        加载CAI相关数据

        Returns:
            (with_penalty_df, no_penalty_df): 带惩罚和不带惩罚的数据
        """
        tables_dir = self.data_dir / "tables"


        cai_comparison_path = tables_dir / "table2_cai_comparison.csv"
        if cai_comparison_path.exists():
            df_cai = pd.read_csv(cai_comparison_path)

        else:

            df_cai = None


        with_penalty_path = tables_dir / "table2_combined_with_penalty.csv"
        no_penalty_path = tables_dir / "table2_combined_no_penalty.csv"

        df_with_penalty = pd.read_csv(with_penalty_path) if with_penalty_path.exists() else None
        df_no_penalty = pd.read_csv(no_penalty_path) if no_penalty_path.exists() else None

        return df_with_penalty, df_no_penalty

    def extract_cai_values(self, df: pd.DataFrame) -> Dict[str, List[float]]:
        """
        从数据框中提取CAI值

        Args:
            df: 包含CAI数据的DataFrame

        Returns:
            按变体分组的CAI值字典
        """
        cai_data = {}


        cai_columns = [col for col in df.columns if 'CAI' in col.upper() and not 'Access' in col]

        if not cai_columns:

            return cai_data


        for col in cai_columns:

            variant = None
            for var in ['L00', 'L01', 'L10', 'L11', 'A00', 'A01', 'A10', 'A11', 'C00', 'C01', 'C10', 'C11']:
                if var in col:
                    variant = var
                    break

            if variant:
                if variant not in cai_data:
                    cai_data[variant] = []


                values = df[col].dropna()
                for val in values:
                    try:
                        if isinstance(val, str):

                            val = float(val.split('(')[0].strip())
                        else:
                            val = float(val)
                        cai_data[variant].append(val)
                    except:
                        continue

        return cai_data

    def plot_cai_distribution(self, save_name: str = "fig_cai_validation"):
        """
        生成CAI分布图表（Panel A: 各变体的最终CAI分布）

        Args:
            save_name: 保存的文件名
        """

        df_with_penalty, df_no_penalty = self.load_cai_data()


        fig, axes = plt.subplots(1, 2, figsize=(15, 6))


        if df_with_penalty is not None:
            self._plot_cai_violin(axes[0], df_with_penalty, "With Penalty")


        if df_no_penalty is not None:
            self._plot_cai_violin(axes[1], df_no_penalty, "No Penalty")


        fig.suptitle('CAI Distribution Across ID3 Variants', fontsize=16, fontweight='bold')


        plt.tight_layout()


        for fmt in ['pdf', 'png']:
            output_path = self.output_dir / f"{save_name}_distribution.{fmt}"
            fig.savefig(output_path, format=fmt, dpi=300 if fmt == 'png' else None,
                       bbox_inches='tight')


        plt.close()

    def _plot_cai_violin(self, ax, df: pd.DataFrame, title: str):
        """
        绘制CAI小提琴图

        Args:
            ax: matplotlib轴对象
            df: 数据框
            title: 标题
        """

        cai_data = self.extract_cai_values(df)

        if not cai_data:
            ax.text(0.5, 0.5, 'No CAI Data Available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(title)
            return


        plot_data = []
        labels = []
        colors = []


        constraint_groups = {
            'Lagrangian': ['L00', 'L01', 'L10', 'L11'],
            'Augmented Max Sampling': ['A00', 'A01', 'A10', 'A11'],
            'Constrained Path Continuation': ['C00', 'C01', 'C10', 'C11']
        }

        positions = []
        pos = 0
        for constraint, variants in constraint_groups.items():
            for var in variants:
                if var in cai_data and cai_data[var]:
                    plot_data.append(cai_data[var])
                    labels.append(self.full_labels[var])
                    colors.append(self.constraint_colors[constraint])
                    positions.append(pos)
                    pos += 1



        if plot_data:
            parts = ax.violinplot(plot_data, positions=positions, widths=0.8,
                                 showmeans=True, showextrema=True, showmedians=True)


            for i, pc in enumerate(parts['bodies']):
                pc.set_facecolor(colors[i])
                pc.set_alpha(0.7)


            ax.axhline(y=0.8, color='red', linestyle='--', linewidth=2, alpha=0.7,
                      label='Target CAI (0.8)')


            ax.set_xticks(positions)
            ax.set_xticklabels([l.split()[-2] + ' ' + l.split()[-1] for l in labels],
                              rotation=45, ha='right')
            ax.set_ylabel('CAI Value', fontsize=12)
            ax.set_ylim([0.75, 0.85])
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right')

        ax.set_title(title, fontsize=14, fontweight='bold')

    def plot_cai_improvement(self, save_name: str = "fig_cai_improvement"):
        """
        生成CAI改进幅度图表（Panel B: CAI改进幅度）

        Args:
            save_name: 保存的文件名
        """

        df_with_penalty, df_no_penalty = self.load_cai_data()


        fig, ax = plt.subplots(1, 1, figsize=(12, 6))


        cai_with = self.extract_cai_values(df_with_penalty) if df_with_penalty is not None else {}
        cai_without = self.extract_cai_values(df_no_penalty) if df_no_penalty is not None else {}


        improvements = []
        labels = []
        colors = []

        for variant in ['L00', 'L01', 'L10', 'L11', 'A00', 'A01', 'A10', 'A11', 'C00', 'C01', 'C10', 'C11']:
            if variant in cai_with and variant in cai_without:
                with_mean = np.mean(cai_with[variant]) if cai_with[variant] else 0
                without_mean = np.mean(cai_without[variant]) if cai_without[variant] else 0


                if without_mean > 0:
                    improvement = ((with_mean - without_mean) / without_mean) * 100
                    improvements.append(improvement)
                    labels.append(self.full_labels[variant])


                    if variant.startswith('L'):
                        colors.append(self.constraint_colors['Lagrangian'])
                    elif variant.startswith('A'):
                        colors.append(self.constraint_colors['Augmented Max Sampling'])
                    else:
                        colors.append(self.constraint_colors['Constrained Path Continuation'])


        if improvements:
            x_pos = np.arange(len(labels))
            bars = ax.bar(x_pos, improvements, color=colors, alpha=0.7, edgecolor='black', linewidth=1)


            for bar, imp in zip(bars, improvements):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{imp:.1f}%', ha='center', va='bottom' if height >= 0 else 'top',
                       fontsize=9)


            ax.set_xticks(x_pos)
            ax.set_xticklabels([l.split()[-2] + ' ' + l.split()[-1] for l in labels],
                              rotation=45, ha='right')
            ax.set_ylabel('CAI Improvement (%)', fontsize=12)
            ax.set_title('CAI Improvement with Penalty vs No Penalty', fontsize=14, fontweight='bold')
            ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            ax.grid(True, alpha=0.3, axis='y')


        plt.tight_layout()
        for fmt in ['pdf', 'png']:
            output_path = self.output_dir / f"{save_name}.{fmt}"
            fig.savefig(output_path, format=fmt, dpi=300 if fmt == 'png' else None,
                       bbox_inches='tight')


        plt.close()

    def generate_all_plots(self):
        """生成所有CAI验证相关图表"""
        logger.info("生成CAI验证图表...")


        self.plot_cai_distribution()


        self.plot_cai_improvement()

        logger.info("CAI验证图表生成完成")
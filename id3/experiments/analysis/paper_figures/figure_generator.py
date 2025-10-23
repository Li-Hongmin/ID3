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
from typing import Dict, List, Optional, Tuple, Any


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PaperFigureGenerator:


    def __init__(self, data_dir: Path = None, output_dir: Path = None):
        """
        初始化图表生成器

        Args:
            data_dir: 数据目录路径
            output_dir: 输出目录路径
        """
        self.data_dir = Path(data_dir) if data_dir else Path("paper_experiment_results")
        self.output_dir = Path(output_dir) if output_dir else self.data_dir / "figures"
        self.output_dir.mkdir(parents=True, exist_ok=True)


        self._setup_style()


        self.colors = {















        }


        self.labels = {

            'lagrangian': 'Lagrangian',
            'ams': 'Augmented Max Sampling',
            'cpc': 'Constrained Path Continuation',


            '00': 'Deterministic Soft',
            '01': 'Deterministic Hard',
            '10': 'Stochastic Soft',
            '11': 'Stochastic Hard',


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

    def _setup_style(self):
        """设置统一的图表样式"""

        plt.style.use('seaborn-v0_8-whitegrid')


        plt.rcParams.update({
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 11,
            'ytick.labelsize': 11,
            'legend.fontsize': 11,
            'figure.titlesize': 18,
            'font.family': 'sans-serif',
            'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
            'axes.unicode_minus': False,
            'figure.dpi': 100,
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.format': 'pdf',
        })


        sns.set_context("paper", font_scale=1.2)
        sns.set_style("whitegrid", {
            'axes.grid': True,
            'grid.alpha': 0.3,
            'axes.edgecolor': '0.2',
            'axes.linewidth': 1.2
        })

    def load_table_data(self) -> Tuple[pd.DataFrame, ...]:
        """


        Returns:

        """
        tables_dir = self.data_dir / "tables"


        access_path = tables_dir / "table1_accessibility_performance.csv"
        df_access = pd.read_csv(access_path) if access_path.exists() else None


        cai_path = tables_dir / "table2_cai_comparison.csv"
        df_cai = pd.read_csv(cai_path) if cai_path.exists() else None


        combined_with_penalty = tables_dir / "table2_combined_with_penalty.csv"
        df_with_penalty = pd.read_csv(combined_with_penalty) if combined_with_penalty.exists() else None

        combined_no_penalty = tables_dir / "table2_combined_no_penalty.csv"
        df_no_penalty = pd.read_csv(combined_no_penalty) if combined_no_penalty.exists() else None

        logger.info(f"加载的数据文件:")
        if df_access is not None:
            logger.info(f"  - Access-only: {df_access.shape}")
        if df_cai is not None:
            logger.info(f"  - CAI comparison: {df_cai.shape}")
        if df_with_penalty is not None:
            logger.info(f"  - With penalty: {df_with_penalty.shape}")
        if df_no_penalty is not None:
            logger.info(f"  - No penalty: {df_no_penalty.shape}")

        return df_access, df_cai, df_with_penalty, df_no_penalty

    def load_experiment_json(self, experiment_type: str, protein: str = None,
                           constraint: str = None, variant: str = None) -> List[Dict]:
        """


        Args:





        Returns:

        """
        exp_dir = self.data_dir / experiment_type
        results = []

        if not exp_dir.exists():
            logger.warning(f"实验目录不存在: {exp_dir}")
            return results


        pattern = "*"
        if protein:
            pattern += f"*{protein}*"
        if constraint:
            pattern += f"*{constraint}*"
        if variant:
            pattern += f"*{variant}*"
        pattern += "*.json"


        for json_file in exp_dir.glob(pattern):
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    results.append(data)
            except Exception as e:
                logger.warning(f"加载文件失败 {json_file}: {e}")

        logger.info(f"从 {experiment_type} 加载了 {len(results)} 个实验")
        return results

    def save_figure(self, fig: plt.Figure, filename: str, formats: List[str] = None):
        """


        Args:



        """
        if formats is None:
            formats = ['pdf', 'png']

        for fmt in formats:
            output_path = self.output_dir / f"{filename}.{fmt}"
            fig.savefig(output_path, format=fmt, dpi=300 if fmt == 'png' else None,
                       bbox_inches='tight', pad_inches=0.1)
            logger.info(f"图表已保存: {output_path}")

    def generate_all_figures(self):




        df_access, df_cai, df_with_penalty, df_no_penalty = self.load_table_data()


        from .cai_validation_plots import CAIValidationPlotter
        from .p00004_case_study import P00004CaseStudyPlotter


        cai_plotter = CAIValidationPlotter(self.data_dir, self.output_dir)
        cai_plotter.generate_all_plots()


        p00004_plotter = P00004CaseStudyPlotter(self.data_dir, self.output_dir)
        p00004_plotter.generate_all_plots()



    def get_variant_label(self, variant_code: str) -> str:
        """
        获取变体的完整标签（不使用缩写）

        Args:
            variant_code: 变体代码（如'L00', 'A01'等）

        Returns:
            完整的标签文本
        """
        return self.labels.get(variant_code, variant_code)
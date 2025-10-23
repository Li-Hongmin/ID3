#!/usr/bin/env python3
"""


"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

def parse_mean_std(value_str):

    if value_str == "N/A":
        return np.nan, np.nan

    if '±' in value_str:
        parts = value_str.split('±')
        return float(parts[0]), float(parts[1])
    else:
        return float(value_str), 0.0

def create_separate_tables(input_csv: Path, output_dir: Path = None):
    """
    从包含mean±std的CSV文件生成分离的均值表和标准差表
    """
    if output_dir is None:
        output_dir = input_csv.parent

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


    df = pd.read_csv(input_csv)


    mean_df = pd.DataFrame()
    std_df = pd.DataFrame()

    mean_df['Protein'] = df['Protein']
    std_df['Protein'] = df['Protein']


    for col in df.columns:
        if col == 'Protein':
            continue

        means = []
        stds = []

        for value in df[col]:
            if pd.isna(value):
                means.append(np.nan)
                stds.append(np.nan)
            else:
                mean_val, std_val = parse_mean_std(str(value))
                means.append(mean_val)
                stds.append(std_val)

        mean_df[col] = means
        std_df[col] = stds


    base_name = input_csv.stem.replace('_with_std', '')

    mean_output = output_dir / f"{base_name}_means.csv"
    std_output = output_dir / f"{base_name}_stds.csv"

    mean_df.to_csv(mean_output, index=False, float_format='%.4f')
    std_df.to_csv(std_output, index=False, float_format='%.4f')









    numeric_cols = [col for col in mean_df.columns if col != 'Protein']
    all_means = mean_df[numeric_cols].values.flatten()
    all_stds = std_df[numeric_cols].values.flatten()


    all_means = all_means[~np.isnan(all_means)]
    all_stds = all_stds[~np.isnan(all_stds)]













    return mean_df, std_df

def main():
    if len(sys.argv) < 2:


        return 1

    input_csv = Path(sys.argv[1])
    if not input_csv.exists():

        return 1

    output_dir = None
    if len(sys.argv) > 2:
        output_dir = Path(sys.argv[2])

    mean_df, std_df = create_separate_tables(input_csv, output_dir)


    return 0

if __name__ == "__main__":
    sys.exit(main())
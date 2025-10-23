#!/usr/bin/env python3
"""

"""
import sys
from pathlib import Path


sys.path.append(str(Path(__file__).parent))


print("=" * 60)
print("🚀 测试双重进度条 - 串行执行")
print("=" * 60)
print("\n运行2个实验，每个100次迭代，展示双重进度条效果\n")

import subprocess


cmd = [
    "python", "run_unified_experiment.py",
    "--proteins", "O15263",
    "--constraints", "lagrangian",
    "--variants", "11",
    "--iterations", "100",
    "--seeds", "2",
    "--enable-cai",
    "--parallel", "1",
    "--force"
]

print(f"执行命令: {' '.join(cmd)}\n")
subprocess.run(cmd)
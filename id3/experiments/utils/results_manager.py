#!/usr/bin/env python3
"""
结果管理模块

负责保存、加载和分析实验结果
"""

import json
import time
from pathlib import Path
from typing import Dict, Any, Optional
import torch
import numpy as np


class ResultsManager:
    """实验结果管理器"""
    
    def __init__(self, base_dir: str = "experiments/results"):
        """
        初始化结果管理器
        
        Args:
            base_dir: 结果保存基础目录
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
    
    def create_experiment_dir(self, experiment_name: str) -> Path:
        """
        创建实验目录
        
        Args:
            experiment_name: 实验名称
            
        Returns:
            实验目录路径
        """
        exp_dir = self.base_dir / experiment_name
        exp_dir.mkdir(parents=True, exist_ok=True)
        return exp_dir
    
    def save_results(
        self,
        results: Dict[str, Any],
        experiment_name: str,
        protein_id: str,
        suffix: Optional[str] = None
    ) -> Path:
        """
        保存实验结果
        
        Args:
            results: 结果数据
            experiment_name: 实验名称
            protein_id: 蛋白质ID
            suffix: 文件名后缀
            
        Returns:
            保存的文件路径
        """
        # 创建实验目录
        exp_dir = self.create_experiment_dir(experiment_name)
        
        # 构建文件名
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        if suffix:
            filename = f"{protein_id}_{suffix}_{timestamp}.json"
        else:
            filename = f"{protein_id}_{timestamp}.json"
        
        filepath = exp_dir / filename
        
        # 转换数据类型
        results_json = self._prepare_for_json(results)
        
        # 保存
        with open(filepath, 'w') as f:
            json.dump(results_json, f, indent=2)
        
        print(f"结果已保存: {filepath}")
        return filepath
    
    def load_results(self, filepath: str) -> Dict[str, Any]:
        """
        加载实验结果
        
        Args:
            filepath: 结果文件路径
            
        Returns:
            结果数据
        """
        with open(filepath, 'r') as f:
            return json.load(f)
    
    def _prepare_for_json(self, obj: Any) -> Any:
        """
        准备数据用于JSON序列化
        
        Args:
            obj: 要序列化的对象
            
        Returns:
            可序列化的对象
        """
        if isinstance(obj, dict):
            return {k: self._prepare_for_json(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._prepare_for_json(v) for v in obj]
        elif isinstance(obj, torch.Tensor):
            return obj.cpu().numpy().tolist()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif hasattr(obj, 'item'):
            return obj.item()
        else:
            return obj
    
    def compare_results(
        self,
        result1: Dict[str, Any],
        result2: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        比较两个实验结果
        
        Args:
            result1: 第一个结果
            result2: 第二个结果
            
        Returns:
            比较结果
        """
        comparison = {}
        
        # 比较关键指标
        metrics = [
            'initial_accessibility',
            'final_accessibility',
            'improvement',
            'optimization_time',
            'amino_acids_correct'
        ]
        
        for metric in metrics:
            if metric in result1 and metric in result2:
                val1 = result1[metric]
                val2 = result2[metric]
                
                if isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                    comparison[metric] = {
                        'result1': val1,
                        'result2': val2,
                        'difference': val1 - val2,
                        'relative_diff': (val1 - val2) / val2 if val2 != 0 else None
                    }
                else:
                    comparison[metric] = {
                        'result1': val1,
                        'result2': val2,
                        'equal': val1 == val2
                    }
        
        return comparison
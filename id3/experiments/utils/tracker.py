#!/usr/bin/env python3
"""
优化跟踪器模块

负责记录和管理优化过程中的指标
"""

import time
from typing import Dict, Optional, List, Any
from collections import defaultdict


class OptimizationTracker:
    """优化过程跟踪器"""
    
    def __init__(
        self, 
        record_interval: int = 10,
        verbose_interval: int = 100
    ):
        """
        初始化跟踪器
        
        Args:
            record_interval: 数据记录间隔
            verbose_interval: 打印输出间隔
        """
        self.record_interval = record_interval
        self.verbose_interval = verbose_interval
        self.data = defaultdict(list)
        self.start_time = None
        self.iteration_count = 0
    
    def start(self):
        """开始跟踪"""
        self.start_time = time.time()
        self.data.clear()
        self.iteration_count = 0
    
    def record(self, iteration: int, **metrics):
        """
        记录一次迭代的数据
        
        Args:
            iteration: 迭代次数
            **metrics: 任意指标键值对
        """
        self.iteration_count = iteration
        
        if iteration % self.record_interval == 0:
            self.data['iterations'].append(iteration)
            self.data['timestamps'].append(time.time() - self.start_time)
            
            for key, value in metrics.items():
                if value is not None:
                    # 处理tensor类型
                    if hasattr(value, 'item'):
                        value = value.item()
                    self.data[key].append(value)
    
    def print_progress(self, iteration: int, **metrics):
        """
        打印进度信息
        
        Args:
            iteration: 迭代次数
            **metrics: 要打印的指标
        """
        if iteration % self.verbose_interval == 0:
            # 构建消息
            parts = [f"Iter {iteration:4d}"]
            
            for key, value in metrics.items():
                if value is not None:
                    # 处理不同类型的值
                    if hasattr(value, 'item'):
                        value = value.item()
                    
                    if isinstance(value, float):
                        parts.append(f"{key}={value:.4f}")
                    elif isinstance(value, bool):
                        parts.append(f"{key}={value}")
                    else:
                        parts.append(f"{key}={value}")
            
            print(" | ".join(parts))
    
    def get_results(self) -> Dict[str, List]:
        """获取跟踪结果"""
        return dict(self.data)
    
    def get_summary(self) -> Dict[str, Any]:
        """获取摘要统计"""
        elapsed_time = time.time() - self.start_time if self.start_time else 0
        
        summary = {
            'total_iterations': self.iteration_count,
            'elapsed_time': elapsed_time,
            'iterations_per_second': self.iteration_count / elapsed_time if elapsed_time > 0 else 0
        }
        
        # 添加各指标的最终值
        for key, values in self.data.items():
            if key not in ['iterations', 'timestamps'] and values:
                summary[f'final_{key}'] = values[-1]
                summary[f'initial_{key}'] = values[0]
                
                # 对于数值类型，计算改善
                if isinstance(values[0], (int, float)):
                    summary[f'{key}_improvement'] = values[0] - values[-1]
        
        return summary


class ProgressTracker:
    """Simple progress tracker for experiment batches"""
    
    def __init__(self, total: int):
        """
        Initialize progress tracker
        
        Args:
            total: Total number of items to track
        """
        self.total = total
        self.current = 0
        self.start_time = time.time()
    
    def update(self, description: str = ""):
        """Update progress"""
        self.current += 1
        elapsed = time.time() - self.start_time
        rate = self.current / elapsed if elapsed > 0 else 0
        
        print(f"   Progress: {self.current}/{self.total} ({self.current/self.total*100:.1f}%) - "
              f"{rate:.2f} items/s - {description}")
    
    def reset(self):
        """Reset tracker"""
        self.current = 0
        self.start_time = time.time()
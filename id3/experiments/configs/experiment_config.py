#!/usr/bin/env python3
"""
实验配置
"""

class ExperimentConfig:
    """12×12实验配置"""
    
    # 12个蛋白质
    PROTEINS = [
        'P00004', 'P01308', 'P01825', 'P04637',
        'P0CG48', 'P0DTC2', 'P0DTC9', 'P31417',
        'P42212', 'P61626', 'P99999', 'O15263'
    ]
    
    # 12个变体（3种约束 × 4种参数组合）
    VARIANTS = [
        # Lagrangian
        'lagrangian-00', 'lagrangian-01', 'lagrangian-10', 'lagrangian-11',
        # CPC  
        'cpc-00', 'cpc-01', 'cpc-10', 'cpc-11',
        # AMS
        'ams-00', 'ams-01', 'ams-10', 'ams-11'
    ]
    
    # 默认参数
    DEFAULT_ITERATIONS = 1000
    DEFAULT_MAX_PARALLEL = 4
    
    @classmethod
    def generate_tasks(cls, num_iterations: int = None):
        """生成所有实验任务"""
        if num_iterations is None:
            num_iterations = cls.DEFAULT_ITERATIONS
            
        tasks = []
        for protein in cls.PROTEINS:
            for variant in cls.VARIANTS:
                tasks.append((protein, variant, num_iterations))
        return tasks
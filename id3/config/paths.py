#!/usr/bin/env python3
"""
路径配置模块

集中管理所有项目路径，避免硬编码
"""

import os
from pathlib import Path


class PathConfig:
    """路径配置管理器"""
    
    # 项目根目录
    PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
    
    # 数据目录
    DATA_DIR = PROJECT_ROOT / "data"
    UTR_TEMPLATES_DIR = DATA_DIR / "utr_templates"
    
    # DeepRaccess模型路径
    DEEPRACCESS_DIR = PROJECT_ROOT / "DeepRaccess"
    DEEPRACCESS_MODEL_PATH = DEEPRACCESS_DIR / "path" / "FCN_structured.pth"
    
    # 实验结果目录
    EXPERIMENTS_DIR = PROJECT_ROOT / "experiments"
    RESULTS_DIR = EXPERIMENTS_DIR / "results"
    
    # 存档目录
    ARCHIVES_DIR = PROJECT_ROOT / "archives"
    
    @classmethod
    def get_deepraccess_model_path(cls) -> str:
        """获取DeepRaccess模型路径"""
        return str(cls.DEEPRACCESS_MODEL_PATH)
    
    @classmethod
    def get_data_dir(cls) -> str:
        """获取数据目录路径"""
        return str(cls.DATA_DIR)
    
    @classmethod
    def get_results_dir(cls) -> str:
        """获取结果目录路径"""
        return str(cls.RESULTS_DIR)
    
    @classmethod
    def ensure_dirs_exist(cls):
        """确保所有必要的目录存在"""
        dirs = [
            cls.DATA_DIR,
            cls.RESULTS_DIR,
            cls.ARCHIVES_DIR,
            cls.UTR_TEMPLATES_DIR
        ]
        for dir_path in dirs:
            dir_path.mkdir(parents=True, exist_ok=True)


# 便利函数
def get_deepraccess_model_path() -> str:
    """获取DeepRaccess模型路径的便利函数"""
    return PathConfig.get_deepraccess_model_path()


def get_data_dir() -> str:
    """获取数据目录路径的便利函数"""
    return PathConfig.get_data_dir()


def get_results_dir() -> str:
    """获取结果目录路径的便利函数"""
    return PathConfig.get_results_dir()
#!/usr/bin/env python3
"""
Path configuration module

Centralized management of all project paths, avoiding hardcoding
"""

import os
from pathlib import Path


class PathConfig:
    """Path configuration manager"""

    # Project root directory
    PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

    # Data directory
    DATA_DIR = PROJECT_ROOT / "data"
    UTR_TEMPLATES_DIR = DATA_DIR / "utr_templates"

    # DeepRaccess model path
    DEEPRACCESS_DIR = PROJECT_ROOT / "DeepRaccess"
    DEEPRACCESS_MODEL_PATH = DEEPRACCESS_DIR / "path" / "FCN_structured.pth"

    # Experiment results directory
    EXPERIMENTS_DIR = PROJECT_ROOT / "experiments"
    RESULTS_DIR = EXPERIMENTS_DIR / "results"

    # Archive directory
    ARCHIVES_DIR = PROJECT_ROOT / "archives"
    
    @classmethod
    def get_deepraccess_model_path(cls) -> str:
        """Get DeepRaccess model path"""
        return str(cls.DEEPRACCESS_MODEL_PATH)

    @classmethod
    def get_data_dir(cls) -> str:
        """Get data directory path"""
        return str(cls.DATA_DIR)

    @classmethod
    def get_results_dir(cls) -> str:
        """Get results directory path"""
        return str(cls.RESULTS_DIR)

    @classmethod
    def ensure_dirs_exist(cls):
        """Ensure all necessary directories exist"""
        dirs = [
            cls.DATA_DIR,
            cls.RESULTS_DIR,
            cls.ARCHIVES_DIR,
            cls.UTR_TEMPLATES_DIR
        ]
        for dir_path in dirs:
            dir_path.mkdir(parents=True, exist_ok=True)


# Convenience functions
def get_deepraccess_model_path() -> str:
    """Convenience function to get DeepRaccess model path"""
    return PathConfig.get_deepraccess_model_path()


def get_data_dir() -> str:
    """Convenience function to get data directory path"""
    return PathConfig.get_data_dir()


def get_results_dir() -> str:
    """Convenience function to get results directory path"""
    return PathConfig.get_results_dir()
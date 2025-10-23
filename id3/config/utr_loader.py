#!/usr/bin/env python3
"""
UTR模板加载器

统一管理UTR序列模板，从data/utr_templates/目录读取
支持缓存和单例模式，提供灵活的UTR序列管理接口
"""

import os
from typing import Dict, Optional
from pathlib import Path
import threading


class UTRLoader:
    """统一的UTR模板加载器
    
    使用单例模式确保全局只有一个实例
    支持从文件加载UTR模板，并提供缓存机制
    """
    
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self):
        """初始化UTR加载器"""
        if UTRLoader._instance is not None:
            raise RuntimeError("UTRLoader是单例类，请使用get_instance()方法获取实例")
        
        self._cache = {}
        self._loaded = False
        self._project_root = self._find_project_root()
    
    @classmethod
    def get_instance(cls) -> 'UTRLoader':
        """获取UTRLoader单例实例"""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance
    
    def _find_project_root(self) -> Path:
        """查找项目根目录"""
        current_dir = Path(__file__).resolve().parent
        
        # 向上查找直到找到包含data/utr_templates的目录
        while current_dir.parent != current_dir:
            utr_templates_path = current_dir / "data" / "utr_templates"
            if utr_templates_path.exists():
                return current_dir
            current_dir = current_dir.parent
        
        # 如果没找到，使用相对路径
        return Path(__file__).resolve().parent.parent.parent
    
    def _parse_fasta_like_file(self, file_path: Path) -> str:
        """解析FASTA格式的UTR文件
        
        Args:
            file_path: UTR模板文件路径
            
        Returns:
            str: UTR序列
        """
        if not file_path.exists():
            raise FileNotFoundError(f"UTR模板文件不存在: {file_path}")
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        sequence_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('>'):
                # 移除空格和换行符，只保留核苷酸序列
                sequence_lines.append(line.replace(' ', '').replace('\t', ''))
        
        if not sequence_lines:
            raise ValueError(f"在文件 {file_path} 中未找到有效的序列")
        
        # 合并所有序列行
        sequence = ''.join(sequence_lines)
        
        # 验证序列只包含有效的核苷酸
        valid_nucleotides = set('ATCGU')
        if not set(sequence.upper()).issubset(valid_nucleotides):
            invalid_chars = set(sequence.upper()) - valid_nucleotides
            raise ValueError(f"序列包含无效字符: {invalid_chars}")
        
        return sequence.upper()
    
    def load_utr_templates(self, 
                          utr5_file: Optional[str] = None,
                          utr3_file: Optional[str] = None,
                          force_reload: bool = False) -> Dict[str, str]:
        """从文件加载UTR模板
        
        Args:
            utr5_file: 5'UTR文件路径，默认为data/utr_templates/5utr_templates.txt
            utr3_file: 3'UTR文件路径，默认为data/utr_templates/3utr_templates.txt
            force_reload: 是否强制重新加载，忽略缓存
            
        Returns:
            Dict[str, str]: 包含UTR序列的字典
        """
        if self._loaded and not force_reload:
            return self._cache
        
        # 设置默认文件路径
        if utr5_file is None:
            utr5_file = self._project_root / "data" / "utr_templates" / "5utr_templates.txt"
        else:
            utr5_file = Path(utr5_file)
            if not utr5_file.is_absolute():
                utr5_file = self._project_root / utr5_file
        
        if utr3_file is None:
            utr3_file = self._project_root / "data" / "utr_templates" / "3utr_templates.txt"
        else:
            utr3_file = Path(utr3_file)
            if not utr3_file.is_absolute():
                utr3_file = self._project_root / utr3_file
        
        try:
            # 加载5'UTR
            utr5_sequence = self._parse_fasta_like_file(utr5_file)
            
            # 加载3'UTR
            utr3_sequence = self._parse_fasta_like_file(utr3_file)
            
            # 更新缓存
            self._cache = {
                'utr5_default': utr5_sequence,
                'utr3_default': utr3_sequence,
                'utr5_file': str(utr5_file),
                'utr3_file': str(utr3_file)
            }
            
            self._loaded = True
            
            return self._cache
            
        except Exception as e:
            raise RuntimeError(f"加载UTR模板失败: {e}")
    
    def get_default_utr5(self) -> str:
        """获取默认5'UTR序列
        
        Returns:
            str: 5'UTR序列
        """
        if not self._loaded:
            self.load_utr_templates()
        
        return self._cache.get('utr5_default', '')
    
    def get_default_utr3(self) -> str:
        """获取默认3'UTR序列
        
        Returns:
            str: 3'UTR序列
        """
        if not self._loaded:
            self.load_utr_templates()
        
        return self._cache.get('utr3_default', '')
    
    def get_utr_info(self) -> Dict[str, any]:
        """获取UTR信息
        
        Returns:
            Dict[str, any]: 包含UTR序列长度和来源文件的信息
        """
        if not self._loaded:
            self.load_utr_templates()
        
        utr5 = self._cache.get('utr5_default', '')
        utr3 = self._cache.get('utr3_default', '')
        
        return {
            'utr5_length': len(utr5),
            'utr3_length': len(utr3),
            'utr5_file': self._cache.get('utr5_file', ''),
            'utr3_file': self._cache.get('utr3_file', ''),
            'total_utr_length': len(utr5) + len(utr3)
        }
    
    def validate_utr_sequences(self) -> bool:
        """验证UTR序列的有效性
        
        Returns:
            bool: 序列是否有效
        """
        try:
            utr5 = self.get_default_utr5()
            utr3 = self.get_default_utr3()
            
            # 检查序列是否为空
            if not utr5 or not utr3:
                return False
            
            # 检查序列长度是否合理
            if len(utr5) < 10 or len(utr5) > 200:
                return False
            
            if len(utr3) < 5 or len(utr3) > 200:
                return False
            
            # 检查序列是否只包含有效核苷酸
            valid_nucleotides = set('ATCGU')
            if not set(utr5).issubset(valid_nucleotides) or not set(utr3).issubset(valid_nucleotides):
                return False
            
            return True
            
        except Exception:
            return False


# 便利函数
def get_utr_loader() -> UTRLoader:
    """获取UTR加载器实例的便利函数"""
    return UTRLoader.get_instance()


def get_default_utrs() -> Dict[str, str]:
    """获取默认UTR序列的便利函数
    
    Returns:
        Dict[str, str]: 包含utr5和utr3的字典
    """
    loader = get_utr_loader()
    return {
        'utr5': loader.get_default_utr5(),
        'utr3': loader.get_default_utr3()
    }


if __name__ == "__main__":
    # 测试代码
    loader = UTRLoader.get_instance()
    
    try:
        loader.load_utr_templates()
        
        print("UTR加载器测试:")
        print(f"5'UTR长度: {len(loader.get_default_utr5())}")
        print(f"3'UTR长度: {len(loader.get_default_utr3())}")
        print(f"5'UTR序列: {loader.get_default_utr5()}")
        print(f"3'UTR序列: {loader.get_default_utr3()}")
        
        info = loader.get_utr_info()
        print(f"UTR信息: {info}")
        
        is_valid = loader.validate_utr_sequences()
        print(f"序列验证: {'通过' if is_valid else '失败'}")
        
    except Exception as e:
        print(f"测试失败: {e}")
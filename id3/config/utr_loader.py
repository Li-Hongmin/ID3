#!/usr/bin/env python3
"""
UTR template loader

Unified management of UTR sequence templates, loaded from data/utr_templates/ directory
Supports caching and singleton pattern, providing flexible UTR sequence management interface
"""

import os
from typing import Dict, Optional
from pathlib import Path
import threading


class UTRLoader:
    """Unified UTR template loader

    Uses singleton pattern to ensure only one global instance
    Supports loading UTR templates from files with caching mechanism
    """
    
    _instance = None
    _lock = threading.Lock()
    
    def __init__(self):
        """Initialize UTR loader"""
        if UTRLoader._instance is not None:
            raise RuntimeError("UTRLoader is a singleton class, please use get_instance() method to get instance")
        
        self._cache = {}
        self._loaded = False
        self._project_root = self._find_project_root()
    
    @classmethod
    def get_instance(cls) -> 'UTRLoader':
        """Get UTRLoader singleton instance"""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance
    
    def _find_project_root(self) -> Path:
        """Find project root directory"""
        current_dir = Path(__file__).resolve().parent

        # Search upward until finding directory containing data/utr_templates
        while current_dir.parent != current_dir:
            utr_templates_path = current_dir / "data" / "utr_templates"
            if utr_templates_path.exists():
                return current_dir
            current_dir = current_dir.parent

        # If not found, use relative path
        return Path(__file__).resolve().parent.parent.parent
    
    def _parse_fasta_like_file(self, file_path: Path) -> str:
        """Parse FASTA format UTR file

        Args:
            file_path: UTR template file path

        Returns:
            str: UTR sequence
        """
        if not file_path.exists():
            raise FileNotFoundError(f"UTR template file does not exist: {file_path}")
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        sequence_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('>'):
                # Remove spaces and tabs, keep only nucleotide sequence
                sequence_lines.append(line.replace(' ', '').replace('\t', ''))

        if not sequence_lines:
            raise ValueError(f"No valid sequence found in file {file_path}")

        # Merge all sequence lines
        sequence = ''.join(sequence_lines)

        # Validate sequence contains only valid nucleotides
        valid_nucleotides = set('ATCGU')
        if not set(sequence.upper()).issubset(valid_nucleotides):
            invalid_chars = set(sequence.upper()) - valid_nucleotides
            raise ValueError(f"Sequence contains invalid characters: {invalid_chars}")
        
        return sequence.upper()
    
    def load_utr_templates(self,
                          utr5_file: Optional[str] = None,
                          utr3_file: Optional[str] = None,
                          force_reload: bool = False) -> Dict[str, str]:
        """Load UTR templates from files

        Args:
            utr5_file: 5'UTR file path, defaults to data/utr_templates/5utr_templates.txt
            utr3_file: 3'UTR file path, defaults to data/utr_templates/3utr_templates.txt
            force_reload: Whether to force reload, ignoring cache

        Returns:
            Dict[str, str]: Dictionary containing UTR sequences
        """
        if self._loaded and not force_reload:
            return self._cache

        # Set default file paths
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
            # Load 5'UTR
            utr5_sequence = self._parse_fasta_like_file(utr5_file)

            # Load 3'UTR
            utr3_sequence = self._parse_fasta_like_file(utr3_file)

            # Update cache
            self._cache = {
                'utr5_default': utr5_sequence,
                'utr3_default': utr3_sequence,
                'utr5_file': str(utr5_file),
                'utr3_file': str(utr3_file)
            }

            self._loaded = True

            return self._cache

        except Exception as e:
            raise RuntimeError(f"Failed to load UTR templates: {e}")
    
    def get_default_utr5(self) -> str:
        """Get default 5'UTR sequence

        Returns:
            str: 5'UTR sequence
        """
        if not self._loaded:
            self.load_utr_templates()
        
        return self._cache.get('utr5_default', '')
    
    def get_default_utr3(self) -> str:
        """Get default 3'UTR sequence

        Returns:
            str: 3'UTR sequence
        """
        if not self._loaded:
            self.load_utr_templates()
        
        return self._cache.get('utr3_default', '')
    
    def get_utr_info(self) -> Dict[str, any]:
        """Get UTR information

        Returns:
            Dict[str, any]: Information containing UTR sequence lengths and source files
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
        """Validate the validity of UTR sequences

        Returns:
            bool: Whether the sequences are valid
        """
        try:
            utr5 = self.get_default_utr5()
            utr3 = self.get_default_utr3()

            # Check if sequences are empty
            if not utr5 or not utr3:
                return False

            # Check if sequence lengths are reasonable
            if len(utr5) < 10 or len(utr5) > 200:
                return False

            if len(utr3) < 5 or len(utr3) > 200:
                return False

            # Check if sequences contain only valid nucleotides
            valid_nucleotides = set('ATCGU')
            if not set(utr5).issubset(valid_nucleotides) or not set(utr3).issubset(valid_nucleotides):
                return False

            return True

        except Exception:
            return False


# Convenience functions
def get_utr_loader() -> UTRLoader:
    """Convenience function to get UTR loader instance"""
    return UTRLoader.get_instance()


def get_default_utrs() -> Dict[str, str]:
    """Convenience function to get default UTR sequences

    Returns:
        Dict[str, str]: Dictionary containing utr5 and utr3
    """
    loader = get_utr_loader()
    return {
        'utr5': loader.get_default_utr5(),
        'utr3': loader.get_default_utr3()
    }


if __name__ == "__main__":
    # Test code
    loader = UTRLoader.get_instance()

    try:
        loader.load_utr_templates()

        print("UTR loader test:")
        print(f"5'UTR length: {len(loader.get_default_utr5())}")
        print(f"3'UTR length: {len(loader.get_default_utr3())}")
        print(f"5'UTR sequence: {loader.get_default_utr5()}")
        print(f"3'UTR sequence: {loader.get_default_utr3()}")

        info = loader.get_utr_info()
        print(f"UTR information: {info}")

        is_valid = loader.validate_utr_sequences()
        print(f"Sequence validation: {'Passed' if is_valid else 'Failed'}")

    except Exception as e:
        print(f"Test failed: {e}")
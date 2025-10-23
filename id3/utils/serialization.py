#!/usr/bin/env python3
"""
Serialization Module

Handles serialization and deserialization of ID3 data structures including
numpy arrays, PyTorch tensors, and step records with metadata.
"""

import numpy as np
import torch
import json
from typing import Any, Dict, List, Union
from pathlib import Path
from datetime import datetime


def serialize_data(obj: Any) -> Any:
    """
    Serialize Python objects to JSON-compatible format.

    Args:
        obj: Object to serialize (supports numpy arrays, PyTorch tensors, dicts, lists)

    Returns:
        serialized_obj: JSON-serializable representation
    """
    if isinstance(obj, np.ndarray):
        return {
            '__type__': 'ndarray',
            'data': obj.tolist(),
            'shape': obj.shape,
            'dtype': str(obj.dtype)
        }
    elif isinstance(obj, torch.Tensor):
        return {
            '__type__': 'tensor',
            'data': obj.detach().cpu().numpy().tolist(),
            'shape': list(obj.shape),
            'dtype': str(obj.dtype),
            'device': str(obj.device)
        }
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, torch.dtype):
        return str(obj)
    elif isinstance(obj, dict):
        return {key: serialize_data(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [serialize_data(item) for item in obj]
    else:
        return obj


def deserialize_data(obj: Any) -> Any:
    """
    Deserialize JSON-compatible format back to Python objects.

    Args:
        obj: Serialized object to deserialize

    Returns:
        deserialized_obj: Original Python object (numpy array, PyTorch tensor, etc.)
    """
    if isinstance(obj, dict):
        if obj.get('__type__') == 'ndarray':
            return np.array(obj['data'], dtype=obj['dtype']).reshape(obj['shape'])
        elif obj.get('__type__') == 'tensor':
            arr = np.array(obj['data'], dtype=obj['dtype']).reshape(obj['shape'])
            return torch.from_numpy(arr)
        else:
            return {key: deserialize_data(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [deserialize_data(item) for item in obj]
    else:
        return obj


def save_step_records_to_json(
    step_records: List[Dict],
    metadata: Dict,
    output_file: Union[str, Path],
    compress: bool = False
) -> Path:
    """
    Save step records and metadata to JSON file (with optional compression).

    Args:
        step_records: List of step record dictionaries
        metadata: Metadata dictionary
        output_file: Output file path
        compress: If True, use gzip compression (.json.gz)

    Returns:
        output_path: Path to saved file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Serialize all data
    serialized_data = {
        'metadata': serialize_data(metadata),
        'step_records': serialize_data(step_records),
        'serialization_info': {
            'timestamp': datetime.now().isoformat(),
            'total_steps': len(step_records),
            'serializer_version': '1.0',
            'compressed': compress
        }
    }

    # Save to file (with optional compression)
    if compress:
        import gzip
        output_path = output_path.with_suffix('.json.gz')
        with gzip.open(output_path, 'wt', encoding='utf-8') as f:
            json.dump(serialized_data, f, indent=2, ensure_ascii=False)
    else:
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(serialized_data, f, indent=2, ensure_ascii=False)
    
    return output_path


def load_step_records_from_json(input_file: Union[str, Path]) -> tuple:
    """
    Load step records and metadata from JSON file (with automatic decompression).

    Args:
        input_file: Input file path (.json or .json.gz)

    Returns:
        (step_records, metadata, serialization_info): Tuple of deserialized data
    """
    input_path = Path(input_file)

    # Load from file (automatic decompression)
    if input_path.suffix == '.gz':
        import gzip
        with gzip.open(input_path, 'rt', encoding='utf-8') as f:
            data = json.load(f)
    else:
        with open(input_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

    # Deserialize data
    step_records = deserialize_data(data['step_records'])
    metadata = deserialize_data(data['metadata'])
    serialization_info = data.get('serialization_info', {})
    
    return step_records, metadata, serialization_info


def generate_step_summary(step_records: List[Dict]) -> Dict:
    """
    Generate summary statistics from step records.

    Args:
        step_records: List of step record dictionaries

    Returns:
        summary: Dictionary with performance metrics and trends
    """
    if not step_records:
        return {
            'total_steps': 0,
            'error': 'No step records provided'
        }

    # Basic statistics
    total_steps = len(step_records)

    # Extract accessibility scores
    accessibility_scores = [
        step.get('true_accessibility', float('inf'))
        for step in step_records
    ]
    valid_scores = [s for s in accessibility_scores if s != float('inf')]

    # Count valid steps
    valid_steps = sum(1 for step in step_records if step.get('is_valid', True))

    # Calculate performance metrics
    if valid_scores:
        best_score = min(valid_scores)
        convergence_step = next(
            (i for i, score in enumerate(accessibility_scores) if score == best_score),
            len(step_records) - 1
        )

        # Calculate stability (standard deviation of recent scores)
        stability_window = max(1, total_steps // 10)
        recent_scores = valid_scores[-stability_window:] if len(valid_scores) >= stability_window else valid_scores
        stability = float(np.std(recent_scores)) if len(recent_scores) > 1 else 0.0
    else:
        best_score = float('inf')
        convergence_step = -1
        stability = float('inf')

    # Extract loss values
    loss_values = [step.get('total_loss', 0) for step in step_records]
    
    return {
        'total_steps': total_steps,
        'valid_steps': valid_steps,
        'validity_rate': valid_steps / total_steps if total_steps > 0 else 0.0,
        
        'performance': {
            'best_accessibility': float(best_score) if best_score != float('inf') else None,
            'convergence_step': convergence_step if convergence_step >= 0 else None,
            'final_accessibility': float(accessibility_scores[-1]) if accessibility_scores else None,
            'stability': float(stability)
        },
        
        'trends': {
            'accessibility_progression': [float(s) for s in accessibility_scores[-20:]],
            'loss_progression': [float(l) for l in loss_values[-20:]],
            'improvement_rate': calculate_improvement_rate(accessibility_scores)
        },
        
        'data_structure': {
            'fields_per_step': list(step_records[0].keys()) if step_records else [],
            'has_cds_logits': 'cds_logits' in step_records[0] if step_records else False,
            'cds_logits_shape': (
                list(step_records[0]['cds_logits'].shape) 
                if step_records and 'cds_logits' in step_records[0] 
                and hasattr(step_records[0]['cds_logits'], 'shape')
                else None
            )
        }
    }


def calculate_improvement_rate(scores: List[float]) -> float:
    """
    Calculate improvement rate from score progression.

    Args:
        scores: List of accessibility scores over time

    Returns:
        improvement_rate: Relative improvement from start to end (negative = worse)
    """
    if len(scores) < 10:
        return 0.0
    
    window_size = max(1, len(scores) // 10)
    initial_avg = np.mean(scores[:window_size])
    final_avg = np.mean(scores[-window_size:])
    
    if initial_avg == 0:
        return 0.0
    
    return float((final_avg - initial_avg) / initial_avg)
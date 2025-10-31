"""
DeepRaccess Wrapper for ID3 Framework

This module provides a wrapper for the DeepRaccess model to integrate with the ID3 framework.
It handles sequence-to-accessibility predictions with support for both discrete and continuous
(soft) probability distributions from ID3 model outputs.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Tuple, Optional, Union


import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'DeepRaccess'))
from DeepRaccess.mymodel import FCN


class DeepRaccessID3Wrapper(nn.Module):
    """
    DeepRaccess wrapper for ID3 framework.

    Wraps the DeepRaccess FCN model to process ID3's probability distributions and predict
    RNA accessibility. Supports both discrete (argmax) and continuous (soft embedding) paths,
    with sliding window processing for long sequences and device management.
    """

    def __init__(self, deepraccess_model_path: str = None, device: str = 'auto'):
        """
        Initialize DeepRaccess wrapper.

        Args:
            deepraccess_model_path: Path to the pretrained DeepRaccess model (.pth file)
            device: Device to use ('auto', 'cuda', 'mps', or 'cpu')
        """
        super().__init__()


        if device == 'auto':
            if torch.cuda.is_available():
                self.device = torch.device('cuda')
                print("🚀 CUDA detected, using GPU acceleration")
            elif torch.backends.mps.is_available():

                print("🍎 MPS (Apple Silicon) detected, but using CPU for stability")
                print("   (To force MPS usage, set device='mps')")
                self.device = torch.device('cpu')
            else:
                self.device = torch.device('cpu')
                print("💻 Using CPU for computation")
        else:

            if device == 'cuda' and not torch.cuda.is_available():
                print("⚠️ CUDA not available, falling back to CPU")
                device = 'cpu'
            elif device == 'mps' and not torch.backends.mps.is_available():
                print("⚠️ MPS not available, falling back to CPU")
                device = 'cpu'
            elif device == 'mps' and torch.backends.mps.is_available():
                print("🚀 Forcing MPS usage (Apple Silicon), will fall back to CPU if issues occur")
            self.device = torch.device(device)
        

        self.deepraccess_model = FCN()
        
        # Attempt to locate model file if path not provided
        if deepraccess_model_path is None:
            # Try common locations for the model file
            possible_paths = [
                # Primary location: project root
                os.path.join(os.path.dirname(__file__), '../../DeepRaccess/path/FCN_structured.pth'),
                'DeepRaccess/path/FCN_structured.pth',
                # Fallback location: scripts directory
                os.path.join(os.path.dirname(__file__), '../../scripts/DeepRaccess/path/FCN_structured.pth'),
                'scripts/DeepRaccess/path/FCN_structured.pth',
                # Alternative model: FCN_uniform
                os.path.join(os.path.dirname(__file__), '../../DeepRaccess/path/FCN_uniform.pth'),
                'DeepRaccess/path/FCN_uniform.pth',
                os.path.join(os.path.dirname(__file__), '../../scripts/DeepRaccess/path/FCN_uniform.pth'),
                'scripts/DeepRaccess/path/FCN_uniform.pth',
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    deepraccess_model_path = path
                    print(f"🔍 Auto-detected model file: {os.path.abspath(path)}")
                    break
        
        if deepraccess_model_path and os.path.exists(deepraccess_model_path):
            print(f"📦 Loading DeepRaccess model: {deepraccess_model_path}")
            checkpoint = torch.load(deepraccess_model_path, map_location=self.device, weights_only=False)

            # Extract state dict from checkpoint
            if 'model_state_dict' in checkpoint:
                state_dict = checkpoint['model_state_dict']
            else:
                state_dict = checkpoint

            # Remove 'module.' prefix if present (from DataParallel)
            new_state_dict = {}
            for key, value in state_dict.items():
                if key.startswith('module.'):
                    new_key = key[7:]
                else:
                    new_key = key
                new_state_dict[new_key] = value

            self.deepraccess_model.load_state_dict(new_state_dict)
            print("✅ DeepRaccess model loaded successfully")
        else:
            print("⚠️  Pretrained model not found, using random initialization")
            print("   Note: For better results, run: bash scripts/setup_deepraccess.sh")
            
        self.deepraccess_model.eval().to(self.device)
        


        

        self.embedding_matrix = self.deepraccess_model.embedding.weight  # [6, 120]
        

        self.max_length = 440
        self.window_size = 440
        self.step_size = 330
        self.overlap = self.window_size - self.step_size  # 110nt
        self.trim_edges = 55
        
        print(f"✅ DeepRaccess wrapper initialized, device: {self.device}")
        
    def id3_to_deepraccess_probs(self, id3_probs: torch.Tensor) -> torch.Tensor:
        """
        Convert ID3 probability distribution to DeepRaccess format.

        ID3 uses [A, C, G, U] encoding (4-dimensional)
        DeepRaccess uses [padding, mask, A, U, G, C] encoding (6-dimensional)

        Args:
            id3_probs: ID3 probability tensor [batch, length, 4]

        Returns:
            deepraccess_probs: DeepRaccess probability tensor [batch, length, 6]
        """
        batch_size, length, _ = id3_probs.shape
        deepraccess_probs = torch.zeros(batch_size, length, 6,
                                      device=id3_probs.device, dtype=id3_probs.dtype)

        # Map ID3 indices to DeepRaccess indices
        deepraccess_probs[:, :, 0] = 0.0                # padding (not used for valid sequences)
        deepraccess_probs[:, :, 1] = 0.0                # mask (not used for valid sequences)
        deepraccess_probs[:, :, 2] = id3_probs[:, :, 0] # A: ID3[0] → DeepRaccess[2]
        deepraccess_probs[:, :, 3] = id3_probs[:, :, 3] # U: ID3[3] → DeepRaccess[3]
        deepraccess_probs[:, :, 4] = id3_probs[:, :, 2] # G: ID3[2] → DeepRaccess[4]
        deepraccess_probs[:, :, 5] = id3_probs[:, :, 1] # C: ID3[1] → DeepRaccess[5]
        
        return deepraccess_probs
    
    def deepraccess_to_id3_probs(self, deepraccess_probs: torch.Tensor) -> torch.Tensor:
        """
        Convert DeepRaccess probability distribution back to ID3 format.

        Args:
            deepraccess_probs: DeepRaccess probability tensor [batch, length, 6]

        Returns:
            id3_probs: ID3 probability tensor [batch, length, 4]
        """
        batch_size, length, _ = deepraccess_probs.shape
        id3_probs = torch.zeros(batch_size, length, 4,
                               device=deepraccess_probs.device, dtype=deepraccess_probs.dtype)

        # Map DeepRaccess indices back to ID3 indices
        id3_probs[:, :, 0] = deepraccess_probs[:, :, 2]  # A: DeepRaccess[2] → ID3[0]
        id3_probs[:, :, 1] = deepraccess_probs[:, :, 5]  # C: DeepRaccess[5] → ID3[1]
        id3_probs[:, :, 2] = deepraccess_probs[:, :, 4]  # G: DeepRaccess[4] → ID3[2]
        id3_probs[:, :, 3] = deepraccess_probs[:, :, 3]  # U: DeepRaccess[3] → ID3[3]
        
        return id3_probs
    
    def compute_soft_embedding(self, probs: torch.Tensor) -> torch.Tensor:
        """
        Compute soft embeddings from probability distributions.

        Args:
            probs: Probability tensor [batch, length, 6] (DeepRaccess format)

        Returns:
            soft_embeddings: Soft embedding tensor [batch, length, 120]
        """
        # Matrix multiplication: [batch, length, 6] @ [6, 120] → [batch, length, 120]
        soft_embeddings = torch.matmul(probs, self.embedding_matrix)
        return soft_embeddings
    
    def sequence_to_indices(self, sequence: str) -> torch.Tensor:
        """
        Convert nucleotide sequence string to DeepRaccess index tensor.

        Args:
            sequence: Nucleotide sequence string (e.g., "AUGC")

        Returns:
            indices: Index tensor [length] with DeepRaccess encoding
        """
        # DeepRaccess encoding: {pad:0, mask:1, A:2, U/T:3, G:4, C:5}
        nucleotide_to_index = {
            'A': 2, 'a': 2,
            'T': 3, 't': 3, 'U': 3, 'u': 3,
            'G': 4, 'g': 4,
            'C': 5, 'c': 5,
            'N': 1, 'n': 1,  # mask for unknown nucleotides
        }
        
        indices = []
        for nt in sequence:
            indices.append(nucleotide_to_index.get(nt, 1))
        
        return torch.tensor(indices, dtype=torch.long)
    
    def process_sliding_window(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Process sequences using sliding window for long sequences.

        Args:
            embeddings: Either indices [batch, length] or soft embeddings [batch, length, 120]

        Returns:
            accessibility: Accessibility predictions [batch, length]
        """
        batch_size, seq_len = embeddings.shape[0], embeddings.shape[1]

        # For short sequences, process directly without sliding window
        if seq_len <= self.max_length:
            if embeddings.dim() == 3:
                return self._forward_with_soft_embeddings(embeddings)
            else:
                padded_embeddings = F.pad(embeddings, (0, self.max_length - seq_len), value=0)
                output = self.deepraccess_model(padded_embeddings)
                return output[:, :seq_len]
        
        # For long sequences, use sliding window with overlap
        num_windows = 1 + (seq_len - self.overlap) // self.step_size
        accessibility_parts = []
        
        for i in range(num_windows):
            start = i * self.step_size
            end = min(start + self.window_size, seq_len)

            # Extract window and pad to max_length
            if embeddings.dim() == 3:
                window_embeddings = embeddings[:, start:end, :]
                # Pad if necessary
                pad_length = self.max_length - (end - start)
                if pad_length > 0:
                    window_embeddings = F.pad(window_embeddings, (0, 0, 0, pad_length), value=0)
                window_output = self._forward_with_soft_embeddings(window_embeddings)
            else:
                window_indices = embeddings[:, start:end]
                # Pad if necessary
                pad_length = self.max_length - (end - start)
                if pad_length > 0:
                    window_indices = F.pad(window_indices, (0, pad_length), value=0)
                window_output = self.deepraccess_model(window_indices)

            # Trim edges to avoid boundary effects, except at sequence ends
            window_length = end - start
            if i == 0:
                valid_output = window_output[:, :window_length - self.trim_edges]
            elif i == num_windows - 1:
                valid_output = window_output[:, self.trim_edges:window_length]
            else:
                valid_output = window_output[:, self.trim_edges:window_length - self.trim_edges]
            
            accessibility_parts.append(valid_output)
        
        # Concatenate all window results
        return torch.cat(accessibility_parts, dim=1)
    
    def _forward_with_soft_embeddings(self, soft_embeddings: torch.Tensor) -> torch.Tensor:
        """
        Forward pass through DeepRaccess model using soft embeddings.

        Args:
            soft_embeddings: Soft embedding tensor [batch, length, 120]

        Returns:
            accessibility: Accessibility predictions [batch, length]
        """
        # Transpose for conv layers: [batch, length, 120] → [batch, 120, length]
        x = soft_embeddings.transpose(1, 2)  # [batch, 120, length]

        # Pass through convolutional layers
        for conv_layer in self.deepraccess_model.convs:
            x = conv_layer(x)

        # DeepRaccess FCN outputs [batch, 1, length], squeeze to [batch, length]
        return x.squeeze(1)  # [batch, length]
    
    def forward(self,
                id3_probs: torch.Tensor,
                discrete: bool = False,
                return_window_scores: bool = False,
                atg_position: Optional[int] = None,
                window_size: int = 35) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """
        Forward pass to compute RNA accessibility from ID3 probabilities.

        Args:
            id3_probs: ID3 probability tensor [batch, length, 4]
            discrete: If True, use discrete (argmax) path; if False, use continuous (soft) path
            return_window_scores: If True, return accessibility at ATG window position
            atg_position: Position of ATG codon for window-based scoring
            window_size: Size of window around ATG (default 35, corresponding to -19 to +15)

        Returns:
            accessibility: Full accessibility predictions [batch, length]
            window_scores: (optional) Accessibility at ATG window position [batch, 1]
        """
        if discrete:
            # Discrete path: argmax → indices → DeepRaccess
            indices = torch.argmax(id3_probs, dim=-1)  # [batch, length]

            # Convert ID3 indices to DeepRaccess indices
            deepraccess_indices = self._convert_id3_indices_to_deepraccess(indices)

            # Process through DeepRaccess model
            accessibility = self.process_sliding_window(deepraccess_indices)

        else:
            # Continuous path: probs → soft embeddings → DeepRaccess
            # This path maintains gradient flow for optimization
            deepraccess_probs = self.id3_to_deepraccess_probs(id3_probs)

            # Compute soft embeddings weighted by probabilities
            soft_embeddings = self.compute_soft_embedding(deepraccess_probs)

            # Process through DeepRaccess model
            accessibility = self.process_sliding_window(soft_embeddings)
        
        # Optionally return accessibility at specific ATG window position
        if return_window_scores and atg_position is not None:
            # Window position is ATG-19 (start of -19 to +15 window)
            target_position = max(0, atg_position - 19)
            if target_position >= accessibility.shape[1]:
                target_position = accessibility.shape[1] - 1
            window_scores = accessibility[:, target_position:target_position+1]
            return accessibility, window_scores
        
        return accessibility
    
    def _convert_id3_indices_to_deepraccess(self, id3_indices: torch.Tensor) -> torch.Tensor:
        """
        Convert ID3 discrete indices to DeepRaccess discrete indices.

        Args:
            id3_indices: ID3 index tensor [batch, length] with values in {0, 1, 2, 3}

        Returns:
            deepraccess_indices: DeepRaccess index tensor [batch, length] with values in {2, 3, 4, 5}
        """
        # ID3: {A:0, C:1, G:2, U:3}
        # DeepRaccess: {pad:0, mask:1, A:2, U:3, G:4, C:5}
        mapping = torch.tensor([2, 5, 4, 3], device=id3_indices.device)  # [A, C, G, U] -> [2, 5, 4, 3]
        deepraccess_indices = mapping[id3_indices]
        return deepraccess_indices
    
    def compute_atg_window_accessibility(self,
                                       id3_probs: torch.Tensor,
                                       atg_position: int,
                                       window_size: int = 35,
                                       discrete: bool = False) -> torch.Tensor:
        """
        Compute accessibility at ATG translation initiation site window.

        DeepRaccess predicts accessibility with a specific window alignment:
        - Each output position corresponds to a 35nt window (-19 to +15 relative to that position)
        - For ATG at position i, we use the prediction at position (i-19) to get the window centered on ATG

        Args:
            id3_probs: ID3 probability tensor [batch, length, 4]
            atg_position: Position of ATG start codon in the sequence
            window_size: Window size (default 35nt, corresponding to -19 to +15)
            discrete: Whether to use discrete or continuous path

        Returns:
            window_accessibility: Accessibility score at ATG window [batch]
        """
        # Get full accessibility predictions
        accessibility = self.forward(id3_probs, discrete=discrete)

        # Calculate target position: ATG position - 19 (window alignment)
        target_position = atg_position - 19

        # Handle edge cases
        if target_position < 0:
            # ATG too close to sequence start
            target_position = 0
            print(f"⚠️ ATG position {atg_position} too close to sequence start, using position 0 output")
        elif target_position >= accessibility.shape[1]:
            # Target position exceeds sequence length
            target_position = accessibility.shape[1] - 1
            print(f"⚠️ Computed position exceeds sequence length, using position {target_position} output")

        # Extract accessibility at window position
        window_accessibility = accessibility[:, target_position]
        
        return window_accessibility


def load_deepraccess_model(model_path: str = None, device: str = 'auto') -> DeepRaccessID3Wrapper:
    """
    Convenience function to load DeepRaccess model wrapper.

    Args:
        model_path: Path to pretrained model file. If None, searches default locations.
        device: Device to use ('auto', 'cuda', 'mps', or 'cpu')

    Returns:
        wrapper: Initialized DeepRaccessID3Wrapper instance
    """
    # Try to locate model automatically if not provided
    if model_path is None:
        # Get the base directory (ID3-github root)
        # This function is in id3/utils/deepraccess_wrapper.py
        # So we go up 3 levels to reach the root
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        default_paths = [
            os.path.join(base_dir, 'DeepRaccess/path/FCN_structured.pth'),
            os.path.join(base_dir, 'DeepRaccess/path/FCN_uniform.pth'),
        ]
        
        for path in default_paths:
            if os.path.exists(path):
                model_path = path
                break
    
    wrapper = DeepRaccessID3Wrapper(model_path, device)
    return wrapper


if __name__ == "__main__":
    # Test script for DeepRaccess wrapper
    print("🧪 Testing DeepRaccess ID3 wrapper...")

    # Load wrapper
    wrapper = load_deepraccess_model()

    # Create random test input
    batch_size, seq_len = 2, 100
    id3_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)

    print(f"Input shape: {id3_probs.shape}")

    # Test continuous (soft) path
    accessibility_soft = wrapper(id3_probs, discrete=False)
    print(f"Continuous path output shape: {accessibility_soft.shape}")
    print(f"Continuous path accessibility range: {accessibility_soft.min():.4f} - {accessibility_soft.max():.4f}")

    # Test discrete (hard) path
    accessibility_hard = wrapper(id3_probs, discrete=True)
    print(f"Discrete path output shape: {accessibility_hard.shape}")
    print(f"Discrete path accessibility range: {accessibility_hard.min():.4f} - {accessibility_hard.max():.4f}")

    # Test ATG window accessibility
    atg_position = 50
    window_accessibility = wrapper.compute_atg_window_accessibility(
        id3_probs, atg_position, window_size=35, discrete=False)
    print(f"ATG window accessibility (-19 to +15): {window_accessibility}")

    print("✅ DeepRaccess wrapper testing complete!")
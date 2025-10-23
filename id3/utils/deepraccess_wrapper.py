"""







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

    





    """
    
    def __init__(self, deepraccess_model_path: str = None, device: str = 'auto'):
        """


        Args:


        """
        super().__init__()


        if device == 'auto':
            if torch.cuda.is_available():
                self.device = torch.device('cuda')
                print("🚀 检测到CUDA，使用GPU加速")
            elif torch.backends.mps.is_available():

                print("🍎 检测到MPS (Apple Silicon)，但使用CPU确保稳定性")
                print("   (如需强制使用MPS，请设置 device='mps')")
                self.device = torch.device('cpu')
            else:
                self.device = torch.device('cpu')
                print("💻 使用CPU计算")
        else:

            if device == 'cuda' and not torch.cuda.is_available():
                print("⚠️ CUDA不可用，回退到CPU")
                device = 'cpu'
            elif device == 'mps' and not torch.backends.mps.is_available():
                print("⚠️ MPS不可用，回退到CPU")
                device = 'cpu'
            elif device == 'mps' and torch.backends.mps.is_available():
                print("🚀 强制使用MPS (Apple Silicon)，如遇问题将回退到CPU")
            self.device = torch.device(device)
        

        self.deepraccess_model = FCN()
        

        if deepraccess_model_path is None:

            possible_paths = [
                os.path.join(os.path.dirname(__file__), '../../DeepRaccess/path/FCN_structured.pth'),
                'DeepRaccess/path/FCN_structured.pth',
                os.path.join(os.path.dirname(__file__), '../../DeepRaccess/path/FCN_uniform.pth'),
                'DeepRaccess/path/FCN_uniform.pth',
            ]
            for path in possible_paths:
                if os.path.exists(path):
                    deepraccess_model_path = path
                    print(f"🔍 自动找到模型文件: {path}")
                    break
        
        if deepraccess_model_path and os.path.exists(deepraccess_model_path):
            print(f"📦 加载DeepRaccess模型: {deepraccess_model_path}")
            checkpoint = torch.load(deepraccess_model_path, map_location=self.device, weights_only=False)
            

            if 'model_state_dict' in checkpoint:
                state_dict = checkpoint['model_state_dict']
            else:
                state_dict = checkpoint
            

            new_state_dict = {}
            for key, value in state_dict.items():
                if key.startswith('module.'):
                    new_key = key[7:]
                else:
                    new_key = key
                new_state_dict[new_key] = value
            
            self.deepraccess_model.load_state_dict(new_state_dict)
            print("✅ DeepRaccess模型加载成功")
        else:
            print("⚠️  未找到预训练模型，使用随机初始化")
            
        self.deepraccess_model.eval().to(self.device)
        


        

        self.embedding_matrix = self.deepraccess_model.embedding.weight  # [6, 120]
        

        self.max_length = 440
        self.window_size = 440
        self.step_size = 330
        self.overlap = self.window_size - self.step_size  # 110nt
        self.trim_edges = 55
        
        print(f"✅ DeepRaccess包装器已初始化，设备: {self.device}")
        
    def id3_to_deepraccess_probs(self, id3_probs: torch.Tensor) -> torch.Tensor:
        """

        


        
        Args:

            
        Returns:

        """
        batch_size, length, _ = id3_probs.shape
        deepraccess_probs = torch.zeros(batch_size, length, 6, 
                                      device=id3_probs.device, dtype=id3_probs.dtype)
        

        deepraccess_probs[:, :, 0] = 0.0                # padding
        deepraccess_probs[:, :, 1] = 0.0                # mask  
        deepraccess_probs[:, :, 2] = id3_probs[:, :, 0] # A: 0→2
        deepraccess_probs[:, :, 3] = id3_probs[:, :, 3] # U: 3→3
        deepraccess_probs[:, :, 4] = id3_probs[:, :, 2] # G: 2→4
        deepraccess_probs[:, :, 5] = id3_probs[:, :, 1] # C: 1→5
        
        return deepraccess_probs
    
    def deepraccess_to_id3_probs(self, deepraccess_probs: torch.Tensor) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

        """
        batch_size, length, _ = deepraccess_probs.shape
        id3_probs = torch.zeros(batch_size, length, 4,
                               device=deepraccess_probs.device, dtype=deepraccess_probs.dtype)
        

        id3_probs[:, :, 0] = deepraccess_probs[:, :, 2]  # A: 2→0
        id3_probs[:, :, 1] = deepraccess_probs[:, :, 5]  # C: 5→1
        id3_probs[:, :, 2] = deepraccess_probs[:, :, 4]  # G: 4→2
        id3_probs[:, :, 3] = deepraccess_probs[:, :, 3]  # U: 3→3
        
        return id3_probs
    
    def compute_soft_embedding(self, probs: torch.Tensor) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

        """

        soft_embeddings = torch.matmul(probs, self.embedding_matrix)
        return soft_embeddings
    
    def sequence_to_indices(self, sequence: str) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

        """

        nucleotide_to_index = {
            'A': 2, 'a': 2,
            'T': 3, 't': 3, 'U': 3, 'u': 3,
            'G': 4, 'g': 4,
            'C': 5, 'c': 5,
            'N': 1, 'n': 1,  # mask for unknown
        }
        
        indices = []
        for nt in sequence:
            indices.append(nucleotide_to_index.get(nt, 1))
        
        return torch.tensor(indices, dtype=torch.long)
    
    def process_sliding_window(self, embeddings: torch.Tensor) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

        """
        batch_size, seq_len = embeddings.shape[0], embeddings.shape[1]
        

        if seq_len <= self.max_length:
            if embeddings.dim() == 3:
                return self._forward_with_soft_embeddings(embeddings)
            else:
                padded_embeddings = F.pad(embeddings, (0, self.max_length - seq_len), value=0)
                output = self.deepraccess_model(padded_embeddings)
                return output[:, :seq_len]
        

        num_windows = 1 + (seq_len - self.overlap) // self.step_size
        accessibility_parts = []
        
        for i in range(num_windows):
            start = i * self.step_size
            end = min(start + self.window_size, seq_len)
            

            if embeddings.dim() == 3:
                window_embeddings = embeddings[:, start:end, :]

                pad_length = self.max_length - (end - start)
                if pad_length > 0:
                    window_embeddings = F.pad(window_embeddings, (0, 0, 0, pad_length), value=0)
                window_output = self._forward_with_soft_embeddings(window_embeddings)
            else:
                window_indices = embeddings[:, start:end]

                pad_length = self.max_length - (end - start)
                if pad_length > 0:
                    window_indices = F.pad(window_indices, (0, pad_length), value=0)
                window_output = self.deepraccess_model(window_indices)
            

            window_length = end - start
            if i == 0:
                valid_output = window_output[:, :window_length - self.trim_edges]
            elif i == num_windows - 1:
                valid_output = window_output[:, self.trim_edges:window_length]
            else:
                valid_output = window_output[:, self.trim_edges:window_length - self.trim_edges]
            
            accessibility_parts.append(valid_output)
        

        return torch.cat(accessibility_parts, dim=1)
    
    def _forward_with_soft_embeddings(self, soft_embeddings: torch.Tensor) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

        """

        x = soft_embeddings.transpose(1, 2)  # [batch, 120, length]
        

        for conv_layer in self.deepraccess_model.convs:
            x = conv_layer(x)
        


        return x.squeeze(1)  # [batch, length]
    
    def forward(self, 
                id3_probs: torch.Tensor, 
                discrete: bool = False,
                return_window_scores: bool = False,
                atg_position: Optional[int] = None,
                window_size: int = 35) -> Union[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        """

        
        Args:





            
        Returns:


        """
        if discrete:

            indices = torch.argmax(id3_probs, dim=-1)  # [batch, length]
            

            deepraccess_indices = self._convert_id3_indices_to_deepraccess(indices)
            

            accessibility = self.process_sliding_window(deepraccess_indices)
            
        else:


            deepraccess_probs = self.id3_to_deepraccess_probs(id3_probs)
            

            soft_embeddings = self.compute_soft_embedding(deepraccess_probs)
            

            accessibility = self.process_sliding_window(soft_embeddings)
        

        if return_window_scores and atg_position is not None:

            target_position = max(0, atg_position - 19)
            if target_position >= accessibility.shape[1]:
                target_position = accessibility.shape[1] - 1
            window_scores = accessibility[:, target_position:target_position+1]
            return accessibility, window_scores
        
        return accessibility
    
    def _convert_id3_indices_to_deepraccess(self, id3_indices: torch.Tensor) -> torch.Tensor:
        """

        
        Args:

            
        Returns:

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


        


        
        Args:




            
        Returns:

        """

        accessibility = self.forward(id3_probs, discrete=discrete)
        

        target_position = atg_position - 19
        

        if target_position < 0:

            target_position = 0
            print(f"⚠️ ATG位置{atg_position}太靠近序列开始，使用位置0的输出")
        elif target_position >= accessibility.shape[1]:

            target_position = accessibility.shape[1] - 1
            print(f"⚠️ 计算位置超出序列长度，使用位置{target_position}的输出")
        

        window_accessibility = accessibility[:, target_position]
        
        return window_accessibility


def load_deepraccess_model(model_path: str = None, device: str = 'auto') -> DeepRaccessID3Wrapper:
    """

    
    Args:


        
    Returns:

    """

    if model_path is None:



        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        default_paths = [
            os.path.join(base_dir, 'DeepRaccess/path/FCN_structured.pth'),
            os.path.join(base_dir, 'DeepRaccess/path/FCN_uniform.pth'),

            '/work/gg53/d58004/Degradation/DeepRaccess/path/FCN_structured.pth',
            '/work/gg53/d58004/Degradation/DeepRaccess/path/FCN_uniform.pth',
        ]
        
        for path in default_paths:
            if os.path.exists(path):
                model_path = path
                break
    
    wrapper = DeepRaccessID3Wrapper(model_path, device)
    return wrapper


if __name__ == "__main__":

    print("🧪 测试DeepRaccess ID3包装器...")
    

    wrapper = load_deepraccess_model()
    

    batch_size, seq_len = 2, 100
    id3_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)
    
    print(f"输入形状: {id3_probs.shape}")
    

    accessibility_soft = wrapper(id3_probs, discrete=False)
    print(f"连续路径输出形状: {accessibility_soft.shape}")
    print(f"连续路径可及性范围: {accessibility_soft.min():.4f} - {accessibility_soft.max():.4f}")
    

    accessibility_hard = wrapper(id3_probs, discrete=True)
    print(f"离散路径输出形状: {accessibility_hard.shape}")
    print(f"离散路径可及性范围: {accessibility_hard.min():.4f} - {accessibility_hard.max():.4f}")
    

    atg_position = 50
    window_accessibility = wrapper.compute_atg_window_accessibility(
        id3_probs, atg_position, window_size=35, discrete=False)
    print(f"ATG窗口可及性（-19到+15）: {window_accessibility}")
    
    print("✅ DeepRaccess包装器测试完成!")
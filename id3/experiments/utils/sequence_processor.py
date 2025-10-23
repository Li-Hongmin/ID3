#!/usr/bin/env python3
"""
序列处理工具模块

负责序列的转换、编码和处理
"""

import torch
from typing import List, Tuple, Optional


class SequenceProcessor:
    """序列处理器"""

    NUCLEOTIDES = ['A', 'C', 'G', 'U']

    def __init__(self, device: Optional[torch.device] = None):
        """
        初始化处理器

        Args:
            device: 计算设备
        """
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device

    def sequence_to_probs(self, sequence: str) -> torch.Tensor:
        """
        将核苷酸序列转换为概率张量

        Args:
            sequence: 核苷酸序列字符串

        Returns:
            概率张量 [L, 4]
        """
        probs = []
        for nt in sequence.upper().replace('T', 'U'):
            prob = torch.zeros(4, dtype=torch.float32, device=self.device)
            if nt in self.NUCLEOTIDES:
                prob[self.NUCLEOTIDES.index(nt)] = 1.0
            else:
                # 未知核苷酸使用均匀分布
                prob[:] = 0.25
            probs.append(prob)
        return torch.stack(probs)

    def probs_to_sequence(self, probs: torch.Tensor) -> str:
        """
        将概率张量转换为核苷酸序列

        Args:
            probs: 概率张量 [L, 4]

        Returns:
            核苷酸序列字符串
        """
        # 获取最大概率的核苷酸
        indices = probs.argmax(dim=-1)
        sequence = ''.join([self.NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
        return sequence

    def prepare_utr_probs(
        self,
        utr5: str,
        utr3: str
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        准备UTR概率张量

        Args:
            utr5: 5'UTR序列
            utr3: 3'UTR序列

        Returns:
            (utr5_probs, utr3_probs)
        """
        utr5_probs = self.sequence_to_probs(utr5)
        utr3_probs = self.sequence_to_probs(utr3)
        return utr5_probs, utr3_probs

    def concat_mrna_probs(
        self,
        utr5_probs: torch.Tensor,
        cds_probs: torch.Tensor,
        utr3_probs: torch.Tensor
    ) -> torch.Tensor:
        """
        拼接完整mRNA概率

        Args:
            utr5_probs: 5'UTR概率 [L5, 4]
            cds_probs: CDS概率 [L_cds, 4]
            utr3_probs: 3'UTR概率 [L3, 4]

        Returns:
            完整mRNA概率 [1, L_total, 4]
        """
        full_probs = torch.cat([utr5_probs, cds_probs, utr3_probs], dim=0)
        return full_probs.unsqueeze(0)  # 添加批次维度

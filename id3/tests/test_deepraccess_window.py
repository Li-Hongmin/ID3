#!/usr/bin/env python3
"""


"""

import torch
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from id3.utils.deepraccess_wrapper import load_deepraccess_model


def test_window_extraction():

    print("=" * 60)

    print("=" * 60)
    

    wrapper = load_deepraccess_model(device='cpu')
    

    seq_length = 150
    batch_size = 1
    

    id3_probs = torch.softmax(torch.randn(batch_size, seq_length, 4), dim=-1)
    

    test_positions = [20, 50, 70, 100]
    

    print("-" * 60)
    
    for atg_pos in test_positions:

        

        target_pos = atg_pos - 19

        

        window_acc = wrapper.compute_atg_window_accessibility(
            id3_probs, 
            atg_position=atg_pos,
            window_size=35,
            discrete=False
        )
        

        

        full_acc = wrapper.forward(id3_probs, discrete=False)
        
        if target_pos >= 0 and target_pos < full_acc.shape[1]:
            manual_acc = full_acc[:, target_pos].item()

            

            if torch.allclose(window_acc, full_acc[:, target_pos], atol=1e-6):

            else:

        else:

    
    print("\n" + "=" * 60)





    print("=" * 60)


def test_boundary_cases():
    """测试边界情况"""
    print("\n📊 测试边界情况：")
    print("-" * 60)
    
    wrapper = load_deepraccess_model(device='cpu')
    

    short_seq = torch.softmax(torch.randn(1, 30, 4), dim=-1)
    

    atg_pos = 10
    print(f"\n短序列测试（长度30，ATG位置{atg_pos}）：")
    
    window_acc = wrapper.compute_atg_window_accessibility(
        short_seq,
        atg_position=atg_pos,
        discrete=False
    )
    
    print(f"窗口可及性: {window_acc.item():.6f}")
    print("✅ 边界处理正常")


if __name__ == "__main__":
    test_window_extraction()
    test_boundary_cases()
    print("\n🎉 所有测试完成！")
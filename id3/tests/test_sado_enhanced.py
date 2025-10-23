#!/usr/bin/env python3
"""


"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import time
from typing import Dict, List
import matplotlib.pyplot as plt

from id3.optimizers.cai.sado import SADOOptimizer


def generate_test_distribution(seq_len: int, num_codons: int = 6, skew: float = 0.2) -> torch.Tensor:

    probs = torch.rand(seq_len, num_codons)

    for i in range(seq_len):
        preferred = np.random.randint(0, num_codons)
        probs[i, preferred] += skew

    probs = probs / probs.sum(dim=1, keepdim=True)
    return probs


def test_enhanced_sado(protein_name: str = 'O15263'):
    """测试增强版SADO的各种配置"""
    
    print(f"\n{'='*80}")
    print(f"测试增强版SADO优化器 - 蛋白质: {protein_name}")
    print(f"{'='*80}\n")
    

    try:
        import json
        data_path = f'data/sequences/{protein_name}/protein_info.json'
        with open(data_path, 'r') as f:
            protein_data = json.load(f)
        amino_sequence = protein_data['amino_acid_sequence']
        print(f"✅ 成功加载蛋白质 {protein_name}")
        print(f"   序列长度: {len(amino_sequence)} AA")
    except:

        amino_sequence = "MKPEITVQYDKENGNVIVTQRNEKHHQLLCEVLKSLPDGTDYKVYVQIQSEKLGYQPPVFDQMDEGKYDDVNGQPEGFSKDLIKKMAEHQDKEKSGEKVLASVKVKEAGAGAGGKKPSMGSKKGKTEAKEKKKKVVKVGDQVQISAVYSITPPDKYELQFQIGDALMEAIALGDLPHNTIYLYTPEGNVFRGFYSAIGECFGSMQRVFDKGTIDQIHLEIEELMQQRLVELARKIERLREKSDHSDKEFKNCKPKKKKRPTFKQFVKNKKAEIDGDIRNKDGKQIVLNDAVKVSPKKEAKEKAEKNRAKKMDDRNAEDFVDISVGASRKSKTKLPTSKKSRPKPKLTKKEEEEEKDKQTQPIFEKIVFDIDSEEKEIKSQLKQFTITKEEKKRKKKSQKRTTADFTKSAKKRKKAKVEPQFEKEAKAKSVKDYFKTDSSKQAKKTLNTKAGKENRGIGKPVNKTYYDKKFKAQKSIKKDSAEGTMKIRFKRSKPKTSVRQLAVSKLKEGTIKVGDRVTVDSNAKDKMPDSKRMLVQDIDVGVRQKYSISKEKFGKIVKMSRIAKKQKPKVDIRFKWWGRTNSQFRRDFVNHFSGKLNEYLKNKMMNLISEIIRWIREDHQDEEKKSEYQFGLPPKKVFRPKRNRSWYTRVMCLFNLKYSEEDWEKDNKRIEAWKEKKEEEKRKLLDEEVAEKLEKVLEKDQRFTDQDVLDAYIEYDKDRYDKKKFAKTNIDLLIEVFGFEVIKREKENAFTEQQKEQINLEKQIFEENRNFVNDVKETVEKNLDLEKTLFGTTKKINLEKTVVKKQVSEKQTPIIKEIVVKVGDVVNISKVISGSEKVDVVSKDSKTMKRDRKKTPNSTVDTTKFLKPGDEISDLSNKGKRIGTKLVNSEPGTKPVQLDPLAFMQTAVGASPMQFVAVSGKPVEIKEQTVGVKPNSKGRKLVDTDSVTPFMQTVGVSSQFYEVSGKQKKITDNTKNSGKKSVGTAMNNPPNTRLLALSGKQKVITSQEKFSGPTKVVDLSTSGDTQNVKKKSAAIKDSLSPVATVDDAKKRNKDDVVKDRKNKKEKMTSIRSKSTQQEEDDEGDHTQLAKKLFQYINRQAKKRRRHHKQKQQQQAQHQKAQKQKQKRHHKEHRKKEKRKKEKQHQREKQREQEKKQKQKQERKQEEKKQRKRQQEKERLKKQERERKRRQEKRHKQQQRKQQELKKKEQEKEKKEQRKRLQGQKQQERQRKAREEEMQREAQERQREREKEREQRQREREMEQEKRRRRERQRQQDKEQKREQEREQRKREREQREREREKEQRKAEQERKEREQREREREREREKRKEEERKERKERKDRERERERDREREREREQEREKREREREQKKEKEEKDRKKERLRQEKENQKQQEKQRQHEREKQRQEEKEKQQRGKQLLKLEKKRERLVRGQNKTDGSGGKHTNTKTQEQRKRERQEKTTVQKESGEASMKRKKKKKQSGFQAHPEKRRLTITKEEIKTKLMFPEGEVIGLGGIIKQTRQQILNHSRVVAIFGGVGGEILNYATANEKPKRALLQVKGCFGGEAIAFYGFGQKESLCGFIGEIVHVDPATGRPYDYPKLTLRQVDPNFQFHVKPAFKQVFAKRNPGLIPWLDKNPKDLIAIFCGGIGGNLLYIDKGSFKVHFPNKDIKGSYLKEFIIRNCKFRGEKIAFVGFGIGLTHQAYINNMVKTGFDSEQLRLSLNWCRRFSRIHNALTVIPWLQKIPAKLSAKMTTMNSWTMDVLRQEAKEEPGEQNEKEKLKLKK"
        print(f"   使用默认O15263序列，长度: {len(amino_sequence)} AA")
    

    seq_len = len(amino_sequence)
    pi_accessibility = generate_test_distribution(seq_len, num_codons=6, skew=0.3)
    

    test_configs = [
        {
            'name': '原始SADO (Gamma初始化)',
            'use_binary_search': False,
            'use_difference_driven': False,
            'gamma': 0.3
        },
        {
            'name': '平衡状态初始化 + 差异驱动',
            'use_binary_search': False,
            'use_difference_driven': True,
            'gamma': 0.3
        },
        {
            'name': '二分查找初始化 + 原始优化',
            'use_binary_search': True,
            'use_difference_driven': False,
            'gamma': 0.3
        },
        {
            'name': '二分查找 + 差异驱动 (完整增强)',
            'use_binary_search': True,
            'use_difference_driven': True,
            'gamma': 0.3
        }
    ]
    
    results = []
    
    for config in test_configs:
        print(f"\n{'='*60}")
        print(f"测试配置: {config['name']}")
        print(f"{'='*60}")
        

        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=torch.device('cuda' if torch.cuda.is_available() else 'cpu'),
            amino_acid_sequence=amino_sequence
        )
        

        optimizer.reset()
        

        cai_values = []
        times = []
        
        for iteration in range(10):
            start_time = time.time()
            

            distribution, metadata = optimizer.optimize(
                pi_accessibility=pi_accessibility.clone(),
                target_cai=0.8,
                amino_acid_sequence=amino_sequence,
                use_binary_search=config['use_binary_search'],
                use_difference_driven=config['use_difference_driven'],
                gamma=config['gamma']
            )
            
            elapsed = time.time() - start_time
            
            cai_values.append(metadata['final_cai'])
            times.append(elapsed)
            

            if iteration == 0 or (iteration + 1) % 5 == 0:
                print(f"  迭代 {iteration+1:2d}: CAI={metadata['final_cai']:.4f}, "
                      f"时间={elapsed*1000:.1f}ms, "
                      f"初始化方法={metadata.get('init_method', 'unknown')}")
        

        avg_cai = np.mean(cai_values)
        std_cai = np.std(cai_values)
        avg_time = np.mean(times) * 1000
        std_time = np.std(times) * 1000
        
        result = {
            'config': config['name'],
            'avg_cai': avg_cai,
            'std_cai': std_cai,
            'avg_time_ms': avg_time,
            'std_time_ms': std_time,
            'cai_values': cai_values,
            'times': times,
            'diversity': optimizer.get_diversity_stats()
        }
        results.append(result)
        
        print(f"\n📊 统计结果:")
        print(f"  平均CAI: {avg_cai:.4f} ± {std_cai:.4f}")
        print(f"  平均时间: {avg_time:.1f} ± {std_time:.1f} ms")
        print(f"  序列多样性: {result['diversity']['unique_ratio']*100:.1f}%")
    
    return results


def compare_long_sequence_performance():

    
    print(f"\n{'='*80}")

    print(f"{'='*80}\n")
    

    test_lengths = [100, 500, 1000, 2000]
    
    for length in test_lengths:

        print("-" * 40)
        

        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        amino_sequence = ''.join(np.random.choice(list(amino_acids), length))
        

        pi_accessibility = generate_test_distribution(length, num_codons=6)
        

        optimizer_standard = SADOOptimizer(
            species='ecoli_bl21de3',
            amino_acid_sequence=amino_sequence
        )
        
        start_time = time.time()
        dist_standard, meta_standard = optimizer_standard.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False
        )
        time_standard = time.time() - start_time
        

        optimizer_enhanced = SADOOptimizer(
            species='ecoli_bl21de3',
            amino_acid_sequence=amino_sequence
        )
        
        start_time = time.time()
        dist_enhanced, meta_enhanced = optimizer_enhanced.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=True,
            use_difference_driven=True
        )
        time_enhanced = time.time() - start_time
        






def main():
    """主测试函数"""
    

    results = test_enhanced_sado('O15263')
    

    print(f"\n{'='*80}")
    print(f"性能比较总结")
    print(f"{'='*80}\n")
    
    print(f"{'配置':<40} {'平均CAI':<12} {'平均时间(ms)':<15} {'多样性'}")
    print("-" * 80)
    
    for result in results:
        print(f"{result['config']:<40} "
              f"{result['avg_cai']:.4f}±{result['std_cai']:.4f}  "
              f"{result['avg_time_ms']:>6.1f}±{result['std_time_ms']:.1f}  "
              f"{result['diversity']['unique_ratio']*100:>6.1f}%")
    

    best_cai = max(results, key=lambda x: x['avg_cai'])
    best_speed = min(results, key=lambda x: x['avg_time_ms'])
    
    print(f"\n🏆 最佳CAI: {best_cai['config']} (CAI={best_cai['avg_cai']:.4f})")
    print(f"⚡ 最快速度: {best_speed['config']} ({best_speed['avg_time_ms']:.1f}ms)")
    

    compare_long_sequence_performance()
    
    print(f"\n{'='*80}")
    print(f"✅ 测试完成!")
    print(f"{'='*80}\n")


if __name__ == '__main__':
    main()
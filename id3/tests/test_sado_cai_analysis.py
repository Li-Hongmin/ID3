"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_cai_analysis')


def test_random_vs_sado():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    


    logger.info("-" * 40)
    
    random_cais = []
    sado_cais = []
    
    for i in range(10):

        pi_random = torch.rand(seq_len, num_codons, device=device)
        pi_random = pi_random * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_random[pos] = pi_random[pos] / pi_random[pos].sum()
        

        _, metadata = optimizer.optimize(
            pi_accessibility=pi_random,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        sado_cais.append(metadata['final_cai'])
        


        selected_indices = torch.argmax(pi_random, dim=1)
        random_cai = compute_cai_from_selection(selected_indices, optimizer)
        random_cais.append(random_cai)
    


    


    logger.info("-" * 40)
    
    uniform_cais = []
    sado_uniform_cais = []
    
    for i in range(10):

        pi_uniform = torch.ones(seq_len, num_codons, device=device) / num_codons
        pi_uniform = pi_uniform * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                valid_count = valid_mask[pos].sum()
                pi_uniform[pos] = valid_mask[pos].float() / valid_count
        

        _, metadata = optimizer.optimize(
            pi_accessibility=pi_uniform,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        sado_uniform_cais.append(metadata['final_cai'])
        

        random_selection = []
        for pos in range(seq_len):
            valid_indices = [i for i in range(num_codons) if valid_mask[pos, i]]
            if valid_indices:
                random_selection.append(np.random.choice(valid_indices))
            else:
                random_selection.append(0)
        random_selection = torch.tensor(random_selection)
        uniform_cai = compute_cai_from_selection(random_selection, optimizer)
        uniform_cais.append(uniform_cai)
    


    


    logger.info("-" * 40)
    
    target_cais = [0.5, 0.6, 0.7, 0.8, 0.9]
    pi_test = torch.rand(seq_len, num_codons, device=device)
    pi_test = pi_test * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
    
    for target in target_cais:
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=target,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        achieved_cai = metadata['final_cai']
        satisfied = "✅" if achieved_cai >= target else "❌"

    


    logger.info("-" * 40)
    
    gammas = [0.0, 0.3, 0.5, 0.6, 0.8, 1.0]
    for gamma in gammas:
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        logger.info(f"gamma={gamma:.1f}: CAI={metadata['final_cai']:.4f}")


def compute_cai_from_selection(indices, optimizer):
    """从密码子选择计算CAI值"""
    cai_product = 1.0
    
    for pos, idx in enumerate(indices):
        if pos < len(optimizer.codon_choices):
            choices = optimizer.codon_choices[pos]

            for choice in choices:
                if choice.get('original_local_index') == idx.item():
                    cai_product *= choice['weight']
                    break
    
    seq_len = len(indices)
    cai = cai_product ** (1.0 / seq_len) if seq_len > 0 else 0.0
    return cai


def analyze_sado_strategy():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    


    logger.info("-" * 40)
    

        if i < len(optimizer.codon_choices):
            choices = optimizer.codon_choices[i]


                logger.info(f"  {choice['codon']}: weight={choice['weight']:.4f}")


def main():
    """主函数"""
    test_random_vs_sado()
    analyze_sado_strategy()
    
    logger.info("\n" + "="*80)
    logger.info("分析总结")
    logger.info("="*80)
    logger.info("SADO倾向于选择高CAI权重的密码子，这是因为：")
    logger.info("1. gamma参数默认0.6，较大程度偏向CAI")
    logger.info("2. 目标CAI设为0.8时，算法会积极优化")
    logger.info("3. E.coli的密码子使用偏好明显，最优密码子权重很高")


if __name__ == "__main__":
    main()
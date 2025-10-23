#!/usr/bin/env python3
"""


"""

import torch
import logging
from id3.constraints.lagrangian import LagrangianConstraint
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
from id3.experiments.utils.data_loader import ProteinDataLoader

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_cai_gradient_flow():

    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    

    data_loader = ProteinDataLoader()
    protein_info = data_loader.get_protein_info('O15263')
    

    constraint = LagrangianConstraint(

        protein_seq=protein_info['sequence'],

        lambda_value=0.01,
        enable_cai=True,
        cai_target=0.8,
        cai_lambda=0.1,
        cai_species='ecoli_bl21de3',
        device=device
    )
    

    deepraccess = DeepRaccessID3Wrapper(device=device)
    

    optimizer = torch.optim.Adam(constraint.parameters(), lr=0.01)
    


    
    for iteration in range(5):

        result = constraint.forward(alpha=1.0, beta=0.0, tau=1.0)
        rna_seq = result['rna_sequence']
        

        if rna_seq.dim() == 2:

            rna_seq_one_hot = torch.zeros(rna_seq.shape[0], rna_seq.shape[1], 4, device=device)
            rna_seq_one_hot.scatter_(2, rna_seq.long().unsqueeze(2), 1.0)
            rna_seq = rna_seq_one_hot
        elif rna_seq.dim() != 3:
            raise ValueError(f"Unexpected rna_seq dimension: {rna_seq.dim()}")
        


        
        # UTR5 one-hot
        utr5_indices = torch.tensor([rna_to_index[n] for n in protein_info['utr5']], 
                                   dtype=torch.long, device=device)
        utr5_one_hot = torch.zeros(1, len(utr5_indices), 4, device=device)
        utr5_one_hot[0].scatter_(1, utr5_indices.unsqueeze(1), 1.0)
        
        # UTR3 one-hot
        utr3_indices = torch.tensor([rna_to_index[n] for n in protein_info['utr3']], 
                                   dtype=torch.long, device=device)
        utr3_one_hot = torch.zeros(1, len(utr3_indices), 4, device=device)
        utr3_one_hot[0].scatter_(1, utr3_indices.unsqueeze(1), 1.0)
        

        full_sequence = torch.cat([utr5_one_hot, rna_seq, utr3_one_hot], dim=1)
        

        accessibility = deepraccess.compute_atg_window_accessibility(
            full_sequence,
            atg_position=len(protein_info['utr5']),
            discrete=False
        )
        

        loss_components = constraint.compute_total_loss(
            accessibility,
            result.get('constraint_penalty', torch.tensor(0.0, device=device)),
            probabilities=result.get('probabilities'),
            cai_metadata=result.get('cai_metadata')
        )
        
        total_loss = loss_components['total_loss']
        


        acc_loss = loss_components.get('accessibility_loss', torch.tensor(0.0))
        cai_loss = loss_components.get('cai_loss', torch.tensor(0.0))
        cons_penalty = loss_components.get('constraint_penalty', torch.tensor(0.0))
        
        logger.info(f"  - Accessibility Loss: {acc_loss.item() if hasattr(acc_loss, 'item') else acc_loss:.4f}")
        logger.info(f"  - CAI Loss: {cai_loss.item() if hasattr(cai_loss, 'item') else cai_loss:.4f}")
        logger.info(f"  - Constraint Penalty: {cons_penalty.item() if hasattr(cons_penalty, 'item') else cons_penalty:.4f}")
        logger.info(f"  - Total Loss: {total_loss.item():.4f}")
        

        if 'cai_loss' in loss_components and loss_components['cai_loss'] > 0:

        else:

        

        optimizer.zero_grad()
        total_loss.backward()
        

        has_gradients = False
        for name, param in constraint.named_parameters():
            if param.grad is not None and param.grad.abs().max() > 0:
                has_gradients = True
                break
        
        if has_gradients:

        else:

        

        optimizer.step()
    

    return True

if __name__ == "__main__":
    test_cai_gradient_flow()
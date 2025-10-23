#!/usr/bin/env python3
"""

"""

import torch
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

sys.path.insert(0, str(Path(__file__).parent))

def test_cai_optimization():

    try:
        from id3.constraints.lagrangian import LagrangianConstraint
        from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper
        
        DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        

        deepraccess = DeepRaccessID3Wrapper()
        

        test_sequence = "ATGGCTAGCTAGCTAGCTAGC"  # 21nt = 7 codons

        L = len(test_sequence)
        

        constraint = LagrangianConstraint(
            L=L,
            amino_acid_sequence=amino_acid_sequence,
            deepraccess=deepraccess,
            device=DEVICE,
            alpha=0.0
        )
        

        cai_metadata = {
            'target': 0.8,
            'lambda_cai': 1.0,
            'species': 'ecoli_bl21de3'
        }
        


        

        with torch.no_grad():
            result = constraint.forward(alpha=0.0, tau=1.0)

            

            prob = result['probabilities']
            disc = result.get('enhanced_sequence', result.get('discrete_sequence'))
            


            

            model_loss = torch.tensor(1.0, device=DEVICE)
            constraint_penalty = torch.tensor(0.0, device=DEVICE)
            
            loss_dict = constraint.compute_total_loss(
                model_loss=model_loss,
                constraint_penalty=constraint_penalty,
                probabilities=prob,
                enhanced_sequence=disc,
                cai_metadata=cai_metadata
            )
            

            for key, value in loss_dict.items():
                if isinstance(value, torch.Tensor):
                    logger.info(f"  {key}: {value.item():.4f}")
                else:
                    logger.info(f"  {key}: {value}")
        

        return True
        
    except Exception as e:

        import traceback
        logger.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    success = test_cai_optimization()
    sys.exit(0 if success else 1)
#!/usr/bin/env python3
"""


"""

import torch
import torch.optim as optim
import numpy as np
import sys
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.lagrangian import LagrangianConstraint
from id3.experiments.analysis.result_parser import ExperimentResultParser

def simulate_experiment_process():

    print("=" * 80)

    print("=" * 80)
    

    protein_seq = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
    



    
    try:

        constraint = LagrangianConstraint(
            amino_acid_sequence=protein_seq,
            batch_size=1,
            initial_lambda=0.01,
            adaptive_lambda=True,
            lambda_lr=0.01,
            lambda_max=10.0,
            device=torch.device('cpu'),
            enable_cai=True,
            cai_target=0.8,
            cai_weight=0.1,  # lambda_cai
            verbose=False
        )
        



        print(f"  Lambda CAI: {constraint.lambda_cai}")
        

        seq_length = len(protein_seq) * 3
        theta = torch.randn(1, seq_length, 4, requires_grad=True)
        

        optimizer = optim.Adam([theta], lr=0.05)
        




        

        iteration_results = []
        

        for iteration in range(10):
            optimizer.zero_grad()
            


            

            probabilities = result.get('probabilities')
            enhanced_sequence = result.get('enhanced_sequence')
            discrete_cai = result.get('discrete_cai_value')
            cai_metadata = result.get('cai_metadata', {})
            


            constraint_penalty = result.get('constraint_penalty', torch.tensor(0.1))
            

            loss_components = constraint.compute_total_loss(
                model_loss=mock_accessibility,
                constraint_penalty=constraint_penalty,
                probabilities=probabilities,
                enhanced_sequence=enhanced_sequence
            )
            
            total_loss = loss_components['total_loss']
            

            iter_result = {
                'iteration': iteration,
                'total_loss': total_loss.item(),
                'accessibility': mock_accessibility.item(),
                'constraint_penalty': constraint_penalty.item() if hasattr(constraint_penalty, 'item') else constraint_penalty,
                'discrete_cai_value': discrete_cai,
                'eval_cai': loss_components.get('eval_cai', None),
                'ecai_value': loss_components.get('ecai_value', torch.tensor(0)).item(),
                'cai_metadata_final': cai_metadata.get('final_cai', None),
                'cai_metadata_original': cai_metadata.get('original_cai', None),
            }
            
            iteration_results.append(iter_result)
            


            print(f"  Total Loss: {total_loss.item():.6f}")
            print(f"  Discrete CAI Value: {discrete_cai:.6f}" if isinstance(discrete_cai, (int, float)) else f"  Discrete CAI Value: {discrete_cai}")
            print(f"  Eval CAI: {iter_result['eval_cai']:.6f}" if isinstance(iter_result['eval_cai'], (int, float)) else f"  Eval CAI: {iter_result['eval_cai']}")
            print(f"  ECAI Value: {iter_result['ecai_value']:.6f}")
            print(f"  Metadata Final CAI: {iter_result['cai_metadata_final']:.6f}" if isinstance(iter_result['cai_metadata_final'], (int, float)) else f"  Metadata Final CAI: {iter_result['cai_metadata_final']}")
            

            if isinstance(iter_result['eval_cai'], (int, float)) and iter_result['eval_cai'] < 1e-6:

            
            if isinstance(discrete_cai, (int, float)) and discrete_cai < 1e-6:

            

            total_loss.backward()
            optimizer.step()
            

            if constraint.adaptive_lambda:
                constraint.update_lambda(constraint_penalty.item() if hasattr(constraint_penalty, 'item') else constraint_penalty)
        

        print(f"\n{'='*60}")

        print(f"{'='*60}")
        

        discrete_cai_values = [r['discrete_cai_value'] for r in iteration_results if isinstance(r['discrete_cai_value'], (int, float))]
        eval_cai_values = [r['eval_cai'] for r in iteration_results if isinstance(r['eval_cai'], (int, float))]
        ecai_values = [r['ecai_value'] for r in iteration_results]
        
        print(f"Discrete CAI Values: {[f'{v:.6f}' for v in discrete_cai_values]}")
        print(f"Eval CAI Values: {[f'{v:.6f}' for v in eval_cai_values]}")
        print(f"ECAI Values: {[f'{v:.6f}' for v in ecai_values]}")
        

        abnormal_eval_cai = [v for v in eval_cai_values if v < 1e-6]
        abnormal_discrete_cai = [v for v in discrete_cai_values if v < 1e-6]
        
        if abnormal_eval_cai:

        
        if abnormal_discrete_cai:

        
        if not abnormal_eval_cai and not abnormal_discrete_cai:


        
        return len(abnormal_eval_cai) > 0 or len(abnormal_discrete_cai) > 0, iteration_results
        
    except Exception as e:

        import traceback

        return False, []

def analyze_original_experiment():
    """分析原始实验数据，寻找模式"""
    print("\n" + "=" * 80)
    print("原始实验数据分析")
    print("=" * 80)
    
    try:

        result_file = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250908_122909_unified_cai_experiments/20250908_134200_P04637_lagrangian_01_seed42.json"
        
        parser = ExperimentResultParser(result_file)
        result = parser.load_experiment()
        
        print(f"原始实验配置:")
        print(f"  蛋白质: {result.get('protein_name', 'N/A')}")
        print(f"  约束类型: {result.get('constraint_type', 'N/A')}")
        print(f"  变体: {result.get('variant', 'N/A')}")
        print(f"  迭代数: {result.get('iterations', 'N/A')}")
        print(f"  学习率: {result.get('learning_rate', 'N/A')}")
        print(f"  CAI启用: {result.get('cai_enabled', 'N/A')}")
        print(f"  CAI目标: {result.get('cai_target', 'N/A')}")
        print(f"  Lambda CAI: {result.get('lambda_cai', 'N/A')}")
        

        trajectory = result.get('trajectory', {})
        if trajectory:
            ecai_values = trajectory.get('ecai_values', [])
            discrete_cai_values = trajectory.get('discrete_cai_values', [])
            
            print(f"\nTrajectory数据:")
            print(f"  ECAI Values: {ecai_values[:5]}... (前5个)")
            print(f"  Discrete CAI Values: {discrete_cai_values[:5]}... (前5个)")
            

            if ecai_values:
                ecai_min, ecai_max = min(ecai_values), max(ecai_values)
                print(f"  ECAI范围: [{ecai_min:.6f}, {ecai_max:.6f}]")
            
            if discrete_cai_values:
                discrete_min, discrete_max = min(discrete_cai_values), max(discrete_cai_values)
                print(f"  Discrete CAI范围: [{discrete_min:.2e}, {discrete_max:.2e}]")
                

                all_small = all(v < 1e-6 for v in discrete_cai_values)
                print(f"  所有discrete_cai_values都 < 1e-6: {all_small}")
                
                if all_small:
                    print(f"  🚨 确认：原始实验中所有discrete_cai_values都异常小")
                    print(f"  🔍 这证实了问题的存在")
        
        return True
        
    except Exception as e:
        print(f"❌ 原始实验分析失败: {str(e)}")
        return False

def main():



    

    reproduced, results = simulate_experiment_process()
    

    analyzed = analyze_original_experiment()
    

    print(f"\n{'='*80}")

    print(f"{'='*80}")
    
    if reproduced:



    else:



    
    if analyzed:



    






if __name__ == "__main__":
    main()
#!/usr/bin/env python3
"""

"""

import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(project_root)

import torch
import matplotlib.pyplot as plt
import numpy as np
from id3.constraints.lagrangian import LagrangianConstraint

def test_subgradient_convergence():

    

    print("=" * 60)
    

    constraint = LagrangianConstraint(
        amino_acid_sequence="MGKR",

        adaptive_lambda=True,

        lambda_max=5.0,
        device='cpu'
    )
    


    print()
    

    constraint_penalties = []
    

    for i in range(20):

        constraint_penalties.append(penalty)
    

    for i in range(30):

        constraint_penalties.append(penalty)
    

    for i in range(50):

        constraint_penalties.append(max(0, penalty))
    


    

    for i, penalty in enumerate(constraint_penalties):
        old_lambda = constraint.lambda_value
        constraint.update_lambda(penalty, tolerance)
        new_lambda = constraint.lambda_value
        
        if i % 10 == 0:
            step_size = constraint.lambda_lr / np.sqrt(constraint.lambda_iteration)
            subgrad = penalty - tolerance
            print(f"Iter {i:3d}: penalty={penalty:.4f}, λ: {old_lambda:.4f}→{new_lambda:.4f}, "
                  f"step={step_size:.4f}, subgrad={subgrad:.4f}")
    


    

    analysis = constraint.get_lambda_convergence_analysis()
    

    print("-" * 40)






    

    rm = analysis['robbins_monro_satisfied']
    print(f"  Σα_t → ∞: {rm['sum_infinite']} ({'✅' if rm['sum_infinite'] else '❌'})")
    print(f"  Σα_t² < ∞: {rm['square_sum_finite']} ({'✅' if rm['square_sum_finite'] else '❌'})")

    



    

    plt.figure(figsize=(15, 5))
    

    plt.subplot(1, 3, 1)
    plt.plot(analysis['lambda_trajectory'], 'b-', linewidth=2)
    plt.title('λ Convergence Trajectory')
    plt.xlabel('Iteration')
    plt.ylabel('λ value')
    plt.grid(True, alpha=0.3)
    

    plt.subplot(1, 3, 2)
    plt.plot(constraint_penalties, 'r-', alpha=0.7, label='Constraint Penalty')
    plt.axhline(y=tolerance, color='g', linestyle='--', label=f'Tolerance={tolerance}')
    plt.title('Constraint Violation Over Time')
    plt.xlabel('Iteration')
    plt.ylabel('Penalty')
    plt.legend()
    plt.grid(True, alpha=0.3)
    

    plt.subplot(1, 3, 3)
    step_sizes = np.array(analysis['step_size_sum']) # This needs to be fixed - should be individual step sizes
    iterations = np.arange(1, len(analysis['lambda_trajectory']) + 1)
    theoretical_steps = constraint.lambda_lr / np.sqrt(iterations)
    plt.plot(theoretical_steps, 'g-', label='Theoretical α_t = α_0/√t')
    plt.title('Step Size Decay')
    plt.xlabel('Iteration')
    plt.ylabel('Step Size')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('subgradient_convergence_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    

    


    if analysis['recent_subgradient_trend'] < 0.05:

    else:

        
    if analysis['lambda_stability'] < 0.1:

    else:


if __name__ == "__main__":
    test_subgradient_convergence()
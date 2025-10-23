#!/usr/bin/env python3
"""
Focused Lagrangian CAI Fix Validation

"""

import torch
import numpy as np
import sys
import time
import json

sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.lagrangian import LagrangianConstraint

def test_lagrangian_cai_fix():


    print("=" * 80)
    

    test_proteins = {
        'short': 'MKHELM',
        'medium': 'MEEPQSDPSVEPPLSQETFSDLWKLL'
    }
    
    results = {}
    
    for name, protein_seq in test_proteins.items():


        
        try:

            constraint = LagrangianConstraint(
                protein_seq,
                enable_cai=True,
                cai_target=0.8,
                cai_lambda=0.1,
                verbose=False
            )
            

            

            test_modes = [



            ]
            
            protein_results = {}
            
            for mode_name, alpha, tau in test_modes:

                
                start_time = time.time()
                result = constraint.forward(alpha=alpha, tau=tau)
                end_time = time.time()
                

                discrete_cai = result.get('discrete_cai_value', None)
                cai_metadata = result.get('cai_metadata', {})
                final_cai = cai_metadata.get('final_cai', None) if cai_metadata else None
                total_loss = result.get('total_loss', None)
                cai_loss = result.get('cai_loss', None)
                

                success = True
                issues = []
                
                if discrete_cai is None:
                    success = False


                    success = False


                    success = False

                
                if final_cai is not None and abs(discrete_cai - final_cai) > 0.2:

                
                if total_loss is None or not np.isfinite(total_loss):
                    success = False

                

                mode_result = {
                    'success': success,
                    'issues': issues,
                    'discrete_cai': discrete_cai,
                    'final_cai': final_cai,
                    'total_loss': total_loss,
                    'cai_loss': cai_loss,
                    'execution_time': end_time - start_time
                }
                
                protein_results[mode_name] = mode_result
                

                if success:

                    print(f"       discrete_cai: {discrete_cai:.6f}")
                    print(f"       final_cai: {final_cai:.6f if final_cai else 'N/A'}")
                    print(f"       total_loss: {total_loss:.6f if total_loss else 'N/A'}")

                else:

                    for issue in issues:
                        print(f"       - {issue}")
            
            results[name] = protein_results
            
        except Exception as e:

            results[name] = {'error': str(e)}
    


    print("=" * 80)
    
    total_tests = 0
    successful_tests = 0
    cai_fix_confirmed = True
    
    for protein_name, protein_results in results.items():
        if 'error' in protein_results:

            cai_fix_confirmed = False
            continue
            

        
        for mode_name, mode_result in protein_results.items():
            total_tests += 1
            if mode_result['success']:
                successful_tests += 1
                status = "✅"
            else:
                status = "❌"

                    cai_fix_confirmed = False
            
            discrete_cai = mode_result.get('discrete_cai', 0)
            print(f"  {status} {mode_name:12}: CAI={discrete_cai:.6f}")
            
            if mode_result.get('issues'):
                for issue in mode_result['issues']:
                    print(f"      - {issue}")
    
    success_rate = successful_tests / total_tests * 100 if total_tests > 0 else 0
    



    
    if cai_fix_confirmed and success_rate >= 80:



    else:

        if success_rate < 80:

        if not cai_fix_confirmed:


    

    detailed_report = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'test_summary': {
            'total_tests': total_tests,
            'successful_tests': successful_tests,
            'success_rate': success_rate,
            'cai_fix_confirmed': cai_fix_confirmed,
            'final_status': final_status
        },
        'detailed_results': results
    }
    
    with open('lagrangian_cai_fix_validation.json', 'w') as f:
        json.dump(detailed_report, f, indent=2, default=str)
    

    
    return cai_fix_confirmed and success_rate >= 80

def test_dynamic_lambda_adjustment():
    """测试动态lambda_cai调整功能"""
    print(f"\n🔧 动态lambda_cai调整验证")
    print("-" * 60)
    
    protein_seq = 'MKHELM'
    
    try:

        constraint = LagrangianConstraint(
            protein_seq,
            enable_cai=True,
            cai_target=0.8,
            adaptive_lambda_cai=True,
            lambda_cai_min=0.01,
            lambda_cai_max=1.0,
            lambda_cai_adjustment_factor=1.5,
            verbose=False
        )
        
        print(f"✅ 启用动态调整的Lagrangian约束创建成功")
        

        lambda_history = []
        cai_history = []
        
        for i in range(5):
            result = constraint.forward(alpha=1.0, tau=1.0)
            
            current_lambda = getattr(constraint, 'cai_lambda', None)
            discrete_cai = result.get('discrete_cai_value', 0)
            
            lambda_history.append(current_lambda)
            cai_history.append(discrete_cai)
            
            print(f"  迭代 {i+1}: λ={current_lambda:.4f if current_lambda else 'None'}, CAI={discrete_cai:.6f}")
            

            if hasattr(constraint, '_update_lambda_cai'):
                constraint._update_lambda_cai(discrete_cai)
        

        lambda_variance = np.var(lambda_history) if len(lambda_history) > 1 else 0
        
        if lambda_variance > 1e-10:
            print(f"✅ 动态调整正常工作 (λ方差: {lambda_variance:.2e})")
            return True
        else:
            print(f"❌ 动态调整未生效 (λ方差: {lambda_variance:.2e})")
            return False
            
    except Exception as e:
        print(f"❌ 动态调整测试失败: {str(e)}")
        return False

def main():



    

    cai_fix_ok = test_lagrangian_cai_fix()
    

    dynamic_ok = test_dynamic_lambda_adjustment()
    


    print("=" * 80)
    
    if cai_fix_ok:

    else:

    
    if dynamic_ok:

    else:

    
    overall_success = cai_fix_ok and dynamic_ok
    
    if overall_success:


        return 0
    else:

        return 1

if __name__ == "__main__":
    exit(main())
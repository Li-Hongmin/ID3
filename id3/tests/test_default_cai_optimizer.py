#!/usr/bin/env python3
"""

"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

from id3.optimizers.cai import DefaultCAIOptimizer, BinarySearchCAIOptimizer, SADOOptimizer


def test_default_optimizer():

    
    print("="*60)

    print("="*60)
    




    

    assert DefaultCAIOptimizer is BinarySearchCAIOptimizer, \

    

    


    optimizer = DefaultCAIOptimizer(species='ecoli_bl21de3')
    


    
    assert isinstance(optimizer, BinarySearchCAIOptimizer), \

    

    




    

    
    return True


if __name__ == "__main__":

    
    success = test_default_optimizer()
    
    if success:
        print("\n" + "="*60)

        print("="*60)
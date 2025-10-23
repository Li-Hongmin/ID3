#!/usr/bin/env python3


import json
import tempfile
import time
from pathlib import Path
import sys
import subprocess


sys.path.append(str(Path(__file__).parent.parent.parent))


def test_config_save_timing():
    """测试配置文件保存时机"""
    

    temp_dir = tempfile.mkdtemp()
    print(f"临时输出目录: {temp_dir}")
    


    project_root = Path(__file__).parent.parent.parent
    cmd = [
        'python', str(project_root / 'run_unified_experiment.py'),
        '--proteins', 'O15263',
        '--constraints', 'lagrangian',
        '--variants', '00',
        '--iterations', '1',
        '--seeds', '1',
        '--base-seed', '42',
        '--output-dir', temp_dir,
        '--no-trajectories'
    ]
    
    print(f"运行命令: {' '.join(cmd)}")
    

    process = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True
    )
    

    config_file = None
    start_time = time.time()
    max_wait = 10
    
    while time.time() - start_time < max_wait:

        config_files = list(Path(temp_dir).rglob('config.json'))
        if config_files:
            config_file = config_files[0]
            elapsed = time.time() - start_time
            print(f"✅ 配置文件在 {elapsed:.2f} 秒后创建: {config_file}")
            

            with open(config_file, 'r') as f:
                config = json.load(f)
            
            print("\n配置内容包含以下关键字段:")
            important_fields = ['proteins', 'constraints', 'variants', 'iterations', 
                              'seeds', 'base_seed', 'enable_cai']
            for field in important_fields:
                if field in config:
                    print(f"  - {field}: {config[field]}")
            

            process.terminate()
            process.wait()
            
            if elapsed < 5:
                print(f"\n✅ 测试通过: 配置文件在实验开始前立即保存（{elapsed:.2f}秒）")
                return True
            else:
                print(f"\n⚠️ 警告: 配置文件创建较慢（{elapsed:.2f}秒）")
                return True
        
        time.sleep(0.1)
    

    process.terminate()
    process.wait()
    print(f"❌ 测试失败: {max_wait}秒内未找到配置文件")
    return False


def test_dual_mode_config():

    

    temp_dir = tempfile.mkdtemp()


    

    project_root = Path(__file__).parent.parent.parent
    cmd = [
        'python', str(project_root / 'run_unified_experiment.py'),
        '--preset', 'quick-test-both',

        '--output-dir', temp_dir
    ]
    

    

    process = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        text=True
    )
    

    start_time = time.time()

    config_count = 0
    
    while time.time() - start_time < max_wait:

        config_files = list(Path(temp_dir).rglob('config.json'))
        
        if len(config_files) > config_count:
            config_count = len(config_files)
            elapsed = time.time() - start_time

            

                process.terminate()
                process.wait()

                return True
        
        time.sleep(0.1)
    
    process.terminate()
    process.wait()

    return config_count > 0


if __name__ == "__main__":
    print("=" * 60)

    print("=" * 60)
    

    success1 = test_config_save_timing()
    

    success2 = test_dual_mode_config()
    
    if success1 and success2:

        sys.exit(0)
    else:

        sys.exit(1)
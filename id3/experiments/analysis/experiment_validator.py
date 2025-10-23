#!/usr/bin/env python3
"""



"""

import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd


class ExperimentValidator:

    
    def __init__(self):
        """初始化验证器"""
        self.expected_experiment_count = 144  # 12x12 experiments
    
    def check_experiment_completeness(self, base_dir: str = "results") -> Dict:
        """

        
        Args:

            
        Returns:

        """
        results_dir = Path(base_dir)
        
        if not results_dir.exists():
            return {
                'error': f"目录不存在: {base_dir}",
                'complete_dirs': [],
                'incomplete_dirs': [],
                'cai_dirs': [],
                'access_dirs': []
            }
        

        complete_144_dirs = []
        incomplete_dirs = []
        cai_dirs = []
        access_dirs = []
        
        for dir_path in sorted(results_dir.glob("*unified*")):
            if dir_path.is_dir():

                json_files = list(dir_path.glob("*.json"))
                experiment_files = [f for f in json_files 
                                  if not f.name.startswith(('config', 'progress', 'summary'))]
                
                count = len(experiment_files)
                dir_name = dir_path.name
                

                if 'cai' in dir_name:
                    cai_dirs.append((dir_name, count))
                else:
                    access_dirs.append((dir_name, count))
                

                if count == self.expected_experiment_count:
                    complete_144_dirs.append(dir_name)
                else:
                    incomplete_dirs.append((dir_name, count))
        
        return {
            'complete_dirs': complete_144_dirs,
            'incomplete_dirs': incomplete_dirs,
            'cai_dirs': cai_dirs,
            'access_dirs': access_dirs,
            'total_complete': len(complete_144_dirs),
            'total_incomplete': len(incomplete_dirs)
        }
    
    def check_cai_status(self, dir_path: Path) -> Dict:
        """

        
        Args:

            
        Returns:

        """
        result = {
            'has_cai_config': False,
            'cai_config': {},
            'has_cai_fields': False,
            'sample_cai_values': [],
            'cai_field_stats': {}
        }
        

        config_file = dir_path / "config.json"
        if config_file.exists():
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)
                    result['has_cai_config'] = config.get('enable_cai', False)
                    result['cai_config'] = {
                        'enable_cai': config.get('enable_cai', 'NOT FOUND'),
                        'cai_target': config.get('cai_target', 'NOT FOUND'),
                        'lambda_cai': config.get('lambda_cai', 'NOT FOUND'),
                        'species': config.get('species', 'NOT FOUND')
                    }
            except Exception as e:
                result['config_error'] = str(e)
        

        json_files = list(dir_path.glob("*.json"))
        experiment_files = [f for f in json_files 
                           if not f.name.startswith(('config', 'progress', 'summary'))]
        
        if experiment_files:

            sample_files = experiment_files[:3]
            cai_field_counts = {
                'final_ecai': 0,
                'initial_ecai': 0,
                'cai_target': 0,
                'cai_optimization': 0,
                'ecai_values': 0,
                'cai_loss': 0
            }
            
            for file_path in sample_files:
                try:
                    with open(file_path, 'r') as f:
                        content = f.read(10000)
                        

                        for field in cai_field_counts:
                            if f'"{field}"' in content:
                                cai_field_counts[field] += 1
                        

                        import re
                        match = re.search(r'"final_ecai":\s*([\d.]+)', content)
                        if match:
                            result['sample_cai_values'].append(float(match.group(1)))
                
                except Exception:
                    pass
            

            result['has_cai_fields'] = any(count > 0 for count in cai_field_counts.values())
            result['cai_field_stats'] = cai_field_counts
        
        return result
    
    def verify_results_full(self) -> Dict:
        """

        
        Returns:

        """
        results_full = Path("results_full")
        
        if not results_full.exists():
            return {
                'exists': False,
                'directories': [],
                'total_size_gb': 0,
                'all_complete': False
            }
        
        directories = []
        total_size = 0
        all_complete = True
        
        for dir_path in sorted(results_full.glob("*")):
            if dir_path.is_dir():

                json_files = list(dir_path.glob("*.json"))
                experiment_files = [f for f in json_files 
                                  if not f.name.startswith(('config', 'progress', 'summary'))]
                
                count = len(experiment_files)
                dir_name = dir_path.name
                

                dir_size = sum(f.stat().st_size for f in experiment_files)
                size_gb = dir_size / (1024 * 1024 * 1024)
                total_size += dir_size
                

                exp_type = 'CAI' if 'cai' in dir_name else 'Access'
                

                if dir_name == "20250910_004126_unified_access_experiments":
                    exp_type = "Access(实际是CAI)"
                
                is_complete = count == self.expected_experiment_count
                if not is_complete:
                    all_complete = False
                
                directories.append({
                    'name': dir_name,
                    'type': exp_type,
                    'count': count,
                    'is_complete': is_complete,
                    'size_gb': size_gb
                })
        
        return {
            'exists': True,
            'directories': directories,
            'total_size_gb': total_size / (1024 * 1024 * 1024),
            'all_complete': all_complete,
            'total_dirs': len(directories)
        }
    
    def print_completeness_report(self, stats: Dict):


        print("-" * 80)
        

        if stats['complete_dirs']:


                print(f"   • {dir_name}")
            if len(stats['complete_dirs']) > 5:

        

        if stats['incomplete_dirs']:

            for dir_name, count in sorted(stats['incomplete_dirs'], key=lambda x: -x[1])[:5]:
                print(f"   • {dir_name}: {count}/144")
        

        cai_complete = sum(1 for _, count in stats['cai_dirs'] if count == 144)
        access_complete = sum(1 for _, count in stats['access_dirs'] if count == 144)
        



        



        
        print("\n" + "=" * 80)

        print("=" * 80)
    
    def print_cai_status_report(self, dir_path: Path, status: Dict):
        """打印CAI状态报告"""
        print(f"\n🔍 检查实验CAI状态: {dir_path.name}")
        print("=" * 80)
        

        if 'cai_config' in status:
            print("\n📋 配置文件:")
            for key, value in status['cai_config'].items():
                print(f"  • {key}: {value}")
        

        if status['has_cai_fields']:
            print("\n✅ 检测到CAI相关字段:")
            for field, count in status['cai_field_stats'].items():
                if count > 0:
                    print(f"  • {field}: {count}/3 个文件包含")
        else:
            print("\n❌ 未检测到CAI相关字段")
        

        if status['sample_cai_values']:
            print(f"\n📊 样本CAI值: {status['sample_cai_values']}")
        

        print("\n📝 结论:")
        if status['has_cai_config'] and status['has_cai_fields']:
            print("  ✅ 这是一个包含CAI优化的实验")
        elif status['has_cai_config'] and not status['has_cai_fields']:
            print("  ⚠️ 配置启用了CAI但结果中没有CAI字段")
        else:
            print("  ❌ 这是一个纯Access优化实验")
    
    def print_results_full_report(self, stats: Dict):

        if not stats['exists']:

            return
        

        print("=" * 80)
        
        for dir_info in stats['directories']:

            print(f"{dir_info['name']:50} {dir_info['type']:15} {status:15} {dir_info['size_gb']:8.2f} GB")
        
        print("=" * 80)

        
        if stats['all_complete']:

        else:


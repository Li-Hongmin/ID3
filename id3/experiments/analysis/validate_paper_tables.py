#!/usr/bin/env python3
"""
Validate paper experiment tables for data completeness and correctness.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from typing import Dict, List, Tuple


PROTEIN_LIST_PAPER = [
    'O15263', 'P00004', 'P01308', 'P01825',
    'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',
    'P31417', 'P42212', 'P61626', 'P99999'
]

class PaperTableValidator:
    """Validate paper experiment tables."""
    
    def __init__(self):
        self.expected_proteins = PROTEIN_LIST_PAPER
        self.expected_constraints = ['L', 'A', 'C']
        self.expected_variants = ['00', '01', '10', '11']
        self.validation_results = {}
        
    def validate_completeness(self, csv_path: str) -> Dict:
        """Validate data completeness in CSV table."""
        df = pd.read_csv(csv_path)
        results = {
            'total_rows': len(df),
            'expected_rows': len(self.expected_proteins),
            'missing_proteins': [],
            'extra_proteins': [],
            'column_count': len(df.columns)
        }
        
        # Check protein coverage
        actual_proteins = set(df.iloc[:, 0].values) if 'Protein' in df.columns else set(df.index)
        expected_proteins = set(self.expected_proteins)
        
        results['missing_proteins'] = list(expected_proteins - actual_proteins)
        results['extra_proteins'] = list(actual_proteins - expected_proteins)
        results['complete'] = len(results['missing_proteins']) == 0
        
        return results
    
    def validate_accessibility_values(self, csv_path: str) -> Dict:
        """Validate accessibility values are in reasonable range."""
        df = pd.read_csv(csv_path)
        results = {
            'min_value': float('inf'),
            'max_value': float('-inf'),
            'mean_value': 0,
            'suspicious_values': [],
            'negative_values': [],
            'extreme_values': []
        }
        
        # Get all numeric columns (skip first column which is protein names)
        numeric_cols = df.columns[1:] if 'Protein' in df.columns else df.columns
        
        all_values = []
        for col in numeric_cols:
            if col in df.columns:
                values = pd.to_numeric(df[col], errors='coerce').dropna()
                all_values.extend(values.tolist())
                
                # Check for suspicious values
                for idx, val in enumerate(values):
                    protein = df.iloc[idx, 0] if 'Protein' in df.columns else df.index[idx]
                    
                    if val < 0:
                        results['negative_values'].append({
                            'protein': protein,
                            'column': col,
                            'value': val
                        })
                    
                    # Accessibility should typically be between 0.5 and 3.0 kcal/mol
                    if val < 0.3 or val > 5.0:
                        results['extreme_values'].append({
                            'protein': protein,
                            'column': col,
                            'value': val
                        })
                    
                    # Values below 0.5 are suspicious (might be loss instead of validation)
                    if val < 0.5:
                        results['suspicious_values'].append({
                            'protein': protein,
                            'column': col,
                            'value': val
                        })
        
        if all_values:
            results['min_value'] = min(all_values)
            results['max_value'] = max(all_values)
            results['mean_value'] = np.mean(all_values)
            results['std_value'] = np.std(all_values)
        
        return results
    
    def validate_cai_values(self, csv_path: str) -> Dict:
        """Validate CAI values are in valid range [0, 1]."""
        df = pd.read_csv(csv_path)
        results = {
            'out_of_range': [],
            'target_met': [],
            'target_missed': [],
            'cai_target': 0.8
        }
        
        # Look for CAI columns
        cai_cols = [col for col in df.columns if 'CAI' in col]
        
        for col in cai_cols:
            if col in df.columns:
                values = pd.to_numeric(df[col], errors='coerce').dropna()
                
                for idx, val in enumerate(values):
                    protein = df.iloc[idx, 0] if 'Protein' in df.columns else df.index[idx]
                    
                    # CAI should be between 0 and 1
                    if val < 0 or val > 1:
                        results['out_of_range'].append({
                            'protein': protein,
                            'column': col,
                            'value': val
                        })
                    
                    # Check if target is met (for With_Penalty_CAI columns)
                    if 'With_Penalty' in col:
                        if abs(val - results['cai_target']) < 0.05:  # Within 5% of target
                            results['target_met'].append({
                                'protein': protein,
                                'column': col,
                                'value': val
                            })
                        else:
                            results['target_missed'].append({
                                'protein': protein,
                                'column': col,
                                'value': val,
                                'deviation': abs(val - results['cai_target'])
                            })
        
        return results
    
    def validate_all_tables(self, table_dir: str) -> Dict:
        """Validate all generated tables."""
        table_dir = Path(table_dir)
        all_results = {}
        
        # Table 1: Accessibility Performance
        table1_path = table_dir / 'table1_accessibility_performance.csv'
        if table1_path.exists():
            print("\n" + "="*80)
            print("VALIDATING TABLE 1: Accessibility Performance")
            print("="*80)
            
            completeness = self.validate_completeness(str(table1_path))
            print(f"\nCompleteness Check:")
            print(f"  - Total rows: {completeness['total_rows']}/{completeness['expected_rows']}")
            print(f"  - Complete: {completeness['complete']}")
            if completeness['missing_proteins']:
                print(f"  - Missing proteins: {completeness['missing_proteins']}")
            
            values = self.validate_accessibility_values(str(table1_path))
            print(f"\nValue Range Check:")
            print(f"  - Range: {values['min_value']:.4f} - {values['max_value']:.4f} kcal/mol")
            print(f"  - Mean ± Std: {values['mean_value']:.4f} ± {values['std_value']:.4f}")
            
            if values['negative_values']:
                print(f"  ⚠️  Found {len(values['negative_values'])} negative values!")
            
            if values['extreme_values']:
                print(f"  ⚠️  Found {len(values['extreme_values'])} extreme values (outside 0.3-5.0 range)")
                for item in values['extreme_values'][:3]:  # Show first 3
                    print(f"      - {item['protein']}, {item['column']}: {item['value']:.4f}")
            
            if values['suspicious_values']:
                print(f"  ⚠️  Found {len(values['suspicious_values'])} suspicious values (<0.5, might be loss not validation)")
                for item in values['suspicious_values'][:3]:  # Show first 3
                    print(f"      - {item['protein']}, {item['column']}: {item['value']:.4f}")
            
            all_results['table1'] = {'completeness': completeness, 'values': values}
        
        # Table 2: CAI Comparison
        table2_path = table_dir / 'table2_cai_comparison.csv'
        if table2_path.exists():
            print("\n" + "="*80)
            print("VALIDATING TABLE 2: CAI Comparison")
            print("="*80)
            
            completeness = self.validate_completeness(str(table2_path))
            print(f"\nCompleteness Check:")
            print(f"  - Total rows: {completeness['total_rows']}/{completeness['expected_rows']}")
            print(f"  - Complete: {completeness['complete']}")
            
            # Check both accessibility and CAI values
            access_values = self.validate_accessibility_values(str(table2_path))
            print(f"\nAccessibility Value Range:")
            print(f"  - Range: {access_values['min_value']:.4f} - {access_values['max_value']:.4f}")
            print(f"  - Mean ± Std: {access_values['mean_value']:.4f} ± {access_values['std_value']:.4f}")
            
            cai_values = self.validate_cai_values(str(table2_path))
            if cai_values['out_of_range']:
                print(f"\n⚠️  CAI values out of [0,1] range: {len(cai_values['out_of_range'])}")
            
            if cai_values['target_met'] or cai_values['target_missed']:
                total_cai = len(cai_values['target_met']) + len(cai_values['target_missed'])
                if total_cai > 0:
                    success_rate = len(cai_values['target_met']) / total_cai * 100
                    print(f"\nCAI Target Achievement (0.8 ± 0.05):")
                    print(f"  - Success rate: {success_rate:.1f}%")
                    print(f"  - Met target: {len(cai_values['target_met'])}/{total_cai}")
                    
                    if cai_values['target_missed']:
                        avg_deviation = np.mean([x['deviation'] for x in cai_values['target_missed']])
                        print(f"  - Average deviation when missed: {avg_deviation:.4f}")
            
            all_results['table2'] = {
                'completeness': completeness,
                'access_values': access_values,
                'cai_values': cai_values
            }
        
        # Table 3: Summary Statistics
        table3_path = table_dir / 'table3_summary_statistics.csv'
        if table3_path.exists():
            print("\n" + "="*80)
            print("VALIDATING TABLE 3: Summary Statistics")
            print("="*80)
            
            df = pd.read_csv(table3_path)
            print(f"\nSummary Categories:")
            for idx, row in df.iterrows():
                category = row['Category'] if 'Category' in df.columns else row[0]
                mean_val = row['Mean'] if 'Mean' in df.columns else row[1]
                count = row['Count'] if 'Count' in df.columns else row[6]
                excellent = row['Excellent (<1.0)'] if 'Excellent (<1.0)' in df.columns else row[7]
                good = row['Good (<1.5)'] if 'Good (<1.5)' in df.columns else row[8]
                
                print(f"\n  {category}:")
                print(f"    - Mean: {mean_val:.4f} kcal/mol")
                print(f"    - Count: {count}")
                print(f"    - Excellent (<1.0): {excellent}/{count} ({excellent/count*100:.1f}%)")
                print(f"    - Good (<1.5): {good}/{count} ({good/count*100:.1f}%)")
            
            all_results['table3'] = {'data': df.to_dict()}
        
        # Overall summary
        print("\n" + "="*80)
        print("VALIDATION SUMMARY")
        print("="*80)
        
        all_complete = all([
            all_results.get('table1', {}).get('completeness', {}).get('complete', False),
            all_results.get('table2', {}).get('completeness', {}).get('complete', False)
        ])
        
        print(f"\n✅ Data Completeness: {'PASS' if all_complete else 'FAIL'}")
        
        # Check for critical issues
        critical_issues = []
        
        if 'table1' in all_results:
            if all_results['table1']['values']['negative_values']:
                critical_issues.append("Negative accessibility values found")
            if all_results['table1']['values']['suspicious_values']:
                critical_issues.append(f"{len(all_results['table1']['values']['suspicious_values'])} suspicious values (<0.5)")
        
        if 'table2' in all_results:
            if all_results['table2']['cai_values']['out_of_range']:
                critical_issues.append("CAI values out of valid range")
        
        if critical_issues:
            print(f"\n⚠️  Critical Issues Found:")
            for issue in critical_issues:
                print(f"  - {issue}")
        else:
            print(f"\n✅ No critical issues found")
        
        # Performance summary
        if 'table1' in all_results:
            values = all_results['table1']['values']
            print(f"\n📊 Performance Summary:")
            print(f"  - Best result: {values['min_value']:.4f} kcal/mol")
            print(f"  - Worst result: {values['max_value']:.4f} kcal/mol")
            print(f"  - Average: {values['mean_value']:.4f} ± {values['std_value']:.4f} kcal/mol")
        
        return all_results


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        table_dir = sys.argv[1]
    else:
        table_dir = "paper_experiment_results/tables"
    
    validator = PaperTableValidator()
    results = validator.validate_all_tables(table_dir)
    
    # Save validation results
    output_path = Path(table_dir) / "validation_results.json"
    with open(output_path, 'w') as f:
        # Convert numpy types to native Python types for JSON serialization
        def convert_types(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_types(v) for v in obj]
            return obj
        
        json.dump(convert_types(results), f, indent=2)
    
    print(f"\n💾 Validation results saved to: {output_path}")
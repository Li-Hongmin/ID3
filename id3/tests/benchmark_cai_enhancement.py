#!/usr/bin/env python3
"""
Unified benchmark for CAI-enhanced discretization methods.

Compares methods over 100 iterations on two protein sequences (shortest and longest)
and reports:
- Time per iteration
- CAI of discrete sequence
- Probability preservation score (sum log pi over chosen codons)
- Repetition rate (unique sequences / total)

Methods benchmarked (configurable):
- binary_search (existing BinarySearchCAIOptimizer via CAIEnhancementOperator)
- sado (existing SADOOptimizer via CAIEnhancementOperator)

This script does NOT change core code; it only evaluates existing components.
"""

import argparse
import json
import os
import time
from pathlib import Path
from typing import List, Dict, Tuple

import numpy as np
import torch

from id3.experiments.utils.data_loader import ProteinDataLoader
from id3.utils.constants import (
    amino_acid_token_map,
    amino_acid_to_codon_matrix,
    codon_to_rna_matrix,
)
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.cai.unified_enhancer import UnifiedCAIEnhancer
from id3.cai.validator import compute_cai_from_sequence
from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import codons as GLOBAL_CODONS


def list_available_protein_ids(data_dir: Path) -> List[str]:
    ids: List[str] = []
    for p in data_dir.iterdir():
        if p.is_file():
            name = p.name
            # Files like P00004.fasta.txt, P0DTC2.fasta, etc.
            if name.endswith('.fasta'):
                ids.append(name.split('.fasta')[0])
            elif name.endswith('.fasta.txt'):
                ids.append(name.split('.fasta.txt')[0])
            elif name.endswith('.txt') and len(name.split('.txt')[0]) > 0:
                # Avoid generic files; keep format similar to above
                base = name.split('.txt')[0]
                if len(base) >= 5 and base[0].isalpha():
                    ids.append(base)
    return sorted(list(set(ids)))


def select_shortest_and_longest(loader: ProteinDataLoader) -> Tuple[str, str]:
    data_dir = loader.data_dir
    protein_ids = list_available_protein_ids(data_dir)
    if not protein_ids:
        # fallback to some defaults defined in loader/test
        candidates = ["O15263", "P0DTC2", "P0CG48", "P01308"]
        lengths = [(pid, len(loader.load_protein_sequence(pid))) for pid in candidates]
    else:
        lengths = []
        for pid in protein_ids:
            try:
                seq = loader.load_protein_sequence(pid)
                lengths.append((pid, len(seq)))
            except Exception:
                continue
    if not lengths:
        raise RuntimeError("No protein sequences found in data/")
    lengths.sort(key=lambda x: x[1])
    shortest = lengths[0][0]
    longest = lengths[-1][0]
    return shortest, longest


def build_constraint_tensors(amino_acid_sequence: str, device: torch.device) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    max_codons = 6
    num_positions = len(amino_acid_sequence)

    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=device)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=device)
    codon_encodings = torch.zeros(num_positions, max_codons, 3, 4, device=device)

    for pos, aa in enumerate(amino_acid_sequence):
        aa_idx = amino_acid_token_map[aa]
        valid = amino_acid_to_codon_matrix[aa_idx].nonzero(as_tuple=True)[0]
        n = min(len(valid), max_codons)
        if n > 0:
            valid_codon_mask[pos, :n] = True
            codon_indices[pos, :n] = valid[:n]
            for i in range(n):
                codon_encodings[pos, i] = codon_to_rna_matrix[valid[i]]
    return valid_codon_mask, codon_indices, codon_encodings


def dirichlet_pi(num_positions: int, valid_codon_mask: torch.Tensor, seed: int = 42) -> torch.Tensor:
    """Create a skewed base π distribution over valid codons per position."""
    rng = np.random.default_rng(seed)
    max_codons = valid_codon_mask.shape[1]
    pi = torch.zeros(num_positions, max_codons)
    for pos in range(num_positions):
        valid = valid_codon_mask[pos].cpu().numpy()
        k = int(valid.sum())
        if k == 0:
            continue
        # Make one heavier component
        alpha = np.ones(k, dtype=np.float64)
        alpha[0] = 5.0
        vec = rng.dirichlet(alpha)
        # Fill into slots
        j = 0
        for slot in range(max_codons):
            if valid[slot]:
                pi[pos, slot] = float(vec[j])
                j += 1
    return pi


def perturb_pi(pi: torch.Tensor, valid_codon_mask: torch.Tensor, noise: float = 0.02, rng: np.random.Generator = None) -> torch.Tensor:
    """Apply small logit noise to π and renormalize over valid slots."""
    if rng is None:
        rng = np.random.default_rng()
    out = pi.clone()
    max_codons = pi.shape[1]
    for pos in range(pi.shape[0]):
        valid = valid_codon_mask[pos]
        if not bool(valid.any()):
            continue
        # Work in logits
        logits = torch.log(out[pos, valid] + 1e-12)
        noise_vec = torch.tensor(rng.normal(0.0, noise, size=(int(valid.sum()),)), dtype=logits.dtype)
        logits = logits + noise_vec
        probs = torch.softmax(logits, dim=-1)
        # write back only valid slots
        out[pos, valid] = probs
        # zero out invalid to keep shape clear
        out[pos, ~valid] = 0.0
    return out


def discrete_to_sequence(discrete_dist: torch.Tensor, codon_indices: torch.Tensor) -> str:
    """Convert one-hot codon distribution [positions, max_codons] into RNA sequence string."""
    if discrete_dist.dim() != 2:
        raise ValueError("discrete_dist must be [positions, max_codons]")
    nucs = ['A', 'C', 'G', 'U']
    parts: List[str] = []
    for pos in range(discrete_dist.shape[0]):
        slot = int(torch.argmax(discrete_dist[pos]).item())
        global_idx = int(codon_indices[pos, slot].item())
        codon_rna = codon_to_rna_matrix[global_idx]  # [3,4]
        for i in range(3):
            nt_idx = int(torch.argmax(codon_rna[i]).item())
            parts.append(nucs[nt_idx])
    return ''.join(parts)


def probs_to_onehot_discrete(probs: torch.Tensor, valid_mask: torch.Tensor) -> torch.Tensor:
    """Convert [positions, max_codons] probs to one-hot discrete over valid slots only."""
    assert probs.dim() == 2
    onehot = torch.zeros_like(probs)
    for pos in range(probs.shape[0]):
        vm = valid_mask[pos]
        if not bool(vm.any()):
            continue
        idx = torch.argmax(probs[pos, vm]).item()
        valid_slots = torch.nonzero(vm).flatten()
        slot = int(valid_slots[idx].item())
        onehot[pos, slot] = 1.0
    return onehot


def log_prob_of_choice(pi: torch.Tensor, discrete_dist: torch.Tensor) -> float:
    """Sum log π over the chosen codons per position (average per position)."""
    total = 0.0
    count = 0
    for pos in range(pi.shape[0]):
        idx = int(torch.argmax(discrete_dist[pos]).item())
        p = float(pi[pos, idx].item())
        if p <= 0:
            total += -30.0  # large negative
        else:
            total += np.log(p)
        count += 1
    return float(np.exp(total / max(1, count)))  # geometric mean as a compact score


class WarmStartSearcher:
    """External warm-start discrete search using UnifiedCAIEnhancer interpolation.

    Does not modify core optimizer code. Maintains last_alpha and searches
    over a small local range to find the minimal gamma achieving target CAI.
    """

    def __init__(self, device: torch.device, species: str = 'ecoli_bl21de3'):
        self.device = device
        self.enhancer = UnifiedCAIEnhancer(device=str(device), species=species, enable_cai=True)
        self.last_alpha = None

    def search(self,
               pi_probs: torch.Tensor,
               amino_seq: str,
               valid_mask: torch.Tensor,
               codon_indices: torch.Tensor,
               target_cai: float) -> Tuple[torch.Tensor, float, float]:
        # Define candidate grid
        if self.last_alpha is None:
            grid = np.linspace(0.0, 1.0, 21)
        else:
            span = 0.15
            low = max(0.0, self.last_alpha - span)
            high = min(1.0, self.last_alpha + span)
            grid = np.linspace(low, high, 9)

        best_alpha = 1.0
        best_probs = None
        best_cai = 0.0

        # Coarse scan and pick smallest alpha that meets target
        for alpha in grid:
            probs = self.enhancer.apply_cai_to_discrete_selection(
                pi_probs, amino_seq, valid_mask, codon_indices, cai_factor=float(alpha)
            )
            onehot = probs_to_onehot_discrete(probs, valid_mask)
            rna = discrete_to_sequence(onehot, codon_indices)
            cai_val = compute_cai_from_sequence(rna)
            if cai_val >= target_cai and alpha <= best_alpha + 1e-12:
                best_alpha = float(alpha)
                best_probs = probs
                best_cai = float(cai_val)
                break

        # If none met target, fallback to max alpha on the grid
        if best_probs is None:
            alpha = float(grid[-1])
            best_probs = self.enhancer.apply_cai_to_discrete_selection(
                pi_probs, amino_seq, valid_mask, codon_indices, cai_factor=alpha
            )
            onehot = probs_to_onehot_discrete(best_probs, valid_mask)
            rna = discrete_to_sequence(onehot, codon_indices)
            best_cai = float(compute_cai_from_sequence(rna))
            best_alpha = alpha

        self.last_alpha = best_alpha
        onehot = probs_to_onehot_discrete(best_probs, valid_mask)
        return onehot, best_alpha, best_cai


def apply_tie_break_jitter(pi: torch.Tensor,
                           onehot: torch.Tensor,
                           valid_mask: torch.Tensor,
                           epsilon: float = 0.02,
                           flip_ratio: float = 0.05,
                           rng: np.random.Generator = None) -> torch.Tensor:
    """Break ties/near-ties by occasionally flipping to the 2nd-best π candidate.

    This is a lightweight post-processing to simulate tie-breaking noise without
    modifying core optimizers. Only flips positions where top-2 π are close.
    """
    if rng is None:
        rng = np.random.default_rng()
    out = onehot.clone()
    max_codons = pi.shape[1]
    pos_indices = list(range(pi.shape[0]))
    rng.shuffle(pos_indices)
    flips_allowed = int(len(pos_indices) * flip_ratio)
    flips_done = 0
    for pos in pos_indices:
        if flips_done >= flips_allowed:
            break
        vm = valid_mask[pos]
        if not bool(vm.any()):
            continue
        # gather valid probs
        valid_probs = pi[pos, vm]
        if valid_probs.numel() < 2:
            continue
        # top2
        top2_vals, top2_idx = torch.topk(valid_probs, k=min(2, valid_probs.numel()))
        if top2_vals.numel() < 2:
            continue
        if float(top2_vals[0] - top2_vals[1]) <= epsilon:
            # flip with 50% chance
            if rng.random() < 0.5:
                valid_slots = torch.nonzero(vm).flatten()
                # current argmax slot in out
                cur_slot = int(torch.argmax(out[pos]).item())
                # target second-best slot index in full space
                second_local = int(top2_idx[1].item())
                second_slot = int(valid_slots[second_local].item())
                if second_slot != cur_slot:
                    out[pos] = torch.zeros(max_codons, dtype=onehot.dtype)
                    out[pos, second_slot] = 1.0
                    flips_done += 1
    return out


def _aa_to_codons(aa: str) -> List[str]:
    from id3.utils.constants import amino_acids_to_codons as AA2C
    return AA2C.get(aa, [])


class SadoFallbackHelper:
    """Helper to run SADO as a duplicate fallback starting from current sequence.

    Translates between our slot order and SADO's local codon order per position.
    """

    def __init__(self, amino_seq: str, valid_mask: torch.Tensor, codon_indices: torch.Tensor, device: torch.device):
        self.amino_seq = amino_seq
        self.valid_mask = valid_mask
        self.codon_indices = codon_indices
        self.device = device
        # pos -> list(local_idx -> slot)
        self.pos_maps: List[List[int]] = []
        for pos, aa in enumerate(amino_seq):
            mapping: List[int] = []
            valid_slots = torch.nonzero(valid_mask[pos]).flatten().tolist()
            for codon_str in _aa_to_codons(aa):
                # map codon string to global index in GLOBAL_CODONS
                try:
                    global_idx = GLOBAL_CODONS.index(codon_str)
                except ValueError:
                    mapping.append(-1)
                    continue
                slot_found = -1
                for j in valid_slots:
                    if int(codon_indices[pos, j].item()) == global_idx:
                        slot_found = int(j)
                        break
                mapping.append(slot_found)
            self.pos_maps.append(mapping)
        # SADO instance
        self.sado = SADOOptimizer(species='ecoli_bl21de3', device=device, amino_acid_sequence=amino_seq)

    def pi_to_sado(self, pi: torch.Tensor) -> torch.Tensor:
        num_positions = pi.shape[0]
        max_local = max(len(m) for m in self.pos_maps) if self.pos_maps else 1
        out = torch.zeros(num_positions, max_local, dtype=pi.dtype, device=self.device)
        for pos in range(num_positions):
            mapping = self.pos_maps[pos]
            for local_idx, slot in enumerate(mapping):
                if 0 <= slot < pi.shape[1]:
                    out[pos, local_idx] = pi[pos, slot]
        return out

    def onehot_sado_to_ours(self, onehot_sado: torch.Tensor) -> torch.Tensor:
        num_positions = onehot_sado.shape[0]
        max_codons = self.codon_indices.shape[1]
        out = torch.zeros(num_positions, max_codons, dtype=onehot_sado.dtype, device=self.device)
        for pos in range(num_positions):
            local_idx = int(torch.argmax(onehot_sado[pos]).item())
            mapping = self.pos_maps[pos]
            slot = mapping[local_idx] if local_idx < len(mapping) else -1
            if slot is not None and slot >= 0:
                out[pos, slot] = 1.0
        return out

    def last_indices_from_onehot(self, onehot_ours: torch.Tensor) -> np.ndarray:
        """Map our discrete choice to SADO's sorted local indices per position.

        SADO codon_choices are sorted by CAI weight and include 'local_index'.
        We must set last_indices[pos] to that sorted local_index for the chosen codon.
        """
        num_positions = onehot_ours.shape[0]
        last = np.zeros(num_positions, dtype=np.int32)
        for pos in range(num_positions):
            slot = int(torch.argmax(onehot_ours[pos]).item())
            global_idx = int(self.codon_indices[pos, slot].item())
            # Convert global index to RNA codon string in our constants order
            codon_str = GLOBAL_CODONS[global_idx]
            # Build map from codon to sorted local_index for this position
            sorted_list = getattr(self.sado, 'codon_choices', [])[pos] if pos < len(getattr(self.sado, 'codon_choices', [])) else []
            local_idx_sorted = 0
            for ch in sorted_list:
                if ch.get('codon') == codon_str:
                    local_idx_sorted = ch.get('local_index', 0)
                    break
            last[pos] = int(local_idx_sorted)
        return last


def run_benchmark_for_sequence(
    amino_acid_sequence: str,
    methods: List[str],
    iterations: int = 100,
    target_cai: float = 0.8,
    seed: int = 42,
    device: str = 'cpu',
    enable_warm_start: bool = False,
    enable_jitter: bool = False,
    enable_sado_fallback: bool = False,
) -> Dict:
    torch_device = torch.device(device)
    valid_mask, codon_indices, _ = build_constraint_tensors(amino_acid_sequence, torch_device)
    base_pi = dirichlet_pi(len(amino_acid_sequence), valid_mask, seed=seed)
    rng = np.random.default_rng(seed + 1)
    warm = WarmStartSearcher(torch_device) if enable_warm_start else None
    sado_helper = SadoFallbackHelper(amino_acid_sequence, valid_mask, codon_indices, torch_device) if enable_sado_fallback else None

    results = {}
    for method in methods:
        op = CAIEnhancementOperator(method=method, species='ecoli_bl21de3', device=torch_device, amino_acid_sequence=amino_acid_sequence)
        times: List[float] = []
        cais: List[float] = []
        probs: List[float] = []
        seen_hashes = set()
        unique_count = 0

        # Reset method state if available (e.g., SADO)
        if hasattr(op, 'reset'):
            try:
                op.reset()
            except Exception:
                pass

        pi = base_pi.clone()
        for it in range(iterations):
            # Small perturbation to simulate iterative changes
            pi = perturb_pi(pi, valid_mask, noise=0.02, rng=rng)
            t0 = time.perf_counter()
            if enable_warm_start and method == 'binary_search':
                # External warm-start path
                discrete_dist, alpha_sel, cai_ws = warm.search(
                    pi, amino_acid_sequence, valid_mask, codon_indices, target_cai
                )
                dt = time.perf_counter() - t0
            else:
                discrete_dist, meta = op.apply_cai_enhancement(
                    pi, amino_acid_sequence, valid_mask, codon_indices, target_cai
                )
                dt = time.perf_counter() - t0

            # Optional jitter to break ties
            if enable_jitter:
                discrete_dist = apply_tie_break_jitter(
                    pi, discrete_dist, valid_mask, epsilon=0.02, flip_ratio=0.05, rng=rng
                )
            times.append(dt)

            # Convert to RNA and evaluate CAI
            rna_seq = discrete_to_sequence(discrete_dist, codon_indices)
            cai_val = compute_cai_from_sequence(rna_seq)
            # Ensure CAI after jitter is not below target (small warm-fix)
            if enable_jitter and cai_val < target_cai:
                tfix0 = time.perf_counter()
                local_warm = WarmStartSearcher(torch_device)
                onehot_fix, alpha_fix, cai_fix = local_warm.search(
                    pi, amino_acid_sequence, valid_mask, codon_indices, target_cai
                )
                discrete_dist = onehot_fix
                rna_seq = discrete_to_sequence(discrete_dist, codon_indices)
                cai_val = compute_cai_from_sequence(rna_seq)
                # account time
                times[-1] += (time.perf_counter() - tfix0)
            cais.append(float(cai_val))

            # Probability preservation score (geometric mean of π over chosen codons)
            prob_score = log_prob_of_choice(pi, discrete_dist)
            probs.append(float(prob_score))

            # Repetition tracking + optional SADO fallback
            h = hash(rna_seq)
            if h in seen_hashes and enable_sado_fallback and method == 'binary_search':
                t1 = time.perf_counter()
                # Initialize SADO from current discrete choice
                last_idx = sado_helper.last_indices_from_onehot(discrete_dist)
                sado_helper.sado.last_indices = last_idx
                sado_helper.sado.history_hashes.add(sado_helper.sado._hash_sequence(last_idx))
                pi_sado = sado_helper.pi_to_sado(pi)
                # Run SADO to get a different sequence with CAI target if possible
                dist_sado, meta_sado = sado_helper.sado.optimize(
                    pi_accessibility=pi_sado,
                    target_cai=target_cai,
                    amino_acid_sequence=amino_acid_sequence,
                    use_binary_search=False,
                    use_difference_driven=True
                )
                discrete_dist = sado_helper.onehot_sado_to_ours(dist_sado)
                rna_seq = discrete_to_sequence(discrete_dist, codon_indices)
                cai_val = compute_cai_from_sequence(rna_seq)
                # If CAI insufficient, do a minimal warm-start push to reach target
                if cai_val < target_cai:
                    local_warm = WarmStartSearcher(self_device := self_device if False else torch_device)
                    onehot_fix, alpha_fix, cai_fix = local_warm.search(
                        pi, amino_acid_sequence, valid_mask, codon_indices, target_cai
                    )
                    discrete_dist = onehot_fix
                    rna_seq = discrete_to_sequence(discrete_dist, codon_indices)
                    cai_val = compute_cai_from_sequence(rna_seq)
                cais[-1] = float(cai_val)
                probs[-1] = log_prob_of_choice(pi, discrete_dist)
                dt += (time.perf_counter() - t1)
                h = hash(rna_seq)
                times[-1] = dt
            if h not in seen_hashes:
                unique_count += 1
                seen_hashes.add(h)

        results[method] = {
            'iterations': iterations,
            'avg_time_sec': float(np.mean(times)),
            'p50_time_sec': float(np.median(times)),
            'total_time_sec': float(np.sum(times)),
            'cai_mean': float(np.mean(cais)),
            'cai_min': float(np.min(cais)),
            'cai_max': float(np.max(cais)),
            'prob_score_mean': float(np.mean(probs)),
            'unique_sequences': int(unique_count),
            'repetition_rate': float(1.0 - unique_count / iterations),
        }

    return results


def main():
    parser = argparse.ArgumentParser(description='Benchmark CAI-enhanced discretization methods')
    parser.add_argument('--methods', type=str, default='binary_search,sado', help='Comma-separated methods: binary_search,sado')
    parser.add_argument('--iterations', type=int, default=100)
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--warm_start', action='store_true', help='Enable external gamma warm-start for binary_search')
    parser.add_argument('--jitter', action='store_true', help='Enable tie-breaking jitter on discrete output')
    parser.add_argument('--target_cai', type=float, default=0.8)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', type=str, default=None, help='Output JSON path')
    parser.add_argument('--sado_fallback', action='store_true', help='Use SADO to avoid duplicates when detected')
    args = parser.parse_args()

    methods = [m.strip() for m in args.methods.split(',') if m.strip()]
    loader = ProteinDataLoader()
    shortest, longest = select_shortest_and_longest(loader)

    report = {
        'config': {
            'methods': methods,
            'iterations': args.iterations,
            'device': args.device,
            'target_cai': args.target_cai,
            'seed': args.seed,
            'sequences': {'shortest': shortest, 'longest': longest},
        },
        'results': {}
    }

    for label, pid in [('shortest', shortest), ('longest', longest)]:
        info = loader.get_protein_info(pid)
        print(f"Benchmarking {label} protein {pid} (len={info['length']})...")
        seq = info['sequence']
        res = run_benchmark_for_sequence(
            seq, methods, iterations=args.iterations, target_cai=args.target_cai,
            seed=args.seed, device=args.device,
            enable_warm_start=args.warm_start,
            enable_jitter=args.jitter,
            enable_sado_fallback=args.sado_fallback,
        )
        report['results'][label] = res

    # Output
    ts = time.strftime('%Y%m%d_%H%M%S')
    out_dir = Path('results') / 'benchmarks'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = Path(args.output) if args.output else (out_dir / f'benchmark_cai_enhancement_{ts}.json')
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"Saved benchmark report to {out_path}")


if __name__ == '__main__':
    main()

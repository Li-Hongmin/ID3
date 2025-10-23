"""








"""

import torch
import torch.nn.functional as F
import numpy as np
from typing import Dict, List, Optional, Tuple
import json
import os


from ..utils.constants import amino_acids_to_codons as rna_amino_acids_to_codons, NUCLEOTIDES as RNA_NUCLEOTIDES, NUCLEOTIDE_MAP


amino_acids_to_codons = {}
for aa, rna_codons in rna_amino_acids_to_codons.items():
    amino_acids_to_codons[aa] = [codon.replace('U', 'T') for codon in rna_codons]


NUCLEOTIDES = ['A', 'C', 'G', 'T']
NUCLEOTIDE_MAP_DNA = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}


amino_acid_token_map = {aa: idx for idx, aa in enumerate(rna_amino_acids_to_codons.keys())}


class CAIModule:
    """
    CAI (Codon Adaptation Index) calculation module.

    Provides methods for computing CAI scores and statistics for DNA/RNA sequences.
    """
    
    def __init__(self,
                 reference_species: str = 'ecoli_bl21de3',
                 reference_sequences: Optional[List[str]] = None,
                 device: torch.device = None):
        """
        Initialize CAI module.

        Args:
            reference_species: Reference organism species
            reference_sequences: Optional reference sequences for custom weight calculation
            device: Computation device (CPU or CUDA)
        """
        self.reference_species = reference_species
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.reference_sequences = reference_sequences
        

        supported_species = ['ecoli_bl21de3', 'ecoli_pET21a']
        if reference_species not in supported_species:
            print(f"⚠️  Warning: Unsupported reference species {reference_species}, using default ecoli_bl21de3")
            self.reference_species = 'ecoli_bl21de3'


        self.codon_weights = self._load_codon_weights()
        self.rare_codons = self._get_rare_codons()
        self.preferred_codons = self._get_preferred_codons()


        self.codon_weights_tensor = self._create_codon_weights_tensor()

        print(f"🧬 CAI module initialization complete")
        print(f"📊 Reference system: {reference_species}")
        print(f"⚠️  Rare codons: {len(self.rare_codons)} codons")
        print(f"✅ Preferred codons: {len(self.preferred_codons)} codons")
    
    def _load_codon_weights(self) -> Dict[str, float]:
        """
        Load codon weights from reference sequences or precomputed data.

        Returns:
            Dictionary mapping codons to their CAI weights
        """
        # Use custom reference sequences if provided
        if self.reference_sequences:
            print(f"📊 Computing RSCU from provided reference sequences (count: {len(self.reference_sequences)})")
            return self._compute_weights_from_sequences(self.reference_sequences)

        # Load precomputed weights
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        weights_file = os.path.join(project_root, 'data/codon_references/ecoli_bl21de3_wi_weights_comparison.json')

        print(f"📊 Loading precomputed CAI weights: {weights_file}")
        return self._load_precomputed_weights(weights_file)
    
    def _load_precomputed_weights(self, weights_file: str) -> Dict[str, float]:
        """
        Load precomputed CAI weights from JSON file.

        Args:
            weights_file: Path to weights file

        Returns:
            Dictionary of codon weights
        """
        with open(weights_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # Extract wi_table
        codon_weights = data['wi_table'].copy()

        print(f"📈 Loading complete: {len(codon_weights)} codon weights")
        print(f"📊 Weight range: {min(codon_weights.values()):.4f} - {max(codon_weights.values()):.4f}")

        return codon_weights
    
    def _compute_weights_from_sequences(self, sequences: List[str]) -> Dict[str, float]:
        """
        Compute codon weights from reference sequences.

        Calculates RSCU (Relative Synonymous Codon Usage) and derives wi weights.

        Returns:
            Dictionary of codon weights
        """
        # Count codon occurrences
        codon_counts = {}

        for seq in sequences:
            # Process complete codons only
            seq_len = len(seq) - (len(seq) % 3)
            for i in range(0, seq_len, 3):
                codon = seq[i:i+3].upper()
                if len(codon) == 3 and all(c in 'ATCG' for c in codon):
                    codon_counts[codon] = codon_counts.get(codon, 0) + 1

        # Calculate RSCU values
        rscu_values = {}

        for aa, rna_codons in amino_acids_to_codons.items():
            if aa == 'X':
                continue

            # Get codons for this amino acid
            aa_codons = []
            aa_counts = []

            for rna_codon in rna_codons:
                dna_codon = rna_codon.replace('U', 'T')
                if dna_codon in codon_counts:
                    aa_codons.append(dna_codon)
                    aa_counts.append(codon_counts[dna_codon])

            if aa_counts:
                total_count = sum(aa_counts)
                num_codons = len(aa_counts)

                # RSCU = (observed count / expected count) * number of synonymous codons
                # = (count / (total/n)) * n = count * n / total
                for codon, count in zip(aa_codons, aa_counts):
                    rscu_values[codon] = count * num_codons / total_count

        # Convert RSCU to CAI weights
        codon_weights = {}

        for aa, rna_codons in amino_acids_to_codons.items():
            if aa == 'X':
                continue

            # Get RSCU values for this amino acid
            aa_rscus = []
            aa_codons = []

            for rna_codon in rna_codons:
                dna_codon = rna_codon.replace('U', 'T')
                if dna_codon in rscu_values:
                    aa_rscus.append(rscu_values[dna_codon])
                    aa_codons.append(dna_codon)

            if aa_rscus:
                max_rscu = max(aa_rscus)
                for codon, rscu in zip(aa_codons, aa_rscus):
                    codon_weights[codon] = rscu / max_rscu if max_rscu > 0 else 1.0

        print(f"📈 Calculation complete: {len(codon_weights)} codon weights")

        # Store RSCU values for reference
        self.rscu_values = rscu_values

        return codon_weights
    
    def _load_fallback_weights(self) -> Dict[str, float]:
        """
        Load fallback E. coli K12 codon frequencies when data files are unavailable.

        Returns:
            Dictionary of codon weights based on K12 frequencies
        """
        # E. coli K12 codon usage frequencies (per 1000 codons)
        ecoli_codon_frequencies = {
            'TTT': 19.7, 'TTC': 15.0, 'TTA': 15.2, 'TTG': 11.9,
            'TCT': 5.7, 'TCC': 5.5, 'TCA': 7.8, 'TCG': 8.0,
            'TAT': 16.8, 'TAC': 14.6, 'TAA': 1.8, 'TAG': 0.0,
            'TGT': 5.9, 'TGC': 8.0, 'TGA': 1.0, 'TGG': 10.7,
            'CTT': 11.9, 'CTC': 10.5, 'CTA': 5.3, 'CTG': 46.9,
            'CCT': 8.4, 'CCC': 6.4, 'CCA': 6.6, 'CCG': 26.7,
            'CAT': 15.8, 'CAC': 13.1, 'CAA': 12.1, 'CAG': 27.7,
            'CGT': 21.1, 'CGC': 26.0, 'CGA': 4.3, 'CGG': 4.1,
            'ATT': 30.5, 'ATC': 18.2, 'ATA': 3.7, 'ATG': 24.8,
            'ACT': 8.0, 'ACC': 22.8, 'ACA': 6.4, 'ACG': 11.5,
            'AAT': 21.9, 'AAC': 24.4, 'AAA': 33.2, 'AAG': 12.1,
            'AGT': 7.2, 'AGC': 16.6, 'AGA': 1.4, 'AGG': 1.6,
            'GTT': 16.8, 'GTC': 11.7, 'GTA': 11.5, 'GTG': 26.4,
            'GCT': 10.7, 'GCC': 31.6, 'GCA': 21.1, 'GCG': 38.5,
            'GAT': 37.9, 'GAC': 20.5, 'GAA': 43.7, 'GAG': 18.4,
            'GGT': 16.8, 'GGC': 31.6, 'GGA': 11.5, 'GGG': 12.1
        }

        # Normalize frequencies to CAI weights
        codon_weights = {}

        # Process each amino acid family
        for aa, rna_codons in amino_acids_to_codons.items():
            # Get frequencies for synonymous codons
            aa_frequencies = []
            aa_codons = []

            for rna_codon in rna_codons:
                # Convert RNA to DNA codon
                dna_codon = rna_codon.replace('U', 'T')
                if dna_codon in ecoli_codon_frequencies:
                    aa_frequencies.append(ecoli_codon_frequencies[dna_codon])
                    aa_codons.append(dna_codon)

            if aa_frequencies:
                max_freq = max(aa_frequencies)
                # Normalize: wi = freq / max_freq
                for codon, freq in zip(aa_codons, aa_frequencies):
                    codon_weights[codon] = freq / max_freq

        return codon_weights
    
    def _get_rare_codons(self) -> List[str]:
        """
        Get list of rare codons for the reference species.

        Rare codons are defined as those with CAI weight < 0.25.

        Returns:
            List of rare codon sequences
        """
        if self.reference_species == 'ecoli_bl21de3':
            # Try to load BL21(DE3) specific rare codon data
            try:
                project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                data_file = os.path.join(project_root, "data/codon_references/build_process/intermediate_data/bl21de3_codon_frequencies.json")
                with open(data_file, 'r', encoding='utf-8') as f:
                    bl21_data = json.load(f)

                rare_codons = list(bl21_data['rare_codons'].keys())
                # Filter out stop codons and non-synonymous codons
                stop_codons = {'TAA', 'TAG', 'TGA'}
                non_synonymous_codons = {'ATG', 'TGG'}
                rare_codons = [codon for codon in rare_codons
                              if codon not in stop_codons and codon not in non_synonymous_codons]
                return rare_codons

            except FileNotFoundError:
                print(f"⚠️  Warning: BL21(DE3) rare codon data not found, using K12 data")

        # Fallback: identify rare codons by weight threshold
        # Rare codons are those with weight < 0.25
        stop_codons = {'TAA', 'TAG', 'TGA'}
        non_synonymous_codons = {'ATG', 'TGG'}
        rare_codons = [codon for codon, weight in self.codon_weights.items()
                      if weight < 0.25 and codon not in stop_codons and codon not in non_synonymous_codons]
        return rare_codons
    
    def _get_preferred_codons(self) -> List[str]:
        """
        Get list of preferred codons for the reference species.

        Preferred codons are defined as those with CAI weight > 0.5.

        Returns:
            List of preferred codon sequences
        """
        if self.reference_species == 'ecoli_bl21de3':
            # Try to load BL21(DE3) specific preferred codon data
            try:
                project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                data_file = os.path.join(project_root, "data/codon_references/build_process/intermediate_data/bl21de3_codon_frequencies.json")
                with open(data_file, 'r', encoding='utf-8') as f:
                    bl21_data = json.load(f)

                preferred_codons = list(bl21_data['preferred_codons'].keys())
                # Filter out stop codons and non-synonymous codons
                stop_codons = {'TAA', 'TAG', 'TGA'}
                non_synonymous_codons = {'ATG', 'TGG'}
                preferred_codons = [codon for codon in preferred_codons
                                   if codon not in stop_codons and codon not in non_synonymous_codons]
                return preferred_codons

            except FileNotFoundError:
                print(f"⚠️  Warning: BL21(DE3) preferred codon data not found, using K12 data")

        # Fallback: identify preferred codons by weight threshold
        # Preferred codons are those with weight > 0.5
        stop_codons = {'TAA', 'TAG', 'TGA'}
        non_synonymous_codons = {'ATG', 'TGG'}
        preferred_codons = [codon for codon, weight in self.codon_weights.items()
                           if weight > 0.5 and codon not in stop_codons and codon not in non_synonymous_codons]
        return preferred_codons
    
    def _create_codon_weights_tensor(self) -> torch.Tensor:
        """
        Create tensor of codon weights for efficient computation.

        Returns:
            Tensor of shape [64] containing CAI weights for all codons
        """
        # Generate all possible DNA codons in standard order
        nucleotides = [nt.replace('U', 'T') for nt in NUCLEOTIDES]
        all_codons = []
        weights = []

        for n1 in nucleotides:
            for n2 in nucleotides:
                for n3 in nucleotides:
                    codon = n1 + n2 + n3
                    all_codons.append(codon)
                    # Default weight 0.1 for unknown codons
                    weight = self.codon_weights.get(codon, 0.1)
                    weights.append(weight)

        weights_tensor = torch.tensor(weights, dtype=torch.float32, device=self.device)
        return weights_tensor
    
    def compute_cai_score(self, rna_sequence: str) -> float:
        """
        Compute CAI score for an RNA sequence.

        Args:
            rna_sequence: RNA sequence string

        Returns:
            CAI score (0.0 - 1.0)
        """
        # Split sequence into codons
        codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]

        # Handle empty sequence
        if not codons:
            return 0.0

        # Collect weights for synonymous codons only
        stop_codons = {'TAA', 'TAG', 'TGA'}
        non_synonymous_codons = {'ATG', 'TGG'}
        weights = []
        for codon in codons:
            if len(codon) == 3:
                # Convert RNA to DNA
                dna_codon = codon.replace('U', 'T')
                # Skip stop codons and non-synonymous codons
                if dna_codon not in stop_codons and dna_codon not in non_synonymous_codons:
                    weight = self.codon_weights.get(dna_codon, 0.1)
                    weights.append(weight)

        if not weights:
            return 0.0

        # CAI = (∏ w_i)^(1/L)
        weights_array = np.array(weights)
        geometric_mean = np.exp(np.mean(np.log(weights_array + 1e-10)))
        # Clip to valid range
        cai_score = np.clip(geometric_mean, 0.0, 1.0)
        return float(cai_score)
    
    def compute_cai_differentiable(self,
                                   codon_probabilities: torch.Tensor,
                                   sequence_length: Optional[int] = None) -> torch.Tensor:
        """
        Compute differentiable CAI score from codon probability distributions.

        Args:
            codon_probabilities: Codon probability tensor [seq_len//3, 64] or [batch, seq_len//3, 64]
            sequence_length: Optional sequence length for validation

        Returns:
            CAI score tensor (differentiable)
        """

        if codon_probabilities.dim() == 2:
            # [seq_len//3, 64] -> [1, seq_len//3, 64]
            codon_probabilities = codon_probabilities.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False
        
        batch_size, num_codons, vocab_size = codon_probabilities.shape

        # Ensure weights tensor is on same device
        from id3.utils.memory_utils import reduce_memory_usage
        if self.codon_weights_tensor.device != codon_probabilities.device:
            self.codon_weights_tensor = self.codon_weights_tensor.to(codon_probabilities.device)

        # Memory optimization
        if not codon_probabilities.requires_grad and codon_probabilities.dtype == torch.float64:
            codon_probabilities = codon_probabilities.to(torch.float32)

        # Create mask to exclude stop codons and non-synonymous codons
        device = codon_probabilities.device
        mask = torch.ones(64, dtype=torch.bool, device=device)

        # Build codon index mapping
        nucleotides = [nt.replace('U', 'T') for nt in NUCLEOTIDES]
        codon_to_idx = {}
        idx = 0
        for n1 in nucleotides:
            for n2 in nucleotides:
                for n3 in nucleotides:
                    codon_to_idx[n1+n2+n3] = idx
                    idx += 1

        # Exclude stop codons and non-synonymous codons
        excluded_codons = {'TAA', 'TAG', 'TGA', 'ATG', 'TGG'}
        for codon in excluded_codons:
            if codon in codon_to_idx:
                mask[codon_to_idx[codon]] = False

        # Apply mask to probabilities and weights
        # Mask out invalid codons (stop codons and non-synonymous)
        valid_probs = codon_probabilities * mask.float()  # [batch, num_codons, 64]
        valid_weights = self.codon_weights_tensor * mask.float()  # [64]

        # Check which positions have synonymous codons
        has_synonymous = torch.sum(valid_probs, dim=-1) > 1e-10  # [batch, num_codons]

        # Compute expected CAI weights per position
        # [batch, num_codons, 64] * [64] -> [batch, num_codons]
        weighted_probs = torch.sum(valid_probs * valid_weights, dim=-1)

        # Extract weights for positions with synonymous codons
        valid_weighted_probs = weighted_probs[has_synonymous]

        if valid_weighted_probs.numel() == 0:
            # All codons are non-synonymous or missing
            result = torch.ones(batch_size, device=device)
        else:
            # Compute geometric mean of CAI weights
            from id3.utils.numerical_stability import safe_geometric_mean
            cai_score = safe_geometric_mean(valid_weighted_probs)
            result = cai_score.unsqueeze(0).expand(batch_size)

        if squeeze_output:
            result = result.squeeze(0)

        return result
    
    def get_codon_penalty(self, codon_probabilities: torch.Tensor) -> torch.Tensor:
        """
        Compute penalty term for rare codon usage.

        Args:
            codon_probabilities: Codon probability tensor

        Returns:
            Average rare codon usage (penalty term)
        """
        # Create mask for rare codons
        rare_codon_mask = torch.zeros(64, device=self.device)
        nucleotides = [nt.replace('U', 'T') for nt in NUCLEOTIDES]

        codon_idx = 0
        for n1 in nucleotides:
            for n2 in nucleotides:
                for n3 in nucleotides:
                    codon = n1 + n2 + n3
                    if codon in self.rare_codons:
                        rare_codon_mask[codon_idx] = 1.0
                    codon_idx += 1

        # Compute rare codon usage
        if codon_probabilities.dim() == 2:
            rare_usage = torch.sum(codon_probabilities * rare_codon_mask, dim=-1)
        else:
            rare_usage = torch.sum(codon_probabilities * rare_codon_mask, dim=-1)

        # Return average rare codon usage as penalty
        return torch.mean(rare_usage)
    
    def integrate_cai_with_logits(self,
                                  logits: torch.Tensor,
                                  integration_method: str = 'additive') -> torch.Tensor:
        """
        Integrate CAI weights with model logits for biased sampling.

        Args:
            logits: Model output logits [batch, seq_len, vocab_size]
            integration_method: Integration method ('additive' or 'multiplicative')

        Returns:
            Enhanced logits with CAI bias
        """
        # Ensure weights on same device
        if self.codon_weights_tensor.device != logits.device:
            self.codon_weights_tensor = self.codon_weights_tensor.to(logits.device)

        if integration_method == 'additive':
            # Add log weights to logits (log-space addition)
            log_weights = torch.log(self.codon_weights_tensor + 1e-10)
            enhanced_logits = logits + log_weights
        elif integration_method == 'multiplicative':
            # Multiply logits by weights
            enhanced_logits = logits * self.codon_weights_tensor
        else:
            raise ValueError(f"Unsupported integration method: {integration_method}")

        return enhanced_logits
    
    def get_cai_statistics(self, rna_sequence: str) -> Dict:
        """
        Compute comprehensive CAI statistics for an RNA sequence.

        Args:
            rna_sequence: RNA sequence string

        Returns:
            Dictionary containing CAI score and codon usage statistics
        """
        codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]

        rare_codon_count = 0
        preferred_codon_count = 0
        total_codons = 0

        codon_weights = []

        for codon in codons:
            if len(codon) == 3:
                total_codons += 1

                # Check if rare codon
                if codon in self.rare_codons:
                    rare_codon_count += 1

                # Check if preferred codon
                if codon in self.preferred_codons:
                    preferred_codon_count += 1

                # Collect weight
                weight = self.codon_weights.get(codon, 0.1)
                codon_weights.append(weight)

        stats = {
            'cai_score': self.compute_cai_score(rna_sequence),
            'total_codons': total_codons,
            'rare_codon_count': rare_codon_count,
            'rare_codon_ratio': rare_codon_count / total_codons if total_codons > 0 else 0,
            'preferred_codon_count': preferred_codon_count,
            'preferred_codon_ratio': preferred_codon_count / total_codons if total_codons > 0 else 0,
            'mean_codon_weight': np.mean(codon_weights) if codon_weights else 0,
            'min_codon_weight': np.min(codon_weights) if codon_weights else 0,
            'max_codon_weight': np.max(codon_weights) if codon_weights else 0
        }

        return stats
    
    def save_codon_weights(self, filepath: str):
        """Save codon weights to JSON file."""
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(self.codon_weights, f, indent=2, ensure_ascii=False)
        print(f"💾 Saved codon weights to: {filepath}")

    
    def export_wi_table(self, filepath: str):
        """
        Export complete wi weight table for differentiable CAI computation.

        Arranged in standard order of 64 codons: AAA, AAC, AAG, AAT, ACA, ACC, ...
        """


        nucleotides = [nt.replace('U', 'T') for nt in NUCLEOTIDES]
        all_codons = []
        for n1 in nucleotides:
            for n2 in nucleotides:
                for n3 in nucleotides:
                    all_codons.append(n1 + n2 + n3)

        # Build wi weight table
        wi_table = {}
        wi_values = []

        stop_codons = {'TAA', 'TAG', 'TGA'}

        for codon in all_codons:
            if codon in stop_codons:
                # Stop codons have weight 0
                weight = 0.0
            else:
                # Use preloaded codon weights, default to 0.1 if not found
                weight = self.codon_weights.get(codon, 0.1)

            wi_table[codon] = weight
            wi_values.append(weight)
        
        # Create export data structure
        export_data = {
            'metadata': {
                'description': 'CAI codon weight table (wi)',
                'reference_species': self.reference_species,
                'calculation_method': 'wi = RSCU_i / max(RSCU)',
                'total_codons': len(all_codons),
                'codon_order': 'AAA, AAC, AAG, AAT, ACA, ACC, ...',
                'stop_codons_weight': 0.0,
                'default_weight': 0.1
            },
            'wi_table': wi_table,
            'wi_values_array': wi_values,
            'codon_order': all_codons,
            'rscu_values': getattr(self, 'rscu_values', {}),
            'statistics': {
                'min_weight': min(wi_values),
                'max_weight': max(wi_values),
                'mean_weight': sum(wi_values) / len(wi_values),
                'non_zero_weights': sum(1 for w in wi_values if w > 0)
            }
        }

        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2, ensure_ascii=False)

        print(f"💾 Exported wi weight table to: {filepath}")
        print(f"📊 Statistics: {len(all_codons)} codons total")
        print(f"   Weight range: [{min(wi_values):.4f}, {max(wi_values):.4f}]")
        print(f"   Mean weight: {sum(wi_values) / len(wi_values):.4f}")
        print(f"   Non-zero weights: {sum(1 for w in wi_values if w > 0)} / {len(wi_values)}")





        
        return export_data
    
    def load_custom_codon_weights(self, filepath: str):
        """Load custom codon weights from file."""
        with open(filepath, 'r', encoding='utf-8') as f:
            self.codon_weights = json.load(f)

        # Rebuild codon weights tensor
        self.codon_weights_tensor = self._create_codon_weights_tensor()
        print(f"📁 Loaded custom codon weights from: {filepath}")



def create_cai_module(reference_species: str = 'ecoli_bl21de3',
                      reference_sequences: Optional[List[str]] = None,
                      device: torch.device = None) -> CAIModule:
    """
    Factory function to create CAI module.

    Args:
        reference_species: Reference species for CAI calculation
        reference_sequences: Optional reference sequences
        device: Computation device (CPU or CUDA)

    Returns:
        Configured CAI module instance
    """
    return CAIModule(reference_species=reference_species, 
                     reference_sequences=reference_sequences, 
                     device=device)


def compute_sequence_cai(rna_sequence: str,
                        reference_species: str = 'ecoli_bl21de3') -> float:
    """
    Convenience function to compute CAI score for a single sequence.

    Args:
        rna_sequence: RNA sequence string
        reference_species: Reference species

    Returns:
        CAI score (0.0 - 1.0)
    """
    cai_module = CAIModule(reference_species=reference_species)
    return cai_module.compute_cai_score(rna_sequence)


if __name__ == "__main__":
    # Test CAI module
    print("🧪 Testing CAI module...")

    # Create CAI module
    cai_module = CAIModule()

    # Test sequences
    test_sequence_good = "ATGCTGGAAAAAGCGCGCGAC"
    test_sequence_bad = "ATGCTACAAAGAGGACCACCC"

    print(f"\n🧬 Test sequence 1 (preferred codons): {test_sequence_good}")
    stats1 = cai_module.get_cai_statistics(test_sequence_good)
    print(f"CAI score: {stats1['cai_score']:.4f}")
    print(f"Rare codon ratio: {stats1['rare_codon_ratio']:.2%}")

    print(f"\n🧬 Test sequence 2 (rare codons): {test_sequence_bad}")
    stats2 = cai_module.get_cai_statistics(test_sequence_bad)
    print(f"CAI score: {stats2['cai_score']:.4f}")
    print(f"Rare codon ratio: {stats2['rare_codon_ratio']:.2%}")

    # Test differentiable CAI calculation
    print(f"\n🔢 Testing differentiable CAI calculation...")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Create codon probability distribution
    num_codons = len(test_sequence_good) // 3
    # Get actual codons
    actual_codons = [test_sequence_good[i:i+3] for i in range(0, len(test_sequence_good), 3)]
    codon_probs = torch.zeros(num_codons, 64, device=device)

    # Set high probability for actual codons, low probability for others
    for i, codon in enumerate(actual_codons):
        if len(codon) == 3:
            # Simple codon index calculation (for demo purposes)
            codon_idx = 0
            try:
                # Simple hash function
                codon_idx = (ord(codon[0]) - ord('A')) * 16 + (ord(codon[1]) - ord('A')) * 4 + (ord(codon[2]) - ord('A'))
                codon_idx = max(0, min(63, codon_idx))
                codon_probs[i, codon_idx] = 0.8
                codon_probs[i, :] += 0.2 / 64
            except (IndexError, ValueError, TypeError) as e:
                # Fallback: uniform distribution
                print(f"Warning: Failed to calculate codon probabilities for codon {codon}: {e}")
                codon_probs[i, :] = 1.0 / 64

    codon_probs = F.softmax(codon_probs, dim=-1)

    cai_score_diff = cai_module.compute_cai_differentiable(codon_probs)
    print(f"Differentiable CAI score: {cai_score_diff.item():.4f}")

    print("✅ CAI module testing complete!")
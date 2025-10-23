import torch


NUCLEOTIDES = ['A', 'C', 'G', 'U']
NUCLEOTIDE_MAP = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}


base_token_map = NUCLEOTIDE_MAP


amino_acids_to_codons = {
    "F": ["UUU", "UUC"],
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "I": ["AUU", "AUC", "AUA"],
    "M": ["AUG"],
    "V": ["GUU", "GUC", "GUA", "GUG"],
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "P": ["CCU", "CCC", "CCA", "CCG"],
    "T": ["ACU", "ACC", "ACA", "ACG"],
    "A": ["GCU", "GCC", "GCA", "GCG"],
    "Y": ["UAU", "UAC"],
    "H": ["CAU", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAU", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAU", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["UGU", "UGC"],
    "W": ["UGG"],
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGU", "GGC", "GGA", "GGG"]
}


codons = [codon for codon_list in amino_acids_to_codons.values() for codon in codon_list]
amino_acid_token_map = {aa: idx for idx, aa in enumerate(amino_acids_to_codons)}
num_codons = len(codons)
num_amino_acids = len(amino_acids_to_codons)


codon_to_rna_matrix = torch.zeros(num_codons, 3, 4)
for i, codon in enumerate(codons):
    for j, base in enumerate(codon):
        codon_to_rna_matrix[i, j, base_token_map[base]] = 1.0


amino_acid_to_codon_matrix = torch.zeros(num_amino_acids, num_codons)
for aa, codon_list in amino_acids_to_codons.items():
    amino_acid_idx = amino_acid_token_map[aa]
    codon_indices = [codons.index(codon) for codon in codon_list]
    amino_acid_to_codon_matrix[amino_acid_idx, codon_indices] = 1.0


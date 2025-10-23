import torch
import torch.nn.functional as F
from .constants import base_token_map, amino_acids_to_codons, codons, amino_acid_token_map, codon_to_rna_matrix, amino_acid_to_codon_matrix
import torch
import torch.nn.functional as F

def gumbel_softmax_hard_rna(logits, tau, amino_acid_sequence, amino_acid_to_codon_matrix, codon_to_rna_matrix):
    """
    Hard Gumbel-Softmax for RNA sequence generation, ensuring valid codons for the target amino acid sequence.
    """
    batch_size, seq_len, num_bases = logits.shape
    # print("............................................")
    # print("batch_size, seq_len, num_bases:", batch_size, seq_len, num_bases)
    assert seq_len % 3 == 0, "RNA sequence length must be a multiple of 3"

    # 1. Gumbel-Softmax Sampling
    gumbel_noise = -torch.log(-torch.log(torch.rand_like(logits) + 1e-20) + 1e-20)
    y = F.softmax((logits + gumbel_noise) / tau, dim=-1)

    # 2. Reshape to (batch_size, num_codons, 3, 4)
    y = y.view(batch_size, seq_len // 3, 3, 4)

    # 3. Retrieve valid codons for the target amino acids
    amino_acid_indices = [amino_acid_token_map[aa] for aa in amino_acid_sequence]
    valid_codon_indices = [amino_acid_to_codon_matrix[idx].nonzero(as_tuple=True)[0] for idx in amino_acid_indices]

    # 4. Create a tensor to store valid codons in one-hot form
    max_codons = max(len(c) for c in valid_codon_indices)  # Ensure uniform tensor size
    valid_codons_padded = torch.zeros(len(valid_codon_indices), max_codons, 3, 4, device=logits.device)

    for i, valid_idx in enumerate(valid_codon_indices):
        valid_codons_padded[i, :len(valid_idx)] = codon_to_rna_matrix[valid_idx]

    # 5. Compute MSE only with valid codons
    y_expanded = y.unsqueeze(4).repeat(1, 1, 1, 1, max_codons)  # Expand y to match valid codon size
    valid_codons_expanded = valid_codons_padded.unsqueeze(0).repeat(batch_size, 1, 1, 1, 1)
    valid_codons_expanded = valid_codons_expanded.permute(0, 1, 3, 4, 2)  # Align dimensions

    mse = ((y_expanded - valid_codons_expanded) ** 2).mean(dim=(2, 3))

    # 6. Mask invalid positions and select the best matching codon
    mask = valid_codons_padded.sum(dim=(2, 3)) > 0
    mask = mask.unsqueeze(0).repeat(batch_size, 1, 1)
    mse[~mask] = 1e10  # Large value for masked positions

    min_indices = mse.argmin(dim=-1)  # Select the best matching valid codon index

    # 7. Convert to One-Hot RNA using selected codons
    y_hard = torch.zeros_like(y)
    for i in range(y.shape[1]):  # Iterate over codon positions
        y_hard[:, i] = valid_codons_padded[i, min_indices[:, i]]

    # 8. Gradient Trick: Keep softmax for gradient propagation
    y_final = (y_hard - y).detach() + y

    return y_final.view(batch_size, seq_len, num_bases)


def gumbel_softmax_sample(logits, temperature=1.0):
    gumbel_noise = -torch.log(-torch.log(torch.rand_like(logits) + 1e-20) + 1e-20)
    return F.softmax((logits + gumbel_noise) / temperature, dim=-1)


def masked_gumbel_softmax(logits, mask, temperature=1.0, hard=False):
    mask = mask.to(logits.device)
    logits = logits.masked_fill(~mask, -float('inf'))
    return F.gumbel_softmax(logits, tau=temperature, hard=hard)


def amino_acids_to_rna_soft_tokens(codon_logits, mask, temperature=0.1):
    soft_codon_encodings = masked_gumbel_softmax(codon_logits, mask, temperature, hard=True)
    rna_soft_tokens = torch.einsum("bi,ijk->bjk", soft_codon_encodings, codon_to_rna_matrix)
    rna_soft_tokens = rna_soft_tokens.reshape(-1, 4)
    return rna_soft_tokens


def rna_to_amino_acids(rna_soft_tokens):
    if (rna_soft_tokens.dim() == 3):
        amino_acids = []
        for i in range(rna_soft_tokens.size(0)):
            amino_acids_i, rna_strings_i = rna_to_amino_acids(rna_soft_tokens[i, :, :])
            amino_acids.append(amino_acids_i)
            if (i > 0 and amino_acids[i] != amino_acids[0]):
                raise ValueError("Inconsistent amino acid sequences across samples")
        return "".join(amino_acids[0]), rna_strings_i

    amino_acids = []
    rna_strings = []
    for i in range(0, len(rna_soft_tokens), 3):
        codon = ""
        for j in range(3):
            base_token = torch.argmax(rna_soft_tokens[i + j]).item()
            base = list(base_token_map.keys())[list(base_token_map.values()).index(base_token)]
            codon += base
        rna_strings.append(codon)
        for amino_acid, codons in amino_acids_to_codons.items():
            if (codon in codons):
                amino_acids.append(amino_acid)
                break
    return "".join(amino_acids), rna_strings


def test_amino_acid_consistency(amino_acid_sequence, num_samples=5):
    amino_acid_indices = [amino_acid_token_map[aa] for aa in amino_acid_sequence]
    codon_logits = amino_acid_to_codon_matrix[amino_acid_indices].clone().detach().requires_grad_(True)
    mask = codon_logits > 0

    consistent = True
    for _ in range(num_samples):
        rna_soft_tokens = amino_acids_to_rna_soft_tokens(codon_logits, mask, temperature=0.5)
        decoded_amino_acids, rna_strings = rna_to_amino_acids(rna_soft_tokens)
        print(f"Sampled RNA -> Amino Acids: {decoded_amino_acids}, RNA strings: {rna_strings}")
        if (decoded_amino_acids != "".join(amino_acid_sequence)):
            consistent = False
            break
    return consistent


def decode_to_codon_level(rna_soft_tokens, codon_to_rna_matrix):
    codon_encoding = torch.einsum('bcjk,jki->bci', rna_soft_tokens, codon_to_rna_matrix)
    return codon_encoding


def rna_to_codon_onehot_unique(rna_onehot, codon_to_rna_matrix, temperature=0.1):
    batch_size, seq_len, num_bases = rna_onehot.shape
    num_codons, codon_len, _ = codon_to_rna_matrix.shape
    codon_to_rna_matrix = codon_to_rna_matrix.to(rna_onehot.device)

    if (seq_len % 3 != 0):
        raise ValueError("RNA sequence length must be a multiple of 3.")

    rna_onehot = rna_onehot.view(batch_size, seq_len // 3, codon_len, num_bases)
    rna_onehot_expanded = rna_onehot.unsqueeze(2)
    codon_to_rna_expanded = codon_to_rna_matrix.unsqueeze(0).unsqueeze(0)
    match_score = -((rna_onehot_expanded - codon_to_rna_expanded) ** 2).sum(dim=(-2, -1))
    codon_onehot = F.softmax(match_score / temperature, dim=-1)
    return codon_onehot


def calculate_rna_constraint_loss(rna_onehot, amino_acid_sequence, amino_acid_to_codon_matrix, codon_to_rna_matrix):
    batch_size, seq_len, num_bases = rna_onehot.shape
    rna_onehot = rna_onehot.reshape(batch_size, seq_len // 3, 3, num_bases)
    num_amino_acids, num_codons = amino_acid_to_codon_matrix.shape

    amino_acid_constraints = []
    amino_acid_sequence = [amino_acid_token_map[aa] for aa in amino_acid_sequence]

    for aa_idx in amino_acid_sequence:
        valid_codon_indices = amino_acid_to_codon_matrix[aa_idx].nonzero(as_tuple=True)[0]
        valid_codons = codon_to_rna_matrix[valid_codon_indices]
        amino_acid_constraints.append(valid_codons)

    max_codons = max(c.shape[0] for c in amino_acid_constraints)
    amino_acid_constraints_padded = torch.zeros(len(amino_acid_constraints), max_codons, 3, 4, device=rna_onehot.device)
    for i, valid_codons in enumerate(amino_acid_constraints):
        amino_acid_constraints_padded[i, :valid_codons.shape[0]] = valid_codons

    rna_expanded = rna_onehot.unsqueeze(4).repeat(1, 1, 1, 1, max_codons)
    amino_acid_constraints_expanded = amino_acid_constraints_padded.unsqueeze(0).repeat(batch_size, 1, 1, 1, 1)
    amino_acid_constraints_expanded = amino_acid_constraints_expanded.permute(0, 1, 3, 4, 2)

    mse = ((rna_expanded - amino_acid_constraints_expanded) ** 2).mean(dim=(2, 3))
    mask = amino_acid_constraints_padded.sum(dim=(2, 3)) > 0
    mask = mask.unsqueeze(0).repeat(batch_size, 1, 1)
    mse[~mask] = 1e10
    min_loss, min_indices = mse.min(dim=-1)
    

    hard_rna_profile = torch.zeros(batch_size, seq_len // 3, 3, 4, device=rna_onehot.device)
    for batch_idx in range(batch_size):
        for pos_idx in range(seq_len // 3):

            best_codon_idx = min_indices[batch_idx, pos_idx]

            hard_rna_profile[batch_idx, pos_idx] = amino_acid_constraints_padded[pos_idx, best_codon_idx]
    

    hard_rna_profile = hard_rna_profile.reshape(batch_size, seq_len, num_bases)
    
    return min_loss.mean(1), hard_rna_profile


def calculate_precision(codon_onehot, valid_codons):
    true_positives = torch.sum(codon_onehot * valid_codons)
    predicted_positives = torch.sum(codon_onehot)
    precision = true_positives / (predicted_positives + 1e-10)
    return precision

if __name__ == '__main__':
    amino_acid_sequence = ["M", "G", "L"]
    is_consistent = test_amino_acid_consistency(amino_acid_sequence)
    print("Is the amino acid sequence consistent across sampled RNA sequences?", is_consistent)

    amino_acid_indices = [amino_acid_token_map[aa] for aa in amino_acid_sequence]
    codon_logits = amino_acid_to_codon_matrix[amino_acid_indices].clone().detach().requires_grad_(True)
    mask = codon_logits > 0
    rna_soft_tokens = amino_acids_to_rna_soft_tokens(codon_logits, mask, temperature=0.5)

    print("RNA soft token encoding sequence:")
    print(rna_soft_tokens)

    loss = rna_soft_tokens.pow(2).sum()
    loss.backward()

    print("Codon logits gradients:")
    print(codon_logits.grad)

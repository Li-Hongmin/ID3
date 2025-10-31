# Protein Sequences Dataset

Test protein sequences used for ID3 framework demonstrations and experiments.

## Data Source

All protein sequences are sourced from **UniProt** (Universal Protein Resource):
- **Database**: Swiss-Prot (reviewed, manually annotated entries)
- **Format**: FASTA with UniProt headers
- **Website**: https://www.uniprot.org/

## Protein List

| UniProt ID | Protein Name | Organism | Length (aa) | Description |
|------------|--------------|----------|-------------|-------------|
| O15263 | Defensin beta 4A | Homo sapiens (Human) | 64 | Antimicrobial peptide |
| P00004 | Cytochrome c | Equus caballus (Horse) | 104 | Electron transfer protein |
| P01308 | Insulin | Homo sapiens (Human) | 110 | Hormone, glucose regulation |
| P01825 | Immunoglobulin heavy variable 4-59 | Homo sapiens (Human) | 117 | Antibody component |
| P0CG48 | Polyubiquitin-C | Homo sapiens (Human) | 228 | Protein degradation signal |
| P31417 | Fatty acid-binding protein 2 | Manduca sexta (Tobacco hornworm) | 122 | Lipid transport |
| P42212 | Green fluorescent protein | Aequorea victoria (Jellyfish) | 238 | Fluorescent marker |
| P61626 | Lysozyme C | Homo sapiens (Human) | 148 | Antibacterial enzyme |
| P99999 | Cytochrome c | Homo sapiens (Human) | 105 | Electron transfer protein |

## Download Instructions

### Download Individual Protein

Visit UniProt and search by ID:
```
https://www.uniprot.org/uniprotkb/{UNIPROT_ID}
```

For example:
- O15263: https://www.uniprot.org/uniprotkb/O15263
- P42212: https://www.uniprot.org/uniprotkb/P42212

Click **"Download"** → **"FASTA (canonical)"** to get the protein sequence.

### Download via Command Line

Using `curl`:
```bash
# Download single protein
curl -o O15263.fasta.txt "https://rest.uniprot.org/uniprotkb/O15263.fasta"

# Download all proteins in this dataset
for id in O15263 P00004 P01308 P01825 P0CG48 P31417 P42212 P61626 P99999; do
    curl -o ${id}.fasta.txt "https://rest.uniprot.org/uniprotkb/${id}.fasta"
done
```

### Download via Python

```python
import requests

def download_uniprot_sequence(uniprot_id, output_file):
    """Download protein sequence from UniProt"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        with open(output_file, 'w') as f:
            f.write(response.text)
        print(f"Downloaded {uniprot_id} to {output_file}")
    else:
        print(f"Failed to download {uniprot_id}: {response.status_code}")

# Example usage
download_uniprot_sequence("O15263", "O15263.fasta.txt")
```

## File Format

Each file follows the standard FASTA format with UniProt header:

```
>sp|UNIPROT_ID|ENTRY_NAME Protein Name OS=Organism OX=TaxID GN=Gene PE=Level SV=Version
AMINO_ACID_SEQUENCE_LINE_1
AMINO_ACID_SEQUENCE_LINE_2
...
```

Header fields:
- `sp`: Swiss-Prot database
- `UNIPROT_ID`: UniProt accession number
- `ENTRY_NAME`: Entry name
- `OS`: Organism species
- `OX`: NCBI taxonomy ID
- `GN`: Gene name
- `PE`: Protein existence level
- `SV`: Sequence version

## Selection Criteria

These proteins were selected for testing the ID3 framework based on:
- **Diversity**: Different protein families and functions
- **Size Range**: 64-238 amino acids (short to medium length)
- **Biological Relevance**: Well-characterized proteins with known functions
- **Species Diversity**: Human, horse, jellyfish, and insect proteins
- **Use Cases**: Therapeutics (insulin), diagnostics (GFP), antimicrobials (defensin)

## Citations

When using these protein sequences, please cite UniProt:

```bibtex
@article{uniprot2023,
  title={UniProt: the Universal Protein Knowledgebase in 2023},
  author={{The UniProt Consortium}},
  journal={Nucleic Acids Research},
  volume={51},
  number={D1},
  pages={D523--D531},
  year={2023},
  doi={10.1093/nar/gkac1052}
}
```

## License

Protein sequences from UniProt are freely available for academic and commercial use under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.

## References

- UniProt: https://www.uniprot.org/
- UniProt REST API: https://www.uniprot.org/help/api
- Swiss-Prot Documentation: https://www.uniprot.org/help/uniprotkb

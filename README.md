# Tiny

Tiny is a powerful terminal-based bioinformatics tool designed for DNA sequence analysis. It provides various features for analyzing, comparing, and discovering patterns in DNA sequences.

## Features

### 1. Basic DNA Analysis
- Sequence validation
- GC content calculation
- Molecular weight calculation
- Base composition analysis
- Complement and reverse complement sequences

### 2. Sequence Comparison
- Pairwise sequence alignment
  - Global alignment (Needleman-Wunsch algorithm)
  - Local alignment (Smith-Waterman algorithm)
  - Semi-global alignment
- Mutation detection
- Sequence identity calculation
- Gap analysis

### 3. Motif Finding
- Variable-length motif detection
- Frequency analysis
- Position tracking
- Consensus score calculation
- Custom minimum frequency thresholds

### 4. Regulatory Element Analysis
- TATA box detection
- GC box detection
- CAAT box detection
- Palindromic sequence identification
- Position information for all elements

## Installation 📦

### Prerequisites
- Python 3.9 or higher
- Poetry (Python package manager)

### Steps

1. To install Poetry on Arch Linux, you can use the following command:
```bash
sudo pacman -S python-poetry
```

Alternatively, if you prefer to install it using the official installer, you can run:
```bash
curl -sSL https://install.python-poetry.org | python3 -
```
Or head over to the official documentaion for poetry:
```bash
https://python-poetry.org/docs/
```

2. Clone the repository:
```bash
git clone https://github.com/Bjorn99/Tiny.git
cd Tiny
```

3. Install dependencies using Poetry:
```bash
poetry install
```
4. Activate the Virtual Environment:
```bash
poetry shell
```

## Usage

For a comprehensive list of examples and use cases, check out the [Examples.md](Examples.md)

### Basic Analysis
```bash
# Analyze single or multiple sequences
poetry run tiny analyze ATCG GCTA

# Analyze sequences from a FASTA file
poetry run tiny analyze --fasta sequences.fasta

# Save analysis results to a file
poetry run tiny analyze ATCG GCTA --output results.json
```

### Sequence Alignment
```bash
# Global alignment
poetry run tiny align ATCGATCG ATCTATCG --mode global

# Local alignment
poetry run tiny align ATCGATCG ATCTATCG --mode local

# Semi-global alignment
poetry run tiny align ATCGATCG ATCTATCG --mode semi-global
```

### Motif Finding
```bash
# Find motifs of length 4 that appear at least twice
poetry run tiny find-motifs ATCGATCG ATCTATCG ATCGAGCG --length 4 --min-freq 2

# Find motifs in sequences from a FASTA file
poetry run tiny find-motifs --fasta sequences.fasta --length 6 --min-freq 3
```

### Regulatory Element Analysis
```bash
# Find regulatory elements in a sequence
poetry run tiny find-regulatory TATAAAAGGCGGGCCAATATCGATCG
```

## Limitations and Considerations ⚠️

1. **Performance Limitations**
   - Not optimized for very long sequences (>10,000 bp)
   - Memory usage increases significantly with sequence length in pairwise alignments
   - Motif finding can be computationally intensive for long sequences

2. **Input Limitations**
   - Only handles DNA sequences (no RNA or protein sequences)
   - Limited to standard nucleotides (A, T, C, G)
   - No support for ambiguous bases (like N, R, Y)
   - Maximum recommended sequence length: 10,000 bp

3. **Analysis Limitations**
   - No support for multiple sequence alignment
   - Limited to basic regulatory element patterns
   - No secondary structure prediction
   - No phylogenetic analysis capabilities
   - No support for genome-scale analyses

4. **File Format Limitations**
   - Only supports FASTA format for file input
   - JSON format for output files
   - No support for other common formats (GenBank, EMBL, etc.)


## Tips for using the tool effectively:

- Always validate your input sequences before analysis
- Use appropriate alignment modes based on your sequences
- Consider sequence length limitations (max 10,000 bp)
- Save results to files for later analysis
- Use the FASTA format for multiple sequence analysis

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License 📄

This project is licensed under the GPL License - see the LICENSE file for details.

## Acknowledgments

- Built with [BioPython](https://biopython.org/)
- Project and dependency management with [Poetry](https://python-poetry.org/)
- CLI interface powered by [Typer](https://typer.tiangolo.com/)
- Terminal formatting by [Rich](https://rich.readthedocs.io/)

## Project Status

Tiny is under active development. Future planned features include:
- Support for RNA sequences
- Multiple sequence alignment
- Phylogenetic analysis
- Secondary structure prediction
- Support for additional file formats
- Performance optimizations for longer sequences
- Advanced statistical analysis
- Integration with external databases

## Support

If you encounter any issues or have questions, please:
1. Check the existing issues on GitHub
2. Create a new issue if your problem isn't already reported
3. Provide as much detail as possible about your problem
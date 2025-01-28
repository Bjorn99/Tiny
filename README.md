# Tiny

Tiny is a powerful terminal-based bioinformatics tool designed for DNA sequence analysis. It provides various features for analyzing, comparing, and discovering patterns in DNA sequences from any organism, including bacterial, fungal, viral, plant, and animal genomes.

## Features

### 1. Basic DNA Analysis
- Sequence validation with IUPAC ambiguous base support
- GC content calculation (handles ambiguous bases)
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

### 5. Enhanced Feature Analysis
- Comprehensive feature overview for GenBank files
- Feature type filtering and counting
- Customizable feature display limits
- Detailed qualifier information
- JSON export for complete feature data

### 6. File Format Support
- FASTA (.fa, .fasta)
- FASTQ (.fq, .fastq)
- GenBank (.gb, .gbk, .genbank)
- EMBL (.embl)
- JSON output format

### 7. Enhanced Visualization
- Progress bars for long operations
- Color-coded output
- Formatted tables
- Summary statistics
- Clear section separators

### 8. Feature Analysis Options
- `--feature-limit`: Control number of features displayed (0 for all)
- `--feature-type`: Filter specific feature types(CDS, gene, tRNA, etc.)
- `--save-features`: Export complete feature data to JSON
- `--format-info`: Show detailed format-specific information

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
Or head over to the official documentation for poetry:
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
tiny analyze ATCG GCTA

# Analyze sequences from files
tiny analyze --input sequence.fasta
tiny analyze --input sequence.gb --format-info

# Control feature display
tiny analyze --input sequence.gb --format-info --feature-limit 10
tiny analyze --input sequence.gb --format-info --feature-type CDS
tiny analyze --input sequence.gb --format-info --save-features

# Save analysis results to a file
tiny analyze ATCG GCTA --output results.json
```

### Sequence Alignment
```bash
# Global alignment
tiny align ATCGATCG ATCTATCG --mode global

# Local alignment
tiny align ATCGATCG ATCTATCG --mode local

# Semi-global alignment
tiny align ATCGATCG ATCTATCG --mode semi-global
```

### Motif Finding
```bash
# Find motifs of length 4 that appear at least twice
tiny find-motifs ATCGATCG ATCTATCG ATCGAGCG --length 4 --min-freq 2

# Find motifs in sequences from a file
tiny find-motifs --fasta sequences.fasta --length 6 --min-freq 3
```

### Regulatory Element Analysis
```bash
# Find regulatory elements in a sequence
tiny find-regulatory TATAAAAGGCGGGCCAATATCGATCG
```

## Limitations and Considerations ⚠️

1. **Performance Limitations**
   - Not optimized for very long sequences (>10,000 bp)
   - Memory usage increases significantly with sequence length in pairwise alignments
   - Motif finding can be computationally intensive for long sequences

2. **Input Capabilities**
   - Handles DNA sequences with IUPAC ambiguous bases
   - Supports multiple file formats (FASTA, FASTQ, GenBank, EMBL)
   - Maximum recommended sequence length: 10,000 bp

3. **Analysis Limitations**
   - No support for multiple sequence alignment
   - No secondary structure prediction
   - No phylogenetic analysis capabilities
   - No support for genome-scale analyses

## Tips for using the tool effectively:

- Validate your input sequences before analysis
- Use appropriate alignment modes based on your sequences
- Consider sequence length limitations (max 10,000 bp)
- Use format-specific information with --format-info flag
- Save results to files for later analysis
- Use file input for multiple sequence analysis

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

## References

This tool implements methods and algorithms from various scientific publications. For a complete list of references, see [REFERENCES.md](REFERENCES.md). Key references include:

- Needleman-Wunsch algorithm: Needleman & Wunsch (1970), Journal of Molecular Biology
- Smith-Waterman algorithm: Smith & Waterman (1981), Journal of Molecular Biology
- IUPAC ambiguous base notation: Cornish-Bowden (1985), Nucleic Acids Research
- Motif finding methods: Bailey & Elkan (1994), ISMB Proceedings
- Regulatory element analysis: Bucher (1990), Journal of Molecular Biology

The tool is built using BioPython (Cock et al., 2009) and other open-source libraries. For implementation details and additional references, please refer to the full references list.
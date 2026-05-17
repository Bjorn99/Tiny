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

### 5. ORF Finding and Translation
- DNA-to-protein translation in all 6 reading frames
- 7 NCBI genetic code tables (Standard, Vertebrate/Yeast/Mold/Invertebrate Mt, Ciliate, Bacterial)
- ORF detection with table-specific start/stop codons
- Nested ORF reporting (matches NCBI ORFfinder convention)
- Strand-specific translation (`--strand +1/-1`)
- `--explain` flag prints algorithm walk-throughs
- Companion Quarto notebook for deep educational dives

### 6. Enhanced Feature Analysis
- Comprehensive feature overview for GenBank files
- Feature type filtering and counting
- Customizable feature display limits
- Detailed qualifier information
- JSON export for complete feature data

### 7. File Format Support
- FASTA (.fa, .fasta)
- FASTQ (.fq, .fastq)
- GenBank (.gb, .gbk, .genbank)
- EMBL (.embl)
- JSON output format

### 8. Enhanced Visualization
- Progress bars for long operations
- Color-coded output
- Formatted tables
- Summary statistics
- Clear section separators

### 9. Feature Analysis Options
- `--feature-limit`: Control number of features displayed (0 for all)
- `--feature-type`: Filter specific feature types(CDS, gene, tRNA, etc.)
- `--save-features`: Export complete feature data to JSON
- `--format-info`: Show detailed format-specific information

## Installation 📦

### Prerequisites
- Python 3.12+ (tested on 3.12, 3.13, 3.14)
- Poetry 2.x (Python package manager)

For detailed install / run / troubleshooting steps see [INSTRUCTIONS.md](INSTRUCTIONS.md).

### Quick start

```bash
git clone https://github.com/Bjorn99/Tiny.git
cd Tiny
poetry install                 # core install (no SAM/BAM support)
poetry install --extras sam    # add pysam for SAM/BAM (Python 3.12/3.13 only)
poetry run tiny --version
poetry run tiny analyze ATCGATCG
```

To work inside the venv (Poetry 2.x):

```bash
# Poetry prints the activation command; run it. Fish example:
source (poetry env info --path)/bin/activate.fish
# Or bash/zsh:
source "$(poetry env info --path)/bin/activate"
```

## Usage

For a comprehensive list of examples and use cases, check out the [Examples.md](Examples.md). For install and troubleshooting see [INSTRUCTIONS.md](INSTRUCTIONS.md).

### Version and help
```bash
tiny --version           # print version (e.g. tiny 0.2.0)
tiny --help              # global help
tiny analyze --help      # per-command help
tiny supported-formats   # list supported file formats
```

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
tiny find-motifs --input sequences.fasta --length 6 --min-freq 3
```

### ORFs and Translation

Translate a DNA sequence to protein (NCBI standard genetic code):

```sh
tiny translate ATGAAATAG
tiny translate --input eg_files/CDKN1B_cds.fasta
tiny translate ATGAGA --table 2          # vertebrate mitochondrial code
tiny translate CTATTTCAT --strand -1     # reverse complement
tiny translate ATG --explain             # print algorithm walk-through
```

Find all open reading frames in all six reading frames:

```sh
tiny find-orfs --input eg_files/CDKN1B_cds.fasta --min-length 100
tiny find-orfs ATGAAATAG --min-length 0
tiny find-orfs ATGCCCTAA --min-length 0 --table 6  # ciliate code
tiny find-orfs ATGAAATAG --min-length 0 --explain
```

`--min-length` is in nucleotides (default `100`). `--table` selects one of
7 NCBI genetic code tables (default `1`, Standard). `--explain` prints a
short walk-through of the algorithm; the full educational write-up lives in
the companion notebook.

### Rendering the notebooks

The notebooks under `notebooks/` are Quarto documents (`.qmd`). Rendering them
requires [Quarto](https://quarto.org) and Jupyter. Install Jupyter via the
optional Poetry group:

```sh
poetry install --with notebooks
cd notebooks && quarto render 01_orf_translation.qmd
```

The rendered HTML is gitignored.

### Regulatory Element Analysis
```bash
# Find regulatory elements in a sequence
tiny find-regulatory TATAAAAGGCGGGCCAATATCGATCG
```

## Limitations and Considerations ⚠️

1. **Performance Limitations**
   - Designed for targeted-panel scale, not whole-genome data
   - Memory usage increases significantly with sequence length in pairwise alignments
   - Motif finding can be computationally intensive for long sequences

2. **Input Capabilities**
   - DNA sequences with full IUPAC ambiguous-base support
   - RNA sequences via `RNASequence` class (programmatic API; CLI exposure planned)
   - Supports multiple file formats (FASTA, FASTQ, GenBank, EMBL, plus SAM/BAM with the `sam` extra)
   - Hard cap on per-sequence length: **10,000 bp** (raises `ResourceLimitError`). Override per-call with `TINY_MAX_SEQUENCE=50000 tiny analyze ...`.

3. **Translation Limitations**
   - 7 of 33 NCBI genetic code tables shipped; remaining tables can be added on request
   - Ambiguous IUPAC codons resolve to single AA when all expansions agree, otherwise `X`
   - Known deviation from BioPython: Tiny uses `X` where BioPython uses `B`/`Z`/`J` for ambiguous amino-acid groups (Asx/Glx/Xle)
   - Selenocysteine (U) and pyrrolysine (O) not supported — these require context-dependent translation
   - ORF definition requires an in-frame stop codon; sequences truncated without a stop are not reported

4. **Analysis Limitations**
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
- Try `--explain` on any command to learn the underlying algorithm
- Select the correct genetic code table with `--table` when translating non-nuclear sequences

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

Tiny is under active revival. **Phase 1 complete** — ORF finding and translation
with 7 NCBI genetic code tables, CLI commands, and a Quarto companion notebook.
**Phase 0** (foundation hardening: typed errors, lazy heavy deps, CI, RNA class,
version flag, resource limits) shipped previously. Future planned features:
- Primer melting temperature (nearest-neighbor)
- Needleman-Wunsch / Smith-Waterman from scratch
- Position weight matrix motifs
- Burrows-Wheeler transform / FM-index
- CLI exposure for RNA sequence analysis

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

# Tiny - Comprehensive Examples Guide

This guide provides detailed examples of using Tiny for various DNA sequence analysis tasks. From basic sequence analysis to complex regulatory element detection, you'll find examples for every feature along with real-world scenarios.

## Table of Contents
- [Basic Sequence Analysis](#basic-sequence-analysis)
- [Working with Different File Formats](#working-with-different-file-formats)
- [Sequence Comparison](#sequence-comparison)
- [Sequence Alignment](#sequence-alignment)
- [Motif Finding](#motif-finding)
- [Regulatory Element Analysis](#regulatory-element-analysis)
- [Real-World Scenarios](#real-world-scenarios)
- [Working with Files](#working-with-files)
- [Tips and Best Practices](#tips-and-best-practices)

## Basic Sequence Analysis

### Single Sequence Analysis
```bash
# Basic analysis of a DNA sequence
tiny analyze ATCGATCGATCGA

# Output:
# Sequence: ATCGATCGATCG
# Length: 12
# GC Content: 50.00%
# Molecular Weight: 3674.4
# Base Composition:
#   A: 3 (25.00%)
#   T: 3 (25.00%)
#   G: 3 (25.00%)
#   C: 3 (25.00%)
```

### Multiple Sequence Analysis
```bash
# Analyze multiple sequences at once
tiny analyze ATCGATCG GCTAGCTA TATATATA

# Save results to JSON
tiny analyze ATCGATCG GCTAGCTA --output results.json
```

## Working with Different File Formats

### FASTA Files
```bash
# Create a FASTA file
cat << EOF > sequences.fasta
>Promoter_Region_1
TATAAAAGGCGGGCCAATATCGATCG
>Promoter_Region_2
CCAATGGCTAGCTAAATATATACGCG
EOF

# Analyze sequences from FASTA file
tiny analyze --input sequences.fasta
```

### GenBank Files
```bash
# Analyze GenBank file with format info
tiny analyze --input sequence.gb --format-info

# Output includes additional GenBank-specific information:
# - Features count
# - Organism information
# - References
# - Annotations
```

### FASTQ Files
```bash
# Analyze FASTQ file with quality scores
tiny analyze --input sequences.fastq --format-info

# Output includes:
# - Sequence quality scores
# - Average quality score
# - Quality statistics
```

### EMBL Files
```bash
# Analyze EMBL format
tiny analyze --input sequence.embl --format-info

# Output includes:
# - EMBL-specific annotations
# - Feature information
# - Cross-references
```

### Working with Ambiguous Bases
```bash
# Analyze sequence with IUPAC ambiguous bases
tiny analyze ATCGRYSWKMBDHVN

# The tool handles these IUPAC codes:
# R (A/G), Y (C/T), S (G/C), W (A/T)
# K (G/T), M (A/C), B (C/G/T), D (A/G/T)
# H (A/C/T), V (A/C/G), N (any base)
```

## Basic Sequence Analysis

### Single Sequence Analysis
```bash
# Basic analysis of a DNA sequence
tiny analyze ATCGATCGATCG

# Output:
# Sequence: ATCGATCGATCG
# Length: 12
# GC Content: 50.00%
# Molecular Weight: 3674.4
# Base Composition:
#   A: 3 (25.00%)
#   T: 3 (25.00%)
#   G: 3 (25.00%)
#   C: 3 (25.00%)
```

### Multiple Sequence Analysis
```bash
# Analyze multiple sequences at once
tiny analyze ATCGATCG GCTAGCTA TATATATA

# Save results to JSON
tiny analyze ATCGATCG GCTAGCTA --output results.json
```

### Working with FASTA Files
```bash
# Create a FASTA file
cat << EOF > sequences.fasta
>Promoter_Region_1
TATAAAAGGCGGGCCAATATCGATCG
>Promoter_Region_2
CCAATGGCTAGCTAAATATATACGCG
EOF

# Analyze sequences from FASTA file
tiny analyze --fasta sequences.fasta
```

## Sequence Comparison

### Basic Comparison
```bash
# Compare two similar sequences
tiny compare ATCGATCG ATCTATCG

# Output:
# Identity: 87.50%
# Mutations:
# Position 3: G → T
# Alignment Score: 14.0
```

### Multiple Mutations
```bash
# Compare sequences with multiple differences
tiny compare ATCGATCGATCG ATCTATGCATCG

# Output:
# Identity: 83.33%
# Mutations:
# Position 3: G → T
# Position 6: C → G
```

## Sequence Alignment

### Global Alignment
```bash
# Align similar length sequences
tiny align ATCGATCG ATCTATCG --mode global

# Output:
# ATCGATCG
# ATC-ATCG
# Score: 14.0
# Identity: 87.50%
# Gaps: 1
```

### Local Alignment
```bash
# Find best matching subsequence
tiny align AAATCGATCGAAA GGGGTCGATGGG --mode local

# Output:
# TCGATC
# TCGATG
# Score: 11.0
# Identity: 83.33%
# Gaps: 0
```

### Semi-Global Alignment
```bash
# Align a shorter sequence against a longer one
tiny align ATCG ATCGATCG --mode semi-global
```

## Motif Finding

### Basic Motif Search
```bash
# Find common patterns
tiny find-motifs ATCGATCGTATA ATCTATCGTATA ATCGAGCGTATA --length 4 --min-freq 2

# Output:
# Motif: ATCG
# Frequency: 3
# Positions: 0, 4, 8
# Consensus Score: 0.875
```

### Advanced Motif Analysis
```bash
# Find longer motifs with higher frequency requirement
tiny find-motifs --fasta sequences.fasta --length 6 --min-freq 3
```

## Regulatory Element Analysis

### Basic Regulatory Analysis
```bash
# Analyze a promoter region
tiny find-regulatory TATAAAAGGCGGGCCAATATCGATCG

# Output:
# TATA Box found at position 0: TATAAAA
# GC Box found at position 7: GGGCGG
```

### Complex Regulatory Search
```bash
# Analyze multiple regulatory elements
tiny find-regulatory CCAATGGCTAGCTAAATATATACGCG

# Output:
# CAAT Box found at position 0: CCAAT
# TATA Box found at position 13: TATATA
```

## Real-World Scenarios

### Scenario 1: Promoter Analysis
```bash
# 1. First, analyze the sequence
tiny analyze TATAAAAGGCGGGCCAATATCGATCG

# 2. Look for regulatory elements
tiny find-regulatory TATAAAAGGCGGGCCAATATCGATCG

# 3. Find conserved motifs
tiny find-motifs TATAAAAGGCGGGCCAATATCGATCG TATAAAAGGCGGGCCAATATCGATCT --length 6 --min-freq 2
```

### Scenario 2: Mutation Analysis
```bash
# 1. Compare wild-type and mutant sequences
tiny compare ATCGATCGTATA ATCTATCGTATA

# 2. Align sequences to visualize mutations
tiny align ATCGATCGTATA ATCTATCGTATA --mode global

# 3. Look for affected regulatory elements
tiny find-regulatory ATCGATCGTATA
tiny find-regulatory ATCTATCGTATA
```

### Scenario 3: Sequence Conservation Study
```bash
# Create a file with multiple sequence variants
cat << EOF > variants.fasta
>Variant1
ATCGATCGTATA
>Variant2
ATCTATCGTATA
>Variant3
ATCGAGCGTATA
EOF

# Find conserved motifs
tiny find-motifs --fasta variants.fasta --length 4 --min-freq 3
```

## Working with Files

### File Input/Output
```bash
# Input from FASTA
tiny analyze --fasta sequences.fasta

# Output to JSON
tiny analyze ATCGATCG --output analysis_results.json

# Batch processing
for file in *.fasta; do
    tiny analyze --fasta "$file" --output "${file%.fasta}_results.json"
done
```

## Tips and Best Practices

1. Working with Different File Formats
```bash
# Check supported formats
tiny supported-formats

# Use format-specific information
tiny analyze --input sequence.gb --format-info
```

2. Handling Ambiguous Bases
```bash
# The tool now handles IUPAC ambiguous bases
tiny analyze ATCGRYSWN
tiny analyze --input sequence.gb  # Often contains ambiguous bases
```

3. Sequence Length Considerations
   - Keep sequences under 10,000 bp for optimal performance
   - Use appropriate alignment modes based on sequence length

4. File Management
   ```bash
   # Organize results by date
   mkdir -p results/$(date +%Y-%m-%d)
   tiny analyze --input sequence.gb --output results/$(date +%Y-%m-%d)/analysis.json
   ```

5. Quality Control
   ```bash
   # View format-specific information
   tiny analyze --input sequence.fastq --format-info
   ```

4. Performance Optimization
   ```bash
   # For multiple sequences, use FASTA files instead of command line input
   # Bad:
   tiny analyze ATCG GCTA CGTA TAGC
   
   # Good:
   tiny analyze --fasta sequences.fasta
   ```

## Common Issues and Solutions

### 1. Sequence Input Issues

#### Invalid Base Characters
```bash
# This will fail
tiny analyze ATCGXTCG
Error: Invalid DNA sequence. Found invalid bases: X.

# This will work (standard bases)
tiny analyze ATCGATCG

# This will also work (IUPAC ambiguous bases)
tiny analyze ATCGRYSWN
```

#### File Format Issues
```bash
# Wrong file extension
tiny analyze --input sequence.txt
Error: Unsupported file format: .txt

# Correct usage with supported formats
tiny analyze --input sequence.fasta
tiny analyze --input sequence.gb
tiny analyze --input sequence.fastq
tiny analyze --input sequence.embl
```

### 2. Memory and Performance Issues

#### Large Sequence Processing
```bash
# For large sequences, split into smaller chunks
# Process 1000 bp at a time
tiny analyze --input subset1.fasta
tiny analyze --input subset2.fasta

# Better alternative: use file input instead of command line
# Bad (memory intensive):
tiny analyze ATCG[...very long sequence...]

# Good:
tiny analyze --input sequence.fasta
```

#### Multiple Sequence Analysis
```bash
# Inefficient way (might cause memory issues):
tiny analyze seq1.fasta seq2.fasta seq3.fasta

# Better approach:
# Combine sequences into one FASTA file
cat seq1.fasta seq2.fasta seq3.fasta > combined.fasta
tiny analyze --input combined.fasta
```

### 3. File Format Specific Issues

#### GenBank Files


##### Basic Feature Overview
```bash
# Missing sequence data
Error: No sequence found in GenBank record

# Solution: Ensure GenBank file contains sequence data
# Check file content and format

# Handling features
tiny analyze --input sequence.gb --format-info

# Output includes:
# - Feature type summary
# - Feature counts
# - Sample features with locations
# - Basic qualifiers
```

##### Controlling Feature Display
```bash
# Show all features (no limit)
tiny analyze --input sequence.gb --format-info --feature-limit 0

# Show first 10 features of each type
tiny analyze --input sequence.gb --format-info --feature-limit 10

# Filter specific feature type
tiny analyze --input sequence.gb --format-info --feature-type CDS
tiny analyze --input sequence.gb --format-info --feature-type gene
tiny analyze --input sequence.gb --format-info --feature-type tRNA

# Save complete feature information
tiny analyze --input sequence.gb --format-info --save-features

# Combined options
tiny analyze --input sequence.gb --format-info --feature-type CDS --feature-limit 10 --save-features
```
##### Feature Analysis Output
```# Feature Type Summary shows:
Feature Type   Count
CDS            1709
gene          1775
rRNA            18
source           1
tRNA            18

# Detailed Features Overview shows:
- Feature locations
- Qualifiers
- Additional annotations
- Cross-references
```

##### Working with Feature Data
```bash
# Export features for further analysis
tiny analyze --input sequence.gb --format-info --save-features
# Creates features.json with complete data

# Analyze specific features
tiny analyze --input sequence.gb --format-info --feature-type CDS > cds_analysis.txt
```

#### FASTQ Quality Scores
```bash
# Invalid quality scores
Error: Invalid quality scores in FASTQ

# Solution: Verify FASTQ format and quality encoding
tiny analyze --input sequence.fastq --format-info
```

### 4. Analysis-Specific Issues

#### Alignment Problems
```bash
# Sequences too different
tiny align ATCG GGGG
# Solution: Use local alignment for divergent sequences
tiny align ATCG GGGG --mode local

# Memory error with long sequences
# Solution: Use smaller sequences or split alignment
```

#### Motif Finding Issues
```bash
# No motifs found
tiny find-motifs --input seqs.fasta --length 6 --min-freq 3
Error: No motifs found with given parameters

# Solutions:
# 1. Decrease motif length
tiny find-motifs --input seqs.fasta --length 4 --min-freq 3

# 2. Decrease minimum frequency
tiny find-motifs --input seqs.fasta --length 6 --min-freq 2
```

### 5. Output and Formatting Issues

#### JSON Output Problems
```bash
# Permission denied
tiny analyze ATCG --output /path/to/results.json
Error: Permission denied

# Solutions:
# 1. Check directory permissions
# 2. Use a different output location
tiny analyze ATCG --output ./results.json
```

#### Display Issues
```bash
# Output too wide for terminal
# Solution: Use a wider terminal or redirect to file
tiny analyze --input sequence.gb --format-info > output.txt
```

### 6. Installation and Dependencies

#### Poetry Installation Issues
```bash
# Poetry installation fails
Error: Permission denied

# Solution 1: Use system package manager
sudo pacman -S python-poetry  # For Arch Linux

# Solution 2: Install for current user only
curl -sSL https://install.python-poetry.org | python3 - --user
```

#### Dependency Conflicts
```bash
# Dependency resolution error
Error: Unable to resolve dependencies

# Solution:
poetry update  # Update all dependencies
poetry install --no-cache  # Clean install
```

### 7. Common Command Usage Mistakes

#### Incorrect Flag Usage
```bash
# Wrong:
tiny analyze -fasta sequence.fasta
Error: Unknown option '-fasta'

# Correct:
tiny analyze --input sequence.fasta
```

#### Mode Specification
```bash
# Wrong:
tiny align ATCG GCTA --mode GLOBAL
Error: Invalid alignment mode

# Correct:
tiny align ATCG GCTA --mode global
```

#### Feature Analysis Issues

##### Large Feature Sets
```bash
# Problem: Too many features displayed
tiny analyze --input sequence.gb --format-info
# Shows limited features by default

# Solution 1: Adjust display limit
tiny analyze --input sequence.gb --format-info --feature-limit 10

# Solution 2: Filter specific features
tiny analyze --input sequence.gb --format-info --feature-type CDS

# Solution 3: Save complete data for external analysis
tiny analyze --input sequence.gb --format-info --save-features
```
##### Feature Type Selection
```bash
# Problem: Unknown feature types
# First, check available types:
tiny analyze --input sequence.gb --format-info
# Shows Feature Type Summary

# Then filter specific type:
tiny analyze --input sequence.gb --format-info --feature-type gene
```

##### Feature Data Export
```bash
# JSON export
tiny analyze --input sequence.gb --format-info --save-features
# Creates features.json

# Combined with filtering
tiny analyze --input sequence.gb --format-info --feature-type CDS --save-features
# Creates filtered features.json
```

Remember to check the help documentation for each command:
```bash
tiny --help                    # General help
tiny analyze --help           # Analysis command help
tiny align --help            # Alignment command help
tiny find-motifs --help      # Motif finding help
tiny find-regulatory --help  # Regulatory element help
tiny supported-formats       # Show supported file formats
```
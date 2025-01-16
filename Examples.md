# Tiny - Comprehensive Examples Guide

This guide provides detailed examples of using Tiny for various DNA sequence analysis tasks. From basic sequence analysis to complex regulatory element detection, you'll find examples for every feature along with real-world scenarios.

## Table of Contents
- [Basic Sequence Analysis](#basic-sequence-analysis)
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

1. Sequence Length Considerations
   - Keep sequences under 10,000 bp for optimal performance
   - Use appropriate alignment modes based on sequence length

2. File Management
   ```bash
   # Organize results in directories
   mkdir -p results/$(date +%Y-%m-%d)
   tiny analyze --fasta sequences.fasta --output results/$(date +%Y-%m-%d)/analysis.json
   ```

3. Quality Control
   ```bash
   # Always validate your sequences first
   tiny analyze ATCGATCG
   
   # Then proceed with detailed analysis
   tiny find-regulatory ATCGATCG
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

1. Invalid Sequences
   ```bash
   # This will fail (invalid base 'X')
   tiny analyze ATCGXTCG
   
   # This will work
   tiny analyze ATCGATCG
   ```

2. Memory Usage
   ```bash
   # For large sequences, split analysis into smaller chunks
   # Process 1000 bp at a time
   tiny analyze --fasta subset1.fasta
   tiny analyze --fasta subset2.fasta
   ```

Remember to check the `--help` option for any command to see all available options:
```bash
tiny --help
tiny analyze --help
tiny align --help
tiny find-motifs --help
tiny find-regulatory --help
```
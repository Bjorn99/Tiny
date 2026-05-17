"""Condensed algorithm walk-throughs for --explain. The notebook carries the deep dive."""

TRANSLATE_EXPLANATION = """\
Translation in one paragraph
============================
A reading frame slides over DNA in steps of three. Each three-letter window
is a **codon**, and a genetic code table maps each codon to one of 20 amino
acids or a stop signal ('*'). Translation reads codons left to right and
emits the amino-acid string.

Tiny ships with 7 NCBI genetic code tables (pass --table N). Table 1 is the
Standard code used by most nuclear genomes. Table 2 (Vertebrate Mitochondrial)
reassigns AGA/AGG to STOP, AUA to Met, and UGA to Trp — mitochondrial genomes
deviate this way because their translation machinery is simpler.

A trailing partial codon (1 or 2 leftover bases after frame offset) is
dropped with a warning. Codons containing ambiguous IUPAC bases (anything
outside ACGT) translate to 'X' — the same convention BioPython uses.

Why six frames? DNA is double-stranded, and the ribosome can begin reading at
offset 0, 1, or 2 on either strand. That gives 3 + 3 = 6 reading frames, any
of which could encode a protein. Use `find-orfs` to search all six at once.
"""

FIND_ORFS_EXPLANATION = """\
ORF finding in one paragraph
============================
An *open reading frame* is a stretch of DNA that starts at a start codon
and ends at the first in-frame stop codon. The algorithm:

  1. Scan the forward strand and the reverse complement.
  2. For each of the 3 offsets (0, 1, 2) on each strand, walk the codons.
  3. When you hit a start codon (table-specific), scan forward until the
     next in-frame stop codon in the same frame.
  4. Record the span [start ... stop] as one ORF.

Nested ORFs (multiple start codons that all end at the same stop) are all
reported, matching the convention used by NCBI ORFfinder. The --min-length
filter removes ORFs shorter than the threshold in nucleotides — the default
of 100 nt corresponds roughly to the smallest biologically interesting CDS.

Start codons vary by genetic code. The Standard code recognises ATG, CTG,
and TTG as starts. Vertebrate mitochondria additionally use ATA, ATC, ATT,
and GTG. Use --table to select the appropriate code for your organism.
"""

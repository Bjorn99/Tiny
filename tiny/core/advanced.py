from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from Bio import motifs
from Bio import Align
from Bio.Seq import Seq
import re
from collections import Counter

@dataclass
class AlignmentResult:
    """Class to hold sequence alignment results."""
    sequence1: str
    sequence2: str
    aligned_seq1: str
    aligned_seq2: str
    alignment_score: float
    identity: float
    gaps: int
    alignment_length: int

@dataclass
class MotifResult:
    """Class to hold motif finding results."""
    motif: str
    positions: List[int]
    frequency: int
    consensus_score: float

class AdvancedAnalysis:
    @staticmethod
    def align_sequences(seq1: str, seq2: str, mode: str = 'global') -> AlignmentResult:
        """
        Align two sequences using different alignment modes.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            mode: Alignment mode ('global', 'local', 'semi-global')
        """
        aligner = Align.PairwiseAligner()
        
        # Set scoring parameters
        aligner.match_score = 2.0
        aligner.mismatch_score = -1.0
        aligner.open_gap_score = -10.0  # Higher penalty for opening gaps
        aligner.extend_gap_score = -0.5  # Lower penalty for extending gaps
        
        # Set mode
        if mode == 'global':
            aligner.mode = 'global'
        elif mode == 'local':
            aligner.mode = 'local'
        else:  # semi-global
            aligner.mode = 'global'
            aligner.query_end_gap_score = 0  # No penalty for gaps at sequence ends
            aligner.target_end_gap_score = 0

        # Get best alignment
        alignments = aligner.align(seq1, seq2)
        if not alignments:
            raise ValueError("No alignment found")

        best_alignment = alignments[0]
        
        # Format aligned sequences
        formatted = best_alignment.format()
        aligned_seqs = formatted.split('\n')
        
        # Extract aligned sequences (removing the middle line with matches)
        if len(aligned_seqs) >= 3:
            aligned_seq1 = aligned_seqs[0]
            aligned_seq2 = aligned_seqs[2]
        else:
            raise ValueError("Unexpected alignment format")

        # Calculate statistics
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
        gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
        alignment_length = len(aligned_seq1)
        identity = (matches / alignment_length) * 100

        return AlignmentResult(
            sequence1=seq1,
            sequence2=seq2,
            aligned_seq1=aligned_seq1,
            aligned_seq2=aligned_seq2,
            alignment_score=best_alignment.score,
            identity=identity,
            gaps=gaps,
            alignment_length=alignment_length
        )

    @staticmethod
    def find_motifs(sequences: List[str], 
                    motif_length: int, 
                    min_frequency: int = 2) -> List[MotifResult]:
        """
        Find common motifs in a set of sequences.
        
        Args:
            sequences: List of DNA sequences
            motif_length: Length of motifs to search for
            min_frequency: Minimum frequency for a motif to be reported
        """
        motif_counts = Counter()
        motif_positions: Dict[str, List[int]] = {}

        # Scan each sequence for potential motifs
        for seq in sequences:
            for i in range(len(seq) - motif_length + 1):
                motif = seq[i:i + motif_length]
                motif_counts[motif] += 1
                if motif not in motif_positions:
                    motif_positions[motif] = []
                motif_positions[motif].append(i)

        # Filter motifs by minimum frequency
        common_motifs = []
        for motif, count in motif_counts.items():
            if count >= min_frequency:
                # Calculate consensus score based on conservation
                consensus_score = AdvancedAnalysis._calculate_consensus_score(motif)
                common_motifs.append(MotifResult(
                    motif=motif,
                    positions=motif_positions[motif],
                    frequency=count,
                    consensus_score=consensus_score
                ))

        # Sort by frequency and consensus score
        return sorted(common_motifs, 
                     key=lambda x: (x.frequency, x.consensus_score), 
                     reverse=True)

    @staticmethod
    def _calculate_consensus_score(motif: str) -> float:
        """Calculate a conservation score for a motif."""
        base_weights = {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}
        score = 0
        for base in motif:
            score += base_weights[base]
        return score / len(motif)

    @staticmethod
    def find_regulatory_elements(sequence: str) -> Dict[str, List[Tuple[int, str]]]:
        """
        Finding potential regulatory elements in a DNA sequence.
        """
        regulatory_elements = {
            'TATA_box': [],
            'GC_box': [],
            'CAAT_box': [],
            'palindromes': []
        }

        # TATA box (common promoter sequence)
        tata_pattern = r'TATA[AT]A[AT]'
        for match in re.finditer(tata_pattern, sequence):
            regulatory_elements['TATA_box'].append((match.start(), match.group()))

        # GC box
        gc_pattern = r'GGGCGG|CCGCCC'
        for match in re.finditer(gc_pattern, sequence):
            regulatory_elements['GC_box'].append((match.start(), match.group()))

        # CAAT box
        caat_pattern = r'CCAAT|ATTGG'
        for match in re.finditer(caat_pattern, sequence):
            regulatory_elements['CAAT_box'].append((match.start(), match.group()))

        # Finding palindromic sequences (min length 6)
        for i in range(len(sequence) - 5):
            for j in range(6, min(11, len(sequence) - i + 1)):
                substr = sequence[i:i+j]
                if AdvancedAnalysis._is_palindrome(substr):
                    regulatory_elements['palindromes'].append((i, substr))

        return regulatory_elements

    @staticmethod
    def _is_palindrome(seq: str) -> bool:
        """Checking if a sequence is palindromic."""
        rev_comp = str(Seq(seq).reverse_complement())
        return seq == rev_comp
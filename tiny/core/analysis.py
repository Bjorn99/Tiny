from dataclasses import dataclass
from typing import List, Optional
from Bio import Align
from tiny.core.sequence import DNASequence

@dataclass
class AnalysisResult:
    """Class to hold sequence analysis results."""
    sequence: DNASequence
    length: int
    gc_content: float
    molecular_weight: float
    base_composition: dict

@dataclass
class Mutation:
    """Class to represent a mutation between sequences."""
    position: int
    original: str
    mutated: str

@dataclass
class ComparisonResult:
    """Class to hold sequence comparison results."""
    identity: float
    alignment_score: float
    mutations: List[Mutation]

def analyze_sequences(sequences: List[DNASequence]) -> List[AnalysisResult]:
    """Analyze a list of DNA sequences."""
    results = []
    for seq in sequences:
        result = AnalysisResult(
            sequence=seq,
            length=len(seq),
            gc_content=seq.gc_content,
            molecular_weight=seq.molecular_weight,
            base_composition=seq.base_composition
        )
        results.append(result)
    return results

def compare_sequences(seq1: DNASequence, seq2: DNASequence) -> ComparisonResult:
    """Compare two DNA sequences."""
    # Calculate sequence identity

    aligner = Align.PairwiseAligner()

    alignments = aligner.align(str(seq1), str(seq2))
    if alignments:
        best_alignment = alignments[0]
        alignment_score = best_alignment.score
        identity = (alignment_score / max(len(seq1), len(seq2))) * 100
    else:
        alignment_score = 0
        identity = 0

    # Find mutations
    mutations = []
    min_length = min(len(seq1), len(seq2))
    for i in range(min_length):
        if seq1.sequence[i] != seq2.sequence[i]:
            mutations.append(Mutation(
                position=i,
                original=seq1.sequence[i],
                mutated=seq2.sequence[i]
            ))

    return ComparisonResult(
        identity=identity,
        alignment_score=alignment_score,
        mutations=mutations
    )
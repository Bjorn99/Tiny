from dataclasses import dataclass
from typing import Dict, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

@dataclass
class DNASequence:
    """Class for representing and analyzing DNA sequences with IUPAC support."""
    
    # IUPAC ambiguous DNA codes and their meanings
    IUPAC_CODES = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],    # Purine
        'Y': ['C', 'T'],    # Pyrimidine
        'M': ['A', 'C'],    # Amino
        'K': ['G', 'T'],    # Keto
        'S': ['C', 'G'],    # Strong
        'W': ['A', 'T'],    # Weak
        'B': ['C', 'G', 'T'],   # not A
        'D': ['A', 'G', 'T'],   # not C
        'H': ['A', 'C', 'T'],   # not G
        'V': ['A', 'C', 'G'],   # not T
        'N': ['A', 'C', 'G', 'T']  # any base
    }

    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self._validate()
        # Create Seq object with molecule type annotation
        self._bio_seq = Seq(self.sequence)
        self._bio_seq.annotations = {"molecule_type": "DNA"}
        self._cached_gc_content: Optional[float] = None
        self._cached_molecular_weight: Optional[float] = None

    def _validate(self):
        """Validate the DNA sequence with IUPAC codes."""
        invalid_bases = set(self.sequence) - set(self.IUPAC_CODES.keys())
        if invalid_bases:
            raise ValueError(
                f"Invalid DNA sequence. Found invalid bases: {', '.join(invalid_bases)}.\n"
                f"Allowed bases: {', '.join(self.IUPAC_CODES.keys())}\n"
                "This includes IUPAC ambiguity codes (R,Y,M,K,S,W,B,D,H,V,N)."
            )

    def _count_gc_weighted(self, base: str) -> float:
        """Calculate weighted GC content for ambiguous bases."""
        possibilities = self.IUPAC_CODES[base]
        gc_count = sum(1 for b in possibilities if b in 'GC')
        return gc_count / len(possibilities)

    @property
    def gc_content(self) -> float:
        """
        Calculate GC content with caching, handling ambiguous bases.
        For ambiguous bases, we take the average probability of G/C content.
        """
        if self._cached_gc_content is None:
            if not self.sequence:
                self._cached_gc_content = 0.0
            else:
                total_gc = sum(self._count_gc_weighted(base) for base in self.sequence)
                self._cached_gc_content = (total_gc / len(self.sequence)) * 100
        return self._cached_gc_content

    @property
    def molecular_weight(self) -> float:
        """
        Calculate approximate molecular weight of the sequence.
        For ambiguous bases, uses average weight.
        """
        # Base weights in g/mol (average for ambiguous bases)
        base_weights = {
            'A': 331.2,
            'C': 307.2,
            'G': 347.2,
            'T': 322.2,
            'N': 327.0,  # Average weight
            'R': 339.2,  # (A+G)/2
            'Y': 314.7,  # (C+T)/2
            'M': 319.2,  # (A+C)/2
            'K': 334.7,  # (G+T)/2
            'S': 327.2,  # (G+C)/2
            'W': 326.7,  # (A+T)/2
            'B': 325.5,  # (C+G+T)/3
            'D': 333.5,  # (A+G+T)/3
            'H': 320.2,  # (A+C+T)/3
            'V': 328.5   # (A+C+G)/3
        }
        
        total_weight = sum(base_weights[base] for base in self.sequence)
        # Add weight of phosphate backbone
        total_weight += (len(self.sequence) - 1) * 174.0
        return total_weight

    @property
    def base_composition(self) -> Dict[str, int]:
        """Get the composition of each base, including ambiguous bases."""
        composition = {}
        for base in self.IUPAC_CODES.keys():
            count = self.sequence.count(base)
            if count > 0:
                composition[base] = count
        return composition

    @property
    def complement(self) -> 'DNASequence':
        """Get the complement sequence."""
        return DNASequence(str(self._bio_seq.complement()))

    @property
    def reverse_complement(self) -> 'DNASequence':
        """Get the reverse complement sequence."""
        return DNASequence(str(self._bio_seq.reverse_complement()))

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def __eq__(self, other: 'DNASequence') -> bool:
        return self.sequence == other.sequence
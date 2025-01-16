from dataclasses import dataclass
from typing import Dict, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

@dataclass
class DNASequence:
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self._validate()
        self._bio_seq = Seq(self.sequence)
        self._cached_gc_content: Optional[float] = None
        self._cached_molecular_weight: Optional[float] = None

    def _validate(self):
        """Validate the DNA sequence."""
        valid_bases = set('ATCG')
        invalid_bases = set(self.sequence) - valid_bases
        if invalid_bases:
            raise ValueError(
                f"Invalid DNA sequence. Found invalid bases: {', '.join(invalid_bases)}. "
                "Only A, T, C, G are allowed."
            )

    @property
    def gc_content(self) -> float:
        """Calculate GC content with caching."""
        if self._cached_gc_content is None:
            if not self.sequence:
                self._cached_gc_content = 0.0
            else:
                gc_count = self.sequence.count('G') + self.sequence.count('C')
                self._cached_gc_content = (gc_count / len(self.sequence)) * 100
        return self._cached_gc_content

    @property
    def molecular_weight(self) -> float:
        """Calculate molecular weight of the sequence."""
        return molecular_weight(self._bio_seq)

    @property
    def base_composition(self) -> Dict[str, int]:
        """Get the composition of each base."""
        return {
            'A': self.sequence.count('A'),
            'T': self.sequence.count('T'),
            'G': self.sequence.count('G'),
            'C': self.sequence.count('C')
        }

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
from dataclasses import dataclass
from typing import Dict
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

@dataclass
class DNASequence:
    """Class for representing and analyzing DNA sequences."""
    sequence: str

    def __post_init__(self):
        self.sequence = self.sequence.upper()
        self._validate()
        self._bio_seq = Seq(self.sequence)

    def _validate(self):
        """Validate the DNA sequence."""
        valid_bases = set('ATCG')
        if not all(base in valid_bases for base in self.sequence):
            raise ValueError("Invalid DNA sequence. Only A, T, C, G are allowed.")

    @property
    def gc_content(self) -> float:
        """Calculate GC content as percentage."""
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / len(self.sequence)) * 100

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
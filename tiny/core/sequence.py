import os
from dataclasses import dataclass
from typing import ClassVar

from Bio.Seq import Seq


@dataclass
class DNASequence:
    """Class for representing and analyzing DNA sequences with IUPAC support."""

    # IUPAC ambiguous DNA codes and their meanings
    IUPAC_CODES: ClassVar[dict[str, list[str]]] = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A", "G"],  # Purine
        "Y": ["C", "T"],  # Pyrimidine
        "M": ["A", "C"],  # Amino
        "K": ["G", "T"],  # Keto
        "S": ["C", "G"],  # Strong
        "W": ["A", "T"],  # Weak
        "B": ["C", "G", "T"],  # not A
        "D": ["A", "G", "T"],  # not C
        "H": ["A", "C", "T"],  # not G
        "V": ["A", "C", "G"],  # not T
        "N": ["A", "C", "G", "T"],  # any base
    }

    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        self._validate()
        # Create Seq object with molecule type annotation
        self._bio_seq = Seq(self.sequence)
        self._bio_seq.annotations = {"molecule_type": "DNA"}
        self._cached_gc_content: float | None = None
        self._cached_molecular_weight: float | None = None

    def _validate(self):
        """Validate the sequence against the class's alphabet and size limit."""
        from tiny.core.errors import InvalidSequenceError, ResourceLimitError

        max_len = int(os.environ.get("TINY_MAX_SEQUENCE", "10000"))
        if len(self.sequence) > max_len:
            raise ResourceLimitError("sequence length", actual=len(self.sequence), limit=max_len)

        invalid_bases = set(self.sequence) - set(self.IUPAC_CODES.keys())
        if invalid_bases:
            raise InvalidSequenceError(
                f"Invalid {self.__class__.__name__} bases: "
                f"{', '.join(sorted(invalid_bases))}. "
                f"Allowed: {', '.join(sorted(self.IUPAC_CODES.keys()))}"
            )

    def _count_gc_weighted(self, base: str) -> float:
        """Calculate weighted GC content for ambiguous bases."""
        possibilities = self.IUPAC_CODES[base]
        gc_count = sum(1 for b in possibilities if b in "GC")
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
            "A": 331.2,
            "C": 307.2,
            "G": 347.2,
            "T": 322.2,
            "U": 324.2,  # uracil monophosphate — lets RNA inherit this method
            "N": 327.0,  # Average weight
            "R": 339.2,  # (A+G)/2
            "Y": 314.7,  # (C+T)/2
            "M": 319.2,  # (A+C)/2
            "K": 334.7,  # (G+T)/2
            "S": 327.2,  # (G+C)/2
            "W": 326.7,  # (A+T)/2
            "B": 325.5,  # (C+G+T)/3
            "D": 333.5,  # (A+G+T)/3
            "H": 320.2,  # (A+C+T)/3
            "V": 328.5,  # (A+C+G)/3
        }

        total_weight = sum(base_weights[base] for base in self.sequence)
        # Add weight of phosphate backbone
        total_weight += (len(self.sequence) - 1) * 174.0
        return total_weight

    @property
    def base_composition(self) -> dict[str, int]:
        """Get the composition of each base, including ambiguous bases."""
        composition = {}
        for base in self.IUPAC_CODES:
            count = self.sequence.count(base)
            if count > 0:
                composition[base] = count
        return composition

    @property
    def complement(self) -> "DNASequence":
        """Get the complement sequence."""
        return DNASequence(str(self._bio_seq.complement()))

    @property
    def reverse_complement(self) -> "DNASequence":
        """Get the reverse complement sequence."""
        return DNASequence(str(self._bio_seq.reverse_complement()))

    def __len__(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return self.sequence

    def __eq__(self, other: "DNASequence") -> bool:
        return self.sequence == other.sequence


class RNASequence(DNASequence):
    """RNA sequence with U replacing T. Inherits DNA logic with overrides."""

    IUPAC_CODES: ClassVar[dict[str, list[str]]] = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "U": ["U"],
        "R": ["A", "G"],
        "Y": ["C", "U"],
        "M": ["A", "C"],
        "K": ["G", "U"],
        "S": ["C", "G"],
        "W": ["A", "U"],
        "B": ["C", "G", "U"],
        "D": ["A", "G", "U"],
        "H": ["A", "C", "U"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "U"],
    }

    _COMPLEMENT = str.maketrans("AUCGRYMKSWBDHVN", "UAGCYRKMSWVHDBN")

    def __init__(self, sequence: str):
        # Skip DNASequence's BioPython Seq construction (it expects T not U).
        self.sequence = sequence.upper()
        self._validate()
        self._cached_gc_content: float | None = None
        self._cached_molecular_weight: float | None = None

    @property
    def complement(self) -> "RNASequence":
        return RNASequence(self.sequence.translate(self._COMPLEMENT))

    @property
    def reverse_complement(self) -> "RNASequence":
        return RNASequence(self.sequence.translate(self._COMPLEMENT)[::-1])

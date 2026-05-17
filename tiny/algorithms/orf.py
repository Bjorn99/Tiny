"""ORF finding and translation. Pure-Python leaf module; no pysam, no Rich.

Each genetic code table is stored as an immutable dict with codons, start_codons,
and stop_codons. Data cross-validated against Bio.Data.CodonTable (BioPython).
"""

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Genetic code tables
# ---------------------------------------------------------------------------
# Keys are NCBI translation table numbers. Each value is a dict with:
#   name          — human-readable name
#   codons        — dict mapping 3-nt uppercase DNA codons to single-letter AA
#   start_codons  — list of codons that can serve as translation start
#   stop_codons   — list of codons that signal termination
# ---------------------------------------------------------------------------

# Per-table codon differences from the Standard code (table 1).
# This avoids duplicating the 64-entry dict for every table.
_STANDARD_CODONS: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

_STANDARD_STARTS = ["ATG", "CTG", "TTG"]
_STANDARD_STOPS = ["TAA", "TAG", "TGA"]

# Codon table overrides keyed by (table_id, codon) -> amino_acid.
# Only entries that DIFFER from the Standard code are listed.
_CODON_OVERRIDES: dict[tuple[int, str], str] = {
    # Table 2 — Vertebrate Mitochondrial
    (2, "AGA"): "*",
    (2, "AGG"): "*",
    (2, "ATA"): "M",
    (2, "TGA"): "W",
    # Table 3 — Yeast Mitochondrial
    (3, "ATA"): "M",
    (3, "TGA"): "W",
    (3, "CTA"): "T",
    (3, "CTC"): "T",
    (3, "CTG"): "T",
    (3, "CTT"): "T",
    # Table 4 — Mold / Protozoan / Coelenterate Mitochondrial
    (4, "TGA"): "W",
    # Table 5 — Invertebrate Mitochondrial
    (5, "AGA"): "S",
    (5, "AGG"): "S",
    (5, "ATA"): "M",
    (5, "TGA"): "W",
    # Table 6 — Ciliate / Dasycladacean / Hexamita Nuclear
    (6, "TAA"): "Q",
    (6, "TAG"): "Q",
    # Table 11 — Bacterial / Archaeal / Plant Plastid (codons same as standard)
}

_TABLE_METADATA: dict[int, dict] = {
    1: {
        "name": "Standard",
        "start_codons": ["ATG", "CTG", "TTG"],
        "stop_codons": ["TAA", "TAG", "TGA"],
    },
    2: {
        "name": "Vertebrate Mitochondrial",
        "start_codons": ["ATA", "ATC", "ATG", "ATT", "GTG"],
        "stop_codons": ["AGA", "AGG", "TAA", "TAG"],
    },
    3: {
        "name": "Yeast Mitochondrial",
        "start_codons": ["ATA", "ATG", "GTG"],
        "stop_codons": ["TAA", "TAG"],
    },
    4: {
        "name": "Mold Mitochondrial",
        "start_codons": ["ATA", "ATC", "ATG", "ATT", "CTG", "GTG", "TTA", "TTG"],
        "stop_codons": ["TAA", "TAG"],
    },
    5: {
        "name": "Invertebrate Mitochondrial",
        "start_codons": ["ATA", "ATC", "ATG", "ATT", "GTG", "TTG"],
        "stop_codons": ["TAA", "TAG"],
    },
    6: {"name": "Ciliate Nuclear", "start_codons": ["ATG"], "stop_codons": ["TGA"]},
    11: {
        "name": "Bacterial",
        "start_codons": ["ATA", "ATC", "ATG", "ATT", "CTG", "GTG", "TTG"],
        "stop_codons": ["TAA", "TAG", "TGA"],
    },
}


def _build_table(table_id: int) -> dict[str, str]:
    """Build the full codon->AA dict for a table by applying overrides to the Standard code."""
    codons = dict(_STANDARD_CODONS)
    for (tid, codon), aa in _CODON_OVERRIDES.items():
        if tid == table_id:
            codons[codon] = aa
    return codons


# Public immutable registry — gen_code[table_id]["codons"]["ATG"]
GENETIC_CODES: dict[int, dict] = {}
for _tid in sorted(_TABLE_METADATA):
    GENETIC_CODES[_tid] = {
        "name": _TABLE_METADATA[_tid]["name"],
        "codons": _build_table(_tid),
        "start_codons": frozenset(_TABLE_METADATA[_tid]["start_codons"]),
        "stop_codons": frozenset(_TABLE_METADATA[_tid]["stop_codons"]),
    }

# Convenience alias — the most commonly used table
STANDARD_CODON_TABLE: dict[str, str] = GENETIC_CODES[1]["codons"]

# Fast reverse-complement via translation table (ACGT only — the codon table
# does not contain ambiguous bases, so revcomp of non-ACGT is irrelevant here).
_REVCOMP_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")

# IUPAC ambiguity expansion — maps each ambiguous base to its possible ACGT bases.
# Used by _resolve_ambiguous_codon() to match BioPython's ambiguous codon handling.
_IUPAC_EXPAND: dict[str, list[str]] = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "S": ["C", "G"],
    "W": ["A", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}


def _resolve_ambiguous_codon(codon: str, table: dict[str, str]) -> str:
    """Translate a codon that may contain IUPAC ambiguity codes.

    Expands ambiguous bases to all possible concrete codons. If every expansion
    translates to the same amino acid, returns that amino acid. Otherwise
    returns 'X'. This matches Bio.Seq.translate behavior for ambiguous codons.
    """
    parts: list[list[str]] = [_IUPAC_EXPAND.get(b, [b]) for b in codon]
    # Cartesian product — for 'ATR': ['A'], ['T'], ['A','G'] → ATA, ATG
    from itertools import product as _product

    aa_set: set[str] = set()
    for concrete in _product(*parts):
        aa_set.add(table.get("".join(concrete), "X"))
        if len(aa_set) > 1:
            return "X"
    return aa_set.pop() if aa_set else "X"


def _reverse_complement(dna: str) -> str:
    """Reverse-complement a DNA string. Pure string operation, no allocation beyond the result."""
    return dna.translate(_REVCOMP_TRANS)[::-1]


# ---------------------------------------------------------------------------
# Translation result
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TranslationResult:
    """Output of translate().

    Warnings are surfaced so callers can decide whether partial codons or
    ambiguous bases are acceptable in their context.
    """

    protein: str
    warnings: list[str] = field(default_factory=list)
    table_id: int = 1
    stop_positions: list[int] = field(default_factory=list)


# ---------------------------------------------------------------------------
# translate()
# ---------------------------------------------------------------------------


def translate(
    dna: str,
    table_id: int = 1,
    frame: int = 0,
    strand: int = 1,
    to_stop: bool = False,
) -> TranslationResult:
    """Translate a DNA sequence to protein.

    Args:
        dna: DNA sequence string (case-insensitive; U is treated as T).
        table_id: NCBI genetic code table number (default 1, Standard).
        frame: Reading frame offset 0/1/2 (0 = first base of input).
        strand: +1 for forward strand, -1 for reverse complement.
        to_stop: If True, stop translation at the first stop codon (inclusive).

    Returns:
        TranslationResult with protein string, warnings, table provenance, and
        stop-codon positions.
    """
    table = GENETIC_CODES[table_id]["codons"]
    stop_codons = GENETIC_CODES[table_id]["stop_codons"]
    warnings: list[str] = []

    # Normalise: uppercase, T for U (so RNA sequences "just work")
    seq = dna.upper().replace("U", "T")

    # Strand
    if strand == -1:
        seq = _reverse_complement(seq)

    # Frame offset
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1, or 2; got {frame}")
    seq = seq[frame:]

    # Detect partial trailing codon
    remainder = len(seq) % 3
    if remainder != 0:
        warnings.append(
            f"Sequence length after frame offset is not a multiple of 3 "
            f"({len(seq)} nt, remainder {remainder}); trailing partial codon dropped."
        )
        seq = seq[: len(seq) - remainder]

    # Translate
    amino_acids: list[str] = []
    stop_positions: list[int] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if set(codon).issubset({"A", "C", "G", "T"}):
            aa = table.get(codon, "X")
        else:
            aa = _resolve_ambiguous_codon(codon, table)
            if aa == "X":
                warnings.append(
                    f"Codon {codon} at position {frame + i} contains ambiguous bases; "
                    f"translated as 'X'."
                )
        amino_acids.append(aa)
        if codon in stop_codons:
            stop_positions.append(len(amino_acids) - 1)
            if to_stop:
                warnings.append(
                    f"Translation halted at stop codon {codon} "
                    f"(AA position {stop_positions[-1]})."
                )
                amino_acids.pop()  # match BioPython: to_stop excludes the *
                break

    return TranslationResult(
        protein="".join(amino_acids),
        warnings=warnings,
        table_id=table_id,
        stop_positions=stop_positions,
    )


# ---------------------------------------------------------------------------
# FrameTranslation + six_frame_translate()
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FrameTranslation:
    """A single reading-frame translation.

    Frames 1-3 are forward (offsets 0, 1, 2); frames 4-6 are the same offsets
    on the reverse complement.
    """

    frame: int
    protein: str


def six_frame_translate(dna: str, table_id: int = 1) -> list[FrameTranslation]:
    """Translate all six reading frames.

    Returns six FrameTranslation objects (frames 1-6), each containing the
    protein that results from translating the given DNA in that frame.
    """
    frames: list[FrameTranslation] = []
    for strand in (1, -1):
        for offset in (0, 1, 2):
            frame_num = offset + 1 if strand == 1 else offset + 4
            result = translate(dna, table_id=table_id, frame=offset, strand=strand)
            frames.append(FrameTranslation(frame=frame_num, protein=result.protein))
    return frames


# ---------------------------------------------------------------------------
# ORF dataclass + find_orfs()
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ORF:
    """An open reading frame: start codon to first in-frame stop, inclusive.

    Positions are 0-indexed on the strand being searched.  ``protein`` includes
    the leading methionine and the trailing ``*`` (stop).
    """

    frame: int
    start: int
    end: int  # one past the stop codon
    protein: str

    @property
    def length_nt(self) -> int:
        return self.end - self.start

    @property
    def strand(self) -> str:
        return "+" if self.frame <= 3 else "-"


def _find_orfs_on_strand(
    seq: str,
    table_id: int,
    frame_base: int,
    start_codons: frozenset[str],
    stop_codons: frozenset[str],
) -> list[ORF]:
    """Find ORFs on one strand.

    Args:
        seq: Uppercase DNA, already revcomp'd if searching the reverse strand.
        table_id: Genetic code table for translating found ORFs.
        frame_base: 0 for forward (frames 1-3), 3 for reverse (frames 4-6).
        start_codons: Set of codons that initiate translation.
        stop_codons: Set of codons that terminate translation.
    """
    out: list[ORF] = []
    for offset in (0, 1, 2):
        # Extract codons for this reading frame
        codons = [seq[i : i + 3] for i in range(offset, len(seq) - 2, 3)]
        i = 0
        while i < len(codons):
            if codons[i] in start_codons:
                for j in range(i + 1, len(codons)):
                    if codons[j] in stop_codons:
                        start_nt = offset + i * 3
                        end_nt = offset + (j + 1) * 3
                        result = translate(seq[start_nt:end_nt], table_id=table_id)
                        out.append(
                            ORF(
                                frame=frame_base + offset + 1,
                                start=start_nt,
                                end=end_nt,
                                protein=result.protein,
                            )
                        )
                        break  # continue scanning from i+1 to catch nested ORFs
                else:
                    # No in-frame stop found — this start codon does not produce
                    # an ORF under our definition.
                    pass
            i += 1
    return out


def find_orfs(
    dna: str,
    table_id: int = 1,
    min_length: int = 100,
) -> list[ORF]:
    """Find all ORFs in all six reading frames.

    Args:
        dna: DNA sequence string.
        table_id: NCBI genetic code table number (default 1).
        min_length: Minimum ORF length in nucleotides including the stop codon
                    (default 100, roughly the smallest biologically interesting CDS).

    Returns:
        List of ORF objects. Nested ORFs (multiple starts sharing a stop) are
        all reported, matching the convention used by NCBI ORFfinder.
    """
    code = GENETIC_CODES[table_id]
    start_codons = code["start_codons"]
    stop_codons = code["stop_codons"]

    seq = dna.upper().replace("U", "T")

    forward = _find_orfs_on_strand(
        seq, table_id, frame_base=0, start_codons=start_codons, stop_codons=stop_codons
    )
    reverse = _find_orfs_on_strand(
        _reverse_complement(seq),
        table_id,
        frame_base=3,
        start_codons=start_codons,
        stop_codons=stop_codons,
    )

    return [orf for orf in forward + reverse if orf.length_nt >= min_length]

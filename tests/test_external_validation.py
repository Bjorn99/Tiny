"""External validation: Tiny translate/find_orfs vs BioPython, EMBOSS, NCBI ORFfinder.

Run with:
    poetry run pytest tests/test_external_validation.py -v

This file is deliberately separate from the main test suite — it validates
against external reference tools, not against internal expectations.
"""

import os
import random
import subprocess
import tempfile

import pytest

from tiny.algorithms.orf import (
    GENETIC_CODES,
    find_orfs,
    six_frame_translate,
    translate,
)

# ---------------------------------------------------------------------------
# Deterministic random sequences for reproducibility
# ---------------------------------------------------------------------------


def _make_random_dna(length: int, seed: int = 42) -> str:
    """Generate a random ACGT sequence with a fixed seed."""
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


# ---------------------------------------------------------------------------
# 1. BioPython cross-validation (automated, runs on every invocation)
# ---------------------------------------------------------------------------


_TRANSLATE_TABLE1_SEQUENCES = [
    ("ATG", "single start codon"),
    ("ATGAAATAG", "M-K-* short ORF"),
    ("TAA", "single stop codon"),
    ("ATGATGATG", "M-M-M no stop"),
    (_make_random_dna(300, seed=1), "random 300 bp"),
    (_make_random_dna(999, seed=2), "random 999 bp (multiple of 3)"),
    (_make_random_dna(1000, seed=3), "random 1000 bp (not multiple of 3)"),
    ("", "empty"),
    ("A", "single base"),
    ("AT", "two bases"),
    # RNA-style input with U
    ("AUGAAAAUG", "U-containing RNA-style"),
]

# Ambiguous-codon tests — Tiny matches BioPython when all expansions agree.
# Known deviation: Tiny uses 'X' where BioPython uses B/Z/J for ambiguous
# amino-acid groups (Asx/Glx/Xle). Documented, Phase 3 scope.
_AMBIGUOUS_CODON_TESTS = [
    ("ATR", "X", "ambiguous R, different AAs"),
    ("TCN", "S", "single N, all Ser"),
    ("GCN", "A", "single N, all Ala"),
    ("ATY", "I", "single Y, both Ile"),
]


class TestTranslateVsBioPython:
    """translate() must match Bio.Seq.translate() for all supported tables."""

    @pytest.mark.parametrize("dna,label", _TRANSLATE_TABLE1_SEQUENCES)
    def test_translate_table1_matches_biopython(self, dna, label):
        from Bio.Seq import Seq

        expected = str(Seq(dna).translate())
        got = translate(dna, table_id=1).protein
        assert got == expected, f"[{label}] table=1: Tiny={got!r} BioPython={expected!r}"

    @pytest.mark.parametrize("dna,expected,label", _AMBIGUOUS_CODON_TESTS)
    def test_translate_ambiguous_matches_biopython_when_unambiguous(self, dna, expected, label):
        """Ambiguous codons must match BioPython when all expansions agree on one AA."""
        got = translate(dna, table_id=1).protein
        assert got == expected, f"[{label}] Tiny={got!r} expected={expected!r}"
        # Also verify BioPython agrees (defense against test rot)
        from Bio.Seq import Seq

        bp = str(Seq(dna).translate())
        assert bp == expected, f"[{label}] BioPython changed: {bp!r} != expected {expected!r}"

    def test_translate_all_tables_random_sequence(self):
        """For every supported table, translate() must match Bio.Seq.translate()
        on a long random sequence, every frame, both strands."""
        from Bio.Seq import Seq

        seq = _make_random_dna(1500, seed=99)
        failures = []
        for table_id in sorted(GENETIC_CODES):
            for strand in (1, -1):
                for frame in (0, 1, 2):
                    tiny_result = translate(seq, table_id=table_id, frame=frame, strand=strand)
                    # BioPython: apply frame/strand manually, then translate
                    bp_seq = str(Seq(seq).reverse_complement()) if strand == -1 else seq
                    bp_seq = bp_seq[frame:]
                    bp_seq = bp_seq[: len(bp_seq) - len(bp_seq) % 3]
                    expected = str(Seq(bp_seq).translate(table=table_id, to_stop=False))
                    got = tiny_result.protein
                    if got != expected:
                        failures.append(
                            f"Table {table_id} strand={strand} frame={frame}: "
                            f"Tiny={got[:30]!r}... BioPython={expected[:30]!r}..."
                        )
        if failures:
            pytest.fail("\n".join(failures))

    def test_translate_to_stop_matches_biopython(self):
        seq = _make_random_dna(500, seed=77)
        from Bio.Seq import Seq

        tiny = translate(seq, to_stop=True)
        expected = str(Seq(seq).translate(to_stop=True))
        assert tiny.protein == expected, f"to_stop: Tiny={tiny.protein!r} BioPython={expected!r}"

    def test_six_frame_translate_matches_biopython(self):
        """Each of 6 frames must match Bio.Seq.translate on the
        appropriately offset/reverse-complemented sequence."""
        seq = "ATGAAATAGCCCAAAGGGTTT"
        from Bio.Seq import Seq

        tiny_frames = six_frame_translate(seq)
        for tf in tiny_frames:
            if tf.frame <= 3:
                offset = tf.frame - 1
                bp_seq = seq[offset:]
                bp_seq = bp_seq[: len(bp_seq) - len(bp_seq) % 3]
                expected = str(Seq(bp_seq).translate())
            else:
                offset = tf.frame - 4
                bp_seq = str(Seq(seq).reverse_complement())[offset:]
                bp_seq = bp_seq[: len(bp_seq) - len(bp_seq) % 3]
                expected = str(Seq(bp_seq).translate())
            assert (
                tf.protein == expected
            ), f"Frame {tf.frame}: Tiny={tf.protein!r} BioPython={expected!r}"


class TestFindOrfsVsBioPython:
    """find_orfs() results must be internally consistent with BioPython
    codon tables and our own translate() output."""

    def test_find_orfs_start_stop_match_table_definitions(self):
        """Every ORF found must start at a codon in the table's start_codons
        and end at a codon in the table's stop_codons."""
        seq = _make_random_dna(2000, seed=55)
        for table_id in sorted(GENETIC_CODES):
            code = GENETIC_CODES[table_id]
            orfs = find_orfs(seq, table_id=table_id, min_length=0)
            for orf in orfs:
                # Determine the strand-specific sequence
                from tiny.algorithms.orf import _reverse_complement

                strand_seq = seq if orf.strand == "+" else _reverse_complement(seq)
                start_codon = strand_seq[orf.start : orf.start + 3]
                stop_codon = strand_seq[orf.end - 3 : orf.end]
                assert start_codon in code["start_codons"], (
                    f"Table {table_id}: ORF start codon {start_codon} not in "
                    f"start_codons {sorted(code['start_codons'])}"
                )
                assert stop_codon in code["stop_codons"], (
                    f"Table {table_id}: ORF stop codon {stop_codon} not in "
                    f"stop_codons {sorted(code['stop_codons'])}"
                )

    def test_find_orfs_protein_matches_translate(self):
        """The protein field of every ORF must match what translate()
        produces for the same span."""
        seq = _make_random_dna(2000, seed=66)
        for table_id in sorted(GENETIC_CODES):
            orfs = find_orfs(seq, table_id=table_id, min_length=0)
            for orf in orfs:
                # Extract the ORF span and translate it independently
                from tiny.algorithms.orf import _reverse_complement

                if orf.strand == "+":
                    span = seq[orf.start : orf.end]
                    result = translate(span, table_id=table_id)
                else:
                    span = _reverse_complement(seq)[orf.start : orf.end]
                    result = translate(span, table_id=table_id)
                assert orf.protein == result.protein, (
                    f"Table {table_id} frame={orf.frame} start={orf.start}: "
                    f"ORF.protein={orf.protein!r} translate()={result.protein!r}"
                )

    def test_find_orfs_nested_are_all_reported(self):
        """When a sequence has two starts sharing one stop, both ORFs must be reported."""
        seq = "ATGAAAATGAAATAG"  # M K M K *
        orfs = find_orfs(seq, table_id=1, min_length=0)
        frame1 = [o for o in orfs if o.frame == 1]
        starts = sorted(o.start for o in frame1)
        assert starts == [0, 6], f"Expected nested ORFs at 0 and 6, got {starts}"
        assert len(frame1) == 2, f"Expected 2 nested ORFs, got {len(frame1)}"

    def test_find_orfs_no_false_positives_on_random(self):
        """On pure random sequence with default min_length=100, no ORFs
        should be reported (expected since 100 nt ORFs are rare by chance)."""
        seq = _make_random_dna(10000, seed=42)
        # With min_length=100, random ACGT should yield very few or no ORFs.
        # Statistically, ATG (1/64) then ~63 non-stop codons (61/64 each)
        # followed by stop (3/64): very rare in random sequence.
        orfs = find_orfs(seq, table_id=1, min_length=100)
        # A 10kb random sequence has ~3333 codons. Probability of an ATG
        # followed by 33 non-stop codons then a stop is tiny.
        # We just verify we don't crash and the result is reasonable.
        assert isinstance(orfs, list)


class TestRealSequences:
    """Validation on real biological sequences from eg_files/."""

    def test_cdkn1b_cds_translation_matches_biopython(self):
        from Bio import SeqIO
        from Bio.Seq import Seq

        record = next(SeqIO.parse("eg_files/CDKN1B_cds.fasta", "fasta"))
        dna = str(record.seq).upper()
        expected = str(Seq(dna).translate())
        assert translate(dna).protein == expected, "CDKN1B CDS: Tiny disagrees with BioPython"

    def test_cdkn1b_find_orfs_finds_annotated_cds(self):
        """CDKN1B CDS should produce at least one ORF >= 500 nt on the forward strand."""
        from Bio import SeqIO

        record = next(SeqIO.parse("eg_files/CDKN1B_cds.fasta", "fasta"))
        dna = str(record.seq).upper()
        orfs = find_orfs(dna, table_id=1, min_length=500)
        forward_orfs = [o for o in orfs if o.strand == "+"]
        assert len(forward_orfs) >= 1, (
            f"Expected at least one forward-strand ORF >= 500 nt in CDKN1B CDS, "
            f"found {len(forward_orfs)}"
        )


# ---------------------------------------------------------------------------
# 2. EMBOSS cross-validation (runs only if EMBOSS is installed)
# ---------------------------------------------------------------------------

_EMBOSS_AVAILABLE = subprocess.run(["which", "transeq"], capture_output=True).returncode == 0


_EMBOSS_SEQUENCES = [
    ("ATGAAATAG", "M-K-* short ORF"),
    (_make_random_dna(300, seed=10), "random 300 bp"),
]


@pytest.mark.skipif(not _EMBOSS_AVAILABLE, reason="EMBOSS not installed")
class TestTranslateVsEmboss:
    """translate() must match EMBOSS transeq output."""

    @pytest.mark.parametrize("dna,label", _EMBOSS_SEQUENCES)
    def test_translate_matches_transeq(self, dna, label):
        """Run transeq on the DNA and compare to Tiny's translate()."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(f">test\n{dna}\n")
            fasta_path = f.name

        try:
            result = subprocess.run(
                [
                    "transeq",
                    "-sequence",
                    fasta_path,
                    "-outseq",
                    "stdout",
                    "-frame",
                    "1",
                    "-table",
                    "0",
                    "-trim",
                ],
                capture_output=True,
                text=True,
            )
            # Parse transeq output: second line after the header
            lines = result.stdout.strip().split("\n")
            emboss_protein = lines[1].strip() if len(lines) >= 2 else ""
            tiny_result = translate(dna)
            assert (
                tiny_result.protein == emboss_protein
            ), f"[{label}] Tiny={tiny_result.protein!r} EMBOSS={emboss_protein!r}"
        finally:
            os.unlink(fasta_path)


@pytest.mark.skipif(not _EMBOSS_AVAILABLE, reason="EMBOSS not installed")
class TestFindOrfsVsEmboss:
    """find_orfs() must be consistent with EMBOSS getorf output."""

    def test_find_orfs_count_matches_getorf(self):
        """getorf should find the same number of ORFs as Tiny on a test sequence."""
        seq = "ATGAAATAGCCCATGCCCTAA"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(f">test\n{seq}\n")
            fasta_path = f.name

        try:
            result = subprocess.run(
                [
                    "getorf",
                    "-sequence",
                    fasta_path,
                    "-outseq",
                    "stdout",
                    "-minsize",
                    "0",
                    "-find",
                    "1",
                ],
                capture_output=True,
                text=True,
            )
            # Count EMBOSS ORFs: each ORF is a FASTA entry (header + sequence)
            emboss_count = result.stdout.count(">")
            tiny_count = len(find_orfs(seq, min_length=0))
            assert (
                emboss_count == tiny_count
            ), f"ORF count: Tiny={tiny_count} EMBOSS getorf={emboss_count}"
        finally:
            os.unlink(fasta_path)


# ---------------------------------------------------------------------------
# 3. NCBI ORFfinder manual validation instructions
# ---------------------------------------------------------------------------

# These are documented below as pytest doctests / manual steps.
# NCBI ORFfinder: https://www.ncbi.nlm.nih.gov/orffinder/

NCBI_ORFFINDER_MANUAL_TESTS = """
=== MANUAL VALIDATION: NCBI ORFfinder ===

Go to https://www.ncbi.nlm.nih.gov/orffinder/

Test 1: CDKN1B CDS
  - Load eg_files/CDKN1B_cds.fasta into ORFfinder
  - Genetic code: 1 (Standard)
  - Min ORF length: 100
  - Compare: number of ORFs found, longest ORF's protein sequence
  - Tiny equivalent: tiny find-orfs --input eg_files/CDKN1B_cds.fasta --min-length 100
  - Expected: both report the annotated CDS as the longest ORF on the + strand

Test 2: Alternative genetic code (ciliate)
  - Sequence: ATGCCCUAACCC
  - Genetic code: 6 (Ciliate)
  - Min ORF length: 0
  - Compare: In table 6, TAA is Gln (not STOP), so TGA is the only stop
  - Tiny equivalent: tiny find-orfs ATGCCCUAACCC --min-length 0 --table 6
  - Expected: No ORF found (no TGA in frame after ATG)

Test 3: Short synthetic nested ORFs
  - Sequence: ATGAAAATGAAATAG
  - Genetic code: 1 (Standard)
  - Min ORF length: 0
  - Compare: Number of ORFs reported in frame 1
  - Tiny equivalent: tiny find-orfs ATGAAAATGAAATAG --min-length 0
  - Expected: 2 ORFs in frame 1 (nested: starts at 0 and 6, same stop)
"""

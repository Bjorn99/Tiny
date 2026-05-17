import pytest

from tiny.algorithms.orf import (
    GENETIC_CODES,
    ORF,
    TranslationResult,
    find_orfs,
    six_frame_translate,
    translate,
)

# ---------------------------------------------------------------------------
# Codon table registry
# ---------------------------------------------------------------------------


def test_all_seven_tables_present():
    assert set(GENETIC_CODES) == {1, 2, 3, 4, 5, 6, 11}


def test_every_table_has_name_codons_start_stop():
    for tid, table in GENETIC_CODES.items():
        assert "name" in table
        assert "codons" in table
        assert "start_codons" in table
        assert "stop_codons" in table
        assert len(table["codons"]) == 64, f"Table {tid} missing codons"
        assert len(table["start_codons"]) >= 1, f"Table {tid} missing start codons"
        assert len(table["stop_codons"]) >= 1, f"Table {tid} missing stop codons"


def test_all_tables_match_biopython():
    """Cross-validate every codon, start, and stop against BioPython."""
    from Bio.Data import CodonTable

    for tid in sorted(GENETIC_CODES):
        t_bp = CodonTable.unambiguous_dna_by_id[tid]
        t_tiny = GENETIC_CODES[tid]
        for codon in t_bp.forward_table:
            assert (
                t_tiny["codons"][codon] == t_bp.forward_table[codon]
            ), f"Table {tid} codon {codon} mismatch"
        assert t_tiny["start_codons"] == frozenset(
            t_bp.start_codons
        ), f"Table {tid} start codons mismatch"
        assert t_tiny["stop_codons"] == frozenset(
            t_bp.stop_codons
        ), f"Table {tid} stop codons mismatch"


# ---------------------------------------------------------------------------
# translate() — happy path
# ---------------------------------------------------------------------------


def test_translate_start_codon():
    assert translate("ATG").protein == "M"


def test_translate_two_codons():
    assert translate("ATGAAA").protein == "MK"


def test_translate_stop_codons():
    assert translate("TAA").protein == "*"
    assert translate("TAG").protein == "*"
    assert translate("TGA").protein == "*"


def test_translate_full_short_orf():
    assert translate("ATGAAATAA").protein == "MK*"


def test_translate_is_case_insensitive():
    assert translate("atgaaa").protein == "MK"


def test_translate_u_is_treated_as_t():
    """RNA sequences with U should translate identically to DNA with T."""
    assert translate("AUGAAA").protein == "MK"


# ---------------------------------------------------------------------------
# translate() — edge cases
# ---------------------------------------------------------------------------


def test_translate_drops_partial_trailing_codon():
    result = translate("ATGAA")  # ATG + partial AA
    assert result.protein == "M"
    assert any("partial codon" in w.lower() for w in result.warnings)


def test_translate_empty_string():
    result = translate("")
    assert result.protein == ""


def test_translate_ambiguous_codon_returns_x():
    result = translate("ATN")
    assert result.protein == "X"
    assert any("ambiguous" in w.lower() for w in result.warnings)


def test_translate_mixed_known_and_ambiguous():
    result = translate("ATGATNAAA")
    assert result.protein == "MXK"


def test_translate_length_one_no_codons():
    result = translate("A")
    assert result.protein == ""
    assert any("partial codon" in w.lower() for w in result.warnings)


def test_translate_length_two_no_codons():
    result = translate("AT")
    assert result.protein == ""
    assert any("partial codon" in w.lower() for w in result.warnings)


# ---------------------------------------------------------------------------
# translate() — frame offset
# ---------------------------------------------------------------------------


def test_translate_frame_1_skips_first_base():
    # G ATG AAA TAG  → frame 1 drops the leading G
    result = translate("GATGAAATAG", frame=1)
    assert result.protein == "MK*"


def test_translate_frame_2_skips_two_bases():
    result = translate("GGATGAAATAG", frame=2)
    assert result.protein == "MK*"


def test_translate_invalid_frame_raises():
    with pytest.raises(ValueError, match="frame must be"):
        translate("ATG", frame=3)


# ---------------------------------------------------------------------------
# translate() — strand
# ---------------------------------------------------------------------------


def test_translate_reverse_strand():
    # CTATTTCAT → revcomp = ATGAAATAG → MK*
    result = translate("CTATTTCAT", strand=-1)
    assert result.protein == "MK*"


# ---------------------------------------------------------------------------
# translate() — to_stop
# ---------------------------------------------------------------------------


def test_translate_to_stop_halts_at_first_stop():
    # ATG GGG TAA ATG AAA  → should stop at TAA, give MG (stop excluded per BioPython convention)
    result = translate("ATGGGGTAAATGAAA", to_stop=True)
    assert result.protein == "MG"


def test_translate_to_stop_no_stop_codon_goes_to_end():
    result = translate("ATGGGGCCC", to_stop=True)
    assert result.protein == "MGP"


# ---------------------------------------------------------------------------
# translate() — alternative genetic codes
# ---------------------------------------------------------------------------


def test_translate_table_2_vertebrate_mitochondrial():
    """AGA/AGG are STOP in vertebrate mt, not Arg. TGA is Trp, not STOP."""
    # ATG AGA AGG TGA → table 1: MRR*, table 2: M**W
    r1 = translate("ATGAGAAGGTGA", table_id=1)
    r2 = translate("ATGAGAAGGTGA", table_id=2)
    assert r1.protein == "MRR*"
    assert r2.protein == "M**W"


def test_translate_table_5_invertebrate_mitochondrial():
    """AGA/AGG are Ser in invertebrate mt."""
    r1 = translate("ATGAGAAGGTGA", table_id=1)
    r5 = translate("ATGAGAAGGTGA", table_id=5)
    assert r1.protein == "MRR*"
    assert r5.protein == "MSSW"


def test_translate_table_6_ciliate():
    """TAA/TAG are Gln in ciliate nuclear code, not STOP."""
    # ATG TAA TAG TGA → table 1: M*** (all three are stops), table 6: MQQ*
    r1 = translate("ATGTAATAGTGA", table_id=1)
    r6 = translate("ATGTAATAGTGA", table_id=6)
    assert r1.protein == "M***"
    assert r6.protein == "MQQ*"


# ---------------------------------------------------------------------------
# translate() — TranslationResult fields
# ---------------------------------------------------------------------------


def test_translation_result_is_dataclass():
    r = translate("ATGAAATAG")
    assert isinstance(r, TranslationResult)
    assert r.table_id == 1
    assert r.protein == "MK*"
    assert isinstance(r.warnings, list)
    assert isinstance(r.stop_positions, list)


def test_translate_records_stop_positions():
    r = translate("ATGAAATAG")
    assert r.stop_positions == [2]  # 0-indexed AA position of stop


def test_translate_default_table_is_standard():
    r = translate("ATG")
    assert r.table_id == 1


# ---------------------------------------------------------------------------
# six_frame_translate()
# ---------------------------------------------------------------------------


def test_six_frame_translate_returns_six_frames():
    result = six_frame_translate("ATGAAATAGCCCAAAGGG")
    assert [f.frame for f in result] == [1, 2, 3, 4, 5, 6]


def test_six_frame_translate_frame_1_matches_translate():
    seq = "ATGAAATAG"
    frames = six_frame_translate(seq)
    frame1 = next(f for f in frames if f.frame == 1)
    assert frame1.protein == translate(seq).protein


def test_six_frame_translate_reverse_frames_use_reverse_complement():
    # revcomp("CTATTTCAT") = "ATGAAATAG" → MK*
    seq = "CTATTTCAT"
    frames = six_frame_translate(seq)
    frame4 = next(f for f in frames if f.frame == 4)
    assert frame4.protein == "MK*"


def test_six_frame_translate_respects_table_id():
    seq = "ATGAGAAGGTGA"
    frames = six_frame_translate(seq, table_id=2)
    frame1 = next(f for f in frames if f.frame == 1)
    assert frame1.protein == "M**W"


# ---------------------------------------------------------------------------
# find_orfs() — basic
# ---------------------------------------------------------------------------


def test_find_orfs_finds_single_forward_orf():
    orfs = find_orfs("ATGAAATAG", min_length=0)
    assert len(orfs) == 1
    orf = orfs[0]
    assert isinstance(orf, ORF)
    assert orf.frame == 1
    assert orf.start == 0
    assert orf.end == 9
    assert orf.protein == "MK*"


def test_find_orfs_returns_empty_when_no_start():
    assert find_orfs("AAACCCGGG", min_length=0) == []


def test_find_orfs_returns_empty_when_start_has_no_in_frame_stop():
    # ATG AAA AAA AAA — no in-frame stop
    assert find_orfs("ATGAAAAAAAAA", min_length=0) == []


# ---------------------------------------------------------------------------
# find_orfs() — min_length filter
# ---------------------------------------------------------------------------


def test_find_orfs_min_length_excludes_short_orf():
    assert find_orfs("ATGAAATAG", min_length=300) == []


def test_find_orfs_min_length_includes_at_threshold():
    orfs = find_orfs("ATGAAATAG", min_length=9)
    assert len(orfs) == 1


def test_find_orfs_default_min_length_is_100():
    assert find_orfs("ATGAAATAG") == []


# ---------------------------------------------------------------------------
# find_orfs() — reverse strand
# ---------------------------------------------------------------------------


def test_find_orfs_detects_reverse_strand_orf():
    seq = "CTATTTCAT"  # revcomp = ATGAAATAG
    orfs = find_orfs(seq, min_length=0)
    assert len(orfs) == 1
    orf = orfs[0]
    assert orf.strand == "-"
    assert orf.frame in (4, 5, 6)
    assert orf.protein == "MK*"


def test_find_orfs_finds_orfs_on_both_strands():
    # Build a sequence with ORFs on both strands
    forward_orf = "ATGAAATAG"  # M K *
    # revcomp of ATGCCCTAA = TTAGGGCTAG (no ATG on forward)
    # We embed both: forward ATGAAA... and reverse complement ATGCCC...
    spacer = "CCCCCC"
    # Both strands should yield ORFs in this construct
    seq2 = forward_orf + spacer + "CTATTTCAT"
    orfs = find_orfs(seq2, min_length=0)
    strands = sorted({o.strand for o in orfs})
    assert "+" in strands
    # The reverse-strand ORF should exist
    rev_orfs = [o for o in orfs if o.strand == "-"]
    assert len(rev_orfs) >= 1


# ---------------------------------------------------------------------------
# find_orfs() — multiple and nested ORFs
# ---------------------------------------------------------------------------


def test_find_orfs_two_orfs_same_frame():
    seq = "ATGAAATAG" + "CCCCCC" + "ATGCCCTAA"
    orfs_frame1 = [o for o in find_orfs(seq, min_length=0) if o.frame == 1]
    assert len(orfs_frame1) == 2
    proteins = sorted(o.protein for o in orfs_frame1)
    assert proteins == ["MK*", "MP*"]


def test_find_orfs_nested_orfs_share_stop():
    # ATG AAA ATG AAA TAG → outer MKMK*, inner MK*
    seq = "ATGAAAATGAAATAG"
    orfs = [o for o in find_orfs(seq, min_length=0) if o.frame == 1]
    assert len(orfs) == 2
    starts = sorted(o.start for o in orfs)
    assert starts == [0, 6]
    assert {o.end for o in orfs} == {15}


# ---------------------------------------------------------------------------
# find_orfs() — alternative tables
# ---------------------------------------------------------------------------


def test_find_orfs_table_6_ciliate_uses_only_tga_stop():
    """In ciliate code only TGA is STOP (TAA/TAG are Gln)."""
    # ATG CCC TAA → table 1 has a short ORF; table 6 has no stop
    seq = "ATGCCCTAA"
    orfs_std = find_orfs(seq, table_id=1, min_length=0)
    orfs_cil = find_orfs(seq, table_id=6, min_length=0)
    assert len(orfs_std) == 1  # TAA is stop
    assert len(orfs_cil) == 0  # TAA is Gln, no stop in frame


# ---------------------------------------------------------------------------
# ORF dataclass properties
# ---------------------------------------------------------------------------


def test_orf_length_nt():
    orf = ORF(frame=1, start=0, end=9, protein="MK*")
    assert orf.length_nt == 9


def test_orf_strand_forward():
    for f in (1, 2, 3):
        orf = ORF(frame=f, start=0, end=9, protein="M*")
        assert orf.strand == "+"


def test_orf_strand_reverse():
    for f in (4, 5, 6):
        orf = ORF(frame=f, start=0, end=9, protein="M*")
        assert orf.strand == "-"


def test_orf_is_frozen():
    orf = ORF(frame=1, start=0, end=9, protein="MK*")
    with pytest.raises(AttributeError):
        orf.frame = 2  # type: ignore[misc]


# ---------------------------------------------------------------------------
# BioPython cross-validation on real CDS
# ---------------------------------------------------------------------------


def test_translate_matches_biopython_on_cdkn1b_cds():
    """Our translate() must agree with Bio.Seq.translate() on a known CDS."""
    from Bio import SeqIO
    from Bio.Seq import Seq

    record = next(SeqIO.parse("eg_files/CDKN1B_cds.fasta", "fasta"))
    dna = str(record.seq).upper()
    expected = str(Seq(dna).translate())
    assert translate(dna).protein == expected

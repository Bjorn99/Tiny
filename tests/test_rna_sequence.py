import pytest

from tiny.core.errors import InvalidSequenceError
from tiny.core.sequence import RNASequence


@pytest.mark.quick
def test_rna_accepts_u_not_t():
    seq = RNASequence("AUCG")
    assert str(seq) == "AUCG"


@pytest.mark.quick
def test_rna_rejects_t():
    with pytest.raises(InvalidSequenceError):
        RNASequence("ATCG")


@pytest.mark.quick
def test_rna_gc_content_matches_dna_equivalent():
    """RNA AUCG should have same GC content as DNA ATCG."""
    rna = RNASequence("AUCG")
    assert rna.gc_content == 50.0


@pytest.mark.quick
def test_rna_complement_uses_u():
    rna = RNASequence("AUCG")
    assert str(rna.complement) == "UAGC"

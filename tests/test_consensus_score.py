import pytest

from tiny.core.advanced import AdvancedAnalysis


@pytest.mark.quick
def test_uniform_motif_has_low_information():
    """A motif of all one base scores high (very conserved)."""
    score_all_A = AdvancedAnalysis._calculate_consensus_score("AAAA")
    assert score_all_A == pytest.approx(2.0)  # max bits per position


@pytest.mark.quick
def test_consensus_score_distinguishes_motifs():
    """Different motifs produce different scores — the bug returned 0.25 for everything."""
    score1 = AdvancedAnalysis._calculate_consensus_score("ATCG")
    score2 = AdvancedAnalysis._calculate_consensus_score("AAAA")
    assert score1 != score2


@pytest.mark.quick
def test_score_is_bounded_zero_to_two():
    """Information content per base is bounded [0, 2] bits."""
    for motif in ["A", "AT", "ATCG", "AAAATTTT"]:
        score = AdvancedAnalysis._calculate_consensus_score(motif)
        assert 0.0 <= score <= 2.0

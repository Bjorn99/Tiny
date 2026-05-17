from typer.testing import CliRunner

from tiny.cli import app

runner = CliRunner()


# ---------------------------------------------------------------------------
# translate command
# ---------------------------------------------------------------------------


def test_translate_inline_sequence_prints_protein():
    result = runner.invoke(app, ["translate", "ATGAAATAG"])
    assert result.exit_code == 0
    assert "MK*" in result.stdout


def test_translate_is_case_insensitive():
    result = runner.invoke(app, ["translate", "atgaaa"])
    assert result.exit_code == 0
    assert "MK" in result.stdout


def test_translate_from_file():
    result = runner.invoke(app, ["translate", "--input", "eg_files/CDKN1B_cds.fasta"])
    assert result.exit_code == 0
    assert "M" in result.stdout


def test_translate_invalid_input_exits_nonzero():
    result = runner.invoke(app, ["translate", "ATGZZZ"])
    assert result.exit_code != 0


def test_translate_missing_input_exits_2():
    result = runner.invoke(app, ["translate"])
    assert result.exit_code == 2


def test_translate_with_table_flag():
    """--table 2 should produce vertebrate mt translation (AGA -> STOP)."""
    # ATG AGA → standard: MR, vertebrate mt: M*
    result = runner.invoke(app, ["translate", "ATGAGA", "--table", "2"])
    assert result.exit_code == 0
    assert "M*" in result.stdout


def test_translate_with_strand_flag():
    """--strand -1 should translate reverse complement."""
    # CTATTTCAT revcomp = ATGAAATAG → MK*
    result = runner.invoke(app, ["translate", "CTATTTCAT", "--strand", "-1"])
    assert result.exit_code == 0
    assert "MK*" in result.stdout


def test_translate_with_to_stop_flag():
    result = runner.invoke(app, ["translate", "ATGGGGTAAATGAAA", "--to-stop"])
    assert result.exit_code == 0
    assert "MG*" in result.stdout


def test_translate_explain_prints_walkthrough():
    result = runner.invoke(app, ["translate", "ATG", "--explain"])
    assert result.exit_code == 0
    text = result.stdout.lower()
    assert "codon" in text
    assert "genetic code" in text or "reading frame" in text


# ---------------------------------------------------------------------------
# find-orfs command
# ---------------------------------------------------------------------------


def test_find_orfs_inline_sequence_lists_orfs():
    result = runner.invoke(app, ["find-orfs", "ATGAAATAG", "--min-length", "0"])
    assert result.exit_code == 0
    assert "MK*" in result.stdout
    assert "1" in result.stdout  # frame column


def test_find_orfs_respects_min_length():
    result = runner.invoke(app, ["find-orfs", "ATGAAATAG"])
    assert result.exit_code == 0
    assert "No ORFs found" in result.stdout


def test_find_orfs_from_file():
    result = runner.invoke(
        app, ["find-orfs", "--input", "eg_files/CDKN1B_cds.fasta", "--min-length", "100"]
    )
    assert result.exit_code == 0
    assert "M" in result.stdout


def test_find_orfs_explain_prints_walkthrough():
    result = runner.invoke(app, ["find-orfs", "ATGAAATAG", "--min-length", "0", "--explain"])
    assert result.exit_code == 0
    text = result.stdout.lower()
    assert "start codon" in text
    assert "stop codon" in text


def test_find_orfs_missing_input_exits_2():
    result = runner.invoke(app, ["find-orfs"])
    assert result.exit_code == 2


def test_find_orfs_with_table_flag():
    """Table 6 (ciliate) uses only TGA as stop — ORF should differ."""
    result = runner.invoke(app, ["find-orfs", "ATGCCCTAA", "--min-length", "0", "--table", "6"])
    assert result.exit_code == 0
    assert "No ORFs found" in result.stdout  # TAA is Gln, not stop

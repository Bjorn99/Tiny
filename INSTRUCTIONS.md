# Tiny — Install, Run, Troubleshoot

A practical guide to getting Tiny working on your machine and handling the errors you'll actually hit.

---

## 1. Install

### 1.1 Prerequisites

- **Python 3.12, 3.13, or 3.14**
- **Poetry 2.x** (Python's dep manager)
- Git

Check what you have:

```bash
python3 --version
poetry --version
git --version
```

### 1.2 Install Poetry

| Platform | Command |
|---|---|
| Arch Linux | `sudo pacman -S python-poetry` |
| Debian/Ubuntu | `sudo apt install python3-poetry` *(may be old; prefer the installer)* |
| macOS | `brew install poetry` |
| Anywhere | `curl -sSL https://install.python-poetry.org \| python3 -` |

If you used the official installer, make sure `~/.local/bin` is on your `PATH`.

### 1.3 Clone and install

```bash
git clone https://github.com/Bjorn99/Tiny.git
cd Tiny

# Core install — no SAM/BAM. Works on Python 3.12, 3.13, 3.14.
poetry install

# OR: install with SAM/BAM support (pysam). Python 3.12 or 3.13 only —
# pysam has no wheels for 3.14 yet.
poetry install --extras sam
```

### 1.4 Verify

```bash
poetry run tiny --version       # → tiny 0.2.0
poetry run tiny analyze ATCG    # should print a small results table
poetry run pytest -q            # → 20 passed
```

---

## 2. Run

There are two ways to invoke `tiny`:

### 2.1 With `poetry run` (recommended for one-offs)

```bash
poetry run tiny analyze ATCGATCG
poetry run tiny --help
```

No setup needed; Poetry handles the venv.

### 2.2 Activate the venv (for many commands)

Poetry 2.x removed `poetry shell`. Use `poetry env activate` — it **prints** the activation command, you must run it.

**Fish:**
```fish
source (poetry env info --path)/bin/activate.fish
```

**Bash / Zsh:**
```bash
source "$(poetry env info --path)/bin/activate"
```

Once activated, just run `tiny ...` directly. `deactivate` exits.

---

## 3. Help

| You want to... | Command |
|---|---|
| See the version | `tiny --version` or `tiny -V` |
| Global help (list commands) | `tiny --help` |
| Per-command help | `tiny analyze --help`, `tiny align --help`, etc. |
| List supported file formats | `tiny supported-formats` |

### Available commands

- `tiny analyze` — sequence-level analysis (length, GC, MW, composition)
- `tiny compare` — pairwise comparison with mutation detection
- `tiny align` — global / local / semi-global alignment
- `tiny find-motifs` — fixed-length motif discovery
- `tiny find-regulatory` — TATA / GC / CAAT boxes + palindromes
- `tiny supported-formats` — format reference

See [Examples.md](Examples.md) for worked examples.

---

## 4. Errors and how to fix them

Tiny raises typed errors that explain the problem. The common ones:

### 4.1 `InvalidSequenceError`

```
Invalid DNASequence bases: X, Z. Allowed: A, B, C, D, G, H, K, M, N, R, S, T, V, W, Y
```

**Cause:** Your input has characters outside the IUPAC DNA alphabet.

**Fix:** Clean the input. For RNA (containing `U`), Tiny has a separate `RNASequence` class in `tiny.core.sequence` — but the CLI commands still target DNA. If you have RNA in a FASTA, convert `U → T` before passing it in:

```bash
sed 's/U/T/g; s/u/t/g' rna.fasta > as_dna.fasta
tiny analyze --input as_dna.fasta
```

### 4.2 `ResourceLimitError`

```
sequence length exceeds limit: got 25000, max allowed 10000.
Override with environment variables if needed.
```

**Cause:** Tiny caps per-sequence length at 10,000 bp by default to keep it laptop-friendly.

**Fix:** Bump the cap for that one call:

```bash
TINY_MAX_SEQUENCE=50000 poetry run tiny analyze --input long.fasta
```

Set it persistently in your shell rc if you always work with larger sequences. Be aware that motif finding and pairwise alignment scale poorly with length.

### 4.3 `OptionalDependencyError`

```
Optional dependency 'pysam' is not installed.
Install with: poetry install -E sam
```

**Cause:** You tried to read a SAM/BAM file without the `sam` extra installed.

**Fix:**

```bash
poetry install --extras sam
```

If that fails on Python 3.14, you'll need to drop to 3.13 — pysam has no 3.14 wheels yet. Easiest path:

```bash
poetry env use python3.13
poetry install --extras sam
```

### 4.4 `UnsupportedFormatError`

**Cause:** The file extension isn't one of `.fa`, `.fasta`, `.fq`, `.fastq`, `.gb`, `.gbk`, `.genbank`, `.embl`, `.sam`, `.bam`.

**Fix:** Either rename the file or convert it. `tiny supported-formats` lists what's accepted.

### 4.5 `poetry install` times out on PyPI

You may see `ConnectionError: Read timed out` mid-install. PyPI's CDN can be flaky.

**Fix:** Just re-run `poetry install`. Poetry skips already-downloaded wheels, so subsequent runs are usually faster. If it keeps failing, try `POETRY_REQUESTS_TIMEOUT=120 poetry install`.

### 4.6 `tiny --help` shows a Click traceback

If you see `TypeError: Parameter.make_metavar() missing 1 required positional argument: 'ctx'`, your `typer` is too old for the installed `click`. Phase 0 pinned `typer ^0.15.0` to fix this — make sure your `pyproject.toml` matches and re-run `poetry install`.

---

## 5. Development workflow

### 5.1 Tests

```bash
poetry run pytest -v            # full suite
poetry run pytest -m quick      # fast subset (used by pre-commit)
poetry run pytest -m "not requires_sam"  # skip SAM/BAM tests
```

### 5.2 Lint and format

```bash
poetry run ruff check .         # lint
poetry run ruff format .        # auto-format
poetry run ruff format --check . # verify formatting (CI does this)
```

### 5.3 Pre-commit hooks

Set up once:

```bash
poetry run pre-commit install
```

After that every `git commit` runs:
- `ruff` (lint + auto-fix)
- `ruff-format`
- Trailing whitespace, end-of-file, YAML/TOML checks
- `pytest -m quick` (must pass)

If a hook auto-fixes files, re-stage them and commit again. Don't bypass with `--no-verify`.

### 5.4 CI

The `.github/workflows/ci.yml` workflow runs on push and PR: a core matrix (no extras) on Python 3.12/3.13, plus a full matrix with the `sam` extra. CI does lint + format check + tests + build.

---

## 6. Where things live

| Path | What |
|---|---|
| `tiny/cli.py` | CLI entry point (lazy-loads core modules) |
| `tiny/core/sequence.py` | `DNASequence`, `RNASequence` |
| `tiny/core/analysis.py` | Basic analysis (GC, MW, compare) |
| `tiny/core/advanced.py` | Alignment, motif finding, regulatory elements |
| `tiny/core/formats.py` | FASTA/FASTQ/GenBank/EMBL/SAM-BAM readers |
| `tiny/core/errors.py` | Typed exception hierarchy |
| `tests/` | Pytest suite (mark `quick` for fast tests) |
| `docs/superpowers/` | Design spec + Phase plans |
| `eg_files/` | Example sequences for trying things out |

---

## 7. Reporting problems

If you hit something not covered here, open an issue with:
- Your Python version (`python3 --version`)
- Your Poetry version (`poetry --version`)
- The exact command you ran
- The full traceback

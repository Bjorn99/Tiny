# Tiny Revival — Design Spec

**Date:** 2026-05-13
**Status:** Draft for user review
**Author:** Brainstormed collaboratively with Claude

---

## Summary

Revive Tiny — an existing terminal-based DNA sequence analysis CLI — as part of a three-artifact portfolio for transitioning from wet-lab cancer genomics into computational biology. **This spec covers Tiny only.** The companion cancer-genomics analysis project and the blog are tracked separately (see "Related Artifacts" at the end).

Tiny's new identity: **a focused sketchbook of classical bioinformatics algorithms, each implemented from scratch with deep explanation, wrapped in a beautiful Rich-powered CLI.** Not a Swiss Army knife. Not a BioPython competitor. A learning artifact that demonstrates fundamentals.

---

## Goals

1. **Learning** — Implement classical bioinformatics algorithms from scratch with full understanding of biology, math, and code.
2. **Portfolio** — Produce a focused, well-documented project that signals comp-bio competence to potential collaborators or hiring managers.
3. **Sustainable** — Achievable at 5-7 hours/week over ~6 months alongside a wet-lab day job.
4. **Architecturally clean** — Strong, stable trunk; pluggable, decoupled leaves.

---

## Non-Goals

These are explicitly **out of scope** for Tiny (some belong to other artifacts):

- Becoming a general-purpose bioinformatics platform competing with BioPython, EMBOSS, or Galaxy.
- NGS pipeline features (FASTQ QC dashboards, variant calling, coverage analysis). These belong in the cancer-genomics analysis project, not in Tiny.
- Multiple sequence alignment, phylogenetics, secondary structure prediction. Tempting but out of scope.
- Production-scale optimization. Tiny prioritizes correctness and clarity over performance.
- New file format support beyond what already exists.

---

## Guiding Principles

### 1. Strong trunk, dynamic leaves
The user's chosen architectural philosophy. Operationalized as:

- `tiny/core/` (sequence, analysis, formats) is the **trunk**. Stable. Tests pin its behavior. Avoid modifying.
- Each new algorithm becomes its own module in `tiny/algorithms/<name>.py` — a **leaf**.
- New CLI commands plug in via Typer without modifying existing commands.
- Heavy dependencies (pysam, future ML libs) MUST be lazy-imported so they cannot break the trunk on `import tiny`.
- When tempted to "improve" trunk code while adding a leaf: stop. Trunk changes need their own commit with explicit reason.

**Validation rule:** `import tiny.core.sequence` and `tiny --help` must succeed even if pysam, BioPython, or other optional deps are broken/missing.

### 2. Every algorithm is also a lesson
For each new feature, three artifacts ship together:
- Clean from-scratch implementation in `tiny/algorithms/`
- CLI command exposing it via Typer
- Companion Jupyter notebook in `notebooks/` covering biology, math, and code

A `--explain` flag on each new command prints a condensed walk-through to the terminal.

### 3. Verify against BioPython, don't replace it
For algorithms BioPython implements (alignment, translation, etc.), keep the existing BioPython-backed command as the production version. Add a new from-scratch command (e.g., `align-dp` alongside `align`). The from-scratch version is for learning; the BioPython version is for actual use.

### 4. Standardize what's inconsistent
- The CLI uses `--input` for file input, but README and Examples.md mention `--fasta` in places. Standardize on `--input` everywhere; update docs to match.

### 5. Yagni ruthlessly
No abstractions for "future flexibility." No configurability beyond what an algorithm naturally requires. Comments only when the *why* is non-obvious. Notebook prose carries the explanatory load.

---

## Current State Assessment

**What's solid and stays:**
- Typer CLI structure (`tiny/cli.py`)
- Rich terminal output (tables, trees, progress bars, panels) — a genuine differentiator
- Module organization: `tiny/cli.py`, `tiny/core/`, `tiny/config/`
- File format support via `tiny/core/formats.py` (FASTA, FASTQ, GenBank, EMBL, SAM/BAM)
- Existing commands: `analyze`, `compare`, `align`, `find-motifs`, `find-regulatory`, `supported-formats`
- BioPython as foundation for file parsing and reference algorithms
- Test scaffold in `tests/`

**What's broken and needs fixing:**
1. **Install/run blocker (Phase 0 priority):** Tiny fails to start because `tiny/core/formats.py` does an eager top-level `import pysam`, and pysam's C extension is incompatible with the installed Python 3.14 (`AttributeError: module 'pysam.libcalignedsegment' has no attribute 'CMATCH'`). This breaks the trunk because of a leaf dependency.
2. **Consensus score bug** (`tiny/core/advanced.py::_calculate_consensus_score`): all bases weighted 0.25 → returns 0.25 for every motif. Meaningless metric. Fix or replace with information content.
3. **CLI flag inconsistency**: docs reference `--fasta`, CLI uses `--input`. Standardize on `--input`.

**What's overengineered and may be simplified:**
- The GenBank feature-tree rendering inside `analyze` in `cli.py` (~80 lines). Move to its own module or simplify. Not urgent.

---

## Target Architecture

```
Tiny/
├── tiny/
│   ├── __init__.py
│   ├── cli.py                    # Thin Typer entrypoint; one sub-command per feature
│   ├── config/
│   │   └── settings.py
│   ├── core/                     # TRUNK — stable, well-tested, rarely changes
│   │   ├── __init__.py
│   │   ├── sequence.py           # DNASequence (+ future RNASequence)
│   │   ├── analysis.py           # Basic analysis (GC, MW, composition)
│   │   ├── formats.py            # File I/O — pysam is now lazy-imported
│   │   └── utils.py
│   └── algorithms/               # LEAVES — each new algorithm is its own module
│       ├── __init__.py
│       ├── orf.py                # ORF finding + translation
│       ├── primer.py             # Primer Tm (nearest-neighbor)
│       ├── alignment_dp.py       # Hand-coded NW + SW
│       ├── motifs_pwm.py         # PWM motif scoring
│       └── bwt.py                # BWT + FM-index
├── notebooks/                    # Companion educational notebooks
│   ├── 01_orf_translation.ipynb
│   ├── 02_primer_tm.ipynb
│   ├── 03_pairwise_alignment_dp.ipynb
│   ├── 04_pwm_motifs.ipynb
│   └── 05_bwt_fm_index.ipynb
├── tests/                        # Pytest, with separate test files per module
│   ├── test_sequence.py
│   ├── test_analysis.py
│   ├── test_orf.py
│   ├── test_primer.py
│   ├── test_alignment_dp.py
│   ├── test_motifs_pwm.py
│   └── test_bwt.py
├── docs/
│   └── superpowers/specs/        # Design docs live here
├── eg_files/                     # Example sequence files
├── pyproject.toml
└── README.md
```

Two new top-level directories: `tiny/algorithms/` and `notebooks/`. Everything else either exists already or is a small fix.

---

## Phase-by-Phase Plan

Each phase is **roughly 2 weeks** at 5-7 hours/week. Reality compresses and expands. The phase ordering matters; the calendar doesn't.

### Phase 0 — Install fix, bug fix, RNA support (Week 1)

**Goal:** A working `tiny --help` and clean foundation for new work.

Concrete tasks:
- Make pysam a **lazy import** in `tiny/core/formats.py` — only imported when `read_sam_bam()` is called. If pysam is unavailable, the SAM/BAM handler raises a clear "install pysam to use this format" error. Trunk loads cleanly without it.
- Fix `_calculate_consensus_score` bug — replace with Shannon information content over base frequencies, or remove the field entirely from `MotifResult` if it's not used elsewhere.
- Add minimal `RNASequence` support — either an `is_rna` flag on `DNASequence` or a sister class. Trivial extension.
- Standardize CLI on `--input` everywhere. Update README.md and Examples.md to remove `--fasta` references.
- Verify `tiny` runs end-to-end on Python 3.14: `tiny analyze ATCG`, `tiny analyze --input eg_files/<some.fasta>`.

Validation:
- `poetry run tiny --help` succeeds.
- `poetry run pytest` passes.
- Removing/uninstalling pysam doesn't break `tiny --help` or `tiny analyze`.

### Phase 1 — ORF finding + translation (Weeks 2-3)

**Goal:** First "algorithm + notebook + CLI" cycle ships, establishing the template.

Concrete tasks:
- New module `tiny/algorithms/orf.py`:
  - Function: scan all 6 reading frames (3 forward, 3 reverse-complement).
  - Function: identify ORFs (start codon `ATG` to stop codon `TAA`/`TAG`/`TGA`).
  - Function: translate DNA → protein using the standard genetic code table.
- New CLI command `tiny find-orfs` and `tiny translate`.
- Add `--explain` flag.
- Notebook `notebooks/01_orf_translation.ipynb`:
  - The biology: what's a reading frame, why six of them, why stop codons matter.
  - The genetic code as a 64×1 lookup table.
  - Worked example with a known gene from `eg_files/`.
  - Comparison against `Bio.Seq.translate` as validation.
- Tests in `tests/test_orf.py` (frame extraction, ORF detection edge cases, translation correctness).
- First blog post drafted from the notebook (optional — can come later).

### Phase 2 — Primer Tm calculator (Weeks 4-5)

**Goal:** The most wet-lab-relevant algorithm. User's existing intuition makes this satisfying.

Concrete tasks:
- New module `tiny/algorithms/primer.py`:
  - Nearest-neighbor Tm calculation (SantaLucia 1998 parameters).
  - GC content sanity check.
  - Basic hairpin detection (self-reverse-complement matches).
  - Basic primer-dimer detection.
- New CLI command `tiny primer-tm <sequence>`.
- Notebook `notebooks/02_primer_tm.ipynb`:
  - The thermodynamics: ΔG, ΔH, ΔS, why Tm depends on neighboring bases.
  - The wet-lab perspective: when does the calculation fail in real PCR?
  - Validate against an online primer Tm calculator (e.g., IDT or BioPython's `Bio.SeqUtils.MeltingTemp`).
- Tests with primers of known Tm from the literature.

### Phase 3 — PWM-based motif finding (Weeks 6-7)

**Goal:** Statistical motif finding done right. Replaces (or supplements) the buggy consensus-score approach.

Concrete tasks:
- New module `tiny/algorithms/motifs_pwm.py`:
  - Build PWM from a set of aligned sequences.
  - Convert frequencies to log-odds scores (with pseudocounts).
  - Score a candidate sequence against a PWM.
  - Compute information content per position.
- Extend `tiny find-motifs` with `--method pwm` flag, alongside the existing sliding-window method.
- Notebook `notebooks/04_pwm_motifs.ipynb`:
  - Information theory primer: entropy, information content, log-odds.
  - PWMs vs simple frequency counting — what's the biological intuition?
  - Example: scan for TATA box using a real PWM (from JASPAR or hand-derived from known TATA sequences).

### Phase 4 — Hand-coded DP alignment (Weeks 8-10)

**Goal:** Implement Needleman-Wunsch and Smith-Waterman from scratch. Verify against BioPython.

Concrete tasks:
- New module `tiny/algorithms/alignment_dp.py`:
  - Needleman-Wunsch (global) with affine gap penalties.
  - Smith-Waterman (local).
  - Traceback function that reconstructs the alignment.
- New CLI command `tiny align-dp` (existing `tiny align` stays as the BioPython-backed version).
- `--explain` flag prints the DP matrix for short sequences.
- Notebook `notebooks/03_pairwise_alignment_dp.ipynb`:
  - The DP recurrence (with diagrams).
  - Visualize the DP matrix in Rich (heatmap-style table).
  - Why affine gaps matter biologically.
  - Cross-check: same inputs, same scores from `tiny align-dp` and `tiny align`.

### Phase 5 — BWT + FM-index (Weeks 11-14)

**Goal:** The crown jewel. Demonstrates understanding of modern NGS alignment internals.

Concrete tasks:
- New module `tiny/algorithms/bwt.py`:
  - Build BWT from a reference sequence (with sentinel character).
  - Build the FM-index data structures: C array, occurrence array, suffix array (sampled).
  - Backward search for exact match.
  - (Stretch: 1-mismatch search via seed-and-extend or backtracking.)
- New CLI command `tiny search <pattern> --reference <file>`.
- Notebook `notebooks/05_bwt_fm_index.ipynb`:
  - Build BWT step-by-step on a small string.
  - Why BWT enables fast string search.
  - How BWA uses this for read alignment.
  - Memory footprint vs naive search.
- This phase intentionally takes longer. It's the hardest and most rewarding.

### Phase 6 — Optional: NCBI/Entrez fetch (Week 15)

**Goal:** Practical connector to real biological databases.

Concrete tasks:
- New module `tiny/algorithms/fetch.py` (lazy-imports `Bio.Entrez`).
- New CLI command `tiny fetch <accession>` saves a FASTA or GenBank file.
- Notebook covers Entrez basics and rate limits.

This is "nice to have" — defer if time-constrained. Not algorithmically deep, but practical.

---

## Trunk Stability Contract

Once Phase 0 is complete, the trunk modules have stable interfaces:

- `DNASequence(seq: str)` — IUPAC validation, GC content, molecular weight, complement, reverse complement.
- `RNASequence(seq: str)` (if added) — equivalent API for RNA.
- `analyze_sequences(seqs) -> List[AnalysisResult]`
- `compare_sequences(s1, s2) -> ComparisonResult`
- `FormatHandler.read_file(path) -> List[SequenceRecord]` — dispatches on extension; SAM/BAM raises a clear error if pysam is missing.

**Rule:** Changing any of these signatures requires updating callers AND tests in the same commit. Treat the trunk as a published API.

---

## Testing Strategy

- **Trunk modules:** comprehensive unit tests pinning current behavior.
- **Each new algorithm:** unit tests covering happy path + 2-3 edge cases.
- **Cross-validation tests:** each from-scratch algorithm tested against its BioPython equivalent on known inputs.
- **No integration tests required** for portfolio purposes — unit tests are sufficient.
- Tests run via `poetry run pytest`.

---

## Risk and Mitigation

| Risk | Mitigation |
|---|---|
| pysam keeps breaking on Python upgrades | Lazy import isolates the breakage. User can downgrade Python or pin pysam later if SAM/BAM matters. |
| Scope creep into "Tiny does everything" | This spec explicitly forbids it. Cancer analysis lives in a separate repo. Review scope at the end of each phase. |
| Wet-lab time pressure stalling the project | 5-7 hr/week cadence is designed to absorb this. Skipping a week is fine. |
| Algorithm too hard (especially BWT/FM-index) | Phase 5 is intentionally long. Plenty of public references (Langmead's notes, "Bioinformatics Algorithms" textbook). |
| Notebooks become a chore | Keep them small. 500-1000 words + working code, not 5000-word essays. Cross-post to blog so they double as content. |

---

## Success Criteria

End-state criteria for "Tiny is done" (i.e., ready as a polished portfolio piece):

1. `poetry install && poetry run tiny --help` works on a fresh checkout.
2. Five algorithm modules under `tiny/algorithms/`, each with: implementation, CLI command, tests, notebook.
3. Trunk modules unchanged in API since Phase 0.
4. README rewritten to reflect the new identity: focused educational tool, not kitchen-sink CLI.
5. At least 3 blog posts published on akarn.netlify.app derived from the notebooks.
6. The codebase is something the user is proud to show.

---

## Related Artifacts (Tracked Separately)

These are part of the broader portfolio strategy but **not** part of this Tiny spec:

- **Cancer Genomics Analysis Project** — separate repo. Project A from brainstorming: mutation signature analysis on TCGA data, reproducing parts of Alexandrov et al. Kicked off around Week 7 of the Tiny timeline, runs in parallel.
- **Blog series** — content lives at akarn.netlify.app (primary), cross-posted to Substack (secondary). Each Tiny notebook becomes ~1 blog post.

These will get their own spec documents when their respective phases begin.

---

## Open Questions (for User Review)

1. Should the BioPython-backed `tiny align` be **removed** once `tiny align-dp` is solid, or kept as a "production" alongside the "learning" version? Current spec says keep both.
2. Should `RNASequence` be a separate class, or a flag on `DNASequence`? Trade-off: separate class is cleaner; flag is less code. Current spec leaves this open.
3. Notebooks vs Quarto: notebooks are simple but Quarto generates better-looking HTML. Defer the choice until first notebook is written.
4. Do you want a `tiny --version` flag added in Phase 0? Trivially useful; not listed above.

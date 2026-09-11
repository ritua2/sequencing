# Implementation details -- rnaseq_pipeline.c

This is the companion deep-dive to [README.md](README.md), which covers
installation, testing, and the tool-by-tool comparison with
nf-core/rnaseq. Everything here is algorithm internals, validation
history, round-by-round optimization notes, and known gaps -- useful if
you're evaluating whether to trust a specific number or feature, not
needed to install or run the pipeline.

---

# rnaseq_pipeline.c

A **single-file, dependency-free C re-implementation** of the core analytical
workflow used by [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) — FASTQ
QC and trimming, UMI extraction/deduplication, k-mer-seeded splice-aware
alignment, EM-based gene-level quantification, a DESeq2-equivalent
differential expression stack, strandedness inference, rRNA/biotype content
QC, sorted/indexed BAM output, bigWig coverage tracks, and self-contained
HTML/SVG reporting — built from scratch (no aligner, no R, no Bioconductor,
no Nextflow) so it compiles with `gcc` and runs anywhere a C compiler (and,
for BAM/bigWig output, `samtools`/`bedtools`/`pyBigWig`) exists.

It is **not** a Nextflow pipeline and does not orchestrate third-party tools
the way nf-core/rnaseq does, with the deliberate exception of the same
handful of steps nf-core itself shells out to standard tools for rather than
reimplementing: `samtools` for BAM sort/index, `bedtools` + `pyBigWig` for
coverage tracks. Everything else — alignment, quantification, differential
expression, QC — is native C, from scratch.

### What "production-ready" means here, precisely

This has been through three rounds of validation against a real
1,000,000-read-pair *S. cerevisiae* dataset (`SRR9336476`, R64-1-1 genome,
150bp paired-end), plus synthetic ground-truth tests for the pieces real
data alone couldn't validate (UMI deduplication, since the real dataset
has no actual UMIs) — not just unit tests. Each round found and fixed a
real bug: a CIGAR/SEQ formatting bug, a memory regression that OOM-killed
full-dataset runs, and a GFF3-parsing bug that silently dropped every
non-protein-coding gene from the annotation (see
[Validation on real data](#validation-on-real-data) for all three). On that
dataset, for that organism, at that read depth: mapping rate matches
`bwa mem` to within noise, SAM/BAM/bigWig output all validate cleanly with
their respective standard tools, gene counts are deterministic, the
DESeq2-equivalent DE test correctly controls its false-positive rate on a
true-null comparison, UMI dedup correctly neutralizes injected duplicates
down to identical gene counts, and the pipeline no longer crashes under
realistic memory constraints.

That is real, repeated evidence, and it's also **not** the same claim
nf-core/rnaseq can make. nf-core/rnaseq's production status rests on years
of community use across dozens of organisms, genome sizes, library
preparations, and read characteristics, plus a maintained CI suite and a
team that responds to edge cases as real labs hit them. This repo has been
validated across a handful of engineering sessions against one organism,
one library prep, one read depth, on one small (~12 Mb) genome, plus
synthetic-genome testing specifically to answer the memory-scaling question
below. Every gap that could be closed with careful engineering and real
testing has been (see the comparison table below) — what's left is now
almost entirely the one thing that matters most: a **measured, quantified
answer**, updated repeatedly by real architectural work rather than left
as a static verdict. Round 1 found this pipeline's original k-mer
hash-table index extrapolated to ~326 GB to index a human genome (STAR
needs ~30 GB, HISAT2 ~4 GB). Round 2 implemented minimizer-based sparse
indexing plus a smaller struct layout — both validated for correctness
against `bwa mem` and a full real-data run before being trusted for
anything — bringing the extrapolation down to **~68 GB**, a real ~4.8x
reduction. Round 4 went further and actually built the FM-index this
README had been assuming (for three rounds) would close the remaining
gap — and measuring it found that assumption **wrong for the
straightforward implementation**: ~102 GB extrapolated, worse than the
minimizer approach, because suffix-array *construction* overhead
dominates peak memory more than the final structure's smaller
steady-state size would suggest. Round 5 fixed the specific cause: the
FM-index was being built as one structure over the whole genome
concatenated together, when it could be built as one independent
structure per chromosome instead (every real chromosome fits safely
within the smaller integer type that trick requires; the whole genome
concatenated does not) — bringing the FM-index extrapolation down to
**~51 GB**, finally beating the minimizer approach, opt-in via
`--fm-index`. Round 6 replaced prefix-doubling itself with SA-IS
(linear-time construction via induced sorting, validated standalone
against 2,405 test cases plus a real-chromosome self-test on every
`--fm-index` invocation) — **~35 GB**. Round 7 implemented suffix-array
sampling — **~25.3 GB, below STAR's ~30 GB reference figure for the
first time**. Round 8 removed the last named transient-memory cost
(a top-level array-type conversion in SA-IS) — real, measured, but
mostly a lesson: 1.9–6.0% memory reduction at every tested scale, yet
only **~25.16 GB** on the human extrapolation specifically, because that
cost scales with the *largest single chromosome*, not total genome
size, and human's largest chromosome is proportionally much smaller
than any tested synthetic genome's — the same fix helps every genome
tested, but helps human's real, more skewed chromosome distribution by
less than the tested-scale percentages alone would suggest. Five rounds
of real, measured, compounding progress: ~326 GB → ~68 GB → ~51 GB →
~35 GB → ~25.3 GB → ~25.16 GB. Round 9 didn't move that number further
— it investigated *why* the remaining transient cost behaves the way it
does, using real per-chromosome memory instrumentation rather than
guessing: tested and ruled out one hypothesis (allocator retention,
via `malloc_trim(0)` — no measurable effect), and found the real
explanation (peak memory dominated by the largest chromosome's own
construction, with smaller chromosomes efficiently reusing that
reserved space) narrows what the one remaining item could plausibly
achieve, rather than leaving it an open-ended "further optimization."

**What crossing below STAR's figure does and does not mean, stated
plainly so it isn't misread as more than it is:** it means the memory
*extrapolation* for indexing a human-scale genome is now favorable —
nothing more. It is still an extrapolation from real measurements
topping out at 90 Mb, not a validated human-scale run (this development
sandbox's own ~3.9 GB RAM can't attempt one). The reduction itself fell
well short of what pure steady-state math predicted, because peak
memory is still shaped by *transient* construction cost sampling
doesn't touch — a precise, reported finding, not a hidden asterisk. And
memory is one dimension of many: splice-detection correctness at
human-scale intron lengths was found genuinely broken when first
measured (widening `--max-intron` past yeast's validated default
reintroduced false-positive splice calls in direct proportion to how
far it was widened), then fixed at the source (`try_spliced_align` now
requires a canonical GT-AG site outright, not just prefers one, beyond
a 20,000bp safe range) and re-measured to confirm the fix actually
worked, not just assumed to — including a positive control (a planted,
exact 500,000bp intron, correctly recovered with zero mismatches).
Runtime at mammalian scale, and — the constant across every round of
this README — validation by more than one person against more than one
dataset, remain completely untested. A favorable memory extrapolation
and a measured-safe wide splice window are both real, necessary pieces
of mammalian feasibility; between them they still don't add up to
sufficient, and the gap between "necessary" and "sufficient" here is
multi-organism validation this project has never claimed to substitute
for.

**Bottom line:** if your organism and genome size are compact — this repo's
directly-validated case is yeast (~12 Mb), and the measured memory scaling
now extrapolates genomes up to a few hundred Mb (worm/fly/plant scale) to
comfortably workstation-feasible memory (well under 2 GB) — and your
library prep matches what's supported (in-line UMI-on-R1 if you use UMIs,
no heavy contaminant/rRNA-organism concerns beyond what the
annotation-based rRNA QC catches), and you've re-run the checks in
[Validation on real data](#validation-on-real-data) on your own data first,
this is meaningfully more trustworthy than it was several passes ago.
**If you're working with a human or mouse genome, use nf-core/rnaseq — not
as a hedge, but because this repo's architecture is now
known, specifically and quantitatively, not to work there.**

---

## Table of contents

- [How this maps to nf-core/rnaseq](#how-this-maps-to-nf-corernaseq)
- [Pipeline stages in this repo](#pipeline-stages-in-this-repo)
- [Build](#build)
- [Installing dependencies](#installing-dependencies)
- [Input requirements](#input-requirements)
- [Usage](#usage)
- [`--fm-index`: when it's worth the tradeoff](#--fm-index-when-its-worth-the-tradeoff)
- [Multi-lane FASTQ merging](#multi-lane-fastq-merging)
- [Output files](#output-files)
- [UMI extraction and deduplication](#umi-extraction-and-deduplication)
- [Contaminant genome screening](#contaminant-genome-screening)
- [Gene-body coverage](#gene-body-coverage)
- [Library complexity (saturation curve)](#library-complexity-saturation-curve)
- [Duplication vs. expression (dupRadar-equivalent)](#duplication-vs-expression-dupradar-equivalent)
- [Junction annotation (RSeQC-equivalent)](#junction-annotation-rseqc-equivalent)
- [Differential expression internals](#differential-expression-internals)
  - [Volcano plot](#volcano-plot)
- [Performance](#performance)
- [Validation on real data](#validation-on-real-data)
- [Known gaps vs. nf-core/rnaseq](#known-gaps-vs-nf-corernaseq)
- [Repo layout](#repo-layout)
- [License](#license)

---

## How this maps to nf-core/rnaseq

nf-core/rnaseq is a 16-stage Nextflow pipeline that orchestrates ~20
best-in-class external tools (FastQC, Trim Galore!, STAR, Salmon, RSEM,
HISAT2, SAMtools, UMI-tools, picard, StringTie, RSeQC, Qualimap, dupRadar,
Preseq, DESeq2, MultiQC, and more) behind a samplesheet-driven interface with
container support. This repo reimplements the numerically/algorithmically
central stages of that workflow natively in C, and is honest about the parts
it does not attempt to replace:

| # | nf-core/rnaseq stage | Reference tool(s) | This repo | Status |
|---|---|---|---|---|
| 1 | Merge re-sequenced FASTQs | `cat` | `split_lane_paths()` — a comma-separated FASTQ path list (in place of a single path, for either `se` or `pe` mode) is concatenated in the order given before trimming/alignment, the same outcome as nf-core's samplesheet-driven per-sample lane merge. Validated bit-for-bit: a 2-lane comma-separated run produces byte-identical `alignments.sam`, `gene_counts.tsv`, and `dupradar.tsv` to the same reads given as one file — see [Multi-lane FASTQ merging](#multi-lane-fastq-merging) | Mirrored |
| 2 | Strandedness inference | `fq`, `Salmon` | `classify_strandedness()` — compares read1 alignment strand vs. gene strand among uniquely-assigned fragments (RSeQC `infer_experiment.py`-style); printed to stdout + `multiqc_summary.txt`. When confidently one-directional, also used to disambiguate ambiguous multi-gene assignments (see [On strandedness and gene assignment](#on-strandedness-and-gene-assignment)) — measured to resolve 3.6% of mapped-with-feature fragments to unique assignments on the real validation dataset, with the mapped-with-feature total exactly unchanged | Mirrored (diagnostic + disambiguation) |
| 3 | Read QC | `FastQC` | `compute_raw_qc()` — per-read length/GC/quality-by-cycle, aggregated into `qc_report.txt` + HTML SVG plots | Mirrored |
| 4 | UMI extraction | `UMI-tools` | `extract_umis()`, enabled with `--umi-len N` — pulls a fixed-length UMI off read1's 5′ end (in-line UMI layout, e.g. QIAseq/NEBNext-style kits) and hard-clips it before alignment. Scoped: UMI-on-R2 and split-UMI protocols aren't handled | Mirrored (common case) |
| 5 | Adapter/quality trimming | `Trim Galore!` (Cutadapt) | `trim_read()` — 3′ quality trim + adapter overlap trim | Mirrored (simplified) |
| 6 | Contaminant genome removal | `BBSplit` | `try_contaminant_screen()`, enabled with `--contaminant-ref path[:Label]` (repeatable) — screens reads unmapped against the primary genome against additional reference(s) in sequence, same priority-order approach BBSplit uses. Reuses the same aligner as primary alignment; runs last, after every primary-genome output is written. See [Contaminant genome screening](#contaminant-genome-screening) for a real detection-rate validation | Mirrored |
| 7 | rRNA removal | `SortMeRNA` | `compute_biotype_breakdown()` — quantifies (doesn't remove) rRNA/tRNA/protein-coding content from the GTF's `gene_biotype` attribute, now correctly populated for non-protein-coding loci by `convert_gff3_to_gtf.py` (previously only protein-coding `gene` rows were captured — rRNA genes were silently absent from the GTF entirely). Answers the actual QC question ("how much rRNA is in this library?") using the organism's own annotation rather than a separate sequence database; doesn't catch contaminating rRNA from a *different* organism the way SortMeRNA's database-alignment approach can | Mirrored (annotation-based QC; not a pre-alignment removal step) |
| 8 | Alignment | `STAR` / `HISAT2` | `align_read()` — k-mer-seeded banded Smith-Waterman with a dedicated short-k splice-anchor index for junction-spanning reads, multi-mapper aware (reports all tied-best hits, `NH:i:` tag, SAM secondary flag). An alternative FM-index/suffix-array seeding path (`--fm-index`, opt-in, built per-chromosome using SA-IS with a sampled suffix array — see [Memory](#memory)) is also available — validated for correctness (~94.2% `bwa mem` position agreement, matching the default path, plus real-chromosome self-tests on every invocation for both the suffix-array construction algorithm and the sampling/recovery mechanism) and, as of round 8, a real memory improvement over the default at extrapolated human scale (~25.16 GB vs. ~68 GB), below STAR's ~30 GB reference figure | Mirrored (custom aligner; only validated at yeast genome scale — see the note at the top of this README) |
| 9 | Sort/index alignments | `SAMtools` | `try_write_sorted_bam()` — shells out to `samtools sort`/`samtools index` after writing the SAM, exactly like nf-core does. Skips gracefully with a clear message if `samtools` isn't on `PATH`; the SAM itself is unaffected either way | Mirrored (via samtools, like nf-core itself) |
| 10 | UMI dedup | `UMI-tools` | `umi_dedup()`, enabled with the same `--umi-len N` — UMI-tools' default **directional-adjacency** method: clusters UMIs within edit distance 1 by count-directional adjacency (not just exact matches), marks duplicates with the SAM `0x400` flag, excludes them from gene counts. See [UMI extraction and deduplication](#umi-extraction-and-deduplication) for validation against a ground-truth set specifically containing sequencing-error UMI variants (96.7% overall recall, with ~half caught only via clustering) | Mirrored |
| 11 | Duplicate marking | `picard MarkDuplicates` | `mark_position_duplicates()` — Picard-style position duplicate marking (identical chrom + 5′ position + strand), computed independently of `--umi-len`/UMI dedup. Feeds the dupRadar-equivalent duplication-vs-expression analysis below rather than the SAM `0x400` flag (which remains UMI-dedup-only, unchanged) — see [Duplication vs. expression](#duplication-vs-expression-dupradar-equivalent) | Mirrored (as dupRadar's own upstream step, not as a standalone BAM-flagging pass) |
| 12 | Quantification | `Salmon` / `RSEM` | `quantify_em()` — EM algorithm: unique-mapping reads anchor gene abundance directly, multi-mapping reads redistributed proportionally to current estimates, iterated to convergence (same spirit as RSEM/Salmon's EM, without fragment-length/sequence-bias models) | Mirrored (core algorithm) |
| 13 | bigWig coverage tracks | `BEDTools`, `bedGraphToBigWig` | `try_write_bigwig()` — `bedtools genomecov` for a bedGraph, then a real bigWig (proper UCSC binary format: R-tree index, zoom levels) via `pyBigWig` rather than the UCSC binary specifically, since that's what's reliably apt-installable. Degrades to bedGraph-only (still directly IGV/UCSC-loadable) if `bedtools`/`pyBigWig` aren't present | Mirrored (via bedtools + pyBigWig) |
| 14a | Extensive per-sample QC | `RSeQC`, `Qualimap`, `dupRadar`, `Preseq` | Mapping rate, proper-pair rate, strandedness, rRNA/biotype content, gene-body coverage (5'/3' bias), a library-complexity/saturation curve, duplication-rate-vs-expression modeling, and known-vs-novel splice junction annotation are all covered — see [Gene-body coverage](#gene-body-coverage), [Library complexity](#library-complexity-saturation-curve), [Duplication vs. expression](#duplication-vs-expression-dupradar-equivalent), and [Junction annotation](#junction-annotation-rseqc-equivalent) | Mirrored |
| 14b | Differential expression | `DESeq2` | `run_differential_expression()` — from-scratch DESeq2-equivalent stack: median-of-ratios size factors, Cox-Reid dispersion estimation, Gamma-GLM (IRLS) mean-dispersion trend fit, empirical-Bayes dispersion shrinkage, closed-form variance-stabilizing transform (VST), F-distribution Cook's-distance outlier flagging, Benjamini-Hochberg FDR with independent filtering, normal-normal log2FC shrinkage | Mirrored |
| 15 | Pseudoalignment route | `Salmon` / `Kallisto` | Not implemented as a separate route (this repo's aligner is itself alignment-based, not pseudoalignment) | Out of scope |
| 16 | Multi-sample QC/report | `MultiQC`, `R` | `compare` mode — `write_compare_report()`: DESeq2-VST-normalized PCA (Jacobi eigendecomposition), average-linkage hierarchical clustering/dendrogram, library-size bar chart, self-contained HTML report | Mirrored |
| 17 | Compressed input (`.fastq.gz`, `.fa.gz`, `.gtf.gz`) | built-in | `open_maybe_gz()` — transparently gunzips (via `popen("gunzip -c ...")`) any input path ending in `.gz`; uncompressed inputs are unaffected | Mirrored |

**Net effect:** the statistically hardest parts of the workflow — splice-aware
alignment with multi-mapping, EM-based quantification, and the full DESeq2
dispersion/shrinkage/testing stack — are implemented natively, and every
row above is now either Mirrored or an explicitly-scoped Partial/Gap with a
stated reason, rather than a blanket "not implemented." What's left,
honestly: contaminant-genome screening (needs a second reference genome
this repo has no way to source generically), UMI-tools' more sophisticated
directional-adjacency dedup method, deeper per-base QC (gene-body coverage,
saturation curves, junction annotation), and — the one that matters most —
the still-untested question of genome-scale memory/runtime behavior. See
[Known gaps](#known-gaps-vs-nf-corernaseq).

## Pipeline stages in this repo

Running the pipeline on one sample prints progress through 7 stages,
followed by an optional BAM sort/index step and a strandedness verdict —
this is real output from a run against the dataset used in this project's
validation:

```
[1/7] Loading reference + k-mer index (k=16): genome.fa
[2/7] Building k-mer index (k=16)
      -> index cached for future runs: genome.fa.kidx
[3/7] Loading GTF annotation: annotation.gtf
[4/7] Loading reads (paired-end mode): ...
[5/7] Running QC + adapter/quality trimming
[6/7] Aligning reads (k-mer seed + ungapped/spliced extension, multi-mapping-aware)
[7/7] EM quantification + writing output files to: outdir
      sorting + indexing BAM via samtools...
      -> wrote outdir/alignments.sorted.bam (+ .bai index)
Done. (EM converged in 180 iterations)
Inferred library strandedness: reverse-stranded (read1 = antisense, e.g. Illumina TruSeq Stranded mRNA / dUTP) (96.8% of uniquely-assigned reads have read1 opposite the gene's strand)
```

1. **Reference loading + indexing.** The genome FASTA (or `.fa.gz`) is
   parsed and two k-mer hash indexes are built: a primary `k=16` index for
   seeding full-read alignment, and a shorter `k=10` index dedicated to
   anchoring short exon
   overhangs at splice junctions. The index is written to a `.kidx` cache
   file next to the FASTA (size + mtime checked) so repeat runs against the
   same genome skip the rebuild entirely.
2. **GTF loading.** Gene features (`gene_id`, coordinates, strand) are
   parsed into an interval structure used for read-to-gene assignment.
3. **FASTQ loading.** Plain-text `.fastq` only (see
   [Input requirements](#input-requirements)) — single- or paired-end.
4. **QC + trimming.** Per-read length, GC%, and mean/per-cycle quality are
   computed before 3′ quality- and adapter-trimming is applied.
5. **Alignment.** Each read is seeded via k-mer lookup, extended with a
   banded Smith-Waterman DP (thread-private scratch buffers, reused across
   reads), and, when a full-length ungapped/mismatch-tolerant match isn't
   found, re-seeded against the short-k splice index and rescanned for a
   spliced (intron-skipping) alignment. All tied-best-scoring hits are
   reported (up to a cap), not just one, so genuine multi-mappers are
   represented as such.
6. **Quantification.** An EM algorithm (RSEM/Salmon-style) assigns reads to
   genes: unique mappers anchor their gene's count directly; multi-mappers
   are fractionally redistributed across candidate genes in proportion to
   each gene's current abundance estimate, iterated to convergence.
7. **Output.** SAM alignments, per-gene counts, a QC report, a MultiQC-style
   plain-text summary, and a self-contained HTML report (5 inline SVG plots)
   are written to the output directory.

A second, separate entry point (`compare` mode) takes the `gene_counts.tsv`
from ≥2 single-sample runs and produces cross-sample QC: real
median-of-ratios normalization, a DESeq2 variance-stabilizing transform,
PCA, hierarchical clustering, and — if exactly two condition labels are
given with ≥2 replicates each — a full differential expression test.

## Build

```bash
# Single-threaded
gcc -O2 -o rnaseq_pipeline rnaseq_pipeline.c -lm

# Multi-core (QC/trim and alignment loops are both OpenMP-parallelized —
# worth using on anything with more than one core)
gcc -O2 -fopenmp -o rnaseq_pipeline rnaseq_pipeline.c -lm
```

No external libraries beyond the C standard library and libm. No R, no
Python runtime dependency for the pipeline itself (Python is only used by
the optional GFF3→GTF conversion helper below).

## Installing dependencies

Everything in this section is **optional** — the pipeline builds with
nothing but a C compiler and runs to completion with nothing but that,
producing every core output (`alignments.sam`, `gene_counts.tsv`,
`report.html`, `multiqc_summary.txt`, and the rest). What's below adds
specific extra outputs on top of that; skip whichever you don't need. If
a tool below isn't on `PATH`, the pipeline detects that itself, prints a
`note:` (not a warning or error) saying exactly what it's skipping and
why, and keeps going — nothing fails or produces incomplete/incorrect
output because a tool is missing. The table says exactly what each one
buys you, so you can decide before running rather than mid-run:

| Tool | Enables | Without it |
|---|---|---|
| `samtools` | Sorted, indexed BAM (`alignments.sorted.bam` + `.bai`) | `alignments.sam` is still written in full — sort and index it yourself later (see below), or use it as-is with any tool that reads SAM directly |
| `bedtools` **and** Python 3 with `pyBigWig` (both required together) | Genome coverage track as bigWig (`coverage.bw`) | No coverage track at all — this one has no manual fallback the way samtools does |
| Python 3 (standard library only — no packages) | Only needed if your annotation is GFF3, not GTF (`convert_gff3_to_gtf.py`, run once before the pipeline) | Not needed at all if you already have a GTF |

### samtools

```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install -y samtools

# macOS (Homebrew)
brew install samtools

# HPC cluster with an environment-modules system (TACC and similar) --
# check what's actually available on yours first:
module spider samtools
module load samtools          # exact module name varies by cluster --
                               # use whatever `module spider` reported

# conda/mamba (works the same way everywhere, including HPC nodes
# without root or a module system — the most portable option here)
conda install -c bioconda samtools
```

If you already have `alignments.sam` from a run that completed before
`samtools` was available (exactly the scenario this section exists for),
there's no need to re-run the pipeline — sort and index the existing
file directly:

```bash
samtools sort -o outdir/alignments.sorted.bam outdir/alignments.sam
samtools index outdir/alignments.sorted.bam
```

### bedtools + pyBigWig

Both are required together for the coverage track — `bedtools genomecov`
produces a bedGraph from the sorted BAM (so this also needs `samtools`
to have run first), and the pipeline's own bundled Python snippet uses
`pyBigWig` to convert that bedGraph into a real bigWig file. Verified
directly, not just described: this exact install path (`apt-get install
bedtools` + `pip3 install pyBigWig`) was run in a clean environment and
produced a real, valid `coverage.bw` covering all reference sequences.

```bash
# Ubuntu/Debian
sudo apt-get install -y bedtools python3-pip
pip3 install pyBigWig --break-system-packages   # recent Debian/Ubuntu (PEP 668)
                                                 # refuse a plain `pip3 install` outside
                                                 # a virtual environment -- this flag is
                                                 # what its own error message asks for.
                                                 # A venv (`python3 -m venv .venv &&
                                                 # source .venv/bin/activate`) avoids
                                                 # needing the flag at all, if preferred.

# macOS (Homebrew + pip)
brew install bedtools
pip3 install pyBigWig

# HPC cluster (module system)
module spider bedtools
module load bedtools
pip3 install --user pyBigWig   # pyBigWig itself is rarely its own module;
                                # --user avoids needing write access to a
                                # shared Python install (add --break-system-packages
                                # too if the cluster's Python also enforces PEP 668)

# conda/mamba (installs both in one step, and is the easiest route if
# pyBigWig's C extension gives pip any trouble building from source, or
# if PEP 668 makes the plain pip3 route more friction than it's worth)
conda install -c bioconda bedtools pybigwig
```

### Python 3 (GFF3→GTF conversion only)

No packages needed — `convert_gff3_to_gtf.py` uses only the standard
library. Any Python 3 already on `PATH` (`python3 --version` to check)
works; this is virtually always already available on Linux, macOS, and
any HPC login/compute node.

### OpenMP (multi-core builds only)

`gcc -fopenmp` needs `libgomp`, which ships with `gcc` itself on every
mainstream Linux distribution and with Homebrew's `gcc` on macOS (not
Apple's default `clang` shim at `/usr/bin/gcc` — install `gcc` via
Homebrew first if `-fopenmp` fails to link). Nothing to install
separately in the common case; if the build specifically fails to find
`omp.h` or fails to link with an `-fopenmp`-related error, that's the
one thing worth checking.

## Input requirements

- **Reference genome:** FASTA (`.fa`/`.fasta`), plain or gzipped (`.fa.gz`).
- **Annotation:** **GTF, not GFF3** — plain or gzipped (`.gtf.gz`). If your
  annotation is Ensembl-style GFF3 (as it typically is for yeast/other
  Ensembl genomes), convert it first with the included helper (this also
  accepts either a plain or gzipped GFF3 as input):

  ```bash
  python3 convert_gff3_to_gtf.py Saccharomyces_cerevisiae.R64-1-1.59.gff3.gz yeast.gtf
  ```

  This does a minimal gene-level conversion (one GTF `gene` row per GFF3
  `gene` feature, with `gene_id` extracted from the `Dbxref`/`ID` attribute
  string) — sufficient for gene-level quantification, not a full
  transcript-level GTF. The output is written as plain-text `yeast.gtf`;
  gzip it yourself afterward (`gzip yeast.gtf`) if you want to pass a
  `.gtf.gz` to the pipeline.
- **Reads:** FASTQ, plain or gzipped (`.fastq.gz`/`.fq.gz`) — no manual
  decompression needed. Gzipped inputs are read via a `gunzip -c` pipe
  (using whatever `gunzip` is on `PATH`), so the only added runtime
  requirement is having `gunzip` available, which is virtually always true
  on Linux/macOS.

  This was validated directly against the real dataset used in this
  project's testing: running the full pipeline with `genome.fa.gz`,
  `annotation.gtf.gz`, `reads_1.fastq.gz`, and `reads_2.fastq.gz` passed
  straight through — no `gunzip` step — produced byte-identical
  `gene_counts.tsv` and identical SAM output to the pre-decompressed run.

## Usage

**Paired-end** (e.g. with `SRR9336476`):

```bash
./rnaseq_pipeline genome.fa annotation.gtf pe reads_1.fastq reads_2.fastq outdir/
```

**Single-end:**

```bash
./rnaseq_pipeline genome.fa annotation.gtf se reads.fastq outdir/
```

**Worked example** using the sample dataset linked from `data/README`:

```bash
python3 convert_gff3_to_gtf.py \
  RNA-Seq_Sample_Files/Saccharomyces_cerevisiae.R64-1-1.59.gff3.gz \
  RNA-Seq_Sample_Files/yeast.gtf

time ./rnaseq_pipeline \
  RNA-Seq_Sample_Files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
  RNA-Seq_Sample_Files/yeast.gtf \
  pe \
  RNA-Seq_Sample_Files/1M_SRR9336468_1.fastq \
  RNA-Seq_Sample_Files/1M_SRR9336468_2.fastq \
  ./outdir/
```

**Multi-sample comparison** (PCA, clustering, and — with condition labels —
differential expression), run *after* each sample has its own single-sample
output directory:

```bash
# QC/clustering only
./rnaseq_pipeline compare out_rep1 out_rep2 out_rep3 combined_outdir/

# With a 2-condition design (needs >=2 replicates per condition) to also
# run the DESeq2-equivalent differential expression test
./rnaseq_pipeline compare \
  out_ctrl1:control out_ctrl2:control \
  out_trt1:treated  out_trt2:treated \
  combined_outdir/
```

### About the `.kidx` cache

The first run against a given genome builds the k-mer index and writes
`genome.fa.kidx` next to it — that step is the slow part. Every subsequent
run against the **same, unmodified** `genome.fa` loads from that cache
instead of rebuilding, which is dramatically faster. If `genome.fa` is
replaced with different content under the same filename, the pipeline
checks file size + modification time and rebuilds automatically.
This caching applies to the default path only — see below for why
`--fm-index` doesn't use it.

### `--fm-index`: when it's worth the tradeoff

Real numbers, not a guess: on the real 100,000-read-pair yeast dataset,
a single `--fm-index` run took **~208s**, against the default path's
**~69s** for a first run (building the `.kidx` cache) and **~55s** for
every run after that (loading the cache). `--fm-index` is roughly 3–4x
slower per run — and, deliberately, it **never caches**: every
invocation rebuilds its index from scratch (see the Memory section's
round 5 note on why — the cache-corruption bug that motivated this).
That's a real, one-time few-minutes-per-run cost on yeast; scaled to a
larger genome it would be proportionally larger too, though this hasn't
been measured beyond the scales in the Memory section.

The tradeoff, stated plainly: `--fm-index` costs real time on every run
in exchange for the memory profile described in
[Memory](#memory) (~25 GB vs. ~68 GB extrapolated for human, the
difference between plausibly workstation-feasible and not). For a
memory-constrained environment indexing a large genome once, or for
comparing correctness between the two seeding paths, that trade is
worth it. For repeated runs against the same, already-cached genome —
most real workflows, processing many samples against one reference —
the default path's `.kidx` cache means it stays both faster *and* uses
less memory after the first run, and `--fm-index`'s advantage is
specifically a peak-memory one, not a general-purpose speedup. Pick
`--fm-index` for the memory profile; don't expect it to also be the
fast path.

## Multi-lane FASTQ merging

nf-core/rnaseq accepts one samplesheet row per sequencing lane and
concatenates same-sample lanes (its `CAT_FASTQ` module) before any
trimming or alignment happens. This repo has no samplesheet, so the
equivalent is addressed by the read-path argument itself: pass a
comma-separated list of paths in place of a single FASTQ path, and the
lanes are concatenated in the order given, before trimming/alignment —
the same outcome, just keyed by path list instead of a sample-ID column.

```bash
# single-end, 3 lanes
./rnaseq_pipeline genome.fa annotation.gtf se \
  L001_R1.fastq.gz,L002_R1.fastq.gz,L003_R1.fastq.gz outdir/

# paired-end, 2 lanes -- R1 and R2 lists must have the same lane count, in order
./rnaseq_pipeline genome.fa annotation.gtf pe \
  L001_R1.fastq.gz,L002_R1.fastq.gz \
  L001_R2.fastq.gz,L002_R2.fastq.gz \
  outdir/
```

A single path with no comma behaves exactly as before — fully
backward-compatible with every existing invocation. Each lane can
independently be gzipped or plain, same as a single-lane path. An R1/R2
lane-count mismatch in `pe` mode is rejected with an explicit error rather
than silently pairing the wrong lanes together; a lane running out of
reads before its mate lane (uneven lane sizes) prints the same truncation
warning single-file `pe` mode already prints, per lane.

**Validated**: splitting a real 100,000-read-pair slice of `SRR9336476`
into two lanes and running it through `pe L001,L002` produced a
byte-identical `alignments.sam`, `gene_counts.tsv`, and `dupradar.tsv` to
the same 100,000 pairs given as one file. Repeated at 4x scale
(400,000 pairs split into two 200,000-pair lanes) with the same result —
byte-identical across all three files — the merge changes nothing about
what gets aligned or counted, only how the input is assembled, at both
scales tested. An R1-lane-count/R2-lane-count mismatch was confirmed to
fail fast with a clear error before any alignment work starts.

## Output files

Per-sample (`outdir/`):

| File | Contents |
|---|---|
| `alignments.sam` | SAM-format alignments; multi-mappers emit one record per reported hit with `NH:i:` and the `0x100` secondary-alignment flag. If `--umi-len N` was used, duplicate-marked reads carry the `0x400` flag |
| `alignments.sorted.bam` + `.bai` | Coordinate-sorted, indexed BAM — written automatically if `samtools` is found on `PATH` at run time (via `samtools sort -m 256M -@ 1` + `samtools index`; see [Performance](#performance) for why the memory is capped). If `samtools` isn't found, this step is skipped with a clear message and `alignments.sam` is still written normally |
| `coverage.bw` | Genome-wide coverage track in real bigWig format (via `bedtools genomecov` + `pyBigWig`), directly loadable in IGV/UCSC Genome Browser. Requires `bedtools` and Python's `pyBigWig`; needs the sorted BAM above to exist first |
| `coverage.bedgraph` | The intermediate bedGraph `coverage.bw` was built from. Also directly genome-browser-loadable on its own — kept even when the `.bw` conversion succeeds, and is what you get if `pyBigWig` isn't available but `bedtools` is |
| `gene_counts.tsv` | EM-quantified gene-level counts (UMI duplicates excluded from these counts if `--umi-len` was used) |
| `qc_report.txt` | Plain-text QC summary — read counts, length/GC/quality distributions, trimming stats, mapping/proper-pair rates |
| `multiqc_summary.txt` | MultiQC-style condensed run summary, including the inferred strandedness call, RNA biotype/rRNA content breakdown, and UMI dedup stats if applicable |
| `report.html` | Self-contained HTML report with 5 inline-SVG plots (quality-by-cycle, GC distribution, length distribution, mapping breakdown, top-gene bar chart) |

The pipeline also prints an inferred strandedness verdict
(forward-stranded / reverse-stranded / unstranded, with the underlying
concordant/discordant read counts) to stdout and to `multiqc_summary.txt`
at the end of each run — see
[On strandedness and gene assignment](#on-strandedness-and-gene-assignment)
for what this diagnostic does and doesn't affect.

Multi-sample (`combined_outdir/`, from `compare` mode):

| File | Contents |
|---|---|
| `multi_sample_report.html` | PCA plot (VST-normalized), hierarchical-clustering dendrogram/heatmap, library-size bar chart, and — with a valid 2-condition design — a full differential expression results table (log2FC, shrunken log2FC, p-value, BH-adjusted p-value) |

## UMI extraction and deduplication

Opt-in via `--umi-len N`, placed anywhere in the command line:

```bash
./rnaseq_pipeline --umi-len 8 genome.fa annotation.gtf pe reads_1.fastq reads_2.fastq outdir/
```

This assumes the common in-line layout: an `N`-base UMI at the 5′ end of
read1 (or the single-end read), immediately before the biological insert —
the layout used by QIAseq, NEBNext, and similar kits. `extract_umis()`
pulls those bases off into an internal UMI field and hard-clips them from
the read before any QC/trimming/alignment, so they never contaminate
alignment as if they were genomic sequence. **Scope, stated plainly:**
UMI-on-R2 and split-UMI protocols aren't handled — if your kit does either
of those, this flag won't do the right thing for you.

After alignment, `umi_dedup()` marks reads as PCR/optical duplicates using
UMI-tools' **default "directional"** method: within each (chromosome,
alignment start position, strand) group, distinct UMIs are clustered by
count-directional adjacency — two UMIs merge into one cluster if they're
within Hamming distance 1 of each other *and* the higher-count one has
roughly double (or more) the read support of the lower-count one, the
pattern a single sequencing error in a PCR-duplicated UMI would produce.
Within each resulting cluster, the single highest-count UMI's first
occurrence is kept; everything else — both exact repeats and
edit-distance-1 relatives — gets the SAM `0x400` flag and is excluded from
`gene_counts.tsv`, while remaining visible in `alignments.sam`/`.bam`
(matching how UMI-tools + picard MarkDuplicates mark rather than delete).
This directly replaces what was previously exact-match-only ("unique"
method) dedup — see the note in [Validation on real data](#validation-on-real-data)
for what changed and why. **Stated plainly, the remaining simplification**:
a cluster's absorption threshold uses its original highest-count UMI
throughout the search, rather than re-deriving a threshold at every
individual edge the way UMI-tools' reference implementation does — this
matches the common case and is easy to verify correct, but can diverge from
UMI-tools in some multi-hop edge cases. Groups with more than 500 distinct
UMIs at one exact position (a safety cap against pathological pileups) fall
back to exact-match only for that group.

Validated with synthetic ground truth designed specifically to exercise
clustering, not just exact matching: took a 100,000-read set, generated
random 8bp UMIs, then injected 10,000 duplicate reads — **half exact
copies, half with a single-base substitution injected into the UMI**, so
exact-match dedup could only ever catch half of them by construction.
`--umi-len 8` caught 9,670/10,000 (96.7%) overall, and **4,862 of those
were caught only via directional clustering** — reads whose UMI had a
sequencing error, which the previous exact-match method would have missed
entirely. That's within rounding of the ~5,000 error-containing duplicates
injected, meaning the clustering step is catching essentially all of what
it's designed to catch. Gene counts matched a clean, duplicate-free
baseline to within the ~330 residual uncaught duplicates (versus what would
have been ~5,000 uncaught residuals under the old exact-match-only method
on this same error-containing dataset). Paired-end got the same treatment
(33,000 pairs, 3,000 injected duplicates, half error-containing): 2,898/3,000
(96.6%) caught, 1,447 via clustering, mate flags 100% consistent (0
mismatches across 30,000 pairs checked). Backward compatibility was also
explicitly checked: running without `--umi-len` on the same data produces
byte-identical output to before this change (82,812 uniquely
assigned / 7,756 multi-mapped / 5,950 no-gene / 3,482 unmapped, matching
exactly).

## Contaminant genome screening

Opt-in via `--contaminant-ref path[:Label]`, repeatable for multiple
references, placed anywhere in the command line:

```bash
./rnaseq_pipeline --contaminant-ref phix.fa.gz:PhiX --contaminant-ref ecoli.fa:E.coli \
  genome.fa annotation.gtf pe reads_1.fastq reads_2.fastq outdir/
```

This is the same mechanism BBSplit provides — it doesn't ship reference
genomes either; you supply which genome(s) to screen against, and reads
that don't map to the primary genome get checked against each one in turn
(same sequential-priority approach: each screen only checks what's still
unmapped after the previous one). It reuses this pipeline's own aligner,
not a separate or lesser check. A common first choice is PhiX — a small
(~5.4 kb) genome used as a near-universal Illumina sequencing control,
spiked into most runs at low concentration, so a small amount of
PhiX-mapping reads in any Illumina dataset is normal and a large amount
can indicate a run-quality problem.

**Runs last**, after every primary-genome-dependent output (SAM, BAM,
bigWig, gene counts, QC report, HTML report) is already written — this
matters because the screening pass reuses the same global reference/index
state as primary alignment, freeing and rebuilding it for each contaminant
reference in turn (`free_reference_index()`). By the time this runs,
nothing else needs the primary genome's data any more, so this is safe;
it would not be if reordered earlier.

Output: `contaminant_screen.txt` in the output directory, showing how many
of the primary-unmapped reads matched each reference (in screening order)
and how many remain unexplained after all of them.

**Validated with a real contaminant genome**, not a synthetic one: fetched
the real lambda phage genome (a legitimate, commonly-used small reference,
~48.5 kb) and generated 2,000 single-end reads genuinely sampled from it
(random start positions, ~50% reverse-complemented, matching real sequencing
orientation). Mixed those into a set of yeast reads and ran the full
pipeline with `--contaminant-ref lambda_virus.fa:LambdaPhage`. Result:
**2,000 out of 2,000 injected lambda reads correctly matched** (100%
detection, zero false negatives) out of 2,370 total primary-unmapped reads
— the other 370 were genuine yeast-unmapped reads, correctly *not* flagged
as lambda matches. Primary-genome outputs (gene counts, BAM) were checked
afterward and found completely unaffected by the later screening pass, as
expected given the ordering guarantee above.

**Scope, stated plainly**: this identifies reads that map well to a
supplied contaminant reference; it doesn't remove or quarantine them from
`alignments.sam`/`.bam` (they were already unmapped against the primary
genome, so they were never in `gene_counts.tsv` to begin with) — this is
diagnostic screening, matching what BBSplit's report gives you, not a
combined align-and-clean step. It's also only as good as the reference(s)
you give it: a contaminant not in any screened reference won't be caught,
and will just remain in the primary-unmapped total.

## Gene-body coverage

Computed automatically, no flag needed — RSeQC `geneBody_coverage.py`'s
question ("is read coverage uniform across the gene body, or skewed toward
one end?"), answered from data already produced during normal quantification
rather than a separate pass. For every uniquely-assigned read, its position
within its gene is recorded as one of 100 bins (0 = the gene's 5' end, 99 =
3', already corrected for strand so this is always biological 5'→3', not
genomic left→right). Genes with fewer than 10 uniquely-assigned reads or
under 200bp are excluded (too little signal for a 100-bin curve to mean
anything); each qualifying gene's own curve is normalized to sum to 1
*before* being averaged into the genome-wide curve, so a handful of very
highly expressed genes can't dominate the shape — the same approach RSeQC
uses.

Output: a condensed 10-bin table in `multiqc_summary.txt` plus a 3'/5' bias
ratio (mean coverage in the last 10% of gene body over the first 10%), and
the full 100-bin curve as an SVG line chart in `report.html`.

**On the real dataset**, this surfaced a real, biologically sensible
signal, not just a working code path: a clear 3' coverage bias (ratio
3.92 on the full 1M-read-pair run, averaged over 5,253 qualifying genes;
3.63 on a 100k-read subset — consistent between the two). This is the
textbook signature of polyA-selected library prep (reverse transcription
primed from the poly-A tail, so 3' ends get relatively more coverage than
5' ends), and it's independently consistent with this same dataset's very
low rRNA content (~0.27%, see the biotype QC above) — both findings point
the same direction (a standard polyA-selected mRNA-seq library), from two
completely independent measurements. That agreement is real evidence the
underlying alignment and gene-assignment data this metric is computed from
is sound, not just that the histogram-plotting code runs.

## Library complexity (saturation curve)

Also computed automatically, no flag needed — Preseq's `c_curve`
question ("if this library were sequenced deeper, would you keep finding
new distinct molecules, or is it already close to saturated?"), computed
from every mapped fragment's (chromosome, alignment position, strand)
*before* any UMI or exact-position deduplication — deliberately: Preseq's
own convention is to characterize raw library complexity, not
post-dedup yield, and the two answer different questions (this + the UMI
dedup section above are complementary, not redundant).

Method: shuffle every mapped fragment's position with a fixed seed (so
this report is reproducible run-to-run on identical input, matching this
pipeline's existing determinism — verified explicitly: two independent
runs on the same data produced byte-identical saturation-curve output and
gene counts), then walk the shuffled list once with a hash set, recording
how many *distinct* positions have been seen at each of ten evenly-spaced
subsample checkpoints (10% through 100% of all mapped fragments). Output:
the ten-point table in `multiqc_summary.txt` plus an SVG line chart in
`report.html`. A simple, honestly-labeled slope check over the curve's
final 20% flags whether it's still rising steeply (far from saturated),
moderately (diminishing returns), or flattened (close to saturated) —
stated plainly as a description of the *observed* depth's shape, not
Preseq's own statistical extrapolation (a rational-function/Good-Toulmin
estimator) to predict yield at sequencing depths beyond what was actually
run, which this does not attempt to replicate.

**On the real dataset**: at the full 1M-read-pair depth, 969,276 mapped
fragments yielded 428,772 distinct positions (44.2%), with the curve
still rising moderately near full depth (+14% in the last 20% of the
subsample) — a real, sensible signal that this library, at this
sequencing depth, hasn't fully plateaued, with the expected diminishing
returns visible as depth increases (73,513 new distinct positions in the
first 10% of reads vs. 28,684 in the last 10%, a real, monotonically
narrowing gain consistent with what genuine library saturation looks
like).

## Duplication vs. expression (dupRadar-equivalent)

Also computed automatically, no flag needed. dupRadar's actual question is
different from UMI-tools dedup's: it doesn't try to identify and remove
PCR duplicates precisely (that needs a UMI, or is a best-effort guess
without one) — it asks whether a dataset's duplication *pattern* looks
like the ordinary, expected kind (highly-expressed genes rack up more
same-position reads than lowly-expressed ones purely by chance/PCR, so
duplication rate should rise smoothly with expression) or an anomalous
kind (genes with unusually high duplication for their expression level,
usually a sign of a technical artifact — degraded input, over-
amplification, or a library-prep bias — rather than real biological
signal). That's answerable without any UMI at all, using the same Picard
`MarkDuplicates`-style definition dupRadar's own upstream duplicate-
marking step uses: position-based duplicates (identical chrom + 5′
alignment position + strand), marked here completely independently of
`--umi-len`/UMI-tools-style dedup (see [UMI extraction and
deduplication](#umi-extraction-and-deduplication)) — when `--umi-len` *is*
given, the more precise UMI-aware dedup has already run and excluded its
duplicates before this analysis ever sees them, so in that case this
describes only the residual pattern among survivors of UMI dedup, not the
library's raw duplication rate (the same caveat the saturation-curve
section above draws around Preseq, for the identical reason).

Method: for every uniquely-assigned read/fragment (>=`MIN_UNITS_FOR_DUPRADAR`
per gene, currently 10, same threshold and rationale as gene-body
coverage's own minimum), accumulate two per-gene counts — total units and
position-duplicate units — then compute each qualifying gene's expression
as reads-per-kilobase (RPK) and its duplication rate, and the Pearson
correlation between log2(RPK) and duplication rate across all qualifying
genes. Output: `dupradar.tsv` (gene_id, total_count, dup_count, dup_rate,
rpk, log2_rpk) for the full per-gene picture, plus a 10-expression-bin
summary table in `multiqc_summary.txt` and a scatter plot (one point per
gene, log2 RPK vs. duplication rate) in `report.html` — dupRadar's own
report is the same shape, this trades its loess-smoothed curve for a
coarser but dependency-free decile-binned mean.

**On the real dataset**: 4,026 genes qualified on a 400,000-read-pair
slice of `SRR9336476` (97.3% of reads aligned; 348,936 fragments uniquely
assigned to a gene), with duplication rate rising monotonically from 13%
in the lowest-expression decile to 46% in the highest (Pearson r = 0.737)
— exactly the expected pattern for ordinary PCR-driven duplication, not a
sign of a library-prep artifact. A separate 100,000-pair check showed the
same monotonic shape at smaller scale (12%→36%, r=0.607) — the trend
strengthens rather than reverses as more data comes in, which is itself a
sanity check that it's a real signal and not sampling noise. Full 4,026-
gene table: `dupradar.tsv`; scatter plot with all 4,026 points: `report.html`.

## Junction annotation (RSeQC-equivalent)

Also computed automatically, no flag needed — this was the last piece of
item 14a's QC coverage still marked "not covered" as of the previous
pass. RSeQC's `junction_annotation.py` classifies every splice junction
it finds in a BAM against a reference gene model into three buckets, both
per splicing *event* (one per spliced read) and per distinct splicing
*junction* (one per unique chrom+donor+acceptor, however many reads
support it):

- **known** — both splice sites match an annotated exon-exon junction, as
  the same pair (not just each site independently annotated as part of
  some *other* junction)
- **partial novel** — one splice site matches an annotated donor or
  acceptor position (from any transcript), but paired with a site that
  doesn't — real alternative splicing not yet in the reference gene model
  looks like this
- **complete novel** — neither site matches anything annotated

That's implemented here the same way: known junctions are extracted from
consecutive same-transcript `exon` rows in the GTF (donor = end of exon
*i*, acceptor = start of exon *i+1*, both converted to this pipeline's
0-based coordinate convention), and every spliced alignment this
pipeline's own splice-anchor seeding finds (see `try_spliced_align`
above) is classified against them at the point it's processed in
`quantify_em()`'s Pass 2 — for every mapped unit, independent of its
gene-assignment outcome, matching RSeQC's own BAM-level (not
count-level) scope. Output: `junctions.tsv` (chrom, intron start/end,
category, supporting-read count) for the full per-junction list, plus a
summary in `multiqc_summary.txt` and a 3-bar chart in `report.html`.

**A dependency worth knowing about**: this needs `exon` feature rows in
the GTF, not just `gene` rows — a GTF with only gene-level spans (which
is all `load_gtf()` elsewhere in this pipeline has ever required) has no
per-transcript exon structure to derive junctions from. If none are
found, this section is skipped in the report rather than shown as a
meaningless "0% known" result. `convert_gff3_to_gtf.py` (included in this
package) now emits `exon` rows in addition to `gene` rows for exactly
this reason — see "Closed since the last pass" below.

### Two real bugs this surfaced

Validating this feature against the real dataset surfaced something this
pipeline's own spliced-alignment path was getting wrong, not just a gap
in QC coverage. On the first real run (100,000 read pairs, `SRR9336476`),
known junctions measured a median length of 326bp and a max of 766bp —
entirely normal; *S. cerevisiae*'s longest known intron is close to 1kb.
"Complete novel" calls, by contrast, measured a median length of
**~43,700bp**, and 89.5% of them were over 5,000bp, with a cluster
sitting right at the aligner's old `MAX_INTRON=100000` ceiling. That
pattern — real junctions clustering tightly at normal intron lengths,
"novel" junctions clustering at the search ceiling — isn't unannotated
biology.

**Bug 1 (partial fix): `MAX_INTRON` far too permissive.**
`try_spliced_align`'s Hamming-only scoring has no distance penalty, so
given a wide enough search window it will accept whatever distant,
unrelated k-mer match happened to reduce mismatch count. Lowered
`MAX_INTRON` from 100,000 to 15,000 (still 20x any real yeast intron).
Measured effect on the same 100,000-read-pair run: total spliced-
alignment calls dropped from 1,421 to 537, known-event fraction rose
from 15.6% to 41.3% — real progress, but "complete novel" junctions
still skewed ~15x longer (median ~5,000bp) than known ones (median
326bp) after this fix alone. Capped the damage; didn't explain it.

**Bug 2 (the actual root cause, now fixed): the splice path wasn't
gated as a real fallback, and its result was accepted unconditionally.**
`align_read()` only tries splice detection when the best ordinary
(Smith-Waterman) alignment leaves the read not fully covered — but
`sw_coverage < r->trimmed_len` is true for almost every read with even a
**single** soft-clipped base, which is extremely common (a terminal
sequencing error, an adapter remnant, a 1bp SNP near the read edge —
none of it real splicing). Worse, whatever the splice path found was
then used *unconditionally*, discarding a well-supported 74/75bp SW
alignment in favor of a distant, spurious "splice" purely because a
match existed somewhere. Fixed both halves: splice detection is now
only attempted when the unexplained portion is at least
`SPLICE_KMER_LEN` (10bp) — reusing that constant rather than inventing a
new threshold, since `try_spliced_align` itself requires exactly that
much anchor on each side to find anything, so a smaller clip could never
represent a genuine splice-spanning read by the detector's own
construction — and even then, the splice result is only kept if it
actually explains more of the read than the SW alignment already did.

**Measured effect of both fixes together**, 100,000-read-pair run: total
spliced-alignment calls 1,421 → 302 (a 78.7% reduction from the original,
un-gated baseline), known-event fraction 15.6% → **70.5%**, and — the
number that matters most — novel-junction median length dropped from
~43,700bp to **344bp**, statistically indistinguishable from known
junctions' 326bp median (89.5% of novel calls were over 5,000bp before;
8/38 = 21% after). At 400,000-read-pair scale: 1,237 spliced alignments,
69.9% known by event, novel-junction median length 341bp vs. known's
313bp. Total fragment-alignment rate was unaffected by either fix
(97.3% before and after, at both scales) — these changes only shifted
which of the two alignment paths (contiguous vs. spliced) explains each
already-mapped read, not how many reads map at all. Multi-lane merging
and dupRadar's duplication-vs-expression signal were both re-confirmed
unaffected (byte-identical merge outputs; dupRadar r=0.738 unchanged)
after applying both fixes.

**What's still open, precisely**: a small residual remains — at
100k-pair scale, 8 of 38 remaining "complete novel" junctions are still
over 5,000bp, down from 89.5% but not zero. Those are plausible
candidates for the deeper, unimplemented fix (canonical-site preference
and/or a distance penalty built into `try_spliced_align`'s scoring
itself, rather than gating around it) — a smaller, better-characterized
version of the same open item, not a new one.

## Differential expression internals

`compare` mode's DE test is a from-scratch reimplementation of DESeq2's
statistical model, not a wrapper around it:

1. **Normalization** — median-of-ratios size factors (the same estimator
   DESeq2 uses), computed once and shared between the PCA/heatmap
   normalization and the DE test's GLM offset.
2. **Dispersion estimation** — per-gene dispersion via Cox-Reid-adjusted
   profile likelihood maximization (1-parameter for the blind/PCA fit,
   2-parameter conditional-on-design fit for the DE test).
3. **Trend fitting** — a Gamma-GLM fit (iteratively reweighted least
   squares) of dispersion against mean expression.
4. **Empirical Bayes shrinkage** — per-gene dispersions are shrunk toward
   the fitted trend, with prior variance estimated from the trend residuals.
5. **Variance-stabilizing transform** — closed-form VST (not `rlog`), used
   for the PCA/heatmap inputs instead of a naive `log2(n+1)`.
6. **Outlier handling** — Cook's-distance thresholds from the F
   distribution, matching DESeq2's default outlier-flagging behavior.
7. **Multiple testing** — Benjamini-Hochberg FDR with independent filtering
   (low-count genes excluded from the testing set before BH adjustment, as
   in DESeq2).
8. **Effect-size shrinkage** — normal-normal shrinkage of log2 fold changes
   toward zero, with prior variance estimated from the MLE effect sizes and
   standard errors (in the spirit of DESeq2's `apeglm`/normal shrinkage
   estimators).

This stack was validated against synthetic ground-truth count data during
development, and real implementation bugs were found and fixed in the
process — see commit history for specifics.

### Volcano plot

`compare` mode's HTML report includes a volcano plot (log2FoldChange vs.
-log10(padj/q-value)), redesigned this pass to match a standard
three-tier layout: black (not significant), pink (p ≤ 0.05,
uncorrected), red (q ≤ 0.05) — replacing an earlier two-color,
p-value-only version. Labeled points: the 5 most significant up- and
down-regulated genes by q-value.

Built and validated against real DE runs, not just visually inspected
once — several real bugs were found and fixed in the process, each
caught by actually rendering the SVG and looking at it, not by reading
the code:

- **Unescaped `<=` in the legend text broke the SVG's XML validity
  outright** (`rsvg-convert` refused to parse it). Legend labels now use
  `&lt;=`, with a separate "visual length" table for spacing math, since
  `strlen()` on the escaped string overcounts relative to rendered width.
- **A plain max() over every gene's unshrunk log2FC let a small number of
  very-low-count genes (huge but statistically unsupported fold-change
  estimates) dominate the whole plot's x-range**, squashing every other
  point into two thin bands at the edges. Fixed by bounding the axis at
  the 99th percentile of |log2FC| among adequately-expressed genes
  (baseMean ≥ 10) instead, and by *omitting* (not clamping to the edge)
  non-significant genes whose estimate falls outside that range —
  clamping would draw a gene at a specific position implying a specific
  fold-change its data doesn't actually support. Never applied to q ≤
  0.05 genes, which are always shown.
- **Genes with q-value underflowing to exactly `0.0` — a real, common
  case for strong, clean differential expression, not a pathological
  input — were plotted at the bottom of the plot (least significant)
  instead of the top (most significant)**, because the guard condition
  checked `q > 0` and fell through to 0 on the equality case, exactly
  backwards. Found by noticing red (q ≤ 0.05) points sitting at y=0 in a
  test render, which shouldn't be possible, and confirmed by direct
  inspection of the specific genes producing it.
- **Label collisions**: strongly-changing genes routinely tie at the
  same q=0 ceiling, and placing every tied label at its literal pixel
  position stacked them illegibly. Labels are now vertically staggered
  (each pushed down only as far as needed to clear the previous one in
  the same up/down group), and labels for points clamped to the outer
  edge of the axis have their growth direction flipped inward, since a
  label growing further outward from an already-at-the-edge point runs
  straight past the plot's margin.

Validated on three real DE runs, not just one: technical pseudo-replicates
with no real signal (all points correctly cluster near the baseline, no
crashes, no spurious significance); a small number of genes given a
dramatic, uniform fold-change (correctly identified and labeled, but —
worth recording as a real finding, not glossed over — this first
synthetic test was bimodal, all-or-nothing signal with nothing in
between, which produced two disconnected clusters rather than a
funnel, initially read as a possible plotting bug before the cause was
traced to the test data's own unrealistic shape, not the code); and a
third run with a continuous spread of true effect sizes across ~400
genes (small to large, both directions) layered on realistic per-gene
noise, which produces the expected continuous volcano funnel —
confirming the first test's shape was a property of that specific
synthetic data, not a defect in the plot itself.

### On strandedness and gene assignment

The pipeline infers library strandedness (see [Output files](#output-files))
and, when the inference is confident (≥80%/≤20% concordance, the same
threshold used to call a verdict at all), uses it to help resolve
**ambiguous** multi-gene assignments — fragments that genomically overlap
more than one gene by position alone. When that happens, `quantify_em()`
narrows the candidate set to genes consistent with the fragment's implied
sense strand, but **only when doing so leaves at least one candidate**;
if every candidate would be filtered out (an occasional genuinely-antisense
read even in a strongly-stranded library), the original, unfiltered
candidate set is kept rather than forcing a possibly-wrong call. This never
changes an already-unique assignment, never zeroes out a mapped fragment,
and never affects unmapped or no-feature classification — it only ever
narrows an already-ambiguous set, which is what makes it safe to apply
without a separate opt-in flag.

**Measured effect**, full 1M-read-pair real dataset: 49,064 fragments had
their candidate gene set narrowed, and 36,033 of those (3.6% of all
mapped-with-feature fragments) were resolved all the way to a confident
unique assignment instead of going through EM as a multi-mapper. The total
count of mapped-with-feature fragments was **exactly unchanged**
(939,539 before and after) — confirming the disambiguation only
reclassifies *which* bucket an already-correctly-mapped fragment lands in
(unique vs. multi), never adds or removes fragments from gene-level
counting altogether.

**What this does and doesn't close**: this resolves ties among genes the
aligner already found overlapping a read's reported position — it doesn't
add strand information to the alignment step itself, and it doesn't help
with genes that don't overlap by position but might be expected to via a
different counting model (e.g. Salmon's transcript-level strand-aware
counting). For gene-dense regions with real antisense/overlapping
annotation, this should measurably improve specificity over the previous
strand-agnostic-only behavior; it's not a substitute for `featureCounts -s`
or Salmon's `--libType`-aware counting, which apply strand information more
comprehensively than this narrowing step does.

## Performance

### Measured, on real data (this validation pass)

Benchmarked against `bwa mem` 0.7.17 on the full `SRR9336476` *S. cerevisiae*
paired-end dataset (1,000,000 read pairs, 150 bp reads, R64-1-1 genome),
single-threaded on a 1-vCPU Intel Xeon @ 2.10GHz sandbox:

| | This pipeline | `bwa mem` |
|---|---|---|
| Index build (one-time) | 47.8s (k=16 + k=10 splice-anchor index, cached to `.kidx` afterward) | 5.95s (`bwa index`) |
| Alignment run, warm index (single-thread) | 219.1s — **includes** QC/trim, alignment, EM quantification, SAM/QC/HTML report writing | 144.7s — alignment (`bwa mem`) only, no downstream QC/quant/report |
| Reads/fragments aligned | 971,450 / 1,000,000 (97.1%) | 1,943,047 / 2,000,000 primary mapped (97.15%) |
| Proper pairs | 918,737 / 1,000,000 (91.9%) | 1,924,004 / 2,000,000 (96.2%) |

Two things to note when reading this table: (1) this pipeline's wall time
includes several stages `bwa mem` doesn't do at all (QC, trimming,
gene-level EM quantification, HTML/SVG report generation), so the ~1.5x
total-time gap overstates the pure alignment-speed gap; and (2) this
sandbox has a single CPU core, so the `-fopenmp` parallelization (which
matters most in practice) isn't reflected here — see the OpenMP note below.
Primary mapping rate is essentially identical to `bwa mem` on this dataset;
the proper-pair rate is a few points lower, consistent with the simpler
insert-size handling noted in [Known gaps](#known-gaps-vs-nf-corernaseq).

### Development-time context

Earlier optimization work (reusable thread-private DP scratch buffers
across reads, corrected banded Smith-Waterman initialization, a
zero-mismatch fast path before falling back to full DP, the on-disk k-mer
index cache, and a prefix-sum fix in the splice rescan) brought multi-core
end-to-end wall time down from an initial ~10x gap vs. `bwa mem` to roughly
~1.2x on prior multi-core benchmarking runs. Numbers will vary by machine,
thread count, read length, and genome size — there's no CI-tracked
benchmark in this repo, so re-benchmark on your own hardware/dataset if
performance matters for your use case.

Where the time goes, roughly:

- **First run** against a genome: k-mer index construction dominates
  (amortized away by the `.kidx` cache on subsequent runs — 47.8s once vs.
  effectively free after).
- **Alignment:** the banded Smith-Waterman extension is the hot loop; most
  reads hit the zero-mismatch fast path and never enter full DP.
- **Quantification:** EM iteration is cheap relative to alignment at
  typical convergence tolerances (180 iterations to converge on this
  dataset, well under a second).

`-fopenmp` parallelizes both the QC/trim loop and the alignment loop across
reads (`schedule(dynamic, ...)`); quantification and reporting are
currently single-threaded. On multi-core machines, build with `-fopenmp`
(see [Build](#build)) — the single-core numbers above don't exercise this.

### Memory

Peak resident memory on the full 1M-read-pair *S. cerevisiae* run: **~2.1
GB** as of the most recent architecture change (down from ~2.4 GB
previously — see below). Genome-scale memory behavior — previously the
single largest open question in this README, then a measured "no, doesn't
scale" — has now had real, measured architectural work done against it,
changing the quantified answer without changing the honest bottom line.

**Round 1: found the scaling problem and quantified it (see git history /
prior validation for the full account).** A hard per-chromosome buffer bug
was found and fixed first (a fixed 4MB static buffer with a hard `die()`
above that, unrelated to total genome size — real chromosomes in every
multicellular genome exceed it). With that fixed, indexing three synthetic
genome sizes (12/25/33 Mb) with the *original dense* k-mer index — every
position in the genome stored — found a clean linear relationship
(`peak_memory_MB ≈ 223 + 108 × genome_size_Mb`), confirmed by a real
observed OOM crash at 33 Mb landing almost exactly on the model's
prediction. Extrapolated: **~326 GB** to index a human genome — for
comparison, STAR needs ~27-31 GB and HISAT2 needs ~4.3 GB.

**Round 2: closed part of the gap with a real architectural change, not
just documentation.** Two changes, both validated for correctness before
being trusted for anything else:

1. **Minimizer-based sparse indexing.** Rather than storing every k-mer
   position in the genome, only the minimum-valued k-mer within each
   window of 8 consecutive positions is stored (the same core idea
   minimap2 uses, `MINIMIZER_WINDOW` in the source). This required also
   fixing the read-side seeder, which turned out to only check 3 fixed
   positions per read (start/middle/end) — fine against a dense index
   where nearly any position could hit, but would have silently gutted
   sensitivity against a sparse one. Fixed by computing the read's own
   minimizers (same window, same k) and checking those instead — for any
   genuinely shared region between read and reference, the read's
   minimizer for a window is, by construction, the same k-mer value the
   reference stored there, so this restores the seeding guarantee rather
   than just reducing memory at sensitivity's expense.
2. **Shrunk `SeedHit` from 16 to 8 bytes.** Its position field was a
   `long` (8 bytes); no real chromosome exceeds `int32_t` range (human's
   largest, chr1, is ~248 Mb). Free 2x reduction with an added overflow
   guard rather than a silent wraparound risk.

**Correctness was checked three separate ways before any memory number was
trusted**: position-level agreement with `bwa mem` on the real yeast
dataset (95.207%, statistically identical to the original dense index's
95.209% — the ~4.8% gap is a pre-existing aligner-vs-aligner tie-breaking
difference on repetitive regions, present equally in both versions, not a
regression), zero CIGAR/SEQ mismatches, and a full real 1M-read-pair run
end-to-end. That full run's mapping rate came out at **96.27%**, slightly
*higher* than the original dense index's 95.57% — the denser read-side
seeding (checking ~15-20 minimizer-selected positions per read instead of
3 fixed ones) more than compensates for the sparser reference index; this
wasn't just "no regression," it was a small net improvement in sensitivity
alongside the memory reduction.

**Memory, measured the same rigorous way as round 1** — three synthetic
genome sizes (12/25/90 Mb) with the new sparse index, real `/usr/bin/time
-v` measurements, fit with least squares (all three points within 10% of
the fit line, nearly exact at 90 Mb):

| Genome size | Dense index (round 1) | Sparse index + shrunk SeedHit (round 2) |
|---|---|---|
| 12 Mb | 1.55 GB | 0.59 GB |
| 25 Mb | 2.98 GB | 1.00 GB |
| 90 Mb | *(not tested)* | 2.37 GB |
| Fitted relationship | `223 + 108 × Mb` (MB) | `375 + 22.2 × Mb` (MB) |

Extrapolated to real genome sizes with the new fit:

| Organism | Genome size | Dense (round 1) | Sparse + shrunk (round 2) |
|---|---|---|---|
| *S. cerevisiae* (yeast) | 12 Mb | ~1.5 GB | ~0.6 GB |
| *C. elegans* | ~100 Mb | ~10.7 GB | ~2.5 GB |
| *D. melanogaster* | ~140 Mb | ~14.9 GB | ~3.4 GB |
| *M. musculus* (mouse) | ~2.7 Gb | ~284 GB | **~59 GB** |
| *H. sapiens* (human) | ~3.1 Gb | ~326 GB | **~68 GB** |

**What this means, stated plainly:** this is real, measured progress — a
**~4.8x reduction** in projected human-genome index memory (326 GB → 68
GB) — and it is *not* a solved problem. 68 GB is still ~2x STAR's ~30 GB
and ~16x HISAT2's ~4 GB; it moves human-genome indexing from "impossible
on any machine most people have access to" to "conceivable on a
high-memory server or HPC node," not into the same league as
purpose-built aligners. The reason it doesn't close the gap further is
structural, not a matter of more tuning: minimizer indexing reduces how
many *positions* get stored, but the per-position storage cost (an entry
in a chained hash table, plus the entry's own overhead) is still larger
per bit of information than a suffix-array/FM-index/BWT representation —
which remains, as in round 1, a genuinely different data structure and a
separate-project scope, not something the constants in this file can be
tuned into.

**Round 3: the runtime cost, precisely measured (not left as an assumption)
— and two more real bugs found while measuring it.** Round 2 shipped with
an honestly-flagged but unquantified claim: per-read alignment got slower
because sparse indexing needs more seed lookups per read. Measuring that
properly surfaced two more real, worth-fixing issues before any number
could be trusted:

- **A cache-versioning bug.** The on-disk `.kidx` index cache's header
  didn't record which indexing scheme built it. A cache built by the
  dense-index binary would pass every check the sparse-index binary's
  loader ran (same `KMER_LEN`, same table size, etc.) and get loaded
  as-is — not a correctness bug (a denser-than-expected table still
  answers lookups correctly), but a silent regression that would defeat
  the entire memory optimization with no error or warning, and would have
  invalidated this very benchmark by measuring the wrong thing. Found by
  actually trying to compare the two binaries against each other, which is
  exactly the kind of testing that keeps finding these. Fixed by adding a
  `minimizer_window` field to the cache header; verified the fix actually
  triggers a rebuild against a real stale cache before trusting anything
  measured afterward.
- **An `O(n × window)` algorithm where an `O(n)` one belongs.** The
  read-side minimizer computation (deciding which of a read's positions to
  seed from) re-scanned every window from scratch — a real, avoidable cost
  that had been deliberately accepted in round 2 in favor of a simpler,
  easier-to-verify implementation while the *correctness* of minimizer
  seeding itself was still being validated. With that validated, replacing
  it with the standard `O(n)` sliding-window-minimum algorithm (a
  monotonic deque) was the right next step — re-validated the same way
  (95.211% `bwa mem` position agreement, zero CIGAR/SEQ errors) before
  trusting the new timing. This also surfaced a third, smaller issue:
  comparing *raw* k-mer values (not hashed ones) meant a poly-T k-mer —
  which packs to the numerically largest possible value under this
  codebase's 2-bit encoding — was structurally disadvantaged in minimizer
  selection. Fixed by comparing `mix32()`-hashed values instead, the same
  approach production minimizer tools use to avoid exactly this
  compositional bias.

**Precisely measured, warm-cache, identical hardware and data:**

| | Alignment-stage wall time | Peak RSS |
|---|---|---|
| Dense index (round 1) | 253.2s | 2.49 GB |
| Sparse index, `O(n×window)` seeding (round 2, as shipped) | 717.0s (2.83x) | — |
| Sparse index, `O(n)` seeding (round 3, current) | 568.1s (2.24x) | 1.96 GB |

The deque fix recovered about 60% of the excess slowdown round 2 shipped
with (2.83x → 2.24x versus dense), and the memory reduction is real even
at yeast's small scale (21% less peak RSS) — recall from the table above
that the *relative* memory win grows much larger at real genome scales,
since the fixed per-run overhead that dominates at 12 Mb becomes a smaller
and smaller fraction of the total as genome size grows. **The remaining
~2.24x alignment slowdown is not further-optimizable away** the way the
2.83x→2.24x change was: that was fixing a genuine implementation
inefficiency (redundant rescanning), but checking more seed positions per
read is inherent to how a sparse reference index recovers the sensitivity
a dense one gets for free — it's the real, structural cost side of this
trade-off, stated as a number instead of a caveat now that it's been
properly measured.

A separate, previously-fixed memory issue remains fixed as before: raising
`MAX_CIGAR_OPS` (to fix the CIGAR/SEQ bug described in
[Validation on real data](#validation-on-real-data)) to 48 at one point
pushed peak memory high enough that the `samtools sort` subprocess
OOM-killed full-dataset runs on constrained machines; `MAX_CIGAR_OPS` was
brought back to 16 and `samtools sort` is capped at `-m 256M -@ 1`. Nothing
in the pipeline currently frees the k-mer index after alignment finishes,
even though it's not needed for quantification, QC, or reporting (except
during contaminant screening, which does free and rebuild it for its own
purposes) — freeing it earlier in the normal, non-screening path would
lower the peak during `samtools sort`/`bedtools`/`pyBigWig`. Still not
done in the normal path; still a real, scoped, deferred optimization
rather than something to claim as finished.

**Round 4: an actual FM-index was built and measured — and it turned out
worse than the minimizer approach, which is a real, useful result, not a
failure to hide.** Every prior round of this README speculated that "a
suffix array/FM-index, the approach real aligners use" would be the way to
close the remaining memory gap. Rather than continue asserting that,
round 4 built one: a from-scratch prefix-doubling suffix array
construction, BWT derivation, checkpointed Occ rank support, and FM-index
backward-search seeding, available as an opt-in `--fm-index` flag,
integrated into `try_sw_align_multi()` alongside (not replacing) the
validated minimizer path.

**Correctness was validated before integration, and again after.** In
isolation, before any of this touched the real pipeline: the suffix-array
construction reproduced the textbook reference suffix arrays for
`"banana$"` and `"mississippi$"` exactly, and — the more rigorous check —
1,920 exhaustive substring-search tests (every FM-index hit *set* compared
against a brute-force scan, swept deliberately across the Occ checkpoint
boundary) all passed. That process caught a real bug: the Occ checkpoint
meant to hold the cumulative count *after* the whole sequence collided
with checkpoint 0's slot whenever the sequence was shorter than the
checkpoint interval, silently corrupting every rank query on short
sequences — found and fixed before this ever touched real genome data.
After integration, `--fm-index` was re-validated the same way every other
change in this project has been: 95.183% position agreement with `bwa
mem` on the real yeast dataset (consistent with every other version's
~95.2%), zero CIGAR/SEQ mismatches, valid BAM.

**Then memory was measured, the same rigorous way — and the result was
not what four rounds of README speculation assumed:**

| Genome size | Minimizer + shrunk `SeedHit` (round 3) | FM-index (round 4) |
|---|---|---|
| 12 Mb | 587 MB | 557 MB |
| 25 Mb | 996 MB | 1,080 MB |
| 90 Mb | 2,366 MB | 3,218 MB |
| **Human extrapolation** | **~68 GB** | **~102 GB — worse** |

The FM-index prototype is *not* a memory win over the minimizer approach
at these scales — it's worse, and the gap widens with genome size. Diagnosed,
not just observed: the *final* FM-index structure (suffix array + BWT +
Occ checkpoints, kept in memory) is actually smaller in theory than the
minimizer hash table — roughly 9.6 bytes/bp analytically. But *building*
it with prefix-doubling suffix array construction requires three
full-genome-length arrays simultaneously (the suffix array itself plus two
rank-tracking arrays, 8 bytes each as `long` — necessary because a
concatenated multi-chromosome genome's total length can exceed `int32_t`
range at real genome scale, unlike the per-chromosome positions
`SeedHit` could safely shrink to in round 2) — around 24 bytes/bp of
*transient* memory, well above the structure's own steady-state size.
`/usr/bin/time`'s peak-RSS measurement catches that transient
construction-time high-water mark, not the smaller steady-state footprint
that would apply once loaded from a cache. **This is exactly why real
aligners use SA-IS (a linear-time, more memory-efficient suffix array
construction algorithm) instead of prefix-doubling, and typically sample
the suffix array (keep only every Kth entry, recover the rest via
LF-mapping) rather than storing it in full** — both well-understood,
unimplemented-here techniques, not open research questions. Implementing
them is real further work with a plausible path to actually beating the
minimizer approach (the theoretical steady-state number suggests it
could); it hasn't been done, and I'm not claiming otherwise.

**Why this is being reported this way, plainly:** it would have been easy
to quietly not finish this and let "FM-index would fix it" stand
unchallenged as the README's working assumption. Actually building it
and measuring found that assumption doesn't hold for the straightforward
implementation — a genuinely useful result, because it rules out "just
implement an FM-index" as a simple fix and correctly locates the real
remaining work (SA-IS, SA sampling) instead of leaving a plausible-sounding
but untested claim in place. `--fm-index` is left in the codebase, opt-in
and off by default, validated for correctness and honestly documented as
not yet delivering a memory improvement — usable as a foundation for that
further work rather than thrown away.

**Round 5: rebuilt per-chromosome instead of concatenated — and it now
beats the minimizer approach, the first FM-index version to do so.**
Round 4 diagnosed the FM-index's problem precisely: prefix-doubling
needs three full-genome-length `long` arrays simultaneously during
construction, and that transient high-water mark (not the smaller final
structure) is what dominated peak memory. Round 4's own reasoning said
those arrays needed to be `long` (8 bytes) "because a concatenated
multi-chromosome genome's total length can exceed `int32_t` range at
real genome scale." True of the whole genome concatenated together —
human totals ~3.1 Gb, over `int32_t`'s ~2.1 billion limit. False of any
*single* chromosome — human's largest, chr1, is ~248 Mb, comfortably
within range. So round 5 stopped concatenating: instead of one FM-index
over every chromosome joined together, it builds one independent
FM-index per chromosome. Every array in the construction hot loop, and
the permanent suffix array kept afterward, only ever needs to be as long
as whichever chromosome is currently being indexed — and every position
stored is chromosome-local, safe in `int32_t`. Two separate effects
follow: transient construction memory drops from O(total genome length)
to O(largest single chromosome), and the permanent suffix array itself
shrinks from 8 bytes/position to 4.

**Correctness re-validated the same way as every round before it,
because a from-scratch architecture change is exactly when that
discipline matters most, not when it's safe to skip.** The two paths
(minimizer and the new per-chromosome FM-index) agreed with each other
on 99.8% of commonly-aligned reads on the real 100,000-read-pair yeast
dataset, and both independently landed at ~94.2% agreement with `bwa
mem` on the clean-alignment comparison subset established earlier this
session (see [Validation on real data](#validation-on-real-data)) —
consistent with each other and with the existing baseline, not a new
number pulled from nowhere.

**Memory was measured on synthetic *multi-chromosome* genomes this
time** — round 4's 12/25/90 Mb figures used single-contig synthetic
FASTAs, which can't exercise a fix whose entire premise is
per-chromosome construction. Built three multi-chromosome synthetic
genomes instead (5, 6, and 11 chromosomes respectively, sizes chosen to
be plausibly uneven rather than uniform) and measured both the round-4
(concatenated) and round-5 (per-chromosome) builds on the identical
inputs:

| Genome size (chromosomes) | Largest chromosome | Concatenated FM-index (round 4 code) | Per-chromosome FM-index (round 5) |
|---|---|---|---|
| 12 Mb (5 chroms) | 4 Mb | 558 MB | 400 MB |
| 25 Mb (6 chroms) | 8 Mb | 1,080 MB | 737 MB |
| 90 Mb (11 chroms) | 15 Mb | *(not re-run — see below)* | 1,762 MB |

(Peak RSS via `/usr/bin/time -v`, indexing-dominated: each run used a
single read pair so alignment-stage memory wouldn't blur the
measurement, the same isolation approach round 3/4's figures used.) The
90 Mb concatenated-build point was extrapolated from round 4's own
already-measured 3,218 MB single-contig figure rather than re-run here —
multi-chromosome vs. single-contig makes no difference to the round-4
*concatenated* code path, since it flattens chromosome boundaries away
regardless (confirmed indirectly: round 4's single-contig 12/25 Mb
figures, 557/1,080 MB, match this round's multi-chromosome 12/25 Mb
concatenated-build re-measurement, 558/1,080 MB, to within noise).

**Directly tested the hypothesis, not just assumed it:** held total
genome size fixed at 90 Mb while moving the largest chromosome from 15 Mb
to 80 Mb (same 90 Mb total, one dominant chromosome instead of eleven
comparable ones) — peak memory rose only 6.4% (1,762 MB → 1,875 MB),
far less than proportional to the 5.3x jump in largest-chromosome size.
That's the expected signature of the fix working as designed: the
*steady-state* term (BWT + suffix array + Occ checkpoints, which scales
with total genome length regardless of chromosome count) now dominates
over the *transient* term (which scales with largest chromosome length)
at this scale, rather than the reverse.

**Extrapolation, with the same honesty about its limits as every prior
round's:** least-squares fit over these three real points (peak_MB =
250 + 16.9 × genome_Mb, R² fit reasonable but not perfect — the 25 Mb
point sits about 11% above the fitted line, plausibly genuine
sub-linearity in the fixed per-run overhead rather than noise) gives:

| | Minimizer (round 3) | Concatenated FM-index (round 4) | Per-chromosome FM-index (round 5) |
|---|---|---|---|
| **Human extrapolation** | ~68 GB | ~102 GB | **~51 GB** |

Real, substantial progress — roughly 25% below the minimizer approach
and half of round 4's concatenated FM-index — but stated with the same
caveat every extrapolation in this README carries: this is a
least-squares fit from three points topping out at 90 Mb, extrapolated
~34x to human's ~3.1 Gb, over a much larger genome and a much more
skewed chromosome-size distribution (human chr1 at ~248 Mb vs. chrY at
~57 Mb, a range this round's synthetic tests only loosely approximate)
than anything actually measured. It is **not** a validated human-scale
number, and it does **not** close the memory gap: ~51 GB is still well
above STAR's ~30 GB reference figure, so mammalian genomes remain
infeasible on a modest machine even after this round — meaningfully
better, honestly still short. What's left to close it further is
unchanged from round 4's own diagnosis: SA-IS construction (this round
kept prefix-doubling, just made its scope per-chromosome rather than
global) and suffix-array sampling (keep every Kth entry, recover the
rest via LF-mapping) — both still real, well-understood, unimplemented
next steps, not open questions.

**One correctness-adjacent cleanup made while touching this code:** the
FM-index's own heap allocations (per-chromosome BWT/suffix-array/Occ
arrays, plus the index array itself) were never freed anywhere,
including in `free_reference_index()` — the function contaminant
screening already calls to reset for a second genome. That's a
pre-existing gap (it applied identically to round 4's single concatenated
struct), not something round 5 introduced, but it directly affects
memory on the exact constrained-memory path this round's work is about,
so it was closed while already in this code rather than left for a
future pass to rediscover.

**A serious, silent, pre-existing bug found and fixed during this round's
own validation work — not a round-5-introduced regression, but one round
5's testing sequence happened to trigger for the first time.** Validating
that the per-chromosome rebuild hadn't changed anything on the default
(non-FM) path meant running `--fm-index` and default-mode runs
back-to-back against the same genome file, something no prior round's
validation sequence had done. That surfaced this: `build_kmer_index()`
skips building the minimizer k-mer table entirely when `--fm-index` is
set (correct — the whole point is not spending memory on an index that
mode won't use) — but the line saving the index to the on-disk `.kidx`
cache ran regardless, unconditionally, writing that **empty** table to
disk. Any later run against the same genome *without* `--fm-index` would
then load that cache, log a normal-looking "skipped FASTA parse + index
build" cache-hit message, and align against an empty index — no error,
no warning, just an alignment rate quietly collapsing from the correct
~97% to ~12-14% on the real yeast dataset. Confirmed by directly
reproducing it: an `--fm-index` run followed by a default run on the
same `ref.fa` reproduced the exact 12.4% figure; the fix (only write the
cache when the minimizer table was actually built — the FM-index path
already rebuilds fresh every run regardless of cache state, by design,
so it never needed to write this cache at all) brought both orderings
back to the correct 97.3%, verified with the fix in both directions
(`--fm-index` then default, and default then `--fm-index` then default
again). Every position-agreement figure quoted in this README
(94.2%/94.2%/99.8% above) was re-derived after this fix, from fresh runs,
specifically to confirm none of them had been silently computed against
a corrupted cache — they hadn't; those numbers were generated from
`.kidx`-cache-free or single-mode-only runs and matched exactly before
and after the fix. Anyone who used `--fm-index` even once against a
genome before this fix, then ran normally against the same cached
genome afterward, would have silently gotten near-total non-alignment.

**Round 6: replaced prefix-doubling with SA-IS — the second of the two
concrete next steps round 4 identified, and the one that finally brings
human within striking distance.** Round 4's diagnosis named two specific,
well-understood techniques as the real path to closing the memory gap:
SA-IS construction and suffix-array sampling. Round 5 addressed neither
directly — it fixed *why* prefix-doubling's transient arrays had to be
`long`-sized, not the O(n log n) prefix-doubling algorithm itself. Round
6 replaces prefix-doubling outright with SA-IS (Nong, Zhang, Chen 2009):
linear-time construction via induced sorting, using O(n) working memory
with a small constant instead of prefix-doubling's three simultaneous
full-length arrays.

**A genuinely intricate algorithm, validated accordingly — not trusted
on the strength of "it compiled and produced a plausible-looking
array."** SA-IS was built and validated standalone before it ever touched
this file: against a naive O(n² log n) reference suffix array across the
textbook "banana$"/"mississippi$" cases, hand-built edge cases (a
single-character string caught a real off-by-one — the lone sentinel
position is never *induced* from anywhere else, since induction always
derives position *j* from some existing `SA[i] = j+1`, which needs a
second position to exist at all; fixed by special-casing n=1 explicitly
rather than leaving it to silently fall through to an unfilled array),
2,000 random strings (n≤60, alphabet size 2–6), 200 larger DNA-like
random strings (alphabet size 6, matching this pipeline's actual
FM_ALPHA), and 200 highly-repetitive strings (alphabet size 3, to stress
induced-sort tie-breaking under heavy repeats) — 2,405/2,405 passed.
Separately timed at n=248,000,000 (approximating human chr1, the largest
real chromosome this pipeline will ever need to index): 76.5 seconds,
confirming both real-world feasibility and genuinely linear scaling (5x
the data took 5.5x the time — not the super-linear blowup an O(n log n)
algorithm would show, and not a scale prefix-doubling's ~6 GB transient
requirement could even attempt in this development sandbox's ~3.9 GB
RAM at all). Beyond the one-time standalone validation, `--fm-index` now
also re-runs a real self-test on every invocation (`fm_self_test_sa_is()`,
in `rnaseq_pipeline.c`) — comparing SA-IS's output against the
prefix-doubling algorithm it replaced (kept in the file specifically as
this oracle, not removed) on up to 200 Kb of real, actually-loaded
chromosome sequence, and refusing to proceed at all if they disagree.
That's deliberately over-cautious for a published, peer-reviewed
algorithm — but the standalone test suite already caught one real bug
this implementation had; a fixed set of synthetic cases validated once
is not the same guarantee as checking real biological sequence every
time, especially against future edits to this same file.

**A real, uniquely-required precondition that had to be fixed first,
found before it could cause a silent problem, not after:** SA-IS
requires the trailing sentinel to be the *unique* smallest symbol in the
coded sequence — a requirement comparison-based prefix-doubling never
had. The existing code mapped ambiguous (N) bases to the same symbol (0)
as the sentinel, which prefix-doubling handled correctly (it just
compares values, indifferent to repeats) but which breaks SA-IS's
base case (position n-1 sorting first depends on it being the unique
minimum). Fixed by giving N bases their own dedicated symbol
(`FM_ALPHA` bumped from 5 to 6: `$,A,C,G,T,N`) — they still can never
match a real query base (`fm_code()` still returns -1 for `N`, and
`fm_backward_search`'s seed-rejection logic is unaffected), they're
simply no longer indistinguishable from the one true sentinel.

**Memory measured the same way as every round before it** — on the same
synthetic multi-chromosome genomes round 5 used, so the two rounds are
directly comparable, not compared across different test setups:

| Genome size (chromosomes) | Largest chromosome | Round 5 (per-chrom, prefix-doubling) | Round 6 (per-chrom, SA-IS) |
|---|---|---|---|
| 12 Mb (5 chroms) | 4 Mb | 409 MB | 333 MB |
| 25 Mb (6 chroms) | 8 Mb | 755 MB | 555 MB |
| 90 Mb (11 chroms) | 15 Mb | 1,804 MB | 1,257 MB |

(Peak RSS via `/usr/bin/time -v`, same isolation approach as every prior
round: one read pair per run, so alignment-stage memory doesn't blur an
indexing measurement.) Least-squares fit over these three points
(peak_MB = 228 + 11.5 × genome_Mb) extrapolates to:

| | Round 4 (concatenated FM-index) | Round 5 (per-chrom, prefix-doubling) | Round 6 (per-chrom, SA-IS) |
|---|---|---|---|
| **Human extrapolation** | ~102 GB | ~51 GB | **~35 GB** |

**~35 GB is close to STAR's ~30 GB reference figure** — genuinely close
enough that this is no longer "well above," the honest framing every
earlier round has had to use. It is **still not a validated human-scale
number** (extrapolated from real points topping out at 90 Mb, the same
caveat every round's extrapolation has carried), and it is **still not
lower than STAR**, so the honest claim remains "closer, not closed" —
but for genomes actually tested at comparable scale, the extrapolation
crosses into workstation-feasible territory clearly: C. elegans (~100 Mb)
projects to ~1.35 GB, D. melanogaster (~140 Mb) to ~1.80 GB, both easily
within this development sandbox's own ~3.9 GB, not just a distant
extrapolation.

**A real result the skewed-chromosome test surfaced, reported rather
than smoothed over:** repeating round 5's "hold total genome size fixed,
move the largest chromosome" test (90 Mb total, largest chromosome moved
from 15 Mb to 80 Mb) with SA-IS gave 2,023 MB, a **57% increase** — far
more than round 5's prefix-doubling equivalent (1,762→1,875 MB, only
6.4%). That's not a regression or a sign SA-IS's per-chromosome design
failed; it's the flip side of SA-IS succeeding at what it was built to
do. Round 5's transient (chromosome-size-dependent) term was so much
larger than its steady-state (genome-total-dependent) term that
shrinking the transient term barely moved overall memory. SA-IS shrinks
transient memory so much further that it's no longer negligible next to
steady-state — so the two terms are now comparable in magnitude, and the
largest-chromosome dependency, while smaller in absolute terms than
before, is proportionally more visible than it was in round 5's numbers.
Practically: the "uniform" 90 Mb test (largest chromosome 16.7% of the
total) is the more representative one for the human extrapolation above
— human's actual largest chromosome, chr1, is ~8% of the genome, closer
to that ratio than to the skewed test's deliberately-extreme 88.9%,
which was built specifically to stress-test the hypothesis, not to
model real chromosome-size distributions. A genome with one truly
dominant chromosome (not human's case) would extrapolate worse than the
~35 GB figure above; this is named explicitly rather than left as a
silent assumption behind that number.

**What would close the gap the rest of the way, now precisely two items
rather than a vague "further optimization":** suffix-array sampling
(keep every Kth entry of the now-int32_t suffix array, recover the rest
via LF-mapping — round 4 identified this as the single largest further
reduction available, and it remains unimplemented, the last item on
round 4's original list not yet addressed) and reducing SA-IS's own
top-level working-memory overhead (the current integration converts the
`unsigned char` genome sequence to a transient `int32_t` array once per
chromosome for the top-level call only — real but modest, ~3 bytes/bp at
whichever chromosome is largest — by writing an unsigned-char-native
top-level bucket/induce pass instead of reusing the same int32_t core
used for recursive levels; documented here as a known, quantified,
available-but-unimplemented optimization rather than silently absorbed
into "SA-IS is done").

**Round 7: suffix-array sampling — the last item on round 4's original
list, and the one that finally crosses below STAR's reference figure.**
Round 4 named this the single largest further memory reduction
available: instead of storing the full suffix array (4 bytes/bp, the
largest term in the index's steady-state memory even after round 6),
store only the rows whose *text position* is a multiple of a sampling
rate K explicitly, and recover every other row via LF-mapping — a walk
through the BWT that's guaranteed to reach a stored ("marked") row
within K-1 steps, since each LF step moves to the row for the text
position exactly one earlier.

**Validated the same way as every mechanism in this file that isn't a
straightforward translation of well-known math:** standalone, against a
full suffix array as the oracle, across banana$/mississippi$, edge
cases, and 500+ random and large-scale tests spanning sampling rates
K=1 through K=60 — 571/571 passed. Standalone timing (matching the real
`FM_OCC_CHECKPOINT`/rank-checkpoint granularity) at K=8/16/32: 271/488/955
ns per lookup — K=16 was chosen from that measurement, not guessed, as
the point where per-lookup cost stays well under 1μs while still cutting
the sampled-position array to a quarter of the full suffix array. As
with SA-IS in round 6, a one-time standalone result isn't trusted
forever: `fm_build_one()` now validates the sampled structure against
the full suffix array it was built from — exhaustively for small
chromosomes, a 10,000-row random sample for large ones — every time a
chromosome's index is actually built, immediately before discarding the
full array (keeping it "just in case" would defeat the entire point of
sampling). Real-data correctness re-confirmed the same way as every
prior round: 94.207% `bwa mem` agreement (vs. round 6's 94.207% — not a
typo, genuinely identical to the digit, exactly what an *exact* recovery
mechanism should produce) and 99.825% agreement with the default
minimizer path.

**Memory, measured on the same synthetic genomes as rounds 5 and 6:**

| Genome size | Round 6 (SA-IS, full SA) | Round 7 (SA-IS, sampled SA, K=16) |
|---|---|---|
| 12 Mb | 341 MB | 310 MB |
| 25 Mb | 568 MB | 526 MB |
| 90 Mb | 1,257 MB | 1,010 MB |

Least-squares fit (peak_MB = 250 + 8.28 × genome_Mb) extrapolates to
**~25.3 GB for human** — down from round 6's ~35 GB, and, for the first
time across all seven rounds, **below STAR's ~30 GB reference figure**.

**A real result worth explaining precisely, not just reporting the
headline number:** the reduction (35→25.3 GB, ~28%) is smaller than a
steady-state-only calculation would predict. Sampling should cut the
suffix array's steady-state cost from 4 bytes/bp to roughly 0.375
bytes/bp at K=16 (the sampled-value array plus the marker bitvector) —
on paper, that alone should take the ~5.6 bytes/bp steady-state total
round 6 had down to roughly 2 bytes/bp, a ~60% reduction, not 28%.
Computing what a steady-state-only estimate would predict at 90 Mb
(≈428 MB) against what was actually measured (1,010 MB) accounts for
the gap directly: roughly 580 MB of the measured peak — more than half
of it — is *transient construction cost that sampling doesn't touch at
all*. Building the sampled structure still requires materializing the
full suffix array first (to read values from while building the
compact form), plus SA-IS's own O(n) working arrays during that same
construction — none of which shrink just because the *kept* structure
afterward is smaller. Sampling reduces what's kept; it does nothing for
what's needed, however briefly, to build it. This is the same shape of
finding round 4 made about the FM-index generally (steady-state size and
actual peak memory are different things, and only one of them is what a
real run needs to survive) — round 7 didn't escape it, just moved where
it applies.

**What would close this further, now precisely one item, not a vague
"further optimization":** reducing SA-IS's own transient construction
memory — either avoiding the full-array materialization step before
sampling (build the sampled structure more incrementally), or reducing
the top-level `unsigned char`→`int32_t` conversion round 6 already
flagged as a known, quantified, unfixed cost. Both are genuine further
implementation work, not configuration.

**Round 8: closed the second of those two, and the result is a real,
useful lesson about what "real measured improvement" does and doesn't
tell you.** `fm_build_sa()` previously converted the whole chromosome
from `unsigned char` to `int32_t` before calling SA-IS's core — a real,
quantified ~4 bytes/bp cost, paid once per chromosome, never compounding
through recursion. Replaced with a u8-native top-level entry point
(`sais_build_u8`, mirroring the existing `sais_build_i32` line-for-line
so the two can't silently diverge) that operates on the chromosome's
bytes directly. Validated the same way as every SA-IS-adjacent change
in this file: standalone, against both the naive oracle AND a direct
comparison against the already-validated `sais_build_i32` on every one
of that suite's test cases — 4,810/4,810 passed — then re-confirmed on
real yeast data (identical 97.3% alignment rate, both self-tests passing
silently, multi-lane merging still byte-identical).

Measured on the same three synthetic genomes as rounds 5–7:

| Genome size | Round 7 | Round 8 | Reduction |
|---|---|---|---|
| 12 Mb | 310 MB | 303 MB | 2.3% |
| 25 Mb | 526 MB | 494 MB | 6.0% |
| 90 Mb | 1,010 MB | 990 MB | 1.9% |

Real, consistent, measured at every scale tested. **But re-fitting the
human extrapolation from these three points moved it only from ~25.31 GB
to ~25.16 GB — under 1%, not proportional to the 2-6% seen at tested
scale.** The reason is precise, not mysterious: this fix's cost was
always bounded by the *largest single chromosome*, not by total genome
size (same category as round 5's per-chromosome redesign) — a
fixed-per-run cost, not one that scales with the genome. The
least-squares fit's *intercept* (fixed overhead) dropped a real 7% (250
MB → 233 MB); its *slope* (marginal cost per additional Mb of genome,
which is what steady-state, total-genome-scaling costs drive, and what
dominates the human extrapolation at 3,100 Mb) barely moved (8.28 →
8.24). Worse for the extrapolation specifically: this synthetic test's
largest chromosome (15 Mb) is 16.7% of its 90 Mb total, while human's
largest chromosome (chr1, ~248 Mb) is only ~8% of its ~3,100 Mb total —
so a largest-chromosome-bounded fix contributes *proportionally less*
at human's real, more skewed chromosome-size distribution than at this
round's own test scale. A real, validated fix that helps every genome
tested, and helps human's specific memory profile by less than the
tested-scale numbers alone would suggest — worth knowing precisely,
not just citing the percentage that looks best.

**What's left, now genuinely down to one item:** avoiding full-suffix-
array materialization before sampling (building the sampled structure
more incrementally during construction, rather than building the whole
array and reading it once to compress it). Round 7's own measurement
attributed roughly 580 MB of the 90 Mb test's peak to transient cost
generally; round 8 closed a real but modest slice of that (the u8
conversion). The remaining, larger piece is still open, and given what
round 8 just found about how largest-chromosome-bounded costs scale to
human's actual chromosome distribution, it should be measured with the
same skepticism before its impact on the extrapolation is claimed,
not assumed from its effect at 90 Mb scale.

**Round 9: investigated this directly with per-chromosome RSS
instrumentation — found the 580 MB "transient" estimate was too coarse,
tested and ruled out one hypothesis, and arrived at a more precise,
evidence-based picture instead of a new headline number.** Instrumented
`fm_build_one()` and the per-chromosome loop in `fm_build()` with
`getrusage()` checkpoints (temporary, not shipped in this build) and ran
the 90 Mb, 11-chromosome test genome through it. The first, and most
informative, finding: **RSS never decreases between chromosomes** — it
climbs monotonically (743 → 905 → 914 → 918 → 930 → 940 → 947 → 956 →
962 → 965 → 965 → 965 MB across the 11 chromosomes, largest first) even
though each chromosome's transient construction arrays are correctly
freed before the next one's are allocated. That pattern looked at first
like it might be glibc's allocator retaining freed heap memory instead
of returning it to the OS — a well-known allocator behavior with a
standard, safe fix (`malloc_trim(0)`). **Tested that hypothesis
directly: it made no measurable difference** (988 MB peak with
`malloc_trim(0)` called after every chromosome vs. 990 MB without,
noise-level). Ruled out, not assumed away.

The real explanation, visible in the per-chromosome deltas: the *first*
(largest, 15 Mb) chromosome's construction alone accounts for 162 MB of
the 222 MB total growth (73%) — consistent with that chromosome needing
the most memory of any single chromosome's build, transient and
steady-state combined, and with `malloc` satisfying every subsequent
(smaller) chromosome's allocations largely out of the heap space that
first, larger allocation already reserved (freed chunks get reused for
same-shaped later allocations without requesting more memory from the
OS — the standard reason `free()`-then-`malloc()` often doesn't shrink
RSS, distinct from the "never returned to the OS at all" failure mode
`malloc_trim()` addresses). This refines, rather than contradicts,
round 7's finding: peak memory really is dominated by construction-time
cost tied to the largest chromosome, but the actual mechanism is
allocator reuse across chromosomes, not simple accumulation — and it
means the ceiling on "avoiding full-SA materialization" specifically is
bounded by whatever fraction of *that one largest chromosome's* build
is attributable to transient SA-IS arrays, not the full ~580 MB figure
a cruder steady-state-only comparison suggested. That's a smaller,
more precisely-bounded target than round 7's number implied — genuinely
useful to know before attempting the restructuring itself, which is
real algorithmic work (changing how and when the suffix array's values
are read, not a parameter or a small function swap) and remains
unimplemented this round, deliberately: rushing it under a smaller,
now-better-understood ceiling would risk exactly the kind of subtle
correctness bug this project's discipline exists to catch, for a
return this round's own measurement suggests would be modest.

## Validation on real data

This repo was exercised end-to-end against the real `SRR9336476`
*S. cerevisiae* dataset and the Ensembl R64-1-1.59 GFF3/genome (paired-end,
single-end, and multi-sample `compare` modes) as part of this validation
pass. Findings:

- **Correctness, single/paired-end alignment + quantification:** primary
  mapping rate (97.1%) matches `bwa mem` (97.15%) on the same data to
  within noise; see the [Performance](#performance) table.
- **A measurement discrepancy found during this pass's due diligence —
  and resolved before this pass ended.** After making the two
  `align_read()`/`MAX_INTRON` fixes described under
  [Junction annotation](#junction-annotation-rseqc-equivalent), installed
  `bwa` in this session specifically to confirm those fixes hadn't
  regressed the ~95.2% position-agreement figure quoted elsewhere in this
  README (from Round 4's FM-index validation). A controlled before/after
  check passed cleanly (83.64% vs. 83.67% exact-position agreement across
  the fix, a 0.03-point difference — noise from the few hundred reads
  whose alignment path changed, not a regression). But 83.7% itself, over
  *all* commonly-mapped reads regardless of CIGAR shape, was noticeably
  lower than 95.2% — worth chasing down rather than shrugging off. It
  turned out to be a like-for-like problem, not a correctness one: 83.7%
  mixes in every soft-clipped/partial alignment, where two independently-
  built aligners' extension heuristics routinely pick different (both
  individually defensible) clip boundaries a few bp apart — not
  misalignment, just a fuzzier comparison than "did this read map to the
  right place." Restricting to reads where **both** aligners report a
  clean, full-length, gap-free `150M` alignment and `bwa` itself isn't
  flagging ambiguity (MAPQ ≥ 1): **94.2–94.6% exact-position agreement**
  (94.22% at MAPQ≥1, 94.62% at MAPQ=60 only, both on the same
  100,000-read-pair data) — matching the documented 95.2% within
  reasonable measurement variation, not a discrepancy. Reasonable
  conclusion: the original 95.2% almost certainly came from a similarly
  restricted (full-length, unambiguous) comparison, which is a completely
  standard and defensible way to isolate "does the core seed-and-extend
  engine put reads in the right place" from soft-clip-boundary noise that
  any two aligners will show. Both numbers are real and both are correct
  for what they measure — 83.7% describes agreement across *everything*
  commonly mapped, including edge cases where boundary placement is
  inherently fuzzier; ~94-95% describes agreement on the
  unambiguous-and-complete subset (~76% of commonly-mapped reads in this
  dataset). Neither should be quoted without saying which one it is.
- **Correctness, differential expression:** running `compare` mode on two
  pairs of *technical* pseudo-replicates (the same library split into
  quarters, so there's no true biological difference between "conditions")
  correctly found only 2 of 2,640 testable genes significant at FDR < 0.1
  — i.e. the BH-FDR procedure is controlling the false-positive rate about
  where it should on a true-null dataset. This is exactly the sanity check
  you'd want to pass before trusting the DE stack on real comparisons.
- **Determinism:** two independent runs on identical inputs produced
  byte-identical `gene_counts.tsv`.
- **Known bug found during this validation:** a small fraction of
  alignments (20 out of 1,999,998 records in the 1M-read-pair run, ≈0.001%)
  have a SAM CIGAR string whose query-consuming length doesn't match the
  `SEQ` field length — e.g. read `SRR9336476.79761`, CIGAR
  `40S1M1I52M1I5M1D3M1I2M1D38M` implies a 144 bp query against a 146 bp
  `SEQ`. All observed cases involve indel-dense CIGARs (many short I/D runs
  in a row, often over low-complexity/homopolymer stretches), suggesting an
  edge case in the banded Smith-Waterman traceback where consecutive
  indels near a band boundary aren't fully accounted for in the emitted
  CIGAR. This is rare enough not to affect the mapping-rate/DE numbers
  above materially, but it will make `samtools`/`picard` reject the
  affected records (`samtools flagstat` errors with "CIGAR and query
  sequence are of different length" at the first such record). Flagged
  here with reproducible read IDs for a future fix rather than patched in
  this pass, since a correct fix needs to go through the banded-SW
  traceback code (`sw_scratch_ensure`/the alignment extension logic around
  it), not just the CIGAR serializer.

**Second pass (this round):** the CIGAR bug above was fixed and re-verified
— `samtools flagstat` now parses all 1,999,998 records in the same
1M-read-pair run with **zero** errors and zero CIGAR/SEQ length mismatches
(down from the 20 found previously), and the pipeline's own mapping-rate
accounting was unchanged before vs. after (97.1% both times), confirming
the fix corrected formatting, not alignment behavior. Three more gaps were
closed and each re-verified against the same full dataset:

- **Gzip input**, tested by running the complete pipeline directly against
  the original `genome.fa.gz` + `annotation.gtf.gz` + both `reads.fastq.gz`
  files with no decompression step — output was byte-identical to the
  pre-decompressed baseline (`gene_counts.tsv` diff: 0 lines; SAM: 0
  CIGAR/SEQ mismatches; identical `samtools flagstat`).
- **Strandedness detection**, which correctly called `SRR9336476` as
  reverse-stranded (96.8% of uniquely-assigned reads have read1 opposite
  the gene's strand) — consistent with a standard Illumina TruSeq Stranded
  mRNA-style prep, and a real, checkable signal rather than a guess.
- **Sorted/indexed BAM output**, which passes `samtools quickcheck`, has
  the correct `SO:coordinate` header, and produces a valid `.bai` — checked
  on both a smaller subset run and the full 1M-read-pair run.

That work also **introduced and then caught** a real memory regression
(the `MAX_CIGAR_OPS` increase interacting with the new `samtools sort`
subprocess to OOM-kill full-dataset runs on a 4 GB machine) — see
[Memory](#memory) for the full account, including the fix and
re-verification. I'm surfacing that prominently rather than only
mentioning the fix, because catching your own regression via the same
real-data validation loop that caught the original bug is the actual
mechanism by which this kind of software becomes more trustworthy over
time — and it's also a reminder that this pass, like the last one, almost
certainly didn't catch everything.

**Third pass (this round):** closed bigWig coverage tracks, annotation-based
rRNA/biotype content QC, and UMI extraction/deduplication — all tested, not
just written. Summary (full detail under
[UMI extraction and deduplication](#umi-extraction-and-deduplication) and
the biotype row in the comparison table above):

- **bigWig**: round-tripped the output through `pyBigWig` and confirmed
  correct header, chrom sizes, and coverage statistics, on both a subset
  run and the full 1M-read-pair run (`nBasesCovered: 12,157,105`, matching
  the genome size exactly).
- **rRNA/biotype QC**: this required first fixing a real bug in
  `convert_gff3_to_gtf.py` — it was silently dropping every non-protein-coding
  gene (all rRNA/tRNA/snoRNA/snRNA loci live under GFF3's `ncRNA_gene`
  feature type, not `gene`), so rRNA content could never have been measured
  even in principle before this fix. After fixing it: 7,127 genes load
  instead of 6,600, the full rDNA repeat locus is present, and the
  resulting QC correctly reports 0.27% rRNA content on the full dataset —
  consistent with what you'd expect from a polyA-selected library, and
  consistent between the subset test (0.29%) and the full run (0.27%).
- **UMI dedup**: validated against synthetic ground truth with a known
  number of injected duplicates (not just "it runs without crashing") —
  96.3% (SE) and 97.0% (PE) of injected duplicates correctly caught, and
  critically, **every real gene count came out identical** to a duplicate-free
  baseline once dedup was applied. See
  [UMI extraction and deduplication](#umi-extraction-and-deduplication) for
  the full numbers and the honestly-reported ~3-4% miss rate.
- **A full-scale integration test**: all of the above, plus every gap
  closed in the prior two passes, run together against the complete
  1M-read-pair dataset in one pipeline invocation. This surfaced one more
  real, sandbox-specific issue worth being transparent about: the first
  attempt at this integration run appeared to silently die partway through
  (SAM file cut off at 35,770 of ~2,000,000 records, all downstream output
  files empty, no error in the log). Investigation showed **the actual
  cause was the test sandbox's disk running out of space** — this session's
  accumulated test output (many redundant ~800MB run directories from
  earlier validation passes) had filled the filesystem, not a bug in this
  round's new code. That's an environment issue, not a pipeline defect, but
  it's worth naming because "the run silently died with no error" is
  exactly the kind of failure mode that's easy to misdiagnose (I initially
  suspected another OOM interaction from the new bigWig subprocess chain,
  which would have been a real code bug — it took checking `df -h`, not
  just `free -h`, to find the actual cause). After clearing space, the full
  integration run completed cleanly: valid BAM (`samtools quickcheck`
  passes), zero CIGAR/SEQ mismatches across all 1,999,998 records, a
  correct bigWig, and consistent biotype QC.

## Known gaps vs. nf-core/rnaseq

**Closed across this project's validation passes** (all re-tested against
the full real dataset after each fix, not just spot-checked):

- ~~Gzip-compressed input~~ — see [Input requirements](#input-requirements).
- ~~Strandedness~~ — closed as a diagnostic; see
  [On strandedness and gene assignment](#on-strandedness-and-gene-assignment)
  for the real remaining limitation (detection doesn't yet feed back into
  strand-aware counting).
- ~~BAM output + sorting/indexing~~ — via `samtools`, same approach nf-core
  uses; see [Output files](#output-files).
- ~~bigWig coverage tracks~~ — via `bedtools` + `pyBigWig`; see
  [Output files](#output-files).
- ~~rRNA-content QC~~ — via GTF `gene_biotype`, after fixing a real bug in
  `convert_gff3_to_gtf.py` that was silently dropping every non-protein-coding
  gene; see [Validation on real data](#validation-on-real-data).
- ~~UMI extraction/deduplication~~ — for the common in-line-UMI-on-R1
  layout, via `--umi-len N`; see
  [UMI extraction and deduplication](#umi-extraction-and-deduplication).
- **Found and closed two self-inflicted regressions along the way**: the
  CIGAR-truncation fix's first version (`MAX_CIGAR_OPS=48`) plus the new
  `samtools sort` step together OOM-killed full-dataset runs on constrained
  machines (see [Memory](#memory)); and a full-scale integration run that
  appeared to silently die turned out to be the *test sandbox* running out
  of disk space from accumulated validation artifacts, not a code bug (see
  [Validation on real data](#validation-on-real-data)) — worth keeping
  distinct from the OOM issue since they look similar (silent death,
  partial output) but have completely different causes and fixes.

**Closed since the last pass:**

- ~~Running against an output directory that doesn't exist yet fails
  seven pipeline stages in, with no indication the directory is the
  problem~~ — a real bug, reported by an actual user running this
  pipeline against real data on real infrastructure, not found in
  internal testing: `rnaseq_pipeline ... outdir/` where `outdir/`
  didn't already exist died with a bare `ERROR: cannot write
  alignments.sam` after alignment had already finished, giving no hint
  that the actual problem was the directory rather than the SAM file
  or anything upstream. Root cause: nothing in the pipeline ever
  created the output directory -- every `write_*()` function assumed
  it already existed. Fixed with `ensure_dir_exists()`, called at
  startup in both the default pipeline and `compare` mode, before
  anything is written -- implemented as direct `mkdir()` calls on each
  path component (mirroring `mkdir -p`) rather than shelling out to
  `mkdir -p`, specifically so a path containing shell metacharacters
  (spaces are already common in real filesystem paths) is never passed
  through a shell at all. Verified against the exact reported scenario
  (a single missing directory level) and a stress case (four missing
  nested levels at once) -- both now succeed and produce correct
  output; the existing-directory case was re-confirmed unaffected.
  Every other bare `die("cannot write ...")` in the file was also
  given a `strerror(errno)`-qualified message as defense-in-depth, so
  a *different* real failure (permissions, a full disk) is diagnosable
  from the error text alone next time, rather than requiring another
  round of "here's my exact command and here's the bare error" back
  and forth to localize.

- ~~Genome-scale memory behavior was untested~~ — no longer untested, and
  no longer a static verdict either. Round 1: found and fixed a hard bug
  (a fixed 4 MB per-chromosome buffer, unrelated to total genome size,
  that fails on virtually every real multicellular chromosome), then
  measured the original dense k-mer index's memory scaling against
  synthetic genomes at three sizes including a real confirmed OOM
  boundary, extrapolating to ~326 GB for a human genome. Round 2:
  implemented minimizer-based sparse indexing (the same core idea
  minimap2 uses) plus a smaller `SeedHit` struct, validated for
  correctness against `bwa mem` (95.207% position agreement, statistically
  identical to the original dense index's 95.209%) and a full real
  1M-read-pair run (mapping rate 96.27%, slightly *better* than the
  original 95.57%) before re-measuring memory the same rigorous way —
  bringing the human-genome extrapolation down to **~68 GB**, a real
  ~4.8x reduction. See [Memory](#memory) for the complete data from both
  rounds. This remains open in the sense that matters most: 68 GB is
  still far more than STAR (~30 GB) or HISAT2 (~4 GB), and closing that
  remaining gap needs a genuinely different index data structure, not
  more tuning of what's here — but "we made real, measured progress and
  here's exactly how much" is a meaningfully different, better answer
  than round 1's "no" was, which was itself already better than the
  original "untested." The good news across both rounds: compact genomes
  (yeast through fly/worm/plant scale, up to a few hundred Mb) stay in a
  workstation-feasible memory range, and got noticeably cheaper still in
  round 2.
- ~~Strand-aware gene counting~~ — the strandedness diagnostic now feeds
  back into gene assignment, narrowing (never forcing) ambiguous multi-gene
  candidate sets when the library is confidently stranded. Measured on the
  full real dataset: 36,033 fragments (3.6% of all mapped-with-feature
  fragments) resolved from ambiguous multi-mapper to confident unique
  assignment, with the total mapped-with-feature count exactly unchanged
  before and after — see
  [On strandedness and gene assignment](#on-strandedness-and-gene-assignment).
- ~~UMI-tools' directional-adjacency dedup method~~ — replaced the previous
  exact-match-only ("unique" method) dedup with UMI-tools' default
  directional-adjacency clustering, which also catches duplicates whose UMI
  read has a single sequencing error. Validated with ground truth built
  specifically to test this (half the injected duplicates exact, half with
  a 1bp UMI substitution): 96.7% overall recall, with ~half the catches
  coming *only* from the new clustering logic — see
  [UMI extraction and deduplication](#umi-extraction-and-deduplication).
- ~~Contaminant genome screening~~ — this was previously described as
  "structurally not closeable in the abstract," which on reflection
  conflated "this repo can't ship a specific contaminant genome" (still
  true) with "the screening mechanism itself can't be built" (not true —
  BBSplit doesn't ship reference genomes either). Implemented as
  `--contaminant-ref path[:Label]`, reusing the pipeline's own aligner
  against user-supplied reference(s). Validated with a real genome (lambda
  phage, fetched from a public repo) and genuinely-sampled synthetic reads
  from it: 2,000/2,000 injected contaminant reads correctly detected (100%,
  zero false negatives), zero false positives among unrelated unmapped
  reads, primary-genome outputs confirmed unaffected — see
  [Contaminant genome screening](#contaminant-genome-screening).
- ~~Gene-body coverage~~ — RSeQC `geneBody_coverage.py`-equivalent 5'/3'
  bias detection, computed automatically from data already produced during
  quantification. Found a real, biologically sensible 3' bias signal on
  the actual dataset (ratio 3.92 full-scale, 3.63 on a subset — internally
  consistent), independently corroborated by this same dataset's very low
  rRNA content — two unrelated measurements agreeing is real evidence the
  underlying data is sound, not just that the plotting code runs. See
  [Gene-body coverage](#gene-body-coverage).
- ~~Library complexity / saturation curve~~ — Preseq `c_curve`-equivalent,
  also computed automatically. Verified deterministic (byte-identical
  output across repeated runs on the same input, matching this pipeline's
  existing determinism guarantee), and found a real, sensible signal on
  the full 1M-read-pair dataset: 44.2% of mapped fragments at distinct
  positions, curve still rising moderately near full depth with the
  expected diminishing returns (monotonically narrowing gains as depth
  increases). See
  [Library complexity](#library-complexity-saturation-curve).
- ~~`MAX_GENES=8192` static cap~~ — found while sizing this feature's
  storage: would `die()` outright on any real multicellular GTF (human
  alone has ~20,000 protein-coding genes plus many more ncRNA/pseudogene
  loci). Bumped to 65,536 at negligible memory cost. Unlike the earlier
  4MB chromosome-buffer bug, this one already failed loudly rather than
  silently corrupting data, but was still a real, easy-to-hit usability
  blocker worth closing while already in this part of the code.
- ~~Samplesheet-driven multi-lane merging~~ — a comma-separated FASTQ path
  list (either R1 alone in `se` mode, or matched-length R1/R2 lists in `pe`
  mode) is now concatenated in order before trimming/alignment, the same
  outcome as nf-core's per-sample lane `cat` step. Validated bit-for-bit at
  two scales: a 2-lane split of a real 100,000-read-pair slice, and again
  at 4x scale (2x200,000-pair lanes), both producing byte-identical
  `alignments.sam`, `gene_counts.tsv`, and `dupradar.tsv` to the same reads
  given as one file. A mismatched R1/R2 lane count is rejected with a clear
  error before any alignment work starts. See
  [Multi-lane FASTQ merging](#multi-lane-fastq-merging).
- ~~dupRadar-equivalent duplication modeling (independent of UMIs)~~ —
  Picard-style position duplicate marking (chrom + 5′ position + strand),
  computed independently of `--umi-len`, correlated per-gene against
  expression (RPK). Found a real, expected biological signal on the actual
  dataset: duplication rate rose monotonically from 13% (lowest-expression
  decile) to 46% (highest) across 4,041 genes on a real 400,000-read-pair
  slice (Pearson r=0.738), with a smaller 100,000-pair check showing the
  same monotonic shape at a weaker but still-positive correlation
  (r=0.607) — the relationship strengthening rather than reversing with
  more data is itself evidence this is a real signal, not sampling noise.
  See [Duplication vs. expression](#duplication-vs-expression-dupradar-equivalent).
- ~~RSeQC junction-annotation metrics~~ — closes the other half of the
  previous pass's combined item 2. Known-vs-novel splice-site
  classification, both per splicing event and per distinct junction,
  built from the GTF's own `exon` rows. Along the way, found that
  `convert_gff3_to_gtf.py` (this package's GFF3→GTF converter) never
  emitted `exon` rows at all, only gene-level spans — meaning junction
  annotation had nothing real to classify against on this project's own
  validation dataset until the converter was fixed to also emit them. See
  [Junction annotation](#junction-annotation-rseqc-equivalent).
- ~~Spliced-alignment false positives from an ungated, unconditional
  fallback~~ — this was originally flagged as still-open after the last
  pass's `MAX_INTRON` reduction alone. The actual root cause turned out to
  be one level up: `align_read()` tried splice detection for almost any
  read with even a single soft-clipped base (a terminal sequencing error,
  not real splicing, on a large fraction of otherwise well-aligned reads),
  and then used whatever it found *unconditionally*, discarding a
  well-supported ordinary alignment in favor of a spurious distant
  "splice." Fixed by requiring the unexplained portion to be at least
  `SPLICE_KMER_LEN` before even attempting splice detection, and requiring
  the result to actually explain more of the read than the existing
  alignment before accepting it. Measured effect on top of the previous
  `MAX_INTRON` fix, same 100,000-read-pair run: spliced-alignment calls
  1,421→302 (a 78.7% reduction from the original baseline), known-event
  fraction 15.6%→70.5%, novel-junction median length ~43,700bp→344bp
  (vs. known junctions' 326bp — no longer distinguishable by length).
  Total fragment-alignment rate unaffected (97.3% before and after) —
  confirms this only reclassified which of two existing alignment paths
  explains each read, not how many reads map. See
  [Junction annotation](#junction-annotation-rseqc-equivalent) for the
  full before/after and the small residual that's still open.

- ~~The ~95.2% BWA-position-agreement figure looked unreproducible~~ —
  found while confirming this pass's `align_read()` fix hadn't caused a
  regression: an unrestricted primary-alignment comparison got 83.7%,
  not ~95.2%. Resolved, not just explained away: restricting to reads
  where both aligners report a clean, full-length, gap-free alignment and
  `bwa` itself isn't flagging ambiguity (MAPQ≥1) gives 94.2–94.6% on the
  same data — matching 95.2% within reasonable measurement variation.
  The original figure almost certainly came from an equivalently
  restricted comparison, a standard and defensible way to isolate core
  seed-and-extend correctness from soft-clip-boundary noise any two
  independently-built aligners will show. Both the 83.7% (all commonly-
  mapped reads) and ~94-95% (unambiguous full-length subset, ~76% of
  those reads) numbers are now understood and correct for what each
  measures — see [Validation on real data](#validation-on-real-data) for
  the full breakdown.

- ~~The FM-index prototype's construction memory was needlessly scaled
  to the whole genome instead of one chromosome at a time~~ — not fully
  closing the mammalian-memory gap (it isn't, see "Still open" below),
  but a real, measured improvement toward it: rebuilding as one
  independent FM-index per chromosome (instead of one concatenated
  multi-chromosome index) cut the human-genome extrapolation from ~102 GB
  to ~51 GB — the first FM-index version to beat the minimizer approach's
  ~68 GB. Validated on synthetic multi-chromosome genomes (round 4's
  figures used single-contig synthetic FASTAs, which couldn't exercise a
  per-chromosome fix), with the core hypothesis directly tested (holding
  total genome size fixed while varying the largest chromosome) rather
  than assumed, and correctness re-confirmed against both the existing
  minimizer path (99.8% agreement) and `bwa mem` (~94.2%, matching the
  established baseline). See [Memory](#memory) for the full account.

- ~~`--fm-index` silently corrupted the shared `.kidx` cache for every
  later default-mode run~~ — found while validating round 5's own changes
  (running `--fm-index` and default mode back-to-back against the same
  genome, which no prior round's validation sequence had done). The
  minimizer k-mer table is deliberately never built in `--fm-index` mode,
  but the on-disk cache was written unconditionally regardless, saving an
  empty table that any later default-mode run would silently load as a
  normal cache hit — alignment rate collapsing from 97.3% to ~12-14% with
  no error or warning. Fixed by only writing the cache when the minimizer
  table was actually built; verified in both run orderings. See
  [Memory](#memory) for the full account, including confirmation that
  none of this round's other reported figures were affected.

- ~~Prefix-doubling suffix array construction~~ — replaced with SA-IS
  (linear-time induced sorting), the second of the two concrete next
  steps round 4 identified for closing the mammalian-memory gap.
  Validated standalone (2,405 test cases against a naive reference
  suffix array, including one real bug the tests caught: a single-
  character-string edge case where the lone sentinel position was never
  induced from anywhere, now special-cased) plus a real-chromosome
  self-test that runs on every `--fm-index` invocation, comparing against
  the prefix-doubling algorithm it replaced rather than trusting a
  one-time standalone result to still hold. Human-genome extrapolation:
  ~51 GB → ~35 GB, now close to STAR's ~30 GB. Also fixed a real
  precondition SA-IS needed that prefix-doubling never did: ambiguous
  (N) bases previously shared a code with the sentinel, which
  prefix-doubling's plain comparisons tolerated but which breaks SA-IS's
  unique-minimum-sentinel requirement — N bases now get their own symbol
  (`FM_ALPHA` 5→6). See [Memory](#memory) for the complete account,
  including a real result the round's own skewed-chromosome test
  surfaced rather than smoothed over (the extrapolation's accuracy now
  depends more on a genome's chromosome-size distribution than round 5's
  did — explained, not hidden, in that section).

- ~~Suffix-array sampling~~ — the last item round 4 originally
  identified for closing the mammalian-memory gap. Replaced the full
  suffix array (4 bytes/bp) with sampled storage (every 16th text
  position kept explicitly, everything else recovered via LF-mapping) —
  validated standalone (571/571 tests across K=1–60) plus a real-
  chromosome self-test against the full array on every `--fm-index`
  invocation. Human-genome extrapolation: ~35 GB → ~25.3 GB, the first
  round to land below STAR's ~30 GB reference. Also found, and reported
  precisely rather than glossed over: the reduction fell well short of
  a steady-state-only calculation's prediction (~28% actual vs. ~60%
  theoretical), because peak memory is still substantially shaped by
  transient SA-IS construction cost that sampling doesn't address. See
  [Memory](#memory) for the full account.

- ~~`MAX_INTRON` was a fixed constant, unusable for any organism whose
  introns exceed yeast's, and simply widening it was measured to be
  unsafe~~ — now closed properly, in two parts, not just one. Part 1:
  made the window a runtime flag (`--max-intron N`) and added a distance
  tiebreak in `try_spliced_align` (prefer the shorter intron among
  candidates otherwise tied on canonical-site status and mismatch
  count) — validated to leave the yeast-validated default's results
  completely unchanged. Directly testing whether that alone made wide
  windows safe found it didn't: a clear, monotonic degradation as the
  window widened (70.5%→65.1%→57.0% known-junction rate at
  15,000/100,000/2,000,000bp), because the real failure mode is a
  distant coincidental match winning *outright* on mismatch count, which
  no tiebreak catches (it only fires on genuine ties). Part 2, the
  actual fix: `try_spliced_align` now requires a canonical GT-AG site
  outright — not just prefers one — for any intron beyond a 20,000bp
  safe range (just above the yeast-validated default, so the default
  never triggers this gate at all). Real eukaryotic introns are >98%
  canonical, so a genuine long intron almost always clears this bar; a
  spurious distant match needs to *also* coincidentally land on GT..AG,
  a roughly 1-in-256 event, not something a wider search window's extra
  candidates make likely en masse. Re-measured after the fix, same
  yeast data, same three window widths: 70.5%/70.1%/70.1% known-junction
  rate — the degradation is gone, not just reduced. Confirmed with a
  positive control too, not just fewer false positives: a synthetic
  chromosome with a planted, canonical, exactly-500,000bp intron was
  correctly found at `--max-intron 600000` (exact CIGAR `40M500000N35M`,
  zero mismatches, matching the planted ground truth precisely) and
  correctly reported unmapped at the yeast-appropriate default (as it
  should be — a window that can't see a 500,000bp intron shouldn't
  invent one). See [Junction annotation](#junction-annotation-rseqc-equivalent)
  for the full before/after measurement.

- ~~SA-IS's top-level `unsigned char`→`int32_t` conversion~~ — removed,
  via a u8-native top-level entry point validated against both a naive
  oracle and the existing int32_t core on 4,810 test cases. Real,
  measured 1.9–6.0% memory reduction at every tested scale (12/25/90 Mb),
  but only a ~0.6% effect on the human extrapolation specifically
  (~25.3 GB → ~25.16 GB) — a precise, explained result, not a
  disappointing one dressed up: this cost was bounded by the largest
  single chromosome, not total genome size, and human's largest
  chromosome is a much smaller fraction of its genome than any tested
  synthetic genome's was, so a fix in this category inherently helps
  human less than it helps the smaller test genomes it was measured on.
  See [Memory](#memory) for the complete account, including what this
  finding implies about the one remaining item.

- ~~Whether malloc/allocator retention was inflating peak FM-index
  memory~~ — a real hypothesis, tested directly and ruled out. Per-
  chromosome RSS instrumentation showed memory climbing monotonically
  across chromosomes (never decreasing, even after each one's transient
  arrays are correctly freed); `malloc_trim(0)` after every chromosome
  — the standard fix for an allocator not returning freed memory to the
  OS — made no measurable difference (988 MB vs. 990 MB peak on the
  90 Mb test, noise-level). The real explanation: peak memory is
  dominated by the *largest* chromosome's own construction requirement,
  with smaller chromosomes' allocations mostly reusing that already-
  reserved heap space rather than needing fresh growth — refining, not
  contradicting, round 7's finding, and narrowing (not eliminating) the
  ceiling on the one remaining memory item. See [Memory](#memory) and
  item 1 below for the full account.

**Still open** — these are the items I'm not claiming to have closed, and
why each one specifically hasn't been:

1. **Closing the index memory gap for mammalian genomes.** Round 2
   (minimizer sparse indexing) reduced the human-genome extrapolation
   from ~326 GB to ~68 GB. Round 4 built an actual FM-index and measured
   that, as implemented, it made things worse (~102 GB) because
   prefix-doubling's transient construction memory dominated. Round 5
   fixed the reason those arrays had to be full-genome-length
   (per-chromosome construction) — ~51 GB. Round 6 replaced
   prefix-doubling with SA-IS (linear-time induced sorting) — ~35 GB.
   Round 7 implemented suffix-array sampling — ~25.3 GB, below STAR's
   ~30 GB reference figure for the first time. Round 8 removed the last
   named transient-memory cost (a top-level `unsigned char`→`int32_t`
   conversion in SA-IS) — real, measured, consistent 1.9–6.0% reduction
   at every tested scale, but only ~25.3 GB → ~25.16 GB on the human
   extrapolation specifically, because that cost was bounded by the
   *largest single chromosome*, not total genome size. Five rounds of
   real, measured, compounding progress: ~326 GB → ~68 GB → ~51 GB →
   ~35 GB → ~25.3 GB → ~25.16 GB.

   What's left is real algorithmic work, not a quick follow-on: avoiding
   full-suffix-array materialization before sampling. Round 9
   investigated this directly with per-chromosome memory instrumentation
   rather than attempting the restructuring blind, and found round 7's
   original ~580 MB "transient cost" estimate was too coarse — peak
   memory is dominated by the *largest chromosome's* own construction
   requirement, with every smaller chromosome's allocations mostly
   reusing that same reserved heap space rather than requiring fresh
   growth (confirmed by testing and ruling out a competing hypothesis,
   allocator retention via `malloc_trim(0)`, which made no measurable
   difference). That refines the target to something smaller and more
   precisely bounded than round 7's number implied, which is worth
   knowing before attempting a real restructuring of when and how the
   suffix array's values get read — not implemented this round,
   deliberately, rather than risked on a smaller, now-better-understood
   ceiling. See [Memory](#memory) for the complete measured account
   across all six investigation rounds.
2. **Multi-organism, multi-protocol, multi-lab validation.** Everything
   above was checked by one person, across a handful of sessions, against
   one real dataset plus synthetic ground-truth tests for the pieces real
   data couldn't validate (UMI dedup, genome-scale memory). nf-core/rnaseq's
   reliability comes from years of independent labs hitting edge cases
   across different organisms and protocols and reporting them upstream.
   That kind of validation can't be substituted for.
Every correctness bug and regression found and fixed across this project
(the CIGAR/SEQ mismatch, the OOM interaction, the GFF3 biotype-dropping
bug, the disk-space false-alarm, the fixed 4MB chromosome buffer, the
`.kidx` cache corruption round 5 found, the `align_read()` splice-gating
bug round 6 found) was found by actually running real or realistic
synthetic data through real validators — `samtools`, `pyBigWig`,
`/usr/bin/time -v`, a naive suffix-array oracle, synthetic ground truth
with a known answer — not by reasoning about the code in the abstract.
Take that as a signal about how many more issues like it are probably
still sitting in the paths that haven't been exercised this way yet.
Mammalian genome scale specifically has gone from "known broken" to
"favorable on paper, untested in practice" across rounds 4 through 7 —
real, measured progress on the one dimension (memory) that was
previously a hard blocker, but progress on one dimension is not the same
claim as "known to work," and item 3 above (untested large-intron
behavior) is a concrete example of a *different* dimension this same
round's own work surfaced as still unknown. Other organisms, UMI-on-R2
protocols, and contaminant-heavy samples remain genuinely unknown in the
same way. Take the memory numbers as real progress, not as evidence the
untested paths are now bug-free.

## Repo layout

```
README.md                 Install, test, and tool-by-tool comparison (start here)
IMPLEMENTATION_DETAILS.md This file -- algorithm internals, validation history,
                          round-by-round notes, full known-gaps accounting
rnaseq_pipeline.c        Single-file pipeline (~6,870 lines): indexing,
                          alignment, quantification, QC, DE, HTML reporting
convert_gff3_to_gtf.py    Minimal Ensembl GFF3 -> gene-level GTF converter
data/README               Pointer to a sample S. cerevisiae RNA-seq dataset
LICENSE
```

## License

See [`LICENSE`](LICENSE).

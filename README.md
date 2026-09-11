A **single-file, C re-implementation** of the core analytical workflow used by [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq)
— FASTQ QC/trimming, UMI extraction/dedup, splice-aware alignment, EM-based gene quantification, a DESeq2-equivalent differential
expression stack, and self-contained HTML/SVG reporting. No aligner, no R, no Bioconductor, no Nextflow required for basic run — it compiles with `gcc` and runs
anywhere a C compiler exists.

**Validated on:** yeast (*S. cerevisiae*, ~12 Mb). **Not validated on** human/mouse-scale genomes or against real biological (non-technical) replicates.

Here's how to build and run the code in this repo:

## 1. Compile
**Assumes:** a C compiler (`gcc`) is available, and you're running the command from the directory containing `rnaseq_pipeline.c`

```bash
gcc -O2 -o rnaseq_pipeline rnaseq_pipeline.c -lm
```

Or with OpenMP for multi-core parallelism (worth using if your machine has more than 1 core - but the QC/trim and alignment loops are both parallelized):

```bash
gcc -O2 -fopenmp -o rnaseq_pipeline rnaseq_pipeline.c -lm
```

That's the whole build — no external libraries beyond the C standard library and libm. Everything below is optional and only adds specific
extra outputs on top of a working default run:

| Tool | Adds | Install | Assumes |
|---|---|---|---|
| `samtools` | Sorted, indexed BAM | `apt install samtools` / `brew install samtools` / `conda install -c bioconda samtools` | `apt` needs `sudo`; `brew` needs [Homebrew](https://brew.sh) already installed; `conda` needs conda/mamba already installed |
| `bedtools` + Python 3 `pyBigWig` | bigWig coverage track | `apt install bedtools` + `pip3 install pyBigWig --break-system-packages` (or `conda install -c bioconda bedtools pybigwig`) | Same as above, plus Python 3 + `pip3` already present |
| Python 3 (stdlib only) | GFF3→GTF conversion, if your annotation isn't already a GTF | Already on virtually any Linux/macOS/HPC system | — |

If a tool above isn't on `PATH`, the pipeline detects that itself, prints a `note:` explaining what it's skipping, and keeps going —
nothing fails because a tool is missing. Full multi-platform install commands (HPC module systems, macOS, troubleshooting `pip`'s
"externally-managed-environment" error) are in [IMPLEMENTATION_DETAILS.md](IMPLEMENTATION_DETAILS.md#installing-dependencies).


## 2. Annotation format — one thing to watch for

The pipeline expects **GTF**, not GFF3. The gene id and gene location information is in the GFF3 file, so it has to be converted before testing. If you are working from the same Ensembl-style GFF3, you'll need to do the same - a minimal gene-level conversion:

```bash
python3 -c "
import gzip, re
gene_re = re.compile(r'gene_id=([^;]+)')
with gzip.open('annotation.gff3.gz','rt') as f, open('annotation.gtf','w') as out:
    for line in f:
        if line.startswith('#'): continue
        p = line.rstrip('\n').split('\t')
        if len(p) < 9 or p[2] != 'gene': continue
        m = gene_re.search(p[8])
        if not m: continue
        out.write(f'{p[0]}\t{p[1]}\tgene\t{p[3]}\t{p[4]}\t.\t{p[6]}\t.\tgene_id \"{m.group(1)}\";\n')
"
```


Or


Use the convert_gff3_to_gtk.py code and run it as below. For the sample dataset, check the link in the the data folder of this repo and download RNA-Seq_Sample_Files. From inside the RNA-Seq_Sample_Files folder, run the following command (assuming that *.py file is in the folder one level up) to convert the *.gff3.gz file into a *.gtf file.

**Assumes:** Python 3 is available; *.gff3.gz is an Ensembl-style GFF3 (this converter is not a general-purpose GFF3 parser — other GFF3 flavors may use different attribute conventions). Skip this step entirely if you already have a GTF.

```
python3 ../convert_gff3_to_gtf.py Saccharomyces_cerevisiae.R64-1-1.59.gff3.gz yeast.gtf
```

## 3. Run

Before running the code, ensure that the input files are available. Create an outdir if one does not exists (mkdir outdir). 

Note: In the latest version of the C code, `outdir/` is created automatically if it doesn't exist.

**Remember to follow the step above to convert your annotation to GTF**, if it's GFF3 (Ensembl-style annotations, e.g. yeast, typically are).

## SERIAL RUN:

**Paired-end** (e.g., with `SRR9336476`):
```bash
./rnaseq_pipeline genome.fa annotation.gtf pe reads_1.fastq reads_2.fastq outdir/
```

**Single-end:**
```bash
./rnaseq_pipeline genome.fa annotation.gtf se reads.fastq outdir/
```

Reads can be plain `.fastq` — the pipeline doesn't read `.gz` directly, so `gunzip` first:

```bash
gunzip -k reads_1.fastq.gz reads_2.fastq.gz
```

Using the sample data set downloaded from the link in the data directory of this repo, and after creating yeast.gtf in step 2, you can test the code as follows:

```
time ./rnaseq_pipeline RNA-Seq_Sample_Files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa RNA-Seq_Sample_Files/yeast.gtf pe RNA-Seq_Sample_Files/1M_SRR9336468_1.fastq RNA-Seq_Sample_Files/1M_SRR9336468_2.fastq ./outdir/
``` 

## PARALLEL RUN
Set the OMP_NUM_THREADS variable to the number of threads that you would like to run the code with to reduce the overall run-time.

```
export OMP_NUM_THREADS=4
```

You can set OMP_NUM_THREADS to the maximum number of cores available on the CPU on which you will be running the code.

Make sure to compile the code with the -fopenmp flag.

```bash
gcc -O2 -fopenmp -o rnaseq_pipeline_omp rnaseq_pipeline.c -lm
```

Then, run the code as follows (similar to the command for the serial run):
```
time ./rnaseq_pipeline_omp RNA-Seq_Sample_Files/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa RNA-Seq_Sample_Files/yeast.gtf pe RNA-Seq_Sample_Files/1M_SRR9336468_1.fastq RNA-Seq_Sample_Files/1M_SRR9336468_2.fastq ./outdir/
```

**Multi-sample comparison** (PCA, clustering, and - with condition labels - differential expression), run *after* each sample has its own single-sample output directory:

```bash
# QC/clustering only
./rnaseq_pipeline compare out_rep1 out_rep2 out_rep3 combined_outdir/

# With a 2-condition design (needs >=2 replicates per condition) to also run the DESeq2-equivalent differential expression test
./rnaseq_pipeline compare \
  out_ctrl1:control out_ctrl2:control \
  out_trt1:treated  out_trt2:treated \
  combined_outdir/
```

## 4. What you will see

```
[1/7] Loading reference + k-mer index (k=16): genome.fa
[2/7] Building k-mer index (k=16)
      -> index cached for future runs: genome.fa.kidx
[3/7] Loading GTF annotation: annotation.gtf
[4/7] Loading reads (paired-end mode): ...
[5/7] Running QC + adapter/quality trimming
[6/7] Aligning reads ...
[7/7] EM quantification + writing output files to: outdir
```

Outputs land in `outdir/`: `alignments.sam`, `gene_counts.tsv`, `qc_report.txt`, `multiqc_summary.txt`.

## 5. About that `.kidx` file

The first run against a given genome builds the k-mer index and writes `genome.fa.kidx` next to it - that step is the slow part. Every subsequent run against the **same, unmodified** `genome.fa` will load from that cache instead of rebuilding, which is dramatically faster. If you ever replace `genome.fa` with different content but keep the same filename, the pipeline checks file size + modification time and rebuilds automatically - you don't need to delete the cache manually.

## 6. Confirm it worked

**Assumes:** step 3 finished without printing `ERROR:` and returned to your shell prompt (if it's still running, wait for it — a full yeast dataset can take a few minutes; nothing below will exist until it's done).

```bash
grep "Fragments aligned" outdir/multiqc_summary.txt
```

A healthy run against real RNA-seq data reads roughly `95-98%` here — if it's near `0%`, something upstream is wrong (wrong genome/annotation pairing, mismatched read files, or similar) before you look any further. Then check the core outputs exist and are non-empty:

```bash
ls -la outdir/gene_counts.tsv outdir/report.html outdir/multiqc_summary.txt
```

Open `outdir/report.html` in a browser for the full QC report (mapping breakdown, quality/GC/length distributions, top genes).

## 7. Multi-sample comparison** (optional)
This step combines results from samples you have *already run individually* through steps 1–3 above; it does not take raw FASTQ files itself.

**Assumes:** you have run step 3 **separately for each sample**, into its own output directory (so if you have 4 samples, you've already run the pipeline 4 times, producing 4 directories, each with its own `gene_counts.tsv` inside it) — `ctrl1`, `ctrl2`, `trt1`, `trt2` below are placeholders for whatever you named *those* directories, not fixed keywords or raw input files. The `:control` / `:treated` condition labels are also names you choose (any two consistently-used labels
work); they are required — with **at least 2 sample directories per label** — only if you want the differential-expression test, not for PCA/clustering alone (omit `:label` from every sample to skip it). `combined_outdir/` is a *new*, separate directory for the combined report — not one of the four input directories, and created automatically if it doesn't exist.

```bash
./rnaseq_pipeline compare ctrl1:control ctrl2:control trt1:treated trt2:treated combined_outdir/
```

| Piece | Meaning |
|---|---|
| `compare` | Combine existing per-sample results — not "align new reads" |
| `ctrl1` | Path to a directory from step 2 — reads `ctrl1/gene_counts.tsv` from inside it |
| `:control` | Your chosen label grouping this sample into the "control" condition |
| `ctrl2:control` | A second replicate, same condition — DE testing needs ≥2 samples per condition |
| `trt1:treated`, `trt2:treated` | The second condition, same rule |
| `combined_outdir/` | New directory for the combined report |

Produces `combined_outdir/multi_sample_report.html`: PCA, hierarchical clustering, and — with the `:label` condition tags shown above — a full differential expression table with a volcano plot.

## 8. Comparison with nf-core/rnaseq and other tools

nf-core/rnaseq orchestrates ~20 external tools (FastQC, Trim Galore!, STAR, Salmon, RSEM, HISAT2, SAMtools, UMI-tools, picard, RSeQC, Qualimap, dupRadar, Preseq, DESeq2, MultiQC...) behind a Nextflow interface. This repo reimplements the analytically central stages natively in C:

| nf-core/rnaseq stage | Reference tool(s) | This repo | Status |
|---|---|---|---|
| Merge re-sequenced FASTQs | `cat` | Comma-separated multi-lane FASTQ paths, concatenated before trimming | Mirrored |
| Strandedness inference | `fq`, `Salmon` | RSeQC `infer_experiment.py`-style classifier | Mirrored |
| Read QC | `FastQC` | Per-read length/GC/quality-by-cycle → HTML report | Mirrored |
| UMI extraction | `UMI-tools` | In-line UMI-on-R1 layouts (e.g. QIAseq/NEBNext) | Mirrored (common case) |
| Adapter/quality trimming | `Trim Galore!` | 3′ quality + adapter overlap trim | Mirrored (simplified) |
| Contaminant genome removal | `BBSplit` | `--contaminant-ref`, repeatable, priority-order screening | Mirrored |
| rRNA content | `SortMeRNA` | Annotation-based rRNA/biotype quantification (not pre-alignment removal) | Mirrored (QC only) |
| Alignment | `STAR` / `HISAT2` | k-mer-seeded banded Smith-Waterman, splice-aware, multi-mapper-aware; optional FM-index path (`--fm-index`) | Mirrored (yeast-scale validated only) |
| Sort/index alignments | `SAMtools` | Shells out to `samtools`, same as nf-core | Mirrored |
| UMI dedup | `UMI-tools` | Directional-adjacency clustering, same default method | Mirrored |
| Duplicate marking | `picard MarkDuplicates` | Position-based marking, feeds dupRadar-equivalent QC | Mirrored |
| Quantification | `Salmon` / `RSEM` | EM algorithm (unique reads anchor, multi-mappers redistributed) | Mirrored (core algorithm) |
| bigWig coverage | `BEDTools` | `bedtools genomecov` + `pyBigWig` | Mirrored |
| Per-sample QC | `RSeQC`, `Qualimap`, `dupRadar`, `Preseq` | Gene-body coverage, saturation curve, duplication-vs-expression, junction annotation | Mirrored |
| Differential expression | `DESeq2` | Median-of-ratios, Cox-Reid dispersion, GLM trend, empirical-Bayes shrinkage, VST, Cook's distance, BH-FDR, volcano plot | Mirrored |
| Pseudoalignment route | `Salmon` / `Kallisto` | Not applicable — this repo is alignment-based | Out of scope |
| Multi-sample QC/report | `MultiQC`, `R` | `compare` mode: PCA, clustering, DE report | Mirrored |
| Compressed input (`.gz`) | built-in | Transparent gunzip on any `.gz` input | Mirrored |

**Correctness, where independently checked:** ~94% position agreement with `bwa mem` on unambiguous reads; UMI dedup 96.7% recall against a
ground-truth set with deliberate sequencing-error variants; DE stack validated against DESeq2-style synthetic ground truth.

**Memory** (the one dimension checked at more than yeast scale): extrapolated to ~25 GB for a human genome via the optional `--fm-index` path — below STAR's ~30 GB reference figure — but this is an extrapolation from real measurements topping out at 90 Mb, not a validated human-scale run.


**Acknowledgement**: Claude/Sonnet 5.0 was leveraged in the process of developing and testing this code.


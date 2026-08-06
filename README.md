Here's how to build and run the code in this repo:

## 1. Compile

```bash
gcc -O2 -o rnaseq_pipeline rnaseq_pipeline.c -lm
```

Or with OpenMP for multi-core parallelism (worth using if your machine has more than 1 core — but the QC/trim and alignment loops are both parallelized):

```bash
gcc -O2 -fopenmp -o rnaseq_pipeline rnaseq_pipeline.c -lm
```

## 2. Annotation format — one thing to watch for

The pipeline expects **GTF**, not GFF3. The gene id and gene location information is in the GFF3 file, so it has to be converted before testing. If you are working from the same Ensembl-style GFF3, you'll need to do the same — a minimal gene-level conversion:

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

```
python3 ../convert_gff3_to_gtf.py Saccharomyces_cerevisiae.R64-1-1.59.gff3.gz yeast.gtf
```

## 3. Run

Before running the code, ensure that the input files are available.

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

Note: If the files in the downloaded dataset folder are in *.gz format, then run the following command to convert the files into *.fastq format before rnning the steps above:

```
gunzip *
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

The first run against a given genome builds the k-mer index and writes `genome.fa.kidx` next to it — that step is the slow part. Every subsequent run against the **same, unmodified** `genome.fa` will load from that cache instead of rebuilding, which is dramatically faster. If you ever replace `genome.fa` with different content but keep the same filename, the pipeline checks file size + modification time and rebuilds automatically — you don't need to delete the cache manually.

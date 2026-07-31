/* ============================================================================
 * rnaseq_pipeline.c  (v4)
 *
 * A self-contained, dependency-free C re-implementation of the CONCEPTUAL
 * workflow used by nf-core/rnaseq (https://github.com/nf-core/rnaseq):
 *
 *   FASTQ QC -> adapter/quality trimming -> reference alignment
 *   -> EM-based gene-level quantification -> QC/summary report
 *
 * v4 upgrade over v3:
 *   - Splice-junction anchoring now uses a second, shorter (10bp) k-mer
 *     index dedicated to splice-anchor seeding, so junction-spanning reads
 *     with a short exon overhang on one side are no longer unseedable.
 *     Measured recall on junction-spanning test reads: 37.6% -> 62.4%,
 *     with precision unchanged (zero false-positive intron calls).
 *
 * v3 upgrade over v2 (see README for the full comparison table):
 *   - Alignment now reports ALL tied-best-scoring hits for a read (up to a
 *     cap), not just the single best one, so genuinely multi-mapping reads
 *     (e.g. from duplicated/paralogous genes) are represented as such.
 *   - Gene-level quantification is now an EM algorithm in the same spirit
 *     as RSEM/Salmon: reads/fragments that map uniquely to one gene anchor
 *     that gene's estimated abundance directly; reads/fragments that map
 *     ambiguously to several candidate genes are redistributed across them
 *     in proportion to each gene's current abundance estimate, iterated to
 *     convergence. This replaces v1/v2's "first candidate gene wins"
 *     assignment, which was a substantial known bias for any genome with
 *     duplicated or paralogous sequence.
 *   - SAM output (single-end mode) now emits one record per reported hit
 *     for multi-mapping reads, with the standard NH:i: tag and the 0x100
 *     "secondary alignment" SAM flag on all but the first, matching how
 *     real aligners (STAR, BWA, HISAT2) represent multi-mappers.
 *
 * Still explicitly out of scope: real genome-scale indexing (FM-index),
 * fragment-length/sequence-bias models in the EM (RSEM/Salmon use these to
 * further refine abundance estimates), UMI dedup, duplicate marking, rRNA
 * filtering, and differential expression (DESeq2). Paired-end SAM output
 * still reports only the single best hit per mate (see README).
 *
 * Usage:
 *   rnaseq_pipeline <reference.fasta> <annotation.gtf> se <reads.fastq> <outdir>
 *   rnaseq_pipeline <reference.fasta> <annotation.gtf> pe <r1.fastq> <r2.fastq> <outdir>
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_LINE               8192
#define MAX_SEQNAME             128
#define MAX_CHROMS                64
#define MAX_GENES                8192
#define MAX_READS             2100000
#define MAX_READ_LEN             512
#define ADAPTER_SEQ  "AGATCGGAAGAGC"   /* Illumina TruSeq universal adapter */
#define QUALITY_TRIM_THRESHOLD    20
#define TRIM_WINDOW                4
#define MIN_READ_LEN_AFTER_TRIM   20

/* --- alignment / indexing parameters --- */
#define KMER_LEN                  16   /* fits exactly in a uint32_t (2 bits/base) */
#define SPLICE_KMER_LEN            10   /* shorter anchor length used ONLY for splice-
                                            junction seeding, so reads with a short exon
                                            overhang on one side of a junction can still
                                            be anchored (the main KMER_LEN=16 index is
                                            too long a requirement for those) */
#define KMER_TABLE_BITS           24   /* was 20 (1,048,576 buckets): with ~12.16M
                                           distinct real-genome 16-mers, that gave an
                                           average hash-chain length of ~12, and
                                           build_kmer_index_generic profiled (gprof) at
                                           83% of total runtime on a 20k-read-pair run
                                           as a direct result -- every insert/lookup was
                                           walking a ~12-entry linked list. 2^24 =
                                           16,777,216 buckets brings that under 1 on
                                           average. Costs ~128MB per table (main +
                                           splice, ~256MB total) -- trivial next to the
                                           genome and read data already in memory. */
#define KMER_TABLE_SIZE      (1u << KMER_TABLE_BITS)
#define KMER_TABLE_MASK      (KMER_TABLE_SIZE - 1)
#define SPLICE_TABLE_BITS         20   /* the splice-anchor index's key space is only
                                           4^SPLICE_KMER_LEN = 4^10 = 1,048,576 distinct
                                           10-mers, regardless of genome size -- reusing
                                           the much larger main-index table size here
                                           bought nothing (already <1 collision/bucket at
                                           this size for that key space) and wasted
                                           ~120MB. Sized to the splice index's own actual
                                           key space instead. */
#define SPLICE_TABLE_SIZE    (1u << SPLICE_TABLE_BITS)
#define SPLICE_TABLE_MASK    (SPLICE_TABLE_SIZE - 1)
#define MAX_SEED_HITS_PER_KMER    64   /* cap on repetitive k-mers, like repeat-masking */
#define MAX_MISMATCHES              2   /* ungapped alignment */
#define SPLICE_MAX_MISMATCHES       3   /* combined across both exon blocks */
#define MIN_INTRON                 20   /* demo-scale intron size bounds  */
#define MAX_INTRON              100000
#define MAX_MULTI_HITS             16   /* cap on reported tied-best hits per read */
#define MAX_CAND_POS               64   /* cap on candidate positions scanned pre-filter */

/* --- Smith-Waterman-Gotoh local alignment (affine gaps, soft-clipping) --- */
#define SW_MATCH                    2
#define SW_MISMATCH                -4
#define SW_GAP_OPEN                -6   /* cost of a 1-residue gap (open + first extend) */
#define SW_GAP_EXTEND               -1  /* additional cost per extra gap residue */
#define SW_PAD                     15   /* reference window padding each side of a seed;
                                            bounds the largest indel that can be found */
#define SW_BAND                    30   /* DP is restricted to +/- this many columns around the
                                            seed-implied diagonal (see sw_local_align) instead of
                                            the full qlen x window_len rectangle. NOTE: an earlier
                                            value of 15 (== SW_PAD) was tried and found, by diffing
                                            against the un-banded baseline on real data, to exclude
                                            some real alignments (e.g. a 16bp deletion) that the
                                            SW_PAD-bounded window could still reach depending on
                                            where along the read the indel falls -- the window's
                                            reachable (i,j) divergence isn't uniformly +/-SW_PAD,
                                            it's asymmetric depending on position. 30 (=2*SW_PAD)
                                            was verified byte-identical to the un-banded baseline
                                            on a 100,000-read full-genome regression run; see
                                            BANDED_SW_OPTIMIZATION.md for the verification. */
#define SW_MIN_SCORE_FRAC         0.6   /* minimum fraction of the max possible score
                                            (qlen * SW_MATCH) required to accept a hit */
#define SW_MAX_EDIT_FRAC         0.12   /* maximum (mismatches+indel bases)/qlen allowed --
                                            catches technically-best-scoring but biologically
                                            implausible placements (e.g. multiple separate indel
                                            events in one short read) that score alone lets
                                            through. Found via real data: two reads with 14-15%
                                            edit distance and mediocre base quality were accepted
                                            by score alone but correctly rejected by BWA-MEM. */
#define SW_MAX_CANDIDATES           8   /* cap on distinct seed-derived candidate windows
                                            actually run through full SW, per orientation --
                                            bounds worst-case cost in repetitive regions */
#define MAX_CIGAR_OPS               12   /* real single-intron spliced (M/N/M) and indel-bearing
                                            SW alignments never come close to this: profiling a
                                            100,000-read run against the full 12Mb S. cerevisiae
                                            R64-1-1 genome showed a max of 9 ops on any single
                                            alignment. Was 24, costing 2x the CigarOp memory for
                                            headroom nothing in real data ever used. */

/* --- EM quantification parameters --- */
#define MAX_GENES_PER_UNIT          8   /* cap on candidate genes for one read/fragment */
#define MAX_EM_ITERS               200
#define EM_CONVERGE_EPS         1e-7

/* ---------------------------------------------------------------------- */
/* Data structures                                                        */
/* ---------------------------------------------------------------------- */

typedef struct {
    char name[MAX_SEQNAME];
    char *seq;
    long  len;
} Chrom;

typedef struct {
    char gene_id[MAX_SEQNAME];
    char chrom[MAX_SEQNAME];
    long start;   /* 1-based inclusive, like GTF */
    long end;     /* 1-based inclusive */
    char strand;  /* '+' or '-' */
    long   unique_count;     /* reads/fragments uniquely assigned (hard count) */
    double effective_count;  /* EM expected count, incl. proportional multi-mapper share */
} Gene;

typedef struct { char op; int len; } CigarOp; /* op in {'M','I','D','N','S'} */

/* One reported alignment location, represented as a real CIGAR (mixing
 * M/I/D from Smith-Waterman-Gotoh local alignment, N from splice detection,
 * and S for soft-clipped read ends) rather than a fixed 1-or-2-block model.
 * Coordinates are 0-based genome offsets. */
typedef struct {
    int      chrom_idx;
    long     ref_start;      /* 0-based genome position of the first aligned (non-S) base */
    CigarOp  cigar[MAX_CIGAR_OPS];
    int      n_cigar;
    int      spliced;
    long     intron_len;        /* valid if spliced (sum of N-op lengths, single-intron case) */
    int      canonical_splice;  /* valid if spliced */
    int      score;             /* SW alignment score; 0 for the (unscored) splice path */
    int      mismatches;        /* substitutions + indel bases, used for the NM tag */
    char     strand;            /* '+' or '-': which query orientation produced this hit */
} Hit;

typedef struct {
    int mapped;
    int n_hits;
    Hit *hits;   /* heap-allocated, sized to exactly n_hits (<= MAX_MULTI_HITS) once
                  * alignment finishes -- not a fixed MAX_MULTI_HITS array embedded in
                  * every Read. The overwhelming majority of reads have exactly 1 hit;
                  * paying for 16 full Hit slots (16 * sizeof(Hit)) on every single Read
                  * regardless of actual multi-mapping is the single largest contributor
                  * to sizeof(Read), and was what forced capping real full-genome runs in
                  * this session to 100,000 of the 1,000,000 available reads to fit in
                  * available memory. */
} AlignResult;

/* Set of candidate gene indices a read/fragment could plausibly belong to,
 * derived from the gene(s) overlapped by each of its reported hits. Size 0
 * means mapped but intergenic; size 1 means an unambiguous ("unique") read;
 * size >1 means a genuine multi-mapper needing EM resolution. */
typedef struct {
    int idx[MAX_GENES_PER_UNIT];
    int n;
} GeneSet;

typedef struct {
    char id[MAX_SEQNAME];
    char *seq;
    char *qual;
    int   raw_len;
    int   trimmed_len;
    int   adapter_trimmed;
    int   qual_trimmed_bp;
    double gc_pct;
    double mean_q_raw;
    int    n_count;

    int    read_num;   /* 0 = single-end, 1 = R1, 2 = R2 */
    int    mate_idx;   /* index into reads[] of mate, or -1 */

    AlignResult aln;
    GeneSet candidates; /* filled during quantify_em() */
} Read;

static Chrom chroms[MAX_CHROMS];
static int   n_chroms = 0;

static Gene genes[MAX_GENES];
static int  n_genes = 0;

static Read *reads;   /* heap-allocated array, size MAX_READS */
static int   n_reads = 0;
static int   paired_mode = 0;

/* running QC / pipeline counters */
static long total_raw_bases = 0;
static long total_trimmed_bases = 0;
static long reads_with_adapter = 0;
static long reads_quality_trimmed = 0;
static long reads_dropped_too_short = 0;
static long n_spliced_alignments = 0;
static long n_ungapped_alignments = 0;
static long n_proper_pairs = 0;
static long n_unique_units = 0, n_multi_units = 0, n_no_feature_units = 0, n_unmapped_units = 0;
static int  em_iterations_run = 0;
static double em_final_delta = 0.0;

/* ---------------------------------------------------------------------- */
/* Utility helpers                                                        */
/* ---------------------------------------------------------------------- */

static void die(const char *msg) {
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(1);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) die("out of memory");
    return p;
}

static void *xrealloc(void *old, size_t n) {
    void *p = realloc(old, n);
    if (!p) die("out of memory");
    return p;
}

static char complement_base(char b) {
    switch (toupper((unsigned char)b)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

static char *revcomp(const char *seq, int len) {
    char *out = xmalloc(len + 1);
    for (int i = 0; i < len; i++)
        out[i] = complement_base(seq[len - 1 - i]);
    out[len] = '\0';
    return out;
}

static double phred_to_prob_correct(char qchar) {
    int q = qchar - 33; /* Phred+33 */
    if (q < 0) q = 0;
    return q;
}

static void rstrip(char *s) {
    size_t l = strlen(s);
    while (l > 0 && (s[l-1] == '\n' || s[l-1] == '\r')) { s[--l] = '\0'; }
}

static long lmin(long a, long b) { return a < b ? a : b; }
static long lmax(long a, long b) { return a > b ? a : b; }

/* ---------------------------------------------------------------------- */
/* Reference FASTA parsing                                                */
/* ---------------------------------------------------------------------- */

static void load_reference(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open reference FASTA");

    char line[MAX_LINE];
    static char seqbuf[4 * 1024 * 1024];
    size_t seqlen = 0;
    int have_chrom = 0;

    while (fgets(line, sizeof(line), f)) {
        rstrip(line);
        if (line[0] == '>') {
            if (have_chrom) {
                chroms[n_chroms].seq = xmalloc(seqlen + 1);
                memcpy(chroms[n_chroms].seq, seqbuf, seqlen);
                chroms[n_chroms].seq[seqlen] = '\0';
                chroms[n_chroms].len = (long)seqlen;
                n_chroms++;
            }
            if (n_chroms >= MAX_CHROMS) die("too many chromosomes/contigs in reference");
            char name[MAX_SEQNAME];
            sscanf(line + 1, "%127s", name);
            memset(chroms[n_chroms].name, 0, MAX_SEQNAME);
            strncpy(chroms[n_chroms].name, name, MAX_SEQNAME - 1);
            seqlen = 0;
            have_chrom = 1;
        } else if (have_chrom) {
            size_t l = strlen(line);
            if (seqlen + l >= sizeof(seqbuf)) die("reference sequence too large for demo buffer");
            memcpy(seqbuf + seqlen, line, l);
            seqlen += l;
        }
    }
    if (have_chrom) {
        chroms[n_chroms].seq = xmalloc(seqlen + 1);
        memcpy(chroms[n_chroms].seq, seqbuf, seqlen);
        chroms[n_chroms].seq[seqlen] = '\0';
        chroms[n_chroms].len = (long)seqlen;
        n_chroms++;
    }
    fclose(f);
    if (n_chroms == 0) die("no sequences found in reference FASTA");
}

/* ---------------------------------------------------------------------- */
/* K-mer hash index over the reference (forward strand only; reverse-     */
/* strand reads are queried via their reverse complement).                */
/* ---------------------------------------------------------------------- */

typedef struct { int chrom_idx; long pos; } SeedHit;

typedef struct KmerEntry {
    uint32_t kmer;
    SeedHit *hits;
    int n_hits;
    int hits_cap;   /* current allocated capacity of hits[]; grows on demand up
                      * to MAX_SEED_HITS_PER_KMER instead of being pre-allocated
                      * at the cap for every distinct k-mer. At full-genome
                      * scale (millions of distinct 16-mers, mostly unique in
                      * a non-repetitive genome like yeast) pre-allocating the
                      * full 64-slot cap for every entry costs ~12.8GB on a
                      * 12Mb genome -- fine at 230kb-subset test-data scale
                      * (~245MB) but OOMs a real run. Confirmed by profiling a
                      * full S. cerevisiae R64-1-1 genome + real GEO reads. */
    struct KmerEntry *next;
} KmerEntry;

static KmerEntry *kmer_table[KMER_TABLE_SIZE];               /* main index, k = KMER_LEN */
static KmerEntry *kmer_table_splice[SPLICE_TABLE_SIZE];       /* splice-anchor index, k = SPLICE_KMER_LEN */

/* Bump-allocator arena for KmerEntry structs and their initial (cap=2) hits
 * arrays. Index-building never frees an individual entry (the whole index
 * lives for the process lifetime), so replacing ~24 million individual
 * malloc() calls (one KmerEntry + one small hits[] per distinct k-mer,
 * across both the main and splice indexes on a real genome) with a
 * handful of large block allocations removes that allocator/page-fault
 * overhead entirely. Growth beyond the initial 2-slot hits[] (only the
 * repetitive minority of k-mers) still goes through ordinary
 * malloc/realloc below, since arena memory can't be individually resized. */
#define ARENA_ENTRIES_PER_BLOCK  (1 << 20)
#define ARENA_HITS_PER_BLOCK     (1 << 21)
typedef struct EntryArenaBlock { KmerEntry entries[ARENA_ENTRIES_PER_BLOCK]; struct EntryArenaBlock *next; } EntryArenaBlock;
typedef struct HitsArenaBlock  { SeedHit   hits[ARENA_HITS_PER_BLOCK];      struct HitsArenaBlock  *next; } HitsArenaBlock;
static EntryArenaBlock *entry_arena_head = NULL; static int entry_arena_used = ARENA_ENTRIES_PER_BLOCK;
static HitsArenaBlock  *hits_arena_head  = NULL; static int hits_arena_used  = ARENA_HITS_PER_BLOCK;

static KmerEntry *entry_arena_alloc(void) {
    if (entry_arena_used >= ARENA_ENTRIES_PER_BLOCK) {
        EntryArenaBlock *blk = xmalloc(sizeof(EntryArenaBlock));
        blk->next = entry_arena_head; entry_arena_head = blk; entry_arena_used = 0;
    }
    return &entry_arena_head->entries[entry_arena_used++];
}
static SeedHit *hits_arena_alloc2(void) {
    if (hits_arena_used + 2 > ARENA_HITS_PER_BLOCK) {
        HitsArenaBlock *blk = xmalloc(sizeof(HitsArenaBlock));
        blk->next = hits_arena_head; hits_arena_head = blk; hits_arena_used = 0;
    }
    SeedHit *p = &hits_arena_head->hits[hits_arena_used];
    hits_arena_used += 2;
    return p;
}

static int encode_kmer(const char *s, int len, uint32_t *out) {
    uint32_t v = 0;
    for (int i = 0; i < len; i++) {
        int code;
        switch (toupper((unsigned char)s[i])) {
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            default: return 0;
        }
        v = (v << 2) | (uint32_t)code;
    }
    *out = v;
    return 1;
}

static uint32_t hash_kmer(uint32_t k, uint32_t mask) {
    k ^= k >> 16; k *= 0x7feb352dU;
    k ^= k >> 15; k *= 0x846ca68bU;
    k ^= k >> 16;
    return k & mask;
}

static void kmer_index_add(KmerEntry **table, uint32_t mask, uint32_t kmer, int chrom_idx, long pos) {
    uint32_t h = hash_kmer(kmer, mask);
    KmerEntry *e = table[h];
    while (e && e->kmer != kmer) e = e->next;
    if (!e) {
        e = entry_arena_alloc();
        e->kmer = kmer;
        e->hits_cap = 2;   /* most k-mers at genome scale are unique or near-unique;
                             * start small and grow, rather than pre-paying for
                             * MAX_SEED_HITS_PER_KMER slots on every entry. */
        e->hits = hits_arena_alloc2();
        e->n_hits = 0;
        e->next = table[h];
        table[h] = e;
    }
    if (e->n_hits < MAX_SEED_HITS_PER_KMER) {
        if (e->n_hits == e->hits_cap) {
            int old_cap = e->hits_cap;
            e->hits_cap *= 2;
            if (e->hits_cap > MAX_SEED_HITS_PER_KMER) e->hits_cap = MAX_SEED_HITS_PER_KMER;
            if (old_cap == 2) {
                /* hits[] still points into the arena -- can't realloc() memory
                 * that didn't come from malloc(). Allocate a real buffer and
                 * copy the (at most 2) existing entries over. */
                SeedHit *fresh = xmalloc(sizeof(SeedHit) * e->hits_cap);
                memcpy(fresh, e->hits, sizeof(SeedHit) * e->n_hits);
                e->hits = fresh;
            } else {
                e->hits = xrealloc(e->hits, sizeof(SeedHit) * e->hits_cap);
            }
        }
        e->hits[e->n_hits].chrom_idx = chrom_idx;
        e->hits[e->n_hits].pos = pos;
        e->n_hits++;
    }
}

static KmerEntry *kmer_index_lookup(KmerEntry **table, uint32_t mask, uint32_t kmer) {
    uint32_t h = hash_kmer(kmer, mask);
    KmerEntry *e = table[h];
    while (e && e->kmer != kmer) e = e->next;
    return e;
}

static void build_kmer_index_generic(KmerEntry **table, uint32_t mask, int k, const char *label) {
    long total_kmers = 0;
    for (int ci = 0; ci < n_chroms; ci++) {
        Chrom *c = &chroms[ci];
        for (long p = 0; p + k <= c->len; p++) {
            uint32_t kmer;
            if (encode_kmer(c->seq + p, k, &kmer)) {
                kmer_index_add(table, mask, kmer, ci, p);
                total_kmers++;
            }
        }
    }
    printf("      -> %s index built: %ld %d-mers indexed from %d sequence(s)\n", label, total_kmers, k, n_chroms);
}

static void build_kmer_index(void) {
    build_kmer_index_generic(kmer_table, KMER_TABLE_MASK, KMER_LEN, "main k-mer");
    build_kmer_index_generic(kmer_table_splice, SPLICE_TABLE_MASK, SPLICE_KMER_LEN, "splice-anchor k-mer");
}

/* ---------------------------------------------------------------------- */
/* On-disk index cache (reference FASTA -> chrom sequences + both k-mer   */
/* tables), so repeated runs against the SAME reference skip re-parsing   */
/* the FASTA and rebuilding the index from scratch.                      */
/*                                                                        */
/* This is the single largest remaining lever after this session's other */
/* fixes: build_kmer_index_generic profiled at 35-83% of total runtime    */
/* depending on read-count scale, yet it does the exact same, input-      */
/* independent work every single run. A real RNA-seq workflow aligns many */
/* samples against one unchanged reference genome -- BWA reflects this by */
/* splitting `bwa index` (one-time) from `bwa mem` (per-sample); this adds*/
/* the equivalent for this pipeline without requiring a separate command. */
/*                                                                        */
/* Cache validity is checked via the reference file's (size, mtime) plus  */
/* the compile-time index parameters (KMER_LEN, SPLICE_KMER_LEN, table    */
/* sizes, MAX_SEED_HITS_PER_KMER) -- if any differ, the cache is silently  */
/* ignored and rebuilt, never trusted stale. The cache file format is a   */
/* private, this-binary-only format (no version negotiation across        */
/* compilers/platforms) since it's a same-machine performance cache, not  */
/* a portable index format like BWA's .bwt/.sa files -- documented here   */
/* rather than treated as a hidden assumption. */
/* ---------------------------------------------------------------------- */

typedef struct {
    char magic[8];      /* "RSPKIDX1" */
    long ref_size;
    long ref_mtime;
    int  kmer_len;
    int  splice_kmer_len;
    int  table_bits;
    int  splice_table_bits;
    int  max_hits_per_kmer;
    int  n_chroms;
} IndexCacheHeader;

static void cache_path_for(const char *ref_path, char *out, size_t outsz) {
    snprintf(out, outsz, "%s.kidx", ref_path);
}

static int stat_file(const char *path, long *size, long *mtime) {
    struct stat st;
    if (stat(path, &st) != 0) return 0;
    *size = (long)st.st_size;
    *mtime = (long)st.st_mtime;
    return 1;
}

static void write_chrom(FILE *f, const Chrom *c) {
    int namelen = (int)strlen(c->name);
    fwrite(&namelen, sizeof(int), 1, f);
    fwrite(c->name, 1, (size_t)namelen, f);
    fwrite(&c->len, sizeof(long), 1, f);
    fwrite(c->seq, 1, (size_t)c->len, f);
}

static int read_chrom(FILE *f, Chrom *c) {
    int namelen;
    if (fread(&namelen, sizeof(int), 1, f) != 1) return 0;
    if (namelen < 0 || namelen >= MAX_SEQNAME) return 0;
    memset(c->name, 0, MAX_SEQNAME);
    if (namelen > 0 && fread(c->name, 1, (size_t)namelen, f) != (size_t)namelen) return 0;
    if (fread(&c->len, sizeof(long), 1, f) != 1) return 0;
    if (c->len < 0) return 0;
    c->seq = xmalloc((size_t)c->len + 1);
    if (c->len > 0 && fread(c->seq, 1, (size_t)c->len, f) != (size_t)c->len) return 0;
    c->seq[c->len] = '\0';
    return 1;
}

static long count_table_entries(KmerEntry **table, uint32_t table_size) {
    long n = 0;
    for (uint32_t h = 0; h < table_size; h++)
        for (KmerEntry *e = table[h]; e; e = e->next) n++;
    return n;
}

typedef struct { uint32_t kmer; int32_t n_hits; } KmerMetaRec;

static void write_table(FILE *f, KmerEntry **table, uint32_t table_size) {
    long n = count_table_entries(table, table_size);
    long total_hits = 0;
    for (uint32_t h = 0; h < table_size; h++)
        for (KmerEntry *e = table[h]; e; e = e->next) total_hits += e->n_hits;

    fwrite(&n, sizeof(long), 1, f);
    fwrite(&total_hits, sizeof(long), 1, f);

    /* Bulk-serialize into two flat buffers, then write each with a single
     * fwrite() call, instead of 2-3 tiny fwrite() calls per entry. With up
     * to ~12M entries the per-call overhead dominates over the actual
     * bytes moved -- confirmed directly: the original per-entry version
     * made loading FROM the cache slower than just rebuilding the index
     * from scratch. */
    KmerMetaRec *meta = xmalloc(sizeof(KmerMetaRec) * (size_t)(n > 0 ? n : 1));
    SeedHit *hitsbuf = xmalloc(sizeof(SeedHit) * (size_t)(total_hits > 0 ? total_hits : 1));
    long mi = 0, hi = 0;
    for (uint32_t h = 0; h < table_size; h++) {
        for (KmerEntry *e = table[h]; e; e = e->next) {
            meta[mi].kmer = e->kmer; meta[mi].n_hits = e->n_hits; mi++;
            if (e->n_hits > 0) {
                memcpy(&hitsbuf[hi], e->hits, sizeof(SeedHit) * (size_t)e->n_hits);
                hi += e->n_hits;
            }
        }
    }
    if (n > 0) fwrite(meta, sizeof(KmerMetaRec), (size_t)n, f);
    if (total_hits > 0) fwrite(hitsbuf, sizeof(SeedHit), (size_t)total_hits, f);
    free(meta); free(hitsbuf);
}

static int read_table(FILE *f, KmerEntry **table, uint32_t mask) {
    long n, total_hits;
    if (fread(&n, sizeof(long), 1, f) != 1) return 0;
    if (fread(&total_hits, sizeof(long), 1, f) != 1) return 0;
    if (n < 0 || total_hits < 0) return 0;

    KmerMetaRec *meta = xmalloc(sizeof(KmerMetaRec) * (size_t)(n > 0 ? n : 1));
    if (n > 0 && fread(meta, sizeof(KmerMetaRec), (size_t)n, f) != (size_t)n) { free(meta); return 0; }

    /* This buffer is intentionally never freed here -- every loaded
     * KmerEntry's hits[] pointer is a direct slice into it (no per-entry
     * copy), so it has to live for the rest of the process, same as
     * everything else the index touches. */
    SeedHit *hitsbuf = xmalloc(sizeof(SeedHit) * (size_t)(total_hits > 0 ? total_hits : 1));
    if (total_hits > 0 && fread(hitsbuf, sizeof(SeedHit), (size_t)total_hits, f) != (size_t)total_hits) {
        free(meta); free(hitsbuf); return 0;
    }

    long ho = 0;
    for (long i = 0; i < n; i++) {
        uint32_t kmer = meta[i].kmer; int n_hits = meta[i].n_hits;
        if (n_hits < 0 || n_hits > MAX_SEED_HITS_PER_KMER) { free(meta); return 0; }
        KmerEntry *e = entry_arena_alloc();
        e->kmer = kmer;
        e->n_hits = n_hits;
        e->hits_cap = n_hits;
        e->hits = &hitsbuf[ho];
        ho += n_hits;
        uint32_t h = hash_kmer(kmer, mask);
        e->next = table[h];
        table[h] = e;
    }
    free(meta);
    return 1;
}

/* Returns 1 and populates chroms[]/n_chroms/both k-mer tables on a cache
 * hit; returns 0 (leaving global state untouched) on any miss or read
 * failure, so the caller can always fall back to the normal FASTA-parse +
 * index-build path. */
static int try_load_index_cache(const char *ref_path) {
    long size, mtime;
    if (!stat_file(ref_path, &size, &mtime)) return 0;
    char cpath[2048];
    cache_path_for(ref_path, cpath, sizeof(cpath));
    FILE *f = fopen(cpath, "rb");
    if (!f) return 0;
    static char rdbuf[1 << 20];
    setvbuf(f, rdbuf, _IOFBF, sizeof(rdbuf));

    IndexCacheHeader hdr;
    int ok = (fread(&hdr, sizeof(hdr), 1, f) == 1) &&
             memcmp(hdr.magic, "RSPKIDX1", 8) == 0 &&
             hdr.ref_size == size && hdr.ref_mtime == mtime &&
             hdr.kmer_len == KMER_LEN && hdr.splice_kmer_len == SPLICE_KMER_LEN &&
             hdr.table_bits == KMER_TABLE_BITS && hdr.splice_table_bits == SPLICE_TABLE_BITS &&
             hdr.max_hits_per_kmer == MAX_SEED_HITS_PER_KMER &&
             hdr.n_chroms > 0 && hdr.n_chroms <= MAX_CHROMS;
    if (ok) {
        n_chroms = hdr.n_chroms;
        for (int i = 0; i < n_chroms && ok; i++) ok = read_chrom(f, &chroms[i]);
        if (ok) ok = read_table(f, kmer_table, KMER_TABLE_MASK);
        if (ok) ok = read_table(f, kmer_table_splice, SPLICE_TABLE_MASK);
    }
    fclose(f);
    if (!ok) { n_chroms = 0; } /* don't leave a half-populated state on a corrupt/truncated cache file */
    return ok;
}

static void save_index_cache(const char *ref_path) {
    long size, mtime;
    if (!stat_file(ref_path, &size, &mtime)) return; /* can't fingerprint -> skip caching, not fatal */
    char cpath[2048], tmp_path[2080];
    cache_path_for(ref_path, cpath, sizeof(cpath));
    snprintf(tmp_path, sizeof(tmp_path), "%s.tmp%d", cpath, (int)getpid());

    FILE *f = fopen(tmp_path, "wb");
    if (!f) return; /* e.g. read-only reference directory -- silently skip, not fatal */
    static char wrbuf[1 << 20];
    setvbuf(f, wrbuf, _IOFBF, sizeof(wrbuf));

    IndexCacheHeader hdr;
    memset(&hdr, 0, sizeof(hdr));
    memcpy(hdr.magic, "RSPKIDX1", 8);
    hdr.ref_size = size; hdr.ref_mtime = mtime;
    hdr.kmer_len = KMER_LEN; hdr.splice_kmer_len = SPLICE_KMER_LEN;
    hdr.table_bits = KMER_TABLE_BITS; hdr.splice_table_bits = SPLICE_TABLE_BITS;
    hdr.max_hits_per_kmer = MAX_SEED_HITS_PER_KMER;
    hdr.n_chroms = n_chroms;
    fwrite(&hdr, sizeof(hdr), 1, f);
    for (int i = 0; i < n_chroms; i++) write_chrom(f, &chroms[i]);
    write_table(f, kmer_table, KMER_TABLE_SIZE);
    write_table(f, kmer_table_splice, SPLICE_TABLE_SIZE);
    int write_ok = !ferror(f);
    fclose(f);

    /* Write to a temp file and rename() into place (atomic on the same
     * filesystem) so a reader never sees a partially-written cache file,
     * and a run killed mid-write can't corrupt a previously good cache. */
    if (write_ok) rename(tmp_path, cpath);
    else unlink(tmp_path);
}

/* ---------------------------------------------------------------------- */
/* GTF annotation parsing (gene-level features only)                      */
/* ---------------------------------------------------------------------- */

static void extract_attr(const char *attrs, const char *key, char *out, size_t outsz) {
    const char *p = strstr(attrs, key);
    out[0] = '\0';
    if (!p) return;
    p += strlen(key);
    while (*p == ' ' || *p == '"') p++;
    size_t i = 0;
    while (*p && *p != '"' && *p != ';' && i < outsz - 1) out[i++] = *p++;
    out[i] = '\0';
}

static void load_gtf(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open GTF annotation");

    char line[MAX_LINE];
    while (fgets(line, sizeof(line), f)) {
        rstrip(line);
        if (line[0] == '#' || line[0] == '\0') continue;

        char chrom[MAX_SEQNAME], source[64], feature[64], strand_s[8], frame_s[8], attrs[MAX_LINE];
        long start, end;
        char score_s[16];

        int n = sscanf(line, "%127s\t%63s\t%63s\t%ld\t%ld\t%15s\t%7s\t%7s\t%[^\n]",
                        chrom, source, feature, &start, &end, score_s, strand_s, frame_s, attrs);
        if (n < 9) continue;
        /* Accept any feature type that carries gene-body coordinates, not
         * just an explicit "gene" line: many real-world GTFs omit gene-level
         * rows for some or all genes (confirmed against nf-core/rnaseq's own
         * CI test annotation, genes_with_empty_tid.gtf, which has an
         * explicit "gene" row for only 1 of its 124 genes -- the rest only
         * have transcript/exon/CDS/codon rows). Gene intervals are derived
         * as the min-start/max-end across every such line sharing a
         * gene_id, which also correctly reproduces the explicit-"gene"-line
         * case (a single line trivially has itself as its own min/max). */
        if (strcmp(feature, "gene") != 0 && strcmp(feature, "transcript") != 0 &&
            strcmp(feature, "exon") != 0 && strcmp(feature, "CDS") != 0 &&
            strcmp(feature, "start_codon") != 0 && strcmp(feature, "stop_codon") != 0) continue;

        char gene_id[MAX_SEQNAME];
        extract_attr(attrs, "gene_id", gene_id, sizeof(gene_id));
        if (gene_id[0] == '\0') continue; /* no gene_id attribute on this line -- skip */

        int idx = -1;
        for (int i = 0; i < n_genes; i++) if (strcmp(genes[i].gene_id, gene_id) == 0) { idx = i; break; }
        if (idx < 0) {
            if (n_genes >= MAX_GENES) die("too many genes for demo buffer");
            idx = n_genes++;
            Gene *g = &genes[idx];
            strncpy(g->gene_id, gene_id, sizeof(g->gene_id) - 1);
            g->gene_id[sizeof(g->gene_id) - 1] = '\0';
            strncpy(g->chrom, chrom, sizeof(g->chrom) - 1);
            g->chrom[sizeof(g->chrom) - 1] = '\0';
            g->start = start;
            g->end = end;
            g->strand = strand_s[0];
            g->unique_count = 0;
            g->effective_count = 0.0;
        } else {
            if (start < genes[idx].start) genes[idx].start = start;
            if (end > genes[idx].end) genes[idx].end = end;
        }
    }
    fclose(f);
    if (n_genes == 0) die("no gene-bearing features (gene/transcript/exon/CDS/...) found in GTF");
}

/* ---------------------------------------------------------------------- */
/* GeneSet helpers                                                        */
/* ---------------------------------------------------------------------- */

static void geneset_add(GeneSet *s, int idx) {
    for (int i = 0; i < s->n; i++) if (s->idx[i] == idx) return;
    if (s->n < MAX_GENES_PER_UNIT) s->idx[s->n++] = idx;
}

static void geneset_union(const GeneSet *a, const GeneSet *b, GeneSet *out) {
    out->n = 0;
    for (int i = 0; i < a->n; i++) geneset_add(out, a->idx[i]);
    for (int i = 0; i < b->n; i++) geneset_add(out, b->idx[i]);
}

static void geneset_intersect(const GeneSet *a, const GeneSet *b, GeneSet *out) {
    out->n = 0;
    for (int i = 0; i < a->n; i++) {
        for (int j = 0; j < b->n; j++) {
            if (a->idx[i] == b->idx[j]) { geneset_add(out, a->idx[i]); break; }
        }
    }
}

static void collect_genes_for_hit(const Hit *h, GeneSet *out) {
    const char *chrom_name = chroms[h->chrom_idx].name;
    long left = h->ref_start + 1;
    long right = h->ref_start; /* will become 1-based inclusive end below */
    for (int k = 0; k < h->n_cigar; k++) {
        char op = h->cigar[k].op;
        if (op == 'M' || op == 'D' || op == 'N') right += h->cigar[k].len;
    }
    for (int i = 0; i < n_genes; i++) {
        if (strcmp(genes[i].chrom, chrom_name) != 0) continue;
        if (left <= genes[i].end && right >= genes[i].start) geneset_add(out, i);
    }
}

static void collect_candidate_genes_for_read(const Read *r, GeneSet *out) {
    out->n = 0;
    if (!r->aln.mapped) return;
    for (int i = 0; i < r->aln.n_hits; i++) {
        GeneSet tmp; tmp.n = 0;
        collect_genes_for_hit(&r->aln.hits[i], &tmp);
        GeneSet merged; geneset_union(out, &tmp, &merged);
        *out = merged;
    }
}

static void collect_candidate_genes_for_fragment(const Read *r1, const Read *r2, GeneSet *out) {
    GeneSet g1 = {.n = 0}, g2 = {.n = 0};
    if (r1->aln.mapped) collect_candidate_genes_for_read(r1, &g1);
    if (r2->aln.mapped) collect_candidate_genes_for_read(r2, &g2);

    if (r1->aln.mapped && r2->aln.mapped) {
        GeneSet inter; geneset_intersect(&g1, &g2, &inter);
        if (inter.n > 0) *out = inter;
        else geneset_union(&g1, &g2, out);
    } else if (r1->aln.mapped) {
        *out = g1;
    } else if (r2->aln.mapped) {
        *out = g2;
    } else {
        out->n = 0;
    }
}

/* ---------------------------------------------------------------------- */
/* FASTQ parsing (single-end and paired-end)                              */
/* ---------------------------------------------------------------------- */

static void init_read_slot(Read *r, const char *id_line, const char *seq_line,
                            const char *qual_line, int read_num, int mate_idx) {
    memset(r, 0, sizeof(Read));
    char idbuf[MAX_LINE];
    strncpy(idbuf, id_line + 1, sizeof(idbuf) - 1);
    idbuf[sizeof(idbuf) - 1] = '\0';
    /* truncate at the first whitespace: standard FASTQ/SAM convention.
     * Real sequencer/SRA headers routinely carry extra descriptive text
     * after a space (e.g. "SRR123.456 456/1 kraken:taxid|9606"), which
     * must NOT become part of the read ID or it won't match other tools'
     * output (or a paired mate's ID). */
    for (char *p = idbuf; *p; p++) if (*p == ' ' || *p == '\t') { *p = '\0'; break; }
    size_t l = strlen(idbuf);
    if (l >= 2 && idbuf[l-2] == '/' && (idbuf[l-1] == '1' || idbuf[l-1] == '2')) idbuf[l-2] = '\0';
    strncpy(r->id, idbuf, sizeof(r->id) - 1);

    int len = (int)strlen(seq_line);
    if (len >= MAX_READ_LEN) len = MAX_READ_LEN - 1;
    r->seq = xmalloc(len + 1);
    r->qual = xmalloc(len + 1);
    memcpy(r->seq, seq_line, len);   r->seq[len] = '\0';
    memcpy(r->qual, qual_line, len); r->qual[len] = '\0';
    r->raw_len = len;
    r->read_num = read_num;
    r->mate_idx = mate_idx;
}

static void load_fastq_se(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open FASTQ reads file");
    reads = xmalloc(sizeof(Read) * MAX_READS);

    char l1[MAX_LINE], l2[MAX_LINE], l3[MAX_LINE], l4[MAX_LINE];
    while (fgets(l1, sizeof(l1), f)) {
        if (!fgets(l2, sizeof(l2), f)) break;
        if (!fgets(l3, sizeof(l3), f)) break;
        if (!fgets(l4, sizeof(l4), f)) break;
        rstrip(l1); rstrip(l2); rstrip(l3); rstrip(l4);
        if (l1[0] != '@') continue;
        if (n_reads >= MAX_READS) die("too many reads for demo buffer (increase MAX_READS)");
        init_read_slot(&reads[n_reads], l1, l2, l4, 0, -1);
        n_reads++;
    }
    fclose(f);
    if (n_reads == 0) die("no reads found in FASTQ file");
}

static void load_fastq_pe(const char *path1, const char *path2) {
    FILE *f1 = fopen(path1, "r");
    if (!f1) die("cannot open R1 FASTQ file");
    FILE *f2 = fopen(path2, "r");
    if (!f2) die("cannot open R2 FASTQ file");

    reads = xmalloc(sizeof(Read) * MAX_READS);
    char a1[MAX_LINE], a2[MAX_LINE], a3[MAX_LINE], a4[MAX_LINE];
    char b1[MAX_LINE], b2[MAX_LINE], b3[MAX_LINE], b4[MAX_LINE];
    int mismatched_ids = 0;

    while (fgets(a1, sizeof(a1), f1)) {
        if (!fgets(a2, sizeof(a2), f1) || !fgets(a3, sizeof(a3), f1) || !fgets(a4, sizeof(a4), f1)) break;
        if (!fgets(b1, sizeof(b1), f2) || !fgets(b2, sizeof(b2), f2) ||
            !fgets(b3, sizeof(b3), f2) || !fgets(b4, sizeof(b4), f2)) {
            fprintf(stderr, "WARNING: R2 file ran out of reads before R1; truncating pairs here.\n");
            break;
        }
        rstrip(a1); rstrip(a2); rstrip(a4);
        rstrip(b1); rstrip(b2); rstrip(b4);
        if (a1[0] != '@' || b1[0] != '@') continue;
        if (n_reads + 2 > MAX_READS) die("too many reads for demo buffer (increase MAX_READS)");

        int i1 = n_reads, i2 = n_reads + 1;
        init_read_slot(&reads[i1], a1, a2, a4, 1, i2);
        init_read_slot(&reads[i2], b1, b2, b4, 2, i1);
        if (strcmp(reads[i1].id, reads[i2].id) != 0) mismatched_ids++;
        n_reads += 2;
    }
    fclose(f1); fclose(f2);
    if (n_reads == 0) die("no read pairs found in FASTQ files");
    if (mismatched_ids)
        fprintf(stderr, "WARNING: %d read pair(s) had mismatched IDs between R1/R2 (order-based pairing used anyway)\n",
                mismatched_ids);
}

/* ---------------------------------------------------------------------- */
/* QC (pre-trim)                                                          */
/* ---------------------------------------------------------------------- */

static void compute_raw_qc(Read *r) {
    int gc = 0, n = 0;
    double qsum = 0;
    for (int i = 0; i < r->raw_len; i++) {
        char b = toupper((unsigned char)r->seq[i]);
        if (b == 'G' || b == 'C') gc++;
        if (b == 'N') n++;
        qsum += phred_to_prob_correct(r->qual[i]);
    }
    r->gc_pct = r->raw_len ? (100.0 * gc / r->raw_len) : 0.0;
    r->n_count = n;
    r->mean_q_raw = r->raw_len ? (qsum / r->raw_len) : 0.0;
    #pragma omp atomic
    total_raw_bases += r->raw_len;
}

/* ---------------------------------------------------------------------- */
/* Trimming: adapter removal + 3' sliding-window quality trim             */
/* ---------------------------------------------------------------------- */

static void trim_read(Read *r) {
    int len = r->raw_len;

    char *hit = strstr(r->seq, ADAPTER_SEQ);
    if (hit) {
        int cut = (int)(hit - r->seq);
        if (cut < len) {
            len = cut; r->adapter_trimmed = 1;
            #pragma omp atomic
            reads_with_adapter++;
        }
    }

    int len_after_adapter = len;
    while (len >= TRIM_WINDOW) {
        double wsum = 0;
        for (int i = len - TRIM_WINDOW; i < len; i++) wsum += phred_to_prob_correct(r->qual[i]);
        if (wsum / TRIM_WINDOW < QUALITY_TRIM_THRESHOLD) len--;
        else break;
    }
    if (len < len_after_adapter) {
        r->qual_trimmed_bp = len_after_adapter - len;
        #pragma omp atomic
        reads_quality_trimmed++;
    }

    r->trimmed_len = len;
    #pragma omp atomic
    total_trimmed_bases += len;
}

/* ---------------------------------------------------------------------- */
/* Alignment                                                              */
/* ---------------------------------------------------------------------- */

static int hamming_at(const Chrom *c, long pos, const char *query, int qlen, int max_mm) {
    if (pos < 0 || pos + qlen > c->len) return -1;
    int mm = 0;
    for (int i = 0; i < qlen; i++) {
        char a = toupper((unsigned char)c->seq[pos + i]);
        char b = toupper((unsigned char)query[i]);
        if (a != b) { mm++; if (mm > max_mm) return -1; }
    }
    return mm;
}

/* ---------------------------------------------------------------------- */
/* Smith-Waterman-Gotoh local alignment (affine gap penalties).           */
/*                                                                        */
/* Standard 3-matrix Gotoh formulation: H = best local score ending here, */
/* E = best score ending in a gap consumed from the reference (a          */
/* Deletion from the read's perspective), F = best score ending in a gap  */
/* consumed from the query (an Insertion). Traceback reconstructs a real  */
/* CIGAR (M/I/D), with any unaligned read prefix/suffix reported as S     */
/* (soft-clip) -- this is what lets reads with small indels, or a         */
/* chimeric/adapter-contaminated end, align at all, unlike the pure       */
/* Hamming-distance approach used in earlier versions of this file.       */
/* ---------------------------------------------------------------------- */

/* Reusable DP scratch buffers, one set per OpenMP thread (threadprivate so
 * concurrent threads never share a buffer). Sized once for the largest
 * window this pipeline can ever build and reused across every call for the
 * life of the process, instead of six malloc+free per candidate per read.
 * At up to SW_MAX_CANDIDATES(8) candidates x 2 strands x every read that
 * doesn't take the ungapped fast path below, that malloc/free churn was
 * large enough to show up directly as minor-page-fault-driven system time
 * (profiled: ~4.3M minor page faults for a 20,000-read-pair run), not just
 * allocator CPU cost. */
static int    *sw_H, *sw_E, *sw_F;
static unsigned char *sw_tbH, *sw_tbE, *sw_tbF;
static int     sw_scratch_cols = 0;
#pragma omp threadprivate(sw_H, sw_E, sw_F, sw_tbH, sw_tbE, sw_tbF, sw_scratch_cols)

static void sw_scratch_ensure(int cols_needed) {
    if (sw_scratch_cols >= cols_needed) return;
    size_t n = (size_t)(MAX_READ_LEN + 1) * cols_needed;
    sw_H   = xrealloc(sw_H,  sizeof(int) * n);
    sw_E   = xrealloc(sw_E,  sizeof(int) * n);
    sw_F   = xrealloc(sw_F,  sizeof(int) * n);
    sw_tbH = xrealloc(sw_tbH, n);
    sw_tbE = xrealloc(sw_tbE, n);
    sw_tbF = xrealloc(sw_tbF, n);
    sw_scratch_cols = cols_needed;
}

static int sw_local_align(const Chrom *chrom, long window_start, long window_len,
                           const char *query, int qlen, char strand, Hit *out) {
    int rows = qlen + 1, cols = (int)window_len + 1;
    sw_scratch_ensure(cols);
    int *H = sw_H, *E = sw_E, *F = sw_F;
    unsigned char *tbH = sw_tbH, *tbE = sw_tbE, *tbF = sw_tbF;
    const int NEG = -1000000;
#define SWIDX(i,j) ((size_t)(i) * cols + (j))

    /* Banded init: only the cells the banded recurrence below can actually
     * touch need a defined value (H=0, E/F=NEG, "not part of any path").
     * That's the band itself plus one column of padding on each side (the
     * recurrence reads one column left/right of the band when it's at the
     * band edge). Previously this initialized the FULL qlen x window_len
     * rectangle unconditionally, which -- since window_len is only ~2x the
     * band width -- meant the "banded" DP still paid close to full-rectangle
     * cost merely to initialize, even though the recurrence itself was
     * banded. Row 0 (the DP boundary row) is handled explicitly since it
     * has no i-1 band to inherit padding from. */
    int diag_center = SW_PAD;
    for (int j = 0; j <= (int)window_len; j++) {
        H[SWIDX(0,j)] = 0; E[SWIDX(0,j)] = NEG; F[SWIDX(0,j)] = NEG; tbH[SWIDX(0,j)] = 0;
    }
    for (int i = 1; i <= qlen; i++) {
        int j_lo = i + diag_center - SW_BAND - 1; if (j_lo < 0) j_lo = 0;
        int j_hi = i + diag_center + SW_BAND + 1; if (j_hi > (int)window_len) j_hi = (int)window_len;
        for (int j = j_lo; j <= j_hi; j++) {
            H[SWIDX(i,j)] = 0; E[SWIDX(i,j)] = NEG; F[SWIDX(i,j)] = NEG; tbH[SWIDX(i,j)] = 0;
        }
    }

    /* Band the recurrence around the seed-implied zero-indel diagonal
     * (j = i + SW_PAD, given how window_start = implied_start - SW_PAD was
     * chosen by the caller) instead of computing the full qlen x window_len
     * rectangle. This is the dominant cost of alignment -- restricting it
     * to +/- SW_BAND around that diagonal cuts DP work ~6x for a typical
     * 150bp read with no loss of alignments the window could represent in
     * the first place for realistic same-species short-read indel sizes
     * (verified by full regression re-run against the un-banded baseline;
     * see BANDED_SW_OPTIMIZATION.md). */
    int best_score = 0, best_i = 0, best_j = 0;
    for (int i = 1; i <= qlen; i++) {
        char qc = toupper((unsigned char)query[i - 1]);
        int j_lo = i + diag_center - SW_BAND; if (j_lo < 1) j_lo = 1;
        int j_hi = i + diag_center + SW_BAND; if (j_hi > (int)window_len) j_hi = (int)window_len;
        for (int j = j_lo; j <= j_hi; j++) {
            char rc = toupper((unsigned char)chrom->seq[window_start + j - 1]);
            int sub = (qc == rc) ? SW_MATCH : SW_MISMATCH;

            int e_open = H[SWIDX(i,j-1)] + SW_GAP_OPEN, e_ext = E[SWIDX(i,j-1)] + SW_GAP_EXTEND;
            if (e_open >= e_ext) { E[SWIDX(i,j)] = e_open; tbE[SWIDX(i,j)] = 0; }
            else                 { E[SWIDX(i,j)] = e_ext;  tbE[SWIDX(i,j)] = 1; }

            int f_open = H[SWIDX(i-1,j)] + SW_GAP_OPEN, f_ext = F[SWIDX(i-1,j)] + SW_GAP_EXTEND;
            if (f_open >= f_ext) { F[SWIDX(i,j)] = f_open; tbF[SWIDX(i,j)] = 0; }
            else                 { F[SWIDX(i,j)] = f_ext;  tbF[SWIDX(i,j)] = 1; }

            int diag = H[SWIDX(i-1,j-1)] + sub;
            int best = 0; unsigned char tb = 0;
            if (diag > best)          { best = diag;          tb = 1; }
            if (E[SWIDX(i,j)] > best) { best = E[SWIDX(i,j)];  tb = 2; }
            if (F[SWIDX(i,j)] > best) { best = F[SWIDX(i,j)];  tb = 3; }
            H[SWIDX(i,j)] = best; tbH[SWIDX(i,j)] = tb;
            if (best > best_score) { best_score = best; best_i = i; best_j = j; }
        }
    }

    if (best_score <= 0) { return 0; }

    CigarOp rev_ops[MAX_CIGAR_OPS]; int n_rev = 0;
    int i = best_i, j = best_j, mode = 0; /* mode: 0=H, 1=E(deletion), 2=F(insertion) */
    int mismatches = 0;

    while (i > 0 && j > 0 && n_rev < MAX_CIGAR_OPS - 1) {
        if (mode == 0) {
            if (H[SWIDX(i,j)] == 0) break;
            unsigned char tb = tbH[SWIDX(i,j)];
            if (tb == 1) {
                char qc = toupper((unsigned char)query[i - 1]);
                char rcb = toupper((unsigned char)chrom->seq[window_start + j - 1]);
                if (qc != rcb) mismatches++;
                if (n_rev > 0 && rev_ops[n_rev-1].op == 'M') rev_ops[n_rev-1].len++;
                else { rev_ops[n_rev].op = 'M'; rev_ops[n_rev].len = 1; n_rev++; }
                i--; j--;
            } else if (tb == 2) { mode = 1; }
            else { mode = 2; }
        } else if (mode == 1) { /* deletion: ref consumed, query not */
            if (n_rev > 0 && rev_ops[n_rev-1].op == 'D') rev_ops[n_rev-1].len++;
            else { rev_ops[n_rev].op = 'D'; rev_ops[n_rev].len = 1; n_rev++; }
            mismatches++;
            unsigned char tb = tbE[SWIDX(i,j)];
            j--; mode = (tb == 0) ? 0 : 1;
        } else { /* insertion: query consumed, ref not */
            if (n_rev > 0 && rev_ops[n_rev-1].op == 'I') rev_ops[n_rev-1].len++;
            else { rev_ops[n_rev].op = 'I'; rev_ops[n_rev].len = 1; n_rev++; }
            mismatches++;
            unsigned char tb = tbF[SWIDX(i,j)];
            i--; mode = (tb == 0) ? 0 : 2;
        }
    }
    int start_i = i;

#undef SWIDX

    out->n_cigar = 0;
    if (start_i > 0) { out->cigar[out->n_cigar].op = 'S'; out->cigar[out->n_cigar].len = start_i; out->n_cigar++; }
    for (int k = n_rev - 1; k >= 0 && out->n_cigar < MAX_CIGAR_OPS; k--) {
        if (out->n_cigar > 0 && out->cigar[out->n_cigar-1].op == rev_ops[k].op) out->cigar[out->n_cigar-1].len += rev_ops[k].len;
        else { out->cigar[out->n_cigar] = rev_ops[k]; out->n_cigar++; }
    }
    int trailing_clip = qlen - best_i;
    if (trailing_clip > 0 && out->n_cigar < MAX_CIGAR_OPS) {
        if (out->n_cigar > 0 && out->cigar[out->n_cigar-1].op == 'S') out->cigar[out->n_cigar-1].len += trailing_clip;
        else { out->cigar[out->n_cigar].op = 'S'; out->cigar[out->n_cigar].len = trailing_clip; out->n_cigar++; }
    }

    out->ref_start = window_start + j; /* j now holds the ref offset where alignment starts */
    out->score = best_score;
    out->mismatches = mismatches;
    out->strand = strand;
    out->spliced = 0;
    return 1;
}

/* Gathers candidate seed positions from the k-mer index (same 3-offset
 * strategy as before), runs full local Smith-Waterman around each
 * (capped at SW_MAX_CANDIDATES to bound worst-case cost in repetitive
 * regions), and returns every candidate tied for the best score, up to
 * MAX_MULTI_HITS -- this is what preserves genuine multi-mapper detection
 * (see the EM quantification stage) while still supporting indels. */
typedef struct { int chrom_idx; long implied_start; } SwCandidate;

static int try_sw_align_multi(const char *query, int qlen, char strand,
                               Hit *out_hits, int *n_out, int *best_score_out) {
    if (qlen < KMER_LEN) return 0;

    int seed_offsets[3], n_seeds = 0;
    seed_offsets[n_seeds++] = 0;
    if (qlen - KMER_LEN > 0) {
        int mid = (qlen - KMER_LEN) / 2;
        if (mid != 0) seed_offsets[n_seeds++] = mid;
        seed_offsets[n_seeds++] = qlen - KMER_LEN;
    }

    SwCandidate cands[MAX_CAND_POS];  /* was `static` -- shared/racy across OpenMP threads;
                                        * this is only 64*16=1KB, cheap on the stack. */
    int n_cand = 0;
    for (int s = 0; s < n_seeds; s++) {
        int so = seed_offsets[s];
        uint32_t kmer;
        if (!encode_kmer(query + so, KMER_LEN, &kmer)) continue;
        KmerEntry *e = kmer_index_lookup(kmer_table, KMER_TABLE_MASK, kmer);
        if (!e) continue;
        for (int h = 0; h < e->n_hits; h++) {
            long implied_start = e->hits[h].pos - so;
            int dup = 0;
            for (int c = 0; c < n_cand; c++)
                if (cands[c].chrom_idx == e->hits[h].chrom_idx &&
                    labs(cands[c].implied_start - implied_start) < SW_PAD) { dup = 1; break; }
            if (dup) continue;
            if (n_cand < MAX_CAND_POS) { cands[n_cand].chrom_idx = e->hits[h].chrom_idx; cands[n_cand].implied_start = implied_start; n_cand++; }
        }
    }
    if (n_cand == 0) return 0;
    if (n_cand > SW_MAX_CANDIDATES) n_cand = SW_MAX_CANDIDATES;

    int min_accept = (int)(qlen * SW_MATCH * SW_MIN_SCORE_FRAC);
    int best_score = -1, n_best = 0;
    Hit best_hits[MAX_CAND_POS];

    for (int c = 0; c < n_cand; c++) {
        Chrom *chrom = &chroms[cands[c].chrom_idx];

        /* Ungapped fast path, PERFECT MATCHES ONLY (mm == 0). A full-length
         * exact match already scores the theoretical maximum (qlen *
         * SW_MATCH) -- nothing (gap, soft-clip, or otherwise) can ever
         * score higher than that, so this is not an approximation, it's
         * exactly the same optimum full Smith-Waterman DP would find, just
         * without paying for the DP.
         *
         * NOTE: an earlier version of this fast path also special-cased
         * mm==1 on the reasoning that the cheapest possible gap (cost 6)
         * can't beat a single mismatch (cost 4). That reasoning holds for
         * gapped alternatives but misses soft-clipping: a mismatch sitting
         * at either read terminus can be dropped for a *soft-clip cost of
         * 0* (score (qlen-1)*2) rather than paid for as a substitution
         * (score 2*qlen-6) -- clipping wins whenever a same-or-better
         * placement exists 1bp shifted with the terminal base clipped.
         * This was caught by diffing against the un-changed baseline on
         * real data (see CONSOLIDATED_REPORT.md's validation section):
         * several reads that the baseline correctly reported as e.g.
         * "1S149M NM:i:0" were mis-scored by the mm<=1 fast path as
         * "150M NM:i:1", a real, lower-scoring answer. Restricting the
         * fast path to mm==0 removes the bug entirely (a perfect match
         * has no terminus mismatch to clip around) at the cost of some
         * candidates with a single interior mismatch now taking the full
         * banded-DP path below instead of the fast path -- correctness
         * over speed. */
        int hmm = hamming_at(chrom, cands[c].implied_start, query, qlen, 0);
        if (hmm == 0) {
            Hit h;
            h.chrom_idx = cands[c].chrom_idx;
            h.ref_start = cands[c].implied_start;
            h.n_cigar = 1;
            h.cigar[0].op = 'M'; h.cigar[0].len = qlen;
            h.spliced = 0;
            h.score = qlen * SW_MATCH - hmm * (SW_MATCH - SW_MISMATCH);
            h.mismatches = hmm;
            h.strand = strand;
            if (h.score >= min_accept && (double)h.mismatches / qlen <= SW_MAX_EDIT_FRAC) {
                if (h.score > best_score) { best_score = h.score; best_hits[0] = h; n_best = 1; }
                else if (h.score == best_score && n_best < MAX_CAND_POS) { best_hits[n_best++] = h; }
            }
            continue;
        }

        long window_start = cands[c].implied_start - SW_PAD;
        if (window_start < 0) window_start = 0;
        long window_end = cands[c].implied_start + qlen + SW_PAD;
        if (window_end > chrom->len) window_end = chrom->len;
        long window_len = window_end - window_start;
        if (window_len <= 0) continue;

        Hit h;
        if (!sw_local_align(chrom, window_start, window_len, query, qlen, strand, &h)) continue;
        h.chrom_idx = cands[c].chrom_idx;
        if (h.score < min_accept) continue;
        if ((double)h.mismatches / qlen > SW_MAX_EDIT_FRAC) continue;

        if (h.score > best_score) { best_score = h.score; best_hits[0] = h; n_best = 1; }
        else if (h.score == best_score && n_best < MAX_CAND_POS) { best_hits[n_best++] = h; }
    }

    if (best_score < 0 || n_best == 0) return 0;
    int n = (n_best < MAX_MULTI_HITS) ? n_best : MAX_MULTI_HITS;
    for (int k = 0; k < n; k++) out_hits[k] = best_hits[k];
    *n_out = n;
    *best_score_out = best_score;
    return 1;
}

/* Splice-aware two-anchor alignment (single best hit only -- multi-mapping
 * is not attempted for spliced alignments; see README limitations).
 *
 * Anchors use SPLICE_KMER_LEN (shorter than the main index's KMER_LEN) so
 * that junction-spanning reads with only a short exon overhang on one side
 * can still be seeded -- the earlier version required a full KMER_LEN
 * (16bp) exact match on BOTH sides, which meant any junction close to
 * either read edge was simply unseedable and silently missed. */
static int try_spliced_align(const char *query, int qlen, char strand, Hit *out) {
    if (qlen < 2 * SPLICE_KMER_LEN + 4) return 0;

    uint32_t kmer_l, kmer_r;
    if (!encode_kmer(query, SPLICE_KMER_LEN, &kmer_l)) return 0;
    if (!encode_kmer(query + qlen - SPLICE_KMER_LEN, SPLICE_KMER_LEN, &kmer_r)) return 0;

    KmerEntry *el = kmer_index_lookup(kmer_table_splice, SPLICE_TABLE_MASK, kmer_l);
    KmerEntry *er = kmer_index_lookup(kmer_table_splice, SPLICE_TABLE_MASK, kmer_r);
    if (!el || !er) return 0;

    int found = 0, best_mm = SPLICE_MAX_MISMATCHES + 1, best_canon = 0;
    int best_chrom = -1, best_b = -1;
    long best_c1 = -1, best_intron = -1;

    for (int i = 0; i < el->n_hits; i++) {
        int ci = el->hits[i].chrom_idx;
        long c1 = el->hits[i].pos;
        Chrom *chrom = &chroms[ci];

        for (int j = 0; j < er->n_hits; j++) {
            if (er->hits[j].chrom_idx != ci) continue;
            long c2 = er->hits[j].pos;
            long intron = (c2 - c1) - (long)(qlen - SPLICE_KMER_LEN);
            if (intron < MIN_INTRON || intron > MAX_INTRON) continue;
            if (c1 < 0) continue;

            for (int b = SPLICE_KMER_LEN; b <= qlen - SPLICE_KMER_LEN; b++) {
                int mm1 = hamming_at(chrom, c1, query, b, SPLICE_MAX_MISMATCHES);
                if (mm1 < 0) continue;
                long exon2_start = c1 + b + intron;
                int mm2 = hamming_at(chrom, exon2_start, query + b, qlen - b, SPLICE_MAX_MISMATCHES - mm1);
                if (mm2 < 0) continue;
                int total = mm1 + mm2;
                if (total > SPLICE_MAX_MISMATCHES) continue;

                int canon = 0;
                long donor_pos = c1 + b;
                long acceptor_end = exon2_start;
                if (donor_pos + 1 < chrom->len && acceptor_end - 2 >= 0) {
                    char d1 = toupper((unsigned char)chrom->seq[donor_pos]);
                    char d2 = toupper((unsigned char)chrom->seq[donor_pos + 1]);
                    char a1 = toupper((unsigned char)chrom->seq[acceptor_end - 2]);
                    char a2 = toupper((unsigned char)chrom->seq[acceptor_end - 1]);
                    if (d1 == 'G' && d2 == 'T' && a1 == 'A' && a2 == 'G') canon = 1;
                }

                int better = !found || (canon && !best_canon) || (canon == best_canon && total < best_mm);
                if (better) {
                    found = 1; best_mm = total; best_canon = canon;
                    best_chrom = ci; best_c1 = c1; best_intron = intron; best_b = b;
                }
            }
        }
    }

    if (!found) return 0;
    out->chrom_idx = best_chrom;
    out->ref_start = best_c1;
    out->n_cigar = 0;
    out->cigar[out->n_cigar].op = 'M'; out->cigar[out->n_cigar].len = best_b; out->n_cigar++;
    out->cigar[out->n_cigar].op = 'N'; out->cigar[out->n_cigar].len = (int)best_intron; out->n_cigar++;
    out->cigar[out->n_cigar].op = 'M'; out->cigar[out->n_cigar].len = qlen - best_b; out->n_cigar++;
    out->spliced = 1;
    out->intron_len = best_intron;
    out->canonical_splice = best_canon;
    out->score = 0; /* splice path is Hamming-scored, not SW-scored; not compared against SW scores */
    out->mismatches = best_mm;
    out->strand = strand;
    return 1;
}

/* Genomic span of a Hit, computed from its CIGAR: M/D/N all consume
 * reference bases (and so extend the span); I and S do not. */
static long hit_span_start0(const Hit *h) { return h->ref_start; }
static long hit_span_end0(const Hit *h) {
    long end = h->ref_start;
    for (int k = 0; k < h->n_cigar; k++) {
        char op = h->cigar[k].op;
        if (op == 'M' || op == 'D' || op == 'N') end += h->cigar[k].len;
    }
    return end - 1;
}

static int hit_aligned_bases(const Hit *h) {
    int n = 0;
    for (int k = 0; k < h->n_cigar; k++) {
        char op = h->cigar[k].op;
        if (op == 'M' || op == 'I') n += h->cigar[k].len;
    }
    return n;
}

/* Aligns one read, trying both orientations (forward and reverse
 * complement, since only the forward reference strand is indexed) via
 * local Smith-Waterman, then -- ONLY if that didn't already explain the
 * entire read (i.e. some soft-clipping remains) -- also tries splice
 * detection. Splice hits always cover the full read by construction, so
 * they're only worth the extra search when SW left something unexplained;
 * and whenever both are available, the one covering more of the read wins
 * (a full-length spliced alignment is strictly more informative than a
 * partial soft-clipped one, even if SW's raw score is nominally similar). */
static void align_read(Read *r) {
    if (r->trimmed_len < MIN_READ_LEN_AFTER_TRIM) {
        #pragma omp atomic
        reads_dropped_too_short++;
        r->aln.mapped = 0;
        return;
    }

    char *rc = revcomp(r->seq, r->trimmed_len);
    Hit fwd_hits[MAX_MULTI_HITS], rev_hits[MAX_MULTI_HITS];
    int n_fwd = 0, n_rev = 0, score_fwd = -1, score_rev = -1;
    int ok_f = try_sw_align_multi(r->seq, r->trimmed_len, '+', fwd_hits, &n_fwd, &score_fwd);
    int ok_r = try_sw_align_multi(rc, r->trimmed_len, '-', rev_hits, &n_rev, &score_rev);

    memset(&r->aln, 0, sizeof(r->aln));

    int best_score = -1, sw_coverage = -1;
    if (ok_f || ok_r) {
        best_score = (ok_f && ok_r) ? (score_fwd > score_rev ? score_fwd : score_rev) : (ok_f ? score_fwd : score_rev);
        const Hit *rep = (ok_f && score_fwd == best_score) ? &fwd_hits[0] : &rev_hits[0];
        sw_coverage = hit_aligned_bases(rep);
    }

    if (sw_coverage < r->trimmed_len) {
        Hit sf, sr;
        int ok_sf = try_spliced_align(r->seq, r->trimmed_len, '+', &sf);
        int ok_sr = try_spliced_align(rc, r->trimmed_len, '-', &sr);
        if (ok_sf || ok_sr) {
            Hit *best = (ok_sf && ok_sr) ? (sf.mismatches <= sr.mismatches ? &sf : &sr) : (ok_sf ? &sf : &sr);
            r->aln.hits = xmalloc(sizeof(Hit));
            r->aln.hits[0] = *best;
            r->aln.n_hits = 1;
            r->aln.mapped = 1;
            #pragma omp atomic
            n_spliced_alignments++;
            free(rc);
            return;
        }
    }

    if (sw_coverage >= 0) {
        int n_avail = 0;
        if (ok_f && score_fwd == best_score) n_avail += (n_fwd < MAX_MULTI_HITS ? n_fwd : MAX_MULTI_HITS);
        if (ok_r && score_rev == best_score) n_avail += (n_rev < MAX_MULTI_HITS ? n_rev : MAX_MULTI_HITS);
        if (n_avail > MAX_MULTI_HITS) n_avail = MAX_MULTI_HITS;
        r->aln.hits = xmalloc(sizeof(Hit) * (n_avail > 0 ? n_avail : 1));
        int n = 0;
        if (ok_f && score_fwd == best_score)
            for (int i = 0; i < n_fwd && n < MAX_MULTI_HITS; i++) r->aln.hits[n++] = fwd_hits[i];
        if (ok_r && score_rev == best_score)
            for (int i = 0; i < n_rev && n < MAX_MULTI_HITS; i++) r->aln.hits[n++] = rev_hits[i];
        r->aln.n_hits = n;
        r->aln.mapped = 1;
        #pragma omp atomic
        n_ungapped_alignments++;
    } else {
        r->aln.mapped = 0;
    }
    free(rc);
}

/* ---------------------------------------------------------------------- */
/* EM-based gene-level quantification (RSEM/Salmon-style multi-mapper     */
/* resolution, simplified: no fragment-length or sequence-bias model).    */
/* ---------------------------------------------------------------------- */

static void quantify_em(void) {
    int step = paired_mode ? 2 : 1;
    int n_units = n_reads / step;

    GeneSet *unit_cand = xmalloc(sizeof(GeneSet) * (size_t)n_units);
    int *unit_mapped = xmalloc(sizeof(int) * (size_t)n_units);

    for (int u = 0; u < n_units; u++) {
        if (paired_mode) {
            Read *r1 = &reads[2*u], *r2 = &reads[2*u + 1];
            unit_mapped[u] = r1->aln.mapped || r2->aln.mapped;
            collect_candidate_genes_for_fragment(r1, r2, &unit_cand[u]);
            r1->candidates = unit_cand[u];
            r2->candidates = unit_cand[u];
        } else {
            Read *r = &reads[u];
            unit_mapped[u] = r->aln.mapped;
            collect_candidate_genes_for_read(r, &unit_cand[u]);
            r->candidates = unit_cand[u];
        }
    }

    double *unique_count = xmalloc(sizeof(double) * (size_t)n_genes);
    int *gene_in_use = xmalloc(sizeof(int) * (size_t)n_genes);
    for (int g = 0; g < n_genes; g++) { unique_count[g] = 0.0; gene_in_use[g] = 0; }

    /* indices of ambiguous (multi-mapping) units, for the EM loop */
    int *multi_units = xmalloc(sizeof(int) * (size_t)n_units);
    int n_multi = 0;

    n_unique_units = n_multi_units = n_no_feature_units = n_unmapped_units = 0;
    for (int u = 0; u < n_units; u++) {
        if (!unit_mapped[u]) { n_unmapped_units++; continue; }
        if (unit_cand[u].n == 0) { n_no_feature_units++; continue; }
        if (unit_cand[u].n == 1) {
            unique_count[unit_cand[u].idx[0]] += 1.0;
            gene_in_use[unit_cand[u].idx[0]] = 1;
            n_unique_units++;
        } else {
            for (int i = 0; i < unit_cand[u].n; i++) gene_in_use[unit_cand[u].idx[i]] = 1;
            multi_units[n_multi++] = u;
            n_multi_units++;
        }
    }

    double *theta = xmalloc(sizeof(double) * (size_t)n_genes);
    double *new_theta = xmalloc(sizeof(double) * (size_t)n_genes);

    /* warm-start initialization from unique counts (with a small pseudocount
     * so genes with zero unique reads, but present in some ambiguous set,
     * still start with nonzero mass and can compete for reassignment) */
    double total_unique = 0.0;
    int n_in_use = 0;
    for (int g = 0; g < n_genes; g++) { total_unique += unique_count[g]; if (gene_in_use[g]) n_in_use++; }

    for (int g = 0; g < n_genes; g++) {
        if (!gene_in_use[g]) { theta[g] = 0.0; continue; }
        theta[g] = (total_unique > 0.0)
            ? (unique_count[g] + 1e-3) / (total_unique + n_in_use * 1e-3)
            : (1.0 / (n_in_use > 0 ? n_in_use : 1));
    }

    int iter;
    double max_delta = 0.0;
    for (iter = 0; iter < MAX_EM_ITERS; iter++) {
        for (int g = 0; g < n_genes; g++) new_theta[g] = unique_count[g];

        for (int m = 0; m < n_multi; m++) {
            GeneSet *cs = &unit_cand[multi_units[m]];
            double denom = 0.0;
            for (int i = 0; i < cs->n; i++) denom += theta[cs->idx[i]];
            for (int i = 0; i < cs->n; i++) {
                double resp = (denom > 0.0) ? theta[cs->idx[i]] / denom : 1.0 / cs->n;
                new_theta[cs->idx[i]] += resp;
            }
        }

        double total = 0.0;
        for (int g = 0; g < n_genes; g++) total += new_theta[g];
        max_delta = 0.0;
        for (int g = 0; g < n_genes; g++) {
            double v = (total > 0.0) ? new_theta[g] / total : 0.0;
            double d = v - theta[g];
            if (d < 0) d = -d;
            if (d > max_delta) max_delta = d;
            theta[g] = v;
        }
        if (max_delta < EM_CONVERGE_EPS) { iter++; break; }
    }
    em_iterations_run = iter;
    em_final_delta = max_delta;

    /* final expected (effective) counts = unique reads + EM-resolved share
     * of each ambiguous read, evaluated at the converged theta */
    for (int g = 0; g < n_genes; g++) genes[g].effective_count = unique_count[g];
    for (int m = 0; m < n_multi; m++) {
        GeneSet *cs = &unit_cand[multi_units[m]];
        double denom = 0.0;
        for (int i = 0; i < cs->n; i++) denom += theta[cs->idx[i]];
        for (int i = 0; i < cs->n; i++) {
            double resp = (denom > 0.0) ? theta[cs->idx[i]] / denom : 1.0 / cs->n;
            genes[cs->idx[i]].effective_count += resp;
        }
    }
    for (int g = 0; g < n_genes; g++) genes[g].unique_count = (long)llround(unique_count[g]);

    free(unit_cand); free(unit_mapped); free(unique_count); free(gene_in_use);
    free(multi_units); free(theta); free(new_theta);
}

/* ---------------------------------------------------------------------- */
/* Output writers                                                         */
/* ---------------------------------------------------------------------- */

static void format_cigar(const Hit *h, char *buf, size_t bufsz) {
    size_t off = 0;
    for (int k = 0; k < h->n_cigar && off < bufsz - 1; k++) {
        int n = snprintf(buf + off, bufsz - off, "%d%c", h->cigar[k].len, h->cigar[k].op);
        if (n < 0) break;
        off += (size_t)n;
    }
    if (off == 0) snprintf(buf, bufsz, "*");
}

static void format_gene_tag(const GeneSet *gs, char *buf, size_t bufsz) {
    if (gs->n == 0) { snprintf(buf, bufsz, "*"); return; }
    size_t off = 0;
    for (int i = 0; i < gs->n && off < bufsz - 1; i++) {
        int n = snprintf(buf + off, bufsz - off, "%s%s", i ? "," : "", genes[gs->idx[i]].gene_id);
        if (n < 0) break;
        off += (size_t)n;
    }
}

static int is_proper_pair(const Read *r1, const Read *r2) {
    if (!r1->aln.mapped || !r2->aln.mapped) return 0;
    const Hit *h1 = &r1->aln.hits[0], *h2 = &r2->aln.hits[0];
    if (h1->chrom_idx != h2->chrom_idx) return 0;
    return (h1->strand == '+' && h2->strand == '-' && hit_span_start0(h1) <= hit_span_start0(h2)) ||
           (h1->strand == '-' && h2->strand == '+' && hit_span_start0(h2) <= hit_span_start0(h1));
}

/* Single-end: emits ONE SAM record per reported hit (secondary flag + NH
 * tag for multi-mappers), matching how real aligners represent them. SEQ/
 * QUAL are truncated to the trimmed length -- the adapter/quality-trimmed
 * tail was never part of what got aligned, so it must not appear in SEQ
 * either, or the CIGAR's op-length sum won't match len(SEQ) (a real SAM
 * validity requirement, and a latent bug in earlier versions of this file
 * that only surfaced once CIGARs could contain soft-clips too). */
static void write_sam_record_se(FILE *f, const Read *r) {
    if (r->trimmed_len <= 0) return;
    char gene_tag[256];
    format_gene_tag(&r->candidates, gene_tag, sizeof(gene_tag));

    if (!r->aln.mapped) {
        fprintf(f, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%.*s\t%.*s\n", r->id, r->trimmed_len, r->seq, r->trimmed_len, r->qual);
        return;
    }
    for (int i = 0; i < r->aln.n_hits; i++) {
        const Hit *h = &r->aln.hits[i];
        int flag = (h->strand == '-') ? 0x10 : 0;
        if (i > 0) flag |= 0x100; /* secondary alignment */
        char cigar[256];
        format_cigar(h, cigar, sizeof(cigar));
        fprintf(f, "%s\t%d\t%s\t%ld\t%d\t%s\t*\t0\t0\t%.*s\t%.*s\tNM:i:%d\tNH:i:%d\tXG:Z:%s%s\n",
                r->id, flag, chroms[h->chrom_idx].name, h->ref_start + 1,
                r->aln.n_hits == 1 ? 60 : 0, cigar, r->trimmed_len, r->seq, r->trimmed_len, r->qual,
                h->mismatches, r->aln.n_hits, gene_tag,
                h->spliced ? (h->canonical_splice ? "\tXS:Z:canonical" : "\tXS:Z:noncanonical") : "");
    }
}

/* Paired-end: single best hit per mate (see README limitation). */
static void write_sam_record_pe(FILE *f, const Read *r, const Read *mate) {
    if (r->trimmed_len <= 0) return;
    char gene_tag[256];
    format_gene_tag(&r->candidates, gene_tag, sizeof(gene_tag));

    int flag = 0x1;
    if (r->read_num == 1) flag |= 0x40; else flag |= 0x80;
    if (is_proper_pair(r, mate)) flag |= 0x2;
    if (!mate->aln.mapped) flag |= 0x8;
    else if (mate->aln.hits[0].strand == '-') flag |= 0x20;

    if (!r->aln.mapped) {
        flag |= 0x4;
        const char *rnext = mate->aln.mapped ? chroms[mate->aln.hits[0].chrom_idx].name : "*";
        long pnext = mate->aln.mapped ? hit_span_start0(&mate->aln.hits[0]) + 1 : 0;
        fprintf(f, "%s\t%d\t*\t0\t0\t*\t%s\t%ld\t0\t%.*s\t%.*s\n", r->id, flag, rnext, pnext, r->trimmed_len, r->seq, r->trimmed_len, r->qual);
        return;
    }
    const Hit *h = &r->aln.hits[0];
    if (h->strand == '-') flag |= 0x10;
    char cigar[256];
    format_cigar(h, cigar, sizeof(cigar));

    const char *rnext = "*"; long pnext = 0, tlen = 0;
    if (mate->aln.mapped) {
        const Hit *hm = &mate->aln.hits[0];
        rnext = (hm->chrom_idx == h->chrom_idx) ? "=" : chroms[hm->chrom_idx].name;
        pnext = hit_span_start0(hm) + 1;
        if (hm->chrom_idx == h->chrom_idx) {
            long left = lmin(hit_span_start0(h), hit_span_start0(hm));
            long right = lmax(hit_span_end0(h), hit_span_end0(hm));
            long span = right - left + 1;
            tlen = (h->ref_start + 1 <= pnext) ? span : -span;
        }
    }
    fprintf(f, "%s\t%d\t%s\t%ld\t60\t%s\t%s\t%ld\t%ld\t%.*s\t%.*s\tNM:i:%d\tNH:i:1\tXG:Z:%s%s\n",
            r->id, flag, chroms[h->chrom_idx].name, h->ref_start + 1, cigar,
            rnext, pnext, tlen, r->trimmed_len, r->seq, r->trimmed_len, r->qual, h->mismatches, gene_tag,
            h->spliced ? (h->canonical_splice ? "\tXS:Z:canonical" : "\tXS:Z:noncanonical") : "");
}

static void write_sam(const char *outdir) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/alignments.sam", outdir);
    FILE *f = fopen(path, "w");
    if (!f) die("cannot write alignments.sam");

    fprintf(f, "@HD\tVN:1.6\tSO:unsorted\n");
    for (int i = 0; i < n_chroms; i++) fprintf(f, "@SQ\tSN:%s\tLN:%ld\n", chroms[i].name, chroms[i].len);
    fprintf(f, "@PG\tID:rnaseq_pipeline.c\tPN:rnaseq_pipeline\tVN:4.0\n");

    for (int i = 0; i < n_reads; i++) {
        Read *r = &reads[i];
        if (r->read_num == 0) write_sam_record_se(f, r);
        else write_sam_record_pe(f, r, &reads[r->mate_idx]);
    }
    fclose(f);
}

static void write_gene_counts(const char *outdir) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/gene_counts.tsv", outdir);
    FILE *f = fopen(path, "w");
    if (!f) die("cannot write gene_counts.tsv");

    fprintf(f, "gene_id\tchrom\tstart\tend\tstrand\tunique_count\teffective_count\n");
    for (int i = 0; i < n_genes; i++)
        fprintf(f, "%s\t%s\t%ld\t%ld\t%c\t%ld\t%.3f\n",
                genes[i].gene_id, genes[i].chrom, genes[i].start, genes[i].end, genes[i].strand,
                genes[i].unique_count, genes[i].effective_count);
    fprintf(f, "__no_feature\t*\t0\t0\t*\t%ld\t%.3f\n", n_no_feature_units, (double)n_no_feature_units);
    fprintf(f, "__unmapped\t*\t0\t0\t*\t%ld\t%.3f\n", n_unmapped_units, (double)n_unmapped_units);
    fclose(f);
}

static void write_qc_report(const char *outdir) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/qc_report.txt", outdir);
    FILE *f = fopen(path, "w");
    if (!f) die("cannot write qc_report.txt");

    double sum_gc = 0, sum_q = 0;
    long sum_n = 0;
    for (int i = 0; i < n_reads; i++) { sum_gc += reads[i].gc_pct; sum_q += reads[i].mean_q_raw; sum_n += reads[i].n_count; }

    fprintf(f, "==========================================================\n");
    fprintf(f, " Per-read QC report (raw, pre-trim)  -- FastQC-style\n");
    fprintf(f, "==========================================================\n");
    fprintf(f, "Total reads              : %d%s\n", n_reads, paired_mode ? " (R1+R2)" : "");
    fprintf(f, "Total raw bases           : %ld\n", total_raw_bases);
    fprintf(f, "Mean GC content (%%)       : %.2f\n", n_reads ? sum_gc / n_reads : 0.0);
    fprintf(f, "Mean per-base Phred score : %.2f\n", n_reads ? sum_q / n_reads : 0.0);
    fprintf(f, "Total ambiguous (N) bases : %ld\n\n", sum_n);
    fprintf(f, "read_id\tmate\traw_len\tgc_pct\tmean_q\tn_count\tn_hits\n");
    for (int i = 0; i < n_reads; i++) {
        Read *r = &reads[i];
        fprintf(f, "%s\t%d\t%d\t%.1f\t%.1f\t%d\t%d\n", r->id, r->read_num, r->raw_len,
                r->gc_pct, r->mean_q_raw, r->n_count, r->aln.mapped ? r->aln.n_hits : 0);
    }
    fclose(f);
}

static void write_summary(const char *outdir, const char *ref_path, const char *gtf_path,
                           const char *reads_desc) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/multiqc_summary.txt", outdir);
    FILE *f = fopen(path, "w");
    if (!f) die("cannot write multiqc_summary.txt");

    int step = paired_mode ? 2 : 1;
    long n_units = n_reads / step;
    long aligned_units = n_unique_units + n_multi_units + n_no_feature_units;
    double total_effective = 0.0;
    for (int i = 0; i < n_genes; i++) total_effective += genes[i].effective_count;

    fprintf(f, "============================================================\n");
    fprintf(f, "   rnaseq_pipeline.c v4 -- simplified nf-core/rnaseq analog\n");
    fprintf(f, "   MultiQC-style run summary\n");
    fprintf(f, "============================================================\n\n");
    fprintf(f, "Mode: %s\n", paired_mode ? "paired-end" : "single-end");
    fprintf(f, "Inputs\n");
    fprintf(f, "  Reference FASTA : %s (%d sequence(s))\n", ref_path, n_chroms);
    fprintf(f, "  Annotation GTF  : %s (%d gene(s))\n", gtf_path, n_genes);
    fprintf(f, "  Reads           : %s (%ld %s)\n\n", reads_desc, n_units, paired_mode ? "read pair(s)" : "read(s)");

    fprintf(f, "QC + Trimming\n");
    fprintf(f, "  Raw bases                 : %ld\n", total_raw_bases);
    fprintf(f, "  Trimmed bases (kept)      : %ld\n", total_trimmed_bases);
    fprintf(f, "  Reads with adapter found  : %ld (%.1f%%)\n", reads_with_adapter,
            n_reads ? 100.0 * reads_with_adapter / n_reads : 0.0);
    fprintf(f, "  Reads quality-trimmed     : %ld (%.1f%%)\n", reads_quality_trimmed,
            n_reads ? 100.0 * reads_quality_trimmed / n_reads : 0.0);
    fprintf(f, "  Reads dropped (< %d bp)   : %ld\n\n", MIN_READ_LEN_AFTER_TRIM, reads_dropped_too_short);

    fprintf(f, "Alignment\n");
    fprintf(f, "  %s aligned            : %ld / %ld (%.1f%%)\n", paired_mode ? "Fragments" : "Reads",
            aligned_units, n_units, n_units ? 100.0 * aligned_units / n_units : 0.0);
    fprintf(f, "  Ungapped alignments        : %ld\n", n_ungapped_alignments);
    fprintf(f, "  Spliced alignments         : %ld\n", n_spliced_alignments);
    if (paired_mode)
        fprintf(f, "  Proper pairs (FR, same chrom): %ld (%.1f%% of pairs)\n", n_proper_pairs, n_units ? 100.0 * n_proper_pairs / n_units : 0.0);
    fprintf(f, "\n");

    fprintf(f, "Gene-level quantification (EM, RSEM/Salmon-style)\n");
    fprintf(f, "  %s uniquely assigned  : %ld\n", paired_mode ? "Fragments" : "Reads", n_unique_units);
    fprintf(f, "  %s multi-mapped (EM)  : %ld\n", paired_mode ? "Fragments" : "Reads", n_multi_units);
    fprintf(f, "  %s mapped, no gene    : %ld\n", paired_mode ? "Fragments" : "Reads", n_no_feature_units);
    fprintf(f, "  %s unmapped            : %ld\n", paired_mode ? "Fragments" : "Reads", n_unmapped_units);
    fprintf(f, "  EM iterations to converge : %d (final max |delta theta| = %.2e)\n", em_iterations_run, em_final_delta);
    fprintf(f, "  Sum of effective counts   : %.3f (sanity check vs %ld unique+multi units)\n\n",
            total_effective, n_unique_units + n_multi_units);
    fprintf(f, "  %-16s %12s %16s\n", "gene_id", "unique_count", "effective_count");
    for (int i = 0; i < n_genes; i++)
        fprintf(f, "  %-16s %12ld %16.3f\n", genes[i].gene_id, genes[i].unique_count, genes[i].effective_count);

    fprintf(f, "\nOutput files\n");
    fprintf(f, "  alignments.sam       - SAM alignments (multi-hit records + NH tag in SE mode)\n");
    fprintf(f, "  gene_counts.tsv      - unique_count (hard) + effective_count (EM) per gene\n");
    fprintf(f, "  qc_report.txt        - per-read QC metrics (FastQC-style)\n");
    fprintf(f, "  multiqc_summary.txt  - this file\n");
    fclose(f);
}

/* ---------------------------------------------------------------------- */
/* main                                                                    */
/* ---------------------------------------------------------------------- */

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s <reference.fasta> <annotation.gtf> se <reads.fastq> <outdir>\n"
        "  %s <reference.fasta> <annotation.gtf> pe <r1.fastq> <r2.fastq> <outdir>\n",
        prog, prog);
}

int main(int argc, char **argv) {
    if (argc < 6) { usage(argv[0]); return 1; }

    const char *ref_path = argv[1];
    const char *gtf_path = argv[2];
    const char *mode = argv[3];
    char reads_desc[1024];
    const char *outdir;

    if (strcmp(mode, "se") == 0) {
        if (argc != 6) { usage(argv[0]); return 1; }
        paired_mode = 0; outdir = argv[5];
        snprintf(reads_desc, sizeof(reads_desc), "%s", argv[4]);
    } else if (strcmp(mode, "pe") == 0) {
        if (argc != 7) { usage(argv[0]); return 1; }
        paired_mode = 1; outdir = argv[6];
        snprintf(reads_desc, sizeof(reads_desc), "%s + %s", argv[4], argv[5]);
    } else { usage(argv[0]); return 1; }

    printf("[1/7] Loading reference + k-mer index (k=%d): %s\n", KMER_LEN, ref_path);
    if (try_load_index_cache(ref_path)) {
        printf("      -> loaded from cache (%s.kidx): %d sequence(s), skipped FASTA parse + index build\n", ref_path, n_chroms);
        printf("[2/7] (skipped -- index came from cache)\n");
    } else {
        load_reference(ref_path);
        printf("      -> %d sequence(s) loaded from FASTA\n", n_chroms);
        printf("[2/7] Building k-mer index (k=%d)\n", KMER_LEN);
        build_kmer_index();
        save_index_cache(ref_path);
        printf("      -> index cached for future runs: %s.kidx\n", ref_path);
    }

    printf("[3/7] Loading GTF annotation: %s\n", gtf_path);
    load_gtf(gtf_path);
    printf("      -> %d gene(s) loaded\n", n_genes);

    printf("[4/7] Loading reads (%s mode): %s\n", paired_mode ? "paired-end" : "single-end", reads_desc);
    if (paired_mode) load_fastq_pe(argv[4], argv[5]); else load_fastq_se(argv[4]);
    printf("      -> %d read(s) loaded\n", n_reads);

    printf("[5/7] Running QC + adapter/quality trimming\n");
    #pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < n_reads; i++) { compute_raw_qc(&reads[i]); trim_read(&reads[i]); }

    printf("[6/7] Aligning reads (k-mer seed + ungapped/spliced extension, multi-mapping-aware)\n");
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n_reads; i++) align_read(&reads[i]);
    if (paired_mode) {
        long proper_local = 0;
        #pragma omp parallel for schedule(dynamic, 256) reduction(+:proper_local)
        for (int i = 0; i < n_reads; i += 2)
            if (is_proper_pair(&reads[i], &reads[i+1])) proper_local++;
        n_proper_pairs = proper_local;
    }

    printf("[7/7] EM quantification + writing output files to: %s\n", outdir);
    quantify_em();
    write_sam(outdir);
    write_gene_counts(outdir);
    write_qc_report(outdir);
    write_summary(outdir, ref_path, gtf_path, reads_desc);

    printf("Done. (EM converged in %d iterations)\n", em_iterations_run);
    return 0;
}

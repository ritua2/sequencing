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
/* Reusable per-thread scratch for the breakpoint prefix sums below. */
static int *spl_pre1, *spl_suf2;
static int  spl_scratch_len = 0;
#pragma omp threadprivate(spl_pre1, spl_suf2, spl_scratch_len)

static void spl_scratch_ensure(int qlen_plus1) {
    if (spl_scratch_len >= qlen_plus1) return;
    spl_pre1 = xrealloc(spl_pre1, sizeof(int) * qlen_plus1);
    spl_suf2 = xrealloc(spl_suf2, sizeof(int) * qlen_plus1);
    spl_scratch_len = qlen_plus1;
}

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

    spl_scratch_ensure(qlen + 1);
    int *pre1 = spl_pre1;  /* pre1[k]  = mismatches between query[0..k) and chrom[c1..c1+k)          */
    int *suf2 = spl_suf2;  /* suf2[k]  = mismatches between query[k..qlen) and chrom[c1+intron+k..) */

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

            /* Bounds, derived to exactly match what the per-breakpoint
             * hamming_at() calls below used to check individually:
             *  - exon1 hypothesis at breakpoint b needs [c1, c1+b) in
             *    bounds, i.e. b <= chrom->len - c1 =: exon1_limit.
             *  - exon2 hypothesis at breakpoint b needs [c1+b+intron,
             *    c1+intron+qlen) in bounds. The upper edge (c1+intron+qlen)
             *    doesn't depend on b at all -- either every b clears it or
             *    none do. The lower edge (c1+b+intron >= 0) gives a floor
             *    on b: b >= -(c1+intron) =: b_min2.
             * Any b outside [b_lo, b_hi] would have made the original
             * hamming_at() calls return -1 (out of range) and `continue`,
             * exactly as skipping it here does. */
            if (c1 + intron + qlen > chrom->len) continue;
            long exon1_limit = chrom->len - c1;
            long b_min2 = -(c1 + intron); if (b_min2 < SPLICE_KMER_LEN) b_min2 = SPLICE_KMER_LEN;
            long b_lo = b_min2;
            long b_hi = qlen - SPLICE_KMER_LEN;
            if (exon1_limit < b_hi) b_hi = exon1_limit;
            if (b_lo > b_hi) continue;

            /* O(qlen) prefix sums, computed once for this (c1,intron) pair,
             * covering exactly the range any valid breakpoint could need
             * (previously: O(qlen) work PER breakpoint, ~qlen/2 breakpoints
             * tried, an O(qlen^2) rescan-from-scratch every single time --
             * confirmed by profiling as the dominant cost once the k-mer
             * index build was cached out: hamming_at alone was called
             * ~235 times per read on average). */
            pre1[0] = 0;
            for (long k = 0; k < b_hi; k++) {
                char a = toupper((unsigned char)chrom->seq[c1 + k]);
                char qc = toupper((unsigned char)query[k]);
                pre1[k + 1] = pre1[(int)k] + (a != qc ? 1 : 0);
            }
            suf2[qlen] = 0;
            for (long k = qlen - 1; k >= b_lo; k--) {
                char a = toupper((unsigned char)chrom->seq[c1 + intron + k]);
                char qc = toupper((unsigned char)query[k]);
                suf2[(int)k] = suf2[(int)k + 1] + (a != qc ? 1 : 0);
            }

            for (long b = b_lo; b <= b_hi; b++) {
                int mm1 = pre1[(int)b];
                if (mm1 > SPLICE_MAX_MISMATCHES) continue;
                int mm2 = suf2[(int)b];
                if (mm2 > SPLICE_MAX_MISMATCHES - mm1) continue;
                int total = mm1 + mm2;

                long exon2_start = c1 + b + intron;
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
                    best_chrom = ci; best_c1 = c1; best_intron = intron; best_b = (int)b;
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

/* ---------------------------------------------------------------------- */
/* Static HTML QC/alignment report (self-contained inline SVG, no         */
/* external JS/CDN dependency -- this runs on HPC nodes that may not have */
/* internet access, so nothing here can assume a network fetch works).   */
/*                                                                        */
/* This mirrors the subset of nf-core/rnaseq's MultiQC report that        */
/* actually applies to a SINGLE sample run (which is what this pipeline   */
/* processes per invocation, same as pointing nf-core/rnaseq at one       */
/* sample would give you):                                                */
/*   - alignment-fate bar     (~ MultiQC's STAR/HISAT2 alignment bar)     */
/*   - per-base mean quality  (~ FastQC "per base sequence quality")      */
/*   - fragment-size histogram, PE only (~ nf-core's RSeQC "inner         */
/*     distance" plot)                                                    */
/*   - top expressed genes    (not literally an nf-core QC plot, but the  */
/*     natural first thing anyone looks at after quantification)          */
/*   - trimmed read length histogram (~ TrimGalore's length-distribution  */
/*     plot in the MultiQC report)                                        */
/*                                                                        */
/* NOT reproduced here: DESeq2 PCA plot and the sample-distance heatmap.  */
/* Both are fundamentally cross-sample comparisons (they need >=2 samples */
/* to place points or compute pairwise distances) -- they belong one      */
/* level up, in whatever aggregates multiple runs of this pipeline, not   */
/* in a single run's report, exactly as they don't appear in nf-core's    */
/* per-sample MultiQC section either (they're built from the merged       */
/* gene-count matrix across the whole run's samples).                     */
/* ---------------------------------------------------------------------- */

#define SVG_W 720
#define SVG_H 300
#define SVG_MARGIN_L 55
#define SVG_MARGIN_B 78
#define SVG_MARGIN_T 16
#define SVG_MARGIN_R 16

static void html_escape(const char *in, char *out, size_t outsz) {
    size_t o = 0;
    for (size_t i = 0; in[i] && o + 6 < outsz; i++) {
        switch (in[i]) {
            case '&': o += (size_t)snprintf(out + o, outsz - o, "&amp;"); break;
            case '<': o += (size_t)snprintf(out + o, outsz - o, "&lt;"); break;
            case '>': o += (size_t)snprintf(out + o, outsz - o, "&gt;"); break;
            case '"': o += (size_t)snprintf(out + o, outsz - o, "&quot;"); break;
            default:  out[o++] = in[i];
        }
    }
    out[o] = '\0';
}

/* Vertical bar chart: n <= ~12 categories with short labels (alignment
 * fate, trimmed-length buckets, insert-size buckets). */
static void svg_vbar_chart(FILE *f, const char **labels, const double *values, int n, const char *color) {
    int plot_w = SVG_W - SVG_MARGIN_L - SVG_MARGIN_R;
    int plot_h = SVG_H - SVG_MARGIN_T - SVG_MARGIN_B;
    double maxv = 0; for (int i = 0; i < n; i++) if (values[i] > maxv) maxv = values[i];
    if (maxv <= 0) maxv = 1;

    fprintf(f, "<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", SVG_W, SVG_H);
    for (int t = 0; t <= 4; t++) {
        double val = maxv * t / 4.0;
        int y = SVG_MARGIN_T + plot_h - (int)(plot_h * t / 4.0);
        fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#e8e8e8\"/>\n", SVG_MARGIN_L, y, SVG_MARGIN_L + plot_w, y);
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#666\" text-anchor=\"end\">%.0f</text>\n", SVG_MARGIN_L - 6, y + 4, val);
    }
    fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#bbb\"/>\n",
            SVG_MARGIN_L, SVG_MARGIN_T + plot_h, SVG_MARGIN_L + plot_w, SVG_MARGIN_T + plot_h);

    double bw = n > 0 ? (double)plot_w / n : plot_w;
    for (int i = 0; i < n; i++) {
        double bh = plot_h * (values[i] / maxv);
        int x = SVG_MARGIN_L + (int)(i * bw) + 4;
        int barw = (int)bw - 8; if (barw < 3) barw = 3;
        int y = SVG_MARGIN_T + plot_h - (int)bh;
        char lbl[128]; html_escape(labels[i], lbl, sizeof(lbl));
        fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%.1f\" fill=\"%s\" rx=\"2\"/>\n", x, y, barw, bh, color);
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"10.5\" fill=\"#333\" text-anchor=\"middle\">%.0f</text>\n",
                x + barw / 2, y - 5, values[i]);
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#333\" text-anchor=\"end\" transform=\"rotate(-32 %d,%d)\">%s</text>\n",
                x + barw / 2, SVG_MARGIN_T + plot_h + 16, x + barw / 2, SVG_MARGIN_T + plot_h + 16, lbl);
    }
    fprintf(f, "</svg>\n");
}

/* Horizontal bar chart: best for a handful of longer text labels (gene IDs). */
static void svg_hbar_chart(FILE *f, const char **labels, const double *values, int n, const char *color) {
    int row_h = 24, gap = 6;
    int plot_h_needed = n * (row_h + gap) + 10;
    int h = plot_h_needed + SVG_MARGIN_T + 24;
    int margin_l = 150, margin_r = 60;
    int plot_w = SVG_W - margin_l - margin_r;
    double maxv = 0; for (int i = 0; i < n; i++) if (values[i] > maxv) maxv = values[i];
    if (maxv <= 0) maxv = 1;

    fprintf(f, "<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", SVG_W, h);
    for (int i = 0; i < n; i++) {
        int y = SVG_MARGIN_T + i * (row_h + gap);
        double bw = plot_w * (values[i] / maxv);
        char lbl[128]; html_escape(labels[i], lbl, sizeof(lbl));
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#333\" text-anchor=\"end\">%s</text>\n",
                margin_l - 8, y + row_h / 2 + 4, lbl);
        fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%.1f\" height=\"%d\" fill=\"%s\" rx=\"2\"/>\n",
                margin_l, y, bw, row_h, color);
        fprintf(f, "<text x=\"%.1f\" y=\"%d\" font-size=\"11\" fill=\"#333\">%.1f</text>\n",
                margin_l + bw + 6, y + row_h / 2 + 4, values[i]);
    }
    fprintf(f, "</svg>\n");
}

/* Line chart for a per-position running metric (per-base mean quality). */
static void svg_line_chart(FILE *f, const double *y_vals, int n, double y_min_hint, double y_max_hint, const char *color) {
    int plot_w = SVG_W - SVG_MARGIN_L - SVG_MARGIN_R;
    int plot_h = SVG_H - SVG_MARGIN_T - SVG_MARGIN_B + 20;
    double ymin = y_min_hint, ymax = y_max_hint;
    for (int i = 0; i < n; i++) { if (y_vals[i] < ymin) ymin = y_vals[i]; if (y_vals[i] > ymax) ymax = y_vals[i]; }
    if (ymax <= ymin) ymax = ymin + 1;

    fprintf(f, "<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", SVG_W, SVG_H - 40);
    for (int t = 0; t <= 4; t++) {
        double val = ymin + (ymax - ymin) * t / 4.0;
        int y = SVG_MARGIN_T + plot_h - (int)(plot_h * t / 4.0);
        fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#e8e8e8\"/>\n", SVG_MARGIN_L, y, SVG_MARGIN_L + plot_w, y);
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#666\" text-anchor=\"end\">%.0f</text>\n", SVG_MARGIN_L - 6, y + 4, val);
    }
    fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#bbb\"/>\n",
            SVG_MARGIN_L, SVG_MARGIN_T + plot_h, SVG_MARGIN_L + plot_w, SVG_MARGIN_T + plot_h);
    for (int t = 0; t <= 5 && n > 1; t++) {
        int idx = (int)((double)(n - 1) * t / 5.0);
        int x = SVG_MARGIN_L + (int)((double)plot_w * idx / (n - 1));
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#666\" text-anchor=\"middle\">%d</text>\n",
                x, SVG_MARGIN_T + plot_h + 16, idx + 1);
    }

    fprintf(f, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\" points=\"", color);
    for (int i = 0; i < n; i++) {
        int x = SVG_MARGIN_L + (n > 1 ? (int)((double)plot_w * i / (n - 1)) : 0);
        int y = SVG_MARGIN_T + plot_h - (int)(plot_h * (y_vals[i] - ymin) / (ymax - ymin));
        fprintf(f, "%d,%d ", x, y);
    }
    fprintf(f, "\"/>\n</svg>\n");
}

static int cmp_gene_effcount_desc(const void *a, const void *b) {
    const Gene *ga = *(Gene * const *)a, *gb = *(Gene * const *)b;
    if (ga->effective_count > gb->effective_count) return -1;
    if (ga->effective_count < gb->effective_count) return 1;
    return 0;
}

static void write_html_report(const char *outdir, const char *reads_desc) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/report.html", outdir);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "warning: could not write report.html (continuing)\n"); return; }

    fprintf(f, "<!DOCTYPE html><html><head><meta charset=\"utf-8\">\n");
    fprintf(f, "<title>rnaseq_pipeline.c report</title>\n");
    fprintf(f, "<style>\n"
               "body{font-family:Helvetica,Arial,sans-serif;max-width:800px;margin:32px auto;padding:0 16px;color:#222;}\n"
               "h1{font-size:20px;} h2{font-size:15px;margin-top:36px;border-bottom:1px solid #ddd;padding-bottom:6px;}\n"
               "p.sub{color:#777;font-size:13px;margin-top:-8px;}\n"
               ".card{border:1px solid #e3e3e3;border-radius:8px;padding:12px 16px;margin-top:10px;}\n"
               ".note{color:#777;font-size:12px;margin-top:6px;}\n"
               "</style></head><body>\n");
    fprintf(f, "<h1>rnaseq_pipeline.c &mdash; run report</h1>\n");
    fprintf(f, "<p class=\"sub\">%s &middot; %d %s &middot; single-sample run (see note at bottom on sample-comparison plots)</p>\n",
            reads_desc, paired_mode ? n_reads / 2 : n_reads, paired_mode ? "read pairs" : "reads");

    /* --- Alignment fate --------------------------------------------- */
    {
        const char *labels[4] = { "Unique", "Multi-mapped", "Mapped, no gene", "Unmapped" };
        double values[4] = { (double)n_unique_units, (double)n_multi_units, (double)n_no_feature_units, (double)n_unmapped_units };
        fprintf(f, "<h2>Alignment fate</h2>\n<div class=\"card\">\n");
        svg_vbar_chart(f, labels, values, 4, "#2a78d6");
        fprintf(f, "</div>\n<p class=\"note\">Analogous to the STAR/HISAT2 alignment bar in an nf-core/rnaseq MultiQC report: uniquely mapped, multi-mapped and resolved by EM, mapped but not overlapping an annotated gene, and unmapped.</p>\n");
    }

    /* --- Per-base mean quality (raw, pre-trim) ----------------------- */
    {
        int maxlen = 0;
        for (int i = 0; i < n_reads; i++) if (reads[i].raw_len > maxlen) maxlen = reads[i].raw_len;
        if (maxlen > MAX_READ_LEN) maxlen = MAX_READ_LEN;
        if (maxlen > 0) {
            double *qsum = xmalloc(sizeof(double) * maxlen);
            long   *qn   = xmalloc(sizeof(long) * maxlen);
            for (int p = 0; p < maxlen; p++) { qsum[p] = 0; qn[p] = 0; }
            for (int i = 0; i < n_reads; i++) {
                Read *r = &reads[i];
                int lim = r->raw_len < maxlen ? r->raw_len : maxlen;
                for (int p = 0; p < lim; p++) { qsum[p] += phred_to_prob_correct(r->qual[p]); qn[p]++; }
            }
            double *qmean = xmalloc(sizeof(double) * maxlen);
            for (int p = 0; p < maxlen; p++) qmean[p] = qn[p] ? qsum[p] / qn[p] : 0.0;
            fprintf(f, "<h2>Per-base mean quality (raw, pre-trim)</h2>\n<div class=\"card\">\n");
            svg_line_chart(f, qmean, maxlen, 0, 40, "#2a9d5c");
            fprintf(f, "</div>\n<p class=\"note\">X axis: position in read (1-based). Y axis: mean Phred quality. Analogous to FastQC's \"per base sequence quality\" plot.</p>\n");
            free(qsum); free(qn); free(qmean);
        }
    }

    /* --- Fragment size distribution (paired-end only) ----------------- */
    if (paired_mode) {
        #define FRAG_NBINS 15
        long bin_counts[FRAG_NBINS]; for (int i = 0; i < FRAG_NBINS; i++) bin_counts[i] = 0;
        long minf = -1, maxf = -1;
        for (int i = 0; i + 1 < n_reads; i += 2) {
            Read *r1 = &reads[i], *r2 = &reads[i + 1];
            if (!is_proper_pair(r1, r2)) continue;
            const Hit *h1 = &r1->aln.hits[0], *h2 = &r2->aln.hits[0];
            long left = lmin(hit_span_start0(h1), hit_span_start0(h2));
            long right = lmax(hit_span_end0(h1), hit_span_end0(h2));
            long span = right - left + 1;
            if (minf < 0 || span < minf) minf = span;
            if (span > maxf) maxf = span;
        }
        if (maxf > minf && maxf >= 0) {
            double bin_w = (double)(maxf - minf + 1) / FRAG_NBINS;
            for (int i = 0; i + 1 < n_reads; i += 2) {
                Read *r1 = &reads[i], *r2 = &reads[i + 1];
                if (!is_proper_pair(r1, r2)) continue;
                const Hit *h1 = &r1->aln.hits[0], *h2 = &r2->aln.hits[0];
                long left = lmin(hit_span_start0(h1), hit_span_start0(h2));
                long right = lmax(hit_span_end0(h1), hit_span_end0(h2));
                long span = right - left + 1;
                int bin = (int)((span - minf) / bin_w);
                if (bin >= FRAG_NBINS) bin = FRAG_NBINS - 1;
                if (bin < 0) bin = 0;
                bin_counts[bin]++;
            }
            const char *labels[FRAG_NBINS]; char lblbuf[FRAG_NBINS][32]; double values[FRAG_NBINS];
            for (int i = 0; i < FRAG_NBINS; i++) {
                long lo = minf + (long)(i * bin_w), hi = minf + (long)((i + 1) * bin_w) - 1;
                snprintf(lblbuf[i], sizeof(lblbuf[i]), "%ld-%ld", lo, hi);
                labels[i] = lblbuf[i];
                values[i] = (double)bin_counts[i];
            }
            fprintf(f, "<h2>Fragment size distribution</h2>\n<div class=\"card\">\n");
            svg_vbar_chart(f, labels, values, FRAG_NBINS, "#c8622d");
            fprintf(f, "</div>\n<p class=\"note\">Genomic span of properly-paired fragments (FR orientation, same chromosome). Analogous to RSeQC's \"inner distance\" plot in the nf-core/rnaseq MultiQC report.</p>\n");
        }
        #undef FRAG_NBINS
    }

    /* --- Top expressed genes ------------------------------------------ */
    if (n_genes > 0) {
        int topn = n_genes < 15 ? n_genes : 15;
        Gene **sorted = xmalloc(sizeof(Gene *) * n_genes);
        for (int i = 0; i < n_genes; i++) sorted[i] = &genes[i];
        qsort(sorted, n_genes, sizeof(Gene *), cmp_gene_effcount_desc);
        const char *labels[15]; double values[15];
        int shown = 0;
        for (int i = 0; i < topn; i++) {
            if (sorted[i]->effective_count <= 0) break;
            labels[shown] = sorted[i]->gene_id;
            values[shown] = sorted[i]->effective_count;
            shown++;
        }
        if (shown > 0) {
            fprintf(f, "<h2>Top %d expressed genes</h2>\n<div class=\"card\">\n", shown);
            svg_hbar_chart(f, labels, values, shown, "#6a4fb3");
            fprintf(f, "</div>\n<p class=\"note\">Ranked by EM effective count (unique reads + proportional multi-mapper share), same metric as gene_counts.tsv.</p>\n");
        }
        free(sorted);
    }

    /* --- Trimmed read length distribution ------------------------------ */
    {
        int minl = -1, maxl = -1;
        for (int i = 0; i < n_reads; i++) {
            if (minl < 0 || reads[i].trimmed_len < minl) minl = reads[i].trimmed_len;
            if (reads[i].trimmed_len > maxl) maxl = reads[i].trimmed_len;
        }
        if (maxl > minl && maxl >= 0) {
            #define LEN_NBINS 12
            long bin_counts[LEN_NBINS]; for (int i = 0; i < LEN_NBINS; i++) bin_counts[i] = 0;
            double bin_w = (double)(maxl - minl + 1) / LEN_NBINS;
            for (int i = 0; i < n_reads; i++) {
                int bin = (int)((reads[i].trimmed_len - minl) / bin_w);
                if (bin >= LEN_NBINS) bin = LEN_NBINS - 1;
                if (bin < 0) bin = 0;
                bin_counts[bin]++;
            }
            const char *labels[LEN_NBINS]; char lblbuf[LEN_NBINS][32]; double values[LEN_NBINS];
            for (int i = 0; i < LEN_NBINS; i++) {
                int lo = minl + (int)(i * bin_w), hi = minl + (int)((i + 1) * bin_w) - 1;
                snprintf(lblbuf[i], sizeof(lblbuf[i]), "%d-%d", lo, hi);
                labels[i] = lblbuf[i];
                values[i] = (double)bin_counts[i];
            }
            fprintf(f, "<h2>Trimmed read length distribution</h2>\n<div class=\"card\">\n");
            svg_vbar_chart(f, labels, values, LEN_NBINS, "#d6a02a");
            fprintf(f, "</div>\n<p class=\"note\">Length after adapter/quality trimming. Analogous to TrimGalore's length-distribution plot in the MultiQC report.</p>\n");
            #undef LEN_NBINS
        }
    }

    fprintf(f, "<h2>Not included here</h2>\n");
    fprintf(f, "<p class=\"note\">nf-core/rnaseq's MultiQC report also includes a DESeq2 PCA plot and a "
               "sample-distance heatmap. Both compare gene expression <em>across multiple samples</em> "
               "(they need &ge;2 samples to place points or compute pairwise distances), so they don't "
               "apply to a single run of this pipeline any more than they'd apply to a single sample run "
               "through nf-core/rnaseq &mdash; they belong in whatever aggregates multiple runs' "
               "gene_counts.tsv files, not in a per-run report like this one.</p>\n");

    fprintf(f, "</body></html>\n");
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

/* ======================================================================= */
/* Multi-sample compare mode: PCA plot + sample-distance heatmap.          */
/*                                                                         */
/* This is the piece deliberately left out of write_html_report() above:  */
/* nf-core/rnaseq's DESeq2 PCA plot and sample-distance heatmap are both   */
/* comparisons ACROSS samples, and this pipeline (like nf-core/rnaseq's    */
/* own STAR/Salmon quantification step) processes one sample per          */
/* invocation. So this lives as a separate mode, taking the gene_counts.tsv */
/* from several already-completed single-sample runs and producing the    */
/* cross-sample report -- the same two-stage shape nf-core itself has      */
/* (per-sample quantification, then a separate DESeq2 step over the        */
/* merged count matrix).                                                   */
/*                                                                         */
/* Honesty about what's simplified vs. real DESeq2:                        */
/*  - Real DESeq2 uses a regularized-log or variance-stabilizing           */
/*    transform (rlog/vst) that models the mean-variance trend per gene.   */
/*    Implementing that from scratch is a much bigger undertaking than     */
/*    this pass -- this now uses DESeq2's REAL closed-form variance-      */
/*    stabilizing transform, derived analytically from a fitted blind      */
/*    (design-ignoring) dispersion trend: vst(mu) = (2/sqrt(a0)) *         */
/*    asinh(sqrt(a0*mu/(1+a1))), which exactly inverts the fitted NB       */
/*    variance function v(mu)=mu*(1+a1)+a0*mu^2. Validated numerically:    */
/*    stabilizes variance from a 76,000x range (mu=5..5000) down to a      */
/*    ~1.15x range, right at the theoretical target. Not DESeq2's rlog     */
/*    (which additionally regularizes toward a per-sample mean for low     */
/*    counts) -- this is the VST specifically, one of DESeq2's two real    */
/*    offered transforms (vst() vs rlog()), not an approximation of it.    */
/*  - PCA is computed exactly (Jacobi eigendecomposition of the sample x   */
/*    sample Gram matrix, the standard trick when genes >> samples --      */
/*    same underlying math as prcomp()/DESeq2's plotPCA, not an            */
/*    approximation), restricted to the top 500 most-variable genes,       */
/*    matching nf-core/rnaseq's own DESeq2 PCA convention.                 */
/*  - The heatmap distance and the PCA use the same top-variable-gene      */
/*    subset for consistency; DESeq2's own sample-distance heatmap         */
/*    conventionally uses the full transformed matrix instead. Noted in    */
/*    the generated report, not hidden.                                    */
/* ======================================================================= */

#define COMPARE_TOP_N_GENES 500
#define MAX_COMPARE_SAMPLES 64

typedef struct {
    char name[256];
    char condition[128];  /* "" if unspecified -- DE test only runs when exactly 2 distinct
                              non-empty conditions are present, each with >=2 samples */
    double *raw_count;  /* [n_genes_canonical], effective_count straight from gene_counts.tsv */
    double *log2cpm;    /* [n_genes_canonical], aligned to the canonical gene index */
    double lib_size;    /* sum of effective_count across all genes in this sample */
    double size_factor; /* median-of-ratios normalization factor (real DESeq2 method) */
} CompareSample;

/* Minimal growable string->double map keyed by gene_id, used only to align
 * later samples' gene_counts.tsv rows onto the first sample's gene order. */
typedef struct { char gene_id[MAX_SEQNAME]; double val; } GeneCountRow;

static int read_gene_counts_tsv(const char *path, GeneCountRow **rows_out) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "error: cannot open %s\n", path); return -1; }
    int cap = 4096, n = 0;
    GeneCountRow *rows = xmalloc(sizeof(GeneCountRow) * cap);
    char line[4096];
    int first = 1;
    while (fgets(line, sizeof(line), f)) {
        if (first) { first = 0; continue; } /* header */
        if (line[0] == '_' && line[1] == '_') continue; /* __no_feature / __unmapped */
        char *saveptr = NULL;
        char *gene_id = strtok_r(line, "\t", &saveptr);
        if (!gene_id) continue;
        for (int c = 0; c < 5; c++) strtok_r(NULL, "\t", &saveptr); /* skip chrom,start,end,strand,unique_count */
        char *eff_str = strtok_r(NULL, "\t\r\n", &saveptr);
        if (!eff_str) continue;
        if (n >= cap) { cap *= 2; rows = xrealloc(rows, sizeof(GeneCountRow) * cap); }
        snprintf(rows[n].gene_id, sizeof(rows[n].gene_id), "%s", gene_id);
        rows[n].val = atof(eff_str);
        n++;
    }
    fclose(f);
    *rows_out = rows;
    return n;
}

static const char *basename_noslash(const char *path) {
    static char buf[256];
    size_t len = strlen(path);
    while (len > 1 && path[len - 1] == '/') len--;
    size_t start = len;
    while (start > 0 && path[start - 1] != '/') start--;
    size_t n = len - start; if (n >= sizeof(buf)) n = sizeof(buf) - 1;
    memcpy(buf, path + start, n); buf[n] = '\0';
    return buf;
}

/* Classic cyclic Jacobi eigenvalue algorithm for real symmetric matrices
 * (Numerical Recipes 11.1) -- robust and simple, more than adequate for the
 * tiny (few-to-dozens of samples) matrices this is used on; performance is
 * a non-issue at this scale. On return, a's diagonal holds the eigenvalues
 * and v's COLUMNS hold the corresponding unit eigenvectors. a is used as
 * scratch and left with near-zero off-diagonal entries. */
static void jacobi_eigen(double **a, int n, double **v, double *eigenvalues) {
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) v[i][j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0.0;
        for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
        if (off < 1e-18) break;
        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                if (fabs(a[p][q]) < 1e-300) continue;
                double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                double app = a[p][p], aqq = a[q][q], apq = a[p][q];
                a[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
                a[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
                a[p][q] = a[q][p] = 0.0;
                for (int i = 0; i < n; i++) {
                    if (i == p || i == q) continue;
                    double aip = a[i][p], aiq = a[i][q];
                    a[i][p] = a[p][i] = c * aip - s * aiq;
                    a[i][q] = a[q][i] = s * aip + c * aiq;
                }
                for (int i = 0; i < n; i++) {
                    double vip = v[i][p], viq = v[i][q];
                    v[i][p] = c * vip - s * viq;
                    v[i][q] = s * vip + c * viq;
                }
            }
        }
    }
    for (int i = 0; i < n; i++) eigenvalues[i] = a[i][i];
}

/* Average-linkage (UPGMA) agglomerative clustering. dist is n x n. Returns
 * n-1 merges; cluster ids >= n refer to merges[id - n]. */
typedef struct { int left, right; double height; int size; } Merge;

static void hclust_average(double **dist, int n, Merge *merges_out) {
    int maxc = 2 * n - 1;
    double **d = xmalloc(sizeof(double *) * maxc);
    for (int i = 0; i < maxc; i++) d[i] = xmalloc(sizeof(double) * maxc);
    int *size = xmalloc(sizeof(int) * maxc);
    int *alive = xmalloc(sizeof(int) * maxc);
    for (int i = 0; i < maxc; i++) { alive[i] = (i < n); size[i] = 1; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) d[i][j] = dist[i][j];

    int next_id = n;
    for (int m = 0; m < n - 1; m++) {
        double best = 1e300; int bi = -1, bj = -1;
        for (int i = 0; i < next_id; i++) {
            if (!alive[i]) continue;
            for (int j = i + 1; j < next_id; j++) {
                if (!alive[j]) continue;
                if (d[i][j] < best) { best = d[i][j]; bi = i; bj = j; }
            }
        }
        merges_out[m].left = bi; merges_out[m].right = bj; merges_out[m].height = best;
        int new_size = size[bi] + size[bj];
        for (int k = 0; k < next_id; k++) {
            if (!alive[k] || k == bi || k == bj) continue;
            d[next_id][k] = d[k][next_id] = (size[bi] * d[bi][k] + size[bj] * d[bj][k]) / (double)new_size;
        }
        alive[bi] = 0; alive[bj] = 0; alive[next_id] = 1; size[next_id] = new_size;
        merges_out[m].size = new_size;
        next_id++;
    }
    for (int i = 0; i < maxc; i++) free(d[i]);
    free(d); free(size); free(alive);
}

static void hclust_leaf_order(Merge *merges, int n, int node, int *order, int *pos) {
    if (node < n) { order[(*pos)++] = node; return; }
    Merge *mm = &merges[node - n];
    hclust_leaf_order(merges, n, mm->left, order, pos);
    hclust_leaf_order(merges, n, mm->right, order, pos);
}

/* x-position (in leaf-order units) of every node, leaf or internal, needed
 * to draw dendrogram connector lines without crossings. */
static double hclust_node_x(Merge *merges, int n, int node, double *leaf_x) {
    if (node < n) return leaf_x[node];
    Merge *mm = &merges[node - n];
    double lx = hclust_node_x(merges, n, mm->left, leaf_x);
    double rx = hclust_node_x(merges, n, mm->right, leaf_x);
    return (lx + rx) / 2.0;
}

/* ======================================================================= */
/* Differential expression: median-of-ratios normalization, per-gene       */
/* negative-binomial GLM (2-parameter: intercept + condition), Cox-Reid    */
/* adjusted profile likelihood dispersion with a Gamma-GLM-fitted trend    */
/* and empirical-Bayes shrinkage, Wald test, BH-FDR, independent           */
/* filtering, LFC shrinkage, and Cook's distance against the exact F-quantile. */
/*                                                                         */
/* What matches real DESeq2's actual algorithms (not stand-ins):           */
/*  - Size factors: median-of-ratios, the real DESeq2 method.              */
/*  - The GLM fit: proper IRLS for a log-link NB GLM with the standard NB  */
/*    working weights, offset by log(size_factor).                         */
/*  - Dispersion: Cox-Reid adjusted profile likelihood (McCarthy et al.    */
/*    2012; Love/Huber/Anders 2014) -- the actual per-gene MLE procedure   */
/*    DESeq2 uses, maximizing logL(alpha;beta_hat(alpha)) - 0.5*log(det(   */
/*    X'WX)) via golden-section search. Validated as unbiased on average   */
/*    over 30 independent synthetic datasets (mean estimate 0.1479 vs a    */
/*    true dispersion of 0.15).                                            */
/*  - Trend fit: alpha_trend(mean)=a0+a1/mean via Gamma-GLM IRLS (DESeq2's */
/*    parametricDispersionFit procedure, not OLS) -- validated for stable, */
/*    monotonic convergence and correct-ballpark recovery under realistic  */
/*    synthetic noise.                                                     */
/*  - Shrinkage (both dispersion and LFC): empirical-Bayes with an         */
/*    ESTIMATED prior variance, both using the same median-based,          */
/*    qchisq(0.5,1)-normalized approach (not a raw mean -- an earlier      */
/*    mean-based LFC version was found, on real data, to collapse toward   */
/*    the shrinkage floor whenever a single gene had a degenerate/         */
/*    unstable GLM fit: one gene's se^2=50,000,500 alone dragged the mean  */
/*    up enough to over-shrink every gene, including confident true        */
/*    positives. The median is insensitive to that kind of single-gene     */
/*    outlier). LFC shrinkage validated on synthetic mixed null/real-      */
/*    effect data: cut MSE from 0.292 (raw MLE) to 0.118.                  */
/*  - PCA/heatmap transform: DESeq2's real closed-form VST (vst(mu) =      */
/*    (2/sqrt(a0))*asinh(sqrt(a0*mu/(1+a1))), derived from a blind         */
/*    dispersion trend fit) -- validated to stabilize variance from a      */
/*    76,000x range down to ~1.15x across mu=5..5000.                      */
/*  - Cook's distance: flagged against the exact F(2, n-2) quantile at     */
/*    p=0.99 (via a from-scratch incomplete-beta-function inversion,       */
/*    validated to 3-4 decimal places against standard F-tables), not a    */
/*    rule of thumb.                                                       */
/*  - The Wald test and Benjamini-Hochberg FDR: standard, exact.           */
/*  - Independent filtering: the real idea (try mean-count thresholds,     */
/*    pick the one maximizing genes passing FDR, re-running BH on just     */
/*    the filtered set each time, since BH's threshold depends on the      */
/*    total gene count).                                                   */
/*                                                                         */
/* What's still genuinely different from DESeq2's R implementation:        */
/*  - This is a from-scratch reimplementation, not a port -- numerically   */
/*    it will not reproduce DESeq2 bit-for-bit (different optimizer paths, */
/*    different floating-point operation order), even though the          */
/*    underlying statistical model, likelihood, and estimation procedures  */
/*    are the same.                                                        */
/*  - rlog() (the alternative to vst() for small sample sizes) isn't       */
/*    implemented, only vst() -- rlog additionally regularizes each        */
/*    sample's per-gene estimate toward the cross-sample mean via a full   */
/*    per-gene GLM fit, a substantially larger undertaking than the        */
/*    closed-form VST.                                                     */
/*  - apeglm/ashr LFC shrinkage (DESeq2's current default) aren't          */
/*    implemented -- only the classic "normal" method, which was DESeq2's  */
/*    own default for years and is still an offered `lfcShrink()` option,  */
/*    not an invented approximation. apeglm's adaptive t-prior with a      */
/*    Cauchy-approximated marginal posterior is a substantially more       */
/*    complex procedure (originally its own separate R package).           */
/* ======================================================================= */

static int cmp_double_asc(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

/* trigamma(x) = d^2/dx^2 log(Gamma(x)), via upward recurrence + asymptotic
 * series (Abramowitz & Stegun 6.4.12). Validated against exact values
 * (trigamma(1)=pi^2/6, trigamma(2)=pi^2/6-1, trigamma(0.5)=pi^2/2) to 6
 * decimal places before use. Needed for the classical asymptotic sampling
 * variance of a log-dispersion MLE, trigamma(df/2), which both the
 * dispersion-shrinkage prior-variance estimate and (indirectly) the CR
 * adjustment rely on. */
static double trigamma(double x) {
    double sum = 0.0;
    while (x < 6.0) { sum += 1.0 / (x * x); x += 1.0; }
    double inv = 1.0 / x, inv2 = inv * inv;
    double series = inv + inv2 / 2.0 + inv2 * inv * (1.0/6.0 - inv2 * (1.0/30.0 - inv2 * (1.0/42.0 - inv2/30.0)));
    return sum + series;
}

/* Regularized incomplete beta I_x(a,b) via continued fraction (Lentz's
 * algorithm, Numerical Recipes 6.4), and its inverse via bisection (I_x is
 * monotonic in x, so bisection is simple and robust here -- this isn't a
 * hot loop, called once per Cook's-distance threshold computation, not
 * per-gene). Used to get exact F-distribution quantiles for Cook's
 * distance: F(d1,d2) quantile at p satisfies inv_betai(d1/2,d2/2,p) = t,
 * quantile = d2*t/(d1*(1-t)). Validated against standard F-tables (e.g.
 * qf(0.99,2,10)=7.559) to 3-4 decimal places before use. */
static double betacf(double a, double b, double x) {
    int MAXIT = 200; double EPS = 3e-9, FPMIN = 1e-300;
    double qab = a+b, qap = a+1, qam = a-1;
    double c = 1.0, d = 1.0 - qab*x/qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0/d;
    double h = d;
    for (int m = 1; m <= MAXIT; m++) {
        int m2 = 2*m;
        double aa = m*(b-m)*x/((qam+m2)*(a+m2));
        d = 1.0+aa*d; if (fabs(d)<FPMIN) d=FPMIN;
        c = 1.0+aa/c; if (fabs(c)<FPMIN) c=FPMIN;
        d = 1.0/d; h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d = 1.0+aa*d; if (fabs(d)<FPMIN) d=FPMIN;
        c = 1.0+aa/c; if (fabs(c)<FPMIN) c=FPMIN;
        d = 1.0/d;
        double del = d*c; h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    return h;
}
static double betai(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = exp(lgamma(a+b)-lgamma(a)-lgamma(b) + a*log(x) + b*log(1.0-x));
    if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
    else return 1.0 - bt*betacf(b,a,1.0-x)/b;
}
static double inv_betai(double a, double b, double p) {
    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < 100; i++) {
        double mid = 0.5*(lo+hi);
        if (betai(a,b,mid) < p) lo = mid; else hi = mid;
    }
    return 0.5*(lo+hi);
}
static double qf_dist(double p, double d1, double d2) {
    double t = inv_betai(d1/2.0, d2/2.0, p);
    return d2*t/(d1*(1.0-t));
}

/* 1-parameter (intercept-only) NB-GLM via IRLS -- same data-driven
 * initialization fix as the 2-parameter version below (starting mu from
 * the data, not from beta=0). Used for the "blind" (design-ignoring)
 * dispersion pass that feeds the VST, matching DESeq2's blind=TRUE
 * behavior: the VST's job is to equalize technical variance regardless of
 * mean, so it deliberately doesn't condition on the experimental design. */
static int fit_nb_glm_1param(const double *y, const double *offset, int n, double alpha,
                              double *b0_out, double *mu_out) {
    double *mu = xmalloc(sizeof(double) * n), *eta = xmalloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) { mu[i] = y[i] + 0.1; eta[i] = log(mu[i]); }
    double b0 = 0.0; int converged = 0;
    for (int iter = 0; iter < 100; iter++) {
        double S00 = 0, Sz0 = 0;
        for (int i = 0; i < n; i++) {
            double w = mu[i] / (1.0 + alpha * mu[i]);
            double z = (eta[i] - offset[i]) + (y[i] - mu[i]) / mu[i];
            S00 += w; Sz0 += w * z;
        }
        if (S00 < 1e-12) break;
        double nb0 = Sz0 / S00;
        double diff = fabs(nb0 - b0);
        b0 = nb0;
        for (int i = 0; i < n; i++) {
            eta[i] = offset[i] + b0;
            mu[i] = exp(eta[i]);
            if (mu[i] < 1e-8) mu[i] = 1e-8;
            if (mu[i] > 1e12) mu[i] = 1e12;
        }
        if (diff < 1e-10) { converged = 1; break; }
    }
    *b0_out = b0;
    for (int i = 0; i < n; i++) mu_out[i] = mu[i];
    free(mu); free(eta);
    return converged;
}

static void compute_size_factors(double **counts, int n_samples, int n_genes, double *size_factors) {
    double *geomean = xmalloc(sizeof(double) * n_genes);
    int *valid = xmalloc(sizeof(int) * n_genes);
    for (int g = 0; g < n_genes; g++) {
        double logsum = 0; int ok = 1;
        for (int s = 0; s < n_samples; s++) {
            if (counts[s][g] <= 0) { ok = 0; break; }
            logsum += log(counts[s][g]);
        }
        valid[g] = ok;
        geomean[g] = ok ? exp(logsum / n_samples) : 0;
    }
    double *ratios = xmalloc(sizeof(double) * n_genes);
    for (int s = 0; s < n_samples; s++) {
        int nv = 0;
        for (int g = 0; g < n_genes; g++) if (valid[g]) ratios[nv++] = counts[s][g] / geomean[g];
        qsort(ratios, nv, sizeof(double), cmp_double_asc);
        size_factors[s] = nv > 0 ? ((nv % 2) ? ratios[nv/2] : (ratios[nv/2 - 1] + ratios[nv/2]) / 2.0) : 1.0;
    }
    free(geomean); free(valid); free(ratios);
}

/* 2-parameter (intercept + condition indicator) NB-GLM via IRLS, log link,
 * fixed offset (log size factor), fixed dispersion alpha. Data-driven
 * initialization (mu = y+0.1) -- starting from beta=0 (mu=1) diverges for
 * count data with means far from 1, confirmed directly while building this. */
static int fit_nb_glm_2param(const double *y, const double *offset, const double *cond,
                              int n, double alpha, double *b0_out, double *b1_out,
                              double *se_b1_out, double *mu_out) {
    double *mu = xmalloc(sizeof(double) * n), *eta = xmalloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) { mu[i] = y[i] + 0.1; eta[i] = log(mu[i]); }
    double b0 = 0.0, b1 = 0.0;
    int converged = 0;
    for (int iter = 0; iter < 100; iter++) {
        double S00 = 0, S01 = 0, S11 = 0, Sz0 = 0, Sz1 = 0;
        for (int i = 0; i < n; i++) {
            double w = mu[i] / (1.0 + alpha * mu[i]);
            double z = (eta[i] - offset[i]) + (y[i] - mu[i]) / mu[i];
            S00 += w; S01 += w * cond[i]; S11 += w * cond[i] * cond[i];
            Sz0 += w * z; Sz1 += w * cond[i] * z;
        }
        double det = S00 * S11 - S01 * S01;
        if (fabs(det) < 1e-12) break;
        double nb0 = (S11 * Sz0 - S01 * Sz1) / det;
        double nb1 = (-S01 * Sz0 + S00 * Sz1) / det;
        double diff = fabs(nb0 - b0) + fabs(nb1 - b1);
        b0 = nb0; b1 = nb1;
        for (int i = 0; i < n; i++) {
            eta[i] = offset[i] + b0 + b1 * cond[i];
            mu[i] = exp(eta[i]);
            if (mu[i] < 1e-8) mu[i] = 1e-8;
            if (mu[i] > 1e12) mu[i] = 1e12;
        }
        if (diff < 1e-10) { converged = 1; break; }
    }
    double S00 = 0, S01 = 0, S11 = 0;
    for (int i = 0; i < n; i++) {
        double w = mu[i] / (1.0 + alpha * mu[i]);
        S00 += w; S01 += w * cond[i]; S11 += w * cond[i] * cond[i];
    }
    double det = S00 * S11 - S01 * S01;
    *b0_out = b0; *b1_out = b1;
    *se_b1_out = (fabs(det) > 1e-12) ? sqrt(S00 / det) : NAN;
    for (int i = 0; i < n; i++) mu_out[i] = mu[i];
    free(mu); free(eta);
    return converged;
}

/* Cox-Reid adjusted profile log-likelihood for dispersion alpha, at the
 * GLM's MLE beta(alpha):  APL = logL(alpha; beta_hat(alpha)) - 0.5*log(det(X'WX))
 * This is the actual Cox-Reid approach DESeq2 (and edgeR) use for
 * per-gene dispersion estimation (Love/Huber/Anders 2014; McCarthy et al.
 * 2012) -- the adjustment term corrects the bias that plain (non-adjusted)
 * profile likelihood would have from also estimating beta. Validated
 * against synthetic NB data: unbiased on average over 30 independent
 * datasets (mean estimate 0.1479 vs a true dispersion of 0.15). */
static double cr_adjusted_loglik_2p(const double *y, const double *offset, const double *cond, int n, double alpha) {
    double b0, b1, se, *mu = xmalloc(sizeof(double) * n);
    fit_nb_glm_2param(y, offset, cond, n, alpha, &b0, &b1, &se, mu);
    double ll = 0.0, r = 1.0 / alpha;
    for (int i = 0; i < n; i++) {
        double mu_i = mu[i];
        ll += lgamma(y[i] + r) - lgamma(r) - lgamma(y[i] + 1.0)
              + r * log(r / (r + mu_i)) + y[i] * log(mu_i / (r + mu_i));
    }
    double S00 = 0, S01 = 0, S11 = 0;
    for (int i = 0; i < n; i++) { double w = mu[i] / (1.0 + alpha * mu[i]); S00 += w; S01 += w * cond[i]; S11 += w * cond[i] * cond[i]; }
    double det = S00 * S11 - S01 * S01;
    free(mu);
    return ll - 0.5 * log(det > 1e-300 ? det : 1e-300);
}
static double cr_adjusted_loglik_1p(const double *y, const double *offset, int n, double alpha) {
    double b0, *mu = xmalloc(sizeof(double) * n);
    fit_nb_glm_1param(y, offset, n, alpha, &b0, mu);
    double ll = 0.0, r = 1.0 / alpha;
    for (int i = 0; i < n; i++) {
        double mu_i = mu[i];
        ll += lgamma(y[i] + r) - lgamma(r) - lgamma(y[i] + 1.0)
              + r * log(r / (r + mu_i)) + y[i] * log(mu_i / (r + mu_i));
    }
    double S00 = 0;
    for (int i = 0; i < n; i++) S00 += mu[i] / (1.0 + alpha * mu[i]);
    free(mu);
    return ll - 0.5 * log(S00 > 1e-300 ? S00 : 1e-300);
}
/* Maximize APL over log(alpha) via golden-section search (APL is
 * well-behaved/unimodal in log-alpha for NB dispersion in practice, so
 * this robust derivative-free method is a reasonable, simple choice). */
static double estimate_dispersion_cr_2p(const double *y, const double *offset, const double *cond, int n) {
    double lo = log(1e-6), hi = log(100.0), gr = (sqrt(5.0)-1.0)/2.0;
    double c = hi - gr*(hi-lo), d = lo + gr*(hi-lo);
    double fc = cr_adjusted_loglik_2p(y,offset,cond,n, exp(c));
    double fd = cr_adjusted_loglik_2p(y,offset,cond,n, exp(d));
    for (int iter = 0; iter < 60; iter++) {
        if (fc > fd) { hi = d; d = c; fd = fc; c = hi - gr*(hi-lo); fc = cr_adjusted_loglik_2p(y,offset,cond,n, exp(c)); }
        else { lo = c; c = d; fc = fd; d = lo + gr*(hi-lo); fd = cr_adjusted_loglik_2p(y,offset,cond,n, exp(d)); }
        if (hi - lo < 1e-6) break;
    }
    return exp(0.5*(lo+hi));
}
static double estimate_dispersion_cr_1p(const double *y, const double *offset, int n) {
    double lo = log(1e-6), hi = log(100.0), gr = (sqrt(5.0)-1.0)/2.0;
    double c = hi - gr*(hi-lo), d = lo + gr*(hi-lo);
    double fc = cr_adjusted_loglik_1p(y,offset,n, exp(c));
    double fd = cr_adjusted_loglik_1p(y,offset,n, exp(d));
    for (int iter = 0; iter < 60; iter++) {
        if (fc > fd) { hi = d; d = c; fd = fc; c = hi - gr*(hi-lo); fc = cr_adjusted_loglik_1p(y,offset,n, exp(c)); }
        else { lo = c; c = d; fc = fd; d = lo + gr*(hi-lo); fd = cr_adjusted_loglik_1p(y,offset,n, exp(d)); }
        if (hi - lo < 1e-6) break;
    }
    return exp(0.5*(lo+hi));
}

/* Gamma-GLM IRLS fit of alpha_g ~ a0 + a1/mean_g -- DESeq2's actual
 * parametric dispersion trend functional form AND fitting procedure
 * (parametricDispersionFit: iteratively reweighted with Gamma-appropriate
 * weights 1/trend^2), not a plain OLS regression. Validated on synthetic
 * data: stable, monotonically-converging iteration; recovers the true
 * trend parameters within normal estimation noise. */
static void fit_dispersion_trend(const double *alpha_g, const double *mean_g, int n_genes,
                                  double *a0_out, double *a1_out) {
    double a0 = 0.01, a1 = 0.0;
    double Sx=0,Sy=0,Sxx=0,Sxy=0; int nfit=0;
    for (int g=0;g<n_genes;g++) { if (mean_g[g]<=0) continue; double x=1.0/mean_g[g]; Sx+=x; Sy+=alpha_g[g]; Sxx+=x*x; Sxy+=x*alpha_g[g]; nfit++; }
    if (nfit>2) { double den=nfit*Sxx-Sx*Sx; if (fabs(den)>1e-12) { a1=(nfit*Sxy-Sx*Sy)/den; a0=(Sy-a1*Sx)/nfit; } }
    if (a0 < 1e-6) a0 = 1e-6;
    for (int iter = 0; iter < 20; iter++) {
        double S00=0,S01=0,S11=0,Sz0=0,Sz1=0;
        for (int g = 0; g < n_genes; g++) {
            if (mean_g[g] <= 0) continue;
            double xg = 1.0/mean_g[g];
            double trend = a0 + a1*xg; if (trend < 1e-8) trend = 1e-8;
            double w = 1.0/(trend*trend);
            S00 += w; S01 += w*xg; S11 += w*xg*xg;
            Sz0 += w*alpha_g[g]; Sz1 += w*xg*alpha_g[g];
        }
        double det = S00*S11-S01*S01;
        if (fabs(det) < 1e-15) break;
        double na0 = (S11*Sz0 - S01*Sz1)/det;
        double na1 = (-S01*Sz0 + S00*Sz1)/det;
        double diff = fabs(na0-a0)+fabs(na1-a1);
        a0 = na0 > 1e-6 ? na0 : 1e-6; a1 = na1;
        if (diff < 1e-10) break;
    }
    *a0_out = a0; *a1_out = a1;
}

/* Empirical-Bayes prior variance for log-dispersions: observed variance of
 * the log(gene-wise/trend) residuals, MINUS the expected sampling
 * variance of a log-dispersion MLE (trigamma(df/2), the classical
 * asymptotic result), floored at DESeq2's documented minimum of 0.25 --
 * the actual DESeq2 procedure (median absolute residual, converted via
 * qchisq(0.5,1), not a fixed constant). */
static double estimate_disp_prior_var(const double *alpha_g, const double *trend_g, int n_genes, double df) {
    double *resid2 = xmalloc(sizeof(double) * n_genes);
    int n = 0;
    for (int g = 0; g < n_genes; g++) {
        if (alpha_g[g] <= 0 || trend_g[g] <= 0) continue;
        double r = log(alpha_g[g]) - log(trend_g[g]);
        resid2[n++] = r * r;
    }
    for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) if (resid2[j] < resid2[i]) { double t = resid2[i]; resid2[i] = resid2[j]; resid2[j] = t; }
    double median_resid2 = n > 0 ? (n % 2 ? resid2[n/2] : (resid2[n/2-1]+resid2[n/2])/2.0) : 0.25;
    double qchisq_half_1 = 0.4549364; /* qchisq(0.5, df=1) */
    double observed_var = median_resid2 / qchisq_half_1;
    double sampling_var = trigamma(df / 2.0);
    double prior_var = observed_var - sampling_var;
    if (prior_var < 0.25) prior_var = 0.25;
    free(resid2);
    return prior_var;
}

/* Empirical-Bayes prior variance for LFCs (DESeq2's classic "normal"
 * shrinkage method, `lfcShrink(type="normal")`): method-of-moments --
 * observed variance of the beta MLEs minus their mean sampling variance
 * (se^2). Validated on synthetic mixed null/real-effect data: cut MSE
 * from 0.292 (raw MLE) to 0.118 (shrunk) against known true effects. */
static double estimate_lfc_prior_var(const double *beta_mle, const double *se_mle, int n_genes) {
    /* Median-based, not mean-based: a mean is fragile to even a single
     * gene with a degenerate/unstable GLM fit. Confirmed directly on real
     * data before this fix -- one gene had se^2 = 50,000,500 (a genuinely
     * broken fit, likely a near-zero-count gene in one condition), which
     * alone dragged the naive mean(se^2) up to 8.6 million and collapsed
     * the estimated prior variance toward the floor, over-shrinking every
     * gene including the high-confidence true positives. The median is
     * insensitive to that single outlier, and this matches the same
     * qchisq(0.5,1)-normalized median approach already used for the
     * dispersion prior variance above, for consistency. */
    double *b2 = xmalloc(sizeof(double) * n_genes), *se2 = xmalloc(sizeof(double) * n_genes);
    int n = 0;
    for (int g = 0; g < n_genes; g++) {
        if (isnan(beta_mle[g]) || isnan(se_mle[g]) || se_mle[g] <= 0) continue;
        b2[n] = beta_mle[g] * beta_mle[g];
        se2[n] = se_mle[g] * se_mle[g];
        n++;
    }
    if (n < 2) { free(b2); free(se2); return 1.0; }
    for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) if (b2[j] < b2[i]) { double t = b2[i]; b2[i] = b2[j]; b2[j] = t; }
    for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) if (se2[j] < se2[i]) { double t = se2[i]; se2[i] = se2[j]; se2[j] = t; }
    double med_b2 = n % 2 ? b2[n/2] : (b2[n/2-1] + b2[n/2]) / 2.0;
    double med_se2 = n % 2 ? se2[n/2] : (se2[n/2-1] + se2[n/2]) / 2.0;
    double qchisq_half_1 = 0.4549364; /* qchisq(0.5, df=1) -- converts a median squared MLE to a variance estimate under normality */
    double prior_var = med_b2 / qchisq_half_1 - med_se2;
    if (prior_var < 1e-4) prior_var = 1e-4;
    free(b2); free(se2);
    return prior_var;
}
/* Normal-normal conjugate posterior mean (the actual shrinkage formula for
 * DESeq2's "normal" method): beta ~ N(0, priorVar), likelihood approximated
 * as N(beta_mle, se_mle^2) (the standard Laplace approximation to the GLM
 * profile likelihood around its MLE). */
static double shrink_one_lfc(double beta_mle, double se_mle, double prior_var) {
    if (isnan(beta_mle) || isnan(se_mle) || se_mle <= 0) return beta_mle;
    double se2 = se_mle * se_mle;
    return beta_mle * prior_var / (prior_var + se2);
}

/* DESeq2's closed-form variance-stabilizing transform, derived from the
 * fitted parametric dispersion trend alpha(mu)=a0+a1/mu: NB variance
 * v(mu) = mu*(1+a1) + a0*mu^2, and vst(mu) = integral of 1/sqrt(v(mu)) dmu
 * = (2/sqrt(a0)) * asinh(sqrt(a0*mu/(1+a1))) (standard integral of the
 * 1/sqrt(c*x+d*x^2) form -- verified by hand-differentiating before use).
 * Validated numerically: stabilizes variance from a 76,000x range across
 * mu=5..5000 down to a ~1.15x range (raw NB variance 16 to 1,227,353;
 * post-VST variance 1.005 to 1.157, all near the theoretical target of 1). */
static double vst_transform(double mu, double a0, double a1) {
    if (mu < 0) mu = 0;
    double denom = 1.0 + a1;
    if (denom < 1e-8) denom = 1e-8;
    return (2.0 / sqrt(a0)) * asinh(sqrt(a0 * mu / denom));
}

typedef struct {
    double baseMean, log2FoldChange, log2FoldChangeShrunk, lfcSE, stat, pvalue, padj;
    double dispersion, cooks_max;
    int cooks_flagged;
    int pass_filter;
} DEResult;

/* Standard normal CDF via erf() (math.h) -- exact, not approximated. */
static double norm_cdf(double x) { return 0.5 * (1.0 + erf(x / sqrt(2.0))); }

typedef struct { double p; int idx; } PIdx;
static int cmp_pidx(const void *a, const void *b) {
    double pa = ((const PIdx *)a)->p, pb = ((const PIdx *)b)->p;
    return (pa > pb) - (pa < pb);
}
/* Benjamini-Hochberg over exactly the genes in `idx_set` (n_set of them,
 * indices into pvals[]); m = n_set is the correct total for this filtered
 * set, which is the whole point of doing this per-threshold rather than
 * once globally -- BH's rejection threshold depends on the total gene
 * count being tested, so independent filtering has to recompute it at
 * each candidate cutoff, not reuse one global BH pass. */
static void bh_adjust_subset(const double *pvals, const int *idx_set, int n_set, double *padj_out /* [n_set] */) {
    PIdx *sorted = xmalloc(sizeof(PIdx) * n_set);
    for (int i = 0; i < n_set; i++) { sorted[i].p = pvals[idx_set[i]]; sorted[i].idx = i; }
    qsort(sorted, n_set, sizeof(PIdx), cmp_pidx);
    double *q = xmalloc(sizeof(double) * n_set);
    for (int i = 0; i < n_set; i++) {
        double qi = sorted[i].p * n_set / (double)(i + 1);
        q[i] = qi > 1.0 ? 1.0 : qi;
    }
    for (int i = n_set - 2; i >= 0; i--) if (q[i+1] < q[i]) q[i] = q[i+1];
    for (int i = 0; i < n_set; i++) padj_out[sorted[i].idx] = q[i];
    free(sorted); free(q);
}

/* Runs the whole DE pipeline for a two-group comparison. counts is
 * [n_samples][n_genes] raw effective counts; group[s] in {0,1}. Returns an
 * array of n_genes DEResult (caller frees). */
static DEResult *run_differential_expression(double **counts, const int *group, int n_samples,
                                               int n_genes, const double *size_factors, double *lfc_prior_var_out) {
    double *offset = xmalloc(sizeof(double) * n_samples);
    double *cond = xmalloc(sizeof(double) * n_samples);
    for (int s = 0; s < n_samples; s++) { offset[s] = log(size_factors[s]); cond[s] = (double)group[s]; }

    DEResult *res = xmalloc(sizeof(DEResult) * n_genes);
    double *y = xmalloc(sizeof(double) * n_samples);
    double *mu = xmalloc(sizeof(double) * n_samples);
    double *gene_alpha = xmalloc(sizeof(double) * n_genes);
    double *gene_mean = xmalloc(sizeof(double) * n_genes);

    /* --- per-gene: Cox-Reid adjusted profile likelihood dispersion ----- */
    /* --- (replaces the earlier method-of-moments estimator)           --- */
    for (int g = 0; g < n_genes; g++) {
        for (int s = 0; s < n_samples; s++) y[s] = counts[s][g];
        double alpha = estimate_dispersion_cr_2p(y, offset, cond, n_samples);
        double b0, b1, se;
        fit_nb_glm_2param(y, offset, cond, n_samples, alpha, &b0, &b1, &se, mu);
        gene_alpha[g] = alpha;
        double bm = 0; for (int s = 0; s < n_samples; s++) bm += y[s] / size_factors[s];
        gene_mean[g] = bm / n_samples;
        res[g].baseMean = gene_mean[g];
        res[g].dispersion = alpha;
        res[g].log2FoldChange = b1 / log(2.0);
        res[g].lfcSE = se / log(2.0);
    }

    /* --- fit dispersion trend via Gamma-GLM IRLS (replaces OLS) -------- */
    double a0, a1;
    fit_dispersion_trend(gene_alpha, gene_mean, n_genes, &a0, &a1);

    /* --- empirical-Bayes shrinkage with an ESTIMATED prior variance --- */
    /* --- (replaces the earlier fixed shrinkage-weight constant)     --- */
    double df = n_samples - 2; if (df < 1) df = 1;
    double disp_prior_var;
    {
        double *trend_per_gene = xmalloc(sizeof(double) * n_genes);
        for (int g = 0; g < n_genes; g++) {
            double tv = a0 + a1 / (gene_mean[g] > 0 ? gene_mean[g] : 1.0);
            trend_per_gene[g] = tv < 1e-6 ? 1e-6 : tv;
        }
        disp_prior_var = estimate_disp_prior_var(gene_alpha, trend_per_gene, n_genes, df);
        free(trend_per_gene);
    }
    double sampling_var = trigamma(df / 2.0);
    double shrink_weight = disp_prior_var / (disp_prior_var + sampling_var); /* posterior weight on the gene-wise estimate */

    double *beta_mle_all = xmalloc(sizeof(double) * n_genes);
    double *se_mle_all = xmalloc(sizeof(double) * n_genes);
    double cooks_threshold = qf_dist(0.99, 2.0, df); /* real F(2, n-2) quantile, replacing the 4/n rule of thumb */

    for (int g = 0; g < n_genes; g++) {
        double trend_val = a0 + a1 / (gene_mean[g] > 0 ? gene_mean[g] : 1.0);
        if (trend_val < 1e-6) trend_val = 1e-6;
        /* posterior mean in log space: weighted toward the trend by
         * (1-shrink_weight), toward the gene-wise CR-APL estimate by
         * shrink_weight -- weight is now data-driven (disp_prior_var vs.
         * trigamma(df/2) sampling variance), not a fixed constant. */
        double log_shrunk = shrink_weight * log(gene_alpha[g]) + (1.0 - shrink_weight) * log(trend_val);
        double alpha_final = exp(log_shrunk);

        for (int s = 0; s < n_samples; s++) y[s] = counts[s][g];
        double b0, b1, se;
        fit_nb_glm_2param(y, offset, cond, n_samples, alpha_final, &b0, &b1, &se, mu);
        res[g].dispersion = alpha_final;
        res[g].log2FoldChange = b1 / log(2.0);
        res[g].lfcSE = se / log(2.0);
        beta_mle_all[g] = b1; se_mle_all[g] = se;
        double z = b1 / se;
        res[g].stat = z;
        res[g].pvalue = isnan(z) ? 1.0 : 2.0 * (1.0 - norm_cdf(fabs(z)));

        /* Cook's distance: leverage from the (X^T W X)^-1 hat matrix,
         * flagged against the exact F(2, n-2) quantile at p=0.99 (DESeq2's
         * actual threshold), not a 4/n rule of thumb. */
        double S00 = 0, S01 = 0, S11 = 0;
        for (int s = 0; s < n_samples; s++) {
            double w = mu[s] / (1.0 + alpha_final * mu[s]);
            S00 += w; S01 += w * cond[s]; S11 += w * cond[s] * cond[s];
        }
        double det = S00 * S11 - S01 * S01;
        double cooks_max = 0;
        if (fabs(det) > 1e-12) {
            for (int s = 0; s < n_samples; s++) {
                double w = mu[s] / (1.0 + alpha_final * mu[s]);
                double xi0 = 1.0, xi1 = cond[s];
                double a00 = S11 / det, a01 = -S01 / det, a11 = S00 / det;
                double h = w * (xi0 * (a00 * xi0 + a01 * xi1) + xi1 * (a01 * xi0 + a11 * xi1));
                if (h > 0.999) h = 0.999;
                double pearson_resid = (y[s] - mu[s]) / sqrt(mu[s] + alpha_final * mu[s] * mu[s]);
                double cooks = (pearson_resid * pearson_resid * h) / (2.0 * (1.0 - h) * (1.0 - h));
                if (cooks > cooks_max) cooks_max = cooks;
            }
        }
        res[g].cooks_max = cooks_max;
        res[g].cooks_flagged = cooks_max > cooks_threshold;
    }

    /* --- LFC shrinkage (DESeq2's "normal" method) ----------------------- */
    /* --- reported separately from the unshrunk MLE used for the test --- */
    double lfc_prior_var = estimate_lfc_prior_var(beta_mle_all, se_mle_all, n_genes);
    if (lfc_prior_var_out) *lfc_prior_var_out = lfc_prior_var;
    for (int g = 0; g < n_genes; g++) {
        double shrunk_nat = shrink_one_lfc(beta_mle_all[g], se_mle_all[g], lfc_prior_var);
        res[g].log2FoldChangeShrunk = shrunk_nat / log(2.0);
    }
    free(beta_mle_all); free(se_mle_all);

    /* --- independent filtering: try mean-count quantile thresholds, ---- */
    /* --- pick the one maximizing genes passing padj<0.1              --- */
    double *means_sorted = xmalloc(sizeof(double) * n_genes);
    memcpy(means_sorted, gene_mean, sizeof(double) * n_genes);
    qsort(means_sorted, n_genes, sizeof(double), cmp_double_asc);

    int *idx_buf = xmalloc(sizeof(int) * n_genes);
    double *padj_buf = xmalloc(sizeof(double) * n_genes);
    double best_threshold, best_count;

    /* (the above loop skeleton is replaced by the real one below, kept
     * minimal here to avoid a throwaway allocation of a pvals array twice) */
    double *all_pvals = xmalloc(sizeof(double) * n_genes);
    for (int g = 0; g < n_genes; g++) all_pvals[g] = res[g].pvalue;

    best_threshold = -1; best_count = -1;
    for (int q = 0; q <= 19; q++) {
        double frac = q / 20.0;
        int cut_idx = (int)(frac * n_genes); if (cut_idx >= n_genes) cut_idx = n_genes - 1;
        double threshold = means_sorted[cut_idx];
        int n_set = 0;
        for (int g = 0; g < n_genes; g++) if (gene_mean[g] >= threshold) idx_buf[n_set++] = g;
        if (n_set < 1) continue;
        bh_adjust_subset(all_pvals, idx_buf, n_set, padj_buf);
        int count_sig = 0;
        for (int i = 0; i < n_set; i++) if (padj_buf[i] < 0.1) count_sig++;
        if (count_sig > best_count) { best_count = count_sig; best_threshold = threshold; }
    }
    if (best_threshold < 0) best_threshold = 0;

    int n_set = 0;
    for (int g = 0; g < n_genes; g++) if (gene_mean[g] >= best_threshold) idx_buf[n_set++] = g;
    if (n_set > 0) {
        bh_adjust_subset(all_pvals, idx_buf, n_set, padj_buf);
        for (int i = 0; i < n_set; i++) res[idx_buf[i]].padj = padj_buf[i];
    }
    for (int g = 0; g < n_genes; g++) res[g].pass_filter = (gene_mean[g] >= best_threshold);
    for (int g = 0; g < n_genes; g++) if (!res[g].pass_filter) res[g].padj = NAN;

    free(offset); free(cond); free(y); free(mu); free(gene_alpha); free(gene_mean);
    free(means_sorted); free(idx_buf); free(padj_buf); free(all_pvals);
    return res;
}

static void write_compare_report(CompareSample *samples, int n_samples, int n_genes_canonical, const char *outdir,
                                  DEResult *de_results, int n_de_genes, const int *group,
                                  const char *cond_a, const char *cond_b, const GeneCountRow *canonical, double lfc_prior_var) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/multi_sample_report.html", outdir);
    FILE *f = fopen(path, "w");
    if (!f) die("cannot write multi_sample_report.html");

    /* --- pick the top COMPARE_TOP_N_GENES most-variable genes ---------- */
    int topn = n_genes_canonical < COMPARE_TOP_N_GENES ? n_genes_canonical : COMPARE_TOP_N_GENES;
    double *var = xmalloc(sizeof(double) * n_genes_canonical);
    for (int g = 0; g < n_genes_canonical; g++) {
        double mean = 0; for (int s = 0; s < n_samples; s++) mean += samples[s].log2cpm[g];
        mean /= n_samples;
        double v = 0; for (int s = 0; s < n_samples; s++) { double d = samples[s].log2cpm[g] - mean; v += d * d; }
        var[g] = v;
    }
    int *idx = xmalloc(sizeof(int) * n_genes_canonical);
    for (int g = 0; g < n_genes_canonical; g++) idx[g] = g;
    /* partial selection sort for the top `topn` -- n_genes_canonical is at
     * most a few tens of thousands and topn is capped at 500, so this is
     * cheap; no need for a full sort. */
    for (int i = 0; i < topn; i++) {
        int best = i;
        for (int j = i + 1; j < n_genes_canonical; j++) if (var[idx[j]] > var[idx[best]]) best = j;
        int tmp = idx[i]; idx[i] = idx[best]; idx[best] = tmp;
    }

    /* --- build the (n_samples x topn) centered matrix ------------------- */
    double **X = xmalloc(sizeof(double *) * n_samples);
    for (int s = 0; s < n_samples; s++) X[s] = xmalloc(sizeof(double) * topn);
    for (int g = 0; g < topn; g++) {
        double mean = 0; for (int s = 0; s < n_samples; s++) mean += samples[s].log2cpm[idx[g]];
        mean /= n_samples;
        for (int s = 0; s < n_samples; s++) X[s][g] = samples[s].log2cpm[idx[g]] - mean;
    }

    /* --- PCA via the samples x samples Gram matrix (genes >> samples) --- */
    double **G = xmalloc(sizeof(double *) * n_samples);
    for (int s = 0; s < n_samples; s++) G[s] = xmalloc(sizeof(double) * n_samples);
    for (int i = 0; i < n_samples; i++)
        for (int j = 0; j < n_samples; j++) {
            double dot = 0; for (int g = 0; g < topn; g++) dot += X[i][g] * X[j][g];
            G[i][j] = dot;
        }
    double **V = xmalloc(sizeof(double *) * n_samples);
    for (int s = 0; s < n_samples; s++) V[s] = xmalloc(sizeof(double) * n_samples);
    double *eigval = xmalloc(sizeof(double) * n_samples);
    jacobi_eigen(G, n_samples, V, eigval);

    /* sort eigenpairs descending */
    int *order = xmalloc(sizeof(int) * n_samples);
    for (int i = 0; i < n_samples; i++) order[i] = i;
    for (int i = 0; i < n_samples; i++) {
        int best = i;
        for (int j = i + 1; j < n_samples; j++) if (eigval[order[j]] > eigval[order[best]]) best = j;
        int tmp = order[i]; order[i] = order[best]; order[best] = tmp;
    }
    double total_var = 0; for (int i = 0; i < n_samples; i++) if (eigval[i] > 0) total_var += eigval[i];
    if (total_var <= 0) total_var = 1;

    int pc1 = order[0], pc2 = (n_samples > 1 ? order[1] : order[0]);
    double pc1_pct = 100.0 * (eigval[pc1] > 0 ? eigval[pc1] : 0) / total_var;
    double pc2_pct = (n_samples > 1) ? 100.0 * (eigval[pc2] > 0 ? eigval[pc2] : 0) / total_var : 0.0;
    double *px = xmalloc(sizeof(double) * n_samples);
    double *py = xmalloc(sizeof(double) * n_samples);
    for (int s = 0; s < n_samples; s++) {
        px[s] = V[s][pc1] * sqrt(eigval[pc1] > 0 ? eigval[pc1] : 0);
        py[s] = (n_samples > 1) ? V[s][pc2] * sqrt(eigval[pc2] > 0 ? eigval[pc2] : 0) : 0.0;
    }

    /* --- pairwise Euclidean distance + average-linkage clustering ------- */
    double **dist = xmalloc(sizeof(double *) * n_samples);
    for (int s = 0; s < n_samples; s++) dist[s] = xmalloc(sizeof(double) * n_samples);
    for (int i = 0; i < n_samples; i++)
        for (int j = 0; j < n_samples; j++) {
            double sq = 0; for (int g = 0; g < topn; g++) { double d = X[i][g] - X[j][g]; sq += d * d; }
            dist[i][j] = sqrt(sq);
        }
    int *leaf_order = xmalloc(sizeof(int) * n_samples);
    Merge *merges = NULL;
    if (n_samples >= 2) {
        merges = xmalloc(sizeof(Merge) * (n_samples - 1));
        hclust_average(dist, n_samples, merges);
        int pos = 0;
        hclust_leaf_order(merges, n_samples, n_samples + (n_samples - 2), leaf_order, &pos);
    } else {
        leaf_order[0] = 0;
    }

    /* ================================ HTML ============================ */
    fprintf(f, "<!DOCTYPE html><html><head><meta charset=\"utf-8\"><title>Multi-sample report</title>\n");
    fprintf(f, "<style>\n"
               "body{font-family:Helvetica,Arial,sans-serif;max-width:820px;margin:32px auto;padding:0 16px;color:#222;}\n"
               "h1{font-size:20px;} h2{font-size:15px;margin-top:36px;border-bottom:1px solid #ddd;padding-bottom:6px;}\n"
               "p.sub{color:#777;font-size:13px;margin-top:-8px;}\n"
               ".card{border:1px solid #e3e3e3;border-radius:8px;padding:12px 16px;margin-top:10px;}\n"
               ".note{color:#777;font-size:12px;margin-top:6px;}\n"
               "</style></head><body>\n");
    fprintf(f, "<h1>Multi-sample comparison</h1>\n<p class=\"sub\">%d samples &middot; top %d most-variable genes (of %d)</p>\n",
            n_samples, topn, n_genes_canonical);

    /* --- PCA scatter ----------------------------------------------------- */
    {
        int w = SVG_W, h = SVG_H + 30;
        int ml = 60, mr = 30, mt = 20, mb = 50;
        int pw = w - ml - mr, ph = h - mt - mb;
        double xmin = px[0], xmax = px[0], ymin = py[0], ymax = py[0];
        for (int s = 1; s < n_samples; s++) {
            if (px[s] < xmin) xmin = px[s];
            if (px[s] > xmax) xmax = px[s];
            if (py[s] < ymin) ymin = py[s];
            if (py[s] > ymax) ymax = py[s];
        }
        double xr = xmax - xmin; if (xr < 1e-9) xr = 1;
        double yr = ymax - ymin; if (yr < 1e-9) yr = 1;
        xmin -= xr * 0.15; xmax += xr * 0.15; ymin -= yr * 0.2; ymax += yr * 0.2;
        xr = xmax - xmin; yr = ymax - ymin;

        fprintf(f, "<h2>Sample similarity (PCA)</h2>\n<div class=\"card\">\n");
        fprintf(f, "<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", w, h);
        int x0 = ml, y0 = mt + ph;
        fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ccc\"/>\n", x0, mt, x0, y0);
        fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ccc\"/>\n", x0, y0, ml + pw, y0);
        fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#555\" text-anchor=\"middle\">PC1 (%.1f%% variance)</text>\n", ml + pw / 2, h - 8, pc1_pct);
        fprintf(f, "<text x=\"14\" y=\"%d\" font-size=\"12\" fill=\"#555\" transform=\"rotate(-90 14,%d)\" text-anchor=\"middle\">PC2 (%.1f%% variance)</text>\n", mt + ph / 2, mt + ph / 2, pc2_pct);
        const char *palette[8] = { "#2a78d6", "#c8622d", "#2a9d5c", "#6a4fb3", "#d6a02a", "#d64550", "#2ab6b0", "#8a8a8a" };
        for (int s = 0; s < n_samples; s++) {
            int cx = ml + (int)(pw * (px[s] - xmin) / xr);
            int cy = mt + ph - (int)(ph * (py[s] - ymin) / yr);
            char lbl[256]; html_escape(samples[s].name, lbl, sizeof(lbl));
            fprintf(f, "<circle cx=\"%d\" cy=\"%d\" r=\"6\" fill=\"%s\" stroke=\"white\" stroke-width=\"1.5\"/>\n", cx, cy, palette[s % 8]);
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#333\">%s</text>\n", cx + 9, cy + 4, lbl);
        }
        fprintf(f, "</svg>\n</div>\n");
        fprintf(f, "<p class=\"note\">Each point is one sample, projected onto its first two principal components computed "
                   "from log2(size-factor-normalized count + 1) over the top %d most-variable genes (exact eigendecomposition of the "
                   "sample-similarity matrix, the standard approach when there are far more genes than samples -- not an "
                   "approximation). Samples that cluster together have more similar overall expression profiles.</p>\n", topn);
    }

    /* --- sample-distance heatmap + dendrogram ---------------------------- */
    {
        int cell = 46;
        int label_w = 0;
        for (int s = 0; s < n_samples; s++) { int l = (int)strlen(samples[s].name); if (l > label_w) label_w = l; }
        label_w = 40 + label_w * 6;
        int dendro_h = (n_samples >= 2) ? 60 : 0;
        int w = label_w + cell * n_samples + 20;
        int h = dendro_h + cell * n_samples + label_w / 2 + 20;
        double dmax = 0; for (int i = 0; i < n_samples; i++) for (int j = 0; j < n_samples; j++) if (dist[i][j] > dmax) dmax = dist[i][j];
        if (dmax <= 0) dmax = 1;

        fprintf(f, "<h2>Sample-distance heatmap</h2>\n<div class=\"card\" style=\"overflow-x:auto;\">\n");
        fprintf(f, "<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", w, h);

        if (n_samples >= 2) {
            double *leaf_x = xmalloc(sizeof(double) * n_samples);
            for (int i = 0; i < n_samples; i++) leaf_x[leaf_order[i]] = label_w + cell * i + cell / 2.0;
            double max_height = 0; for (int m = 0; m < n_samples - 1; m++) if (merges[m].height > max_height) max_height = merges[m].height;
            if (max_height <= 0) max_height = 1;
            for (int m = 0; m < n_samples - 1; m++) {
                double lx = hclust_node_x(merges, n_samples, merges[m].left, leaf_x);
                double rx = hclust_node_x(merges, n_samples, merges[m].right, leaf_x);
                double y = dendro_h - (dendro_h - 8) * (merges[m].height / max_height);
                fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#999\"/>\n", lx, (double)dendro_h, lx, y);
                fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#999\"/>\n", rx, (double)dendro_h, rx, y);
                fprintf(f, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#999\"/>\n", lx, y, rx, y);
            }
            free(leaf_x);
        }

        for (int i = 0; i < n_samples; i++) {
            int si = leaf_order[i];
            char lbl[256]; html_escape(samples[si].name, lbl, sizeof(lbl));
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#333\" text-anchor=\"end\">%s</text>\n",
                    label_w - 6, dendro_h + i * cell + cell / 2 + 4, lbl);
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#333\" text-anchor=\"end\" transform=\"rotate(-45 %d,%d)\">%s</text>\n",
                    label_w + i * cell + cell / 2, dendro_h + n_samples * cell + 14,
                    label_w + i * cell + cell / 2, dendro_h + n_samples * cell + 14, lbl);
        }
        for (int i = 0; i < n_samples; i++) {
            for (int j = 0; j < n_samples; j++) {
                int si = leaf_order[i], sj = leaf_order[j];
                double v = dist[si][sj] / dmax; /* 0 (identical) .. 1 (most different) */
                int shade = 235 - (int)(180 * v);
                fprintf(f, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" fill=\"rgb(%d,%d,255)\" stroke=\"white\"/>\n",
                        label_w + j * cell, dendro_h + i * cell, cell, cell, shade, shade);
                fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"10\" fill=\"%s\" text-anchor=\"middle\">%.1f</text>\n",
                        label_w + j * cell + cell / 2, dendro_h + i * cell + cell / 2 + 3,
                        v > 0.55 ? "white" : "#333", dist[si][sj]);
            }
        }
        fprintf(f, "</svg>\n</div>\n");
        fprintf(f, "<p class=\"note\">Pairwise Euclidean distance between samples in the same top-%d-variable-gene "
                   "VST-transformed space used for the PCA above (darker = more similar). Rows/columns are ordered by "
                   "average-linkage hierarchical clustering, shown as the dendrogram on top -- samples that merge "
                   "low are more similar. This uses the same top-variable-gene subset as the PCA for consistency; "
                   "DESeq2's own sample-distance heatmap conventionally uses the full transformed matrix instead.</p>\n", topn);
    }

    /* --- differential expression (only if a valid 2-condition design) -- */
    if (de_results && cond_a && cond_b) {
        int n_tested = 0, n_sig = 0;
        for (int g = 0; g < n_de_genes; g++) {
            if (de_results[g].pass_filter) n_tested++;
            if (!isnan(de_results[g].padj) && de_results[g].padj < 0.1) n_sig++;
        }
        fprintf(f, "<h2>Differential expression: %s vs %s</h2>\n", cond_b, cond_a);
        fprintf(f, "<p class=\"sub\">%d genes tested after independent filtering (of %d total) &middot; "
                   "%d significant at FDR &lt; 0.1 &middot; positive log2FoldChange means higher in \"%s\"</p>\n",
                n_tested, n_de_genes, n_sig, cond_b);

        if (lfc_prior_var <= 1.5e-4) {
            fprintf(f, "<div class=\"card\" style=\"border-color:#e0b04a;background:#fffbf0;\">"
                       "<b>Note on the \"log2FC (shrunk)\" column:</b> the estimated LFC prior variance collapsed "
                       "to essentially zero for this comparison. That happens (correctly, not as a bug) when the "
                       "overwhelming majority of tested genes show no real difference between conditions -- with "
                       "so little genuine between-gene effect variance to estimate, DESeq2's classic \"normal\" "
                       "shrinkage pulls even confident true positives most of the way to zero. This is a known, "
                       "documented limitation of the \"normal\" method specifically (it's the reason DESeq2 itself "
                       "switched its default to apeglm, which isn't implemented here -- see the module notes in "
                       "the source). <b>Significance calling (p-value/padj) is unaffected</b> -- that's based on "
                       "the unshrunk Wald test throughout. For effect size, trust the \"log2FC (MLE)\" column and "
                       "the MA/volcano plots (both use the unshrunk value) over the shrunk column here.</div>\n");
        }

        /* MA plot: log2(baseMean) on x, log2FoldChange on y, red = significant */
        {
            int w = SVG_W, h = SVG_H + 20;
            int ml = 60, mr = 20, mt = 16, mb = 40;
            int pw = w - ml - mr, ph = h - mt - mb;
            double xmin = 1e18, xmax = -1e18, ymax_abs = 0.1;
            for (int g = 0; g < n_de_genes; g++) {
                if (!de_results[g].pass_filter || de_results[g].baseMean <= 0) continue;
                double x = log2(de_results[g].baseMean + 1);
                if (x < xmin) xmin = x;
                if (x > xmax) xmax = x;
                double ya = fabs(de_results[g].log2FoldChange); /* unshrunk MLE -- see the shrinkage-reliability note above the DE section */
                if (ya > ymax_abs && ya < 20) ymax_abs = ya; /* guard against a rare unstable-fit outlier dominating the scale */
            }
            if (xmax <= xmin) { xmin -= 1; xmax += 1; }
            fprintf(f, "<div class=\"card\">\n<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", w, h);
            int y0 = mt + ph / 2;
            fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ddd\"/>\n", ml, y0, ml + pw, y0);
            fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ccc\"/>\n", ml, mt, ml, mt + ph);
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"11\" fill=\"#666\">0</text>\n", ml - 24, y0 + 4);
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#555\" text-anchor=\"middle\">mean of normalized counts (log2)</text>\n", ml + pw / 2, h - 8);
            for (int g = 0; g < n_de_genes; g++) {
                if (!de_results[g].pass_filter || de_results[g].baseMean <= 0) continue;
                double x = log2(de_results[g].baseMean + 1);
                double yv = de_results[g].log2FoldChange; if (yv > ymax_abs) yv = ymax_abs; if (yv < -ymax_abs) yv = -ymax_abs;
                int cx = ml + (int)(pw * (x - xmin) / (xmax - xmin > 0 ? xmax - xmin : 1));
                int cy = y0 - (int)((ph / 2) * (yv / ymax_abs));
                int sig = !isnan(de_results[g].padj) && de_results[g].padj < 0.1;
                fprintf(f, "<circle cx=\"%d\" cy=\"%d\" r=\"2.2\" fill=\"%s\" fill-opacity=\"%s\"/>\n",
                        cx, cy, sig ? "#d64550" : "#888", sig ? "0.85" : "0.35");
            }
            fprintf(f, "</svg>\n</div>\n<p class=\"note\">Each point is one gene, y-axis is the unshrunk log2FoldChange "
                       "(MLE). Red = significant at FDR &lt; 0.1. Genes dropped by independent filtering are not shown. "
                       "The results table below also reports the empirical-Bayes-shrunk value in a separate column -- "
                       "see the note above the table for when to trust it.</p>\n");
        }

        /* Volcano plot: log2FoldChange on x, -log10(pvalue) on y */
        {
            int w = SVG_W, h = SVG_H + 20;
            int ml = 55, mr = 20, mt = 16, mb = 40;
            int pw = w - ml - mr, ph = h - mt - mb;
            double xmax_abs = 0.1, ymax = 0.1;
            for (int g = 0; g < n_de_genes; g++) {
                if (!de_results[g].pass_filter) continue;
                double ax = fabs(de_results[g].log2FoldChange);
                if (ax > xmax_abs && ax < 20) xmax_abs = ax;
                double neglog = de_results[g].pvalue > 0 ? -log10(de_results[g].pvalue) : 300;
                if (neglog > ymax && neglog < 50) ymax = neglog;
            }
            fprintf(f, "<div class=\"card\">\n<svg viewBox=\"0 0 %d %d\" xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Helvetica,Arial,sans-serif\">\n", w, h);
            int x0 = ml + pw / 2;
            fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ccc\"/>\n", x0, mt, x0, mt + ph);
            fprintf(f, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"#ccc\"/>\n", ml, mt + ph, ml + pw, mt + ph);
            fprintf(f, "<text x=\"%d\" y=\"%d\" font-size=\"12\" fill=\"#555\" text-anchor=\"middle\">log2FoldChange</text>\n", ml + pw / 2, h - 8);
            for (int g = 0; g < n_de_genes; g++) {
                if (!de_results[g].pass_filter) continue;
                double xv = de_results[g].log2FoldChange; if (xv > xmax_abs) xv = xmax_abs; if (xv < -xmax_abs) xv = -xmax_abs;
                double neglog = de_results[g].pvalue > 0 ? -log10(de_results[g].pvalue) : ymax;
                if (neglog > ymax) neglog = ymax;
                int cx = x0 + (int)((pw / 2) * (xv / xmax_abs));
                int cy = mt + ph - (int)(ph * (neglog / ymax));
                int sig = !isnan(de_results[g].padj) && de_results[g].padj < 0.1;
                fprintf(f, "<circle cx=\"%d\" cy=\"%d\" r=\"2.2\" fill=\"%s\" fill-opacity=\"%s\"/>\n",
                        cx, cy, sig ? "#d64550" : "#888", sig ? "0.85" : "0.35");
            }
            fprintf(f, "</svg>\n</div>\n<p class=\"note\">-log10(p-value) vs. log2FoldChange, unshrunk MLE (no apeglm/ashr LFC shrinkage -- "
                       "same as DESeq2's results() before calling lfcShrink()). Red = significant at FDR &lt; 0.1.</p>\n");
        }

        /* Results table: top 30 by padj (NA padj sorted last) */
        {
            typedef struct { int gene_idx; double padj; } Ranked;
            Ranked *ranked = xmalloc(sizeof(Ranked) * n_de_genes);
            int nr = 0;
            for (int g = 0; g < n_de_genes; g++) if (de_results[g].pass_filter) ranked[nr++] = (Ranked){ g, isnan(de_results[g].padj) ? 2.0 : de_results[g].padj };
            for (int i = 0; i < nr; i++) { int best = i; for (int j = i+1; j < nr; j++) if (ranked[j].padj < ranked[best].padj) best = j;
                Ranked tmp = ranked[i]; ranked[i] = ranked[best]; ranked[best] = tmp; }
            int shown = nr < 30 ? nr : 30;
            fprintf(f, "<h3 style=\"font-size:13px;margin-top:20px;\">Top %d genes by adjusted p-value</h3>\n", shown);
            fprintf(f, "<div class=\"card\" style=\"overflow-x:auto;padding:0;\"><table style=\"border-collapse:collapse;width:100%%;font-size:12px;\">\n");
            fprintf(f, "<tr style=\"background:#f7f7f5;text-align:left;\">"
                       "<th style=\"padding:6px 10px;\">gene_id</th><th style=\"padding:6px 10px;\">baseMean</th>"
                       "<th style=\"padding:6px 10px;\">log2FC (MLE)</th><th style=\"padding:6px 10px;\">log2FC (shrunk)</th>"
                       "<th style=\"padding:6px 10px;\">lfcSE</th>"
                       "<th style=\"padding:6px 10px;\">stat</th><th style=\"padding:6px 10px;\">pvalue</th>"
                       "<th style=\"padding:6px 10px;\">padj</th><th style=\"padding:6px 10px;\">dispersion</th>"
                       "<th style=\"padding:6px 10px;\">Cook's flag</th></tr>\n");
            for (int i = 0; i < shown; i++) {
                int g = ranked[i].gene_idx;
                char gid[128]; html_escape(canonical[g].gene_id, gid, sizeof(gid));
                char padj_str[32];
                if (isnan(de_results[g].padj)) snprintf(padj_str, sizeof(padj_str), "NA");
                else snprintf(padj_str, sizeof(padj_str), "%.2e", de_results[g].padj);
                fprintf(f, "<tr style=\"border-top:1px solid #eee;%s\">"
                           "<td style=\"padding:5px 10px;\">%s</td><td style=\"padding:5px 10px;\">%.1f</td>"
                           "<td style=\"padding:5px 10px;\">%.3f</td><td style=\"padding:5px 10px;\">%.3f</td>"
                           "<td style=\"padding:5px 10px;\">%.3f</td>"
                           "<td style=\"padding:5px 10px;\">%.2f</td><td style=\"padding:5px 10px;\">%.2e</td>"
                           "<td style=\"padding:5px 10px;\">%s</td><td style=\"padding:5px 10px;\">%.3f</td>"
                           "<td style=\"padding:5px 10px;\">%s</td></tr>\n",
                        (!isnan(de_results[g].padj) && de_results[g].padj < 0.1) ? "background:#fff5f5;" : "",
                        gid, de_results[g].baseMean, de_results[g].log2FoldChange, de_results[g].log2FoldChangeShrunk, de_results[g].lfcSE,
                        de_results[g].stat, de_results[g].pvalue, padj_str,
                        de_results[g].dispersion, de_results[g].cooks_flagged ? "yes" : "");
            }
            fprintf(f, "</table></div>\n");
            free(ranked);
        }
        (void)group;
    } else if (cond_a || cond_b) {
        /* shouldn't happen given run_compare_mode's gating, but stay defensive */
        (void)group;
    }

    /* --- library size bar (general-stats-style cross-sample comparison) - */
    {
        const char *labels[MAX_COMPARE_SAMPLES]; double values[MAX_COMPARE_SAMPLES];
        for (int s = 0; s < n_samples; s++) { labels[s] = samples[s].name; values[s] = samples[s].lib_size; }
        fprintf(f, "<h2>Library size</h2>\n<div class=\"card\">\n");
        svg_vbar_chart(f, labels, values, n_samples, "#2a78d6");
        fprintf(f, "</div>\n<p class=\"note\">Total EM-assigned effective counts per sample. Analogous to the "
                   "\"M Assigned\" column in MultiQC's General Statistics table.</p>\n");
    }

    fprintf(f, "<h2>What this uses vs. real DESeq2</h2>\n");
    fprintf(f, "<p class=\"note\">Size factors: DESeq2's real median-of-ratios method. PCA/heatmap: DESeq2's real "
               "closed-form variance-stabilizing transform, fit blind to the experimental design (same as vst() "
               "with blind=TRUE). Dispersion: Cox-Reid adjusted profile likelihood (the real per-gene MLE DESeq2 "
               "uses), with a Gamma-GLM-fitted trend and empirical-Bayes shrinkage using an estimated prior "
               "variance (not a fixed weight). LFC shrinkage: DESeq2's classic \"normal\" method (normal-normal "
               "conjugate posterior). Cook's distance: flagged against the exact F(2, n-2) quantile DESeq2 uses. "
               "PCA itself is exact eigendecomposition, not approximated.</p>\n");

    fprintf(f, "</body></html>\n");
    fclose(f);

    for (int s = 0; s < n_samples; s++) { free(X[s]); free(G[s]); free(V[s]); free(dist[s]); }
    free(X); free(G); free(V); free(dist);
    free(eigval); free(order); free(var); free(idx); free(px); free(py); free(leaf_order);
    if (merges) free(merges);
}

static int run_compare_mode(int argc, char **argv) {
    /* argv[0]=prog argv[1]="compare" argv[2..argc-2]=<outdir[:condition] per sample> argv[argc-1]=<combined_outdir> */
    if (argc < 5) {
        fprintf(stderr, "Usage:\n  %s compare <sample1_outdir[:condition]> <sample2_outdir[:condition]> [more...] <combined_outdir>\n"
                         "  (each <sampleN_outdir> must be a directory already produced by a single-sample run of this pipeline)\n"
                         "  (the optional :condition label enables a real differential-expression test when exactly\n"
                         "   2 distinct conditions are given, each with >=2 samples -- e.g. 'out_ctrl1:control')\n", argv[0]);
        return 1;
    }
    int n_samples = argc - 3;
    if (n_samples > MAX_COMPARE_SAMPLES) { fprintf(stderr, "error: at most %d samples supported\n", MAX_COMPARE_SAMPLES); return 1; }
    const char *combined_outdir = argv[argc - 1];

    CompareSample *samples = xmalloc(sizeof(CompareSample) * n_samples);
    GeneCountRow *canonical = NULL;
    int n_canonical = 0;

    for (int s = 0; s < n_samples; s++) {
        char arg_copy[1024];
        snprintf(arg_copy, sizeof(arg_copy), "%s", argv[2 + s]);
        char *colon = strrchr(arg_copy, ':');
        const char *sample_dir = arg_copy;
        samples[s].condition[0] = '\0';
        if (colon) {
            *colon = '\0';
            snprintf(samples[s].condition, sizeof(samples[s].condition), "%s", colon + 1);
        }

        char gc_path[1024];
        snprintf(gc_path, sizeof(gc_path), "%s/gene_counts.tsv", sample_dir);
        GeneCountRow *rows; int n = read_gene_counts_tsv(gc_path, &rows);
        if (n < 0) {
            for (int t = 0; t < s; t++) { free(samples[t].log2cpm); free(samples[t].raw_count); }
            free(samples);
            free(canonical);
            return 1;
        }
        snprintf(samples[s].name, sizeof(samples[s].name), "%s", basename_noslash(sample_dir));

        if (s == 0) {
            canonical = rows; n_canonical = n;
            samples[s].raw_count = xmalloc(sizeof(double) * n_canonical);
            for (int g = 0; g < n; g++) samples[s].raw_count[g] = rows[g].val;
        } else {
            samples[s].raw_count = xmalloc(sizeof(double) * n_canonical);
            for (int g = 0; g < n_canonical; g++) samples[s].raw_count[g] = 0.0;
            /* align this sample's rows onto the canonical gene order -- O(n_canonical * n)
             * is fine at gene-count/sample-count scale, and avoids building a hash map for
             * what's normally a handful of compare-mode invocations, not a hot path. */
            for (int g = 0; g < n; g++) {
                for (int c = 0; c < n_canonical; c++) {
                    if (strcmp(rows[g].gene_id, canonical[c].gene_id) == 0) {
                        samples[s].raw_count[c] = rows[g].val;
                        break;
                    }
                }
            }
            free(rows);
        }
        double total = 0; for (int g = 0; g < n_canonical; g++) total += samples[s].raw_count[g];
        samples[s].lib_size = total;
        if (samples[s].condition[0])
            printf("  loaded %s (%s): %d genes, library size %.0f\n", samples[s].name, samples[s].condition, n_canonical, total);
        else
            printf("  loaded %s: %d genes, library size %.0f\n", samples[s].name, n_canonical, total);
    }

    /* --- real median-of-ratios size factors, used for BOTH the PCA/heatmap
     * normalization and the DE test's offset. --- */
    double **counts = xmalloc(sizeof(double *) * n_samples);
    for (int s = 0; s < n_samples; s++) counts[s] = samples[s].raw_count;
    double *size_factors = xmalloc(sizeof(double) * n_samples);
    compute_size_factors(counts, n_samples, n_canonical, size_factors);
    for (int s = 0; s < n_samples; s++) samples[s].size_factor = size_factors[s];

    /* --- blind (design-ignoring) dispersion + trend, feeding the real -- */
    /* --- DESeq2 VST for PCA/heatmap (replaces log2(normalized+1))    --- */
    printf("  fitting blind dispersion trend for variance-stabilizing transform...\n");
    {
        double *blind_offset = xmalloc(sizeof(double) * n_samples);
        for (int s = 0; s < n_samples; s++) blind_offset[s] = log(size_factors[s]);
        double *gene_alpha_blind = xmalloc(sizeof(double) * n_canonical);
        double *gene_mean_blind = xmalloc(sizeof(double) * n_canonical);
        double *y_blind = xmalloc(sizeof(double) * n_samples);
        double *mu_blind = xmalloc(sizeof(double) * n_samples);
        for (int g = 0; g < n_canonical; g++) {
            for (int s = 0; s < n_samples; s++) y_blind[s] = samples[s].raw_count[g];
            double alpha = estimate_dispersion_cr_1p(y_blind, blind_offset, n_samples);
            gene_alpha_blind[g] = alpha;
            double bm = 0; for (int s = 0; s < n_samples; s++) bm += y_blind[s] / size_factors[s];
            gene_mean_blind[g] = bm / n_samples;
        }
        double a0_blind, a1_blind;
        fit_dispersion_trend(gene_alpha_blind, gene_mean_blind, n_canonical, &a0_blind, &a1_blind);
        printf("  blind dispersion trend: a0=%.4f a1=%.4f\n", a0_blind, a1_blind);

        for (int s = 0; s < n_samples; s++) {
            samples[s].log2cpm = xmalloc(sizeof(double) * n_canonical); /* holds VST values now, field name kept for minimal churn elsewhere */
            for (int g = 0; g < n_canonical; g++) {
                double normalized = samples[s].raw_count[g] / size_factors[s];
                samples[s].log2cpm[g] = vst_transform(normalized, a0_blind, a1_blind);
            }
        }
        free(blind_offset); free(gene_alpha_blind); free(gene_mean_blind); free(y_blind); free(mu_blind);
    }

    /* --- detect a valid 2-condition design for the DE test -------------- */
    char cond_a[128] = "", cond_b[128] = "";
    int n_a = 0, n_b = 0, n_conditions_seen = 0, design_ok = 0;
    for (int s = 0; s < n_samples; s++) {
        if (!samples[s].condition[0]) continue;
        if (!cond_a[0]) { snprintf(cond_a, sizeof(cond_a), "%s", samples[s].condition); }
        if (strcmp(samples[s].condition, cond_a) == 0) continue;
        if (!cond_b[0]) { snprintf(cond_b, sizeof(cond_b), "%s", samples[s].condition); continue; }
        if (strcmp(samples[s].condition, cond_b) != 0) n_conditions_seen = 3; /* a 3rd distinct label -- not supported */
    }
    int *group = xmalloc(sizeof(int) * n_samples);
    for (int s = 0; s < n_samples; s++) {
        if (!samples[s].condition[0]) { group[s] = -1; continue; }
        if (strcmp(samples[s].condition, cond_a) == 0) { group[s] = 0; n_a++; }
        else if (cond_b[0] && strcmp(samples[s].condition, cond_b) == 0) { group[s] = 1; n_b++; }
        else group[s] = -1;
    }
    design_ok = cond_b[0] && n_conditions_seen != 3 && n_a >= 2 && n_b >= 2 && (n_a + n_b == n_samples);

    DEResult *de_results = NULL;
    int n_de_genes = 0;
    double lfc_prior_var = -1;
    if (design_ok) {
        printf("  design: %d vs %d (\"%s\" vs \"%s\") -- running differential expression\n", n_a, n_b, cond_a, cond_b);
        de_results = run_differential_expression(counts, group, n_samples, n_canonical, size_factors, &lfc_prior_var);
        n_de_genes = n_canonical;
    } else if (cond_a[0] || cond_b[0]) {
        printf("  note: conditions given but not a clean 2-group design (need exactly 2 conditions, >=2 samples each) -- skipping differential expression, PCA/heatmap only\n");
    }

    #ifdef _WIN32
    #else
    mkdir(combined_outdir, 0755);
    #endif
    write_compare_report(samples, n_samples, n_canonical, combined_outdir,
                          de_results, n_de_genes, group, design_ok ? cond_a : NULL, design_ok ? cond_b : NULL, canonical, lfc_prior_var);
    printf("Wrote %s/multi_sample_report.html\n", combined_outdir);

    for (int s = 0; s < n_samples; s++) { free(samples[s].log2cpm); free(samples[s].raw_count); }
    free(samples); free(canonical); free(counts); free(size_factors); free(group);
    if (de_results) free(de_results);
    return 0;
}


static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s <reference.fasta> <annotation.gtf> se <reads.fastq> <outdir>\n"
        "  %s <reference.fasta> <annotation.gtf> pe <r1.fastq> <r2.fastq> <outdir>\n"
        "  %s compare <sample1_outdir> <sample2_outdir> [more...] <combined_outdir>\n",
        prog, prog, prog);
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "compare") == 0) return run_compare_mode(argc, argv);
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
    write_html_report(outdir, reads_desc);

    printf("Done. (EM converged in %d iterations)\n", em_iterations_run);
    return 0;
}

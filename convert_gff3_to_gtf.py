import gzip, re, sys

# Ensembl-style GFF3 top-level gene loci. 'gene' covers protein-coding genes;
# everything else that also carries its own gene_id (rRNA/tRNA/snoRNA/snRNA
# loci under ncRNA_gene, plus pseudogenes and TE genes) is a real, separate
# gene the previous version of this script silently dropped -- which meant
# rRNA genes specifically (needed for rRNA-content QC) were never in the
# GTF this pipeline loads at all. All of these feature types have exactly
# the same GFF3 shape (a top-level row with gene_id + biotype), so one
# converter handles all of them.
GENE_FEATURE_TYPES = {'gene', 'ncRNA_gene', 'pseudogene', 'transposable_element_gene'}

# Second-level feature types that sit between a gene and its exons in
# Ensembl-style GFF3 (mRNA for protein-coding genes; tRNA/rRNA/ncRNA/snRNA/
# snoRNA for the corresponding ncRNA_gene loci; pseudogenic_transcript for
# pseudogenes). Each carries Parent=gene:<gene_id> and ID=transcript:<tx_id>
# (or an equivalent without the "gene:"/"transcript:" prefix on some
# non-Ensembl GFF3 sources -- both are handled below).
TRANSCRIPT_FEATURE_TYPES = {
    'mRNA', 'tRNA', 'rRNA', 'ncRNA', 'snRNA', 'snoRNA',
    'pseudogenic_transcript', 'transcript',
}

def _strip_prefix(value):
    """Ensembl GFF3 IDs are namespaced ('gene:YAL069W', 'transcript:YAL069W_mRNA');
    plenty of other GFF3 sources use the bare ID with no namespace at all.
    Handle both by taking whatever follows the last colon, if any."""
    return value.rsplit(':', 1)[-1]

def convert(gff3_path, gtf_path):
    gene_re = re.compile(r'gene_id=([^;]+)')
    biotype_re = re.compile(r'biotype=([^;]+)')
    id_re = re.compile(r'(?:^|;)ID=([^;]+)')
    parent_re = re.compile(r'(?:^|;)Parent=([^;]+)')
    opener = gzip.open if gff3_path.endswith('.gz') else open

    # Pass 1: build transcript_id -> gene_id, from every second-level
    # feature row (mRNA/tRNA/etc.) that links a transcript back to its
    # parent gene. A separate pass rather than deferred single-pass
    # resolution because GFF3 doesn't strictly guarantee parents appear
    # before children in file order, and a second full read is cheap next
    # to the one-time cost of building the k-mer index downstream anyway.
    transcript_to_gene = {}
    with opener(gff3_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] not in TRANSCRIPT_FEATURE_TYPES:
                continue
            id_m = id_re.search(p[8])
            parent_m = parent_re.search(p[8])
            if not id_m or not parent_m:
                continue
            transcript_to_gene[_strip_prefix(id_m.group(1))] = _strip_prefix(parent_m.group(1))

    with opener(gff3_path, 'rt') as f, open(gtf_path, 'w') as out:
        n = 0
        by_biotype = {}
        n_exons = 0
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue

            if p[2] in GENE_FEATURE_TYPES:
                m = gene_re.search(p[8])
                if not m:
                    continue
                bt_m = biotype_re.search(p[8])
                biotype = bt_m.group(1) if bt_m else 'unknown'
                out.write(f'{p[0]}\t{p[1]}\tgene\t{p[3]}\t{p[4]}\t.\t{p[6]}\t.\t'
                          f'gene_id "{m.group(1)}"; gene_biotype "{biotype}";\n')
                n += 1
                by_biotype[biotype] = by_biotype.get(biotype, 0) + 1

            elif p[2] == 'exon':
                # Needed for junction-annotation QC (known-vs-novel splice
                # sites): the pipeline derives known exon-exon junctions
                # from consecutive same-transcript exon rows, so without
                # these every observed splice would trivially classify as
                # novel regardless of whether it's real, annotated
                # splicing. Emitting exon rows changes nothing about gene
                # quantification, which has always used gene-level (not
                # exon-level) intervals -- see load_gtf() in
                # rnaseq_pipeline.c.
                parent_m = parent_re.search(p[8])
                if not parent_m:
                    continue
                tx_id = _strip_prefix(parent_m.group(1))
                gene_id = transcript_to_gene.get(tx_id)
                if not gene_id:
                    continue  # exon's parent transcript has no resolvable parent gene -- skip
                out.write(f'{p[0]}\t{p[1]}\texon\t{p[3]}\t{p[4]}\t.\t{p[6]}\t.\t'
                          f'gene_id "{gene_id}"; transcript_id "{tx_id}";\n')
                n_exons += 1

    print(f'Wrote {n} genes and {n_exons} exon rows to {gtf_path}')
    for bt, count in sorted(by_biotype.items(), key=lambda kv: -kv[1]):
        print(f'  {bt}: {count}')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python3 {sys.argv[0]} <input.gff3[.gz]> <output.gtf>')
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])

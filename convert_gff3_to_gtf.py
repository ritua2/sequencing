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

def convert(gff3_path, gtf_path):
    gene_re = re.compile(r'gene_id=([^;]+)')
    biotype_re = re.compile(r'biotype=([^;]+)')
    opener = gzip.open if gff3_path.endswith('.gz') else open
    with opener(gff3_path, 'rt') as f, open(gtf_path, 'w') as out:
        n = 0
        by_biotype = {}
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] not in GENE_FEATURE_TYPES:
                continue
            m = gene_re.search(p[8])
            if not m:
                continue
            bt_m = biotype_re.search(p[8])
            biotype = bt_m.group(1) if bt_m else 'unknown'
            out.write(f'{p[0]}\t{p[1]}\tgene\t{p[3]}\t{p[4]}\t.\t{p[6]}\t.\t'
                      f'gene_id "{m.group(1)}"; gene_biotype "{biotype}";\n')
            n += 1
            by_biotype[biotype] = by_biotype.get(biotype, 0) + 1
    print(f'Wrote {n} genes to {gtf_path}')
    for bt, count in sorted(by_biotype.items(), key=lambda kv: -kv[1]):
        print(f'  {bt}: {count}')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python3 {sys.argv[0]} <input.gff3[.gz]> <output.gtf>')
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])

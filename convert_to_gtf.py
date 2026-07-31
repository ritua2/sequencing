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

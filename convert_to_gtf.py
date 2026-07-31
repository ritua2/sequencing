import gzip, re, sys

def convert(gff3_path, gtf_path):
    gene_re = re.compile(r'gene_id=([^;]+)')
    with gzip.open(gff3_path, 'rt') as f, open(gtf_path, 'w') as out:
        n = 0
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            m = gene_re.search(p[8])
            if not m:
                continue
            out.write(f'{p[0]}\t{p[1]}\tgene\t{p[3]}\t{p[4]}\t.\t{p[6]}\t.\tgene_id "{m.group(1)}";\n')
            n += 1
    print(f'Wrote {n} genes to {gtf_path}')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python3 {sys.argv[0]} <input.gff3.gz> <output.gtf>')
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])

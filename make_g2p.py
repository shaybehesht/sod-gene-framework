"""Cache the HPO gene-phenotype file as a pandas pickle after checking its MD5."""
import hashlib, pandas as pd
P = 'data/genes_to_phenotype.txt'
MD5 = '6dac74395b2a9192dc0e71594093594f'   # file downloaded November 2, 2026
md5 = hashlib.md5(open(P, 'rb').read()).hexdigest()
assert md5 == MD5, f'genes_to_phenotype.txt MD5 is {md5}, expected {MD5}'
g = pd.read_csv(P, sep='\t', dtype=str)
g.to_pickle('data/g2p.pkl')
print('cached', g.shape, 'MD5 ok')

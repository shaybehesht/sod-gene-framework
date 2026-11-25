"""GO Biological Process over-representation test.

GO gene sets (data/go_sets_<tag>.csv.gz: one row per GO term and gene) and backgrounds
(data/<tag>_background.csv) are exported by R/export_go_annotations.R. Terms are already
restricted to 10-500 background genes, with ancestor terms included.
Test: one-sided hypergeometric; Benjamini-Hochberg FDR across all tested terms.
"""
import pandas as pd
from scipy.stats import hypergeom, false_discovery_control

D = 'data/'
NAMES = dict(pd.read_csv(D + 'go_term_names.csv').values)
_g = pd.read_pickle(D + 'g2p.pkl')
SYM = _g.drop_duplicates('ncbi_gene_id').set_index('ncbi_gene_id').gene_symbol
_cache = {}


def load_sets(tag):
    if tag not in _cache:
        s = pd.read_csv(D + f'go_sets_{tag}.csv.gz', dtype=str)
        _cache[tag] = {go: set(g) for go, g in s.groupby('GO_ID').ncbi_gene_id}
    return _cache[tag]


def go_ora(genes, tag='nonIEA_HPO'):
    """Return (results table, n query genes tested, background size, query genes not in background)."""
    sets = load_sets(tag)
    bg = set(pd.read_csv(D + f'{tag}_background.csv', dtype=str).ncbi_gene_id)
    q = set(genes) & bg
    N, n = len(bg), len(q)
    if n == 0:
        raise ValueError('no query genes in the background')
    rows = []
    for go, s in sets.items():
        ov = sorted(q & s); k, K = len(ov), len(s)
        rows.append((go, NAMES.get(go, ''), k, K, (k / n) / (K / N), hypergeom.sf(k - 1, N, K, n),
                     ';'.join(SYM.get(x, x) for x in ov)))
    df = pd.DataFrame(rows, columns=['GO_ID', 'Description', 'Count', 'BgTermSize', 'FoldEnrichment', 'pvalue', 'genes'])
    df['p_adjust_BH'] = false_discovery_control(df.pvalue.values, method='bh')
    return df.sort_values(['p_adjust_BH', 'pvalue']), n, N, sorted(set(genes) - bg)

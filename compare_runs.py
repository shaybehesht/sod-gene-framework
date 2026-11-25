"""Write the HPO term table and the Figure 2 change table by comparing the
original and corrected runs (outputs/out_orig vs outputs/out_corr)."""
import re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from hpo import load_obo

# HPO terms used, with gene counts before and after the zero-frequency exclusion
t = pd.read_csv('TableS3_curated_HPO_terms.csv')
g = pd.read_pickle('data/g2p.pkl')
def zf(f):
    f = str(f).strip(); m = re.match(r'^(\d+)/(\d+)$', f)
    if m: return int(m.group(1)) == 0
    m = re.match(r'^([\d.]+)%$', f); return bool(m) and float(m.group(1)) == 0
g['z'] = g.frequency.map(zf)
T, C, ALT = load_obo('data/hp.obo'); g['hpo_id'] = g.hpo_id.map(lambda x: ALT.get(x, x))
n_all = g.groupby('hpo_id').ncbi_gene_id.nunique(); n_nz = g[~g.z].groupby('hpo_id').ncbi_gene_id.nunique()
rows = []
for r in t.itertuples():
    sod = r.hpo_id == 'HP:0100842'
    if sod and r.domain == 'optic': continue
    feat = {'midline': 'Midline defect', 'optic': 'ONH or optic nerve aplasia', 'pituitary': 'Pituitary involvement'}[r.domain]
    rows.append({'Feature': 'Inclusion flag (separate)' if sod else feat, 'HPO ID': r.hpo_id, 'HPO term': r.hpo_name,
                 'Genes annotated (all entries)': int(n_all.get(r.hpo_id, 0)),
                 'Genes annotated (zero-frequency entries excluded)': int(n_nz.get(r.hpo_id, 0)),
                 'Note': 'No longer counts toward optic or midline' if sod else ''})
pd.DataFrame(rows).to_csv('outputs/curated_terms_used.csv', index=False)

# Figure 2 cell changes
o = pd.read_csv('outputs/out_orig/GO_heatmap_data.csv'); c = pd.read_csv('outputs/out_corr/GO_heatmap_data.csv')
cols = ['Reported SOD genes', 'Full framework', 'All three features', 'ONH + midline', 'Pituitary + midline']
lab = {'Full framework': 'SOD candidate genes'}
def b(f):
    if f >= 0.05: return 'grey'
    v = -np.log10(f)
    for hi, nm in [(2, '1.3–2'), (3, '2–3'), (4, '3–4'), (5, '4–5'), (6, '5–6')]:
        if v < hi: return nm
    return '≥6'
def ncol(run, col):
    d = 'outputs/in_' + run
    if col == 'Reported SOD genes': return len(pd.read_csv(f'{d}/lit.csv'))
    if col == 'Full framework': return len(pd.read_csv(f'{d}/framework_genes.csv'))
    gg = pd.read_csv(f'{d}/groups.csv').group
    return int((gg == {'All three features': 'All three'}.get(col, col)).sum())
out = []
for col in cols:
    a, z = ncol('orig', col), ncol('corr', col)
    if a != z: out.append({'Panel': 'A (column header)', 'GO term': '—', 'Column': lab.get(col, col), 'Previous': f'n = {a}', 'Corrected': f'n = {z}', 'What to change': 'Header n'})
for i in range(len(o)):
    for col in cols:
        fo, fc, no, nc = o.loc[i, col], c.loc[i, col], int(o.loc[i, col + '_n']), int(c.loc[i, col + '_n'])
        bo, bc = b(fo), b(fc); what = []
        if bo == 'grey' and bc != 'grey': what.append('becomes significant (add colour and number)')
        elif bo != 'grey' and bc == 'grey': what.append('no longer significant (grey, remove number)')
        else:
            if bc != 'grey' and no != nc: what.append(f'number {no} → {nc}')
            if bc != 'grey' and bo != bc: what.append(f'colour bin {bo} → {bc}')
        if what:
            out.append({'Panel': 'A', 'GO term': o.loc[i, 'term'], 'Column': lab.get(col, col), 'Previous': f'FDR {fo:.2e}, n = {no}',
                        'Corrected': f'FDR {fc:.2e}, n = {nc}', 'What to change': '; '.join(what)})
mo = pd.read_csv('outputs/out_orig/GO_module_summary.csv').set_index('module'); mc = pd.read_csv('outputs/out_corr/GO_module_summary.csv').set_index('module')
fmt = lambda x: f'{x:.1e}'.replace('e-0', 'e-')
for m in mc.index:
    a, z = mo.loc[m], mc.loc[m]
    if m != 'Dystroglycan O-mannosylation' and (fmt(a.FDR) != fmt(z.FDR) or a.n_exclusive != z.n_exclusive or a.n_total != z.n_total or a.robust != z.robust):
        out.append({'Panel': 'B (bar label)', 'GO term': '—', 'Column': m, 'Previous': f'FDR {fmt(a.FDR)}', 'Corrected': f'FDR {fmt(z.FDR)}', 'What to change': 'label FDR; bars unchanged'})
Fo = pd.read_csv('outputs/out_orig/GO_Full_framework.csv'); Fc = pd.read_csv('outputs/out_corr/GO_Full_framework.csv')
po = dict(zip(Fo.Description, Fo.p_adjust_BH))['protein O-linked glycosylation']; pc = dict(zip(Fc.Description, Fc.p_adjust_BH))['protein O-linked glycosylation']
out.append({'Panel': 'B (bar label)', 'GO term': '—', 'Column': 'Dystroglycan O-mannosylation', 'Previous': f'FDR {fmt(po)} (parent term)',
            'Corrected': f'FDR {fmt(mc.loc["Dystroglycan O-mannosylation"].FDR)} (lowest term, as Methods state) or {fmt(pc)} (parent term, as before)',
            'What to change': 'pre-existing: figure used the parent term; Methods say lowest FDR'})
pd.DataFrame(out).to_csv('outputs/figure2_changes.csv', index=False)
print(pd.DataFrame(out).to_string())

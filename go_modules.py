"""GO enrichment, term clustering and module assignment for the SOD framework.

Uses go_ora() from go_ora.py unchanged (one-sided hypergeometric, BH FDR,
non-IEA GO BP sets restricted to the HPO background, 10-500 genes).

Clustering of significant terms follows Methods: Jaccard distance between the
framework genes in each term, average linkage, cut height 0.8.

Modules are named by their defining terms (MODULES below), taken from the
original analysis. A gene's primary module is the module whose most specific
significant term (fewest background genes) contains it.

Usage: python3 go_modules.py <framework_genes.csv> <groups.csv> <lit.csv> <outdir> [same_disease_genes.csv]
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from go_ora import go_ora

# Module definitions. Jaccard clustering of the 62 significant terms gives nine clusters; the seven
# modules below are defined by the specific GO terms within them. Two clusters share the axis/left-right
# module, and one cluster (limb and appendage development, cell fate commitment, dopaminergic neuron
# differentiation) contains only broad developmental terms and is not a module. Terms not listed here
# are reported as broad terms.
MODULES = {
    'Hedgehog signaling': ['smoothened signaling pathway', 'regulation of smoothened signaling pathway'],
    'Forebrain and pituitary development': ['forebrain development', 'pituitary gland development', 'diencephalon development'],
    'Axis and left-right patterning': ['axis specification', 'anterior/posterior axis specification',
        'embryonic pattern specification', 'anterior/posterior pattern specification', 'embryonic axis specification',
        'heart looping', 'embryonic heart tube morphogenesis', 'determination of heart left/right asymmetry',
        'determination of left/right symmetry', 'left/right pattern formation', 'embryonic heart tube development',
        'specification of symmetry', 'determination of bilateral symmetry', 'nodal signaling pathway'],
    'Dystroglycan O-mannosylation': ['protein O-linked glycosylation via mannose', 'protein O-linked glycosylation'],
    'Axon guidance': ['neuron projection guidance'],
    'Microtubule / chromatid cohesion': ['sister chromatid cohesion', 'microtubule cytoskeleton organization'],
    'mRNA splicing regulation': ['regulation of mRNA splicing, via spliceosome', 'regulation of mRNA processing'],
}
HEAT = [('Forebrain and pituitary development', t) for t in ['pituitary gland development', 'diencephalon development', 'forebrain development', 'endocrine system development']] + \
       [('Hedgehog signaling', t) for t in ['smoothened signaling pathway', 'regulation of smoothened signaling pathway']] + \
       [('Axis and left-right patterning', t) for t in ['axis specification', 'nodal signaling pathway', 'determination of left/right symmetry', 'regionalization']] + \
       [('Dystroglycan O-mannosylation', t) for t in ['protein O-linked glycosylation via mannose', 'protein O-linked glycosylation']] + \
       [('Axon guidance', 'neuron projection guidance')] + \
       [('Microtubule / chromatid cohesion', t) for t in ['sister chromatid cohesion', 'microtubule cytoskeleton organization']] + \
       [('mRNA splicing regulation', 'regulation of mRNA splicing, via spliceosome')]


def cluster_terms(sig, cut=0.8):
    gs = [set(x.split(';')) for x in sig.genes]
    m = len(gs); S = np.ones((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            S[i, j] = S[j, i] = len(gs[i] & gs[j]) / len(gs[i] | gs[j])
    D = 1 - S; np.fill_diagonal(D, 0)
    cl = fcluster(linkage(squareform(D, checks=False), 'average'), cut, 'distance')
    return cl, pd.DataFrame(S, index=sig.Description, columns=sig.Description)


def main(fw_csv, groups_csv, lit_csv, out, same_csv=None):
    os.makedirs(out, exist_ok=True)
    fw = pd.read_csv(fw_csv, dtype=str)            # columns: gene, ncbi_gene_id
    grp = pd.read_csv(groups_csv, dtype=str)       # columns: ncbi_gene_id, group
    lit = pd.read_csv(lit_csv, dtype=str)          # columns: gene, ncbi_gene_id
    sets = {'Reported SOD genes': list(lit.ncbi_gene_id), 'Full framework': list(fw.ncbi_gene_id)}
    for gname in ['All three', 'ONH + midline', 'Pituitary + midline']:
        sets[{'All three': 'All three features'}.get(gname, gname)] = list(grp[grp.group == gname].ncbi_gene_id)
    res, log = {}, {}
    for nm, ids in sets.items():
        df, n, N, miss = go_ora(ids)
        res[nm] = df
        df.to_csv(os.path.join(out, 'GO_' + nm.replace(' ', '_').replace('+', 'and') + '.csv'), index=False)
        log[nm] = dict(n_genes=len(ids), n_tested=n, background=N, n_sig=int((df.p_adjust_BH < 0.05).sum()),
                       not_testable=len(miss))
    sens = {}
    if same_csv:
        s = pd.read_csv(same_csv, dtype=str)
        sens['same'], _, _, _ = go_ora(list(s.ncbi_gene_id))
        sens['same'].to_csv(os.path.join(out, 'GO_Full_framework_same_disease.csv'), index=False)
    sens['allev'], _, _, _ = go_ora(list(fw.ncbi_gene_id), tag='allEvidence_HPO')
    sens['allev'].to_csv(os.path.join(out, 'GO_Full_framework_all_evidence.csv'), index=False)

    F = res['Full framework']
    sig = F[F.p_adjust_BH < 0.05].copy().reset_index(drop=True)
    sig['cluster'], sim = cluster_terms(sig)
    t2m = {t: m for m, ts in MODULES.items() for t in ts}
    sig['module'] = sig.Description.map(t2m).fillna('')
    sig.to_csv(os.path.join(out, 'GO_significant_terms_with_clusters.csv'), index=False)
    sim.to_csv(os.path.join(out, 'GO_term_similarity_jaccard.csv'))

    fdr = dict(zip(F.Description, F.p_adjust_BH)); genes = dict(zip(F.Description, F.genes))
    K = dict(zip(F.Description, F.BgTermSize))
    fsame = dict(zip(sens['same'].Description, sens['same'].p_adjust_BH)) if 'same' in sens else {}
    fall = dict(zip(sens['allev'].Description, sens['allev'].p_adjust_BH))
    mg, rows = {}, []
    for m, ts in MODULES.items():
        live = [t for t in ts if fdr.get(t, 1) < 0.05]
        mg[m] = set().union(*[set(genes[t].split(';')) for t in live]) if live else set()
        lead = min(ts, key=lambda t: fdr.get(t, 1))
        rows.append(dict(module=m, label_term=lead, n_terms_significant=len(live), n_terms_defined=len(ts),
                         terms_lost='; '.join(t for t in ts if t not in live), n_total=len(mg[m]),
                         FDR=min(fdr.get(t, 1) for t in ts), FDR_same=min(fsame.get(t, 1) for t in ts) if fsame else np.nan,
                         FDR_allev=min(fall.get(t, 1) for t in ts)))
    assign, allm = {}, {}
    for g in sorted(set().union(*mg.values())):
        c = [(int(K[t]), fdr[t], m) for m, ts in MODULES.items() for t in ts
             if fdr.get(t, 1) < 0.05 and g in genes[t].split(';')]
        assign[g] = min(c)[2]; allm[g] = ';'.join(sorted({x[2] for x in c}))
    for r in rows:
        m = r['module']
        r['n_assigned'] = sum(1 for v in assign.values() if v == m)
        r['assigned_genes'] = ', '.join(sorted(g for g, v in assign.items() if v == m))
        excl = sorted(g for g in mg[m] if ';' not in allm[g])
        r['n_exclusive'] = len(excl); r['exclusive_genes'] = ', '.join(excl)
        r['shared_genes'] = ', '.join(sorted(g for g in mg[m] if ';' in allm[g]))
        r['robust'] = bool(r['FDR'] < 0.05 and r['FDR_same'] < 0.05)
    M = pd.DataFrame(rows); M.to_csv(os.path.join(out, 'GO_module_summary.csv'), index=False)
    pd.DataFrame([(g, assign[g], allm[g]) for g in assign], columns=['gene', 'assigned_module', 'all_modules']) \
        .to_csv(os.path.join(out, 'GO_gene_module_assignment.csv'), index=False)

    H = []
    for m, t in HEAT:
        row = dict(module=m, term=t)
        for nm, df in res.items():
            h = df[df.Description == t]
            row[nm] = float(h.p_adjust_BH.iloc[0]) if len(h) else 1.0
            row[nm + '_n'] = int(h.Count.iloc[0]) if len(h) else 0
        H.append(row)
    pd.DataFrame(H).to_csv(os.path.join(out, 'GO_heatmap_data.csv'), index=False)

    # gene status: in a module / broad terms only / no enriched term / not tested
    sym_in_sig = set(';'.join(sig.genes).split(';'))
    bg = set(pd.read_csv('data/nonIEA_HPO_background.csv', dtype=str).ncbi_gene_id)
    untested = sorted(fw[~fw.ncbi_gene_id.isin(bg)].gene)
    tested = set(fw.gene) - set(untested)
    status = dict(module=sorted(assign), broad=sorted(sym_in_sig - set(assign)),
                  none=sorted(tested - sym_in_sig), untested=untested)
    json.dump(status, open(os.path.join(out, 'gene_status.json'), 'w'), indent=1)
    log['n_sig_terms'] = len(sig); log['n_clusters'] = int(sig.cluster.nunique())
    log.update({k: len(v) for k, v in status.items()})
    json.dump(log, open(os.path.join(out, 'go_summary.json'), 'w'), indent=1)
    print(json.dumps(log, indent=1))
    print(M[['module', 'n_total', 'n_assigned', 'FDR', 'FDR_same', 'robust', 'terms_lost']].to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])

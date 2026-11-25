"""Write the gene-set inputs for go_modules.py from a build_sets.py output.

Usage: python3 prep_inputs.py <S1c_curated*_genes.csv> <original|corrected> <outdir>
"""
import sys, os
import pandas as pd

flags, mode, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
df = pd.read_csv(flags, dtype={'ncbi_gene_id': str})
lit = pd.read_csv('literature_genes.csv', dtype=str)
lit = lit[lit.in_original_27 == 'yes'] if mode == 'original' else lit[lit.include_corrected == 'yes']

sel = df[df.selected].copy()

def group(r):
    if r.n_domains == 3: return 'All three'
    if r.optic and r.midline: return 'ONH + midline'
    if r.pituitary and r.midline: return 'Pituitary + midline'
    if r.optic and r.pituitary: return 'ONH + pituitary'
    return 'SOD term only'

sel['group'] = sel.apply(group, axis=1)
sel[['ncbi_gene_id', 'gene_symbol', 'group']].rename(columns={'gene_symbol': 'gene'}).to_csv(f'{out}/groups.csv', index=False)

fw = pd.concat([sel[['gene_symbol', 'ncbi_gene_id']].rename(columns={'gene_symbol': 'gene'}),
                lit[['gene', 'ncbi_gene_id']]]).drop_duplicates('gene')
fw.to_csv(f'{out}/framework_genes.csv', index=False)
lit[['gene', 'ncbi_gene_id']].to_csv(f'{out}/lit.csv', index=False)

# sensitivity set (as in the original analysis): HPO genes only, no literature genes.
# original: >=2 features in one disease, or the SOD term (125 genes, as published);
# corrected: >=2 features under one disease ID; the SOD term no longer supplies features.
same = df[df.selected_same_disease] if mode == 'original' else df[df.selected & (df.max_domains_same_disease >= 2)]
s = same[['gene_symbol', 'ncbi_gene_id']].rename(columns={'gene_symbol': 'gene'})
s.to_csv(f'{out}/same_disease_genes.csv', index=False)
print(f"{mode}: HPO={len(sel)} lit={len(lit)} framework={len(fw)} same-disease set={len(s)} | groups",
      sel.group.value_counts().to_dict())

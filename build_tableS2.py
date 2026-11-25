"""Assemble Table S2 (corrected) from the pipeline outputs."""
import json, re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
from go_modules import MODULES, HEAT

OUT = 'outputs/TableS2_SOD_gene_framework.xlsx'
flags = pd.read_csv('outputs/'+'S1c_curated_corr_genes.csv', dtype={'ncbi_gene_id': str}).set_index('gene_symbol')
flags_o = pd.read_csv('outputs/'+'S1c_curated_orig_genes.csv', dtype={'ncbi_gene_id': str}).set_index('gene_symbol')
fw = pd.read_csv('outputs/in_corr/framework_genes.csv', dtype=str)
grp = pd.read_csv('outputs/in_corr/groups.csv', dtype=str).set_index('gene').group
lit = pd.read_csv('literature_genes.csv', dtype=str).fillna('')
lit = lit[lit.include_corrected == 'yes'].set_index('gene')
st = json.load(open('outputs/out_corr/gene_status.json'))
asg = pd.read_csv('outputs/out_corr/GO_gene_module_assignment.csv').set_index('gene')
M = pd.read_csv('outputs/out_corr/GO_module_summary.csv')
Mo = pd.read_csv('outputs/out_orig/GO_module_summary.csv')
sig = pd.read_csv('outputs/out_corr/GO_significant_terms_with_clusters.csv').fillna({'module': ''})
sim = pd.read_csv('outputs/out_corr/GO_term_similarity_jaccard.csv', index_col=0)
F = pd.read_csv('outputs/out_corr/GO_Full_framework.csv').fillna({'genes': ''})
S = pd.read_csv('outputs/out_corr/GO_Full_framework_same_disease.csv')
A = pd.read_csv('outputs/out_corr/GO_Full_framework_all_evidence.csv')
fd = dict(zip(F.Description, F.p_adjust_BH)); fs = dict(zip(S.Description, S.p_adjust_BH)); fa = dict(zip(A.Description, A.p_adjust_BH))
robust = dict(zip(M.module, M.robust))
yn = lambda b: 'Yes' if bool(b) else 'No'

# ---------------- gene sheet
rows = []
for g in sorted(fw.gene, key=str.upper):
    in_hpo_sel = g in grp.index
    in_lit = g in lit.index
    src = 'HPO + published SOD report' if (in_hpo_sel and in_lit) else ('HPO' if in_hpo_sel else 'Published SOD report only')
    r = {'Gene': g, 'NCBI Gene ID': int(fw.set_index('gene').ncbi_gene_id[g]), 'Source': src}
    if g in flags.index:
        x = flags.loc[g]
        r.update({'Midline defect (HPO)': yn(x.midline), 'ONH or optic nerve aplasia (HPO)': yn(x.optic),
                  'Pituitary involvement (HPO)': yn(x.pituitary), 'No. of SOD features (HPO)': int(x.n_domains),
                  'Annotated to Septo-optic dysplasia (HP:0100842)': yn(x.SOD_term),
                  '≥2 features under the same disease ID': yn(x.max_domains_same_disease >= 2)})
    else:
        r.update({k: 'NA' for k in ['Midline defect (HPO)', 'ONH or optic nerve aplasia (HPO)', 'Pituitary involvement (HPO)',
                                    'No. of SOD features (HPO)', 'Annotated to Septo-optic dysplasia (HP:0100842)',
                                    '≥2 features under the same disease ID']})
    r['SOD feature group'] = grp.get(g, 'Not grouped (published report only)').replace('SOD term only', 'Not grouped (SOD term only)')
    r['Published case features'] = lit.published_SOD_features.get(g, '') if in_lit else ''
    r['Primary citation'] = lit.primary_citation.get(g, '') if in_lit else ''
    r['Evidence note'] = lit.evidence_note.get(g, '') if in_lit else ''
    status = 'Not tested' if g in st['untested'] else 'Module' if g in st['module'] else 'Broad terms only' if g in st['broad'] else 'No enriched term'
    r['Tested in GO analysis'] = 'No' if status == 'Not tested' else 'Yes'
    r['GO result'] = status
    r['Primary module'] = asg.assigned_module.get(g, '')
    r['All modules'] = asg.all_modules.get(g, '').replace(';', '; ')
    r['Robust module'] = yn(robust[r['Primary module']]) if r['Primary module'] else ''
    if status == 'Broad terms only':
        r['Significant broad GO terms (top 3)'] = '; '.join(sig[sig.genes.str.split(';').map(lambda s: g in s)].Description.head(3))
    else:
        r['Significant broad GO terms (top 3)'] = ''
    if g in flags.index:
        x = flags.loc[g]
        r['Midline terms matched'] = x.midline_terms if isinstance(x.midline_terms, str) else ''
        r['Optic terms matched'] = x.optic_terms if isinstance(x.optic_terms, str) else ''
        r['Pituitary terms matched'] = x.pituitary_terms if isinstance(x.pituitary_terms, str) else ''
        r['Disease IDs (midline)'] = x.midline_diseases if isinstance(x.midline_diseases, str) else ''
        r['Disease IDs (optic)'] = x.optic_diseases if isinstance(x.optic_diseases, str) else ''
        r['Disease IDs (pituitary)'] = x.pituitary_diseases if isinstance(x.pituitary_diseases, str) else ''
        r['Disease IDs (HP:0100842)'] = x.SOD_term_diseases if isinstance(x.SOD_term_diseases, str) else ''
        r['Disease IDs with ≥2 features'] = x.same_disease_ids if isinstance(x.same_disease_ids, str) else ''
    ch = []
    if g in flags_o.index and g in flags.index:
        a, b = flags_o.loc[g], flags.loc[g]
        for f, lab in [('midline', 'midline'), ('optic', 'optic'), ('pituitary', 'pituitary')]:
            if a[f] and not b[f]: ch.append(f'{lab} flag removed')
        if (a.max_domains_same_disease >= 2) and not (b.max_domains_same_disease >= 2): ch.append('no longer ≥2 features in one disease')
        go = {3: 'All three'}.get(int(a.n_domains), None)
        if int(a.n_domains) != int(b.n_domains): ch.append(f'features {int(a.n_domains)} → {int(b.n_domains)}')
    if g in ('VAX1', 'COL4A1'): ch.append('added from published reports')
    r['Change from previous version'] = '; '.join(ch)
    rows.append(r)
T = pd.DataFrame(rows)

# ---------------- other sheets
terms = pd.read_csv('outputs/curated_terms_used.csv')
figchg = pd.read_csv('outputs/figure2_changes.csv')
removed = pd.DataFrame([
    ('KDM6A', 'Removed', 'Midline flag rested on a corpus callosum term recorded in 0 of 6 individuals (OMIM:147920); one feature left.'),
    ('KMT2D', 'Removed', 'Midline flag rested on the same zero-frequency entry in OMIM:147920; one feature left.'),
    ('VAX1', 'Added', 'Published individual with small optic nerves and corpus callosum agenesis (Slavotinek et al., 2012); described in Pasca et al. 2025 but missing from the original 27.'),
    ('COL4A1', 'Added', 'Published individual with ONH and almost complete corpus callosum agenesis (Grego et al., 2019).'),
], columns=['Gene', 'Change', 'Reason'])
summary_changes = [
    ('Rule', 'Annotations recorded in 0 of n individuals are not counted as positive matches (3,083 such rows in the file).'),
    ('Rule', 'Septo-optic dysplasia (HP:0100842) is its own inclusion flag and no longer supplies the optic and midline features.'),
    ('Rule', 'The co-occurrence column is now ≥2 features annotated under the same disease ID, without the SOD-term pass. Disease IDs are listed for every match.'),
    ('Rule', 'Published-report genes require a published individual with ≥2 documented SOD features; VAX1 and COL4A1 added.'),
    ('Count', 'HPO-derived genes 144 → 142; all three features 32 → 25; ≥2 features in one disease ID 125 → 119; published-report genes 27 → 29; framework 155 → 155.'),
    ('Count', 'Groups: all three 25, ONH + midline 66, pituitary + midline 45, ONH + pituitary 1, SOD term only 5 (GPR161, RSPO2, TBX4, WDR11, WNT3; not grouped).'),
    ('Count', 'GO: 149 genes tested, 62 significant terms (was 59), 9 clusters, same seven modules with identical gene membership, same five robust modules.'),
    ('Count', 'Gene status: 64 in a module, 36 broad terms only, 49 no enriched term, 6 not tested.'),
    ('Note', 'Module FDR is the lowest FDR among the module\'s terms, as stated in Methods. The previous Table S2 reported dystroglycan at 1.7e-3, which is its parent term; the lowest is 3.4e-4.'),
]

MOD = []
for r in M.itertuples():
    o = Mo.set_index('module').loc[r.module]
    MOD.append({'Module': r.module, 'Defining GO terms': '; '.join(MODULES[r.module]), 'Lowest-FDR term': r.label_term,
                'Genes in module': r.n_total, 'Genes assigned': r.n_assigned, 'Only in this module': r.n_exclusive,
                'Shared with another module': r.n_total - r.n_exclusive, 'Lowest FDR': r.FDR, 'Previous lowest FDR': o.FDR,
                'Lowest FDR, ≥2 features in one disease ID': r.FDR_same, 'Lowest FDR, all evidence codes': r.FDR_allev,
                'Robust': yn(r.robust), 'Genes (assigned)': r.assigned_genes, 'Genes shared with other modules': r.shared_genes})
MOD = pd.DataFrame(MOD)
GO = sig[['GO_ID', 'Description', 'Count', 'BgTermSize', 'FoldEnrichment', 'pvalue', 'p_adjust_BH', 'cluster', 'module', 'genes']].copy()
GO.columns = ['GO ID', 'GO term', 'Framework genes in term', 'Background genes in term', 'Fold enrichment', 'P value', 'FDR',
              'Cluster', 'Module (defining term)', 'Genes']
GO['Genes'] = GO['Genes'].str.replace(';', ', ')
GO['In previous analysis'] = GO['GO term'].map(lambda t: 'Yes' if t in set(pd.read_csv('outputs/out_orig/GO_significant_terms_with_clusters.csv').Description) else 'New')
LIT = lit.reset_index()[['gene', 'ncbi_gene_id', 'report_group', 'primary_citation', 'citation_verified', 'published_SOD_features', 'evidence_note', 'in_original_27']]
LIT.columns = ['Gene', 'NCBI Gene ID', 'Report group (Pasca et al. 2025)', 'Primary citation', 'Citation checked', 'Published SOD features', 'Evidence note', 'In original 27']
LIT['NCBI Gene ID'] = LIT['NCBI Gene ID'].astype(int)

guide = [
    ('Gene', 'HGNC symbol.'), ('NCBI Gene ID', 'Entrez Gene identifier.'),
    ('Source', 'HPO: annotated to terms from ≥2 of the three SOD features, or to Septo-optic dysplasia (HP:0100842). Published SOD report: a published individual with ≥2 documented SOD features (sheet Literature genes).'),
    ('Midline / ONH or optic nerve aplasia / Pituitary (HPO)', 'Annotated in genes_to_phenotype.txt (downloaded November 2, 2026; MD5 6dac74395b2a9192dc0e71594093594f; HPO release 2026-09-01) to a curated term for that feature (sheet HPO terms used). Annotations recorded in 0 of n individuals are excluded. NA: gene absent from the HPO file.'),
    ('No. of SOD features (HPO)', 'Number of the three features annotated (0-3). HP:0100842 does not count toward any feature.'),
    ('Annotated to Septo-optic dysplasia (HP:0100842)', 'Gene annotated directly to the composite SOD term; an inclusion flag on its own.'),
    ('≥2 features under the same disease ID', 'At least two features annotated under one OMIM or Orphanet disease ID (sensitivity definition). A shared disease ID can combine findings from different patients, so this is not patient-level co-occurrence.'),
    ('SOD feature group', 'Feature combination used for the group analyses. Genes included only through HP:0100842, and published-report genes, are not grouped.'),
    ('Published case features / Primary citation / Evidence note', 'For published-report genes: the SOD features in the cited report and qualifiers on the evidence. Citations marked "to confirm" or "to add" still need the primary source.'),
    ('Tested in GO analysis', 'No if the gene is absent from the HPO background or lacks a qualifying non-IEA GO Biological Process annotation.'),
    ('GO result', 'Module: in a significant term that defines one of the seven modules. Broad terms only: only in significant broad terms. No enriched term: tested, in no significant term.'),
    ('Primary module / All modules', 'Primary module: the module whose most specific significant term (fewest background genes) contains the gene. All modules: every module whose significant terms contain it.'),
    ('Robust module', 'Primary module stays significant when the analysis is restricted to genes with ≥2 features under one disease ID.'),
    ('Terms matched / Disease IDs', 'The HPO terms and supporting disease IDs behind each feature flag, so every match can be traced.'),
    ('Change from previous version', 'What changed for this gene relative to the November 9, 2026 version.'),
]

# ---------------- write
wb = Workbook()
ARIAL = 'Arial'
hdr_font = Font(name=ARIAL, bold=True, size=10, color='FFFFFF'); hdr_fill = PatternFill('solid', fgColor='1F3A5F')
body = Font(name=ARIAL, size=10); title_font = Font(name=ARIAL, bold=True, size=12)
thin = Side(style='thin', color='D9D9D9')

def write(ws, df, title, start=3, widths=None, wrap_cols=()):
    ws.cell(1, 1, title).font = title_font
    for j, c in enumerate(df.columns, 1):
        cell = ws.cell(start, j, c); cell.font = hdr_font; cell.fill = hdr_fill
        cell.alignment = Alignment(wrap_text=True, vertical='center')
    for i, row in enumerate(df.itertuples(index=False), start + 1):
        for j, v in enumerate(row, 1):
            if isinstance(v, float) and pd.isna(v): v = None
            cell = ws.cell(i, j, v); cell.font = body
            cell.border = Border(bottom=thin)
            if df.columns[j - 1] in wrap_cols: cell.alignment = Alignment(wrap_text=True, vertical='top')
            else: cell.alignment = Alignment(vertical='top')
            if isinstance(v, float) and v < 0.01: cell.number_format = '0.00E+00'
            elif isinstance(v, float): cell.number_format = '0.000'
    ws.freeze_panes = ws.cell(start + 1, 2)
    for j, c in enumerate(df.columns, 1):
        w = (widths or {}).get(c, min(max(len(str(c)), *(len(str(x)) for x in df[c].head(200))) + 2, 45))
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[start].height = 42

ws = wb.active; ws.title = 'Table S2'
write(ws, T, 'Table S2. SOD gene framework: 155 genes linked to SOD features, with HPO annotations, supporting disease IDs, and GO module assignments (corrected, November 12, 2026)',
      wrap_cols=('Published case features', 'Evidence note', 'Midline terms matched', 'Optic terms matched', 'Pituitary terms matched', 'Change from previous version'))

ws = wb.create_sheet('Modules')
write(ws, MOD, 'GO modules in the full framework (non-IEA GO BP, HPO background; FDR < 0.05)', wrap_cols=('Defining GO terms', 'Genes (assigned)', 'Genes shared with other modules'),
      widths={'Defining GO terms': 50, 'Genes (assigned)': 45, 'Genes shared with other modules': 40})

ws = wb.create_sheet('GO terms')
write(ws, GO, 'Significant GO Biological Process terms in the full framework (62 terms; clustered by Jaccard distance, average linkage, cut height 0.8)',
      wrap_cols=('Genes',), widths={'Genes': 60, 'GO term': 42})

ws = wb.create_sheet('GO term similarity')
ws.cell(1, 1, 'Jaccard similarity between the framework genes in each pair of significant GO terms (1 − distance used for clustering)').font = title_font
ws.cell(3, 1, 'GO term').font = hdr_font; ws.cell(3, 1).fill = hdr_fill
for j, t in enumerate(sim.columns, 2):
    c = ws.cell(3, j, t); c.font = Font(name=ARIAL, bold=True, size=8, color='FFFFFF'); c.fill = hdr_fill
    c.alignment = Alignment(text_rotation=90, vertical='bottom'); ws.column_dimensions[get_column_letter(j)].width = 5
ws.row_dimensions[3].height = 230
for i, t in enumerate(sim.index, 4):
    ws.cell(i, 1, t).font = Font(name=ARIAL, size=9)
    for j, v in enumerate(sim.loc[t], 2):
        c = ws.cell(i, j, round(float(v), 3)); c.font = Font(name=ARIAL, size=8); c.number_format = '0.00'
        if v >= 0.2:
            shade = int(255 - min(v, 1) * 150); c.fill = PatternFill('solid', fgColor=f'{shade:02X}{shade + 20 if shade < 235 else 255:02X}FF')
ws.column_dimensions['A'].width = 48; ws.freeze_panes = 'B4'

ws = wb.create_sheet('HPO terms used')
write(ws, terms, 'Curated HPO terms defining each SOD feature (exact matches, no child terms). HP:0100842 is a separate inclusion flag.',
      widths={'HPO term': 44})

ws = wb.create_sheet('Literature genes')
write(ws, LIT, 'Published-report genes (29): primary citations and evidence', wrap_cols=('Primary citation', 'Published SOD features', 'Evidence note'),
      widths={'Primary citation': 42, 'Published SOD features': 45, 'Evidence note': 50})

ws = wb.create_sheet('Figure 2 changes')
write(ws, figchg, 'Figure 2: every cell that changes with the corrections (FDR < 0.05 shown in colour in the figure)', widths={'GO term': 40})

ws = wb.create_sheet('Changes')
ws.cell(1, 1, 'Changes from the November 9, 2026 version').font = title_font
r0 = 3
for k, (kind, txt) in enumerate(summary_changes):
    ws.cell(r0 + k, 1, kind).font = Font(name=ARIAL, bold=True, size=10)
    c = ws.cell(r0 + k, 2, txt); c.font = body; c.alignment = Alignment(wrap_text=True, vertical='top')
r1 = r0 + len(summary_changes) + 2
for j, c in enumerate(removed.columns, 1):
    x = ws.cell(r1, j, c); x.font = hdr_font; x.fill = hdr_fill
for i, row in enumerate(removed.itertuples(index=False), r1 + 1):
    for j, v in enumerate(row, 1):
        c = ws.cell(i, j, v); c.font = body; c.alignment = Alignment(wrap_text=True, vertical='top')
ws.column_dimensions['A'].width = 10; ws.column_dimensions['B'].width = 110; ws.column_dimensions['C'].width = 90

ws = wb.create_sheet('Column guide')
ws.cell(1, 1, 'Column guide').font = title_font
for i, (a, b) in enumerate(guide, 3):
    ws.cell(i, 1, a).font = Font(name=ARIAL, bold=True, size=10)
    c = ws.cell(i, 2, b); c.font = body; c.alignment = Alignment(wrap_text=True, vertical='top')
n = len(T) + 3
s0 = len(guide) + 5
ws.cell(s0, 1, 'Summary (formulas on sheet Table S2)').font = Font(name=ARIAL, bold=True, size=11)
summ = [('Total genes', f"=COUNTA('Table S2'!A4:A{n})"),
        ('Source: HPO only', f"=COUNTIF('Table S2'!C4:C{n},\"HPO\")"),
        ('Source: HPO + published SOD report', f"=COUNTIF('Table S2'!C4:C{n},\"HPO + published*\")"),
        ('Source: published SOD report only', f"=COUNTIF('Table S2'!C4:C{n},\"Published*\")"),
        ('All three features', f"=COUNTIF('Table S2'!J4:J{n},\"All three\")"),
        ('ONH + midline', f"=COUNTIF('Table S2'!J4:J{n},\"ONH + midline\")"),
        ('Pituitary + midline', f"=COUNTIF('Table S2'!J4:J{n},\"Pituitary + midline\")"),
        ('In a module', f"=COUNTIF('Table S2'!O4:O{n},\"Module\")"),
        ('Broad terms only', f"=COUNTIF('Table S2'!O4:O{n},\"Broad terms only\")"),
        ('No enriched term', f"=COUNTIF('Table S2'!O4:O{n},\"No enriched term\")"),
        ('Not tested', f"=COUNTIF('Table S2'!O4:O{n},\"Not tested\")")]
for i, (a, f) in enumerate(summ, s0 + 1):
    ws.cell(i, 1, a).font = body; ws.cell(i, 2, f).font = body
ws.column_dimensions['A'].width = 48; ws.column_dimensions['B'].width = 110
wb.save(OUT)
print('wrote', OUT, T.shape)
print(T.columns.tolist())

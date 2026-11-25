# SOD gene framework

Code for the septo-optic dysplasia (SOD) gene framework and Gene Ontology (GO) enrichment
analysis.

SOD is diagnosed when at least two of three features are present: optic nerve hypoplasia
(ONH), pituitary involvement, and midline brain defects. This repository builds a framework
of genes linked to these features from the Human Phenotype Ontology (HPO) and published
reports, tests which biological processes they share, and groups enriched GO terms into
pathway modules. It reproduces Table S2 and the data behind Figure 2.

## Quick start

```bash
git clone https://github.com/shaybehesht/sod-gene-framework.git
cd sod-gene-framework
pip install -r requirements.txt

# Add the HPO annotation file used in the paper (downloaded November 2, 2026).
# make_g2p.py checks its MD5: 6dac74395b2a9192dc0e71594093594f
cp /path/to/genes_to_phenotype.txt data/

bash run_all.sh
```

Results are written to `outputs/`. Reference outputs from the published run are in
`outputs_reference/`; a fresh run reproduces `outputs_reference/out_corr/` exactly.

The GO gene sets are included, so R is not needed to run the pipeline. To regenerate them
(R with Bioconductor org.Hs.eg.db and GO.db 3.23.1):

```bash
Rscript R/export_go_annotations.R data/genes_to_phenotype.txt
```

## Methods

### SOD gene framework

Gene-phenotype relationships were taken from the HPO annotation file `genes_to_phenotype.txt`
(https://purl.obolibrary.org/obo/hp/genes_to_phenotype.txt; downloaded November 2, 2026; MD5
`6dac74395b2a9192dc0e71594093594f`), which links genes to HPO terms through their OMIM and
Orphanet disease annotations. Term names and hierarchy come from HPO release 2026-09-01
(`data/hp.obo`), and obsolete term IDs were mapped to current IDs.

Each SOD feature was defined by a curated set of developmental HPO terms
(`curated.py`, listed in `TableS3_curated_HPO_terms.csv`):

| Feature | Terms | Examples |
| --- | --- | --- |
| Midline defects | 12 | corpus callosum and septum pellucidum abnormalities |
| ONH | 4 | optic nerve hypoplasia, optic disc hypoplasia, optic nerve aplasia |
| Pituitary involvement | 19 | structural pituitary and stalk abnormalities, hypopituitarism, panhypopituitarism |

Rules:

* Terms are matched exactly, without child terms.
* Annotations recorded in zero individuals (`0/n` or `0%`) are not counted.
* Broad terms that do not specify a developmental pituitary abnormality, such as
  Abnormality of the hypothalamus-pituitary axis (HP:0000864), are excluded.
* A gene is included if it is annotated to terms from at least two of the three features.
* Septo-optic dysplasia (HP:0100842) is used only as an inclusion flag. A gene annotated
  directly to it is included, but the term does not count toward any single feature.

This gives 142 HPO-derived genes. We also compiled 29 genes from published SOD reports
(`literature_genes.csv`: 27 from Pasca et al., 2025, plus VAX1 and COL4A1 from primary
reports). Sixteen were already HPO-derived, and the remaining 13 were added, for a
155-gene framework.

HPO-derived genes are grouped by feature combination: all three features (25 genes), ONH
with midline defects (66), pituitary involvement with midline defects (45), ONH with
pituitary involvement (1), and SOD term only (5). Literature-only genes are not grouped.

### GO enrichment

* GO Biological Process annotations from Bioconductor org.Hs.eg.db 3.23.1 and GO.db 3.23.1
  (GO source dated January 23, 2026), including ancestor terms and excluding
  inferred-from-electronic-annotation (IEA) and no-data (ND) evidence. The gene sets are
  exported by `R/export_go_annotations.R` and stored in `data/go_sets_*.csv.gz`.
* Background: all genes in the HPO file with a qualifying GO annotation (5,083 genes).
* Terms with 10 to 500 background genes were tested (3,717 terms).
* One-sided hypergeometric test, Benjamini-Hochberg FDR across all tested terms;
  FDR < 0.05 is significant (`go_ora.py`).
* Analyses were run on the 29 reported genes, the full framework, and each feature group.

### Modules

Significant terms from the full framework were clustered by Jaccard distance between the
framework genes they contain (average linkage, cut height 0.8), which gives nine clusters.
Seven modules were defined from the specific GO terms within these clusters (`MODULES` in
`go_modules.py`): two clusters form the axis and left-right patterning module, and one
cluster contains only broad developmental terms (limb and appendage development, cell fate
commitment) and is not a module. The seven modules are: Hedgehog signaling; forebrain and
pituitary development; axis and left-right patterning; dystroglycan O-mannosylation; axon
guidance; microtubule organization and sister chromatid cohesion; and mRNA splicing
regulation. Module FDR is the lowest FDR among its terms.

A gene can belong to several modules. For Table S2, each gene is given one primary module:
the module of the most specific significant term (fewest background genes) that contains it.

### Sensitivity analyses

* All GO evidence codes, including IEA.
* A stricter framework keeping only HPO genes with at least two SOD features annotated under
  the same disease ID (119 genes). Modules that stay significant here are called robust.

### Expected results

| | Value |
| --- | --- |
| HPO-derived genes | 142 |
| Framework genes | 155 (149 testable) |
| Significant GO terms (full framework) | 62 |
| Genes in a module | 64 |
| Genes with broad terms only | 36 |
| Genes with no enriched term | 49 |
| Genes not testable | 6 |

## Repository layout

| File | Purpose |
| --- | --- |
| `run_all.sh` | Runs the full pipeline |
| `make_g2p.py` | Checks the HPO file MD5 and caches it as `data/g2p.pkl` |
| `hpo.py` | Reads `hp.obo` (hierarchy, obsolete IDs) |
| `curated.py` | Curated HPO term lists for the three features |
| `build_sets.py` | Matches genes to features; per-gene flags, matched terms, disease IDs |
| `prep_inputs.py` | Framework, feature-group, literature, and sensitivity gene sets |
| `go_ora.py` | Hypergeometric test and BH correction |
| `go_modules.py` | Enrichment runs, term clustering, modules, primary-module assignment, Figure 2 data |
| `compare_runs.py` | HPO term table and comparison with the earlier version |
| `build_tableS2.py` | Writes `outputs/TableS2_SOD_gene_framework.xlsx` |
| `R/export_go_annotations.R` | Exports GO gene sets, backgrounds, and term names from org.Hs.eg.db / GO.db |
| `literature_genes.csv` | Published SOD genes with citations |
| `TableS3_curated_HPO_terms.csv` | Curated HPO terms (Table S3) |
| `data/` | `hp.obo` (HPO 2026-09-01), GO gene sets (`go_sets_*.csv.gz`), GO backgrounds, GO term names |
| `outputs_reference/` | Outputs of the published run |

## Earlier version

`build_sets.py` has two switches that reproduce an earlier version of the framework when
set to 0: `DROP_ZERO` (count 0/n annotations) and `SOD_SEP` (let HP:0100842 supply the optic
and midline features). The published analysis uses `DROP_ZERO=1 SOD_SEP=1`. `run_all.sh`
runs both and `compare_runs.py` lists the differences.

## Notes

* `genes_to_phenotype.txt` is not included. Download it from HPO; the MD5 check confirms
  the version.
* `data/go_sets_*.csv.gz` hold every tested GO term with all of its background genes, so
  the enrichment can be checked without R. Newer Bioconductor releases will give slightly
  different GO sets; use version 3.23.1 to reproduce the published numbers exactly.
* Tested with Python 3.13, pandas 2, scipy 1.18.

## Acknowledgments

Analysis code was written and edited with the help of Claude (Anthropic) and reviewed by the
authors.

## Contact

Shaghayegh (Shay) Beheshti · GitHub [@shaybehesht](https://github.com/shaybehesht)

#!/usr/bin/env bash
# Rebuild the SOD gene framework, the GO analysis, and Table S2, original and corrected.
# Run from this folder. Needs Python 3 with pandas, numpy, scipy, openpyxl.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p outputs

# GO gene sets in data/ were exported with: Rscript R/export_go_annotations.R data/genes_to_phenotype.txt
# (needs R with org.Hs.eg.db and GO.db 3.23.1; rerun only to regenerate them)

# 0. cache the HPO file (checks the MD5 of the November 2, 2026 download)
[ -f data/g2p.pkl ] || python3 make_g2p.py

# 1. gene framework: original rules, then both corrections
DROP_ZERO=0 SOD_SEP=0 ONLY=S1c_curated TAG=_orig python3 build_sets.py
DROP_ZERO=1 SOD_SEP=1 ONLY=S1c_curated TAG=_corr python3 build_sets.py

# 2. gene sets for the GO step
python3 prep_inputs.py outputs/S1c_curated_orig_genes.csv original outputs/in_orig
python3 prep_inputs.py outputs/S1c_curated_corr_genes.csv corrected outputs/in_corr

# 3. GO enrichment, term clustering, modules
for m in orig corr; do
  python3 go_modules.py outputs/in_$m/framework_genes.csv outputs/in_$m/groups.csv \
    outputs/in_$m/lit.csv outputs/out_$m outputs/in_$m/same_disease_genes.csv > outputs/go_$m.log
done

# 4. comparison tables and Table S2
python3 compare_runs.py > outputs/figure2_changes.txt
python3 build_tableS2.py
echo "done: outputs/TableS2_SOD_gene_framework.xlsx"

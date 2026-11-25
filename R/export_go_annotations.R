# Export the GO Biological Process gene sets used by go_ora.py.
#
# Usage (from the repository root):
#   Rscript R/export_go_annotations.R data/genes_to_phenotype.txt
#
# Writes, for each evidence setting (nonIEA_HPO, allEvidence_HPO):
#   data/<tag>_background.csv    background genes (HPO genes with >=1 qualifying GO BP annotation)
#   data/go_sets_<tag>.csv.gz    one row per GO term and background gene (terms with 10-500 genes)
#   data/go_term_names.csv       GO ID and term name
#   data/R_sessionInfo.txt       package versions
#
# Published run: org.Hs.eg.db 3.23.1 and GO.db 3.23.1 (GO source dated January 23, 2026).
# Annotations include ancestor terms (GOALL). nonIEA excludes IEA and ND evidence;
# allEvidence keeps all evidence codes. The BP root term (GO:0008150) is excluded.
# This reproduces the gene sets of the original analysis script (saved there as
# <tag>_tested_GO_gene_sets.rds).

suppressPackageStartupMessages({
  library(AnnotationDbi); library(org.Hs.eg.db); library(GO.db)
})
args <- commandArgs(trailingOnly = TRUE)
hpo_file <- if (length(args)) args[1] else "data/genes_to_phenotype.txt"
MIN_SIZE <- 10L; MAX_SIZE <- 500L

hpo <- read.delim(hpo_file, sep = "\t", quote = "", comment.char = "", colClasses = "character")
names(hpo) <- sub("^X\\.", "", sub("^#", "", names(hpo)))
hpo_genes <- unique(hpo$ncbi_gene_id[nzchar(hpo$ncbi_gene_id)])

ann <- AnnotationDbi::select(org.Hs.eg.db, keys = hpo_genes, keytype = "ENTREZID",
                             columns = c("GOALL", "EVIDENCEALL", "ONTOLOGYALL"))
ann <- unique(ann[!is.na(ann$GOALL) & !is.na(ann$ONTOLOGYALL) &
                   ann$ONTOLOGYALL == "BP" & ann$GOALL != "GO:0008150", ])   # drop the BP root

export <- function(a, tag) {
  a <- unique(a[, c("GOALL", "ENTREZID")])
  bg <- sort(unique(a$ENTREZID))
  sizes <- table(a$GOALL)
  keep <- names(sizes)[sizes >= MIN_SIZE & sizes <= MAX_SIZE]
  a <- a[a$GOALL %in% keep, ]
  write.csv(data.frame(ncbi_gene_id = bg), sprintf("data/%s_background.csv", tag), row.names = FALSE)
  out <- gzfile(sprintf("data/go_sets_%s.csv.gz", tag), "w")
  write.csv(data.frame(GO_ID = a$GOALL, ncbi_gene_id = a$ENTREZID), out, row.names = FALSE)
  close(out)
  message(tag, ": ", length(bg), " background genes, ", length(keep), " terms")
  keep
}
t1 <- export(ann[!is.na(ann$EVIDENCEALL) & !ann$EVIDENCEALL %in% c("IEA", "ND"), ], "nonIEA_HPO")
t2 <- export(ann, "allEvidence_HPO")

ids <- sort(unique(c(t1, t2)))
write.csv(data.frame(GO_ID = ids, Description = Term(GOTERM[ids])), "data/go_term_names.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "data/R_sessionInfo.txt")

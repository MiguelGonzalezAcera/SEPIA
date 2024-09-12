# 05-2 - Gene Ontology analysis of common and unique genes
# This script runs a GO ORA of the obtained lists for the venn diagrams

# Required:
# 00 - DifferentialExpressionMouseModels
# 05-1 - VennDiagrams

suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(GOstats))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(enrichplot))


go_venn <- function(universe, genelist) {
    # Load the universe
    load(universe)

    # Select annotation
    database <- org.Mm.eg.db

    # Read the gene file
    genes <- readLines(genelist)

    # Transform the ensembl names into gene symbol. NOTE that the name of the variable must change.
    entrezgeneids <- as.character(mapIds(database, as.character(genes), "ENTREZID", "ENSEMBL"))

    # Obtain universe ids
    universeids <- unique(as.character(mapIds(database, as.character(rownames(counts(dss))), "ENTREZID", "ENSEMBL")))

    # Do the ORA
    hgCutoff <- 0.1
    x <- enrichGO(entrezgeneids, database, ont="BP", pvalueCutoff = hgCutoff, readable = TRUE,
                    pAdjustMethod = "BH", universe = universeids, minGSSize = 1, maxGSSize = 1000)

    # Transform the result into a data frame
    GOtable <- as.data.frame(x)

    outname <- gsub("VennGeneLists/", "", genelist, fixed = TRUE)
    outname <- gsub("_ensembl.txt", "", outname, fixed = TRUE)

    # Save the data frame
    write.table(GOtable, file = sprintf("VennGOresult/%s_BP.tsv", outname),
                sep = "\t", row.names = FALSE)

    # Save the objects
    save(x, file = sprintf("VennGOresult/%s_BP.Rda", outname))
}

# Get the universes of each tissue
univ_col <- "detables/AcDSS_universe.Rda"
univ_ile <- "detables/TNFdAREIle_universe.Rda"

# Get the lists of genes
genefiles_ile <- list.files("VennGeneLists", pattern = "ileum.*.txt", full.names = TRUE)
genefiles_col <- list.files("VennGeneLists", pattern = "colon.*.txt", full.names = TRUE)

# Run through the genelists and enrich each and every one of them
for (gene_ile in genefiles_ile) {
    go_venn(univ_ile, gene_ile)
}

for (gene_col in genefiles_col) {
    go_venn(univ_col, gene_col)
}
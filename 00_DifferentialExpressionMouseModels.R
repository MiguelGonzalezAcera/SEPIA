# 00 - Differential expression analyses of mouse models
# This script performs the differential expression for each mouse model separately.
# It must be run initially, as it generates the neccesary tables for posterior analyses

# Raw counts files are obtained by mapping the raw fastq files to the mouse genome and
# extracting the counts per gene with featureCounts

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(org.Mm.eg.db))
# suppressPackageStartupMessages(library(data.table))
# suppressPackageStartupMessages(library(gsubfn))

# Read the table with the location of all raw data
SEPIAtable <- read.table('/DATA/Thesis_proj/SEPIA_code_repo/SEPIA/raw_counts/table_reference.tsv', header=TRUE)

# Iter through the models for the differential expression analyses
models <- SEPIAtable$Model

for (model in models) {
  # Load the counts table
  Counts_tab <- read.table(SEPIAtable[SEPIAtable$Model == model,][['Raw_counts']], fileEncoding = "UTF8", header = TRUE)

  # Load the design table
  sampleTableSingle <- read.table(SEPIAtable[SEPIAtable$Model == model,][['Design']], fileEncoding = "UTF8")

  # Get the control and sample names
  control <- SEPIAtable[SEPIAtable$Model == model,][['Control']]
  sample <- SEPIAtable[SEPIAtable$Model == model,][['Sample']]

  # Move the gene IDs as row names
  row.names(Counts_tab) <- Counts_tab$Geneid
  Counts_tab$Geneid <- NULL

  # Select the columns specified in the provided design and resort the genes
  Counts_tab <- Counts_tab[, row.names(sampleTableSingle)]
  Counts_tab <- Counts_tab[order(row.names(Counts_tab)), ]

  # Add row names as column and subset control samples
  sampleTableSingle$rn <- row.names(sampleTableSingle)
  control_samples <- sampleTableSingle[sampleTableSingle$Tr1 == control,][['rn']]

  # Design model matrix, including batch effect correction
  Tr1 <- relevel(factor(sampleTableSingle$Tr1), control)
  design <- model.matrix(~ Tr1)

  # --------------------------------------------------------------

  # Create the experiment from a SummarizedExperiment object
  dss <- DESeqDataSetFromMatrix(countData = Counts_tab,
                                colData = sampleTableSingle,
                                design = design)

  # Save the universe (of genes). This is important for downstream analyses
  save(dss, file = sprintf("detables/%s_universe.Rda", model))

  # filter the counts by minimum
  keep <- rowSums(counts(dss)) >= 15
  dss <- dss[keep, ]

  # Run the analysis
  dds <- DESeq(dss, betaPrior = FALSE)

  # --------------------------------------------------------------

  # Save the normalized counts
  # Get the table from the result
  norm_counts <- counts(estimateSizeFactors(dds), normalized = TRUE)

  # Get the names of the columns
  norm_counts_colnames <- colnames(norm_counts)
  # Add the gene names as a new column
  norm_counts <- cbind(norm_counts, as.character(mapIds(org.Mm.eg.db, as.character(rownames(norm_counts)), 'SYMBOL', 'ENSEMBL')))
  # Rename the columns with the new name
  colnames(norm_counts) <- c(norm_counts_colnames, "Genename")

  # Save the object both as table and as R object
  write.table(norm_counts, file = sprintf("detables/%s_norm_counts.tsv", model), sep="\t")
  df_norm <- as.data.frame(norm_counts)
  save(df_norm, file = sprintf("detables/%s_norm_counts.Rda", model))

  # --------------------------------------------------------------

  # Get the EnsemblIDs from the counts table
  df_norm$EnsGenes <- rownames(df_norm)

  # Get the result of the differential expression analyses
  res <- results(dds, name = paste("Tr1", sample, sep = ""), cooksCutoff = FALSE)

  # Save the full result object
  # Contrast name will be replaced by the sample and controls
  res_name <- paste("detables/", sample, "_", control, ".Rda", sep = "")
  save(res, file = res_name)

  # A simple helper function that makes a so-called "MA-plot", i.e. a scatter plot of
  # log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis).
  MA_name <- paste("plots/", sample, "_", control, "_MA.png")
  png(file = MA_name)
  plotMA(res)
  dev.off()

  # Transform result into data frame
  resdf <- data.frame(res)[complete.cases(data.frame(res)), ]

  # Transform row ensembl IDs into column
  resdf$EnsGenes <- rownames(resdf)

  # Add also gene symbols
  resdf$Genes <- as.character(mapIds(org.Mm.eg.db, as.character(rownames(resdf)), "SYMBOL", "ENSEMBL"))

  # Save table with all the new names. Replace contrast
  res_tab_name <- paste("detables/", sample, "_", control, ".tsv", sep = "")
  write.table(resdf, file = res_tab_name, sep = "\t", row.names = FALSE)

  # Subset design table with sample to get vector
  sample_samples <- sampleTableSingle[sampleTableSingle$Tr1 == sample, ][["rn"]]

  # Merge res table with the counts of its samples and controls
  resdf_wcounts <- merge(resdf, df_norm[,c("EnsGenes", control_samples, sample_samples)], by="EnsGenes", all.x=TRUE)

  # filter by normalized counts in order to remove false positives
  # Oder of stuff: Select samples or control columns, transform to numeric with the function up,
  # transform to a data matrix, get the medians, Boolean on who's under 25, select rows
  resdf_wcounts$FLAG <- ifelse((rowMedians(data.matrix(sapply(resdf_wcounts[control_samples], as.numeric))) > 25) | (rowMedians(data.matrix(sapply(resdf_wcounts[sample_samples], as.numeric))) > 25), 'OK', 'WARN: Inconsinstent Counts')

  #Save new table
  res_exp_tab_name = paste("detables/", sample, "_", control, "_expanded.tsv", sep = "")
  write.table(resdf_wcounts, res_exp_tab_name, sep = "\t", row.names = FALSE)

  # --------------------------------------------------------------

  # Transform the counts using Variance Stabilizing Transformation
  vsd <- vst(dds)

  # Get the counts from the transformed object
  tr_counts <- assay(vsd)

  # Get the names of the columns
  tr_counts_colnames <- colnames(tr_counts)
  # Add the gene names as a new column
  tr_counts <- cbind(tr_counts, as.character(mapIds(org.Mm.eg.db, as.character(rownames(tr_counts)), 'SYMBOL', 'ENSEMBL')))
  # Rename the columns with the new name
  colnames(tr_counts) <- c(tr_counts_colnames, "Genename")

  # Save the object both as table and as R object
  write.table(tr_counts, file = sprintf("detables/%s_tr_counts.tsv", model), sep="\t")
  df_norm <- as.data.frame(tr_counts)
  save(df_norm, file = sprintf("detables/%s_tr_counts.Rda", model))

}

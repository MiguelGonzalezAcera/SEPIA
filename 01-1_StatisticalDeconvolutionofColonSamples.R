# 01-1 - Statistical deconvolution of colon samples - Part 1
# This script performs the statistical deconvolution of the bulk RNAseq tables.

# Signature matrix is obtained processing the single cell dataset and
# running the DWLS sing matrix generation functions on them

library(DWLS)

# Determine the files with the normalized counts
# list the file names for ileum and colon
files_col <-  c(
    "detables/AcDSS_norm_counts.Rda",
    "detables/cDSS_norm_counts.Rda",
    "detables/Casp8dIECCol_norm_counts.Rda",
    "detables/TC_norm_counts.Rda",
    "detables/OxC_norm_counts.Rda",
    "detables/AcTNBS_norm_counts.Rda",
    "detables/cTNBS_norm_counts.Rda",
    "detables/TNFdARECol_norm_counts.Rda",
    "detables/Hhepa_norm_counts.Rda",
    "detables/Crode_norm_counts.Rda"
)

# Load the signature matrix from GSE168033
load("deconvolution/signature.Rda")

# loop through each file and execute de deconv.
for (file in files_col) {
    # Load the normalized counts obtained in the differential expression assay
    load(file)

    # Remove NAs
    df_norm <- df_norm[complete.cases(df_norm), ]

    # Remove unnecessary columns
    df_norm$Genename <- NULL

    # Turn all number columns to numeric (security measure)
    rows_df <- rownames(df_norm)
    df_norm <- sapply(df_norm, function(x) as.numeric(as.character(x)))
    rownames(df_norm) <- rows_df

    # Create the deconvolution table object and names
    deconv <- NULL
    deconv_names <- c()

    # Iter through each sample to obtain the proportions of cell types
    for (i in seq_len(ncol(df_norm))) {
    
        # select the column in question
        df_norm_i <- df_norm[, i]
        
        # Trim the data (required for dwls)
        tr <- trimData(signature_matrix, df_norm_i)
        
        # Deconvolve
        deconv_i <- solveDampenedWLS(tr$sig, tr$bulk)
        
        # Add new column and names to table
        deconv <- cbind(deconv, deconv_i)
        deconv_names <- c(deconv_names, colnames(df_norm)[i])
    }

    # Name the thew table with the chosen names
    colnames(deconv) <- deconv_names

    # get the name of the experiment in question
    tname <- gsub("detables/", "", file)
    tname <- gsub("_norm_counts.Rda", "", tname)

    write.table(deconv, file = sprintf("deconvolution/%s_deconv.tsv", tname), sep = "\t")
}
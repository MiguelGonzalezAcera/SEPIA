# 03 - Ridge plots of cell markers
# This script generates ridge plots of cell markers for the models of the dataset.

# Required:
# 00 - DifferentialExpressionMouseModels

library(ggplot2)
library(ggridges)

files_col <- c(
    "detables/AcDSS_Healthy.Rda",
    "detables/cDSS_Healthy.Rda",
    "detables/Casp8dIECCol_Casp8floxCol.Rda",
    "detables/TC_RKO.Rda",
    "detables/OxC_Healthy.Rda",
    "detables/AcTNBS_Healthy.Rda",
    "detables/cTNBS_Healthy.Rda",
    "detables/TNFdARECol_WTCol.Rda",
    "detables/Hhepa_Healthy.Rda",
    "detables/Crode_Healthy.Rda"
)

files_ile <- c(
    "detables/Casp8dIECIle_Casp8floxIle.Rda",
    "detables/TNFdAREIle_WTIle.Rda",
    "detables/Everm_Healthy.Rda"
)

# Establish the names for ileum and colon
tabnames_col <- c(
    "AcDSS", "cDSS", "Casp8ΔIEC", "TC", "OxC",
    "AcTNBS", "cTNBS", "TNFΔARE", "Hhepa", "Crode"
)

tabnames_ile <- c("Casp8ΔIEC", "TNFΔARE", "Everm")

# Get the lists of genes
genefiles <- c(
    "cell_type_markers/InfMacrophages_ensembl.txt",
    "cell_type_markers/T_cells_ensembl.txt",
    "cell_type_markers/B_cells_ensembl.txt",
    "cell_type_markers/Goblet_ensembl.txt",
    "cell_type_markers/TAProg_ensembl.txt",
    "cell_type_markers/EntericNeuron_ensembl.txt"
)

# Loop through each of the lists producing the plots for colon and ileum
for (genelist in genefiles) {
    # Read the list of genes
    genes <- readLines(genelist)

    # Init the resulting object
    ridge_df <- NULL

    # loop through each model to get the relevant genes
    for (model in files_col) {
        # Ignore B cells in the case of TC
        if (genelist == "cell_type_markers/B_cells_ensembl.txt" && model == "detables/TC_RKO.Rda") next

        # Read the model in question and transform to dataframe
        load(model)
        full_df <- as.data.frame(res)

        # Get the fold change for the selected genes
        FC_df <- full_df[genes, "log2FoldChange", drop = FALSE]

        # Add a column with the name of the model in question
        FC_df$model <- tabnames_col[match(model, files_col)]

        # Add a column with the gene names
        FC_df$EnsGenes <- rownames(FC_df)

        # Merge with the resulting object
        if (is.null(ridge_df) == TRUE) {
            # If the final df is empty, fill it with one column
            ridge_df <- FC_df
        } else {
            # If not, add the column to the df
            ridge_df <- rbind(ridge_df, FC_df)
        }
    }
    
    # Relevel the model column
    ridge_df$model <- factor(ridge_df$model, levels = rev(unique(ridge_df$model)))

    # Plot the ridge plot
    ridgeplot <- ggplot(ridge_df, aes(x = log2FoldChange, y = model, fill = stat(x))) + 
        geom_density_ridges_gradient(scale = 2) +
        scale_fill_gradientn(colors = c("#7F00FF", "#7F00FF", "#7F00FF","#7F00FF", "#7F00FF","#7F00FF","#7F00FF","#7F00FF","#7F00FF","#7F00FF","#BD76EB", "#faebd7", "#FDB66C","#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000"), limits = c(-22, 22))

    # Save the plot
    tname <- gsub("cell_type_markers/", "", genelist)
    tname <- gsub("_ensembl.txt", "", tname)

    ggsave(sprintf("plots/%s_colon_ridge.svg", tname), plot = ridgeplot, width = 12, height = 12, units = 'cm', dpi = 600)

    # reinit for the ileal samples
    ridge_df <- NULL

    # loop through each model to get the relevant genes
    for (model in files_ile) {
        # Ignore B cells in the case of TC
        if (genelist == "cell_type_markers/B_cells_ensembl.txt" && model == "detables/TC_RKO.Rda") next

        # Read the model in question and transform to dataframe
        load(model)
        full_df <- as.data.frame(res)

        # Get the fold change for the selected genes
        FC_df <- full_df[genes, "log2FoldChange", drop = FALSE]

        # Add a column with the name of the model in question
        FC_df$model <- tabnames_ile[match(model, files_ile)]

        # Add a column with the gene names
        FC_df$EnsGenes <- rownames(FC_df)

        # Merge with the resulting object
        if (is.null(ridge_df) == TRUE) {
            # If the final df is empty, fill it with one column
            ridge_df <- FC_df
        } else {
            # If not, add the column to the df
            ridge_df <- rbind(ridge_df, FC_df)
        }
    }
    
    # Relevel the model column
    ridge_df$model <- factor(ridge_df$model, levels = rev(unique(ridge_df$model)))

    # Plot the ridge plot
    ridgeplot <- ggplot(ridge_df, aes(x = log2FoldChange, y = model, fill = stat(x))) + 
        geom_density_ridges_gradient(scale = 2) +
        scale_fill_gradientn(colors = c("#7F00FF", "#7F00FF", "#7F00FF","#7F00FF", "#7F00FF","#7F00FF","#7F00FF","#7F00FF","#7F00FF","#7F00FF","#BD76EB", "#faebd7", "#FDB66C","#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000", "#FF8000"), limits = c(-22, 22))

    # Save the plot
    tname <- gsub("cell_type_markers/", "", genelist)
    tname <- gsub("_ensembl.txt", "", tname)

    ggsave(sprintf("plots/%s_ileum_ridge.svg", tname), plot = ridgeplot, width = 12, height = 12, units = 'cm', dpi = 600)

}
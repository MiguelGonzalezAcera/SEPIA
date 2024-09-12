# 05-4 - Term Similarity Network of GO terms
# This script generates a table with a network of the term similarity of GO results
# Table then can be loaded in Cytoscape

# Required:
# 00 - DifferentialExpressionMouseModels
# 05-1 - VennDiagrams
# 05-2 - GOofVennLists

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(simplifyEnrichment))

# Get the GO results where we will perform the networks
filelist <- c(
    "VennGOresult/colon_common_genes_up_BP.Rda",
    "VennGOresult/ileum_common_genes_dw_BP.Rda",
    "VennGOresult/ileum_common_genes_up_BP.Rda"
)

# Iter through the tables for each experiment
for (res in filelist) {
    # Load the GO experiment and get it as table
    load(res)
    GOtable <- as.data.frame(x)

    # Run the term similarity for the matrix of comparisons
    mat <- term_similarity_from_enrichResult(x, method = "kappa")

    # Filter the matrix for only our elements
    mat_filt <- mat[GOtable$ID, GOtable$ID]

    # Cluster the matrix and save the plot in question
    df <- simplifyEnrichment(mat[GOtable$ID, GOtable$ID], plot = FALSE)

    # Get the 20 biggest clusters and obtain the terms belonging to them
    if (length(names(sort(table(df$cluster), decreasing = TRUE))) < 20) {
        maxlen <- length(names(sort(table(df$cluster), decreasing = TRUE)))
    } else {
        maxlen <- 20
    }

    sel_clusters <- names(sort(table(df$cluster), decreasing = TRUE))[1:maxlen]
    df_biggest <- subset(df, cluster %in% sel_clusters)

    # Filter that to get only the best 10 elements per cluster to get the nodes
    final_sel <- c()
    for (i in sel_clusters) {
        cl_list <- df_biggest[df_biggest$cluster == i,]$id
        cl_df <- subset(GOtable, ID %in% cl_list)
        
        if (length(cl_df[order(cl_df$p.adjust),]$ID) < 10) {
            maxlen2 <- length(cl_df[order(cl_df$p.adjust), ]$ID)
        } else {
            maxlen2 <- 10
        }
        
        final_sel <- c(final_sel, cl_df[order(cl_df$p.adjust), ][1:maxlen2, ]$ID)
    }

    # Refilter the interaction matrix with the selected terms
    mat3 <- mat[final_sel, final_sel]

    # Iter through the list and produce the network table
    net_df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("Source", "Target", "Kappa", "Name_S", "Cluster",  "Padj_S", "Count_S", "Name_T", "Cluster",  "Padj_T", "Count_T"))

    # Make a vector to avoid repeated edges
    done_v <- c()

    # Loop
    for (i in final_sel) {
        # Loop again to iter through the matrix
        for (j in  final_sel) {
            # Get the kappa value for checling
            kappa_val <- mat3[i, j]
            
            # Get the values fro the relations
            if (i == j) {
                next
            } else if (kappa_val < 0.3) {
                next
            } else if (j %in% done_v) {
                next
            } else {
                net_df[nrow(net_df) + 1, ] <- c(i, j, kappa_val, 
                    GOtable[GOtable$ID == i,]$Description, df[df$id == i,]$cluster, GOtable[GOtable$ID == i,]$p.adjust, GOtable[GOtable$ID == i,]$Count,
                    GOtable[GOtable$ID == j,]$Description, df[df$id == j,]$cluster, GOtable[GOtable$ID == j,]$p.adjust, GOtable[GOtable$ID == j,]$Count
                )
            }
        }
    
        # Add i to done terms
        done_v <- c(done_v, c(i))
    }

    # Change columns to numeric classes where it applies
    net_df$Kappa <- as.numeric(net_df$Kappa)
    net_df$Padj_S <- as.numeric(net_df$Padj_S)
    net_df$Padj_T <- as.numeric(net_df$Padj_T)
    net_df$Count_S <- as.integer(net_df$Count_S)
    net_df$Count_T <- as.integer(net_df$Count_T)

    # Rename the columns
    colnames(net_df) <- c("Source", "Target", "Kappa", "Name", "Cluster",  "Padj", "Count", "Name", "Cluster",  "Padj", "Count")

    # Write the network table
    write.table(net_df, file = gsub(".Rda", "_network.tsv", res, fixed = TRUE), sep="\t", row.names = FALSE)
}
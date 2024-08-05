## Required:
# 00 - DifferentialExpressionMouseModels
# This script generates the bubble plot for seelcted cell types

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
import math
import matplotlib
from statistics import mean
import numpy as np

def main():
    # Generate a bubble plot of established gene markers for intestinal cell types.

    # list the file names for ileum and colon
    files_col = [
        "detables/AcDSS_Healthy_expanded.tsv",
        "detables/cDSS_Healthy_expanded.tsv",
        "detables/Casp8dIECCol_Casp8floxCol_expanded.tsv",
        "detables/TC_RKO_expanded.tsv",
        "detables/OxC_Healthy_expanded.tsv",
        "detables/AcTNBS_Healthy_expanded.tsv",
        "detables/cTNBS_Healthy_expanded.tsv",
        "detables/TNFdARECol_WTCol_expanded.tsv",
        "detables/Hhepa_Healthy_expanded.tsv",
        "detables/Crode_Healthy_expanded.tsv"
          ]
    
    files_ile = [
        "detables/Casp8dIECIle_Casp8floxIle_expanded.tsv",
        "detables/TNFdAREIle_WTIle_expanded.tsv",
        "detables/Everm_Healthy_expanded.tsv"
    ]

    # Establish the names for ileum and colon
    tabnames_col = ["AcDSS","cDSS","Casp8ΔIEC", "TC","OxC", "AcTNBS","cTNBS","TNFΔARE", "Hhepa","Crode"]
    
    tabnames_ile = ["Casp8ΔIEC","TNFΔARE","Everm"]

    # Get the lists of genes
    genefiles = [
        "cell_type_markers/InfMacrophages_ensembl.txt",
        "cell_type_markers/T_cells_ensembl.txt",
        "cell_type_markers/B_cells_ensembl.txt",
        "cell_type_markers/Goblet_ensembl.txt",
        "cell_type_markers/TAProg_ensembl.txt",
        "cell_type_markers/EntericNeuron_ensembl.txt"
    ]

    # Craft the complete genelist for the plot
    genelist = []

    for genefile in genefiles:
        with open(genefile, 'r') as filehandle:
            genes = [i.rstrip() for i in filehandle.readlines()]

        genelist += genes

    # generate color map
    cvals  = [-2, 0, 2]
    colors = ['purple','#faebd7','orange']

    norm = plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm,cvals), colors))
    cm = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)

    # ------------------------------------------------------------
    # Colon

    # Initiate the tables for the values represented in the plot
    fc = list()
    pv = list()

    for file in files_col:
        # read dataframe
        df = pd.read_csv(file, sep='\t', index_col=None)
        #print(file)
        
        # filter by genenames
        df = df[df['EnsGenes'].isin(genelist)]

        # keep order
        df = df.iloc[pd.Index(df['EnsGenes']).get_indexer(genelist)]

        #df = df.sort_values(["EnsGenes"])

        # Ger genenames
        genenames = df['Genes'].values.tolist()
        #print(genenames)
        
        #QC
        #print(df['EnsGenes'].values.tolist())
        #print(len(df['log2FoldChange'].values.tolist()))
        
        fc.append(df['log2FoldChange'].values.tolist())
        
        # fix limits on adjusted p values
        pvList = [-math.log10(float(i)+1e-148) for i in df['padj'].tolist()]
        pvList = [5 if i > 5 else i for i in pvList]
        
        pv.append(pvList)

    # transform into numpy array
    fcArr = np.array(fc)
    pvArr = np.array(pv)

    # Establish dimensions of the grid
    N = len(files_col)
    M = len(genenames)

    x, y = np.meshgrid(np.arange(M), np.arange(N))

    # Create figure
    fig, ax = plt.subplots(figsize=(15,3))

    # Get the radius
    R = (pvArr/pvArr.max()/2)*0.95

    # Generate the circles with coordinates
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]

    # Color the circles with a scale
    col = PatchCollection(circles, array=fcArr.flatten(), cmap=cm)
    col.set_clim([-2, 2])

    # Add to the ax
    ax.add_collection(col)

    # Format ticks
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
        yticklabels=tabnames_col)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.set_xticklabels(genenames, rotation=90, style='italic')
    #ax.grid(which='minor')

    fig.colorbar(col)
    plt.tight_layout()
    plt.savefig('plots/CellMarker_BubblePlot_Colon.png', dpi=600)
    plt.close()

    # ------------------------------------------------------------
    # Ileum

    # Initiate the tables for the values represented in the plot
    fc = list()
    pv = list()

    for file in files_ile:
        # read dataframe
        df = pd.read_csv(file, sep='\t', index_col=None)
        #print(file)
        
        # filter by genenames
        df = df[df['EnsGenes'].isin(genelist)]

        # keep order
        df = df.iloc[pd.Index(df['EnsGenes']).get_indexer(genelist)]

        #df = df.sort_values(["EnsGenes"])

        # Ger genenames
        genenames = df['Genes'].values.tolist()
        #print(genenames)
        
        #QC
        #print(df['EnsGenes'].values.tolist())
        #print(len(df['log2FoldChange'].values.tolist()))
        
        fc.append(df['log2FoldChange'].values.tolist())
        
        # fix limits on adjusted p values
        pvList = [-math.log10(float(i)+1e-148) for i in df['padj'].tolist()]
        pvList = [5 if i > 5 else i for i in pvList]
        
        pv.append(pvList)

    # transform into numpy array
    fcArr = np.array(fc)
    pvArr = np.array(pv)

    # Establish dimensions of the grid
    N = len(files_ile)
    M = len(genenames)

    x, y = np.meshgrid(np.arange(M), np.arange(N))

    # Create figure
    fig, ax = plt.subplots(figsize=(15,1.7))

    # Get the radius
    R = (pvArr/pvArr.max()/2)*0.95

    # Generate the circles with coordinates
    circles = [plt.Circle((j,i), radius=r) for r, j, i in zip(R.flat, x.flat, y.flat)]

    # Color the circles with a scale
    col = PatchCollection(circles, array=fcArr.flatten(), cmap=cm)
    col.set_clim([-2, 2])

    # Add to the ax
    ax.add_collection(col)

    # Format ticks
    ax.set(xticks=np.arange(M), yticks=np.arange(N),
        yticklabels=tabnames_ile)
    ax.set_xticks(np.arange(M+1)-0.5, minor=True)
    ax.set_yticks(np.arange(N+1)-0.5, minor=True)
    ax.set_xticklabels(genenames, rotation=90, style='italic')
    #ax.grid(which='minor')

    fig.colorbar(col)
    plt.tight_layout()
    plt.savefig('plots/CellMarker_BubblePlot_Ileum.png', dpi=600)
    plt.close()

if __name__ == '__main__':
    main()
# 01-2 - Statistical deconvolution of colon samples - Part 2
# This script plots the deconvolution results from the colonic mouse model samples.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

def main():
    # Generate a stacked barplot from deconvolution results

    # Get the resulting files from the deconvolution assay
    filelist = [
        # Contains AcDSS, cDSS and TC
        'deconvolution/AcDSS_deconv.tsv',
        'deconvolution/Casp8dIECCol_deconv.tsv',
        # Contains OxC and Crode
        'deconvolution/OxC_deconv.tsv',
        # contains AcTNBS and cTNBS
        'deconvolution/AcTNBS_deconv.tsv',
        'deconvolution/TNFdARECol_deconv.tsv',
        'deconvolution/Hhepa_deconv.tsv'
    ]

    # Merge the tables into one

    resdf = pd.DataFrame()
    for file in filelist:
        # Transpose the table
        df = pd.read_csv(file, sep='\t').T

        if resdf.empty:
            resdf = df
        else:
            resdf = pd.concat([resdf,df])

    # Place the samples in order and assign to a model name
    samples_dict = {
        'Ctrl_DSS': ['CErldcM1', 'CErldcM2', 'CErldcM3', 'CErldcM4', 'CErldcM5'],
        'AcDSS': ['DSSdcM2', 'DSSdcM3', 'DSSdcM5'],
        'cDSS': ['cDSSdcM1', 'cDSSdcM2', 'cDSSdcM3', 'cDSSdcM4', 'cDSSdcM5'],
        'Ctrl_Casp8dIEC': ['Col_182', 'Col_184', 'Col_195', 'Col_326'], 
        'Casp8dIEC': ['Col_100', 'Col_131', 'Col_186', 'Col_87', 'Col_9'],
        'Rag1KO': ['RKOdcM1', 'RKOdcM2', 'RKOdcM3', 'RKOdcM4', 'RKOdcM5'],
        'TC': ['TCdcM1','TCdcM2', 'TCdcM3', 'TCdcM5'], 
        'Ctrl_OxC': ['Control_4', 'Control_5_1', 'Control_6_1', 'Control_7_1', 'Control_8_1'],
        'OxC': ['Ox_1_1', 'Ox_3_1', 'Ox_4', 'Ox_7', 'Ox_8'],
        'Ctrl_TNBS': ['Bl6colon71', 'Bl6colon72', 'Bl6colon73', 'Bl6colon74', 'Bl6colon75', 'Bl6colon76'],
        'AcTNBS': ['Co13V581a', 'Co14V581a', 'Co15V581a', 'Co53V582a', 'Co55V582a'],
        'cTNBS': ['Co08V381c', 'Co10V381c', 'Co11V381c'],
        'Ctrl_TNFdARE': ['C10TRR241ARE', 'C11TRR241ARE', 'C12TRR241ARE', 'C4TRR241ARE', 'C9TRR241ARE'],
        'TNFdARE': ['C13TRR241ARE', 'C14TRR241ARE', 'C15TRR241ARE', 'C5TRR241ARE', 'C6TRR241ARE'],
        'Ctrl_Crode': ['Control_4', 'Control_5_1', 'Control_6_1', 'Control_7_1', 'Control_8_1'],
        'Crode': ['Cb_1', 'Cb_3_1', 'Cb_4_1', 'Cb_5_1', 'Cb_6'],
        'Ctrl_Hhepa': ['a19', 'a20', 'a9', 'a10'], 
        'Hhepa': ['a15', 'a16', 'a17', 'a18']
    }

    # Fix the names of the cell types, select the necessary ones and rename ISlC more appropriately
    resdf = resdf[['Colonocytes','Colonocytes-EEC','Colonocytes-Secretory','Keratinocytes',
               'T_cells','B_cells','Macrophages_APCS','Neutrophils',
               'RBCs','Fibroblasts','Lymphatics','Endothelial_cells','Muscle-like_cells']]
    resdf = resdf.rename(columns={'Keratinocytes': 'Intermediate squamous-like cells'})

    # Establish the colors for the plot
    colors = [
        '#fdd9b4','#fda762','#f3701b','#c54102',
        '#e5ccff','#b266ff','#9933ff','#7f00ff',
        '#dbf1d5','#aedea7','#74c476','#37a055','#0c7734'
        ]
    
    # Get the cell types to a list
    cell_types = resdf.columns.tolist()

    # Loop to get the mean and error per model and cell type
    means_list = []
    error_list = []

    for modelID in samples_dict:
        # Start tmp lists
        means_row = []
        error_row = []

        # Get the row per category
        resdf_tmp = resdf.loc[samples_dict[modelID]]

        # Loop through cell types and get mean and error
        for ctype in cell_types:
            means_row.append(resdf_tmp[ctype].mean())
            error_row.append(resdf_tmp[ctype].sem())

        # Add model name
        means_row = means_row + [modelID]
        error_row = error_row

        # Add to general list
        means_list.append(means_row)
        error_list.append(error_row)

    # Make into dataframe
    means_df = pd.DataFrame(means_list)
    means_df.columns = cell_types + ['Model']
    error_df = pd.DataFrame(error_list)
    error_df.columns = cell_types

    # Arrange the positions for the bars and size
    N = len(means_df.index)
    ind = np.arange(N)
    ind = [0,1,2, 3.5,4.5, 6,7, 8.5,9.5, 11,12,13, 14.5,15.5, 17,18, 19.5,20.5]
    width = 0.9

    # Create the figure canvas
    fig, ax = plt.subplots(figsize=(15, 6))

    # Set x tick labels position, text and rotation
    ax.xaxis.set_major_locator(ticker.FixedLocator((ind)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((list(samples_dict.keys()))))

    plt.setp(ax.get_xticklabels(), rotation=90)

    # plot the first set of stacked bars
    p = plt.bar(ind, means_df[cell_types[0]].tolist(), width, yerr = error_df[cell_types[0]].tolist(), error_kw=dict(ecolor='black', lw=1, capsize=3, capthick=1), color=colors[0])

    bottom = means_df[cell_types[0]].tolist()

    # Iter through the cell types to plot the rest of the stacks
    for i in range(1,len(cell_types)):
        p = plt.bar(ind, means_df[cell_types[i]].tolist(), width, bottom=bottom, yerr=error_df[cell_types[i]].tolist(), error_kw=dict(ecolor='black', lw=1, capsize=3, capthick=1), color=colors[i])

        bottom = [sum(x) for x in zip(bottom, means_df[cell_types[i]].tolist())]

    # Get arial for the legend labels
    font = font_manager.FontProperties(family='Arial', style='normal')

    # Craft a custom legend and add it to the plot
    legend_elements = [
        mpatches.Patch(color='#fdd9b4', alpha=1, label='Colonocytes'),
        mpatches.Patch(color='#fda762', alpha=1, label='Colonocytes-EEC'),
        mpatches.Patch(color='#f3701b', alpha=1, label='Colonocytes-Secretory'),
        mpatches.Patch(color='#c54102', alpha=1, label='Intermediate squamous-like cells'),
        mpatches.Patch(color='#e5ccff', alpha=1, label='T_cells'),
        mpatches.Patch(color='#b266ff', alpha=1, label='B_cells'),
        mpatches.Patch(color='#9933ff', alpha=1, label='Macrophages_APCS'),
        mpatches.Patch(color='#7f00ff', alpha=1, label='Neutrophils'),
        mpatches.Patch(color='#dbf1d5', alpha=1, label='RBCs'),
        mpatches.Patch(color='#aedea7', alpha=1, label='Fibroblasts'),
        mpatches.Patch(color='#74c476', alpha=1, label='Lymphatics'),
        mpatches.Patch(color='#37a055', alpha=1, label='Endothelial_cells'),
        mpatches.Patch(color='#0c7734', alpha=1, label='Muscle-like_cells')
        ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1, 1), fontsize="10", prop=font)

    # Save as svg and png
    plt.tight_layout()
    fig.savefig('deconvolution/Deconv_colon.svg', dpi=600)
    fig.savefig('deconvolution/Deconv_colon.png', dpi=600)


if __name__ == '__main__':
    main()
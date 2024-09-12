## Required:
# 00 - DifferentialExpressionMouseModels
# This script generates the venn diagrams and exports the genelists of common and unique gene sets

import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from venn import generate_petal_labels, draw_venn, generate_colors

def vennPlot(venn_dict, keys, label, sec_lab):
    # Fixture with only the common and exclusive numbers
    petal_labels = generate_petal_labels(venn_dict.values(), fmt="{size}")

    petal_labels_filt = {key: petal_labels[key] for key in keys}

    axUp = draw_venn(
        petal_labels=petal_labels_filt, dataset_labels=venn_dict.keys(),
        hint_hidden=False, colors=generate_colors(cmap='plasma', n_colors=5), ax=None, 
        fontsize=20, legend_loc=None, 
        figsize = (7, 7)
    )
    axUp.set_title(f'Venn_{sec_lab}', fontdict = {'fontsize' : 18})
    figUp = axUp.get_figure()
    plt.tight_layout()
    figUp.savefig(f"plots/{label}_{sec_lab}_Venn.svg", dpi=600)
    figUp.savefig(f"plots/{label}_{sec_lab}_Venn.png", dpi=600)

    petal_labels = generate_petal_labels(venn_dict.values(), fmt="{size}")

def vennDiag(models, label, keys):
    # Create the dictionaries for the lists and the venn plots

    # Venn source
    venn_dict_up = {}
    venn_dict_dw = {}

    # Common genes
    common_genes_up = []
    common_genes_dw = []

    # Exclusive genes
    excl_genes_up = {}
    excl_genes_dw = {}

    #iter through the experiments to fill the lists up with stuff

    for model in models:
        # Read the file and filter by significance
        df = pd.read_csv(models[model], sep='\t', index_col=None)
        df = df[(df['padj'] < 0.05) & (df['FLAG'] == 'OK')]
        
        # Add the genes up and down to the exclusive list of genes for later processing
        excl_genes_up[model] = df[df['log2FoldChange'] > 0]['EnsGenes'].tolist()
        excl_genes_dw[model] = df[df['log2FoldChange'] < 0]['EnsGenes'].tolist()
        
        # Same thing for the venn dictionaries, as 
        venn_dict_up[model] = set(df[df['log2FoldChange'] > 0]['EnsGenes'].tolist())
        venn_dict_dw[model] = set(df[df['log2FoldChange'] < 0]['EnsGenes'].tolist())
        
        # Set up the common genes up n down
        if not common_genes_up:
            common_genes_up = df[df['log2FoldChange'] > 0]['EnsGenes'].tolist()
            common_genes_dw = df[df['log2FoldChange'] < -0]['EnsGenes'].tolist()
        else:
            common_genes_up = list(set(common_genes_up).intersection(df[df['log2FoldChange'] > 0]['EnsGenes'].tolist()))
            common_genes_dw = list(set(common_genes_dw).intersection(df[df['log2FoldChange'] < 0]['EnsGenes'].tolist()))
            
    # Save the common genes
    with open(f'VennGeneLists/{label}_common_genes_up_ensembl.txt', 'w') as com_up:
        com_up.write('\n'.join(common_genes_up))
    com_up.close()
    with open(f'VennGeneLists/{label}_common_genes_dw_ensembl.txt', 'w') as com_dw:
        com_dw.write('\n'.join(common_genes_dw))
    com_dw.close()

    # Iter through the exclusive genes lists to get the ones that appear only in each experiment
    for base_model in excl_genes_up:
        # Get list of current genes
        excl_genes_up_model = excl_genes_up[base_model]
        
        # Iter again to compare one against another
        for problem_model in excl_genes_up:
            # If the model is the same, continue, if not do the thingy
            if base_model == problem_model:
                continue
            else:
                # As seen i https://stackoverflow.com/questions/41125909/find-elements-in-one-list-that-are-not-in-the-other
                excl_genes_up_model = list(set(excl_genes_up_model) - set(excl_genes_up[problem_model]))
            
        # Save the generated list
        with open(f'VennGeneLists/{label}_{base_model}_excl_up_ensembl.txt', 'w') as excl_up:
            excl_up.write('\n'.join(excl_genes_up_model))
        excl_up.close()

    # Repeat for the down genes
    for base_model in excl_genes_dw:
        # Get list of current genes
        excl_genes_dw_model = excl_genes_dw[base_model]
        
        # Iter again to compare one against another
        for problem_model in excl_genes_dw:
            # If the model is the same, continue, if not do the thingy
            if base_model == problem_model:
                continue
            else:
                # As seen i https://stackoverflow.com/questions/41125909/find-elements-in-one-list-that-are-not-in-the-other
                excl_genes_dw_model = list(set(excl_genes_dw_model) - set(excl_genes_dw[problem_model]))
            
        # Save the generated list
        with open(f'VennGeneLists/{label}_{base_model}_excl_dw_ensembl.txt', 'w') as excl_dw:
            excl_dw.write('\n'.join(excl_genes_dw_model))
        excl_dw.close()

    # Generate the common diagrams
    vennPlot(venn_dict_up, keys['common'], label, 'common_up')
    vennPlot(venn_dict_dw, keys['common'], label, 'common_dw')

    # Generate the exclusive venn diagrams
    vennPlot(venn_dict_up, keys['excl'], label, 'excl_up')
    vennPlot(venn_dict_dw, keys['excl'], label, 'excl_dw')
    

def main():
    # Get the dicts of models
    models_col = {
        'AcDSS':'detables/AcDSS_Healthy_expanded.tsv',
        'OxC':'detables/OxC_Healthy_expanded.tsv',
        'TC':'detables/TC_RKO_expanded.tsv',
        'Hhepa':'detables/Hhepa_Healthy_expanded.tsv',
        'Crode':'detables/Crode_Healthy_expanded.tsv'
    }

    models_ile = {
        'Casp8ΔIEC':'detables/Casp8dIECIle_Casp8floxIle_expanded.tsv',
        'TNFΔARE':'detables/TNFdAREIle_WTIle_expanded.tsv',
        'Everm':'detables/Everm_Healthy_expanded.tsv'
    }

    keys_col = {
        'common': ['11111'],
        'excl': ['00001', '00010', '00100', '01000', '10000']
    }

    keys_ile = {
        'common': ['111'],
        'excl': ['100','010','001']
    }

    vennDiag(models_col, 'colon', keys_col)
    vennDiag(models_ile, 'ileum', keys_ile)


if __name__ == '__main__':
    main()
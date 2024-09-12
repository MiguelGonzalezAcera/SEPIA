## Required:
# 00 - DifferentialExpressionMouseModels
# 05-1 - VennDiagrams
# 05-2 - GOofVennLists
# This script generates the Lollipop plots from the results of GO ORA

import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
import glob
import textwrap

def main():
    # Get the file lists
    filelists = glob.glob('VennGOresult/*_BP.tsv')

    # Iter through
    for file in filelists:
        # Read the file
        df = pd.read_csv(file, sep='\t')

        # Select only significant lines
        df = df[df['p.adjust'] < 0.05]

        # Wrapt text labels for them to fit
        df['Description'] = ['\n'.join(textwrap.wrap(i, 20)) for i in df['Description'].tolist()]

        # Select only the first three
        df = df.head(3)

        # Create the figure
        fig, ax = plt.subplots(figsize=(4.5, 3.5))

        # Run under a try-except in case the file is empty
        try:
            # Plot the barplots (base of the lolipop)
            ax = df.plot.barh(x='Description', y='Count', width=0.05, ax=ax, color='black', alpha=0.3, legend=False)

            # Select the color for the circles
            if "_up_" in file:
                color_candy = 'orange'
            else:
                color_candy = 'purple'

            # Draw the circles (head of the lolipop)
            ax = ax.scatter(df['Count'], df['Description'], s=df['Count']*5, edgecolor='black', color=color_candy)

            # Title the plot
            plt.title(file.split('/')[-1].replace('_BP.tsv', ''), fontdict = {'fontsize' : 12})

            # Manage ticks, fonts and margins
            plt.autoscale(enable=True, axis='y', tight=False)
            plt.xlabel('')
            plt.ylabel('')
            plt.yticks(fontsize=14)
            fig.tight_layout()
            plt.margins(0.2,0.25)

            # Save the figure
            fig.savefig(file.replace('.tsv','.png'))
    
        except:
            # In case a file is empty, state which
            print(file)
            print('failed')

if __name__ == '__main__':
    main()
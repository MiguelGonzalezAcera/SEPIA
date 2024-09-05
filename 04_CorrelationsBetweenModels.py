## Required:
# 00 - DifferentialExpressionMouseModels
# This script generates the correlation plot for all models, segregated by ileum and colon

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math

import matplotlib as mp
import numpy as np


# functions for the correlation plot in the lower half
# from https://stackoverflow.com/questions/27164114/show-confidence-limits-and-prediction-limits-in-scatter-plot
def plot_ci_manual(t, s_err, n, x, x2, y2, color_used, ax=None):
    """Return an axes of confidence bands using a simple approach.
    
    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}
    
    References
    ----------
    .. [1] M. Duarte.  "Curve fitting," Jupyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb
    
    """
    if ax is None:
        ax = plt.gca()
    
    ci = t * s_err * np.sqrt(1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    ax.fill_between(x=x2, y1=y2 + ci, y2=y2 - ci, color=color_used)

    return ax

# Computations ----------------------------------------------------------------    
# Modeling with Numpy
def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b) 

#---------------------------------------

# function to plot the fold change and the regression, IBDome style baybee
def IBDomelinreg(x, y, sign_genes,  **kwargs):
    # establish cmap
    cmap = mp.colors.LinearSegmentedColormap.from_list("", ["purple","#faebd7","#FF8000"])
    norm = mp.colors.Normalize(vmin=-1.5, vmax=1.5)
    
    # Get list of exclusive genes
    sign_genes_list = sign_genes[x.name][y.name]
    
    # Filter by significant genes
    x = x[x.index.isin(sign_genes_list)]
    y = y[y.index.isin(sign_genes_list)]
    
    # perform the linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(x),y)
    #print(x.name, y.name, slope)
    
    # get the linear regression
    x = np.array(x)
    y = np.array(y)
    p, cov = np.polyfit(x, y, 1, cov=True)                     # parameters and covariance from of the fit of 1-D polynom.
    y_model = equation(p, x)                                   # model using the fit parameters; NOTE: parameters here are coefficients

    # Statistics
    n = x.size                                           # number of observations
    m = p.size                                                 # number of parameters
    dof = n - m                                                # degrees of freedom
    t = stats.t.ppf(0.975, n - m)                              # t-statistic; used for CI and PI bands

    # Estimates of Error in Data/Model
    resid = y - y_model                                        # residuals; diff. actual data from predicted values
    chi2 = np.sum((resid / y_model)**2)                        # chi-squared; estimates error in data
    chi2_red = chi2 / dof                                      # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid**2) / dof)                    # standard deviation of the error

    # Plotting --------------------------------------------------------------------
    ax = plt.gca()
    plt.ylim(min(y)-1,max(y)+1)
    diff_val = np.max(x)-np.min(x)
    plt.xlim(min(x)-1,max(x)+1)

    
    # Draw axis
    ax.axvline(x=0, c='#D1D1D1', alpha=0.5)
    ax.axhline(y=0, c='#D1D1D1', alpha=0.5)
    
    # Data
    ax.plot(
        x, y, ".", color=cmap(norm(slope), alpha=0.5), markersize=6, 
        markeredgewidth=1, markeredgecolor=cmap(norm(slope)), markerfacecolor=cmap(norm(slope), alpha=0.5)
    )

    # Fit
    ax.plot(x, y_model, "-", color="0.1", linewidth=2.5, alpha=0.5, label="Fit")  

    x2 = np.linspace(-30,30, 100)
    y2 = equation(p, x2)

    # Confidence Interval (select one)
    plot_ci_manual(t, s_err, n, x, x2, y2, cmap(norm(slope), alpha=0.3), ax=ax)

    # Prediction Interval
    #pi = t * s_err * np.sqrt(1 + 1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))   
    #ax.fill_between(x2, y2 + pi, y2 - pi, color="None", linestyle="--")
    #ax.plot(x2, y2 - pi, "--", color="0.5", label="95% Prediction Limits")
    #ax.plot(x2, y2 + pi, "--", color="0.5")
    # ax.text(min(tmp_mt_ibdome['log_tpm']),max(tmp_mt_ibdome['median_modified_riley_score']), f"y = {slope:.2f}x + {intercept:.2f}. Pval = {p_value:.2f}. R2 = {r_value**2:.2f}")

# Function to plot the R squared of each regression
def corrdot(x, y, sign_genes, **kwargs):
    #print(args[1])
    
    # establish cmap
    cmap = mp.colors.LinearSegmentedColormap.from_list("", ["#faebd7","purple"])
                                                            #"#FF8000"])
    norm = mp.colors.Normalize(vmin=-1.5, vmax=1.5)
    
    # Get list of exclusive genes
    sign_genes_list = sign_genes[x.name][y.name]
    
    # Filter by significant genes
    x = x[x.index.isin(sign_genes_list)]
    y = y[y.index.isin(sign_genes_list)]
    
    # perform the linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(x),y)
    
    # get the R value as text
    corr_text = f"{r_value**2:2.2f}".replace("0.", ".")
    
    # Make the actual plot
    ax = plt.gca()
    ax.set_axis_off()
    marker_size = abs(r_value**2) * 5000
    ax.scatter([.5], [.5], marker_size, [r_value**2], alpha=0.75, cmap=cmap,
               vmin=0, vmax=1, transform=ax.transAxes)
    font_size = abs(r_value**2) * 30 + 5
    ax.annotate(corr_text, [.5, .5,],  xycoords="axes fraction",
                ha='center', va='center', fontsize=font_size)
    

# Function to put the name of the model in the diagonal
def find_the_x(x, sign_genes, absmax, **kwargs):
    
    # Get the center of the plot
    xlist = x[x.index.isin(sign_genes[x.name][x.name])].tolist()
    c = (max(xlist) - min(xlist))/2 + min(xlist)
    
    # proportion of size against max
    p = sign_genes[x.name]['pie']['size']/absmax
    
    # Get the radius
    # 1st get the area of the max circle in the allocated area, so maxarea = math.pi*((max(xlist) - c)**2)
    # 2nd Proportionate the area to the size we want using the p, usedarea = maxarea*p
    # 3rd Get the radius from the area using the reverse equation, r = sqrt(usedarea/math.pi)
    r = np.sqrt(((math.pi*((max(xlist) - c)**2))*p)/math.pi)
    #r = (max(xlist) - c)*p
    # print(x.name, sign_genes[x.name]['pie']['size'])
    # print(f'radius: {r}')
    
    #print(x.name, max(xlist), min(xlist), c, sign_genes[x.name]['pie']['size'])
    # Set the name of the model
    model_name = x.name
    
    # Make the plot
    ax = plt.gca()

    ax.pie(
        sign_genes[x.name]['pie']['genes'],
        center=(c, c), radius=r,
        colors=['#FF8000','purple'], frame = True
    )
    
    plt.xlim(min(xlist)-1,max(xlist)+1)
    plt.ylim(min(xlist)-1,max(xlist)+1)
    
    #ax.set_axis_off()

def corrplot(mod_list, label, dim=7):
    # Step 1, get a nested dict of the lists of genes that are significant between two models
    sign_genes_dict = {}

    # Start the absolute max value
    absmax = 0

    # Loop 1
    for model1 in mod_list:
        # Add first model to list
        sign_genes_dict[model1] = {}
        
        # Read the table for model1 and filter by significance and flag
        model1df = pd.read_csv(mod_list[model1], sep='\t', index_col=None)
        #model1df = model1df[(model1df['padj'] < 0.05) & (model1df['FLAG'] == "OK")]
        model1df = model1df[model1df['padj'] < 0.05]
        
        # Get the number of genes up n down and size
        sign_genes_dict[model1]['pie'] = {
            'genes': [
                len(model1df[(model1df['padj'] < 0.05) & (model1df['log2FoldChange'] > 0)]['padj'].tolist()),
                len(model1df[(model1df['padj'] < 0.05) & (model1df['log2FoldChange'] < 0)]['padj'].tolist()),
            ],
            'size': len(model1df[model1df['padj'] < 0.05]['padj'].tolist())
        }
        
        # Update the absmax
        if len(model1df[model1df['padj'] < 0.05]['padj'].tolist()) > absmax:
            absmax = len(model1df[model1df['padj'] < 0.05]['padj'].tolist())
        
        # Loop 2
        for model2 in mod_list:        
            # Read the table for model2 and filter by significance and flag
            model2df = pd.read_csv(mod_list[model2], sep='\t', index_col=None)
            #model2df = model2df[(model2df['padj'] < 0.05) & (model2df['FLAG'] == "OK")]
            model2df = model2df[model2df['padj'] < 0.05]
            
            # Get a list of the common genes contained in both tables
            comgenes = list(set(model1df['EnsGenes'].tolist()).intersection(model2df['EnsGenes'].tolist()))
            #print(model1, model2, len(comgenes))
            
            # Save the list in the dictionary
            sign_genes_dict[model1][model2] = comgenes
            
    # Add some chunk to absmax to avoid collision
    absmax = absmax * 1.1

    # Step 2: Get the fold changes of each model merged in a single data frame, putting zeros in place of NAs
    # Init the result
    FCdf = pd.DataFrame()

    # Loop through
    for model in mod_list:
        # Read the model and get the gene id and fold change
        modelDf = pd.read_csv(mod_list[model], sep='\t', index_col=None)
            
        modelDf = modelDf[['EnsGenes', 'log2FoldChange']]
        
        # change column names to model name
        modelDf.columns = ['EnsGenes',model]
        
        # Merge into result
        if FCdf.empty:
            FCdf = modelDf
        else:
            FCdf = FCdf.merge(modelDf, on='EnsGenes', how='outer').fillna(0)
            
    # Set ensGenes as index
    FCdf.set_index('EnsGenes', inplace=True)

    # Step 3: Define the functions and make the plot
    # Set the aesthetics in seaborn
    sns.set_theme(font='Arial', style='white', font_scale=1.6)

    # Make the frid in question
    g = sns.PairGrid(FCdf, aspect=1, diag_sharey=False)
    g.map_lower(IBDomelinreg, sign_genes = sign_genes_dict)
    g.map_upper(corrdot, sign_genes = sign_genes_dict)
    g.map_diag(find_the_x, sign_genes = sign_genes_dict, absmax = absmax)
    g.figure.set_size_inches(dim,dim)

    g.savefig(f'plots/{label}_corr_plot.svg', dpi=600)
    g.savefig(f'plots/{label}_corr_plot.png', dpi=600)


def main():
    # Step 0: Define a dictionary with our input files and their location
    models_mouse_col = {
        'AcDSS':'detables/AcDSS_Healthy_expanded.tsv',
        'cDSS':'detables/cDSS_Healthy_expanded.tsv',
        'Casp8ΔIEC':'detables/Casp8dIECCol_Casp8floxCol_expanded.tsv',
        'OxC':'detables/OxC_Healthy_expanded.tsv',
        'AcTNBS':'detables/AcTNBS_Healthy_expanded.tsv',
        'cTNBS':'detables/cTNBS_Healthy_expanded.tsv',
        'TC':'detables/TC_RKO_expanded.tsv',
        'TNFΔARE':'detables/TNFdARECol_WTCol_expanded.tsv',
        'Hhepa':'detables/Hhepa_Healthy_expanded.tsv',
        'Crode':'detables/Crode_Healthy_expanded.tsv'
    }

    models_mouse_ile = {
        'Casp8ΔIEC':'detables/Casp8dIECIle_Casp8floxIle_expanded.tsv',
        'TNFΔARE':'detables/TNFdAREIle_WTIle_expanded.tsv',
        'Everm':'detables/Everm_Healthy_expanded.tsv'
    }

    # Run each throuth the function
    corrplot(models_mouse_col, 'colon', dim=14)
    corrplot(models_mouse_ile, 'ileum', dim=7)


if __name__ == '__main__':
    main()
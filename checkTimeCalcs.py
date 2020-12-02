import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import linregress

# generate summary information for raw data frame
def generate_summary(df):
    d = []
    # get information
    genotypes = df.Genotype.unique()
    times = df.Firing_time.unique()
    # process data
    for g in genotypes:
        da = df[df.Genotype == g]
        da['Ratio'] = da['Leading']/da['Lagging']
        data = da.groupby(['Firing_time'], as_index=False).agg({'Leading':'sum', \
                'Lagging':'sum', 'Ratio':'median'}).reset_index()
        data['Genotype'] = g
        data['Libraries'] = len(da.Library.unique())
        data['MLE_ratio'] = data['Leading']/data['Lagging']
        data['log_MLE_ratio'] = np.log(data['MLE_ratio'])
        data['log_ratio'] = np.log(data['Ratio'])
        d.append(data)
    data = pd.concat(d)
    return data

# draw scatter plot for ratio
def draw_ratio_scatter(df, genotypes,output=None, use_MLE_ratio=True, logrithm=True, use_efficiency=False):
    palette = sns.hls_palette(16, l=0.5, s=1)
    palette = [palette[0], palette[14], '#003f3f','#000000']
    df = df[df.Genotype.isin(genotypes)]
    sns.set(style='ticks')
    fig, ax = plt.subplots(figsize=(10,8))
    if logrithm:
        feature = 'log_MLE_ratio' if use_MLE_ratio else 'log_ratio'
    else:
        feature = 'MLE_ratio' if use_MLE_ratio else 'Ratio'
    palette = palette[:len(df.Genotype.unique())]
    sns.scatterplot(x='Firing_time', y=feature, hue='Genotype', hue_order=genotypes, data=df, palette=palette, ax=ax)
    if use_efficiency:
        xlim = [0,1]
    else:
        xlim = [15,40]
    # calculate regression line
    ws = {}
    bs = {}
    r2s = {}
    c=0
    title = 'Leading/lagging ratio with ARS firing time'
    for g in genotypes:
        da = df[df.Genotype==g]
        times = da.Firing_time.values
        ratios = da[feature].values
        # remove nan
        times = times[~np.isnan(ratios)]
        ratios = ratios[~np.isnan(ratios)]
        # remove inf
        ratios[ratios == np.inf] = np.max(ratios[ratios!=np.inf])
        ratios[ratios == 0] = np.min(ratios[ratios!=0])
        w,b,r2,_,_ = linregress(times, ratios)
        plt.plot(xlim,[xlim[0]*w+b, xlim[1]*w+b], c=palette[c], linewidth=3)
        # store infomation
        title += '\n{}: Coefficient={:.4f}, bias={:.4f}, R-square={:.4f}'.format(g, w,b, r2)
        c += 1
        ws[g] = w
        bs[g] = b
        r2s[g] = r2
    # formatting for publication
    sns.despine()
    plt.subplots_adjust(left=0.12, right=0.96, bottom=0.06)
    ax.set_xlim(xlim)
    if logrithm:
        ax.set_ylim([-1, 1])
    else:
        ax.set_ylim([0, 2])
    plt.locator_params(axis='y', nbins=5)
    plt.suptitle(title)
    plt.ylabel('')
    plt.xlabel('')
    ax.get_legend().remove()
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    fig.savefig(output)


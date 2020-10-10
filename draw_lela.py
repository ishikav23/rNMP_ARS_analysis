#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import linear_model
from matplotlib.ticker import FuncFormatter


def main():

    # argparse
    parser = argparse.ArgumentParser(description='Draw leading/lagging comparison plot for ars comparison')
    parser.add_argument('ars', type=argparse.FileType('r'), help='Ars region frequency file')
    parser.add_argument('-m', type=int, default=100, help='Minimum number of ribose as threshold, default=100')
    parser.add_argument('-o', default='', help='Output file basename')
    args = parser.parse_args()

    if args.o == '':
        args.o = args.ars.name.split('.')[0]

    # parameters
    FREQ_COL_NUM = 9

    # get information for bed file
    data = []
    columns = args.ars.readline().rstrip('\n').split('\t')

    for l in args.ars:
        ws = l.rstrip('\n').split()
        data.append(ws[:5] + [int(ws[5]), ws[6], int(ws[7]),float(ws[8])])
    df = pd.DataFrame(data, columns = columns[:FREQ_COL_NUM])

    # set Categorical data
    lib_params = {'Genotype':['WT','Rrnh201','EMrnh201','HYrnh201'],\
            'RE':['RE1','RE2','RE3'],\
            'Time': ['early', 'medium', 'late'],\
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C'] }
    df = df[df.Genotype.isin(lib_params['Genotype'])]
    for w in ['Genotype','RE','Time','String']:
        categorical_order = []
        unique_values = df[w].unique()
        for v in lib_params[w]:
            if v in unique_values:
                categorical_order.append(v)
        for v in unique_values:
            if v not in categorical_order:
                categorical_order.append(v)
        df[w] = pd.Categorical(df[w], categorical_order)

    # get params
    times = df.Time.unique()
    flanks = df.Flank.unique()
    la,le = sorted(df.Strand.unique())

    # set param needed for plot
    genotype_needed =  ['WT','Rrnh201','EMrnh201','HYrnh201']
    flanks_needed = [15000]

    # set percentage formatter
    f_percentage = FuncFormatter(lambda y, _: '{:.0%}'.format(y))
    # formatter for scientific
    def my_scientific_format(y):
        a = f'{y:.1E}'
        m, e = a.split('E')
        e = e.replace('+0','')
        e = e.replace('-0','-')
        return f'{m}E{e}'
    f_scientific = FuncFormatter(lambda y, _: my_scientific_format(y))

    # set color
    clist = sns.hls_palette(8, l=0.5, s=1)
    clist2 = ['#a570f3','#1bce77']

    for geno in genotype_needed:
        # scatter plot for time
        for f in flanks_needed:
            da = df[(df.Genotype == geno) & (df.Flank == f) & (df.Strand == le)]
            db = df[(df.Genotype == geno) & (df.Flank == f) & (df.Strand == la)]
            dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype','RE','Time','Flank'])

            if len(dc) == 0:
                continue

            color = {}
            for i in range(len(times)):
                color[times[i]] = clist2[i]

            ds = {}
            regrs = {}
            for t in times:
                ds[t]=dc[dc.Time==t]
                # linear regression
                regrs[t]=linear_model.LinearRegression(fit_intercept=False)
                regrs[t].fit(ds[t].Sum_lagging.values.reshape(-1,1), ds[t].Sum_leading.values)

            fig, ax = plt.subplots(figsize=(8,7))
            sns.scatterplot(x='Sum_lagging', y='Sum_leading', hue='Time', data=dc, palette=color, s=100)
            # set same lim
            lim_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
            ax.set_ylim([0,lim_max])
            ax.set_xlim([0,lim_max])
            # draw regr lines
            for t in times:
                plt.plot([0, lim_max], regrs[t].predict([[0], [lim_max]]),color=color[t], linewidth=1.5)
            # draw reference line
            plt.plot([0, lim_max], [0,lim_max], linewidth='1.5', color='black', linestyle='--')
            # set y tick lables
            ax.yaxis.set_major_formatter(f_scientific)
            ax.xaxis.set_major_formatter(f_scientific)
            # remove right and top line
            sns.despine()
            # legend
            ax.legend(frameon=False, prop={'size':15})
            # fontsize
            ax.set_ylabel('Counts of leading strand', fontsize=24)
            ax.set_xlabel('Counts of lagging strand', fontsize=24)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(18)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(18)
                tick.label.set_rotation(45)
            plt.suptitle('Ribose incorporation in leading/lagging strand ({}, flank = {:d}bp, N={:d})\n'.format(geno, f, len(dc.Library.unique())) + \
                    ', '.join([' slope_{}={:.4f}'.format(t, regrs[t].coef_[0]) for t in times]))
            # for publication
            plt.subplots_adjust(left=0.11, right=0.98, bottom=0.11)
            plt.xlabel('')
            plt.ylabel('')
            ax.get_legend().remove()
            plt.savefig(args.o + '_time_{}_{}_scatter.png'.format(geno, f))
            plt.close('all')
    print('Scatter plots for time generated!')

    # bar plot for average
    sns.set(style='ticks', font_scale=3)
    da = df[(df.Genotype.isin(genotype_needed)) & (df.Flank.isin(flanks_needed)) & (df.Strand == le)]
    db = df[(df.Genotype.isin(genotype_needed)) & (df.Flank.isin(flanks_needed)) & (df.Strand == la)]
    dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype','RE','Time','Flank'])
    dc['Ratio'] = dc['Sum_leading']/dc['Sum_lagging']
    fig, ax = plt.subplots(figsize=(10,6))
    plt.subplots_adjust(left=0.1, top=0.98, right=0.95, bottom=0.1)
    sns.barplot(x='Genotype', y='Ratio', hue='Time', data=dc, ci='sd', palette=clist2,\
             errwidth=3, capsize=0.18, edgecolor="white", ax=ax)
    # add data points
    sns.swarmplot(x='Genotype', y='Ratio', hue='Time', data=dc, dodge=True, color='black', size=7, ax=ax)
    # formatting
    labels = ax.get_xticklabels()
    ax.set_xticklabels(['WT','R','EM','HY'])
    sns.despine()
    plt.xlabel('')
    plt.ylabel('')
    ax.get_legend().remove()
    fig.savefig(f'{args.o}_bar_15000.png')
    plt.close('all')
    print('Bar plot generated!')


    print('Done!')

if __name__ == '__main__':
    main()




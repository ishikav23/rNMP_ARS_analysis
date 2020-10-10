#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu

# read data
def read_data(leading_file, lagging_file):
    data = []
    header = leading_file.readline().rstrip().split('\t')
    header_check = lagging_file.readline().rstrip().split('\t')
    assert header == header_check, 'Leading and lagging file should have same header'
    for fr, st in [[leading_file, 'Leading'],[lagging_file, 'Lagging']]:
        for l in fr:
            ws = l.rstrip().split('\t')
            if len(ws) != len(header):
                continue
            for i in range(1, len(ws)):
                data.append(ws[0].split('-') + [st, header[i], float(ws[i])])
    df = pd.DataFrame(data)
    df.columns = ['Library','String','Genotype','RE','Strand','Feature', 'Value']
    return df


# draw box plot
def draw(name, plotname, df):
    sns.set(style='ticks')
    clist = sns.hls_palette(50, l=0.5, s=1)
    cpalette = [clist[48], clist[26]]
    # check di or mono
    di = False
    if len(df.Feature.unique()) == 16:
        di = True
    if di:
        labels = ['AA','CA','GA','TA',
                'AC','CC','GC','TC',
                'AG','CG','GG','TG',
                'AU','CU','GU','TU']
        fig, ax = plt.subplots(figsize=(10,6))
        fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.1)
    else:
        labels = ['A','C','G','U']
        fig, ax = plt.subplots(figsize=(6,6))
        fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.1)
    sns.boxplot(x='Feature', y='Value', hue='Strand', data=df, ax=ax, palette=cpalette)
    sns.swarmplot(x='Feature', y='Value', hue='Strand', data=df, ax=ax, color='.25', dodge=True, size=3)
    # add statistical annotation
    if di:
        if name == 'pol3':
            box_pairs = [((x, 'Leading'), (x, 'Lagging')) for x in df.Feature.unique()[1::4]]
        elif name == 'pol2':
            box_pairs = [((x, 'Leading'), (x, 'Lagging')) for x in df.Feature.unique()[::4]]
        else:
            box_pairs = [((x, 'Leading'), (x, 'Lagging')) for x in df.Feature.unique()[::4]] + \
                    [((x, 'Leading'), (x, 'Lagging')) for x in df.Feature.unique()[1::4]]
    else:
        box_pairs = [((x, 'Leading'), (x, 'Lagging')) for x in df.Feature.unique()]
    # One-sided MWW
    pvalues = []
    for pair in box_pairs:
        data1 = df[df.Feature == pair[0][0]].groupby('Strand')['Value'].get_group(pair[0][1])
        data2 = df[df.Feature == pair[0][0]].groupby('Strand')['Value'].get_group(pair[1][1])
        if data1.median() < data2.median():
            alt = 'less'
        else:
            alt = 'greater'
        stats, p = mannwhitneyu(data1, data2, alternative=alt, use_continuity=False)
        pvalues.append(p)
    ax, results = add_stat_annotation(ax, data=df, x='Feature', y='Value', hue='Strand',\
            box_pairs=box_pairs, perform_stat_test=False, pvalues=pvalues, loc='outside',\
            fontsize='xx-large', linewidth=2, verbose=2)

    sns.despine()

    # formatting
    if di:
        if ax.get_ylim()[1] <= 0.7:
            plt.ylim((0,0.7))
        else:
            plt.ylim((0,0.9))
    else:
        plt.ylim((0,0.7))
    ax.get_legend().remove()
    plt.xlabel('')
    plt.ylabel('')
    ax.set_xticklabels(labels, fontsize=24)
    plt.yticks(fontsize=24)
    plt.savefig(plotname)
    plt.close('all')


def main():
    parser = argparse.ArgumentParser(description='Generate boxplot for heatmap')
    parser.add_argument('leading', type=argparse.FileType('r'), help='Normalized leading file')
    parser.add_argument('lagging', type=argparse.FileType('r'), help='Normalized lagging file')
    parser.add_argument('-e', nargs='+', default=[], help='Libraries to be excluded')
    parser.add_argument('-o', help='Output basename')
    args = parser.parse_args()

    if not args.o:
        args.o = 'box_plot'

    # read leading and lagging
    df = read_data(args.leading, args.lagging)

    # groups
    groups = {'WT':['WT'],
            'PolWT':['rnh201', 'PolWT'],
            'pol1':['Pol1L868M','Pol1Y869A'],
            'pol2': ['Pol2M644G'],
            'pol3':['Pol3L612M', 'Pol3L612G']}

    # plot
    for name, genotypes in groups.items():
        subset = df[df.Genotype.isin(genotypes) & ~df.Library.isin(args.e)]
        if len(subset) == 0:
            continue
        plot_name = f'{args.o}_{name}.png'
        draw(name, plot_name, subset)

    print('Done!')


if __name__ == '__main__':
    main()


#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import linear_model
from matplotlib.ticker import FuncFormatter

# replace T to U
def replace_columns(columns, FREQ_COL_NUM):
    # mono
    if len(columns) == FREQ_COL_NUM + 12:
        for i in range(FREQ_COL_NUM, len(columns)):
            columns[i] = columns[i].replace('T', 'U')
    # dinuc
    elif len(columns) == FREQ_COL_NUM + 48:
        # NR
        if columns[FREQ_COL_NUM + 12] == 'AT':
            for i in range(FREQ_COL_NUM, len(columns)):
                if (i-FREQ_COL_NUM) % 16 >= 12:
                    columns[i]= columns[i][0] + 'U' + columns[i][2:]
        # RN
        elif columns[FREQ_COL_NUM + 12] == 'TA':
            for i in range(FREQ_COL_NUM, len(columns)):
                if (i-FREQ_COL_NUM) % 16 >= 12:
                    columns[i]= 'U' + columns[i][1:]
    return columns


def main():

    # argparse
    parser = argparse.ArgumentParser(description='Draw box plot for ars comparison')
    parser.add_argument('ars', type=argparse.FileType('r'), help='Ars region frequency file')
    parser.add_argument('-m', type=int, default=100, help='Minimum number of ribose as threshold, default=100')
    parser.add_argument('-b', action='store_true', help='Draw barplot and scatterplot')
    parser.add_argument('--notime', action='store_true', help='No time information in input file')
    parser.add_argument('--noflank', action='store_true', help='No flank information in input file')
    parser.add_argument('--ylim', type=float, default=0, help='Set y axis range of raw data, default = maximum')
    parser.add_argument('-o', default='', help='Output file basename')
    args = parser.parse_args()

    if args.o == '':
        args.o = args.ars.name.split('.')[0]

    # parameters
    FREQ_COL_NUM = 9
    if args.notime:
        FREQ_COL_NUM -= 1
    if args.noflank:
        FREQ_COL_NUM -= 1

    # get information for bed file
    data = []
    columns = args.ars.readline().rstrip('\n').split('\t')
    columns = replace_columns(columns, FREQ_COL_NUM)

    for l in args.ars:
        ws = l.rstrip('\n').split()
        for i in range(FREQ_COL_NUM,len(ws)):
            data.append(ws[:5] + [int(ws[5]), ws[6], int(ws[7]),float(ws[8]), columns[i],float(ws[i])])
    df = pd.DataFrame(data, columns = columns[:FREQ_COL_NUM]+['Composition', 'Value'])

    # set Categorical data
    lib_params = {'Genotype':['WT', 'pip', 'rnh1', 'rnh201', 'RED', 'PolWT','Pol2M644G','Pol3L612M', 'Pol3L612G','Pol1L868M','Pol1Y869A'],\
            'RE':['RE1','RE2','RE3'],\
            'Time': ['early', 'medium', 'late'],\
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C'] }
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
    try:
        la,le = sorted(df.Strand.unique())
    except:
        print('Not only two strands in strand column of input. The values are:')
        print(df.Strand.unique.sort_values())
        sys.exit(1)

    # set genotype needed for plot
    genotype_needed =  list(df.Genotype.unique())

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
    clist = sns.hls_palette(50, l=0.5, s=1)
    clist3 = [clist[17], clist[33], clist[42]]
    cpalette = [clist[48], clist[26]]

    if args.b:
        for geno in genotype_needed:
            # scatter plot for flank
            for t in times:
                da = df[(df.Genotype == geno) & (df.Time == t) & (df.Strand == le) & (df.Composition == 'A')]
                db = df[(df.Genotype == geno) & (df.Time == t) & (df.Strand == la) & (df.Composition == 'A')]
                dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype','RE','Time','Flank'])
                if len(dc) == 0:
                    continue

                color = {}
                for i in range(len(flanks)):
                    color[flanks[i]] = clist3[i]

                ds = {}
                regrs = {}
                for l in sorted(flanks):
                    ds[l]=dc[dc.Flank==l]
                    # linear regression
                    regrs[l]=linear_model.LinearRegression()
                    regrs[l].fit(ds[l].Sum_lagging.values.reshape(-1,1), ds[l].Sum_leading.values)

                fig, ax = plt.subplots(figsize=(8,7))
                plt.subplots_adjust(left=0.2, right=0.95, bottom=0.2)
                sns.scatterplot(x='Sum_lagging', y='Sum_leading', hue='Flank', data=dc, palette=color, s=100)
                # set same lim
                lim_max = max(ax.get_xlim()[1], ax.get_ylim()[1])
                ax.set_ylim([0,lim_max])
                ax.set_xlim([0,lim_max])
                # draw regr lines
                for l in flanks:
                    plt.plot([0, lim_max], regrs[l].predict([[0], [lim_max]]),color=color[l], linewidth=1.5)
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
                plt.suptitle('Ribose incorporation in leading/lagging strand ({}, {} ARS, N={:d})\n'.format(geno, t, len(dc.Library.unique())) + ', '.join([' slope_{:d}={:.4f}'.format(l, regrs[l].coef_[0]) for l in flanks]))
                plt.savefig(args.o + '_flank_{}_{}_scatter.png'.format(geno, t))
                plt.close('all')
            print('Scatter plots for flank generated!')


            # scatter plot for time
            for f in flanks:
                da = df[(df.Genotype == geno) & (df.Flank == f) & (df.Strand == le) & (df.Composition == 'A')]
                db = df[(df.Genotype == geno) & (df.Flank == f) & (df.Strand == la) & (df.Composition == 'A')]
                dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype','RE','Time','Flank'])
                if len(dc) == 0:
                    continue

                color = {}
                for i in range(len(times)):
                    color[times[i]] = clist3[i]

                ds = {}
                regrs = {}
                for t in times:
                    ds[t]=dc[dc.Time==t]
                    # linear regression
                    regrs[t]=linear_model.LinearRegression()
                    regrs[t].fit(ds[t].Sum_lagging.values.reshape(-1,1), ds[t].Sum_leading.values)

                fig, ax = plt.subplots(figsize=(8,7))
                plt.subplots_adjust(left=0.2, right=0.95, bottom=0.2)
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
                plt.suptitle('Ribose incorporation in leading/lagging strand ({}, flank = {:d}bp, N={:d})\n'.format(geno, f, len(dc.Library.unique())) + ', '.join([' slope_{}={:.4f}'.format(t, regrs[t].coef_[0]) for t in times]))
                plt.savefig(args.o + '_time_{}_{}_scatter.png'.format(geno, f))
                plt.close('all')
            print('Scatter plots for time generated!')


        # barplot count
        sns.set(style='ticks')
        for t in times:
            for f in flanks:
                # create data
                da = df[(df.Time == t) & (df.Flank == f)]
                da = da.sort_values(by=['Genotype', 'String', 'RE', 'Library'])
                if len(da) == 0:
                    continue
                fig, ax = plt.subplots(figsize=(len(da)*0.07+2,7))
                plt.subplots_adjust(bottom=0.25)
                sns.barplot(x='Library', y='Sum', hue='Strand', hue_order=[le,la], data=da, palette=cpalette, ci=None)
                if args.ylim:
                    ax.set_ylim(0, args.ylim)
                # remove box
                sns.despine()
                # fontsize size
                ax.set_xlabel(ax.get_xlabel(), fontsize=24)
                ax.set_ylabel('Count', fontsize=24)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(24)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(24)
                plt.suptitle('Ribose incorporation in leading/lagging strand\n flank = {:d}bp, {} ARS'.format(f, t))
                plt.savefig(args.o + '_strand_{:d}_{}_bar.png'.format(f,t))
                plt.close('all')
        print('Count bar plots generated!')


        # barplot percentage
        for t in times:
            for f in flanks:
                # create data
                da = df[(df.Time == t) & (df.Strand == le) & (df.Flank == f) & (df.Composition == 'A')]
                db = df[(df.Time == t) & (df.Strand == la) & (df.Flank == f) & (df.Composition == 'A')]
                dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype','RE','Time','Flank'])
                if len(dc) == 0:
                    continue
                dc = dc.sort_values(by=['Genotype', 'String', 'RE', 'Library'])

                dc['total'] = dc.Sum_lagging + dc.Sum_leading
                dc['Lagging'] = dc.Sum_lagging/dc.total
                dc['Leading'] = dc.Sum_leading/dc.total
                dc['Ratio'] = dc.Sum_leading/dc.Sum_lagging
                fig, ax = plt.subplots(figsize=(len(dc)*0.7+9,10))
                plt.subplots_adjust(left=0.05, bottom=0.45, right=0.9, top=0.92)
                sns.barplot(x=dc.Library, y=dc.Lagging + dc.Leading, color=clist[26], label='Lagging', ax=ax)
                sns.barplot(x=dc.Library, y=dc.Leading, color=clist[48],label='Leading', ax=ax)
                plt.suptitle('Ribose incorporation in leading/lagging strand\n flank = {:d}bp, {} ARS'.format(f, t))
                # legend
                plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1,frameon=False, prop={'size': 24})
                ax.set_ylim([0,1])
                # xtick label
                lib_full = dc.apply(lambda x: '-'.join([str(x.String), str(x.Genotype), str(x.RE), str(x.Library)]), axis=1)
                ax.set_xticklabels(lib_full, rotation=90)
                # add a line
                plt.axhline(y=0.5, linewidth=2, color='black', linestyle='--')
                # percentage y-axis
                ax.yaxis.set_major_formatter(f_percentage)
                # remove box
                sns.despine()
                # fontsize size
                ax.set_xlabel(ax.get_xlabel(), fontsize=24)
                ax.set_ylabel('Percentage', fontsize=24)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(24)
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(20)
                # add labels
                for nle, nla, rect in zip(dc.Leading, dc.Lagging,ax.patches):
                    ax.text(rect.get_x() + rect.get_width()/2, 0.05, '{:.0%}'.format(nle), ha='center', color='black', fontsize=18, fontweight='bold')
                    ax.text(rect.get_x() + rect.get_width()/2, 0.93, '{:.0%}'.format(nla), ha='center', color='black', fontsize=18, fontweight='bold')
                plt.savefig(args.o + '_strand_percent_{:d}_{}_bar.png'.format(f,t))
                plt.close('all')
        print('Percentage bar plots generated!')
        sys.exit()

    # boxplots for various comparison
    sns.set(style='ticks')

    # compare leading/lagging
    for t in times:
        for f in flanks:
            for s in genotype_needed:
                if s == 'rnh201':
                    draw = sns.boxplot
                else:
                    draw = sns.swarmplot
                color = 'Set1'

                # compare leading/lagging
                da = df[(df.Genotype == s) & (df.Time == t) & (df.Flank == f) & (df.Sum >= args.m) & ~ df.Composition.str.contains('%') & ~ df.Composition.str.contains('n')]
                db = df[(df.Genotype == s) & (df.Time == t) & (df.Flank == f) & (df.Sum >= args.m) & (df.Composition.str.contains('%'))]
                dc = df[(df.Genotype == s) & (df.Time == t) & (df.Flank == f) & (df.Sum >= args.m) & (df.Composition.str.contains('n'))]
                # skip if none
                if len(da) == 0:
                    continue

                # draw plot compare leading/lagging wt
                for norm, data in zip(['raw','denorm', 'norm'], [da, db,dc]):
                    fig, ax =plt.subplots(figsize=(len(columns)*0.2+4,7))
                    plt.subplots_adjust(left=0.25, bottom=0.15, top=0.9)
                    draw(x='Composition', y='Value', hue='Strand', hue_order=[le,la],data=data, palette=cpalette)
                    if args.ylim and norm == 'raw':
                        ax. set_ylim(0, args.ylim)
                    plt.suptitle('Comparison for leading/lagging \n{}, {}, flank = {:d}bp, {} ARS, threshold={:d}'.format(norm, s, f, t, args.m))
                    # remove top and right line
                    sns.despine()
                    # set x tick labels
                    ax.set_xticklabels(list(map(lambda x: x.get_text().replace('n','').replace('%',''), ax.get_xticklabels())), fontsize=24)
                    # labels:
                    ax.set_xlabel('Ribonucleotide', fontsize=24)
                    if norm == 'raw':
                        ax.set_ylabel('Count', fontsize=24)
                    elif norm == 'denorm':
                        ax.set_ylabel('Percentage', fontsize=24)
                        ax.yaxis.set_major_formatter(f_percentage)
                    elif norm == 'norm':
                        ax.set_ylabel('Z-score', fontsize=24)
                    # tick size
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    # set legend
                    ax.legend(frameon=False, prop={'size':20})
                    plt.savefig(args.o + '_strand_{}_{}_{}_{}.png'.format(norm, s,t, f))
                    plt.close('all')


    # compare flank
    for t in times:
        for s in [le,la]:
            for c in genotype_needed:
                if s == le:
                    color = 'Set1'
                else:
                    color = 'Set2'

                if c == 'rnh201':
                    draw = sns.boxplot
                else:
                    draw = sns.swarmplot
                    color = 'Set1'

                da = df[(df.Genotype == c) & (df.Time == t) & (df.Strand == s) & (df.Sum >= args.m) & ~ df.Composition.str.contains('%') & ~ df.Composition.str.contains('n')]
                db = df[(df.Genotype == c) & (df.Time == t) & (df.Strand == s) & (df.Sum >= args.m) & (df.Composition.str.contains('%'))]
                dc = df[(df.Genotype == c) & (df.Time == t) & (df.Strand == s) & (df.Sum >= args.m) & (df.Composition.str.contains('n'))]

                # skip if none
                if len(da) == 0:
                    continue

                # draw plot
                for norm, data in zip(['raw','denorm', 'norm'], [da, db,dc]):
                    fig, ax =plt.subplots(figsize=(len(columns)*0.2+4,7))
                    plt.subplots_adjust(left=0.25, bottom=0.15, top=0.9)
                    draw(x='Composition', y='Value', hue='Flank', data=data, palette=color)
                    if args.ylim and norm == 'raw':
                        ax. set_ylim(0, args.ylim)
                    # remove top and right line
                    sns.despine()
                    plt.suptitle('Comparison for different flank length \n{}, {}, {} strand, {} ARS, threshold={:d}'.format(norm, c, s, t, args.m))
                    # set x tick labels
                    ax.set_xticklabels(list(map(lambda x: x.get_text().replace('n','').replace('%',''), ax.get_xticklabels())), fontsize=24)
                    # labels:
                    ax.set_xlabel('Ribonucleotide', fontsize=24)
                    if norm == 'raw':
                        ax.set_ylabel('Count', fontsize=24)
                    elif norm == 'denorm':
                        ax.set_ylabel('Percentage', fontsize=24)
                        ax.yaxis.set_major_formatter(f_percentage)
                    elif norm == 'norm':
                        ax.set_ylabel('Z-score', fontsize=24)
                    # tick size
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    # set legend
                    ax.legend(frameon=False, prop={'size':20})
                    plt.savefig(args.o + '_flank_{}_{}_{}_{}.png'.format(c, s, t, norm))
                    plt.close('all')


    # box plot compare time
    for f in flanks:
        for s in [le,la]:
            for c in genotype_needed:
                if s == le:
                    color = 'Set1'
                else:
                    color = 'Set2'

                if c == 'rnh201':
                    draw = sns.boxplot
                else:
                    draw = sns.swarmplot
                    color = 'Set1'

                da = df[(df.Genotype == c) & (df.Flank == f) & (df.Strand == s) & (df.Sum >= args.m) & ~ df.Composition.str.contains('%') & ~ df.Composition.str.contains('n')]
                db = df[(df.Genotype == c) & (df.Flank == f) & (df.Strand == s) & (df.Sum >= args.m) & (df.Composition.str.contains('%'))]
                dc = df[(df.Genotype == c) & (df.Flank == f) & (df.Strand == s) & (df.Sum >= args.m) & (df.Composition.str.contains('n'))]

                # skip if none
                if len(da) == 0:
                    continue

                # draw plot
                for norm, data in zip(['raw','denorm', 'norm'], [da, db,dc]):
                    fig, ax =plt.subplots(figsize=(len(columns)*0.2+4,7))
                    plt.subplots_adjust(left=0.25, bottom=0.15, top=0.9)
                    draw(x='Composition', y='Value', hue='Time',data=data, palette=color)
                    if args.ylim and norm == 'raw':
                        ax. set_ylim(0, args.ylim)
                    plt.suptitle('Comparison for different firing time \n{}, {}, {} strand, flank = {:d}bp, threshold={:d}'.format(norm, c, s, f, args.m))
                    # remove top and right line
                    sns.despine()
                    # set x tick labels
                    ax.set_xticklabels(list(map(lambda x: x.get_text().replace('n','').replace('%',''), ax.get_xticklabels())), fontsize=24)
                    # labels:
                    ax.set_xlabel('Ribonucleotide', fontsize=24)
                    if norm == 'raw':
                        ax.set_ylabel('Count', fontsize=24)
                    elif norm == 'denorm':
                        ax.set_ylabel('Percentage', fontsize=24)
                        ax.yaxis.set_major_formatter(f_percentage)
                    elif norm == 'norm':
                        ax.set_ylabel('Z-score', fontsize=24)
                    # tick size
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    # set legend
                    ax.legend(frameon=False, prop={'size':20})
                    plt.savefig(args.o + '_time_{}_{}_{:d}_{}.png'.format(c, s, f, norm))
                    plt.close('all')

    print('Done!')

if __name__ == '__main__':
    main()




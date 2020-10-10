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
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C','RS'] }
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
            fig, ax = plt.subplots(figsize=(len(dc)*0.7+9,6))
            sns.barplot(x=dc.Library, y=dc.Lagging + dc.Leading, color=clist[26], label='Lagging', ax=ax)
            sns.barplot(x=dc.Library, y=dc.Leading, color=clist[48],label='Leading', ax=ax)
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
            # setting for publication
            plt.subplots_adjust(left=0.036, bottom=0.28, right=0.935, top=0.97)
            # tick label
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(36)
            label_texts = [tick.label.get_text().split('-')[-1] for tick in ax.xaxis.get_major_ticks()]
            ax.set_xticklabels(label_texts, fontsize=36)
            ax.set_xlabel('')
            ax.set_ylabel('')
            # add labels
            for nle, nla, rect in zip(dc.Leading, dc.Lagging,ax.patches):
                inte, deci = f'{nle*100:.1f}'.split('.')
                nle_label = f'{inte}.\n{deci}%'
                inte, deci = f'{nla*100:.1f}'.split('.')
                nla_label = f'{inte}.\n{deci}%'
                ax.text(rect.get_x() + rect.get_width()/2, 0.05, nle_label, ha='center', fontsize=28, fontweight='bold')
                ax.text(rect.get_x() + rect.get_width()/2, 0.8, nla_label, ha='center', fontsize=28, fontweight='bold')
            plt.savefig(args.o + '_strand_percent_{:d}_{}_bar.png'.format(f,t))
            plt.close('all')
    print('Percentage bar plots generated!')

    print('Done!')

if __name__ == '__main__':
    main()




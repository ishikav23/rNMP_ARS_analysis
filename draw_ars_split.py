#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
import sys

# turn off warning
pd.options.mode.chained_assignment = None

# set legend to ACGU
def set_legend(ax):
    # ax: axes
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) != 5:
        print('Number of legends is not 4! Cannot set to A,C,G,U')
    else:
        ax.legend(handles=handles, labels=['Type','A','C','G','U'])


# sort features
def categorize_df(df, orders):
    for w in orders.keys():
        if w not in df.columns:
            continue
        categorical_order = []
        unique_values = df[w].unique()
        for v in orders[w]:
            if v in unique_values:
                categorical_order.append(v)
        for v in unique_values:
            if v not in categorical_order:
                categorical_order.append(v)
        df[w] = pd.Categorical(df[w], categorical_order)


# hide label and titles
def hide_label(fig, ax):
    # hide title
    fig.suptitle("")
    # hide axis labels
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)
    # hide legend
    ax.get_legend().remove()
    # increase plot size
    fig.subplots_adjust(top=0.99, left=0.2, right=0.99, bottom=0.08)
    # tick size
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(36)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(36)


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw box plot for ars comparison')
    parser.add_argument('ars', type=argparse.FileType('r'), help='Ars region frequency file')
    parser.add_argument('-m', type=int, default=100, help='Minimum number of ribose as threshold, default=100')
    parser.add_argument('-o', default='', help='Output basename')
    parser.add_argument('--bar', action='store_true', help='Draw bar plot of each library')
    parser.add_argument('--box', action='store_true', help='Draw box plot for ratios')
    parser.add_argument('--line', action='store_true', help='Draw  plot line chart for ratio trend')
    parser.add_argument('--no-label', action='store_true', help='Hide label, header and legend in the plot')
    args = parser.parse_args()
    if not any([args.bar, args.box, args.line]):
        args.line = True
    if args.o == '':
        args.o = args.ars.name.split('.')[0]


    # get information for bed file
    df = pd.read_csv(args.ars, sep='\t')
    # convert to kbp
    df.Position = df.Position / 1000
    # set Categorical data
    lib_params = {'Genotype':['WT', 'pip', 'rnh1', 'rnh201', 'RED', 'PolWT','Pol2M644G',\
         'Pol3L612M', 'Pol3L612G','Pol1L868M','Pol1Y869A'],\
            'RE':['RE1','RE2','RE3'],\
            'Time': ['early', 'medium', 'late'],\
            'Strand': ['leading', 'lagging'],\
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C'] }
    categorize_df(df, lib_params)


    # set genotype needed for plot
    genotype_needed =  df.Genotype.unique()
    # column number for first rnmp column
    RNMP_COL_NUM = 9

    # formatter for scientific
    def my_scientific_format(y):
        a = f'{y:.1E}'
        m, e = a.split('E')
        e = e.replace('+0','')
        e = e.replace('-0','-')
        return f'{m}E{e}'
    f_scientific = FuncFormatter(lambda y, _: my_scientific_format(y))

    # get available categories
    libs = df.Library.unique()
    times = df.Time.unique().sort_values()
    treatments = df.Genotype.unique()
    p = sorted(df.Position.unique())
    len_part = p[1] - p[0]
    flank_start = df.Position.min()
    len_flank = df.Position.max()
    print('Finish reading data!')
    print('Available replication peak times: {}.\n'.format(', '.join(times)) + \
            'Available treatments: {}.\n'.format(', '.join(treatments)) + \
            'Flank:{} - {}kbp, Part length: {}kbp\n'.format(flank_start, len_flank, len_part) + \
            'Libraries: {}.'.format(', '.join([str(x) for x in libs])))

    # set color parameters
    sns.set(style='whitegrid')
    colorlist = sns.hls_palette(8, l=.5, s=1)
    # color for bases
    colorb = [colorlist[i] for i in [6,4,1,2]]
    # barplot count
    co = {'leading':colorlist[0], 'lagging':colorlist[1]}
    # build color dict for time
    # remove yellow color which is not print friendly
    colorlist = colorlist[:1] + colorlist[3:]
    colort = {}
    for i in range(len(times)):
        colort[times[i]] = colorlist[i]

    # draw barplot for libs
    # array to store ratios:
    ratios = []
    # array to store percentage of leading and lagging
    # lib_percentage = []
    lib_needed = {}
    for t in times:
        lib_needed[t] = []
        for l in libs:
            da = df[(df.Time == t) & (df.Library == l) & (df.Strand == 'leading')].sort_values(by='Position')
            db = df[(df.Time == t) & (df.Library == l) & (df.Strand == 'lagging')].sort_values(by='Position')
            dc = da.merge(db, suffixes=['_leading','_lagging'],on=['Library','String','Genotype', 'Time', 'RE','Position'])
            dc['Total'] = dc.Sum_lagging + dc.Sum_leading
            # filter for threshold
            total_ribo = dc.Total.sum()
            if total_ribo >= args.m:
                lib_needed[t].append(l)
            ratios.append(dc)

            if args.bar:
                fig, ax = plt.subplots(figsize=(len(da)/2+2,6))
                sns.barplot(x='Position', y='Sum', hue='Strand', hue_order=['leading','lagging'], data=da, palette=co, ci=None)
                plt.suptitle('Ribos incorporation in leading/lagging Strand\n {}, Time < {} min'.format(l, t))
                plt.savefig(args.o + '_{}_{}_bar.png'.format(l,t))
                plt.close('all')

                dc['Lagging'] = dc.Sum_lagging/dc.Total
                dc['Leading'] = dc.Sum_leading/dc.Total
                dc['Ratio'] = dc.Sum_leading/dc.Sum_lagging
                fig, ax = plt.subplots(figsize=(len(dc)/2+2,6))
                sns.barplot(x=dc.Position, y=dc.Leading, color=colorlist[0], ci=None, label='Leading')
                sns.barplot(x=dc.Position, y=dc.Lagging, bottom=dc.Leading, color=colorlist[1], ci=None, label='Lagging')
                plt.suptitle('Ribos incorporation percentage in leading/lagging Strand\n {}, Time < {} min'.format(l, t))
                ax.set_ylabel('Leading Percentage')
                plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
                plt.savefig(args.o + '_{}_{}_bar_percent.png'.format(l,t))
                plt.close('all')
    if args.bar:
        print('bar plots generated!')
    # generate data frame
    dmerge = pd.concat(ratios, sort=True)


    # box plots for ratios
    # box plot for ratios
    if args.box:
        for tr in genotype_needed:
            db = dmerge[dmerge.Genotype == tr]
            if db.empty:
                continue
            fig, ax = plt.subplots(figsize=(len(db.Position.unique())*0.8+2,8))
            sns.boxplot(x='Position', y='Ratio', hue='Time', data=db, palette=colort)
            plt.suptitle('Leading/lagging ratio with increasing distance to ARS ({})'.format(tr))
            # set maximum y lim to 4
            ax.set_ylim([0, min(4, ax.get_ylim()[1])])
            ax.set_ylabel('Leading/Lagging')
            ax.set_xlabel('Distance')
            plt.savefig(args.o + '_ratio_{}_box.png'.format(tr))
            plt.close('all')
        print('Box plot for ratio generated!')


    if args.line:
        # feature
        features= ['Sum']
        for i in df.columns[RNMP_COL_NUM:]:
            if i.find('%') == -1:
                features.append(i)
        # calculate the estimators for ratio
        flanks = df.Position.unique()
        # store average
        dsummary = []
        dleading = []
        dlagging = []
        for tr in genotype_needed:
            for t in times:
                for f in flanks:
                    db = dmerge[(dmerge.Genotype == tr) & (dmerge.Time == t) & (dmerge.Position == f) & (dmerge.Library.isin(lib_needed[t]))]
                    if db.empty:
                        continue
                    db['Ratio'] = db['Sum_leading'] / db['Sum_lagging']
                    db['Length'] = db['Sum_leading'] /db['RPB_leading']
                    dsummary.append([tr,t,f, db.Total.sum(), db.Ratio.median()] + [db['{}_leading'.format(fe)].sum()/db['{}_lagging'.format(fe)].sum() for fe in features])
                    dleading.append([tr,t,f, db.Sum_leading.sum(), db.RPB_leading.mean()] + [db['{}_leading'.format(fe)].mean() for fe in features])
                    dlagging.append([tr,t,f, db.Sum_lagging.sum(), db.RPB_lagging.mean()] + [db['{}_lagging'.format(fe)].mean() for fe in features])
        df_summary = pd.DataFrame(dsummary, columns = ['Genotype', 'Time', 'Position', 'Total', 'Median'] + ['{}'.format(fe) for fe in features])
        df_leading = pd.DataFrame(dleading, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB'] + ['{}'.format(fe) for fe in features])
        df_lagging = pd.DataFrame(dlagging, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB'] + ['{}'.format(fe) for fe in features])
        categorize_df(df_summary, lib_params)
        categorize_df(df_leading, lib_params)
        categorize_df(df_lagging, lib_params)

        # draw line chart for small window summaries
        feature_groups = {'Raw':[],'Normalized':[]}
        columns = list(df_summary.columns)
        for i in columns[6:]:
            if i.find('n') == -1:
                feature_groups['Raw'].append(i)
            else:
                feature_groups['Normalized'].append(i)
        # start plotting
        for tr in genotype_needed:
            db = df_summary[(df_summary.Genotype == tr)].sort_values(by='Time')
            if db.empty:
                continue
            # Median, MLE_Sum
            for r in ['Median','Sum']:
                fig, ax = plt.subplots(figsize=(9,9))
                sns.lineplot(x='Position', y=r, hue='Time', data=db, palette=colort, linewidth=4)
                plt.suptitle('Leading/lagging ratio with increasing distance to ARS\n ({}, {})'.format(tr, r))
                ax.set_ylabel('Leading/Lagging ratio')
                ax.set_xlabel('Distance to ARS (kb)')
                if args.no_label:
                    hide_label(fig, ax)
                plt.savefig(args.o + '_ratio_{}_{}.png'.format(tr,r).lower())
                plt.close('all')
            # total
            fig, ax = plt.subplots(figsize=(9,9))
            sns.lineplot(x='Position', y='Total', hue='Time', data=db, palette=colort, linewidth=4)
            plt.suptitle('Ribonucleotides incorporation with increasing distance to ARS\n ({})'.format(tr))
            ax.set_ylabel('Counts')
            ax.set_xlabel('Distance to ARS (kb)')
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(24)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(24)
            if args.no_label:
                hide_label(fig, ax)
            plt.savefig(args.o + '_total_{}.png'.format(tr).lower())
            plt.close('all')
            # Compositions for leading/lagging ratio
            # build dataframe
            for k,v in feature_groups.items():
                dcomp = []
                for index, l in db.iterrows():
                    for fe in v:
                        dcomp.append([l.Time, l.Position, fe, l[fe]])
                df_comp = pd.DataFrame(dcomp, columns=['Time','Position','Feature','Value'])
                for w in ['Time']:
                    for v in df_summary[w].unique():
                        if v not in lib_params[w]:
                            lib_params[w].append(v)
                    df_comp[w] = pd.Categorical(df_comp[w], lib_params[w])
                # draw split
                for t in times:
                    da = df_comp[df_comp.Time == t]
                    fig, ax = plt.subplots(figsize=(9,9))
                    sns.lineplot(x='Position', y='Value', hue='Feature', data=da, palette=colorb, linewidth=4)
                    plt.suptitle('Leading/lagging ratio for rNMP with increasing distance to ARS \n({}, {}, {})'.format(tr, t, k))
                    ax.set_ylabel('Leading/Lagging ratio')
                    ax.set_xlabel('Distance to ARS (kb)')
                    set_legend(ax)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    if args.no_label:
                        hide_label(fig, ax)
                    plt.savefig(args.o + '_rNMP' +'_{}_{}_{}.png'.format(tr,t,k).lower())
                    plt.close('all')

        # leading and lagging
        for df_lela, slela in [[df_leading, 'leading'], [df_lagging, 'lagging']]:
            for tr in genotype_needed:
                dlela = df_lela[(df_lela.Genotype == tr)].sort_values(by='Time')
                if dlela.empty:
                    continue
                # total
                for ydata in ['Total','RPB']:
                    fig, ax = plt.subplots(figsize=(9,9))
                    sns.lineplot(x='Position', y=ydata, hue='Time', data=dlela, palette=colort, linewidth=4)
                    plt.suptitle(f'Ribonucleotides incorporation with increasing distance to ARS \n({slela} strand, {tr})')
                    if ydata == 'Total':
                        ax.set_ylabel('Counts')
                        ax.yaxis.set_major_formatter(f_scientific)
                    else:
                        if tr == 'WT':
                            ax.yaxis.set_major_formatter(f_scientific)
                        ax.set_ylabel('RPB')
                    ax.set_xlabel('Distance to ARS (kb)')
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    if args.no_label:
                        hide_label(fig, ax)
                    plt.savefig(args.o + f'_{ydata}_{slela}_{tr}.png'.lower())
                    plt.close('all')
                # Compositions
                # build dataframe
                for k,v in feature_groups.items():
                    dcomp = []
                    for index, l in dlela.iterrows():
                        for fe in v:
                            dcomp.append([l.Time, l.Position, fe, l[fe]])
                    df_comp = pd.DataFrame(dcomp, columns=['Time','Position','Feature','Value'])
                    # draw split
                    for t in times:
                        da = df_comp[df_comp.Time == t]
                        fig, ax = plt.subplots(figsize=(9,9))
                        sns.lineplot(x='Position', y='Value', hue='Feature', data=da, palette=colorb, linewidth=4)
                        plt.suptitle(f'Frequency of rNMP with increasing distance to ARS\n ({slela} strand, {tr}, {t}, {k})')
                        ax.yaxis.set_major_formatter(f_scientific)
                        ax.set_ylabel('Frequency')
                        ax.set_xlabel('Distance to ARS (kb)')
                        set_legend(ax)
                        for tick in ax.yaxis.get_major_ticks():
                            tick.label.set_fontsize(24)
                        for tick in ax.xaxis.get_major_ticks():
                            tick.label.set_fontsize(24)
                        if args.no_label:
                            hide_label(fig, ax)
                        plt.savefig(args.o + f'_rNMP_{slela}_{tr}_{t}_{k}.png'.lower())
                        plt.close('all')
                # stacked area plot
                for t in times:
                    da = dlela[dlela.Time == t].sort_values(by='Position')
                    da['total'] = 0
                    for fe in feature_groups['Normalized']:
                        da['total'] += da[fe]
                    data = []
                    for fe in feature_groups['Normalized']:
                        data.append((da[fe]/da['total']).tolist())
                    fig, ax = plt.subplots(figsize=(9,9))
                    plt.stackplot(da.Position.tolist(), data, labels=['A','C','G','U'], colors=colorb)
                    plt.suptitle(f'Relative rNMP incorporation rate during DNA replication process\n ({slela} strand, {tr}, {t}, {k})')
                    plt.legend()
                    ax.set_ylabel('Rate')
                    ax.set_xlabel('Distance to ARS (kb)')
                    plt.ylim([0,1])
                    plt.xlim([da.Position.min(), da.Position.max()])
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    if args.no_label:
                        hide_label(fig, ax)
                    plt.savefig(args.o + f'_rNMP_rate_{slela}_{tr}_{t}_{k}.png'.lower())
                    plt.close('all')


        print('line chart generated!')
        print('Done!')

if __name__ == '__main__':
    main()




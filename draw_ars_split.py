#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter, ScalarFormatter
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
    return df


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
    parser.add_argument('--ppb', action='store_true', help='PPB data is available')
    parser.add_argument('--bar', action='store_true', help='Draw bar plot of each library')
    parser.add_argument('--box', action='store_true', help='Draw box plot for ratios')
    parser.add_argument('--line', action='store_true', help='Draw  plot line chart for ratio trend')
    parser.add_argument('--no-label', action='store_true', help='Hide label, header and legend in the plot')
    parser.add_argument('--same_scale', action='store_true', help='Use same scale for all plots')
    parser.add_argument('--hy', action='store_true', help='HydEn-seq libraries')
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
            'Time': ['early', 'medium', 'late', 'high','low'],\
            'Strand': ['leading', 'lagging'],\
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C'] }
    categorize_df(df, lib_params)


    # set genotype needed for plot
    genotype_needed =  ['WT', 'rnh201','PolWT','Pol2M644G','Pol3L612M', 'Pol3L612G','Pol1L868M','Pol1Y869A']
    # column number for first rnmp column
    RNMP_COL_NUM = 9 if not args.ppb else 10

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
    colorlist = ['#a570f3','#1bce77','#d95f02']
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
    # generate data frame
    dmerge = pd.concat(ratios, sort=True)

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
                    dsummary.append([tr,t,f, db.Total.sum(), db.Ratio.median(), db.Ratio.std()] + [db['{}_leading'.format(fe)].sum()/db['{}_lagging'.format(fe)].sum() for fe in features])
                    if args.ppb:
                        dleading.append([tr,t,f, db.Sum_leading.sum(), db.RPB_leading.mean(), db.RPB_leading.std(), db.PPB_leading.mean(), db.PPB_leading.std()] + [db['{}_leading'.format(fe)].mean() for fe in features])
                        dlagging.append([tr,t,f, db.Sum_lagging.sum(), db.RPB_lagging.mean(), db.RPB_lagging.std(), db.PPB_lagging.mean(), db.PPB_lagging.std()] + [db['{}_lagging'.format(fe)].mean() for fe in features])
                    else:
                        dleading.append([tr,t,f, db.Sum_leading.sum(), db.RPB_leading.mean(), db.RPB_leading.std()] + [db['{}_leading'.format(fe)].mean() for fe in features])
                        dlagging.append([tr,t,f, db.Sum_lagging.sum(), db.RPB_lagging.mean(), db.RPB_lagging.std()] + [db['{}_lagging'.format(fe)].mean() for fe in features])
        df_summary = pd.DataFrame(dsummary, columns = ['Genotype', 'Time', 'Position', 'Total', 'Median', 'Std'] + ['{}'.format(fe) for fe in features])
        if args.ppb:
            df_leading = pd.DataFrame(dleading, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB', 'RPB_std', 'PPB', 'PPB_std'] + ['{}'.format(fe) for fe in features])
            df_lagging = pd.DataFrame(dlagging, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB', 'RPB_std', 'PPB', 'PPB_std'] + ['{}'.format(fe) for fe in features])
        else:
            df_leading = pd.DataFrame(dleading, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB', 'RPB_std'] + ['{}'.format(fe) for fe in features])
            df_lagging = pd.DataFrame(dlagging, columns = ['Genotype', 'Time', 'Position', 'Total', 'RPB', 'RPB_std'] + ['{}'.format(fe) for fe in features])
        df_summary = categorize_df(df_summary, lib_params)
        df_leading = categorize_df(df_leading, lib_params)
        df_lagging = categorize_df(df_lagging, lib_params)
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
            for r in ['Sum']:
                fig, ax = plt.subplots(figsize=(9,9))
                sns.lineplot(x='Position', y=r, hue='Time', data=db, palette=colort, linewidth=4)
                # add std
                facecolors = ['#cc99cc', '#99cc99']
                c = 0
                for t in db.Time.unique():
                    dc = db[db.Time == t].sort_values('Position')
                    ax.fill_between(dc.Position, dc.Sum + dc.Std, dc.Sum - dc.Std, facecolor = facecolors[c], alpha=0.4)
                    c += 1
                if not args.hy:
                    plt.ylim((max(0, ax.get_ylim()[0]), min(ax.get_ylim()[1], 6)))
                plt.suptitle('Leading/lagging ratio with increasing distance to ARS\n ({}, {})'.format(tr, r))
                ax.set_ylabel('Leading/Lagging ratio')
                ax.set_xlabel('Distance to ARS (kb)')
                if args.no_label:
                    hide_label(fig, ax)
                fig.subplots_adjust(top=0.99, left=0.15, right=0.99, bottom=0.08)
                if args.same_scale:
                    plt.ylim([0.2,6])
                    plt.yscale('log')
                    ax.yaxis.set_major_formatter(ScalarFormatter())
                    ax.yaxis.set_ticks([0.2,0.5,1,2,5])
                plt.savefig(args.o + '_ratio_{}_{}.png'.format(tr,r).lower())
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
        # leading and lagging
        for tr in genotype_needed:
            for df_lela, slela in [[df_leading, 'leading'], [df_lagging, 'lagging']]:
                dlela = df_lela[(df_lela.Genotype == tr)].sort_values(by='Time')
                if dlela.empty:
                    continue
                # total
                for ydata in ['RPB', 'PPB']:
                    fig, ax = plt.subplots(figsize=(9,9))
                    sns.lineplot(x='Position', y=ydata, hue='Time', data=dlela, palette=colort, linewidth=4)
                    ylims = ax.get_ylim()
                    # add std
                    facecolors = ['#cc99cc', '#99cc99']
                    c = 0
                    for t in dlela.Time.unique():
                        dc = dlela[dlela.Time == t].sort_values('Position')
                        ax.fill_between(dc.Position, dc[ydata] + dc[f'{ydata}_std'], dc[ydata] - dc[f'{ydata}_std'], facecolor = facecolors[c], alpha=0.4)
                        c += 1
                    if not args.hy:
                        plt.ylim((max(0, ax.get_ylim()[0]), min(ax.get_ylim()[1], 3)))
                    plt.suptitle(f'Ribonucleotides incorporation with increasing distance to ARS \n({slela} strand, {tr})')
                    # ax.ticklabel_format(useOffset=False)
                    if ydata in ['Total','RPB']:
                        ax.set_ylabel('Counts')
                    else:
                        ax.set_ylabel(ydata)
#                         ax.set_ylim([2.2e-8,7.8e-8])
#                         ax.yaxis.set_ticks([3e-8,5e-8,7e-8])
                        ax.yaxis.set_major_formatter(f_scientific)
                    ax.set_xlabel('Distance to ARS (kb)')
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(24)
                    if args.no_label:
                        hide_label(fig, ax)
                    # plt.ylim((0.03, 0.082))
                    plt.savefig(args.o + f'_{ydata}_{slela}_{tr}.png'.lower())
                    plt.close('all')
        print('line chart generated!')
        print('Done!')

if __name__ == '__main__':
    main()




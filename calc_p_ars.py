#!/usr/bin/python3

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from sklearn import linear_model
from matplotlib.ticker import FuncFormatter

def main():

    # argparse
    parser = argparse.ArgumentParser(description='Run statistical test for ARS region comparison\n Strand: Wilcoxon, Time and flank: ')
    parser.add_argument('ars', type=argparse.FileType('r'), help='Ars region frequency file')
    parser.add_argument('-l', nargs='*', default=[5000,10000, 15000], help='Length of flank region to check, default=[5000,10000]')
    parser.add_argument('-t', nargs='*', default=['25','30', 'all'], help='Replication peak time threshold to check, default=[25.30]')
    parser.add_argument('-m', type=int, default=100, help='Minimum number of ribose as threshold, default=100')
    parser.add_argument('-s', nargs='*',default=['leading', 'lagging'], help='[leading,lagging]')
    parser.add_argument('--ttest', action='store_true', help='Use paired t-test instead of Wilcoxon')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    le,la = args.s

    # get information for bed file
    df = pd.read_csv(args.ars, sep='\t')
    df['Genotype'] = pd.Categorical(df['Genotype'], ['WT', 'pip', 'rnh1', 'rnh201', 'RED'])
    df['RE'] = pd.Categorical(df['RE'], ['RE1', 'RE2', 'RE3'])
    df['String'] = pd.Categorical(df['String'], ['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C'])
    df['Strand'] = pd.Categorical(df['Strand'], args.s)
    df['Time'] = pd.Categorical(df['Time'], args.t)
    df['Flank'] = pd.Categorical(df['Flank'], args.l)
    df = df.sort_values(by=['Genotype', 'String', 'RE', 'Time', 'Flank', 'Strand'])

    # group
    group = [['WT'], ['rnh201'], ['WT', 'pip', 'rnh1'], ['rnh201', 'RED']]

   # fix t, f, compare s
    cols = df.columns
    feature = cols[7:]
    args.o.write('Time\tFlank\tGenotype\tNum\t' + '\t'.join(feature) + '\n')
    for t in args.t:
        for f in args.l:
            for g in group:
                da = df[(df.Time == t) & (df.Strand == le) & (df.Flank == f) & (df.Genotype.isin(g))]
                db = df[(df.Time == t) & (df.Strand == la) & (df.Flank == f) & (df.Genotype.isin(g))]
                if args.ttest:
                    p = [stats.ttest_rel(da[x], db[x])[1] for x in feature]
                else:
                    p = [stats.wilcoxon(da[x], db[x])[1] for x in feature]
                args.o.write('\t'.join([t,str(f),'&'.join(g), str(len(da))] + [str(x) for x in p]) +'\n')
    print('Done!')

if __name__ == '__main__':
    main()




#!/usr/bin/env python3

import pandas as pd
import argparse
import sys


# read files
def read_files(frs, lib_info):
    data = []
    for fr in frs:
        columns = fr.readline().rstrip().split('\t')
        columns = ['Library', 'String', 'Genotype','RE'] + columns[1:]
        for l in fr:
            ws = l.rstrip().split()
            info = ws[0].split('-')
            if len(info) != 3:
                continue
            # skip unwilling libraries
            if info[0] not in lib_info:
                continue
            data.append(info + [lib_info[info[0]][2]]+ [float(x) for x in ws[1:]])
    df = pd.DataFrame(data)
    df.columns = columns

    return df


def main():
    parser = argparse.ArgumentParser(description='Sort data for figure 1')
    parser.add_argument('info', type=argparse.FileType('r'), help='Information of libraries')
    parser.add_argument('tsv', nargs='+', type=argparse.FileType('r'), help='Input files')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # read informations
    lib_info = {}
    for l in args.info:
        ws = l.rstrip('\n').split('\t')
        if len(ws) != 4:
            continue
        lib_info[ws[0]] = ws[1:]

    # read csv
    df = read_files(args.tsv, lib_info)

    # set Categorical data
    lib_params = {'Genotype':['WT', 'pip', 'rnh1', 'rnh201', 'RED', 'PolWT','Pol2M644G','Pol3L612M', 'Pol3L612G','Pol1L868M','Pol1Y869A'],\
            'RE':['RE1','RE2','RE3'],\
            'String':['E134', 'BY4741', 'BY4742', 'YFP17', 'W303', 'S288C','RS'] }
    for w in ['Genotype','RE','String']:
        categorical_order = []
        unique_values = df[w].unique()
        for v in lib_params[w]:
            if v in unique_values:
                categorical_order.append(v)
        for v in unique_values:
            if v not in categorical_order:
                categorical_order.append(v)
        df[w] = pd.Categorical(df[w], categorical_order)

    # sort
    df = df.sort_values(['Genotype', 'String', 'RE', 'Library'])

    # merge columns
    df['chrom'] = df['Library'].astype(str) + '-' + df['String'].astype(str) +\
        '-' + df['Genotype'].astype(str) + '-' + df['RE'].astype('str')
    df = df[['chrom'] + df.columns[4:-1].tolist()]

    # out
    df.to_csv(args.o, sep='\t', index=False)

if __name__ == '__main__':
    main()


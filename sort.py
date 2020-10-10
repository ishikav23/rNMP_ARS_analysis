#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Sort data for figure 1')
    parser.add_argument('tsv', type=argparse.FileType('r'), help='Input file')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # read csv
    df = pd.read_csv(args.tsv, sep='\t')

    # change orders
    orders = {'Genotype': ['WT','rnh201', 'PolWT','Pol2M644G','Pol3L612M','Pol3L612G', 'Pol1L868M', 'Pol1Y869A'],\
              'String': ['RS','EM', 'HY']}
    for k,v in orders.items():
        for va in df[k].unique():
            if va not in v:
                v.append(va)
        df[k] = pd.Categorical(df[k], v)

    # sort
    df = df.sort_values(list(orders.keys()) + ['Library'])

    # out
    df.to_csv(args.o, sep='\t', index=False)

if __name__ == '__main__':
    main()


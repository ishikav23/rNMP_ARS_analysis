#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
from checkTimeInputs import *
from checkTimeCalcs import *

def main():

    # argparse
    parser = argparse.ArgumentParser(description='check whether ARS firing time could affect the ribonucleotide incorporation')
    parser.add_argument('ars', type=argparse.FileType('r'), help='Bed file for ars region with time')
    parser.add_argument('list', type=argparse.FileType('r'), help='List for bed file with genotype')
    parser.add_argument('-csv', help='Start from a generated dataframe csv file, skip data reading')
    parser.add_argument('-bed', default='.', help='Folder of bed file, default=\'.\'')
    parser.add_argument('-l', type=int, default=15000, help='Length of flank region, default=15000')
    parser.add_argument('-o', default='Output', help='Output file basename')
    parser.add_argument('--block_ribosomal', action='store_false',  help='Do not block ribosomal DNA')
    parser.add_argument('--efficiency', action='store_true', help='Use efficiency instead of time')
    args = parser.parse_args()

    # read lib info
    libinfo = read_libinfo(args.list)
    libs = list(libinfo.keys())
    print('Libraries:' + ','.join(libs))

    # read data
    if not args.csv:
        # read ars
        ars = read_ars(args.ars)
        print('ARS information read!')
        # extend position
        windows = generate_windows(ars, args.l, args.block_ribosomal)
        # add data
        data = read_data(windows, libs, args.bed)
        df = generate_df(data, libinfo)
        df.to_csv(args.o + '_data.csv', index=False)
    else:
        df = pd.read_csv(args.csv)
    print('Data read!')

    # generate summary
    df_summary = generate_summary(df)
    genotypes=df_summary.Genotype.unique()
    genotypes_possible = ['Rrnh201','EMrnh201','rnh201','WT']
    genotypes_used = [x for x in genotypes_possible if x in genotypes]
    # plot
    draw_ratio_scatter(df_summary, genotypes_used, output=args.o+'_MLE_scatter.png', use_efficiency=args.efficiency)
    draw_ratio_scatter(df_summary, genotypes_used, output=args.o+'_mean_scatter.png', use_MLE_ratio=False, use_efficiency=args.efficiency)

    print('Done!')


if __name__ == '__main__':
    main()


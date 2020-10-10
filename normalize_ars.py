#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from collections import defaultdict

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Normalize the ars region ')
    parser.add_argument('raw', type=argparse.FileType('r'), help='ARS ribos frequency file needed to be normalized, library information should be add so that frequency start at 9th column')
    parser.add_argument('bg', type=argparse.FileType('r'), help='Background frequency')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('--norm', default='zscore', choices=['zscore', 'prob', 'sum1'], help='Representation of normalized frequency, zscore, probabilty or sum1, default=zscore')
    parser.add_argument('--name', default='', help='Prefix of the input file, default = prefix of input')
    parser.add_argument('--notime', action='store_true', help='No time information in input file')
    parser.add_argument('--nopos', action='store_true', help='No pos information in input file')

    args = parser.parse_args()

    if args.name == '':
        args.name = args.raw.name.split('/')[-1].split('_')[0]

    # load bg frequency
    bg = defaultdict(dict)
    wsh = args.bg.readline().rstrip('\n').split('\t')
    for l in args.bg:
        ws = l.split('\t')
        for i in range(1,len(ws)):
            bg[ws[0]][wsh[i]] = float(ws[i])

    # load freqs
    freq_start = 8
    if args.notime:
        freq_start -= 1
    if args.nopos:
        freq_start -= 1

    di = args.raw.readline().rstrip('\n').split('\t')
    args.o.write('\t'.join(di[:freq_start] + ['RPB'] + di[freq_start:] + [i+'%' for i in di[freq_start:]]) + '\t')
    args.o.write('\t'.join([i+'n' for i in di[freq_start:]]) + '\n')
    for l in args.raw:
        ws = l.rstrip('\n').split('\t')

        # get background name
        norm = '_'.join([args.name] + ws[4:freq_start-1])
        if norm not in bg:
            sys.exit(f'Cannot find background information for {norm}')

        args.o.write('\t'.join(ws[:freq_start]))
        # Probability for ribos incor = count/divided by length
        p = float(ws[freq_start-1])/sum(list(bg[norm].values()))
        args.o.write(f'\t{p}\t')
        args.o.write('\t'.join(ws[freq_start:]))

        # deal with freq
        freq = list(map(float, ws[freq_start:]))
        args.o.write('\t' + '\t'.join([str(i/sum(freq)) for i in freq]))

        # calc norm freq
        freq_norm = [float(ws[i])/bg[norm][di[i]] for i in range(freq_start,len(ws))]
        if args.norm == 'prob':
            args.o.write('\t' + '\t'.join([str(i) for i in freq_norm]))
        elif args.norm == 'sum1':
            try:
                args.o.write('\t' + '\t'.join([str(i/sum(freq_norm)) for i in freq_norm]))
            except ZeroDivisionError:
                args.o.write('\t' + '\t'.join(['0.0' for i in freq_norm]))
        elif args.norm == 'zscore':
            freq_mean = np.mean(freq_norm)
            if freq_mean == 0:
                args.o.write('\t' + '\t'.join(['0.0' for i in freq_norm]))
            else:
                freq_std = np.std(freq_norm)
                args.o.write('\t' + '\t'.join([str((i-freq_mean)/freq_std) for i in freq_norm]))
        args.o.write('\n')


    print('Done!')


if __name__ == '__main__':
    main()

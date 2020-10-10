#!/usr/bin/env python3
import argparse
import sys
from collections import OrderedDict

def main():
    parser = argparse.ArgumentParser(description='Sum up bg file to generate background for ARS heatmaps')
    parser.add_argument('info', type=argparse.FileType('r'), help='ARS info file')
    parser.add_argument('-s', type=int, default=0, help='Start postion, exclude. (0)')
    parser.add_argument('-e', type=int, default=2**32, help='End position, include. (2**32)')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # header
    header = args.info.readline().rstrip('\n').split('\t')
    args.o.write('chrom\t' + '\t'.join(header[1:]) + '\n')

    # get data
    data = OrderedDict()
    for l in args.info:
        ws = l.rstrip('\n').split('\t')
        # skip unselected entries
        features = ws[0].split('_')
        pos = float(features[2])
        if pos <= args.s or pos > args.e:
            continue
        # sum up all selected entries
        name = '-'.join(features[:2] + [features[3]])
        if name not in data:
            data[name] = [float(x) for x in ws[1:]]
        else:
            for i in range(len(ws[1:])):
                data[name][i] += float(ws[1+i])

    # output
    for k,v in data.items():
        args.o.write('\t'.join([k] + [str(x) for x in v]) + '\n')

    print('Done!')





if __name__ == '__main__':
    main()
